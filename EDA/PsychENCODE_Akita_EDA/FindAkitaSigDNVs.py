import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import SubplotSpec
from collections import defaultdict
from pybedtools import BedTool
import argparse
import shutil
from datetime import date
import sys

parser = argparse.ArgumentParser(description='Compare Akita predictions for PsychENCODE GWAS againts a null distribution to find extreme variants.', prog='FindAkitaSigVariants.py')
parser.add_argument('--scoring_system', required=True, type=str,  help='Either `msd` or `max`; one of two scoring systems used to calculate variant effect scores along a sequence')
parser.add_argument('--disease', nargs = 1, required=True, type=str,  help='Name of disease(s) to score here')

args = parser.parse_args()
scoring_system = str(args.scoring_system) #either "msd" or "max"
disease = args.disease #this should be just one disease
print(disease)
assert len(disease) == 1, "This code expects to only take 1 disease as input and parallelize them over many SGE array tasks!"
disease = disease[0] #get the name of the disease out of the list
print(disease)
sys.stdout.flush()

def get_file_birthtime(filename):
    unix_timestamp = os.stat(filename).st_birthtime
    file_birthdate = datetime.fromtimestamp(unix_timestamp)
    year = file_birthdate.year
    month = file_birthdate.month
    day = file_birthdate.day
    return f"{year}-{month}-{day}"

#Load Akita Functions

import os
import json
import subprocess
import shutil
os.environ["CUDA_VISIBLE_DEVICES"] = '-1' ### run on CPU

import tensorflow as tf
print(tf.__version__)
if tf.__version__[0] == '1':
    tf.compat.v1.enable_eager_execution()

import numpy as np
import pandas as pd
import pysam
import matplotlib.pyplot as plt
from cooltools.lib.numutils import set_diag
import sys
sys.path.append("/pollard/data/projects/sdrusinsky/pollard_lab/basenji")
from basenji import dataset, dna_io, seqnn

import cooler
import cooltools
from cooltools.lib.numutils import observed_over_expected
from cooltools.lib.numutils import adaptive_coarsegrain
from cooltools.lib.numutils import interpolate_bad_singletons
from cooltools.lib.numutils import interp_nan, set_diag
from cooltools.lib.plotting import *


import kipoiseq
from kipoiseq import Interval
import pyfaidx
from collections import Counter

fasta_file = '/pollard/data/projects/sdrusinsky/pollard_lab/data/hg38_genome.fa'



### load params, specify model ###

model_dir = '/pollard/data/projects/sdrusinsky/pollard_lab/'
params_file = model_dir+'params.json'
model_file  = model_dir+'model_best.h5'
with open(params_file) as params_open:
    params = json.load(params_open)
    params_model = params['model']
    params_train = params['train']

seqnn_model = seqnn.SeqNN(params_model)


### restore model ###
# note: run %%bash get_model.sh 
# if you have not already downloaded the model
seqnn_model.restore(model_file)
print('successfully loaded')


### names of targets ###
data_dir =   '/pollard/data/projects/sdrusinsky/pollard_lab/basenji/manuscripts/akita/data/'

hic_targets = pd.read_csv(data_dir+'/targets.txt',sep='\t')
hic_file_dict_num = dict(zip(hic_targets['index'].values, hic_targets['file'].values) )
hic_file_dict     = dict(zip(hic_targets['identifier'].values, hic_targets['file'].values) )
hic_num_to_name_dict = dict(zip(hic_targets['index'].values, hic_targets['identifier'].values) )

# read data parameters
data_stats_file = '%s/statistics.json' % data_dir
with open(data_stats_file) as data_stats_open:
    data_stats = json.load(data_stats_open)
seq_length = data_stats['seq_length']
target_length = data_stats['target_length']
hic_diags =  data_stats['diagonal_offset']
target_crop = data_stats['crop_bp'] // data_stats['pool_width']
target_length1 = data_stats['seq_length'] // data_stats['pool_width']

### for converting from flattened upper-triangluar vector to symmetric matrix  ###

def from_upper_triu(vector_repr, matrix_len, num_diags):
    z = np.zeros((matrix_len,matrix_len))
    triu_tup = np.triu_indices(matrix_len,num_diags)
    z[triu_tup] = vector_repr
    for i in range(-num_diags+1,num_diags):
        set_diag(z, np.nan, i)
    return z + z.T


def preprocess_from_cool(myseq_str, genome_hic_cool):
    print("Seq-str: ", myseq_str)
    num_counts= np.sum(genome_hic_cool.matrix(balance=False).fetch(myseq_str))
    seq_hic_obs = genome_hic_cool.matrix(balance=True).fetch(myseq_str)
    seq_hic_smoothed =  adaptive_coarsegrain(
                     seq_hic_obs,  
                     genome_hic_cool.matrix(balance=False).fetch(myseq_str),  
                     cutoff=3, max_levels=8)
    seq_hic_nan = np.isnan(seq_hic_smoothed)
    seq_hic_obsexp = observed_over_expected(seq_hic_smoothed, ~seq_hic_nan)[0]
    seq_hic_obsexp = np.log(seq_hic_obsexp)
    seq_hic_obsexp = np.clip(seq_hic_obsexp,-2,2)
    seq_hic_obsexp_init = np.copy(seq_hic_obsexp)
    seq_hic_obsexp = interp_nan(seq_hic_obsexp)
    seq_hic_obsexp = np.nan_to_num(seq_hic_obsexp)
    seq_hic = np.clip(seq_hic_obsexp,-2,2)
    for i in [-1,0,1]: set_diag(seq_hic, 0,i)
        
    from astropy.convolution import Gaussian2DKernel
    from astropy.convolution import convolve
    kernel = Gaussian2DKernel(x_stddev=1,x_size=5)

    seq_hic = convolve(seq_hic, kernel)
    return seq_hic, num_counts, seq_hic_obs


def get_expt(region_chr, region_start, region_stop,genome_hic_cool,target_length_cropped):
    myseq_str = "{}:{}-{}".format(region_chr, region_start, region_stop)
    expt, num_counts, expt_obs = preprocess_from_cool(myseq_str, genome_hic_cool)
    new_start = int((target_length - target_length_cropped)/2)
    new_end = int(target_length-new_start)
    expt = expt[new_start:target_length-new_start, new_start:target_length-new_start]
    return(expt)


# @title `variant_centered_sequences`

class FastaStringExtractor:
    
    def __init__(self, fasta_file):
        self.fasta = pyfaidx.Fasta(fasta_file)
        self._chromosome_sizes = {k: len(v) for k, v in self.fasta.items()}

    def extract(self, interval: Interval, **kwargs) -> str:
        # Truncate interval if it extends beyond the chromosome lengths.
        chromosome_length = self._chromosome_sizes[interval.chrom]
        trimmed_interval = Interval(interval.chrom,
                                    max(interval.start, 0),
                                    min(interval.end, chromosome_length),
                                    )
        # pyfaidx wants a 1-based interval
        sequence = str(self.fasta.get_seq(trimmed_interval.chrom,
                                          trimmed_interval.start + 1,
                                          trimmed_interval.stop).seq).upper()
        # Fill truncated values with N's.
        pad_upstream = 'N' * max(-interval.start, 0)
        pad_downstream = 'N' * max(interval.end - chromosome_length, 0)
        return pad_upstream + sequence + pad_downstream

    def close(self):
        return self.fasta.close()


def variant_generator(vcf_file, gzipped=False):
  """Yields a kipoiseq.dataclasses.Variant for each row in VCF file."""
  def _open(file):
    return gzip.open(vcf_file, 'rt') if gzipped else open(vcf_file)
    
  with _open(vcf_file) as f:
    for line in f:
      if line.startswith('#') or line.startswith('CHROM'):
        continue
      chrom, pos, id, ref, alt_list = line.split('\t')[:5]
      # Split ALT alleles and return individual variants as output.
      for alt in alt_list.split(','):
        yield kipoiseq.dataclasses.Variant(chrom=chrom, pos=pos,
                                           ref=ref, alt=alt, id=id)


def one_hot_encode(sequence):
  return kipoiseq.transforms.functional.one_hot_dna(sequence).astype(np.float32)


def get_N_composition(seq: str):
    """
    Get % of N's in input sequence
    
    Input: 
        seq: string of sequence
    Returns: % of Ns in input sequence
    """
    count = Counter(seq)
    
    for key, value in count.items():
        count[key] = round(value/len(seq)*100,2)
#     if 'N' in count.keys():
    if count['N'] > 0:
        return count['N']
    else:
        return 0 

def variant_centered_sequences(vcf_file, sequence_length, gzipped=False,
                               chr_prefix=''):
  seq_extractor = kipoiseq.extractors.VariantSeqExtractor(
    reference_sequence=FastaStringExtractor(fasta_file))

  for variant in variant_generator(vcf_file, gzipped=gzipped):
    interval = Interval(chr_prefix + variant.chrom,
                        variant.pos, variant.pos)
    interval = interval.resize(sequence_length)
    center = interval.center() - interval.start

    reference = seq_extractor.extract(interval, [], anchor=center)
    ref_N_composition = get_N_composition(reference)
    
    alternate = seq_extractor.extract(interval, [variant], anchor=center)
    alt_N_composition = get_N_composition(alternate)
    yield {'inputs': {'ref': dna_io.dna_1hot(reference),
                      'alt': dna_io.dna_1hot(alternate)},
           'metadata': {'chrom': chr_prefix + variant.chrom,
                        'pos': variant.pos,
                        'id': variant.id,
                        'ref': variant.ref,
                        'alt': variant.alt,
                        'ref_N_composition':ref_N_composition,
                        'alt_N_composition':alt_N_composition}}
    
def msd(alternate_prediction, reference_prediction):
    
    #returns Mean squared difference between alt and ref predictions for each cell line
    return np.nanmean(np.square(alternate_prediction - reference_prediction),axis = 1).reshape(-1)

def max_diff(alternate_prediction, reference_prediction):
    #returns max difference between absolute value of alt and ref predictions for each cell line


    return np.max(abs(alternate_prediction - reference_prediction),axis = 1).reshape(-1)



def plot_variant(cell_type,hg38_chr,hg38_pos,hg38_wt,hg38_alt):
    """
    cell_type (int): 0-5. Same convenction as `hic_num_to_name_dict`
    hg38_wt (str): Nucleotide in reference
    hg38_alt: altered nucleotide
    """
    if cell_type == 0:
        cooler_file = '/pollard/data/projects/sdrusinsky/pollard_lab/basenji/manuscripts/akita/data/coolers/HFF_hg38_4DNFIP5EUOFX.mapq_30.2048.cool'
    elif cell_type == 1:
        cooler_file = '/pollard/data/projects/sdrusinsky/pollard_lab/basenji/manuscripts/akita/data/coolers/H1hESC_hg38_4DNFI1O6IL1Q.mapq_30.2048.cool'
    elif cell_type == 2:
        cooler_file = '/pollard/data/projects/sdrusinsky/pollard_lab/basenji/manuscripts/akita/data/coolers/GM12878_inSitu_MboI_all.hg38.2048.cool'

    elif cell_type == 3:
        cooler_file = '/pollard/data/projects/sdrusinsky/pollard_lab/basenji/manuscripts/akita/data/coolers/IMR90_inSitu_MboI_all.hg38.2048.cool'
    elif cell_type == 4:
        pass
    
    genome_hic_cool = cooler.Cooler(cooler_file)
    myseq_str = hg38_chr+':'+str(hg38_pos)
    
    seq_extractor = kipoiseq.extractors.VariantSeqExtractor(
    reference_sequence=FastaStringExtractor(fasta_file))
    
    #get predictions
    interval = Interval(hg38_chr,hg38_pos, hg38_pos)
    interval = interval.resize(seq_length)
    center = interval.center() - interval.start
    variant = kipoiseq.dataclasses.Variant(chrom = hg38_chr, pos=hg38_pos,
                                           ref=hg38_wt, alt=hg38_alt)

    reference = seq_extractor.extract(interval, [], anchor=center)
    alternate = seq_extractor.extract(interval, [variant], anchor=center)
    
    ref_1hot = dna_io.dna_1hot(reference)
    alt_1hot = dna_io.dna_1hot(alternate)
    
    ref_pred = seqnn_model.model.predict(np.expand_dims(ref_1hot,0))
    alt_pred = seqnn_model.model.predict(np.expand_dims(alt_1hot,0))
    
    print(f'ln MSD: {np.log(msd(alt_pred,ref_pred))}')
    print(f'ln Max Difference: {np.log(max_diff(alt_pred,ref_pred))}')
    

    #experimental map for this region
    hic_params = params['model']['head_hic']
    cropping = hic_params[5]['cropping']
    target_length = params_model['target_length']
    print(target_length)
    target_length_cropped = target_length - 2 * cropping
    target = get_expt(hg38_chr, hg38_pos - (seq_length//2), hg38_pos + (seq_length//2),genome_hic_cool,target_length_cropped) # get experimental data
    print(target)
    
    
    target_index = cell_type

    print(' ')
    print(myseq_str)

    


    plt.figure(figsize=(8,4))
    vmin=-2; vmax=2

    # plot pred
    plt.subplot(121) 
    mat = from_upper_triu(ref_pred[:,:,target_index], target_length_cropped, hic_diags)
    im = plt.matshow(mat, fignum=False, cmap= 'RdBu_r', vmax=vmax, vmin=vmin)
    plt.colorbar(im, fraction=.04, pad = 0.05, ticks=[-2,-1, 0, 1,2]);
    plt.title('Predicted Reference-'+str(hic_num_to_name_dict[target_index]),y=1.15 )
    plt.ylabel(myseq_str)

    # # plot target 
    plt.figure(figsize=(8,4))
    plt.subplot(122) 
#     mat = from_upper_triu(test_target, target_length1_cropped, hic_diags)
    im = plt.matshow(target, fignum=False, cmap= 'RdBu_r', vmax=vmax, vmin=vmin)
    plt.colorbar(im, fraction=.04, pad = 0.05, ticks=[-2,-1, 0, 1,2]);
    plt.title( 'Experimental Reference-'+str(hic_num_to_name_dict[target_index]),y=1.15)
    plt.tight_layout()
    plt.show()
    plt.close()

    #plot mutation
    plt.figure(figsize=(8,4))
    plt.subplot(221)
    mat = from_upper_triu(alt_pred[:,:,target_index], target_length_cropped, hic_diags)
    im = plt.matshow(mat, fignum=False, cmap= 'RdBu_r', vmax=vmax, vmin=vmin)
    plt.colorbar(im, fraction=.04, pad = 0.05, ticks=[-2,-1, 0, 1,2]);
    plt.title( f'{hg38_chr} {hg38_pos} {hg38_wt}>{hg38_alt}-'+str(hic_num_to_name_dict[target_index]),y=1.15)

    
    
    #plot difference between ref and alt
    
    plt.figure(figsize=(8,4))
    plt.subplot(222)
    mat = from_upper_triu(alt_pred[:,:,target_index] - ref_pred[:,:,target_index], target_length_cropped, hic_diags)
    im = plt.matshow(mat, fignum=False, cmap= 'RdBu_r', vmax=vmax, vmin=vmin)
    plt.colorbar(im, fraction=.04, pad = 0.05, ticks=[-2,-1, 0, 1,2]);
    plt.title( f'{hg38_chr} {hg38_pos} {hg38_wt}>{hg38_alt} ALT - REF -'+str(hic_num_to_name_dict[target_index]),y=1.15)

    
    
    plt.tight_layout()
    plt.show()

    
def msd(alternate_prediction, reference_prediction):

    #returns Mean squared difference between alt and ref predictions for each cell line
    return np.nanmean(np.square(alternate_prediction - reference_prediction),axis = 1).reshape(-1)

def max_diff(alternate_prediction, reference_prediction):
    #returns max difference between absolute value of alt and ref predictions for each cell line


    return np.max(abs(alternate_prediction - reference_prediction),axis = 1).reshape(-1)




########## Load null distribution ############
null_path = "/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/Akita/Rare1KGVariants"
null_dist = pd.DataFrame()
for idx, file in enumerate(os.listdir(null_path)):
  if scoring_system in file and file.endswith('.vcf.csv'):
    split_null_df = pd.read_csv(os.path.join(null_path,file))
    null_dist = null_dist.append(split_null_df)

score_cols = ['0_HFF', '1_H1hESC', '2_GM12878', '3_IMR90','4_HCT116']
######## Load PsychENCODe GWAS variant scores #######
akita_dir = "/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/Akita"
disease_akita_scores = pd.DataFrame()

PsychENCODE_disease_akita_path = os.path.join(akita_dir,disease)
for file in os.listdir(PsychENCODE_disease_akita_path):
    if file.endswith('.csv') and f'{scoring_system}' in file: #only open the csv files with the variant score (either max or summed, depending on the input argument to this file)
        PsychENCODE_disease_akita_var_df = pd.read_csv(os.path.join(akita_dir,disease,file))
                                                                   
        PsychENCODE_disease_akita_var_df['disease'] = disease
        disease_akita_scores = disease_akita_scores.append(PsychENCODE_disease_akita_var_df)
disease_akita_scores = disease_akita_scores[disease_akita_scores['chrom']!= 'chrom'] #some rows have the header inside. remove them
disease_akita_scores[score_cols] = disease_akita_scores[score_cols].astype(float) #coerce the akita scores to be float, since they were strings before due to header being inside the rows
disease_akita_scores[score_cols] = np.log(disease_akita_scores[score_cols]) #upper percentiles will be returned ln transformed, so apply that to the PsychENCODE Akita scores, else everything will appear significant
def get_quantiles(threshold, null_distribution,score_cols):

    
    upper_percentiles = defaultdict() #only looking for big differences between ref and alt. No need for lower percentiles; not a matter of upregulation or downregulation, but is the contact highly different or not?
    for score_col in score_cols:
        null_distribution[score_col] = np.log(null_distribution[score_col])
        upper_percentiles[score_col] = null_distribution[score_col].quantile(threshold)
    return upper_percentiles 
    
upper_percentiles = get_quantiles(0.999,null_dist,score_cols)



def get_sig_variants(GWAS_variant_predictions,upper_percentiles,score_cols):
    """
    GWAS_variant_predictions: df with variant effect scores for each lead and tag SNP in the GWAS of interest
    upper_percentiles: Dictionary with upper threshold for each null distribution above which a variant effect is significant
    score_cols: columns in GWAS_variant_predictions corresponding to variant effect scores
    
    """

    sig_variants = pd.DataFrame()
    variant_info = ['chrom','pos','ref','alt','disease']
    for score_col in score_cols:
        #append variants to df if they are more extreme than the upper or lower threshold
        sig_variants = sig_variants.append(GWAS_variant_predictions[
                                    (GWAS_variant_predictions[score_col] >= upper_percentiles[score_col])][variant_info + [score_col]].copy())

    #group variants that are extreme against different null distributions back to the same row
    sig_variants = sig_variants.groupby(['chrom','pos','ref','alt']).max().reset_index()
    
    return sig_variants

sig_variants = get_sig_variants(disease_akita_scores,upper_percentiles,score_cols)


filename = f"/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/PsychENCODE_GWAS_scripts/EDA/PsychENCODE_Akita_EDA/AkitaSig_{scoring_system}_{disease}_DNVs.csv"
if os.path.exists(filename): #if there is already data, back it up before re-doing this analysis
    #make backup of existing sig_variants
    backup_dir = f"/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/PsychENCODE_GWAS_scripts/EDA/PsychENCODE_Akita_EDA/backups/{disease}"
    if not os.path.exists(backup_dir):
        os.mkdir(backup_dir)
    backup_path = os.path.join(backup_dir,f"AkitaSig_{scoring_system}_{disease}_DNVs_{get_file_birthtime(filename)}.csv")
    shutil.copy2(filename,backup_path)

sig_variants.to_csv(filename)





def get_pvals(sig_variants,null_distribution,score_cols):
    """
    take variants already shown to be significant (to save computation time) and calculate p value for each
    
    """

    sig_variants_pval = sig_variants.copy()

    for idx,row in sig_variants.iterrows():
        #get seq classes this sig variant was significant (i.e not null) in 
        sig_score_cols = row[score_cols][row[score_cols].notnull()].index #the index here is actually the columns of the dataframe because it's one series
        values = row[score_cols][row[score_cols].notnull()].values


        #iterate through all seq classes this sig variant was significant in and get actual p value using the null distribution
        for idx2,score_col in enumerate(sig_score_cols):


            #get p value, calculated as the number of variants in the null distribution at least as extreme as this one
            #only interested in one tail (right tail)
            pval = sum(null_distribution[score_col] >= values[idx2]) / null_distribution.shape[0]


            sig_variants_pval.loc[idx,score_col] = pval #record pval for this variant
    return sig_variants_pval

sig_variants_pval = get_pvals(sig_variants,null_dist,score_cols)

filename = f"/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/PsychENCODE_GWAS_scripts/EDA/PsychENCODE_Akita_EDA/AkitaSigPVals_{scoring_system}_{disease}_DNVs.csv"
if os.path.exists(filename): #if there is already data, back it up before re-doing this analysis
    #make backup of existing sig_variants
    backup_dir = f"/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/PsychENCODE_GWAS_scripts/EDA/PsychENCODE_Akita_EDA/backups/{disease}"
    if not os.path.exists(backup_dir):
        os.mkdir(backup_dir)
    backup_path = os.path.join(backup_dir,f"AkitaSigPVals_{scoring_system}_{disease}_DNVs_{get_file_birthtime(filename)}.csv")
    shutil.copy2(filename,backup_path)


sig_variants_pval.to_csv(filename)
