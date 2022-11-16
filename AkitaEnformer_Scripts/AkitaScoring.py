#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import json
import subprocess
import shutil

import tensorflow as tf
print(tf.__version__)
if tf.__version__[0] == '1':
    tf.compat.v1.enable_eager_execution()
from tensorflow.python.client import device_lib
import numpy as np
import pandas as pd
import pysam
import matplotlib.pyplot as plt
from cooltools.lib.numutils import set_diag
import sys
sys.path.append("./models/Akita/basenji")
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
import argparse
import sys
from pathlib import Path

# for GPUs
# print(device_lib.list_local_devices())

# assert tf.test.is_built_with_cuda(),"The installed version of Tensorflow does not include GPU support"
# print(tf.config.list_physical_devices())
# assert tf.config.list_physical_devices('GPU'), "Connect to GPU!" #make sure GPU is enabled

#make sure tensorflow doesn't use all the memory on the GPU



parser = argparse.ArgumentParser(description='Run Akita on variants associated with diseases (via the GWAS catalog) relevant to the PsychENCODE project.', prog='AkitaPsychENCODE_GWAS_predictions.py')
parser.add_argument('--vcf_path', required=True, type=str,  help='Full path to the vcf file containing the mutations to be run by Akita')
parser.add_argument('--reference_genome', required=True, type=str,  help='Either hg38 or hg19')
parser.add_argument('--experiment_name', required=True, type=str,  help='Name of the experiment')

args = parser.parse_args()
vcf_path = str(args.vcf_path)
experiment_name = str(args.experiment_name)
reference_genome = str(args.reference_genome)
reference_genome = reference_genome.lower()
print(f"Reference genome: {reference_genome}")
assert reference_genome in ['hg38','hg19'], "Desired reference genome is not supported!"

### load params, specify model ###
model_dir = '../models/Akita/basenji/manuscripts/akita'
params_file = model_dir+'params.json'
model_file  = model_dir+'model_best.h5'
with open(params_file) as params_open:
    params = json.load(params_open)
    params_model = params['model']
    params_train = params['train']

seqnn_model = seqnn.SeqNN(params_model)


# In[3]:


### restore model ###
# note: run %%bash get_model.sh 
# if you have not already downloaded the model
seqnn_model.restore(model_file)
print('successfully loaded')


# In[4]:


### names of targets ###
data_dir =  os.path.join(model_dir,"basenji/manuscripts/akita/data/")

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


if reference_genome == 'hg38':
  fasta_file = '/pollard/data/projects/sdrusinsky/pollard_lab/data/hg38_genome.fa'
elif reference_genome == 'hg19':
  fasta_file = "/pollard/data/projects/sdrusinsky/pollard_lab/data/hg19.fa"




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


def get_expt(region_chr, region_start, region_stop):
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
    return gzip.open(vcf_file, 'rt') if gzipped else open(vcf_file,encoding='utf-8')
    
  with _open(vcf_file) as f:
    for line in f:
      if line.startswith('#') or line.startswith('CHROM') or line.startswith('    '): # header lines start wth either '#', "CHROM" to denote the colnames row, or '    ' because of a formatting error when generating the headers for the vcf
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


# In[7]:
#take vvariants and get model(alt) - model(ref) predictions
#take MSD or max along sequence axis to get variant score for each track. Save these scores + variant position and allele metadata
output_dir = f"../model_outputs/Akita/{experiment_name}"
if not os.path.exists(output_dir):
  os.mkdir(output_dir)
vcf_basename = vcf_path.split('/')[-1]
it = variant_centered_sequences(vcf_path, sequence_length=seq_length,
                            gzipped=False, chr_prefix='')
msd_scores_list = []
maxed_scores_list = []
for idx, example in enumerate(it):
    reference_prediction = seqnn_model.model.predict({k: np.expand_dims(v,0) for k,v in example['inputs'].items()}['ref'])
    alternate_prediction = seqnn_model.model.predict({k: np.expand_dims(v,0) for k,v in example['inputs'].items()}['alt'])
    msd_vals = msd(alternate_prediction, reference_prediction)
    max_diffs = max_diff(alternate_prediction, reference_prediction)
    
  

    msd_scores = {f'{i}_{name[:20]}': score for i, (name, score) in enumerate(zip(hic_num_to_name_dict.values(),msd_vals))}
    msd_scores_list.append({**example['metadata'],
                               **msd_scores})
    maxed_scores ={f'{i}_{name[:20]}': score for i, (name, score) in enumerate(zip(hic_num_to_name_dict.values(),max_diffs))}
    maxed_scores_list.append({**example['metadata'],
                               **maxed_scores})
    

    
    
    if idx == 0:
        #save first variant with header to csv
        max_df = pd.DataFrame(maxed_scores_list)
        max_df.to_csv(f"{scratch_path}/Akita_max_predictions_PsychENCODEGWASLeadTagVariants_{vcf_basename}.csv",index=False,header=True)
        maxed_scores_list = []


        msd_df = pd.DataFrame(msd_scores_list)
        msd_df.to_csv(f"{scratch_path}/Akita_msd_predictions_PsychENCODEGWASLeadTagVariants_{vcf_basename}.csv",index=False,header=True)
        msd_scores_list = []

    elif idx != 0 and idx % 500 == 0:
        #append scores every 500 (the lists are being appended to in the for loop). Do not add header as row name
        max_df = pd.DataFrame(maxed_scores_list)
        max_df.to_csv(f"{scratch_path}/Akita_max_predictions_PsychENCODEGWASLeadTagVariants_{vcf_basename}.csv",mode = 'a',index=False,header=False)
        maxed_scores_list = []


        msd_df = pd.DataFrame(msd_scores_list)
        msd_df.to_csv(f"{scratch_path}/Akita_msd_predictions_PsychENCODEGWASLeadTagVariants_{vcf_basename}.csv",mode = 'a',index=False,header=False)
        msd_scores_list = []
        print(idx)

global_path = f"/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/Akita/{disease}" #where final data will be moved to
Path(global_path).mkdir(parents=True, exist_ok=True) #just in case the folder where the data will be stored on the Pollard server does not exist yet, create it

#save the last set of SNPs when idx % 500 != 0. i.e, when total number of variants isn't divisable by 500.
max_df = pd.DataFrame(maxed_scores_list)
max_df.to_csv(f"{scratch_path}/Akita_max_predictions_PsychENCODEGWASLeadTagVariants_{vcf_basename}.csv",mode = 'a',index=False,header=False)
shutil.move(f"{scratch_path}/Akita_max_predictions_PsychENCODEGWASLeadTagVariants_{vcf_basename}.csv",f"{global_path}/Akita_max_predictions_PsychENCODEGWASLeadTagVariants_{vcf_basename}.csv") #move from scratch to global location on Pollard NAS


msd_df = pd.DataFrame(msd_scores_list)
msd_df.to_csv(f"{scratch_path}/Akita_msd_predictions_PsychENCODEGWASLeadTagVariants_{vcf_basename}.csv",mode = 'a',index=False,header=False)
shutil.move(f"{scratch_path}/Akita_msd_predictions_PsychENCODEGWASLeadTagVariants_{vcf_basename}.csv",f"{global_path}/Akita_msd_predictions_PsychENCODEGWASLeadTagVariants_{vcf_basename}.csv") #move from scratch to global location on Pollard NAS


print('done')




