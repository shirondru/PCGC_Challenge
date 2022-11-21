import pandas as pd
import os
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
import tensorflow as tf

import tensorflow_hub as hub
import joblib
import gzip
import kipoiseq
from kipoiseq import Interval
import pyfaidx
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import argparse
import sys
from scipy.stats import zscore
from datetime import date, datetime
import pickle
import shutil


parser = argparse.ArgumentParser(description='Compare Enformer predictions for PsychENCODE GWAS againts a null distribution to find extreme variants.', prog='FindEnformerSigDNVs.py')
parser.add_argument('--scoring_system', required=True, type=str,  help='Either `summed` or `max`; one of two scoring systems used to calculate variant effect scores along a sequence')
parser.add_argument('--relevant_cols', required=True, type=str,  help='Path to text file containing relevant Enformer features for use with this disease')
parser.add_argument('--disease', nargs = 1, required=True, type=str,  help='Name of disease(s) to score here')

args = parser.parse_args()
scoring_system = str(args.scoring_system) #either "summed" or "max"
disease = args.disease #this should be just one disease
print(disease)
assert len(disease) == 1, "This code expects to only take 1 disease as input and parallelize them over many SGE array tasks!"
disease = disease[0] #get the name of the disease out of the list
print(disease)

relevant_cols_file = str(args.relevant_cols)
print(relevant_cols_file)
relevant_cols_path="/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/PsychENCODE_GWAS_scripts/EDA/PsychENCODE_Enformer_EDA"
with open(os.path.join(relevant_cols_path,relevant_cols_file),'r') as f:
  file = f.read()
  relevant_cols = file.split('\n')
relevant_cols = [x for x in relevant_cols if x!= ''] #remove emptry strings
relevant_cols = [x.lower() for x in relevant_cols] #coerce all strings to be lowercase
sys.stdout.flush()

def get_file_birthtime(filename):
    unix_timestamp = os.stat(filename).st_mtime
    file_birthdate = datetime.fromtimestamp(unix_timestamp)
    year = file_birthdate.year
    month = file_birthdate.month
    day = file_birthdate.day
    return f"{year}-{month}-{day}"


model_path = "/pollard/data/projects/sdrusinsky/pollard_lab/enformer"
fasta_file = '/pollard/data/projects/sdrusinsky/pollard_lab/data/hg38_genome.fa'
# fasta_file = "/pollard/data/projects/sdrusinsky/pollard_lab/data/hg19.fa"

clinvar_vcf = '/pollard/data/projects/sdrusinsky/pollard_lab/data/clinvar.vcf.gz'

# cols = ['chrom','txStart','txEnd','ENST','strand',,'cdsStart','cdsEnd','exonCount','exonStarts','exonEnds','ENST_y','A','B','C','D','Gene','E','F']
gene_annotations = pd.read_csv("/pollard/data/projects/sdrusinsky/pollard_lab/data/knownGene.tsv",sep = '\t') #from ucsc genome browser hg38


# Download targets from Basenji2 dataset 
# Cite: Kelley et al Cross-species regulatory sequence activity prediction. PLoS Comput. Biol. 16, e1008050 (2020).
df_targets = pd.read_csv("/pollard/data/projects/sdrusinsky/pollard_lab/data/enformer_df_targets.csv")
df_targets = df_targets.assign(assay = df_targets['description'].str.split(':',expand=True)[0])

#### LOAD ENFORMER #####

# @title `Enformer`, `EnformerScoreVariantsNormalized`, `EnformerScoreVariantsPCANormalized`,
SEQUENCE_LENGTH = 393216

class Enformer:

  def __init__(self, tfhub_url):
    self._model = hub.load(tfhub_url).model

  def predict_on_batch(self, inputs):
    predictions = self._model.predict_on_batch(inputs)
    return {k: v.numpy() for k, v in predictions.items()}

  @tf.function
  def contribution_input_grad(self, input_sequence,
                              target_mask, output_head='human'):
    input_sequence = input_sequence[tf.newaxis]

    target_mask_mass = tf.reduce_sum(target_mask)
    with tf.GradientTape() as tape:
      tape.watch(input_sequence)
      prediction = tf.reduce_sum(
          target_mask[tf.newaxis] *
          self._model.predict_on_batch(input_sequence)[output_head]) / target_mask_mass

    input_grad = tape.gradient(prediction, input_sequence) * input_sequence
    input_grad = tf.squeeze(input_grad, axis=0)
    return tf.reduce_sum(input_grad, axis=-1)


class EnformerScoreVariantsRaw:

  def __init__(self, tfhub_url, organism='human'):
    self._model = Enformer(tfhub_url)
    self._organism = organism
  
  def predict_on_batch(self, inputs):
    ref_prediction = self._model.predict_on_batch(inputs['ref'])[self._organism]
    alt_prediction = self._model.predict_on_batch(inputs['alt'])[self._organism]

    return alt_prediction.mean(axis=1) - ref_prediction.mean(axis=1)


class EnformerScoreVariantsNormalized:

  def __init__(self, tfhub_url, transform_model,
               organism='human'):
    assert organism == 'human', 'Transforms only compatible with organism=human'
    self._model = EnformerScoreVariantsRaw(tfhub_url, organism)
    
    self._transform = transform_model.steps[0][1]  # StandardScaler.
    
  def predict_on_batch(self, inputs):
    scores = self._model.predict_on_batch(inputs)
    return self._transform.transform(scores)


class EnformerScoreVariantsPCANormalized:

  def __init__(self, tfhub_url, transform_model,
               organism='human', num_top_features=500):
    self._model = EnformerScoreVariantsRaw(tfhub_url, organism)
    self.transform_model = transform_model
    self._num_top_features = num_top_features
    
  def predict_on_batch(self, inputs):
    scores = self._model.predict_on_batch(inputs)
    return self.transform_model.transform(scores)[:, :self._num_top_features]


# TODO(avsec): Add feature description: Either PCX, or full names.



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
    alternate = seq_extractor.extract(interval, [variant], anchor=center)

    yield {'inputs': {'ref': one_hot_encode(reference),
                      'alt': one_hot_encode(alternate)},
           'metadata': {'chrom': chr_prefix + variant.chrom,
                        'pos': variant.pos,
                        'id': variant.id,
                        'ref': variant.ref,
                        'alt': variant.alt}}
    
    
    
# @title `plot_tracks`

def plot_tracks(tracks, interval, height=1.5):
    fig, axes = plt.subplots(len(tracks)+1, 1, figsize=(20, height * len(tracks)), sharex=True)
    for idx, (ax, (title, y)) in enumerate(zip(axes, tracks.items())):
        ax.fill_between(np.linspace(interval.start, interval.end, num=len(y)), y)
        ax.set_title(title)
        sns.despine(top=True, right=True, bottom=True)
    if idx == len(tracks):
        ax.bar(np.linspace(interval.start, interval.end, num=len(y)), y)

    ax.set_xlabel(str(interval))
    plt.tight_layout()


model = Enformer(model_path)

fasta_extractor = FastaStringExtractor(fasta_file)



#The column names in the Enformer predictions were truncated. Harmonize this
enformer_relevant_column_idxs = df_targets[df_targets['description'].str.lower().isin(relevant_cols)]['index'].tolist()
null_path = "/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/Enformer/Rare1KGVariants"
null = pd.DataFrame()
file_num = 0
for file in os.listdir(null_path):
  if scoring_system in file and file.endswith('.vcf.csv'):
    if file_num == 0: #open just the first few rows of the first results file so you can get the colnames you want, and open the full files faster
      split_null_df_small = pd.read_csv(os.path.join(null_path,file),nrows=5) #just use 5 rows to open the file quickly and get the info for the columns you need
      enformer_relevant_cols = list(split_null_df_small.columns[5:][enformer_relevant_column_idxs])
      file_num+=1 
    split_null_df = pd.read_csv(os.path.join(null_path,file),usecols = ['chrom','pos','ref','alt'] + enformer_relevant_cols)
    null = null.append(split_null_df)
# get 99.9th percentile for each sequence class and plot value per sequence class

def get_quantiles(threshold, null_distribution,score_cols):
    upper_percentiles = defaultdict()
    lower_percentiles = defaultdict()
    for score_col in score_cols:
        upper_percentiles[score_col] = null_distribution[score_col].quantile(threshold)
        lower_percentiles[score_col] = null_distribution[score_col].quantile(1 - threshold)
    return upper_percentiles, lower_percentiles
    
    
upper_percentiles, lower_percentiles = get_quantiles(0.999,null,enformer_relevant_cols)


def get_sig_variants(GWAS_variant_predictions,upper_percentiles,lower_percentiles,score_cols):
    """
    GWAS_variant_predictions: df with variant effect scores for each lead and tag SNP in the GWAS of interest
    upper_percentiles: Dictionary with upper threshold for each null distribution above which a variant effect is significant
    lower_percentiles: Dictionary with lower threshold for each null distribution below which a variant effect is significant
    score_cols: columns in GWAS_variant_predictions corresponding to variant effect scores
    """
    sig_variants = pd.DataFrame()
    variant_info = ['chrom','pos','ref','alt','disease']
    for score_col in score_cols:
        #append variants to df if they are more extreme than the upper or lower threshold
        sig_variants = sig_variants.append(GWAS_variant_predictions[
                                    (GWAS_variant_predictions[score_col] >= upper_percentiles[score_col]) | 
                                    (GWAS_variant_predictions[score_col] <= lower_percentiles[score_col])][variant_info + [score_col]].copy())
    #group variants that are extreme against different null distributions back to the same row
    sig_variants = sig_variants.groupby(['chrom','pos','ref','alt']).max().reset_index()
    return sig_variants




#load PsychENCODE disease variant scores (for all diseases) into one dataframe. Then return significant variants
enformer_dir = "/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/Enformer"
disease_enformer_scores = pd.DataFrame()
PsychENCODE_disease_enformer_path = os.path.join(enformer_dir,disease)
for file in os.listdir(PsychENCODE_disease_enformer_path):
    if file.endswith('.csv') and f'{scoring_system}' in file: #only open the csv files with the variant score (either max or summed, depending on the input argument to this file)
        PsychENCODE_disease_enformer_var_df = pd.read_csv(os.path.join(enformer_dir,disease,file),
                                                                   usecols = ['chrom','pos','ref','alt'] + enformer_relevant_cols)
        PsychENCODE_disease_enformer_var_df['disease'] = disease
        disease_enformer_scores = disease_enformer_scores.append(PsychENCODE_disease_enformer_var_df)
sig_variants = get_sig_variants(disease_enformer_scores,upper_percentiles,lower_percentiles,enformer_relevant_cols)                                   

filename = f"/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/PsychENCODE_GWAS_scripts/EDA/PsychENCODE_Enformer_EDA/EnformerSig_{scoring_system}_{disease}_DNVs.csv"
if os.path.exists(filename): #if there is already data, back it up before re-doing this analysis
  #make backup of existing sig_variants
  backup_dir = f"/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/PsychENCODE_GWAS_scripts/EDA/PsychENCODE_Enformer_EDA/backups/{disease}"
  if not os.path.exists(backup_dir):
    os.mkdir(backup_dir)
  backup_path = os.path.join(backup_dir,f"EnformerSig_{scoring_system}_{disease}_DNVs_{get_file_birthtime(filename)}.csv")
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
            #check if the variant effect score for this variant in this class is positive or negative first to determine the tail the comparison should be made with
            if values[idx2] > 0: 
                pval = sum(null_distribution[score_col] >= values[idx2]) / null_distribution.shape[0]

            elif values[idx2] < 0:
                pval = sum(null_distribution[score_col] <= values[idx2]) / null_distribution.shape[0]

            sig_variants_pval.loc[idx,score_col] = pval #record pval for this variant
    return sig_variants_pval






sig_variants_pval = get_pvals(sig_variants,null,enformer_relevant_cols)

filename = f"/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/PsychENCODE_GWAS_scripts/EDA/PsychENCODE_Enformer_EDA/EnformerSigPVals_{scoring_system}_{disease}_DNVs.csv"
if os.path.exists(filename): #if there is already data, back it up before re-doing this analysis
  #make backup of existing sig_variants
  backup_dir = f"/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/PsychENCODE_GWAS_scripts/EDA/PsychENCODE_Enformer_EDA/backups/{disease}"
  if not os.path.exists(backup_dir):
    os.mkdir(backup_dir)
  backup_path = os.path.join(backup_dir,f"EnformerSigPVals_{scoring_system}_{disease}_DNVs_{get_file_birthtime(filename)}.csv")
  shutil.copy2(filename,backup_path)


sig_variants_pval.to_csv(filename)

