#!/usr/bin/env python
# coding: utf-8

# In[1]:
import os
import tensorflow as tf
import gzip
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
from scipy.stats import zscore
from pandas import HDFStore
import h5py
import itertools
import argparse
import sys
import datetime
import shutil
from pathlib import Path
from tensorflow.python.client import device_lib
import git


## For using GPU
# print(device_lib.list_local_devices())
# assert tf.test.is_built_with_cuda(),"The installed version of Tensorflow does not include GPU support"
# print(tf.config.list_physical_devices())
# assert tf.config.list_physical_devices('GPU'), "Connect to GPU!" #make sure GPU is enabled



parser = argparse.ArgumentParser(description='Run Enformer on variants associated with diseases (via the GWAS catalog) relevant to the PsychENCODE project.', prog='EnformerPsychENCODE_GWAS_predictions.py')
parser.add_argument('--vcf_path', required=True, type=str,  help='Full path to the vcf file containing the mutations to be run by Enformer')
parser.add_argument('--reference_genome', required=True, type=str,  help='Either hg38 or hg19')
parser.add_argument('--experiment_name', required=True, type=str,  help='Name of the experiment')

args = parser.parse_args()
vcf_path = str(args.vcf_path)
experiment_name = str(args.experiment_name)
reference_genome = str(args.reference_genome)
reference_genome = reference_genome.lower()
assert reference_genome in ['hg38','hg19'], "Desired reference genome is not supported!"
#get top level directory of the git repository:
git_repo = git.Repo(os.getcwd(), search_parent_directories=True)
git_root = git_repo.git.rev_parse("--show-toplevel") #this is path/to/PCGC_Challenge
fasta_file = f'{git_root}/models/{reference_genome}_genome.fa'

model_path = 'https://tfhub.dev/deepmind/enformer/1'








# Download targets from Basenji2 dataset 
# Cite: Kelley et al Cross-species regulatory sequence activity prediction. PLoS Comput. Biol. 16, e1008050 (2020).
targets_txt = 'https://raw.githubusercontent.com/calico/basenji/master/manuscripts/cross2020/targets_human.txt'
df_targets = pd.read_csv(targets_txt, sep='\t')



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

  def __init__(self, tfhub_url, transform_pkl_path,
               organism='human'):
    assert organism == 'human', 'Transforms only compatible with organism=human'
    self._model = EnformerScoreVariantsRaw(tfhub_url, organism)
    with tf.io.gfile.GFile(transform_pkl_path, 'rb') as f:
      transform_pipeline = joblib.load(f)
    self._transform = transform_pipeline.steps[0][1]  # StandardScaler.
    
  def predict_on_batch(self, inputs):
    scores = self._model.predict_on_batch(inputs)
    return self._transform.transform(scores)


class EnformerScoreVariantsPCANormalized:

  def __init__(self, tfhub_url, transform_pkl_path,
               organism='human', num_top_features=500):
    self._model = EnformerScoreVariantsRaw(tfhub_url, organism)
    with tf.io.gfile.GFile(transform_pkl_path, 'rb') as f:
      self._transform = joblib.load(f)
    self._num_top_features = num_top_features
    
  def predict_on_batch(self, inputs):
    scores = self._model.predict_on_batch(inputs)
    return self._transform.transform(scores)[:, :self._num_top_features]


# TODO(avsec): Add feature description: Either PCX, or full names.


# In[51]:


# @title `variant_centered_sequences`

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


def variant_generator(vcf_file, gzipped=False,skip_lines = 0,max_lines = np.inf):
  """Yields a kipoiseq.dataclasses.Variant for each row in VCF file."""
  def _open(file):
    return gzip.open(vcf_file, 'rt') if gzipped else open(vcf_file,encoding='utf-8')
    
  with _open(vcf_file) as f:
    for idx,line in enumerate(f):
      if line.startswith('#') or line.startswith('CHROM') or line.startswith('    '): # header lines start wth either '#', "CHROM" to denote the colnames row, or '    ' because of a formatting error when generating the headers for the vcf
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
# In[63]:


model = Enformer(model_path)
fasta_extractor = FastaStringExtractor(fasta_file)


# In[55]:



#take vvariants and get model(alt) - model(ref) predictions
#take sum or max along sequence axis to get variant score for each track. Save these scores + variant position and allele metadata
output_dir = f"{git_repo}/model_outputs/Enformer/{experiment_name}"
if not os.path.exists(output_dir):
  os.mkdir(output_dir)
vcf_basename = vcf_path.split('/')[-1]

it = variant_centered_sequences(vcf_path, sequence_length=SEQUENCE_LENGTH,
                            gzipped=False, chr_prefix='')
summed_scores_list = []
maxed_scores_list = []
for idx, example in enumerate(it):
  reference_prediction = model.predict_on_batch({k: v[tf.newaxis] for k,v in example['inputs'].items()}['ref'])['human'][0]
  alternate_prediction = model.predict_on_batch({k: v[tf.newaxis] for k,v in example['inputs'].items()}['alt'])['human'][0]
  pred = alternate_prediction - reference_prediction
  summed_scores = np.sum(pred,axis = 0) #sum allelic differences along sequence axis
  maxed_scores = np.max(abs(pred),axis = 0) #take max absolute difference along sequence axis
  
  
  summed_scores = {f'{i}_{name[:20]}': score for i, (name, score) in enumerate(zip(df_targets.description, summed_scores))}
  summed_scores_list.append({**example['metadata'],
                             **summed_scores})
  maxed_scores = {f'{i}_{name[:20]}': score for i, (name, score) in enumerate(zip(df_targets.description, maxed_scores))}
  maxed_scores_list.append({**example['metadata'],
                               **maxed_scores})
   


#save variant effect scores to csv file
max_df = pd.DataFrame(maxed_scores_list)
max_df.to_csv(f"{output_dir}/Enformer_max_predictions_{vcf_basename}.csv",index=False,header=True)


summed_df = pd.DataFrame(summed_scores_list)
summed_df.to_csv(f"{output_dir}/Enformer_summed_predictions_{vcf_basename}.csv",index=False,header=True)






print('done')






