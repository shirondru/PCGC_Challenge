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


import git
#get top level directory of the git repository:
git_repo = git.Repo(os.getcwd(), search_parent_directories=True)
git_root = git_repo.git.rev_parse("--show-toplevel") #this is path/to/PCGC_Challenge


sys.path.append(f"{git_root}/models/Akita/basenji")
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



parser = argparse.ArgumentParser(description='Run Akita to score variants.', prog='AkitaScoring.py')
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
fasta_file = f'{git_root}/models/{reference_genome}_genome.fa'



### load params, specify model ###
model_dir = f'{git_root}/models/Akita/basenji/manuscripts/akita/'
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




### load more Akita info to targets ###
data_dir =  os.path.join(model_dir,"data")

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


#for parsing the reference genome, in order to find the DNA nucleotides that flank each variant
#these flanking nucleotides will be used to form a sequence of the expected length for input into Akita, with the variant at the center
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


#for parsing VCF file, and then scoring each variant one by one
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

#N's might be added to pad the input sequence. Akita's output will tell you if that was the case for any variant so you know.
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

#calls the variant_generator, converts each variant to a DNA sequence of the expected length with the variant at the center, 
#and one hot encodes that sequence. This is done to the reference allele and alternate allele. This will be passed into akita to form predictions
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
#take variants and get model(alt) - model(ref) predictions
#take MSD or max along sequence axis to get variant score for each track. Save these scores + variant position and allele metadata
output_dir = f"{git_root}/model_outputs/Akita/{experiment_name}"
if not os.path.exists(output_dir):
  os.mkdir(output_dir)
vcf_basename = vcf_path.split('/')[-1]
it = variant_centered_sequences(vcf_path, sequence_length=seq_length,
                            gzipped=False, chr_prefix='')
msd_scores_list = []
maxed_scores_list = []
#loop through each variant in the VCF file, form akita predictions for each allele of each variant, and score the mean-squared and max differences 
# between each allele, in each of Akita's 5 predicted cell types. Save results for each variant in two different dataframes (one per scoring method). 
# Each dataframe also contains metadata about the variant position (from the original VCF) and % of N's in the input sequence.
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
    

    
    

max_df = pd.DataFrame(maxed_scores_list)
max_df.to_csv(f"{output_dir}/Akita_max_predictions_{vcf_basename}.csv",index=False,header=True)


msd_df = pd.DataFrame(msd_scores_list)
msd_df.to_csv(f"{output_dir}/Akita_msd_predictions_{vcf_basename}.csv",index=False,header=True)




