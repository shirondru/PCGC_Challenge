import pandas as pd
import os
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
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

parser = argparse.ArgumentParser(description='Compare Enformer predictions for variants againts a null distribution to find extreme variants.', prog='FindEnformerSigDNVs.py')
parser.add_argument('--scoring_system', required=True, type=str,  help='Either `summed` or `max`; one of two scoring systems used to calculate variant effect scores along a sequence')
parser.add_argument('--relevant_cols', required=True, type=str,  help='Path to text file containing relevant Enformer features for use with this experiment')
parser.add_argument('--experiment_name', required=True, type=str,  help='Name of experiment to find significant variants. This name should match the corresponding directory in git_root/model_outputs/Enformer/')
parser.add_argument('--git_root', required=True, type=str,  help='the path to the top level of the git repository. i.e., path/to/PCGC_Challenge')

args = parser.parse_args()
scoring_system = str(args.scoring_system) #either "summed" or "max"
experiment_name = args.experiment_name
git_root=args.git_root

#get the list of enformer tracks that are chosen to be used/relevant for this experiment. For the Ricther 2020 DNVs, this includes tracks from tissues in the heart, blood and blood vessels, and smooth muscle
relevant_cols_file = str(args.relevant_cols)
relevant_cols_path=f"{git_root}/EDA/Enformer_EDA/"
with open(os.path.join(relevant_cols_path,relevant_cols_file),'r') as f:
  file = f.read()
  relevant_cols = file.split('\n')
relevant_cols = [x for x in relevant_cols if x!= ''] #remove emptry strings
relevant_cols = [x.lower() for x in relevant_cols] #coerce all strings to be lowercase



# Download targets from Basenji2 dataset 
# Cite: Kelley et al Cross-species regulatory sequence activity prediction. PLoS Comput. Biol. 16, e1008050 (2020).
df_targets = pd.read_csv(f"{git_root}/EDA/Enformer_EDA/enformer_df_targets.csv", sep='\t')




#The column names in the Enformer predictions were truncated. Harmonize this
enformer_relevant_column_idxs = df_targets[df_targets['description'].str.lower().isin(relevant_cols)]['index'].tolist()
null_path = f"{git_root}/model_outputs/Enformer/Rare1KGVariants" #this directory will be empty until populated with the Enformer variant scores from https://ucsf.box.com/s/gqz9cg3ncj08svu3z5ij3w21xb83nlzw
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


def get_sig_variants(variant_predictions,upper_percentiles,lower_percentiles,score_cols):
    """
    variant_predictions: df with variant effect scores for each ref and alt allele of each variant
    upper_percentiles: Dictionary with upper threshold for each null distribution above which a variant effect is significant
    lower_percentiles: Dictionary with lower threshold for each null distribution below which a variant effect is significant
    score_cols: columns in variant_preditions corresponding to variant effect scores
    """
    sig_variants = pd.DataFrame()
    variant_info = ['chrom','pos','ref','alt','disease']
    for score_col in score_cols:
        #append variants to df if they are more extreme than the upper or lower threshold
        sig_variants = sig_variants.append(variant_predictions[
                                    (variant_predictions[score_col] >= upper_percentiles[score_col]) | 
                                    (variant_predictions[score_col] <= lower_percentiles[score_col])][variant_info + [score_col]].copy())
    #group variants that are extreme against different null distributions back to the same row
    sig_variants = sig_variants.groupby(['chrom','pos','ref','alt']).max().reset_index()
    return sig_variants




#load Enformer variant scores for this experiment. These are separated into many files, so they are appeneded together
enformer_dir = f"{git_root}/model_outputs/Enformer/"
disease_enformer_scores = pd.DataFrame()
disease_enformer_path = os.path.join(enformer_dir,experiment_name)
for file in os.listdir(disease_enformer_path):
    if file.endswith('.csv') and f'{scoring_system}' in file: #only open the csv files with the variant score (either max or summed, depending on the input argument to this file)
        disease_enformer_var_df = pd.read_csv(os.path.join(enformer_dir,experiment_name,file),
                                                                   usecols = ['chrom','pos','ref','alt'] + enformer_relevant_cols)
        disease_enformer_var_df['disease'] = experiment_name
        disease_enformer_scores = disease_enformer_scores.append(disease_enformer_var_df)
sig_variants = get_sig_variants(disease_enformer_scores,upper_percentiles,lower_percentiles,enformer_relevant_cols)                                   

#save file containing enformer scores for variant<>genomic track combinations that were significant relative to the null distribution. All other entries will be empty
filename = f"{git_root}/EDA/Enformer_EDA/EnformerSig_{scoring_system}_{experiment_name}_DNVs.csv"
sig_variants.to_csv(filename)


## for variant<>genomic track combinations that were more extreme than expected under the null distribution, compute the exact p value. This was not done in the first step to save computation time
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

#save file containing p values for variant<>genomic track combinations that were significant. All other entries will be empty
filename = f"{git_root}/EDA/Enformer_EDA/EnformerSigPVals_{scoring_system}_{experiment_name}_DNVs.csv"
sig_variants_pval.to_csv(filename)

