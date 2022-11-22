import pandas as pd
import os
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from datetime import date
from datetime import datetime
import argparse
import sys
import shutil


parser = argparse.ArgumentParser(description='Compare Sei predictions for variants againts a null distribution to find extreme variants.', prog='FindSeiSigDNVs.py')
parser.add_argument('--experiment_name', required=True, type=str,  help='Name of experiment to find significant variants. This name should match the corresponding directory in git_root/model_outputs/Sei/')
parser.add_argument('--git_root', required=True, type=str,  help='the path to the top level of the git repository. i.e., path/to/PCGC_Challenge')

args = parser.parse_args()
experiment_name = args.experiment_name 
git_root = args.git_root

print(git_root)


#load sei sequence class scores for the experiment
sei_disease_dir = f"{git_root}/model_outputs/Sei/{experiment_name}" #for the DNVs this should correspond to "{git_root}/model_outputs/Sei/CHD_DNVs_Richter2020"
all_sei_scores = pd.DataFrame()
for folder in os.listdir(sei_disease_dir):
    if folder.endswith('.vcf'): 
        disease_sei_var_df = pd.read_csv(os.path.join(sei_disease_dir,folder,'sequence_class_scores.tsv'),sep = '\t')
        disease_sei_var_df['disease'] = experiment_name
        all_sei_scores = all_sei_scores.append(disease_sei_var_df)




##### load null distribution #####
#null distribution 
null_path = f"{git_root}/model_outputs/Sei/Rare1KGVariants" #this directory will be empty until populated with the Sei variant scores from https://ucsf.box.com/s/gqz9cg3ncj08svu3z5ij3w21xb83nlzw
null_class_scores = pd.DataFrame()
for subdir in os.listdir(null_path):
    if subdir.endswith('.vcf'):
        data_dir = os.path.join(null_path,subdir)
        null_class_scores = null_class_scores.append(pd.read_csv(os.path.join(data_dir,"sequence_class_scores.tsv"),sep='\t'))

score_cols = ['PC1 Polycomb / Heterochromatin', 'L1 Low signal', 'TN1 Transcription',
       'TN2 Transcription', 'L2 Low signal', 'E1 Stem cell', 'E2 Multi-tissue',
       'E3 Brain / Melanocyte', 'L3 Low signal', 'E4 Multi-tissue',
       'TF1 NANOG / FOXA1', 'HET1 Heterochromatin', 'E5 B-cell-like',
       'E6 Weak epithelial', 'TF2 CEBPB', 'PC2 Weak Polycomb',
       'E7 Monocyte / Macrophage', 'E8 Weak multi-tissue', 'L4 Low signal',
       'TF3 FOXA1 / AR / ESR1', 'PC3 Polycomb', 'TN3 Transcription',
       'L5 Low signal', 'HET2 Heterochromatin', 'L6 Low signal', 'P Promoter',
       'E9 Liver / Intestine', 'CTCF CTCF-Cohesin', 'TN4 Transcription',
       'HET3 Heterochromatin', 'E10 Brain', 'TF4 OTX2', 'HET4 Heterochromatin',
       'L7 Low signal', 'PC4 Polycomb / Bivalent stem cell Enh',
       'HET5 Centromere', 'E11 T-cell', 'TF5 AR', 'E12 Erythroblast-like',
       'HET6 Centromere']
# get 99.9th percentile for each sequence class

def get_quantiles(threshold, null_distribution,score_cols):

    
    upper_percentiles = defaultdict()
    lower_percentiles = defaultdict()
    for score_col in score_cols:
        upper_percentiles[score_col] = null_distribution[score_col].quantile(threshold)
        lower_percentiles[score_col] = null_distribution[score_col].quantile(1 - threshold)
    
    return upper_percentiles, lower_percentiles
    
    
upper_percentiles, lower_percentiles = get_quantiles(0.999,null_class_scores,score_cols)


def get_sig_variants(variant_predictions,upper_percentiles,lower_percentiles,score_cols):

    """
    variant_predictions: df with variant effect scores for each alt and ref alleles of variants of interest
    upper_percentiles: Dictionary with upper threshold for each null distribution above which a variant effect is significant
    lower_percentiles: Dictionary with lower threshold for each null distribution below which a variant effect is significant
    score_cols: columns in variant_predictions corresponding to variant effect scores
    
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

sig_variants = get_sig_variants(all_sei_scores,upper_percentiles,lower_percentiles,score_cols)

#save file containing sei sequence class scores for variant<>sequence class combinations that were significant relative to the null distribution. All other entries will be empty
filename = f"{git_root}/EDA/Sei_EDA/SeiSig_{experiment_name}_DNVs.csv"
sig_variants.to_csv(filename)


## for variant<>sequence class combinations that were more extreme than expected under the null distribution, compute the exact p value. This was not done in the first step to save computation time
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

sig_variants_pval = get_pvals(sig_variants,null_class_scores,score_cols)

#save file containing p values for variant<> sequence class combinations that were significant. All other entries will be empty
filename = f"{git_root}/EDA/Sei_EDA/SeiSigPVals_{experiment_name}_DNVs.csv"
sig_variants_pval.to_csv(filename)