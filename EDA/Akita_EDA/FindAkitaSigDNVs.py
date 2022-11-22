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



parser = argparse.ArgumentParser(description='Compare Akita predictions for variants againts a null distribution to find extreme variants.', prog='FindAkitaDNVs.py')
parser.add_argument('--scoring_system', required=True, type=str,  help='Either `msd` or `max`; one of two scoring systems used to calculate variant effect scores along a sequence')
parser.add_argument('--experiment_name', required=True, type=str,  help='Name of experiment to find significant variants. This name should match the corresponding directory in git_root/model_outputs/Akita/')
parser.add_argument('--git_root', required=True, type=str,  help='the path to the top level of the git repository. i.e., path/to/PCGC_Challenge')

args = parser.parse_args()
scoring_system = str(args.scoring_system) #either "msd" or "max"
experiment_name = args.experiment_name 
git_root = args.git_root




########## Load null distribution ############
null_path = f"{git_root}/model_outputs/Akita/Rare1KGVariants"
null_dist = pd.DataFrame()
for idx, file in enumerate(os.listdir(null_path)):
  if scoring_system in file and file.endswith('.vcf.csv'):
    split_null_df = pd.read_csv(os.path.join(null_path,file))
    null_dist = null_dist.append(split_null_df)

score_cols = ['0_HFF', '1_H1hESC', '2_GM12878', '3_IMR90','4_HCT116']



######## Load Akita variant scores for this experiment #######
akita_dir = f"{git_root}/model_outputs/Akita/CHD_DNVs_Richter2020"
disease_akita_scores = pd.DataFrame()
disease_akita_path = os.path.join(akita_dir,experiment_name)
for file in os.listdir(disease_akita_path):
    if file.endswith('.csv') and f'{scoring_system}' in file: #only open the csv files with the variant score (either max or summed, depending on the input argument to this file)
        disease_akita_var_df = pd.read_csv(os.path.join(akita_dir,experiment_name,file))
                                                                   
        disease_akita_var_df['disease'] = experiment_name
        disease_akita_scores = disease_akita_scores.append(disease_akita_var_df)
disease_akita_scores = disease_akita_scores[disease_akita_scores['chrom']!= 'chrom'] #some rows have the header inside. remove them
disease_akita_scores[score_cols] = disease_akita_scores[score_cols].astype(float) #coerce the akita scores to be float, since they were strings before due to header being inside the rows
disease_akita_scores[score_cols] = np.log(disease_akita_scores[score_cols]) #upper percentiles will be returned ln transformed, so apply that to the  Akita scores, else everything will appear significant
def get_quantiles(threshold, null_distribution,score_cols):

    
    upper_percentiles = defaultdict() #only looking for big differences between ref and alt. No need for lower percentiles; not a matter of upregulation or downregulation, but is the contact highly different or not?
    for score_col in score_cols:
        null_distribution[score_col] = np.log(null_distribution[score_col])
        upper_percentiles[score_col] = null_distribution[score_col].quantile(threshold)
    return upper_percentiles 
    
upper_percentiles = get_quantiles(0.999,null_dist,score_cols)



def get_sig_variants(variant_predictions,upper_percentiles,score_cols):
    """
    variant_predictions: df with variant effect scores for each ref and alt allele in the variant of interest
    upper_percentiles: Dictionary with upper threshold for each null distribution above which a variant effect is significant
    score_cols: columns in variant_predictions corresponding to variant effect scores
    
    """

    sig_variants = pd.DataFrame()
    variant_info = ['chrom','pos','ref','alt','disease']
    for score_col in score_cols:
        #append variants to df if they are more extreme than the upper or lower threshold
        sig_variants = sig_variants.append(variant_predictions[
                                    (variant_predictions[score_col] >= upper_percentiles[score_col])][variant_info + [score_col]].copy())

    #group variants that are extreme against different null distributions back to the same row
    sig_variants = sig_variants.groupby(['chrom','pos','ref','alt']).max().reset_index()
    
    return sig_variants

sig_variants = get_sig_variants(disease_akita_scores,upper_percentiles,score_cols)

## for variant<>cell_type combinations that were more extreme than expected under the null distribution, compute the exact p value. This was not done in the first step to save computation time
filename = f"{git_root}/EDA/Akita_EDA/AkitaSig_{scoring_system}_{experiment_name}_DNVs.csv"
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


#save file containing p values for variant<>cell typecombinations that were significant. All other entries will be empty
filename = f"{git_root}/EDA/Enformer_EDA/AkitaSigPVals_{scoring_system}_{experiment_name}_DNVs.csv"
sig_variants_pval.to_csv(filename)
