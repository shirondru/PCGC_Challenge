import os
import pandas as pd

test_dir="/wynton/home/hernandez/shirondru/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/test/testEnformerPredictions/"
ground_truth_df = pd.read_csv(test_dir+'ground_truth_summed_Enformer_preds.csv')

ground_truth_df = ground_truth_df.rename(columns = {'chrom':'CHROM',
	                                                 'pos':'POS',
	                                                 'ref':'REF',
	                                                 'id':'ID',
	                                                  'alt':'ALT'})
#fill empty columns
ground_truth_df['QUAL'] = ""
ground_truth_df['INFO'] = ""
ground_truth_df['FILTER'] = ""
#force column order
ground_truth_df = ground_truth_df[['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']] 

#save VCF with header for Enformer test
header = """##fileformat=VCFv4.1
##fileDate=20220722
##source=/wynton/home/hernandez/shirondru/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/test/testEnformerPredictions/gen_testEnformer_vcfs_1.py
##reference=/wynton/home/hernandez/shirondru/pollard_lab/data/hg38_genome.fa
#CHROM POS ID REF ALT QUAL FILTER INFO
"""
#write VCF header
#save chunks of df to later parallelize model predictions on each chunk


output_VCF = os.path.join(test_dir,f"ground_truth_Enformer_SNPs.vcf")
with open(output_VCF, 'w') as vcf:
    vcf.write(header)

ground_truth_df.to_csv(output_VCF, sep="\t",mode = 'a',index=False) 	                                                  