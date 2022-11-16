"""
This script finds variants that were identified via OpenTargets Genetics to be a lead or tag variant associated with some phenotype,
but did not align with the reference genome
Example:
The variant 6_26745970 C>T

If i follow the rsid for this variant (rs13213200) supplied by OpenTargets on gnomAD I see the following:

gnomAD 2.1: 6-26755915-C-T. Here, the chrom, ref and alt alleles are the same as in OpenTargets Genetics, but the position is different (26755915 vs 26745970).
gnomAD 3.1: 6-26745970-G-A. Here, the chrom and position are the same as in Open Targets Genetics, but the ref and alt alleles are different (G>A vs C>T)

Further, if I search for position chr6:26745970 on the UCSC genome browser using both GRCh37 and GRCh38 assemblies, the reference allele is G in both cases, never C.

This suggests 6_26745970 C>T is not a real variant. 

This script identified other such variants so I can remove them from analysis
"""

import os
import pandas as pd
from Bio.Seq import Seq

sei_dir="/wynton/home/hernandez/shirondru/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/Sei"

#### get all sequence class scores from all diseases into one dataframe to find where there was no ref_match
main_seq_class_df = pd.DataFrame()
for disease in os.listdir(sei_dir):
	if disease != 'test':
		disease_path = os.path.join(sei_dir,disease)
		for _ in os.listdir(disease_path):
			if _.endswith('.vcf'):
				sei_pred_folder = os.path.join(disease_path,_)
				main_seq_class_df = main_seq_class_df.append(pd.read_csv(os.path.join(sei_pred_folder,'sequence_class_scores.tsv'),sep='\t')) #append sequence class scores from this split_vcf into main_seq_class_df


no_match_df = main_seq_class_df[main_seq_class_df['ref_match']!=True]
output_dir = "/wynton/home/hernandez/shirondru/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/PsychENCODE_GWAS_scripts/test/data_for_tests/bad_variants"
no_match_df.to_csv(os.path.join(output_dir,"no_ref_match_variants.csv"))


#save the mismatched variants as a vcf as well, so I can form prediction with the reverse complement, as it seems some of the variants had issues because the sequence was flipped during liftover from hg19 to hg38
