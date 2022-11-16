import os
import pandas as pd
import numpy as np
from collections import defaultdict



############## THESE TESTS CHECK THAT PREDICTIONS FOR ALL VARIANTS IN EACH MODEL WAS FORMED ###################

def test_all_Enformer_max_predictions_finished_properly():
	'''
	This test checks that every vcf file that was created to hold PsychENCODE disease variants generated a corresponding csv file with the model's predictions representing the max differences along the sequence axis between the reference and alternate allele for all SNPs
	The reason multiple vcf's were generated was to parallelize the process -- with multiple vcf's each containing max 5000 variants for each disease.
	This allowed for multiple jobs to be submitted, speeding up the process of forming predictions.
	'''
	enformer_dir = "/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/Enformer"
	vcf_list_dir = "/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/getGWASVariants/GWAS_psychENCODE_LeadTagVariants" #Subdirectories in this directory contain the split vcfs (len <=5000 each)
	for disease in os.listdir(enformer_dir):
		if disease != 'test':#predictions for unit tests are in here too, ignore these because there is a test to make sure they had the expected results. Here we are testing all the disease predictions finished.
			for vcf in os.listdir(os.path.join(vcf_list_dir,disease,'AkitaEnformer')): ##this is where all the split vcfs (len <=5000 each) were stored. I should get a csv with model predictions for each of these
				if vcf.endswith(".vcf"):
					assert os.path.exists(os.path.join(enformer_dir,disease,f"Enformer_max_predictions_PsychENCODEGWASLeadTagVariants_{vcf}.csv")),"You don't have Enformer Max predictions for every vcf file!"

def test_all_Enformer_summed_predictions_finished_properly():
	'''
	This test checks that every vcf file that was created to hold PsychENCODE disease variants generated a corresponding csv file with the model's predictions representing the summed differences along the sequence axis between the reference and alternate allele for all SNPs
	The reason multiple vcf's were generated was to parallelize the process -- with multiple vcf's each containing max 5000 variants for each disease.
	This allowed for multiple jobs to be submitted, speeding up the process of forming predictions.
	'''
	enformer_dir = "/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/Enformer"
	vcf_list_dir = "/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/getGWASVariants/GWAS_psychENCODE_LeadTagVariants" #Subdirectories in this directory contain the split vcfs (len <=5000 each)
	for disease in os.listdir(enformer_dir):
		if disease != 'test':#predictions for unit tests are in here too, ignore these because there is a test to make sure they had the expected results. Here we are testing all the disease predictions finished.
			for vcf in os.listdir(os.path.join(vcf_list_dir,disease,'AkitaEnformer')): ##this is where all the split vcfs (len <=5000 each) were stored. I should get a csv with model predictions for each of these
				if vcf.endswith(".vcf"):
					assert os.path.exists(os.path.join(enformer_dir,disease,f"Enformer_summed_predictions_PsychENCODEGWASLeadTagVariants_{vcf}.csv")),"You don't have Enformer summed predictions for every vcf file!"





def test_all_max_Akita_predictions_finished_properly():
	'''
	This test checks that every vcf file that was created to hold PsychENCODE disease variants generated a corresponding csv file with the model's predictions representing the max differences between the reference and alternate allele Hi-C contact map for all SNPs
	The reason multiple vcf's were generated was to parallelize the process -- with multiple vcf's each containing max 5000 variants for each disease.
	This allowed for multiple jobs to be submitted, speeding up the process of forming predictions.
	'''
	akita_dir = "/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/Akita"
	vcf_list_dir = "/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/getGWASVariants/GWAS_psychENCODE_LeadTagVariants" #Subdirectories in this directory contain the split vcfs (len <=5000 each)
	for disease in os.listdir(akita_dir):
		if disease != 'test': #predictions for unit tests are in here too, ignore these because there is a test to make sure they had the expected results. Here we are testing all the disease predictions finished.
			for vcf in os.listdir(os.path.join(vcf_list_dir,disease,'AkitaEnformer')): ##this is where all the split vcfs (len <=5000 each) were stored. I should get a csv with model predictions for each of these
				if vcf.endswith(".vcf"):
					assert os.path.exists(os.path.join(akita_dir,disease,f"Akita_max_predictions_PsychENCODEGWASLeadTagVariants_{vcf}.csv")),"You don't have Akita Max predictions for every vcf file!"



def test_all_msd_Akita_predictions_finished_properly():
	'''
	This test checks that every vcf file that was created to hold PsychENCODE disease variants generated a corresponding csv file with the model's predictions representing the mean squared differences between the reference and alternate allele Hi-C contact map for all SNPs
	The reason multiple vcf's were generated was to parallelize the process -- with multiple vcf's each containing max 5000 variants for each disease.
	This allowed for multiple jobs to be submitted, speeding up the process of forming predictions.
	'''
	akita_dir = "/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/Akita"
	vcf_list_dir = "/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/getGWASVariants/GWAS_psychENCODE_LeadTagVariants" #Subdirectories in this directory contain the split vcfs (len <=5000 each)
	for disease in os.listdir(akita_dir):
		if disease != 'test': #predictions for unit tests are in here too, ignore these because there is a test to make sure they had the expected results. Here we are testing all the disease predictions finished.
			for vcf in os.listdir(os.path.join(vcf_list_dir,disease,'AkitaEnformer')): ##this is where all the split vcfs (len <=5000 each) were stored. I should get a csv with model predictions for each of these
				if vcf.endswith(".vcf"):
					assert os.path.exists(os.path.join(akita_dir,disease,f"Akita_msd_predictions_PsychENCODEGWASLeadTagVariants_{vcf}.csv")),"You don't have Akita mean squared difference predictions for every vcf file!"



def test_all_Sei_predictions_finished_properly():
	'''
	This test checks that every vcf file that was created to hold PsychENCODE disease variants generated a corresponding sequence_class_scores.tsv file
	The reason multiple vcf's were generated was to parallelize the process -- with multiple vcf's each containing max 5000 variants for each disease.
	This allowed for multiple jobs to be submitted, speeding up the process of forming predictions.
	'''
	sei_dir = "/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/Sei"
	vcf_list_dir = "/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/getGWASVariants/GWAS_psychENCODE_LeadTagVariants" #this is where all the split vcfs (len <=5000 each) were stored. I should get a sequence_class_scores.tsv for each of these
	for disease in os.listdir(sei_dir):
		for vcf_or_folder in os.listdir(os.path.join(vcf_list_dir,disease)):
			if vcf_or_folder.endswith(".vcf"):
				vcf = vcf_or_folder
				assert os.path.exists(os.path.join(sei_dir,disease,vcf,"sequence_class_scores.tsv")),"You don't have Sei predictions for every vcf file!"


def test_expected_SNP_predictions_per_disease_Sei():
	'''
	Take the original disease vcf containing all the SNPs for each disease, the `sequence_class_scores.tsv` output file containing predictions from Sei for each disease, and either doing an inner join on chrom, pos, ref, alt columns
	see the number of rows doesn't decrease. If the number of rows doesn't change, that means every variant in the original disease vcf had predictions formed by Sei for this disease, and there were no SNPs in the output that were not in the original disease vcf

	This will also test that the original disease vcf was split properly into smaller chunks -- because the model's output for each chunk covered all variants
	And this will also test that the expected number of variants had predictions formed for each disease, and none were missing or added.
	'''	
	

	path_to_vcfs = "/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/getGWASVariants/GWAS_psychENCODE_LeadTagVariants"
	sei_dir = "/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/Sei"
	colnames = ['chrom','pos','id','ref','alt','qual','filter','info']
	for disease in os.listdir(sei_dir):
		if disease != 'test':
			original_vcf_df = pd.read_csv(os.path.join(path_to_vcfs,f"PsychENCODE_GWASVariants_{disease}_SeiNoHeader.vcf"),sep='\t',names=colnames) #number of rows in the main disease vcf that was split into smaller vcfs to parallelize jobs
			original_vcf_df['variant'] = original_vcf_df['chrom'].astype(str) + ':' + original_vcf_df['pos'].astype(str) + ' ' + original_vcf_df['ref'].astype(str) + '>' + original_vcf_df['alt'].astype(str)
			original_num_variants = original_vcf_df.shape[0]
			#loop through the Sei predictions for this disease
			disease_pred_df = pd.DataFrame()
			for folder in os.listdir(os.path.join(sei_dir,disease)):
				if folder.endswith('.vcf'): 
					disease_pred_df= disease_pred_df.append(pd.read_csv(os.path.join(sei_dir,disease,folder,'sequence_class_scores.tsv'),sep='\t',usecols = ['chrom','pos','ref','alt'])) #add number of rows in the split_vcf to tally
			disease_pred_df['variant'] = disease_pred_df['chrom'].astype(str) + ':' + disease_pred_df['pos'].astype(str) + ' ' + disease_pred_df['ref'].astype(str) + '>' + disease_pred_df['alt'].astype(str)
			merge_result = original_vcf_df.merge(disease_pred_df,on='variant')

			assert original_num_variants == merge_result.shape[0],f"The SNPs in the original {disease} VCF are not exactly the same as those that came out of the Sei predictions!"

def test_expected_SNP_predictions_per_disease_Enformer_max():
	'''
	Take the original disease vcf containing all the SNPs for each disease, the output file containing predictions from Enformer (max) for each disease, and either doing an inner join on chrom, pos, ref, alt columns
	see the number of rows doesn't decrease. If the number of rows doesn't change, that means every variant in the original disease vcf had predictions formed by Enformer for this disease, and there were no SNPs in the output that were not in the original disease vcf

	This will also test that the original disease vcf was split properly into smaller chunks -- because the model's output for each chunk covered all variants
	And this will also test that the expected number of variants had predictions formed for each disease, and none were missing or added.
	'''	
	path_to_vcfs = "/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/getGWASVariants/GWAS_psychENCODE_LeadTagVariants"
	enformer_dir = "/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/Enformer"
	colnames = ['chrom','pos','id','ref','alt','qual','filter','info']
	for disease in os.listdir(enformer_dir):
		if disease != 'test':
			original_vcf_df = pd.read_csv(os.path.join(path_to_vcfs,f"PsychENCODE_GWASVariants_{disease}.vcf"),sep='\t',names=colnames,skiprows=6) #number of rows in the main disease vcf that was split into smaller vcfs to parallelize jobs
			original_vcf_df['variant'] = original_vcf_df['chrom'].astype(str) + ':' + original_vcf_df['pos'].astype(str) + ' ' + original_vcf_df['ref'].astype(str) + '>' + original_vcf_df['alt'].astype(str)
			original_num_variants = original_vcf_df.shape[0]
			
			#loop through the Enformer predictions for this disease
			disease_pred_df = pd.DataFrame()
			for file in os.listdir(os.path.join(enformer_dir,disease)):
				if file.endswith('.csv') and 'max' in file: 
					disease_pred_df= disease_pred_df.append(pd.read_csv(os.path.join(enformer_dir,disease,file),usecols = ['chrom','pos','ref','alt'])) #add number of rows in the split_vcf to tally
			disease_pred_df['variant'] = disease_pred_df['chrom'].astype(str) + ':' + disease_pred_df['pos'].astype(str) + ' ' + disease_pred_df['ref'].astype(str) + '>' + disease_pred_df['alt'].astype(str)
			merge_result = original_vcf_df.merge(disease_pred_df,on='variant')

			assert original_num_variants == merge_result.shape[0],f"The SNPs in the original {disease} VCF are not exactly the same as those that came out of the Enformer max predictions!"

def test_expected_SNP_predictions_per_disease_Enformer_summed():
	'''
	Take the original disease vcf containing all the SNPs for each disease, the output file containing predictions from Enformer (summed) for each disease, and either doing an inner join on chrom, pos, ref, alt columns
	see the number of rows doesn't decrease. If the number of rows doesn't change, that means every variant in the original disease vcf had predictions formed by Enformer for this disease, and there were no SNPs in the output that were not in the original disease vcf

	This will also test that the original disease vcf was split properly into smaller chunks -- because the model's output for each chunk covered all variants
	And this will also test that the expected number of variants had predictions formed for each disease, and none were missing or added.
	'''	
	

	path_to_vcfs = "/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/getGWASVariants/GWAS_psychENCODE_LeadTagVariants"
	enformer_dir = "/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/Enformer"
	colnames = ['chrom','pos','id','ref','alt','qual','filter','info']
	for disease in os.listdir(enformer_dir):
		if disease != 'test':
			original_vcf_df = pd.read_csv(os.path.join(path_to_vcfs,f"PsychENCODE_GWASVariants_{disease}.vcf"),sep='\t',names=colnames,skiprows=6) #number of rows in the main disease vcf that was split into smaller vcfs to parallelize jobs
			original_vcf_df['variant'] = original_vcf_df['chrom'].astype(str) + ':' + original_vcf_df['pos'].astype(str) + ' ' + original_vcf_df['ref'].astype(str) + '>' + original_vcf_df['alt'].astype(str)
			original_num_variants = original_vcf_df.shape[0]

			#loop through the Enformer predictions for this disease
			disease_pred_df = pd.DataFrame()
			for file in os.listdir(os.path.join(enformer_dir,disease)):
				if file.endswith('.csv') and 'summed' in file: 
					disease_pred_df= disease_pred_df.append(pd.read_csv(os.path.join(enformer_dir,disease,file),usecols = ['chrom','pos','ref','alt'])) #add number of rows in the split_vcf to tally
			disease_pred_df['variant'] = disease_pred_df['chrom'].astype(str) + ':' + disease_pred_df['pos'].astype(str) + ' ' + disease_pred_df['ref'].astype(str) + '>' + disease_pred_df['alt'].astype(str)
			merge_result = original_vcf_df.merge(disease_pred_df,on='variant')

			assert original_num_variants == merge_result.shape[0],f"The SNPs in the original {disease} VCF are not exactly the same as those that came out of the Enformer summed predictions!"


def test_expected_SNP_predictions_per_disease_Akita_max():
	'''
	Take the original disease vcf containing all the SNPs for each disease, the output file containing predictions from Akita (max) for each disease, and either doing an inner join on chrom, pos, ref, alt columns
	see the number of rows doesn't decrease. If the number of rows doesn't change, that means every variant in the original disease vcf had predictions formed by Akita for this disease, and there were no SNPs in the output that were not in the original disease vcf

	This will also test that the original disease vcf was split properly into smaller chunks -- because the model's output for each chunk covered all variants
	And this will also test that the expected number of variants had predictions formed for each disease, and none were missing or added.
	'''	
	path_to_vcfs = "/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/getGWASVariants/GWAS_psychENCODE_LeadTagVariants"
	akita_dir = "/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/Akita"
	colnames = ['chrom','pos','id','ref','alt','qual','filter','info']
	for disease in os.listdir(akita_dir):
		if disease != 'test':
			original_vcf_df = pd.read_csv(os.path.join(path_to_vcfs,f"PsychENCODE_GWASVariants_{disease}.vcf"),sep='\t',names=colnames,skiprows=6) #number of rows in the main disease vcf that was split into smaller vcfs to parallelize jobs
			original_vcf_df['variant'] = original_vcf_df['chrom'].astype(str) + ':' + original_vcf_df['pos'].astype(str) + ' ' + original_vcf_df['ref'].astype(str) + '>' + original_vcf_df['alt'].astype(str)
			original_num_variants = original_vcf_df.shape[0]
			
			#loop through the Enformer predictions for this disease
			disease_pred_df = pd.DataFrame()
			for file in os.listdir(os.path.join(akita_dir,disease)):
				if file.endswith('.csv') and 'max' in file: 
					disease_pred_df= disease_pred_df.append(pd.read_csv(os.path.join(akita_dir,disease,file),usecols = ['chrom','pos','ref','alt'])) #add number of rows in the split_vcf to tally
			disease_pred_df['variant'] = disease_pred_df['chrom'].astype(str) + ':' + disease_pred_df['pos'].astype(str) + ' ' + disease_pred_df['ref'].astype(str) + '>' + disease_pred_df['alt'].astype(str)
			merge_result = original_vcf_df.merge(disease_pred_df,on='variant')

			assert original_num_variants == merge_result.shape[0],f"The SNPs in the original {disease} VCF are not exactly the same as those that came out of the Akita max predictions!"

def test_expected_SNP_predictions_per_disease_Akita_msd():
	'''
	Take the original disease vcf containing all the SNPs for each disease, the output file containing predictions from Akita (mean squared difference) for each disease, and either doing an inner join on chrom, pos, ref, alt columns
	see the number of rows doesn't decrease. If the number of rows doesn't change, that means every variant in the original disease vcf had predictions formed by Akita for this disease, and there were no SNPs in the output that were not in the original disease vcf

	This will also test that the original disease vcf was split properly into smaller chunks -- because the model's output for each chunk covered all variants
	And this will also test that the expected number of variants had predictions formed for each disease, and none were missing or added.
	'''	
	path_to_vcfs = "/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/getGWASVariants/GWAS_psychENCODE_LeadTagVariants"
	akita_dir = "/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/Akita"
	colnames = ['chrom','pos','id','ref','alt','qual','filter','info']
	for disease in os.listdir(akita_dir):
		if disease != 'test':
			original_vcf_df = pd.read_csv(os.path.join(path_to_vcfs,f"PsychENCODE_GWASVariants_{disease}.vcf"),sep='\t',names=colnames,skiprows=6) #number of rows in the main disease vcf that was split into smaller vcfs to parallelize jobs
			original_vcf_df['variant'] = original_vcf_df['chrom'].astype(str) + ':' + original_vcf_df['pos'].astype(str) + ' ' + original_vcf_df['ref'].astype(str) + '>' + original_vcf_df['alt'].astype(str)
			original_num_variants = original_vcf_df.shape[0]
			
			#loop through the Enformer predictions for this disease
			disease_pred_df = pd.DataFrame()
			for file in os.listdir(os.path.join(akita_dir,disease)):
				if file.endswith('.csv') and 'msd' in file: 
					disease_pred_df= disease_pred_df.append(pd.read_csv(os.path.join(akita_dir,disease,file),usecols = ['chrom','pos','ref','alt'])) #add number of rows in the split_vcf to tally
			disease_pred_df['variant'] = disease_pred_df['chrom'].astype(str) + ':' + disease_pred_df['pos'].astype(str) + ' ' + disease_pred_df['ref'].astype(str) + '>' + disease_pred_df['alt'].astype(str)
			merge_result = original_vcf_df.merge(disease_pred_df,on='variant')

			assert original_num_variants == merge_result.shape[0],f"The SNPs in the original {disease} VCF are not exactly the same as those that came out of the Akita msd predictions!"


############### THESE TESTS CHECK THAT DEEP LEARNING MODELS FORM EXPECTED PREDICTIONS AND ARE BEING USED PROPERLY #############

def test_Sei_predictions_PsychENCODE_variants():
	"""
	Many variants from the GWAS catalog for diseases relevant to the PsychENCODE project had predictions formed with Sei using my code. 
	I took 15 random variants from this pool and also formed predictions using the website hosing the Sei framework, called those predictions for those 15 variants the ground truth, 
	and I test here that the ground truth is approximately equal to the predictions formed using my code
	"""

	def get_variant(chrom,pos,ref,alt):
		"""
		Helper function to query 3 df's -- with only one of which containing the variant of interest -- and returning the row corresponding to that variant
		"""
		def query_df_for_variant(disease_sei_pred_df,chrom,pos,ref,alt):
			"""
			Extra helper function to generally query any df for a row corresponding to a variant
			"""
			return disease_sei_pred_df[(disease_sei_pred_df['chrom']==chrom) & (disease_sei_pred_df['pos'] == pos) & (disease_sei_pred_df['ref'] == ref) & (disease_sei_pred_df['alt'] == alt)]

		bipolar_variant = query_df_for_variant(bipolar_disorder_preds,chrom,pos,ref,alt)
		adhd_variant = query_df_for_variant(adhd_disorder_preds,chrom,pos,ref,alt)
		alz_variant = query_df_for_variant(alzheimers_disorder_preds,chrom,pos,ref,alt)

		#return whichever row holds the variant that was queried
		if bipolar_variant.shape[0] == 1:
			return bipolar_variant
		elif adhd_variant.shape[0] == 1:
			return adhd_variant
		elif alz_variant.shape[0] == 1:
			return alz_variant

	sei_dir = "/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/Sei"
	sei_framework_preds = pd.read_csv("/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/PsychENCODE_GWAS_scripts/test/data_for_tests/sei_website_preds_for_random_PsychENCODE_variants.tsv",sep='\t') #Holds the ground truth predictions from the Sei website for the 15 variants 
	## my predictions for the 15 variants are in these files (plus a bunch of other variants)
	bipolar_disorder_preds = pd.read_csv(os.path.join(sei_dir,"BipolarDisorder","PsychENCODE_GWASVariants_BipolarDisorder_SeiNoHeader.vcf00.vcf","sequence_class_scores.tsv"),sep='\t') #5 of the variants were the first 5 in the BipolarDisorder vcf, so I open the vcf00 file for this disease
	adhd_disorder_preds = pd.read_csv(os.path.join(sei_dir,"ADHD","PsychENCODE_GWASVariants_ADHD_SeiNoHeader.vcf14.vcf","sequence_class_scores.tsv"),sep='\t')  #5 of the variants were the last 5 for ADHD, so I opened the vcf14 file for this disease
	alzheimers_disorder_preds = pd.read_csv(os.path.join(sei_dir,"Alzheimer","PsychENCODE_GWASVariants_Alzheimer_SeiNoHeader.vcf00.vcf","sequence_class_scores.tsv"),sep='\t') # 5 of the variants were the first 5 for Alzheimer, so I opened the vcf00 file for this disease


	#find my predictions for the same variants
	#They are the first 5 predictions for Bipolar Disorder, the last 5 variants for ADHD, and the first 5 variants for Alzheimer 
	#because Sei predictions come out sorted by the largest seqclass_max_abs_diff, I need to search each file for the correct variant. 
	my_preds = pd.DataFrame() #instantiate df to hold predictions formed from Sei with my code
	for idx, row in sei_framework_preds.iterrows():
		chrom = row['chrom']
		pos = row['pos']
		ref = row['ref']
		alt = row['alt']

		#Because the Sei predictions in sei_framework_preds are sorted by seqclass_max_abs_diff,
		#check each df to see if the row holding my predictions for the variant is within using `get_variant`
		correct_variant_row = get_variant(chrom,pos,ref,alt)
		my_preds = my_preds.append(correct_variant_row)
	feature_cols = my_preds.columns[-40:]


	assert np.allclose(my_preds[feature_cols],sei_framework_preds[feature_cols],rtol=1e-3, atol=1e-3), "My Sei predictions for at least 1/15 PsychENCODE GWAS variants is different than obtained using the Sei website!"


def test_Akita_predictions():
	"""
	Test I am properly implementing Akita to form predictions for variants by checking to see I get the same values for variants as Katie Gjoni from the Pollard lab.
	Here I test that we both get the same mean squared difference values for 3 SNPs
	"""
	my_preds = pd.read_csv("/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/Akita/test/Akita_msd_predictions_PsychENCODEGWASLeadTagVariants_Variants_KatieG_formed_AkitaPreds.vcf.csv")
	katie_preds = pd.read_csv("/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/PsychENCODE_GWAS_scripts/test/data_for_tests/AkitaPreds_by_KatieG/AkitaPreds_by_KatieG.tsv",sep='\t')

	#Katie G and I handle indels differently, so it is not expected we get the same mean square difference values for the indels. 
	#therefore, remove all variants that aren't SNPs (2 indels), and only test the SNPs are the same. 3 SNPs will be left
	my_preds = my_preds[(my_preds['ref'].str.len() == 1) & (my_preds['alt'].str.len() == 1)]
	katie_preds = katie_preds[(katie_preds['Ref'].str.len() == 1) & (katie_preds['Alt'].str.len() == 1)]

	assert np.allclose(my_preds['0_HFF'].astype(float),katie_preds['MSE']),'Your Akita predictions for SNPs are different than a different member of the Pollard Lab!'

def test_Enformer_predictions():
	"""

	This will test that Enformer's predictions using my code for 5 variants taken from a ClinVar vcf file, supplied by the authors, is extremely close (each prediction within approx 2 decimal places) to the predictions obtained running the authors enformer-usage.ipynb notebook on Google Colab. 
	Those predictions are considered to be ground truth. The expected/ground truth predictions that enformer-usage.ipynb was used to create are stored in ground_truth_summed_Enformer_preds.csv and a copy of the notebok used to generate them is in  enformer_usage_get_ground_truth_summed_Enformer_preds.ipynb

	`deepmind_output` generated in:
	/Users/shirondrusinsky/Documents/Pollard_Lab/PsychENCODE_GWAS_scripts/test/data_for_tests/testEnformerPredictions/enformer_usage_get_ground_truth_summed_Enformer_preds.ipynb
	"""
	test_output = pd.read_csv("/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/Enformer/test/Enformer_msd_predictions_PsychENCODEGWASLeadTagVariants_ground_truth_Enformer_SNPs.vcf.csv")
	deepmind_output = pd.read_csv("/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/PsychENCODE_GWAS_scripts/test/data_for_tests/test_Enformer_predictions_data/ground_truth_summed_Enformer_preds.csv")

	assert np.allclose(test_output[test_output.columns[5:]].to_numpy(),deepmind_output[deepmind_output.columns[5:]].to_numpy(),rtol=1e-2,atol=1e-2), f"Your Enformer predictions don't match the authors prediction somewhere in the first 5 ClinVar vcf file variants!!"


####################### Test that code used to find variants with many more significant tracks than other variants in the same disease-associated locus worked as expected ###########

def test_consistent_std_val_per_identifier_col_Enformer():
	"""
	This tests that all std_val columns have only 1 unique value per `lead_vars` value. Having multiple would indicate an error in its calculation and could affect the results in `test_sig_loci_calculation_Enformer`
	"""
	for enformer_scoring_system in ['max','summed']:
		sig_variants_pval = pd.read_csv(f"/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/PsychENCODE_GWAS_scripts/EDA/PsychENCODE_Enformer_EDA/EnformerSigPVals_{scoring_system}_Variants_PsychENCODE_GWAS.csv")
		std_val_cols = [x for x in list(sig_variants_pval.columns) if 'std_val' in x]
		for std_val_col in std_val_cols:
			assert all(sig_variants_pval.groupby('lead_var')[std_val_col].nunique() == 1), f"You have multiple std_val values in Enformer {std_val_col}"

def test_sig_loci_calculation_Enformer():
	"""
	Tests that the identification of variants that had much more significnat tracks than others in the same disease associated locus did not lead to unexpected results
	"""

	for enformer_scoring_system in ['max','summed']:
		sig_variants_pval = pd.read_csv(f"/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/PsychENCODE_GWAS_scripts/EDA/PsychENCODE_Enformer_EDA/EnformerSigPVals_{scoring_system}_Variants_PsychENCODE_GWAS.csv")
		#dict containing the columns to identify most extreme variants per locus (keys) and the cols used to help calculate them (values) 
		identifier_cols = {'sig_in_locus_all_brain_feats_all_vars':'num_sig_all_brain_feats',
							'sig_in_locus_CAGE_feats_all_vars':'num_sig_CAGE_brain_feats',
							'sig_in_locus_all_brain_feats_SNP_only': 'num_sig_all_brain_feats',
							'sig_in_locus_CAGE_feats_SNP_only':'num_sig_CAGE_brain_feats'}
		sig_loci_cols = [x for x in list(sig_variants_pval.columns) if 'locus' in x and 'std_val' not in x]
		unique_loci = list(sig_variants_pval['lead_var'].unique())
		for unique_locus in unique_loci: #test all disease associated loci
			for sig_locus_col in sig_loci_cols:
				if 'SNP_only' in sig_locus_col:
					random_test_df = sig_variants_pval[(sig_variants_pval['lead_var'] == unique_locus) & sig_variants_pval['is_SNP'] == True]
				else:
					random_test_df = sig_variants_pval[(sig_variants_pval['lead_var'] == unique_locus)]
				identifier_calculating_col = identifier_cols[sig_locus_col]
				mean_sig = random_test_df[identifier_calculating_col].mean()
				std_sig = random_test_df[identifier_calculating_col].std()
				sig_var_rows = random_test_df[random_test_df[sig_locus_col]==True]
				std_val = sig_var_rows[f"std_val_{sig_locus_col}"].values[0] # test_consistent_std_val_per_identifier_col ensures there is only one value in this column per locus 

				# if there are -1 vals in the corresponding std_val col the tests are different and reflect the meaning behind the -1 value
				if random_test_df[random_test_df[f"std_val_{sig_locus_col}"] == -1].shape[0] >=1: 
					assert sig_var_rows.shape[0] == 1, f"if std_val_{sig_locus_col} has value -1 there should only be one putative causal variant in {sig_locus_col}!"
					assert sig_var_rows[identifier_calculating_col].values[0] == random_test_df[identifier_calculating_col].max(),f"if std_val_{sig_locus_col} has value -1 the value for the putative causal variant in this locus should be the largest value in {identifier_calculating_col}!"
				# if all the identifiger calculating column vals equal 0, the mean and std both equal 0 and the test will fail
					if not all(random_test_df[identifier_calculating_col] == 0):
						assert sig_var_rows[identifier_calculating_col].values[0] < mean_sig + std_sig, f"if std_val_{sig_locus_col} has value -1 the value for the putative causal variant in this locus should be the less than the mean + std of significant tracks across all variants in this locus, else the std_Val would be >=1!"
				else:
					assert all(sig_var_rows[identifier_calculating_col] >= (std_val*std_sig) + mean_sig), f"Putative causal variants in {sig_locus_col} do not all have {std_val}*std + mean {identifier_calculating_col}!"
					assert not all(sig_var_rows[identifier_calculating_col] >= ((std_val+1)*std_sig) + mean_sig), f"Putative causal variants in {sig_locus_col} have (std_val+1)*std + mean {identifier_calculating_col}. If that is the case then std_val should equal std_val+1! Maybe std_val == 30, which is the max std_val considered?"
					assert not all(random_test_df[identifier_calculating_col] == 0), f"All {identifier_calculating_col} == 0, but std_val != -1. std_val should equal -1!"


def test_sigDNV_Akita():
	"""
	This tests that the variants identified as significant truly are significant under the null hypothesis used in the Akita sig DNV script
	"""
	pass
	# score_cols = ['0_HFF', '1_H1hESC', '2_GM12878', '3_IMR90','4_HCT116']

	# def _get_quantiles(threshold, null_distribution,score_cols):
	# 	"""
	# 	Helper function used to get significance threshold 
	# 	"""
	    
	#     upper_percentiles = defaultdict() #only looking for big differences between ref and alt. No need for lower percentiles; not a matter of upregulation or downregulation, but is the contact highly different or not?
	#     for score_col in score_cols:
	#         null_distribution[score_col] = np.log(null_distribution[score_col])
	#         upper_percentiles[score_col] = null_distribution[score_col].quantile(threshold)
	#     return upper_percentiles 
	    
	# upper_percentiles = get_quantiles(0.999,null_dist,score_cols)


	# #load null
	# null_path = "/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/Akita/Rare1KGVariants"
	# null_dist = pd.DataFrame()
	# for idx, file in enumerate(os.listdir(null_path)):
	#   if scoring_system in file and file.endswith('.vcf.csv'):
	#     split_null_df = pd.read_csv(os.path.join(null_path,file))
	#     null_dist = null_dist.append(split_null_df)

	# #load DNV results
	# AkitaEDA_path = "/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/PsychENCODE_GWAS_scripts/EDA/PsychENCODE_Akita_EDA"
	# for file in os.listdir(AkitaEDA_path):
	# 	filepath = os.path.join(AkitaEDA_path,file)
	# 	if 'DNV' in filepath and filepath.endswith('.csv') and 'PVals' in filepath:
	# 		DNV_result_df = pd.read_csv(filepath)

	# 		#get 5 random significant variants
	# 		random_DNV_result_df = DNV_result_df.sample(5)

	# 		#get the results that hve the original Akita scores, not the PValues
	# 		non_pvals_result_file = filepath.replace('PVals','') 
	# 		non_pvals_DNV_result_df = pd.read_csv(non_pvals_result_file)

	# 		non_pvals_DNV_result_df = non_pvals_DNV_result_df.merge(random_DNV_result_df,on = ['chrom','pos','ref','alt']) #get the same 5 random variants for this df
	# 		for score_col in score_cols:
	# 			assert non_pvals_DNV_result_df

	# 		elif 'PVals' in filepath:
	# 			pass




def test_sigDNV_Sei():
	"""
	This tests that the variants identified as significant truly are significant under the null hypothesis used in the Sei sig DNV script
	"""
	pass

def test_sigDNV_Enformer():
	"""
	This tests that the variants identified as significant truly are significant under the null hypothesis used in the Enformer sig DNV script
	"""
	pass


def test_results_have_all_columns():
	"""
	Do this for all null distributions, variant predictions, and EDA outputs for all models. Make sure all columns exist
	"""

##### BELOW CODE IS TO CREATE A FILE CALLED `missing_sequence_class_scores.sh` so I can run `sh missing_sequence_class_scores.sh` to re-run sei predictions for those that are missing

		# sei_dir = "/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/Sei"
		# job_script="/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/PsychENCODE_GWAS_scripts/SeiScripts/SeiPsychENCODEGWASPredictions_job_script.sh"
		# vcf_list_dir = "/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/getGWASVariants/GWAS_psychENCODE_LeadTagVariants" #this is where all the split vcfs (len <=5000 each) were stored. I should get a sequence_class_scores.tsv for each of these
		# missing_sequence_class_scores =[]
		# for disease in os.listdir(sei_dir):
		# 	for vcf_or_folder in os.listdir(os.path.join(vcf_list_dir,disease)):
		# 		if vcf_or_folder.endswith(".vcf"):
		# 			vcf = vcf_or_folder
		# 			if not os.path.exists(os.path.join(sei_dir,disease,vcf,"sequence_class_scores.tsv")):
		# 				missing_sequence_class_scores.append(f"{disease};{vcf}")

		# with open("missing_sequence_class_scores.sh",'w') as f:
		# 	for missing_sequence_class_score in missing_sequence_class_scores:
		# 		disease = missing_sequence_class_score.split(';')[0]
		# 		vcf = missing_sequence_class_score.split(';')[1]
		# 		f.write(f"cd {sei_dir}")
		# 		f.write('\n')
		# 		f.write(f"qsub -cwd {job_script} {vcf} {os.path.join(vcf_list_dir,disease,vcf)} {disease}")
		# 		f.write('\n')
		# 		f.write('\n')



# ### BELOW CODE IS TO CREATE A FILE CALLED `missing_enformer_scores.sh` so I can run `sh missing_enformer_scores.sh` to re-run Enformer predictions for those that are missing
# import os
# enformer_dir = "/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/Enformer"
# job_script="/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/PsychENCODE_GWAS_scripts/AkitaEnformer_Scripts/AkitaEnformerPsychENCODEGWASPredictions_job_script.sh"
# vcf_list_dir = "/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/getGWASVariants/GWAS_psychENCODE_LeadTagVariants" #this is where all the split vcfs (len <=5000 each) were stored. I should get a sequence_class_scores.tsv for each of these
# missing_max_scores =[]
# missing_summed_scores = []
# for disease in os.listdir(enformer_dir):
# 	if disease != 'test': #files for use in unit test exist in here. Ignore these
# 		for vcf in os.listdir(os.path.join(vcf_list_dir,disease,'AkitaEnformer')): #get vcf's meant for use with Akita and Enformer (these have no header)
# 			if vcf.endswith(".vcf"):
# 				if not os.path.exists(os.path.join(enformer_dir,disease,f"Enformer_max_predictions_PsychENCODEGWASLeadTagVariants_{vcf}.csv")):
# 					missing_max_scores.append(f"{disease};{vcf}")
# 				if not os.path.exists(os.path.join(enformer_dir,disease,f"Enformer_summed_predictions_PsychENCODEGWASLeadTagVariants_{vcf}.csv")):
# 					missing_summed_scores.append(f"{disease};{vcf}")
# missing_data_list = list(set(missing_max_scores+missing_summed_scores)) #get unique entries from both lists. This new list stores disease/vcf for which there is either/both missing max or summed scores
# with open("/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/PsychENCODE_GWAS_scripts/test/missing_enformer_scores.sh",'w') as f:
# 	for missing_data in missing_data_list:
# 		disease = missing_data.split(';')[0]
# 		vcf = missing_data.split(';')[1]
# 		f.write(f"cd {os.path.join(enformer_dir,disease)}")
# 		f.write('\n')
# 		f.write(f"echo Redo of {disease}-{vcf} job >> AkitaPsychENCODEGWASPredictions_job.o")
# 		f.write('\n')
# 		f.write(f"qsub -cwd {job_script} Enformer {os.path.join(vcf_list_dir,disease,'AkitaEnformer',vcf)} {disease}")
# 		f.write('\n')
