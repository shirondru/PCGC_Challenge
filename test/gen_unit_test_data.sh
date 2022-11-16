"""
************ Generate data for test.py/test_Enformer_predictions *********************

This does not contain code to generate `ground_truth_summed_Enformer_preds.csv`, which was generated using:
	/Users/shirondrusinsky/Documents/Pollard_Lab/PsychENCODE_GWAS_scripts/test/data_for_tests/testEnformerPredictions/enformer_usage_get_ground_truth_summed_Enformer_preds.ipynb

"""
cd data_for_tests/test_Enformer_predictions_data
python3 gen_testEnformer_vcfs_1.py #generate ground_truth_Enformer_SNPs.vcf
qsub -cwd run_EnformerTestPredictions_job_2.sh #generate Enformer predictions using my code for SNPs in ground_truth_Enformer_SNPs.vcf




"""
************** Generate data for test.py/test_Akita_predictions ****************
This creates Akita predictions using my code for Variants_KatieG_formed_AkitaPreds.vcf, 
which are the same variants for which predictions were made by another member of the Pollard Lab for a unit test. 
Those predictions are stored in AkitaPreds_by_KatieG.tsv
"""
cd /Users/shirondrusinsky/Documents/Pollard_Lab/PsychENCODE_GWAS_scripts/test/data_for_tests/AkitaPreds_by_KatieG
sh gen_myAkitaPreds_for_KatieG_variants.sh