sei_website_preds_for_random_PsychENCODE_variants.tsv:
Predictions for the first 5 variants in PsychENCODE_GWASVariants_BipolarDisorder_SeiNoHeader.vcf, last 5 variants in PsychENCODE_GWASVariants_ADHD_SeiNoHeader.vcf and first 5 variants in PsychENCODE_GWASVariants_Alzheimer_SeiNoHeader.vcf were computed using the website hosting the Sei framework (https://hb.flatironinstitute.org/sei). The predictions from the website were treated as the ground truth, proper way to use Sei. So the 15 resulting predictions were compared against the predictions received using my code in a unit test in test.py to ensure proper predictions were formed.

sei_manuscripts_HGMD_preds.xlsx:
This file comes from Supplementary File 6 of the Sei manuscript. It contains the sei predictions for some HGMD mutations. I performed predictions using my code with Sei to confirm the same results were achieved and Sei is being used properly

AkitaPreds_by_KatieG:
	AkitaPreds_by_KatieG.tsv:
	This file contains Akita variant predictions made by Katie Gjoni from the Pollard Lab. I form Akita predictions for the same variants using my code and test the results are the same.

	Variants_KatieG_formed_AkitaPreds.vcf:
	The same variants as in AkitaPreds_by_KatieG.tsv, put into vcf formatting to run with Akita

test_Enformer_predictions_data:
Code and data necessary to run unit test to ensure Enformer predictions were formed properly
