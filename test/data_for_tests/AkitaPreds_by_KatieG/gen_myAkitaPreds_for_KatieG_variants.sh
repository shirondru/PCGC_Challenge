. /wynton/home/hernandez/shirondru/pollard_lab/pollard_environment/bin/activate

model=Akita
data_dir=/wynton/home/hernandez/shirondru/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/PsychENCODE_GWAS_scripts/test/data_for_tests/AkitaPreds_by_KatieG
model_pred_dir=/wynton/home/hernandez/shirondru/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/Akita
cd $data_dir
job_script=/wynton/home/hernandez/shirondru/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/PsychENCODE_GWAS_scripts/AkitaEnformer_Scripts/AkitaEnformerPsychENCODEGWASPredictions_job_script.sh
disease=test

	

#define the directory the Akita predictions will go into
model_outdir=$model_pred_dir/$disease
mkdir $model_outdir
vcf_path=$data_dir/Variants_KatieG_formed_AkitaPreds.vcf #Variants KatieG has formed predictions with
outfile=$outdir/$file
cd $model_outdir #so the job script output txt file is saved in the appropriate folder
	
qsub -cwd $job_script $model $vcf_path $disease #run the job that calls Akita/Enformer




# sh AkitaEnformer_split_vcfs_run_job.sh Akita
