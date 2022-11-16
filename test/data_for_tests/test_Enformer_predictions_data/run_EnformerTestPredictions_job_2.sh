#$ -o test_EnformerPredictions.o
#$ -e test_EnformerPredictions.e                       
#$ -cwd                            
#$ -r y                            
#$ -j y                            
#$ -l mem_free=10G                  
#$ -l scratch=25G                   
#$ -l h_rt=03:00:00                
##$ -t 1-10   



. /wynton/home/hernandez/shirondru/pollard_lab/pollard_environment/bin/activate

model=Enformer
vcf_path=/wynton/home/hernandez/shirondru/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/test/testEnformerPredictions/ground_truth_Enformer_SNPs.vcf
disease=test

python3 /wynton/home/hernandez/shirondru/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/AkitaEnformer_Scripts/${model}PsychENCODE_GWAS_predictions.py --vcf_path $vcf_path --disease $disease



qstat -j $JOB_ID


#qsub -cwd run_EnformerTestPredictions_job_2.sh