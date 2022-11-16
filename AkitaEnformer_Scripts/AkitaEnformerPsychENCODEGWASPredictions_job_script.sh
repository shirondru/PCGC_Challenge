#!/bin/bash
#$ -o /pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/PsychENCODE_GWAS_scripts/AkitaEnformer_Scripts/logs/stdout/$JOB_ID.o
#$ -e /pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/PsychENCODE_GWAS_scripts/AkitaEnformer_Scripts/logs/stderr/$JOB_ID.e
#$ -r y                                                        
#$ -l gpu_mem=20000M                 
#$ -l scratch=15G                  
#$ -l h_rt=07:00:00                
                                   
. /pollard/data/projects/sdrusinsky/pollard_lab/variant_modeling/bin/activate


	
model=$1 #either Akita or Enformer. Case sensitive
reference_genome=$2 #either hg19 or hg38
split_vcf_list_path=$3 
split_vcf_file=$(cat $split_vcf_list_path) #open txt file containing split vcf file paths and save as string
read -a split_vcf_list <<< $split_vcf_file #convert into a parse-able array

vcf_path="${split_vcf_list[$SGE_TASK_ID - 1]}" #Get input vcf data file corresponding to this array job task
disease=$(echo $vcf_path | cut -d'/' -f10) #Get name of disease corresponding to the input vcf data file for this array job task
vcf_basename=$(echo $vcf_path | cut -d'/' -f12)
echo VCF_PATH: $vcf_path
echo DISEASE: $disease
echo VCF_BASENAME: $vcf_basename

output_log_dir=/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/$model/$disease/logs
mkdir $output_log_dir #ensure output directory exists
cd $output_log_dir
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/pollard/data/projects/sdrusinsky/pollard_lab/usr #add to LD_LIBRARY_PATH libcudnn and libcusolver.so
module load cuda/11.0
export CUDA_VISIBLE_DEVICES=$SGE_GPU

# GPU use monitoring
gpuprof=$(dcgmi group -c mygpus -a $SGE_GPU | awk '{print $10}')
dcgmi stats -g $gpuprof -e
dcgmi stats -g $gpuprof -s $JOB_ID

python3 /pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/PsychENCODE_GWAS_scripts/AkitaEnformer_Scripts/${model}PsychENCODE_GWAS_predictions.py --vcf_path $vcf_path --reference_genome $reference_genome --disease $disease  1>>$output_log_dir/$vcf_basename.$JOB_ID.$SGE_TASK_ID.o 2>>$output_log_dir/$vcf_basename.$JOB_ID.$SGE_TASK_ID.e


# GPU use monitoring
dcgmi stats -g $gpuprof -x $JOB_ID
dcgmi stats -g $gpuprof -v -j $JOB_ID
dcgmi group -d $gpuprof
qstat -j $JOB_ID
