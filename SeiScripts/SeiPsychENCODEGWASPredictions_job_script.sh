#!/bin/bash
#$ -o /pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/PsychENCODE_GWAS_scripts/SeiScripts/logs/stdout/$JOB_ID.o
#$ -e /pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/PsychENCODE_GWAS_scripts/SeiScripts/logs/stderr/$JOB_ID.e
#$ -r y                                                        
#$ -l gpu_mem=20000M                 
#$ -l scratch=15G                   
#$ -l h_rt=07:00:00                
                                        
module load Sali
module load anaconda
eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)" #tell conda using bash so no need to conda init
conda activate /pollard/data/projects/sdrusinsky/pollard_lab/conda_variant_modeling #activate conda environment with proper pytorch and cuda installed
module load cuda/11.0
export CUDA_VISIBLE_DEVICES=$SGE_GPU

reference_genome=$1 #either hg19 or hg38
split_vcf_list_path=$2
split_vcf_file=$(cat $split_vcf_list_path) #open txt file containing split vcf file paths and save as string
read -a split_vcf_list <<< $split_vcf_file #convert into a parse-able array

vcf_path="${split_vcf_list[$SGE_TASK_ID - 1]}" #Get input vcf data file corresponding to this array job task
disease=$(echo $vcf_path | cut -d'/' -f10) #Get name of disease corresponding to the input vcf data file for this array job task
vcf_basename=$(echo $vcf_path | cut -d'/' -f11) #-f11 for Sei but -f12 for AkitaEnformer because of the extra dir
echo VCF_BASENAME
echo $vcf_basename
output_log_dir=/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/Sei/$disease/logs
mkdir -p $output_log_dir  #ensure output directory exists
cd $output_log_dir

tmp_outdir=$TMPDIR/$JOB_ID/$SGE_TASK_ID/Sei/$disease/$vcf_basename #location in scratch output will be temporarily stored
echo $TMPDIR
echo $tmp_outdir
mkdir -p $tmp_outdir
global_outdir=/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/Sei/$disease/$vcf_basename #global location data will be stored outside of scratch
mkdir -p $global_outdir




cd /pollard/data/projects/sdrusinsky/pollard_lab/sei-framework

# GPU use monitoring
gpuprof=$(dcgmi group -c mygpus -a $SGE_GPU | awk '{print $10}')
dcgmi stats -g $gpuprof -e
dcgmi stats -g $gpuprof -s $JOB_ID

sh run_pipeline.sh $vcf_path $reference_genome $tmp_outdir --cuda 1>>$output_log_dir/$vcf_basename.$JOB_ID.$SGE_TASK_ID.o 2>>$output_log_dir/$vcf_basename.$JOB_ID.$SGE_TASK_ID.e

mv $tmp_outdir/* $global_outdir #move data out of scratch to permanent storage

# GPU use monitoring
dcgmi stats -g $gpuprof -x $JOB_ID
dcgmi stats -g $gpuprof -v -j $JOB_ID
dcgmi group -d $gpuprof
qstat -j $JOB_ID

# cd $sei_outdir
# rm -r chromatin-profiles-hdf5 #remove this directory because not needed for our purpose and requires a lot of storage
# rm chromatin_profile_diffs.tsv #this also requires a lot of storage
