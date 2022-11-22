#!/bin/bash
#$ -o logs/FindSeiSigVariants_$JOB_ID_$TASK_ID.o
#$ -e logs/FindSeiSigVariants_$JOB_ID_$TASK_ID.e                       
#$ -cwd                            
#$ -r y                            
#$ -l mem_free=150G                  
#$ -l scratch=25G                   
#$ -l h_rt=133:00:00                
                                   
. /pollard/data/projects/sdrusinsky/pollard_lab/variant_modeling/bin/activate #activate env
top_level_dir=$1
shift
diseases=("$@") #save  arguments as an array
echo $SGE_TASK_ID
task_disease="${diseases[$SGE_TASK_ID - 1]}"
echo $task_disease


python3 $top_level_dir/EDA/Sei_EDA/FindSeiSigDNVs.py --experiment_name $task_disease --git_root $top_level_dir

qstat -j $JOB_ID
