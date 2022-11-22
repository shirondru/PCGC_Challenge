#!/bin/bash
#$ -o logs/FindAkitaSigVariants_$JOB_ID_$TASK_ID.o
#$ -e logs/FindAkitaSigVariants_$JOB_ID_$TASK_ID.e                       
#$ -cwd                            
#$ -r y                            
#$ -l mem_free=100G                  
#$ -l scratch=25G                   
#$ -l h_rt=133:00:00                
##$ -t 1-10                        
                                   
. /pollard/data/projects/sdrusinsky/pollard_lab/variant_modeling/bin/activate
top_level_dir=$(git rev-parse --show-toplevel)

top_level_dir=$1
shift
diseases=("$@") #save  arguments as an array

#### Take diseases array and copy each diseaes inside of it twice and save into a new array called `duplicated_disease_arr` so you can form max and msd
#### predictions with each disease based on the $SGE_TASK_ID
duplicated_disease_arr=()
echo "${diseases[@]}"
for disease in "${diseases[@]}"; do
	duplicated_disease_arr+=($disease)
	duplicated_disease_arr+=($disease)
done
echo $SGE_TASK_ID
echo "${duplicated_disease_arr[@]}"
task_disease="${duplicated_disease_arr[$SGE_TASK_ID - 1]}"
echo $task_disease
#form prediction with max scoring system if SGE_TASK_ID  is even, else use msd scoring system
if [ `expr  $SGE_TASK_ID % 2` == 0 ]; then 
	scoring_system=max
fi

if [ `expr  $SGE_TASK_ID % 2` != 0 ]; then 
	scoring_system=msd	
fi


python3 $top_level_dir/EDA/Akita_EDA/FindAkitaSigDNVs.py --scoring_system $scoring_system --experiment_name $task_disease --git_root $top_level_dir



qstat -j $JOB_ID

