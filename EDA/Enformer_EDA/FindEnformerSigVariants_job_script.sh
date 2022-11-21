#!/bin/bash
#$ -o logs/FindEnformerSigVariants_$JOB_ID_$TASK_ID.o
#$ -e logs/FindEnformerSigVariants_$JOB_ID_$TASK_ID.e                       
#$ -cwd                            
#$ -r y                            
#$ -l mem_free=150G                  
#$ -l scratch=50G                   
#$ -l h_rt=133:00:00                
                                   
. /pollard/data/projects/sdrusinsky/pollard_lab/variant_modeling/bin/activate
top_level_dir=$(git rev-parse --show-toplevel)

var_type=$1
shift
relevant_cols=$1
shift
diseases=("$@") #save  arguments as an array

#### Take diseases array and copy each diseaes inside of it twice and save into a new array called `duplicated_disease_arr` so you can form max and summed
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
#form prediction with max scoring system if SGE_TASK_ID  is even, else use summed scoring system
if [ `expr  $SGE_TASK_ID % 2` == 0 ]; then 
	scoring_system=max
fi

if [ `expr  $SGE_TASK_ID % 2` != 0 ]; then 
	scoring_system=summed	
fi

echo $scoring_system
echo $relevant_cols
echo $task_disease

if [ "$var_type" == "GWAS" ]; then 
	echo Finding Sig GWAS SNPs 
	python3 $top_level_dir/EDA/Enformer_EDA/FindEnformerSigVariants.py --scoring_system $scoring_system --relevant_cols $relevant_cols --disease $task_disease

fi

if [ "$var_type" == "DNV" ]; then 
	echo Finding Sig DNVs
	python3 $top_level_dir/EDA/Enformer_EDA/FindEnformerSigDNVs.py --scoring_system $scoring_system --relevant_cols $relevant_cols --disease $task_disease
fi


qstat -j $JOB_ID

