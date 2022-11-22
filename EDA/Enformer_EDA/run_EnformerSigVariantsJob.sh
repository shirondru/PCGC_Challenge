#!/bin/bash

. /pollard/data/projects/sdrusinsky/pollard_lab/variant_modeling/bin/activate #activate env

top_level_dir=$(git rev-parse --show-toplevel)
job_script=$top_level_dir/EDA/Enformer_EDA/FindEnformerSigVariants_job_script.sh


relevant_cols_txt=$1
echo $relevant_cols_txt
shift #shift all arguments to the left so you can retrieve the diseases array
diseases=("$@") #save diseases as an array
num_diseases=$(echo "${diseases[@]}" | wc -w)
total_num_tasks=$((2*$num_diseases)) #two tasks for each disease; one for summed and one for max Enformer Analysis

if [[ ! -d "$top_level_dir/EDA/Enformer_EDA/logs" ]]; then
	mkdir logs
fi


qsub -cwd -t 1-$total_num_tasks -tc 10 -N EnformerFindSigVariants $job_script $relevant_cols_txt $top_level_dir ${diseases[*]}

# example command:
# diseases=(PTSD Autism)
# relevant_cols_txt=NeuroEnformerCols.txt
# sh run_EnformerSigVariantsJob.sh $relevant_cols_txt ${diseases[*]}