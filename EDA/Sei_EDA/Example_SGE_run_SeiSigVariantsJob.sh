#!/bin/bash

. /pollard/data/projects/sdrusinsky/pollard_lab/variant_modeling/bin/activate #activate env
top_level_dir=$(git rev-parse --show-toplevel)
job_script=$top_level_dir/EDA/Sei_EDA/FindSeiSigVariants_job_script.sh



diseases=("$@") #save diseases as an array
total_num_tasks=$(echo "${diseases[@]}" | wc -w) #one task per disease

if [[ ! -d "$top_level_dir/EDA/Sei_EDA/logs" ]]; then
	mkdir logs
fi

qsub -cwd -t 1-$total_num_tasks -tc 10 -N SeiFindSigVariants $job_script $top_level_dir ${diseases[*]} 



# Example command
# diseases=(PTDS Autism)

# sh run_SeiSigVariantsJob.sh  ${diseases[*]}