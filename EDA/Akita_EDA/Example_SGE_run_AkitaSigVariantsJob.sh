#!/bin/bash

. /pollard/data/projects/sdrusinsky/pollard_lab/variant_modeling/bin/activate #activate env
top_level_dir=$(git rev-parse --show-toplevel)
job_script=$top_level_dir/EDA/Akita_EDA/FindAkitaSigVariants_job_script.sh


diseases=("$@") #save diseases as an array
num_diseases=$(echo "${diseases[@]}" | wc -w)
total_num_tasks=$((2*$num_diseases)) #two tasks for each disease; one for msd and one for max Akita Analysis


if [[ ! -d "$top_level_dir/EDA/Akita_EDA/logs" ]]; then
	mkdir logs
fi

qsub -cwd -t 1-$total_num_tasks -tc 10 -N AkitaFindSigVariants $job_script $top_level_dir ${diseases[*]}



# example command:
# diseases=(PTSD Autism)
# sh Example_SGE_run_AkitaSigVariantsJob.sh ${diseases[*]}