#!/bin/bash

. /pollard/data/projects/sdrusinsky/pollard_lab/variant_modeling/bin/activate #activate env

job_script=/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/PsychENCODE_GWAS_scripts/EDA/PsychENCODE_Akita_EDA/FindAkitaSigVariants_job_script.sh

var_type=$1 #either DNV or GWAS. Case sensitive
shift #shift all arguments to the left so you can retrieve the rest as an array
diseases=("$@") #save diseases as an array
num_diseases=$(echo "${diseases[@]}" | wc -w)
total_num_tasks=$((2*$num_diseases)) #two tasks for each disease; one for msd and one for max Akita Analysis

qsub -cwd -t 1-$total_num_tasks -tc 10 -N AkitaFindSigVariants $job_script $var_type ${diseases[*]}



# example command:
# diseases=(PTSD Autism)
# sh run_AkitaSigVariantsJob.sh DNV ${diseases[*]}