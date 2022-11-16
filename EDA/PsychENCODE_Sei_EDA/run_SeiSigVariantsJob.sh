#!/bin/bash

. /pollard/data/projects/sdrusinsky/pollard_lab/variant_modeling/bin/activate #activate env

job_script=/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/PsychENCODE_GWAS_scripts/EDA/PsychENCODE_Sei_EDA/FindSeiSigVariants_job_script.sh


var_type=$1 #either DNV or GWAS. Case sensitive
shift #shift all arguments to the left so you can retrieve the rest as an array
diseases=("$@") #save diseases as an array
total_num_tasks=$(echo "${diseases[@]}" | wc -w) #one task per disease

qsub -cwd -t 1-$total_num_tasks -tc 10 -N SeiFindSigVariants $job_script $var_type ${diseases[*]}



# Example command
# diseases=(PTDS Autism)
# var_type=DNV
# sh run_SeiSigDNVsJob.sh $var_type ${diseases[*]}