#!/bin/bash

. /pollard/data/projects/sdrusinsky/pollard_lab/variant_modeling/bin/activate #activate env

job_script=/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/PsychENCODE_GWAS_scripts/EDA/PsychENCODE_Enformer_EDA/FindEnformerSigVariants_job_script.sh

var_type=$1 #either DNV or GWAS. Case sensitive
echo $var_type
shift #shift all arguments to the left 
relevant_cols_txt=$1
echo $relevant_cols_txt
shift #shift all arguments to the left so you can retrieve the diseases array
diseases=("$@") #save diseases as an array
num_diseases=$(echo "${diseases[@]}" | wc -w)
total_num_tasks=$((2*$num_diseases)) #two tasks for each disease; one for summed and one for max Enformer Analysis




qsub -cwd -t 1-$total_num_tasks -tc 10 -N EnformerFindSigVariants $job_script $var_type $relevant_cols_txt ${diseases[*]}

# example command:
# diseases=(PTSD Autism)
# relevant_cols_txt=NeuroEnformerCols.txt
# sh run_EnformerSigVariantsJob.sh DNV $relevant_cols_txt ${diseases[*]}