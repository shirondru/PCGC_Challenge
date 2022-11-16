#!/bin/bash
                    
#$ -o get_rare_1kg_variants.o
#$ -e get_rare_1kg_variants.e                       
#$ -cwd                            
#$ -r y                            
#$ -l mem_free=100G                  
#$ -l scratch=25G                   
#$ -l h_rt=8:00:00                
##$ -t 1-10                        
                                   
. /pollard/data/projects/sdrusinsky/pollard_lab/variant_modeling/bin/activate

python3 /pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/PsychENCODE_GWAS_scripts/DNVs/get_rare_1kg_vars.py
qstat -j $JOB_ID
