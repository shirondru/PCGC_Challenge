reference_genome=$1 #either hg19 or hg38
vcf_path=$2
experiment_name=$3

outdir=./model_outputs/Sei/$experiment_name #where model outputs will be stored

cd ./models/Sei/sei-framework

#run Sei. re-route stderr and stdout
sh run_pipeline.sh $vcf_path $reference_genome $outdir  1>>$oudir/log.o 2>>$oudir/log.e
