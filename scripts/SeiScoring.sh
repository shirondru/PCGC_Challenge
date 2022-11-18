vcf_path=$1
reference_genome=$2 #either hg19 or hg38
experiment_name=$3
vcf_basename=$(echo $vcf_path | rev | cut -d'/' -f 1 | rev | cut -d'.' -f1) #name of vcf file, minus ".vcf"
top_level_dir=$(git rev-parse --show-toplevel)
cd $top_level_dir #cd into top level directory (`PCGC_Challenge`)

outdir=$top_level_dir/model_outputs/Sei/$experiment_name #where model outputs will be stored

if [[ ! -e $outdir ]]; then #make outdir if it doesn't exist
    mkdir $outdir
fi

cd $top_level_dir/models/Sei/sei-framework

#run Sei. re-route stderr and stdout
sh 1_variant_effect_prediction.sh $vcf_path $reference_genome $outdir  1>>$outdir/log.o 2>>$outdir/log.e

ref_fp=$outdir/chromatin-profiles-hdf5/${vcf_basename}.ref_predictions.h5
alt_fp=$outdir/chromatin-profiles-hdf5/${vcf_basename}.alt_predictions.h5
sh 2_varianteffect_sc_score.sh $ref_fp $alt_fp $outdir 
