reference_genome=$1 #either hg19 or hg38
vcf_path=$2
experiment_name=$3

top_level_dir=$(git rev-parse --show-toplevel)

#run Enformer and re-route stderr and stdout
python3 $top_level_dir/scripts/EnformerScoring.py $vcf_path $reference_genome $experiment_name 1>>$top_level_dir/model_outputs/Enformer/$experiment_name/log.o 2>>$top_level_dir/model_outputs/Enformer/$experiment_name/log.e

#run Akita and re-route stderr and stdout
python3 $top_level_dir/scripts/AkitaScoring.py $vcf_path $reference_genome $experiment_name 1>>$top_level_dir/model_outputs/Akita/$experiment_name/log.o 2>>$top_level_dir/model_outputs/Akita/$experiment_name/log.e

#run Sei. Sei Script re-routes stderr and stdout already
sh $top_level_dir/scripts/SeiScoring.sh $vcf_path $reference_genome $experiment_name 