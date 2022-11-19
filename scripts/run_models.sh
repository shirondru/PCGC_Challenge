reference_genome=$1 #either hg19 or hg38
akita_enformer_vcf_path=$2 #path to vcf file formatted for input into Akita/Enformer
sei_vcf_path=$3 #path to vcf file formatted for input into Sei
experiment_name=$4

top_level_dir=$(git rev-parse --show-toplevel) #path to the top level of the git hub repo (i.e, */PCGC_Challenge)

## download the reference genome (either hg19, or hg38) if it doesn't already exist. Define path to the fasta file
fasta_file=$top_level_dir/models/${reference_genome}_genome.fa
if [[ ! -e $fasta_file ]]; then #download reference genome fasta file if it doesn't exist
	wget -O - http://hgdownload.cse.ucsc.edu/goldenPath/${reference_genome}/bigZips/${reference_genome}.fa.gz | gunzip -c > ${fasta_file}
fi

#run Enformer and re-route stderr and stdout
python3 $top_level_dir/scripts/EnformerScoring.py --vcf_path $akita_enformer_vcf_path --reference_genome $reference_genome --experiment_name $experiment_name 1>>$top_level_dir/model_outputs/Enformer/$experiment_name/log.o 2>>$top_level_dir/model_outputs/Enformer/$experiment_name/log.e

#run Akita and re-route stderr and stdout
python3 $top_level_dir/scripts/AkitaScoring.py --vcf_path $akita_enformer_vcf_path --reference_genome $reference_genome --experiment_name $experiment_name 1>>$top_level_dir/model_outputs/Akita/$experiment_name/log.o 2>>$top_level_dir/model_outputs/Akita/$experiment_name/log.e

#run Sei. Sei Script re-routes stderr and stdout already
sh $top_level_dir/scripts/SeiScoring.sh $sei_vcf_path $reference_genome $experiment_name 
