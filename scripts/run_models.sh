reference_genome=$1 #either hg19 or hg38
akita_enformer_vcf_path=$2 #path to vcf file formatted for input into Akita/Enformer
sei_vcf_path=$3 #path to vcf file formatted for input into Sei
experiment_name=$4

top_level_dir=$(git rev-parse --show-toplevel)


if [[ $reference_genome = 'hg38' ]]; then #download hg38 if file doesn't exist
	fasta_file=$top_level_dir/models/hg38_genome.fa
    !wget -O - http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz | gunzip -c > {fasta_file}
elif [[ $reference_genome = 'hg19' ]]; then
	fasta_file=$top_level_dir/models/hg19_genome.fa
    !wget -O - http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz | gunzip -c > {fasta_file}
fi

#run Enformer and re-route stderr and stdout
python3 $top_level_dir/scripts/EnformerScoring.py $akita_enformer_vcf_path $reference_genome $experiment_name 1>>$top_level_dir/model_outputs/Enformer/$experiment_name/log.o 2>>$top_level_dir/model_outputs/Enformer/$experiment_name/log.e

#run Akita and re-route stderr and stdout
python3 $top_level_dir/scripts/AkitaScoring.py $akita_enformer_vcf_path $reference_genome $experiment_name 1>>$top_level_dir/model_outputs/Akita/$experiment_name/log.o 2>>$top_level_dir/model_outputs/Akita/$experiment_name/log.e

#run Sei. Sei Script re-routes stderr and stdout already
sh $top_level_dir/scripts/SeiScoring.sh $sei_vcf_path $reference_genome $experiment_name 