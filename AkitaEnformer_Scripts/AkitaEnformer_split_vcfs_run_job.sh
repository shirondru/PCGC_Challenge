#! /bin/bash
. /pollard/data/projects/sdrusinsky/pollard_lab/variant_modeling/bin/activate

model=$1 #either Akita or Enformer. Case sensitive
shift #shift all arguments to the left so you can retrieve the rest as an array
reference_genome=$1 #either hg19 or hg38
shift
diseases=("$@") #save remaining arguments as an array

data_dir=/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/getGWASVariants/GWAS_psychENCODE_LeadTagVariants
model_pred_dir=/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/${model}
job_script=/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/PsychENCODE_GWAS_scripts/AkitaEnformer_Scripts/AkitaEnformerPsychENCODEGWASPredictions_job_script.sh

# diseases=(ADHD TouretteSyndrome BipolarDisorder UnipolarDepression AnorexiaNervosa ObsessiveCompulsiveDisorder PTSD Schizophrenia Anxiety Alzheimer Autism)
# diseases=(RandomGWASVariants)
# diseases=(CHD_DNVs_Richter2020)
cd $data_dir
for disease in ${diseases[*]}; do
	#define directory the split vcf files will go into. Keep vcfs for use with Akita/Enformer separate from Sei, because the variants in each file are the same but Sei's lack a header. Keeping them together would result in both sets of file being run, leading to double the jobs with redundant variants being modeled
	outdir=$data_dir/$disease/AkitaEnformer 
	mkdir -p $outdir

	#define the directory the Akita/Enformer predictions will go into
	model_outdir=$model_pred_dir/$disease
	mkdir $model_outdir

	file=PsychENCODE_GWASVariants_${disease}.vcf #Use files with header
	outfile=$outdir/$file
	split -l 15000 -d --additional-suffix=.vcf  $file $outfile #split files into 5000 rows each to parallelize Akita/Enformer operations as different jobs
	cd $model_outdir #so the job script output txt file is saved in the appropriate folder
	for split_file_path in $outdir/*.vcf;do
		split_file_basename=$(basename -- "$split_file_path")
		# rm AkitaPsychENCODEGWASPredictions_job.o #remove output file so its clear what error and output messages come from this most recent submission of jobs, in case this gets run multiple times
	done
done


total_num_tasks=0 #Will contain number of SGE tasks. This is the number of split vcf files for which a task needs to be created
#make an array with each split vcf path
for disease in ${diseases[*]}; do
	for vcf in  $outdir/*vcf; do #count number of split vcf files for this disease, add to running tally
		total_num_tasks=$((total_num_tasks + 1))
		split_vcf_list+=(${vcf}) # create list of all vcf file paths
	done
done

split_vcf_list_path=$outdir/PsychENCODE_GWASVariants_${disease}_split_vcf_list.txt
echo "${split_vcf_list[@]}" > $split_vcf_list_path #save list of split vcf files as txt file that will be opened in job script and used to form predictions based on SGE task ID
qsub -q gpu.q -t 1-$total_num_tasks -tc 30 -N ${model}PsychENCODEGWASPredictions $job_script $model $reference_genome $split_vcf_list_path #run the job that calls Akita/Enformer


# sh AkitaEnformer_split_vcfs_run_job.sh Akita hg38 $diseases
