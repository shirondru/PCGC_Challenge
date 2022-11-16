. /pollard/data/projects/sdrusinsky/pollard_lab/variant_modeling/bin/activate

data_dir=/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/getGWASVariants/GWAS_psychENCODE_LeadTagVariants
sei_pred_dir=/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/Sei
cd $data_dir
job_script=/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/PsychENCODE_GWAS_scripts/SeiScripts/SeiPsychENCODEGWASPredictions_job_script.sh
reference_genome=$1 #either hg19 or hg38
shift
diseases=("$@") #save remaining arguments as an array
# diseases=(ADHD TouretteSyndrome BipolarDisorder UnipolarDepression AnorexiaNervosa ObsessiveCompulsiveDisorder PTSD Schizophrenia Anxiety Alzheimer Autism)
# diseases=(RandomGWASVariants)
# diseases=(PTSD)
for disease in ${diseases[*]}; do
	#define directory the split vcf files will go into
	outdir=$data_dir/$disease
	mkdir $outdir

	#define the directory the sei predictions will go into
	sei_outdir=$sei_pred_dir/$disease
	mkdir $sei_outdir

	file=PsychENCODE_GWASVariants_${disease}_SeiNoHeader.vcf #only consider files without Sei header
	outfile=$outdir/$file
	split -l 15000 -d --additional-suffix=.vcf  $file $outfile #split files into 5000 rows each to parallelize Sei operations as different jobs
	cd $sei_outdir #so the job script output txt file is saved in the appropriate folder
	for split_file_path in $outdir/*.vcf;do
		split_file_basename=$(basename -- "$split_file_path")
		# rm SeiPsychENCODEGWASPredictions_job.o #remove output file so its clear what error and output messages come from thismost recent submission of jobs, in case this gets run multiple times
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

qsub -q gpu.q -t 1-$total_num_tasks -tc 30 -N SeiPsychENCODEGWASPredictions $job_script $reference_genome $split_vcf_list_path #run the job that calls Sei


#example command
# diseases=(PTSD)
# sh Sei_split_vcfs_run_job.sh hg38 $diseases


##### After this ran, I manually removed chromatin-profiles-hdf5 folders and chromatin_profile_diffs.tsv files to save storage using the following code####
# sei_pred_dir=/pollard/data/projects/sdrusinsky/pollard_lab/GWASPredictions/PsychENCODE_GWAS_Predictions/Sei
# cd $sei_pred_dir
# for disease in ADHD TouretteSyndrome BipolarDisorder UnipolarDepression AnorexiaNervosa ObsessiveCompulsiveDisorder PTSD Schizophrenia Anxiety Alzheimer Autism; do	
# 	for vcf_file in *.vcf; do
# 		rm $disease/$vcf_file/chromatin_profile_diffs.tsv
# 		rm -r $disease/$vcf_file/chromatin-profiles-hdf5
# 	done
# done
