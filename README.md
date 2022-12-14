# PCGC Challenge Scripts
 Scripts to run deep learning models that take DNA sequence as input and predict gene regulatory profiles, to prioritize causal variants and mechanisms

Models include Akita, Sei, an Enformer


### Requirements:
1. Clone repository containing Akita:

	``` 
	cd ./models/Akita
	git clone https://github.com/calico/basenji.git
	cd basenji/manuscripts/akita
	sh get_model.sh
	```
2. Clone repository containing Sei:
	``` 
	cd ./models/Sei
	git clone https://github.com/FunctionLab/sei-framework.git
	cd sei-framework
	sh ./download_data.sh
	```

(OPTIONAL). Download Enformer: Only necessary if running code from a machine without internet access. Otherwise, Enformer can be loaded from tfhub. Code in this repo loads Enformer from tfhub, so downloading the model is optional. If downloading Enformer is desired, the path to Enformer must be updated in all scripts that call Enformer.
* Use either link:
	* https://tfhub.dev/deepmind/enformer/1
	* https://github.com/deepmind/deepmind-research/tree/master/enformer

3. Create a virtual environment
	```
	virtualenv -p python3 variant_modeling #create a virtual environment called variant_modeling
	```
4. Activate environment and install requirements
	```
	cd variant_modeling
	. bin/activate
	python3 -m pip install -r requirements.txt
	```



## Usage:
* Define reference genome (either hg38 or hg19)
* Define the name of the experiment, so output files are recognizable
* Define path to the vcf files containing the variants to be scored. Files must be in VCF format. Two VCF files are required per set of variants. One with the header removed, for compatibility with Sei, and one with the header intact, for compatibility with the parser used for Akita and Sei.
* Run script that scores the variants with the three models.
```
reference_genome=hg38 #either hg38 or hg19
experiment_name=test #your choice. This will help determine the names of output files
top_level_dir=$(git rev-parse --show-toplevel) #path for top level of this repo
akita_enformer_vcf_paths=$top_level_dir/DNVs/example.vcf #path to vcf file formatted for input into Akita/Enformer
sei_vcf_paths=$top_level_dir/DNVs/example_sei.vcf #path to vcf file formatted for input into Sei. Contains same variants as akita_enformer_vcf_paths

## run models
sh $top_level_dir/scripts/run_models.sh $reference_genome $akita_enformer_vcf_paths $sei_vcf_paths $experiment_name 


```

However, you will probably need to use a cluster if scoring more than just a few variants, which complicates things as different clusters run differently. Here is an example of a script that could be used to submit a job using SGE (with CPUs, more code will need to be added to use GPUs):

```
#!/bin/bash

top_level_dir=$(git rev-parse --show-toplevel) #path for top level of this repo


#$ -o $top_level_dir/scripts/$JOB_ID.o #path to save stdout
#$ -e $top_level_dir/scripts$JOB_ID.e #path to save stderr
#$ -r y  # if job crashes, it should be restarted                                                        
#$ -l mem_free=15G                 
#$ -l scratch=15G                  
#$ -l h_rt=24:00:00  #max 24hr run time           
                                   
. $top_level_dir/variant_modeling/bin/activate #activate environment



reference_genome=hg38
akita_enformer_vcf_paths=$top_level_dir/DNVs/example.vcf #path to vcf file formatted for input into Akita/Enformer
sei_vcf_paths=$top_level_dir/DNVs/example_sei.vcf #path to vcf file formatted for input into Sei. Contains same variants as akita_enformer_vcf_paths
experiment_name=test

qsub -cwd $top_level_dir/scripts/run_models.sh $reference_genome $akita_enformer_vcf_paths $sei_vcf_paths $experiment_name #submit the job

```
This is a simple job script that submits run_models.sh as a job that uses CPUs. For faster jobs, you should consider using GPUs, if available. Additionally, you should consider parallelizing this task. For example, splitting each VCF file into multiple, smaller VCFs, and/or running separate jobs for Akita, Sei, and Enformer, instead of submitting run_models.sh (which sequentially scores varaints with each model in one job). The top of this script contains come parameters, including the path to output stdout/stderr, as well as the desired maximum memory, scratch, and CPU time alottment. These should be changed depending on the size of the job being submitted.



#### Next, find which genomic profiles from each model were more extreme than expected under the null hypothesis of random rare variants:
* run ./EDA/Akita_EDA/FindAkitaSigDNVs.py
* run ./EDA/Sei_EDA/FindSeiSigDNVs.py
* run ./EDA/Enformer_EDA/FindEnformerSigDNVs.py
```
top_level_dir=$(git rev-parse --show-toplevel) #path for top level of this repo
experiment_name=test #this must match the experiment name from above, to ensure variant scores from the correct path will be used

#### Find significant predictions from Akita ####
scoring_system=msd #for variants scored using mean squared allelic difference
python3 $top_level_dir/EDA/Akita_EDA/FindAkitaSigDNVs.py --scoring_system $scoring_system --experiment_name $experiment_name --git_root $top_level_dir

scoring_system=max #for variants scored using maximum absolute allelic difference
python3 $top_level_dir/EDA/Akita_EDA/FindAkitaSigDNVs.py --scoring_system $scoring_system --experiment_name $experiment_name --git_root $top_level_dir

#### Find significant predictions from Enformer ####
relevant_cols=HeartEnformerCols.txt #text file containing the Enformer tracks relevant to the Heart. Only these tracks will be tested for significance.

scoring_system=summed #for variants scored using summed allelic difference
python3 $top_level_dir/EDA/Enformer_EDA/FindEnformerSigDNVs.py --scoring_system $scoring_system --relevant_cols $relevant_cols --experiment_name $experiment_name --git_root $top_level_dir

scoring_system=max #for variants scored using max absolute allelic difference
python3 $top_level_dir/EDA/Enformer_EDA/FindEnformerSigDNVs.py --scoring_system $scoring_system --relevant_cols $relevant_cols --experiment_name $experiment_name --git_root $top_level_dir

#### Find significant predictions from Sei ####
python3 $top_level_dir/EDA/Sei_EDA/FindSeiSigDNVs.py --experiment_name $experiment_name --git_root $top_level_dir
```
This will yield two sets of files for each experiment and for each scoring system. One ({model_name}Sig_{experiment_name})includes the raw variant effect scores for only the significant genomic profiles for each variant. All other entries will be empty. The other file ({model_name}SigPVals_{experiment_name}) contains the p values for these significant genomic profiles. All other entires will be empty.

Alternatively, an example SGE script is provided so this may be submitted as a job. This is recommended: 
```
top_level_dir=$(git rev-parse --show-toplevel)
experiment_names=(test) #can include more than one experiment here if desired.


sh $top_level_dir/EDA/Akita_EDA/Example_SGE_run_AkitaSigVariantsJob.sh ${experiment_names[*]} #This will find significant variants from Akita using both MSD and Max scoring systems

relevant_cols=NeuroEnformerCols.txt
sh $top_level_dir/EDA/Enformer_EDA/Example_SGE_run_EnformerSigVariantsJob.sh $relevant_cols ${experiment_names[*]} #This will find significant variants from Enformer using both Summed and Max scoring systems

sh $top_level_dir/EDA/Sei_EDA/Example_SGE_run_SeiSigVariantsJob.sh  ${experiment_names[*]} #This will find Sei significant variants.

```

# Data visualizations of Akita and Enformer predictions for variant with extreme effects, and code to find variants with extreme predicted effects in all models, can be found in ./EDA/PCGC_EDA.ipynb


## Data Availability:
* Variant effect prediction scores for 58,090 CHD DNVs from Richter et al. (2020):
	* Enformer prediction scores can be found here (not hosted on GitHub because of storage limitations) https://ucsf.box.com/s/uopnbwg6rpul9ba0s1wu15zq4u4cpmmg
	* Akita prediction scores can be found in this repo under ./model_outputs/Akita/CHD_DNVs_Richter2020
	* Sei prediction scores can be found in this repo under ./model_outputs/Sei/CHD_DNVs_Richter2020
* Variant effect predictions for rare 1000 Genomes variants used for the null distributions:
	* Predictions for Akita, Sei, and Enformer can be found here: https://ucsf.box.com/s/uopnbwg6rpul9ba0s1wu15zq4u4cpmmg
* List of 456 DNVs whose predicted effects were significant in Akita, Sei, and Enformer can be found here: ./EDA/DNVs_Sig_in_ALL_models.csv
	* This csv contains chromosome, position (hg38), reference allele, and alternate allele for these 456 DNVs, as well as Akita and Sei predictions for each. This file also contains summary Enformer results (e.g., number of CAGE tracks significantly altered), rather than each individual predicted effect, because of storage considerations. Each variant appears in multiple rows to account for different scoring methods from each model. 
* Bootstrap results and visualization, as well as a list of genes within 100kb of prioritized CHD DNVs can be found in ./EDA/Enrichment_analysis.ipynb
