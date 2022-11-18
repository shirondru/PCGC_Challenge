# PCGC Challenge Scripts
 Scripts to run deep learning models that take DNA sequence as input and predict gene regulatory profiles, to prioritize causal variants and mechanisms

Models include Akita, Sei, an Enformer


Requirements:
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

(OPTIONAL). Download Enformer: Only necessary if running code from a machine without internet access. Otherwise, Enformer can be loaded from tfhub.
* Use either link:
	* https://tfhub.dev/deepmind/enformer/1
	* https://github.com/deepmind/deepmind-research/tree/master/enformer



Usage:
```
reference_genome=hg38 #either hg38 or hg19
experiment_name=test #your choice. This will help determine the names of output files
top_level_dir=$(git rev-parse --show-toplevel) #path for top level of this repo
akita_enformer_vcf_paths=$top_level_dir/DNVs/example.vcf #path to vcf file formatted for input into Akita/Enformer
sei_vcf_paths=$top_level_dir/DNVs/example_sei.vcf #path to vcf file formatted for input into Sei. Contains same variants as akita_enformer_vcf_paths

## run models
sh $top_level_dir/scripts/run_models.sh $reference_genome $akita_enformer_vcf_paths $sei_vcf_paths $experiment_name 


```






