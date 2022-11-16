# PCGC Challenge Scripts
 Scripts to run deep learning models that take DNA sequence as input and predict gene regulatory profiles, to prioritize causal variants and mechanisms

Models include Akita, Sei, an Enformer


Usage:
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







