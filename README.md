# PCGC Challenge Scripts
 Scripts to run deep learning models that take DNA sequence as input and predict gene regulatory profiles, to prioritize causal variants and mechanisms

Models include Akita, Sei, an Enformer


Usage:
1. Clone repository containing Akita:
	* https://github.com/calico/basenji
	* Place this repository under ./models/Akita
	``` 
	cd ./models/Akita/basenji-master/manuscripts/akita
	sh get_model.sh
	sh get_data.sh
	```
2. Clone repository containing Sei:
	* https://github.com/FunctionLab/sei-framework
	* Place this repository under ./models/Sei

(OPTIONAL). Download Enformer: Only necessary if running code from a machine without internet access. Otherwise, Enformer can be loaded from tfhub.
* Use either link:
	* https://tfhub.dev/deepmind/enformer/1
	* https://github.com/deepmind/deepmind-research/tree/master/enformer







