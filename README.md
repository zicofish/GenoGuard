This repository contains all the code that was written while we produce the following two papers:

[1] Z. Huang, E. Ayday, J. Fellay, J.-P. Hubaux and A. Juels, "GenoGuard: Protecting Genomic Data against Brute-Force Attacks," *36th IEEE Symposium on Security and Privacy (S&P 2015)*, San Jose, CA, USA, 2015.

[2] S. Samani, Z. Huang, E. Ayday, M. Elliot, J. Fellay, J.-P. Hubaux, and Z. Kutalik, "Quantifying Genomic Privacy via Inference Attack with High-Order SNV Correlations," in *2nd International Workshop on Genome Privacy and Security (in conjunction with IEEE S&P; 2015)*, 2015.

# Dataset
This repository contains a small dataset (part of chromosome 22) in ./genoguard/hapmap/chr22 for you to run the code successfully.
The data is produced from the [Hapmap dataset](http://hapmap.ncbi.nlm.nih.gov/downloads/index.html.en). We apply handy pre-processing steps (in the two files `dat.dataPreprocess` and `dat.dataConverter`) to transform the original Hapmap data into more simplified formats, so that it is convenient for us to use.
More explanation about the format can be found in the code.

# GenoGuard
To test GenoGuard on encrypting and decrypting a genome sequence, simply run the python code `honey_encryption/GenoGuard.py`, which is the implementation of honey encryption with the proposed DTE models on genomic data:

`python GenoGuard.py`

# Inference of hidden SNVs
The pipeline for our experiment of inference attack is described as follows. 

To produce data with hidden SNVs, run the function `hideSNVs` in the code `dat/dataConverter.py` :

`python dataConverter.py`

To predict the hidden SNVs with the recombination model, run the code `inf/inference.py`:

`python inference.py`

Of course, you can predict the hidden SNVs with other models, which is a simple extension I have not implemented in the current version. But you can find instructions on how to do this in `inference.py`. 

# Contact
Zhicong Huang (zhicong.huang@epfl.ch)