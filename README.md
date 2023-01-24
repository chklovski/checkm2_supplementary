This repository contains all the workflows necessary to update CheckM2 and re-train the models with new high-quality genomes. 

The workflow is provided for reference and transparency, but updating CheckM2 and its models will be done centrally by the CheckM2 team and new releases will be generated on the CheckM2 github (https://github.com/chklovski/CheckM2). 

To update CheckM2 models, you will need the CheckM2 environment as well as have bbmap and zenodo_backpack installed (both available from conda). 

To run this workflow, you will need to provide your new, high-quality genomes in .faa (protein) format. Please make sure to manually curate them to ensure only high-quality, complete, uncontaminated genomes are provided when retraining CheckM2 models. 

This workflow consists of downloading synthetic genomes generated in the previous CheckM2 release version from zenodo via zenodobackpack, generating new synthetic genomes from provided complete genomes, merging old and new input vectors in a sparse format, and then retraining three new models on the new merged input vectors - a completeness prediction model using gradient boost, a contamination prediction model using gradient boost, and a completeness prediction model using neural networks. 

To run the workflow, git clone the repository, then place all new genomes in the 'new_genomes' folder. 

After that, run `python retrainer.py <threads> <release_version> <old_DOI>` (positional arguments)

Where `<release_version>` is tied to the GTDB version the new complete genomes are from, e.g. r202, r207 etc. The old DOI is the DOI containing the synthetic genomes from Zenodo (complete list follows). 

Current DOI's are:

|DOI | GTDB Release Version | Upload Date |
|  :---:   |  :---:   |  :---:  |
| 10.5281/zenodo.7563512 | r202 | 24/01/2023 |

More threads will make training faster - please keep in mind the training can consume a substantial (> 250GB) amount of RAM and some time, especially if not using GPU for neural network training. 

The script will automatically go through the workflow and generate several outputs: 

a `new_models` folder containing the newly trained models for incorporation into a CheckM2 update, as well as a new min_ref_rsdata_<version>.npz incorporating the new genomes for switching between specific/general model when predicting completeness. These can be added to a new CheckM2 release. 

a `<new_version>_database.zb.tar.gz` containing combined old and new synthetic genome vectors to uploaded to Zenodo. When uploading to Zenodo, you MUST specify the zenodo version to be exactly the same as the version you provided to retrainer.py.

For general information, the full_feature_vector_list.tsv consists of all input vectors that are used in CheckM2 as well as their arrangement. Neural networks only use the first 20021 of these vectors, while the gradient boost models use 21241. 
