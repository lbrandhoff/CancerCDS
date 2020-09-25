# CancerCDS

This clinical decision support system is used to prioritise somatic mutations in VCF format by their cancer driverness.\
The tool is designed to be used for the GrCH37 genome.

## System requirements
The tool is written in Python 3 (3.6)
It requires a local installation of the Variant Effect Predictor Tool (release 98.3)

### Setup requirements for the VEP tool
* Local installation + cache (Homo sapiens, release 98, GrCH37)
* dbNSFP plugin (>=v3.5)
* gnomAD Genome variants, tabix indexed
* GrCH37 fasta, tabix indexed

### The following packages are required
* numpy (>=v1.17.4)
* pandas (>=v0.25.3)
* pyfaidx (>=v0.5.9.1)


## Run
First, download the classifiers (link) and to put them into the `./CancerCDS/classifiers` directory.\
Then, run the run_project.py file, required inputs are the following:
* *-i*    input VCf file
* *-o*    output directory
* *-c*    cancer t
acronym (optional), a list of available cancer types can be found in `supported_cTypes.txt`
* *-p* minimum decision probability of each variant in output of being cancer driving (min 0.5)
* *-r* path to GrCH37 reference genome
* *-vt* path to VEP executable
* *-vc* path to VEP `cache` directory
* *-gn* path to gnomAD genome file
* *-db* path to dbNSFP genome file

## Output

If a viable input cancer type is given, the tool will analyse the input file with the respective type-specific and pan-cancer classifier.
If no cancer type is given, onyl the pan-cancer classification is performed.\
Output files can be found in the supplied output directory. Each file contains different information for each sufficient driver mutation, namely HGVSg notation, affected gene symbol, decision probability of being cancer driving and treatment recommendations.
