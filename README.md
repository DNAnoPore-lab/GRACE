
# GRACE
GRACE (**G**astric cancer detection with f**RA**gmentomi**C** and **E**pigenetic features) is a multidimensional diagnostic model for gastic cancer. The fast workflow of GRACE that can realize accurate detection of gastric cancer within 6.5 hours of blood draw.

## Requirements
#### Hardware requirements
- MinION or GridION device (Oxford Nanopore Technology, UK)
- R10 flow cell (FLO-MIN114 R.10.4.1)
- Server with GPU for basecalling
#### Software requirements
- dorado (0.8.0, 1.0.2)  
- samtools (1.13)  
- awk (5.1.0)  
- wc (8.32)  
- samtools (1.13)  
- modkit (0.4.0)  
- bedtools (2.30.0)  
- R (4.3.0)  
- Python (3.13) dependencies: sys, joblib, shutil, subprocess, numpy, pandas  

**Note that:**   
Both 0.8.0 and 1.0.2 version of dorado are needed for analysis, please use different suffix to distinguish between two versions (e.g. dorado0.8.0, dorado1.0.2).

## Installation
Download the folder of this project to your home:
```
wget -O ~/GRACE_pipeline https://github.com/DNAnoPore-lab/GRACE.git
```

## Demo
A mapped BAM `aligned.bam` is available as an input demo for GRACE fast diagnositc pipeline.

## Instruction for use

### Basecalling
Please download dorado basecalling model `dna_r10.4.1_e8.2_400bps_sup@v5.0.0` at [nanoporetech/dorado: Oxford Nanopore's Basecaller](https://github.com/nanoporetech/dorado/) before start.  
Basecalling with dorado (0.8.0):
```
dorado0.8.0 basecaller [path_to_dorado_model]/dna_r10.4.1_e8.2_400bps_sup@v5.0.0 [path_to_pod5_directory]/pod5/ -r --modified-bases 5mCG_5hmCG --no-trim > demux.bam
```
### Demultiplex
Demultiplex with dorado (1.0.2):
```
dorado1.0.2 demux -t 200 --kit-name SQK-NBD114-24 --barcode-both-ends --output-dir demux demux.bam
```
Enter the directory `demux/` and check for BAM of correct barcode for alignment.
### Alignment
Please download FASTA of reference GRChg38 at https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz before start, and make sure they are placed in `Fragmentomics/Utility/` directory under the name of `Homo_sapiens_assembly38.fasta`.  
Generate index for reference:
```
samtools faidx Homo_sapiens_assembly38.fasta
```
Align with dorado (0.8.0):
```
dorado0.8.0 aligner ~/GRACE_pipeline/Fragmentomics/Utility/Homo_sapiens_assembly38.fasta [path_to_BAM_of_correct_barcode] > ~/GRACE_pipeline/aligned.bam
```
### GRACE fast diagnostic pipeline
Please make sure `aligned.bam` is placed in the same directory as `GRACE_pipeline.py`.  
Please do not change any name or structure of the directory in this project.  
Run the following command line:
```
python GRACE_pipeline.py
```

