
# GRACE
GRACE (**G**astric cancer detection with f**RA**gmentomi**C** and **E**pigenetic features) is a multidimensional diagnostic model for gastic cancer. The fast workflow of GRACE can realize accurate detection of gastric cancer within 6.5 hours of blood draw.  

The GRACE fast diagnostic pipeline integrates the data filtration process[1], the extraction process of fragmentomic (fragment ratio, end motif[1]) and epigenetic (genome-level, bin-level) features, and the GRACE model prediction process. The output of this pipeline is the features extracted and the G-score predicted by GRACE.

## Requirements
#### Hardware requirements
- MinION or GridION device (Oxford Nanopore Technology, UK)
- R10 flow cell (FLO-MIN114 R.10.4.1)
- Server with GPU for basecalling
#### Software requirements
- dorado (0.8.0, 1.0.2)  
- samtools (1.13)  
- perl (5.34.0)  
- awk (5.1.0)  
- wc (8.32)  
- samtools (1.13)  
- modkit (0.4.0)  
- bedtools (2.30.0)  
- R (4.3.0) packages: vroom, data.table, tidyr, readr, tidyverse, dplyr 
- Python (3.13) dependencies: sys, joblib, shutil, subprocess, numpy, pandas  

**Note that:**   
Both 0.8.0 and 1.0.2 version of dorado are needed for analysis, please use different suffix to distinguish between two versions (e.g. dorado0.8.0, dorado1.0.2).

## Installation
Clone this project to your home (typically takes 10s to install) :
```
git clone https://github.com/DNAnoPore-lab/GRACE.git
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
Please download FASTA of reference GRChg38 at https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz before start, and make sure it is placed in the same directory as `GRACE_pipeline.py` under the name of `Homo_sapiens_assembly38.fasta`.  
Generate index for reference:
```
samtools faidx Homo_sapiens_assembly38.fasta
```
Align with dorado (0.8.0):
```
dorado0.8.0 aligner Homo_sapiens_assembly38.fasta [path_to_BAM_of_correct_barcode] > aligned.bam
```
### GRACE fast diagnostic pipeline
Please make sure `aligned.bam` is placed in the same directory as `GRACE_pipeline.py`.  
Please do not change any name or structure of the directory in this project.  
Run the following command line:
```
python GRACE_pipeline.py
```
The expected run time for demo is 2 minutes.  
The expected output contains the features extracted and the G-score:
```
short_fragment_ratio: [[0.18833348 0.59683794 0.16401111]]
end_motif_feature: [['0.83321152' '0.33620816' '0.46776787' '0.38006139' '0.38737027'
  '0.40929689' '0.46776787' '0.46776787' '0.06577986' '0.29235492']]
gene_level_feature: [[79.91704943]]
bin_level_feature: [['79.56989' '76' '81.39535' '84.25197' '88' '72.22222' '82.25806'
  '68.14159' '81.70732' '83.60656' '79.01235' '80.70175' '88.73874'
  '84.74576' '83.63636' '83' '85.9375' '70.74468' '87.95181' '73.33333'
  '80' '80.89888' '71.42857' '76.31579' '85.05747' '82.8125' '73.52941'
  '79.13043' '80.43478' '76.92308' '86.53846' '69.76744' '76' '91.66667'
  '77.9661' '68.35443' '91.34615' '92.85714' '83.33333' '81.01266'
  '83.78378' '84.93151' '90' '81.81818' '76.47059' '74.02597' '74.19355'
  '74.5098' '76.66667' '71.26437' nan '62.96296' '75.4386' '75' '83.75'
  '92.53731' '70.4918' '88.50575' '88.46154' '90.625' '80' '72'
  '80.50847' '83.92857' '66.92308' '78.23529' '73.17073' '60.37736'
  '83.15789' '79.26829' '78.31325' '73.03371' '76' '71.27197' '76.92308'
  '85.96491' '71.875' '76.21359' '73.33333' '83.82353' '82.78689'
  '77.02703' '75.24752' '81.41026' '89.41176' '88.28125' '79.16667'
  '86.66667' '77.10843' '86.95652' '75' '61.53846' '82.92683' '89.24051'
  '87.2093' '76.42276' '84.88372' '73.40426' '73.35766' '72.25806'
  '85.49618' '71.34831' '74.44444' '95' '100' '92' '74.54545' '90'
  '95.65217' '87.83784' '68.96552' '84.76562' '89.69072' '71.79487'
  '94.3662' '83.09859' '85.36585' '75.81967' '86.30137' '92.1875'
  '76.85185' '78.57143' '79.48718' '85.25641' '70.4918' '80.82192'
  '53.33333' '67.34694' '81.5534' '87.58865' '68.18182' '67.30769'
  '63.75' '82.8125' '82.07547' '82.85714' '77.13178' '71.31783'
  '78.18182' '76.74419' '92.85714' '63.49206' '85' '91.42857' '89.47368'
  '78.84615' '65.51724' '77.02703' '75' '76.47059' '72.72727' '75'
  '87.09677' '85' '100' '88.99083' '71.69811' '83.60656' '76.74419'
  '85.96491' '78.88889' '97.36842' '84.90566' '86.66667' '88.67925'
  '96.9697' '76.5625' '78.94737' '86.84211' '62.5' '75.45455' '63.51351'
  '88.88889' '70.76923' '76.74419' '80.70175' '84.74576' '72.88136'
  '74.35897' '75' '73.33333' '82.14286' '73.40426' '70.65217' '66.66667'
  '73.95833' '68.68132' '76.19048' '81.25' '87.01299' '82.23684'
  '86.66667' '83.33333' '74.46809' '68.69919' '79.5082' '85.2459' '75'
  '74.24242' '84.66258' '91.48936' '85.29412' '80' '78.37838' '66.66667'
  '76.74419' '78.125' '85.45455' '63.63636' '91.13924' '90.84507'
  '72.72727' '78.40909' '75' '74.4186' '94.59459' '83.63636' '78.57143'
  '54.31034' '85.71429' '65.90909' '65.87302' '83.33333' '78.66667'
  '82.08955' '76.36364' '81.48148' '71.73913' '75' '89.28571' '84.84848'
  '89.62264' '90.27778' '87.23404' '81.57895' '75' '84.81013' '90.66667'
  '75.5102' '75.60976' '59.18367' '89.55224' '68.51852' '81.08108'
  '88.17204' '89.47368' '87.7551' '74.28571' '84.14634' '86.25'
  '87.27273' '70.21277' '83.33333' '81.7757' '82.22222' '84.84848'
  '88.88889' '84.61538' '82.25806' '84.31373' '81.73913' '73.68421'
  '68.96552' '82.17822' '89.3617' '86.66667' '90.21739' '79.80769'
  '80.45977' '80' '79.6875' '80.14706' '75.75758' '74.4186' '85.60606'
  '72.91667' '81.69014' '75' '83.33333' '77.77778' '89.79592' '67.92453'
  '85' '77.41935' '81.66667' '66.66667' '82.40741' '63.7931' '81.69014'
  '86.04651' '82.92683' '87.32394' '72.72727' '32.95455' '64.40678'
  '73.8806' '85.41667']]
  
  G-score: [0.65013023]
```



[1] The data filtration process and the end motif feature extraction process in this pipeline were based on codes deposited by Katsman et al ([Puputnik/Fragmentomics_GenomBiol: Pipeline to replicate Katsman et. al fragmentomics results](https://github.com/Puputnik/Fragmentomics_GenomBiol)).
