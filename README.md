![Snakemake](https://img.shields.io/badge/snakemake- >=4.8.0-brightgreen.svg?style=flat-square)
# khanlab_pipeline
This is the implementation of [KhanLab](https://ccr.cancer.gov/Genetics-Branch/javed-khan) NGS Pipeline using Snakemake.

## Installation

The easiest way to get this pipeline is to clone the repository.

```
git clone git@github.com:hsienchao/khanlab_pipeline.git
```
This pipeline is available on NIH biowulf cluster, contact me if you would like to do a test run. The data from this pipeline could directly be ported in [OncoGenomics-DB](https://clinomics.ncifcrf.gov/production/public/), an application created to visualize NGS data available to NIH users.
## Requirements
[snakemake 5.13.0](https://snakemake.readthedocs.io/en/stable/)  
[mutt](http://www.mutt.org/)  
[gnu parallel](http://www.gnu.org/software/parallel/)  
SLURM or PBS for resource management  
Bioinformatics Tools Listed in [config files](config/config_common.json)  

### HiC:
- [HiCPro] (https://github.com/nservant/HiC-Pro)
- [JuiceBox] (https://github.com/aidenlab/Juicebox)
- [Bowtie2] (https://github.com/BenLangmead/bowtie2)


## Conventions

- Sample names cannot have "/" or "." in them
- Fastq files end in ".fastq.gz"

## Run pipeline

### HiC

#### Input 

- Sample sheet in YAML format
- Sample sheet can be generated using script sampleToYaml.py. Example:
```
python scripts/sampleToYaml.py /data/khanlab/projects/HiC/manage_samples/HiC_sample_sheet.xlsx RH4_Ent6_H3K27ac_HiChIP_HH3JVBGX7 > RH4_Ent6_H3K27ac_HiChIP_HH3JVBGX7.hic.yaml
```

##### Required columns: 
1. Amplified_Sample_Library_Name1. 
2. FlowCell_GSE.
3. Genome
4. SampleFiles

##### Optional columns
1. SpikeIn
2. SpikeInGenome

- Example:
```
samples:
  RH4_D6_H3K27ac_HiChIP_HKJ22BGX7:
    Amplified_Sample_Library_Name: RH4_D6_H3K27ac_HiChIP
    FlowCell_GSE: HKJ22BGX7
    Genome: hg19
    SampleFiles: Sample_RH4_D6_H3K27ac_HiChIP_HKJ22BGX7
    SpikeIn: 'yes'
    SpikeInGenome: mm10
```
#### Output
1. Juicebox hic file: [output dir]/[sample_id].allValidPairs.hic
2. HiCpro pairs:
  * Pairs for the reference genome:
    * [output dir]/HiCproOUTPUT/hic_results/data/[sample_id]/[sample_id]/[sample_id].allValidPair
  * Mergestate for the reference genome:
    * [output dir]/HiCproOUTPUT/hic_results/stats/[sample_id]/[sample_id]/[sample_id].allValidPair.mergestat
  * Pairs for the spikeIn genome:
    * [output dir]/HiCproAQuAOUTPUT/hic_results/data/[sample_id]/[sample_id]/[sample_id].allValidPair
  * Mergestate for the spikeIn genome:
    * [output dir]/HiCproAQuAOUTPUT/hic_results/stats/[sample_id]/[sample_id]/[sample_id].allValidPair.mergestat

#### Launching the pipeline
##### 1. General launch script:

    launch [options]
	
    required options:

            -type|t         <string>        Pipeline type (available options: hic,chipseq,ranseq,dnaseq)
            -workdir|w      <string>        Working directory where all the results will be stored.
            -sheet|s        <string>        Sample sheet in YAML format

    optional options: 
            -datadir|d <string> FASTQ file location (default: /data/khanlab/projects/DATA) 
            -dryrun Dryrun only 
            -dag Generate DAG PDF

  * Example
```
      launch -type hic -w /data/khanlab/projects/HiC/processed_DATA/RH4_D6_H3K27ac_HiChIP_HKJ22BGX7 -s /data/khanlab/projects/HiC/processed_DATA/sample_sheets/RH4_D6_H3K27ac_HiChIP_HKJ22BGX7.hic.yaml
```
#####  2. Launch by sample ID (Khanlab automation):
```
    scripts/automate_hic.sh [sampleID]
```    
    This script will parse the Khanlab HiC sample sheet and check if FASTQ files are ready. Then it will lauch the HiC pipeline automatically. 
    
    The Khanlab data location:
    
    -- Sample sheet: /data/khanlab/projects/HiC/manage_samples/HiC_sample_sheet.xlsx
    -- FASTQ files: /data/khanlab/projects/DATA
    -- Processed data: /data/khanlab/projects/HiC/processed_DATA

  * Example 1: process specific sample
```
   scripts/automate_hic.sh RH4_Ent6_H3K27ac_HiChIP_HH3JVBGX7
```
  * Example 2: process all samples
```
   scripts/automate_hic.sh all
```

### ChIPseq
Not implemented yet.
### RNAseq
Not implemented yet.
### DNAseq
Not implemented yet.
### Methylseq
Not implemented yet.
## Rulegraph
DAG for HiC example
![alt tag](dag.hic.pdf)

For questions or comments, please contact: Hsienchao Chou