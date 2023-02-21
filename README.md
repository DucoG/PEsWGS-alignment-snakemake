# Paired-end Shallow Whole Genome Sequencing Snakemake Pipeline

This is a snakemake pipeline for the analysis of paired-end shallow whole genome sequencing data. The pipeline is designed to work with raw fastq files generated from Illumina HiSeq or NovaSeq. It includes the following steps:

- Align reads to a reference genome using BWA
- Index the resulting BAM file
- Mark duplicates using Picard
- Index the marked BAM file
- Generate fastqc report for raw, aligned, and marked reads
- Generate multiqc report for all three reports

## Setup

To use this pipeline, first clone the repository and navigate to the project directory:

```bash
git clone https://github.com/DucoG/PEsWGS-alignment-snakemake.git
cd PEsWGS-alignment-snakemake
```

### Conda Environment

The pipeline requires a conda environment which can be installed using the `env.yml` file in the repository. To install the environment, run the following command:
```bash
conda env create --file env.yml
```

This will create a conda environment named `wgs` which is required for the pipeline.

### Config File

The `config_snake.yaml` file contains variables and settings required for the pipeline. Please update this file with your own values before running the pipeline.

## Usage
Once the environment and config file have been set up, the pipeline can be run using snakemake. Before running the entire pipeline for all your samples, it is advisable to first check if all settings are correct using a dry run using the `-n` option for a dry-run and `-p` for printing the ran commands:
```bash
snakemake -np --use-conda
```

To run the entire pipeline, use the following command:
```bash
snakemake --use-conda -j <NUMBER_OF_JOBS>
```
This will generate marked BAM files and QC reports for the data in `data/raw`. Your directory will look as follows:
```bash
.
├── data                        <- contains all input and output files for the pipeline
│   ├── aligned
│   ├── marked
│   ├── raw
├── logs                        <- contains log files for each step of the pipeline
│   ├── index_bam
│   └── mark_duplicates
└── qc_outputs                  <- contains output files from quality control steps
    ├── marked
    │   ├── fastqc_output
    │   └── multiqc_output
    │       └── multiqc_data
    └── raw
        ├── fastqc_output
        └── multiqc_output
            └── multiqc_data
```


To run a specific rule, for example, to run only the `fastq2bam` rule, use the following command:

```bash
snakemake data/aligned/<SAMPLE_ID>_hg38.bam --use-conda
```
This will generate a BAM file for the specified sample in `data/aligned`.

## Note

This pipeline does not trim the sequences automatically. To determine if trimming is necessary, start by generating the raw multiqc report. Based on that, trimming can be done on the raw data files before aligning.

