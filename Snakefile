
import os
from pathlib import Path
from datetime import datetime

configfile: "config_snake.yaml"

hg_path_dict = config['hg_path_dict']

IDS, = glob_wildcards("data/raw/{id}.fastq.gz")
IDS = list(set(['_'.join(x.split('_')[:-1]) for x in IDS]))

fastq_files, = glob_wildcards("data/raw/{fastq_file}.fastq.gz")

def get_ext(step_folder):
    ext_map = {
        'raw': '.fastq.gz',
        'aligned': '.bam',
        'marked': '.bam'
    }
    return ext_map.get(step_folder)


rule all:
    input:
        # generate_output_files('data/raw'),
        "qc_outputs/raw/multiqc_output/multiqc_report.html",
        "qc_outputs/marked_hg19/multiqc_output/multiqc_report.html",
        # "qc_outputs/marked_hg38noalt/multiqc_output/multiqc_report.html"

rule bwa_mem2_index:
    input:
        "{genome}",
    output:
        "{genome}.0123",
        "{genome}.amb",
        "{genome}.ann",
        "{genome}.bwt.2bit.64",
        "{genome}.pac",
    log:
        "logs/bwa-mem2_index{genome}.log",
    wrapper:
        "v1.23.4/bio/bwa-mem2/index"

ref_path_list = lambda refname: hg_path_dict[refname]
rule bwa_mem2_mem:
    input:
        reads=["data/raw/{id}_R1.fastq.gz", "data/raw/{id}_R2.fastq.gz"],
        # idx=lambda wildcards: ref_path_list(wildcards.refname)
        # Index can be a list of (all) files created by bwa, or one of them
        idx=multiext("/home/d.gaillard/source/PEsWGS-alignment-snakemake/ref_genome/hg19.fa", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
    output:
        bam = 'data/aligned_{refname}/{id}_{refname}.bam',
    log:
        "logs/bwa_mem2/{id}_{refname}.log",
    params:
        extra="",
        sort="samtools",  # Can be 'none', 'samtools' or 'picard'.
        sort_order="coordinate",  # Can be 'coordinate' (default) or 'queryname'.
        sort_extra="",  # Extra args for samtools/picard.
    threads: 8
    wrapper:
        "v1.23.4/bio/bwa-mem2/mem"

rule samtools_index:
    input:
        bam = 'data/{step_folder}/{sample}.bam'
    output:
        bai = 'data/{step_folder}/{sample}.bam.bai'
    log:
        "logs/samtools_index/{step_folder}_{sample}.log",
    params:
        extra="",  # optional params string
    threads: 2  # This value - 1 will be sent to -@
    wrapper:
        "v1.23.4/bio/samtools/index"

rule mark_duplicates:
    input:
        bams='data/aligned_{refname}/{id}_{refname}.bam',
    # optional to specify a list of BAMs; this has the same effect
    # of marking duplicates on separate read groups for a sample
    # and then merging
    output:
        bam="data/marked_{refname}/{id}_{refname}_mrk.bam",
        metrics="qc_outputs/marked_{refname}/{id}_{refname}.metrics.txt",
    log:
        "logs/picard/marked/{id}_{refname}.log",
    params:
        extra="--REMOVE_DUPLICATES false",
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=1024,
    wrapper:
        "v1.23.4/bio/picard/markduplicates"

rule fastqc_raw:
    input: "data/raw/{sample}.fastq.gz"
    output:
        html='qc_outputs/raw/fastqc_output/{sample}_fastqc.html',
        zip='qc_outputs/raw/fastqc_output/{sample}_fastqc.zip' # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: "--quiet"
    log:
        "logs/fastqc_raw/{sample}.log"
    threads: 1
    wrapper:
        "v1.23.4/bio/fastqc"
        
rule fastqc_marked:
    input:
        lambda wildcards: f'data/{wildcards.step_folder}_{wildcards.refname}/{wildcards.id}' + get_ext(wildcards.step_folder)
    output:
        html='qc_outputs/{step_folder}_{refname}/fastqc_output/{id}_fastqc.html',
        zip='qc_outputs/{step_folder}_{refname}/fastqc_output/{id}_fastqc.zip' # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: "--quiet"
    log:
        "logs/fastqc_{step_folder}/{id}_{refname}.log"
    threads: 1
    wrapper:
        "v1.23.4/bio/fastqc"

def generate_neccesary_fastqcs_raw(fastq_files):
    base_path = Path('qc_outputs').joinpath("raw").joinpath('fastqc_output')
    filenames = [fastq_file + "_fastqc.zip" for fastq_file in fastq_files]
    full_paths = [ str(base_path.joinpath(filename)) for filename in filenames]
    return full_paths

# rule multiqc_raw:
#     input:
#         infiles = generate_neccesary_fastqcs_raw(fastq_files)
#     output:
#         multiqc_report = "qc_outputs/raw/multiqc_output/multiqc_report.html"
#     conda:
#         config['wgs_env']
#     shell:
#         """
#         multiqc --outdir qc_outputs/raw/multiqc_output qc_outputs/raw
#         """ 

rule multiqc_dir_raw:
    input:
        generate_neccesary_fastqcs_raw(fastq_files)
    output:
        "qc_outputs/raw/multiqc_output/multiqc_report.html"
    params:
        extra=""  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc_raw.log"
    wrapper:
        "v1.23.4/bio/multiqc"


def generate_neccesary_fastqcs_marked(IDS, refname):
    base_path = Path('qc_outputs').joinpath(f"marked_{refname}").joinpath('fastqc_output')
    filenames = [ID + f"_{refname}_mrk_fastqc.zip" for ID in IDS]
    full_paths = [ str(base_path.joinpath(filename)) for filename in filenames]
    return full_paths

# rule multiqc_marked:
#     input:
#         infiles = lambda wildcards: generate_neccesary_fastqcs_marked(IDS, wildcards.refname)
#     output:
#         multiqc_report = "qc_outputs/marked_{refname}/multiqc_output/multiqc_report.html"
#     conda:
#         config['wgs_env']
#     shell:
#         """
#         multiqc --outdir qc_outputs/marked_{wildcards.refname}/multiqc_output qc_outputs/marked_{wildcards.refname}
#         """ 

rule multiqc_dir_marked:
    input:
        lambda wildcards: generate_neccesary_fastqcs_marked(IDS, wildcards.refname)
    output:
        "qc_outputs/marked_{refname}/multiqc_output/multiqc_report.html"
    params:
        extra=""  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc_{refname}.log"
    wrapper:
        "v1.23.4/bio/multiqc"