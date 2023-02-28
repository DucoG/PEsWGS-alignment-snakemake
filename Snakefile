
import os
from pathlib import Path
from datetime import datetime

configfile: "config_snake.yaml"

hg_path_dict = {
    "hg19": "/home/d.gaillard/source/reference_genomes/hg19/hg19.fa",
    "hg38noalt": "/home/d.gaillard/source/reference_genomes/hg38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
    "hg38full": "/home/d.gaillard/source/reference_genomes/hg38_full/hg38.fa"
}

def generate_qc_outputs(step, refname):
    data_dir = Path("data")
    raw_files = list(data_dir.joinpath("raw").glob("*.fastq.gz"))
    qc_dir = Path("qc_outputs")

    step_to_suffix = {
        "raw": "_fastqc.html",
        f"aligned_{refname}": f"_{refname}_fastqc.html",
        f"marked_{refname}": f"_{refname}_mrk_fastqc.html",
    }

    output_qc_paths = []
    if step in ['raw', 'aligned', 'marked']:
        qc_suffix = step_to_suffix[f"{step}_{refname}"]
        if step in ["aligned", "marked"]:
            raw_files = [f for f in raw_files if "R1" in str(f)]
        for file in raw_files:
            new_file_name = f"{Path(file.stem).stem}{qc_suffix}".replace("_R1", "")
            new_file_path = qc_dir.joinpath(f"{step}_{refname}", "fastqc_output" , new_file_name)
            output_qc_paths.append(str(new_file_path))
    else:
        raise ValueError(f"Unknown step {step}")

    return output_qc_paths

def get_ext(step_folder):
    ext_map = {
        'raw': '.fastq.gz',
        'aligned': '.bam',
        'marked': '.bam'
    }
    return ext_map.get(step_folder)

def get_index_ext(step_folder):
    ext_map = {
        'raw': '.fastq.gz',
        'aligned': '.bam.bai',
        'marked': '.bam.bai'
    }
    return ext_map.get(step_folder)

def get_timestamp_string():
	return datetime.fromtimestamp(datetime.timestamp(datetime.now())).strftime("%y-%m-%d-%H:%M")


rule all:
    input:
        # generate_output_files('data/raw'),
        # "qc_outputs/raw/multiqc_output/multiqc_report.html",
        "qc_outputs/marked_hg19/multiqc_output/multiqc_report.html",
        # "qc_outputs/marked_hg38noalt/multiqc_output/multiqc_report.html"

rule fastq2bam:
    input: 
        fastq1 = 'data/raw/{id}_R1.fastq.gz',
        fastq2 = 'data/raw/{id}_R2.fastq.gz',
        ref_path = lambda wildcards: hg_path_dict[wildcards.refname]
    output: 
        bam = 'data/aligned_{refname}/{id}_{refname}.bam'
    threads: 4
    resources:
        mem_mb=10000
    conda:
        config['wgs_env']
    log:
        "logs/fastq2bam/{id}_{refname}.log"
    shell: 
        """
        (bwa mem -t {threads} {input.ref_path} {input.fastq1} {input.fastq2} | samtools sort > {output.bam}) 2> {log}
        """

rule index_bam:
    input:
        bam = 'data/{step_folder}/{sample}.bam'
    output:
        bai = 'data/{step_folder}/{sample}.bam.bai'
    conda:
        config['wgs_env']
    log:
        "logs/index_bam/{step_folder}_{sample}.log"
    threads: 2
    shell:
        """
        (samtools index {input.bam}) 2> {log}
        """

rule mark_duplicates:
    input:
        in_bam = 'data/aligned_{refname}/{id}_{refname}.bam',
        in_bai = 'data/aligned_{refname}/{id}_{refname}.bam.bai'
    output:
        out_bam = 'data/marked_{refname}/{id}_{refname}_mrk.bam',
        metrics_file = 'qc_outputs/marked_{refname}/{id}_{refname}_mrk_stats.txt'
    conda:
        config['wgs_env']
    log:
        "logs/mark_duplicates/{id}_{refname}.log"
    threads: 2
    resources:
        mem_mb=1024
    shell:
        """
        (picard MarkDuplicates -Xmx{resources.mem_mb}m INPUT={input.in_bam} OUTPUT={output.out_bam} M={output.metrics_file}) 2> {log}
        """


# rule fastqc
rule fastqc:
    input:
        file = lambda wildcards: f'data/{wildcards.step_folder}_{wildcards.refname}/{wildcards.id}' + get_ext(wildcards.step_folder),
        index_file = lambda wildcards: f'data/{wildcards.step_folder}_{wildcards.refname}/{wildcards.id}' + get_index_ext(wildcards.step_folder)
    output:
        fastqc_report = 'qc_outputs/{step_folder}_{refname}/fastqc_output/{id}_fastqc.html'
    conda:
        config['wgs_env']
    threads: 4
    resources:
        mem_mb=2024
    shell:
        """
        fastqc -t {threads} -o qc_outputs/{wildcards.step_folder}_{wildcards.refname}/fastqc_output {input.file}
        """

rule multiqc:
    input:
        infiles = lambda wildcards: generate_qc_outputs(wildcards.step_folder, wildcards.refname)
    output:
        multiqc_report = "qc_outputs/{step_folder}_{refname}/multiqc_output/multiqc_report.html"
    conda:
        config['wgs_env']
    threads: 2
    shell:
        """
        multiqc --outdir qc_outputs/{wildcards.step_folder}_{wildcards.refname}/multiqc_output qc_outputs/{wildcards.step_folder}_{wildcards.refname}
        """ 