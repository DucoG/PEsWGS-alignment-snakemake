
import os
from pathlib import Path
from datetime import datetime



configfile: "config_snake.yaml"

def generate_output_files(raw_folder):
    out_folder = raw_folder.replace('raw', 'marked')
    output_files = []
    for file in os.listdir(raw_folder):
        if file.endswith("R1.fastq.gz"):
            output_files.append(
                os.path.join(out_folder,file.replace('_R1.fastq.gz', '_hg38_mrk.bam.bai')))
    return output_files

def generate_qc_outputs(step):
    data_dir = Path("data")
    raw_files = list(data_dir.joinpath("raw").glob("*.fastq.gz"))
    qc_dir = Path("qc_outputs")

    step_to_suffix = {
        "raw": "_fastqc.html",
        "aligned": "_hg38_fastqc.html",
        "marked": "_hg38_mrk_fastqc.html",
    }

    output_qc_paths = []
    if step in step_to_suffix:
        qc_suffix = step_to_suffix[step]
        if step in ["aligned", "marked"]:
            raw_files = [f for f in raw_files if "R1" in str(f)]
        for file in raw_files:
            new_file_name = f"{Path(file.stem).stem}{qc_suffix}".replace("_R1", "")
            new_file_path = qc_dir.joinpath(step, "fastqc_output" , new_file_name)
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

def get_timestamp_string():
	return datetime.fromtimestamp(datetime.timestamp(datetime.now())).strftime("%y-%m-%d-%H:%M")


rule all:
    input:
        generate_output_files('data/raw'),
        "qc_outputs/raw/multiqc_output/multiqc_report.html",
        "qc_outputs/marked/multiqc_output/multiqc_report.html"

# rule from fastq to bam
rule fastq2bam:
    input: 
        fastq1 = 'data/raw/{id}_R1.fastq.gz',
        fastq2 = 'data/raw/{id}_R2.fastq.gz',
        reference = config['reference_genome']
    output: 
        bam = 'data/aligned/{id}_' + config['reference_name'] + '.bam'
    threads: 4
    conda:
        config['wgs_env']
    log:
        "logs/fastq2bam/{id}.log"
    shell: 
        """
        (bwa mem -t {threads} {input.reference} {input.fastq1} {input.fastq2} | samtools sort > {output.bam}) 2> {log}
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
    shell:
        """
        (samtools index {input.bam}) 2> {log}
        """

rule mark_duplicates:
    input:
        in_bam = 'data/aligned/{sample}.bam',
        in_bai = 'data/aligned/{sample}.bam.bai'
    output:
        out_bam = 'data/marked/{sample}_mrk.bam',
        metrics_file = 'qc_outputs/marked/{sample}_mrk_stats.txt'
    conda:
        config['wgs_env']
    log:
        "logs/mark_duplicates/{sample}.log"
    shell:
        """
        (picard MarkDuplicates INPUT={input.in_bam} OUTPUT={output.out_bam} M={output.metrics_file}) 2> {log}
        """


# rule fastqc
rule fastqc:
    input:
        file = lambda wildcards: f'data/{wildcards.step_folder}/{wildcards.sample}' + get_ext(wildcards.step_folder)
    output:
        fastqc_report = 'qc_outputs/{step_folder}/fastqc_output/{sample}_fastqc.html'
    conda:
        config['wgs_env']
    shell:
        """
        fastqc -o qc_outputs/{wildcards.step_folder}/fastqc_output {input.file}
        """

rule multiqc:
    input:
        infiles = lambda wildcards: generate_qc_outputs(wildcards.step_folder)
    output:
        multiqc_report = "qc_outputs/{step_folder}/multiqc_output/multiqc_report.html"
    conda:
        config['wgs_env']
    shell:
        """
        multiqc --outdir qc_outputs/{wildcards.step_folder}/multiqc_output qc_outputs/{wildcards.step_folder}
        """