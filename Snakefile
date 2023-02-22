
import os
from pathlib import Path
from datetime import datetime

configfile: "config_snake.yaml"

def generate_output_files(raw_folder):
    refname = config['reference_name']
    out_folder = raw_folder.replace('raw', f'marked_{refname}')
    output_files = []
    for file in os.listdir(raw_folder):
        if file.endswith("R1.fastq.gz"):
            output_files.append(
                os.path.join(out_folder,file.replace('_R1.fastq.gz', f'_{refname}_mrk.bam.bai')))
    return output_files

def generate_qc_outputs(step):
    refname = config['reference_name']
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

def get_timestamp_string():
	return datetime.fromtimestamp(datetime.timestamp(datetime.now())).strftime("%y-%m-%d-%H:%M")


rule all:
    input:
        generate_output_files('data/raw'),
        "qc_outputs/raw/multiqc_output/multiqc_report.html",
        f"qc_outputs/marked_{config['reference_name']}/multiqc_output/multiqc_report.html"

rule fastq2bam:
    params:
        refname = config['reference_genome']
    input: 
        fastq1 = 'data/raw/{id}_R1.fastq.gz',
        fastq2 = 'data/raw/{id}_R2.fastq.gz'
    output: 
        bam = 'data/aligned_{refname}/{id}_{refname}.bam'
    threads: 4
    conda:
        config['wgs_env']
    log:
        "logs/fastq2bam/{id}_{refname}.log"
    shell: 
        """
        (bwa mem -t {threads} {params.refname} {input.fastq1} {input.fastq2} | samtools sort > {output.bam}) 2> {log}
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
        in_bam = 'data/aligned_{refname}/{id}_{refname}.bam',
        in_bai = 'data/aligned_{refname}/{id}_{refname}.bam.bai'
    output:
        out_bam = 'data/marked_{refname}/{id}_{refname}_mrk.bam',
        metrics_file = 'qc_outputs/marked_{refname}/{id}_{refname}_mrk_stats.txt'
    conda:
        config['wgs_env']
    log:
        "logs/mark_duplicates/{id}_{refname}.log"
    shell:
        """
        (picard MarkDuplicates INPUT={input.in_bam} OUTPUT={output.out_bam} M={output.metrics_file}) 2> {log}
        """


# rule fastqc
rule fastqc:
    params:
        refname = config['reference_name']
    input:
        file = lambda wildcards: f'data/{wildcards.step_folder}_{wildcards.refname}/{wildcards.id}' + get_ext(wildcards.step_folder)
    output:
        fastqc_report = 'qc_outputs/{step_folder}_{refname}/fastqc_output/{id}_fastqc.html'
    conda:
        config['wgs_env']
    shell:
        """
        fastqc -o qc_outputs/{wildcards.step_folder}/fastqc_output {input.file}
        """

rule multiqc:
    input:
        infiles = lambda wildcards: generate_qc_outputs(wildcards.step_folder),
        refname = config['reference_genome']
    output:
        multiqc_report = "qc_outputs/{step_folder}_{refname}/multiqc_output/multiqc_report.html"
    conda:
        config['wgs_env']
    shell:
        """
        multiqc --outdir qc_outputs/{wildcards.step_folder}_{wildcards.refname}/multiqc_output qc_outputs/{wildcards.step_folder}_{wildcards.refname}
        """ 