from pathlib import Path

configfile: "config_snake.yaml"

<<<<<<< HEAD
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
=======
IDS, = glob_wildcards("data/raw/{id}.fastq.gz")
IDS = list(set(['_'.join(x.split('_')[:-1]) for x in IDS]))
>>>>>>> origin/wrapper_version

# fastq_files, = glob_wildcards("data/raw/{fastq_file}.fastq.gz")

hg_path_dict = config['hg_path_dict']
def refname_2_bwamem2index(refname):
    return multiext(hg_path_dict[refname], ".amb", ".ann", ".bwt.2bit.64", ".pac")

rule all:
    input:
        # generate_output_files('data/raw'),
<<<<<<< HEAD
        # "qc_outputs/raw/multiqc_output/multiqc_report.html",
        "qc_outputs/marked_hg19/multiqc_output/multiqc_report.html",
        # "qc_outputs/marked_hg38noalt/multiqc_output/multiqc_report.html"
=======
        "qc_outputs/raw/multiqc_output/multiqc_report.html",
        # "qc_outputs/marked_hg19/multiqc_output/multiqc_report.html",
        "qc_outputs/marked_hg38noalt/multiqc_output/multiqc_report.html"
>>>>>>> origin/wrapper_version


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

rule bwa_mem2_mem:
    input:
        reads=["data/raw/{id}_R1.fastq.gz", "data/raw/{id}_R2.fastq.gz"],
        # Index can be a list of (all) files created by bwa, or one of them
        idx=lambda wildcards: refname_2_bwamem2index(wildcards.refname)
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

def fastqc_input(wildcards):
    return f"data/{wildcards.sample_type}/{wildcards.sample}.fastq.gz" if wildcards.sample_type == "raw" else f"data/{wildcards.sample_type}/{wildcards.sample}.bam"

rule fastqc:
    input:
        lambda wildcards: fastqc_input(wildcards)
    output:
        html='qc_outputs/{sample_type}/fastqc_output/{sample}_fastqc.html',
        zip='qc_outputs/{sample_type}/fastqc_output/{sample}_fastqc.zip'
    params: "--quiet"
    log:
        "logs/fastqc_{sample_type}/{sample}.log"
    threads: 1
    wrapper:
        "v1.23.4/bio/fastqc"

def generate_neccesary_fastqcs(sample_type):
    base_path = Path('qc_outputs').joinpath(sample_type).joinpath('fastqc_output')
    if sample_type == "raw":
        filenames = [ID + "_fastqc.zip" for ID in IDS]
    else:
        refname = sample_type.split('_')[-1]
        filenames = [ID + f"_{refname}_mrk_fastqc.zip" for ID in IDS]
    full_paths = [str(base_path.joinpath(filename)) for filename in filenames]
    return full_paths

rule multiqc_dir:
    input:
        lambda wildcards: generate_neccesary_fastqcs(wildcards.sample_type)
    output:
        "qc_outputs/{sample_type}/multiqc_output/multiqc_report.html"
    params:
        extra=""  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc_{sample_type}.log"
    wrapper:
        "v1.23.4/bio/multiqc"
