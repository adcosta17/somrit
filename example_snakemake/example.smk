#
# Data & Fastq related rules
#

def get_ref(wildcards):
    return config["reference"]

def get_base_dir(wildcards):
    return config["base_dir"]

def get_fastq(wildcards):
    return config["fastq"]

def get_repbase(wildcards):
    return config["repbase"]

def get_centromeres(wildcards):
    return config["centromeres"]

def get_telomeres(wildcards):
    return config["telomeres"]

def get_bam(wildcards):
    return config["base_dir"]+"/"+wildcards.sample+"/mapped/"+wildcards.sample+".bam"

def get_fastq(wildcards):
    return config["base_dir"]+"/"+wildcards.sample+"/fastq/"+wildcards.sample+".fastq.gz"

rule all:
    input:
        expand("{s}/Realigned_classified_filtered_{s}.tsv", s=config["samples"])

rule extract_inserts:
    output:
        tsv="{sample}/{sample}.tsv",
        merged=temp("{sample}/{sample}_merged.txt")
    threads: 10
    params:
        memory_per_thread="8G",
        base_dir=get_base_dir,
        bam=get_bam,
        fastq=get_fastq
    shell:
        """
        cd {config[somrit_dir]}
        python somrit.py extract --bam {params.bam} --output-merged {params.base_dir}/{output.merged} --output-tsv {params.base_dir}/{output.tsv} --fastq-file {params.fastq} --threads {threads}
        cd {params.base_dir}
        """

rule realign_inserts:
    input:
        tsv="{sample}/{sample}.tsv",
    output:
        bam="{sample}/Realigned_{sample}.bam",
        tsv="{sample}/Realigned_{sample}.tsv"
    threads: 10
    params:
        memory_per_thread="20G",
        base_dir=get_base_dir,
        ref=get_ref,
        bam=get_bam,
        fastq=get_fastq
    shell:
        """
        cd {config[somrit_dir]}
        python somrit.py realign --bam-list {params.bam} --tsv-list {params.base_dir}/{input.tsv} --fastq-list {params.fastq} --output-dir {params.base_dir} --tsv-prefix {wildcards.sample}/Realigned_{wildcards.sample} --bam-prefix {wildcards.sample}/Realigned_{wildcards.sample} --reference-genome {params.ref} --threads {threads} --filter-depth --max-insert-size 10000 --max-depth 500 
        cd {params.base_dir}
        """


rule classify_inserts:
    input:
        tsv="{sample}/{sample}.tsv",
        bam="{sample}/Realigned_{sample}.bam",
        bam_index="{sample}/Realigned_{sample}.bam.bai",
        realign_tsv="{sample}/Realigned_{sample}.tsv"
    output:
        tsv="{sample}/Realigned_classified_{sample}.tsv"
    threads: 1
    params:
        memory_per_thread="24G",
        base_dir=get_base_dir,
        fastq=get_fastq,
        repbase=get_repbase
    shell:
        """
        cd {config[somrit_dir]}
        python somrit.py classify --bam-list {params.base_dir}/{input.bam} --sample-list {wildcards.sample} --tsv-list {params.base_dir}/{input.tsv} --realign-tsv {params.base_dir}/{input.realign_tsv} --annotation-file {params.repbase} --fastq-list {params.fastq} --output-tsv {params.base_dir}/{output.tsv} 
        cd {params.base_dir}
        """

rule filter_inserts:
    input:
        tsv="{sample}/Realigned_classified_{sample}.tsv",
        bam="{sample}/Realigned_{sample}.bam",
        bam_index="{sample}/Realigned_{sample}.bam.bai"
    output:
        tsv="{sample}/Realigned_classified_filtered_{sample}.tsv"
    threads: 10
    params:
        memory_per_thread="20G",
        base_dir=get_base_dir,
        ref=get_ref,
        repbase=get_repbase,
        centromeres=get_centromeres,
        telomeres=get_telomeres,
    shell:
        """
        cd {config[somrit_dir]}
        python somrit.py filter --threads {threads} --input-tsv {params.base_dir}/{input.tsv} --bam {params.base_dir}/{input.bam} --fastq-list NA --reference-genome {params.ref} --centromeres {params.centromeres} --telomeres {params.telomeres} --output-tsv {params.base_dir}/{output.tsv}
        cd {params.base_dir}
        """


rule make_bam_index:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    params:
        memory_per_thread="16G"
    threads: 1
    shell:
        "samtools index {input}"



