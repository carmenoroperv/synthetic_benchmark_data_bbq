import sys
import os.path
import json

configfile: "config/lower_HiFi_SNVs_VAF_in_Illumina.yaml"
workdir: "../"

# input
reference_name = config["reference"]["name"]               # reference name
fai_file = config["reference"]["fai_file"]                 # reference fai file

sample = config["sample"]["name"]                          # sample name
bam_file = config["sample"]["bam_file"]                    # sample bam file

snv_calls_file = config["variants"]["snv_calls_file"]           # SNV calls for which the VAF will be changed
alt_keep_prob = config["variants"]["alt_keep_prob"]        # 1 - probability that the alternative base is changed to the reference base at the given position and read



# select only contigs that are longer than 10**6
contig_dict = {}
with open(fai_file) as contigs:
    for line in contigs:
        contig_dict[line.split()[0]] = line.split()[1]
contig_list = []
for contig_name, length in contig_dict.items():
    if int(length) > 10**6:
        contig_list.append(contig_name)

MODIFIED_BAM_FILES = expand(f"results/{sample}_{{contig}}.aligned.sorted.bam", contig = contig_list)
MODIFIED_BAM_MERGED = f"results/{sample}.aligned.sorted.bam"

rule all:
    input: MODIFIED_BAM_MERGED

rule change_VAF:
    input: 
        bam_file = bam_file,
        snvs_file = snv_calls_file
    output: 
        bam_file = temp(f"results/{sample}_{{contig}}.aligned.bam")
    resources: 
        mem_mb = 50000,
        time = "24:00:00" #
    params:
        log = f"logs/{reference_name}/{sample}/alt_keep_prob{alt_keep_prob}/change_VAF/{{contig}}.out",
        alt_keep_prob = f"{alt_keep_prob}",
        contig = "{contig}"
    conda: "../../envs/change_illumina_VAF.yaml"
    script: 
        "../scripts/change_VAF_in_Illumina.py"

rule sort_modified:
    input: 
        bam_file = rules.change_VAF.output.bam_file,
    output: 
        sorted_bam_file = temp(f"results/{sample}_{{contig}}.aligned.sorted.bam")
    resources:
        mem_mb = 10000,
        time = "04:00:00"
    threads: 2
    params:
        log = f"logs/{reference_name}/{sample}/alt_keep_prob{alt_keep_prob}/sort_modified/{{contig}}.out"
    conda:
        "../../envs/samtools.yaml"
    shell:
        "samtools sort -o {output.sorted_bam_file} -@ {threads} {input.bam_file}"
        
rule merge_modified: 
    input: 
        bam_files = MODIFIED_BAM_FILES
    output: 
        merged_bam_file = temp(f"results/{sample}.aligned.bam")
    resources:
        mem_mb = 50000,
        time = "12:00:00"
    threads: 1
    params:
        log = f"logs/{reference_name}/{sample}/alt_keep_prob{alt_keep_prob}/merge_modified.out", 
        bam_files_str = ' '.join(MODIFIED_BAM_FILES)
    conda:
        "../../envs/samtools.yaml"
    shell:
        "samtools merge -O BAM {output.merged_bam_file} {params.bam_files_str}"

rule sort_merged:
    input: 
        bam_file = rules.merge_modified.output.merged_bam_file
    output: 
        sorted_bam_file = f"results/{sample}.aligned.sorted.bam"
    resources:
        mem_mb = 80000,
        time = "12:00:00"
    threads: 2
    params:
        log = f"logs/{reference_name}/{sample}/alt_keep_prob{alt_keep_prob}/sort_merged.out"
    conda:
        "../../envs/samtools.yaml"
    shell:
        "samtools sort -o {output.sorted_bam_file} -@ {threads} {input.bam_file}"