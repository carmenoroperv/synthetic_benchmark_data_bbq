import sys
import os
import json

configfile: "config/call_SNVs_HiFi.yaml"
workdir: "../"

# input
sample_name = config["hifi_sample"]["name"]
bam_file = config["hifi_sample"]["bam_file"]

fasta_file = config["reference"]["fasta_file"]
fai_file = config["reference"]["fai_file"]
twobit_file = config["reference"]["twobit_file"]
splits_bed_file = config["reference"]["split_bed"]

min_BQ = config["var_calling_params"]["min_BQ"]
min_MQ = config["var_calling_params"]["min_MQ"]

# select only contigs that are longer than 10**6
contig_dict = {}
with open(fai_file) as contigs:
    for line in contigs:
        contig_dict[line.split()[0]] = line.split()[1]
long_contigs = []
for contig_name, length in contig_dict.items():
    if int(length) > 10**6:
        long_contigs.append(contig_name)

contigs = []
starts = []
ends = []
with open(splits_bed_file) as f:
    for line in f:
        chrom, start, end = line.split()
        if str(chrom) in long_contigs:
            contigs.append(chrom)
            starts.append(start)
            ends.append(end)

SNV_CALLS = expand(f"results/split/contig_{{contig}}_{{start}}_{{end}}.txt", zip, contig = contigs, start = starts, end = ends)
SNV_CALLS_CONCAT = f"results/snv_calls_minBQ{min_BQ}_minMQ{min_MQ}.txt"

SNV_CALLS_ALL_QUAL = expand(f"results/split_all/contig_{{contig}}_{{start}}_{{end}}.txt", zip, contig = contigs, start = starts, end = ends)
SNV_CALLS_ALL_QUAL_CONCAT = f"results/snv_calls_allQual_noFilt.txt"

rule all:
    input: SNV_CALLS_CONCAT, SNV_CALLS_ALL_QUAL_CONCAT
    
rule call_snvs: 
    input:
        bam = bam_file,
        fasta = fasta_file
    output: 
        snv_calls_contig = temp(f"results/split/contig_{{contig}}_{{start}}_{{end}}.txt")
    resources: 
        mem_mb = 10000, 
        time = "12:00:00"
    params: 
        log = f"logs/{sample_name}/calling_hifi_snvs/call_snvs/contig_{{contig}}_{{start}}_{{end}}.out",
        chrom = '{contig}',
        start = '{start}',
        end = '{end}',
        min_BQ = min_BQ,
        min_MQ = min_MQ,
    conda: "../../envs/call_hifi_snvs.yaml"
    script: 
        "../scripts/call_snvs_from_HiFi.py"

rule concat_snv_calls:
    input: 
        snp_calls_contigs = SNV_CALLS
    output: 
        combined_calls = f"results/snv_calls_minBQ{min_BQ}_minMQ{min_MQ}.txt"
    resources: 
        mem_mb = 8000, 
        time = "04:00:00"
    params: 
        log = f"logs/{sample_name}/calling_hifi_snvs/concat_snv_calls.out",
        input_list = ' '.join(SNV_CALLS)
    conda: "../../envs/call_hifi_snvs.yaml"
    shell:
        "cat {params.input_list} > {output.combined_calls}"

rule call_snvs_all_qualities_no_filters:
    input:
        bam = bam_file,
        fasta = fasta_file, 
        splitfile = splits_bed_file
    output: 
        snv_calls_contig = temp(f"results/split_all/contig_{{contig}}_{{start}}_{{end}}.txt"),
    resources: 
        mem_mb = 20000, 
        time = "12:00:00"
    params: 
        log = f"logs/{sample_name}/calling_hifi_snvs/call_snvs_all_qualities_no_filters/contig_{{contig}}_{{start}}_{{end}}.out",
        chrom = '{contig}',
        start = '{start}',
        end = '{end}',
        min_BQ = 1,
        min_MQ = 1,
    conda: "../../envs/call_hifi_snvs.yaml"
    script: 
        "../scripts/call_snvs_from_HiFi_all_qualities_no_filters.py"

rule concat_snvs_all_qualities_no_filters:
    input: 
        snp_calls_contigs = SNV_CALLS_ALL_QUAL
    output: 
        combined_calls = f"results/snv_calls_allQual_noFilt.txt"
    resources: 
        mem_mb = 8000, 
        time = "04:00:00"
    params: 
        log = f"logs/{sample_name}/calling_hifi_snvs/concat_snvs_all_qualities_no_filters.out",
        input_list = ' '.join(SNV_CALLS_ALL_QUAL)
    conda: "../../envs/call_hifi_snvs.yaml"
    shell:
        "cat {params.input_list} > {output.combined_calls}"

