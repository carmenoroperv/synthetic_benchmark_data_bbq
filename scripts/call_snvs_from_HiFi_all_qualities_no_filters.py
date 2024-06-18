import pysam
import os
import sys
import py2bit
import statistics

####################################################################################################################################################################################
# flag_filter=3844 removes:
    # - read unmapped (0x4)
    # - not primary alignment (0x100)
    # - read fails platform/vendor quality checks (0x200)
    # - read is PCR or optical duplicate (0x400)
    # - supplementary alignment (0x800)
####################################################################################################################################################################################

def get_alleles_w_quals(pileupcolumn):
    base_quals = {'A':[], 'C':[], 'G':[], 'T':[]}
    for pileup_read in pileupcolumn.pileups:
        # test for deletion at pileup
        if pileup_read.is_del or pileup_read.is_refskip:
            continue

        # fetch read information
        allele = pileup_read.alignment.query_sequence[pileup_read.query_position].upper() 
        base_qual = pileup_read.alignment.query_qualities[pileup_read.query_position]
        if (allele not in "ATGC"):
            continue
        base_quals[allele].append(base_qual)
    return base_quals

def snv_call(bam_in, fa_in, chrom, start, end, min_BQ, min_MQ):
    ref_fa = fa_in.fetch(chrom)
    
    for pileupcolumn in bam_in.pileup(chrom, start, end, min_base_quality = min_BQ, min_mapping_quality = min_MQ, truncate=True):
        chrom = pileupcolumn.reference_name
        ref_pos = pileupcolumn.pos
        ref = ref_fa[ref_pos].upper()
        n_total = pileupcolumn.n     
        
        if ref not in 'ACGT':
            continue
        
        base_quals = get_alleles_w_quals(pileupcolumn)
        n = sum(len(base_quals[x]) for x in base_quals)
        n_ref = len(base_quals[ref])
        
        for A in [x for x in ['A','C','G','T'] if x != ref]:
            if len(base_quals[A]) == 0:
                continue
            median_BQ = statistics.median(base_quals[A])
            snv_calls.write(f'{chrom}\t{ref_pos+1}\t{ref}\t{A}\t{n_ref}\t{len(base_quals[A])}\t{n}\t{n_total}\t{median_BQ}\t{len(base_quals[A])/n_total}\n')

####################################################################################################################################################################################
####################################################################################################################################################################################

# log and printing
os.makedirs(os.path.dirname(snakemake.params.log), exist_ok=True)
sys.stderr = open(snakemake.params.log, "a+")
sys.stdout = sys.stderr

# input files
bam_in = snakemake.input.bam
bamfile = pysam.AlignmentFile(bam_in, "rb")
fasta_in = snakemake.input.fasta
fastafile = pysam.FastaFile(fasta_in)

# input parameters
chrom = snakemake.params.chrom
start = int(snakemake.params.start)
end = int(snakemake.params.end)
min_BQ = int(snakemake.params.min_BQ)
min_MQ = int(snakemake.params.min_MQ)

# output
os.makedirs(os.path.dirname(snakemake.output.snv_calls_contig), exist_ok=True)
with open(snakemake.output.snv_calls_contig, 'w') as fp:
    pass
snv_calls = open(snakemake.output.snv_calls_contig, 'a')

# call SNVs
snv_call(bamfile, fastafile, chrom, start, end, min_BQ, min_MQ)

    