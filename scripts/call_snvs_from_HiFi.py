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

        # fetch sequence and quality information across reads in the given position
        allele = pileup_read.alignment.query_sequence[pileup_read.query_position].upper() 
        base_qual = pileup_read.alignment.query_qualities[pileup_read.query_position]
        if (allele not in "ATGC"):
            continue
        base_quals[allele].append(base_qual)
    return base_quals

def snv_call(bam_in, fa_in, chrom, start, end, min_BQ, min_MQ):
    n_snvs = 0                          # counting the number of called SNVs
    indel_positions = []                # registering the indel positions
    snv_positions_list = []             # registering the snv positions
    snv_positions = {}                  # registering the called snvs  
    
    ref_fa = fa_in.fetch(chrom)
    for pileupcolumn in bam_in.pileup(chrom, start, end, min_base_quality = min_BQ, min_mapping_quality = min_MQ, truncate=True, flag_filter = 3844):
        chrom = pileupcolumn.reference_name                                                 # chromosome name
        ref_pos = pileupcolumn.pos                                                          # position in chromosome
        ref =  ref_fa[ref_pos].upper()                                                      # reference base
        position = [x.upper() for x in pileupcolumn.get_query_sequences(add_indels = True)] # all bases in the current position
        base_quals = get_alleles_w_quals(pileupcolumn)                                      # dictionary with base qualities
        n = sum(len(base_quals[x]) for x in base_quals)                                     # coverage
        n_unique_query_seq = len(set(position))                                             # number of unique bases (including deletions)
        n_unqiue_bases = len(set([i for i in position if i in 'ACGT']))                     # number of unique bases (base in [ACGT])
        n_bases = len([i for i in position if i in 'ACGT'])                                 # total number of bases (base in [ACGT])
        n_total = pileupcolumn.n                                                            # total coverage (including indels)

        # indel position is defined as a position: 
            # with more than 1 unique sequences, 
            # where at least 1 of the sequences is a nucleotide base, 
            # coverage is above 20,
            # nucleotide bases make up less than 90% of the coverage.
        if n_unique_query_seq > 1 and n_unqiue_bases >= 1 and n_bases / n_total < 0.9 and n_total > 20: 
            indel_positions.append(ref_pos+1)

        # check that all nucleotide bases have base qualities reported 
        
        # require that reference base is in [ACGT]
        if ref not in 'ACGT':
            continue
        
        alternatives = []
        for A in [x for x in ['A','C','G','T'] if x != ref]:   # iterate over possible alternative alleles
            if len(base_quals[A]) != 0:                        # check and save, if the alternative is present
                alternatives.append(A)
        
        # require that only 1 alternative is present, less than 10% of reads contain indels and coverage is equal or larger to 10
        if len(alternatives) > 1 or len(alternatives) == 0 or n/n_total < 0.9 or n < 10:
            continue
        
        alt = alternatives[0]
        n_alt = len(base_quals[alt])
        n_ref = len(base_quals[ref])
        repeat_block = ref_fa[ref_pos-16:ref_pos+16]
        # require that 
            # both alt and ref allele have allele frequency > 0.1
            # variant is not in a repeat block
                # Repeat block def.: considering 4-mers in +-16 bp window around position. If the number of unique 4-mers is smaller than 16, region is defined as repeat block
        if n_alt/n_total < 0.1 or n_ref/n_total < 0.1 or len(set([repeat_block[i:i+4] for i in range(len(repeat_block))])) < 16:
            continue

        # save snv information in a dictionary and snv position in addition to a list
        snv_positions[str(ref_pos+1)] = [str(chrom), str(ref_pos+1), ref, alt, str(n_ref), str(n_alt), str(n), str(n_total), str(statistics.median(base_quals[alt])), str(n_alt/n_total)]
        snv_positions_list.append(ref_pos+1)

    # require that variant is not closer than 25 bp to an indel position
    for i in indel_positions: 
        for j in snv_positions_list:
            if abs(i-j) <= 25:
                if str(j) in snv_positions:
                    del snv_positions[str(j)]
            if j-i > 25: 
                break
	
    for snv in snv_positions:
        snv_str = "\t".join(snv_positions[snv])
        snv_calls.write(f'{snv_str}\n')
        n_snvs += 1
    print(f"Found {n_snvs} SNVs")
    #print(f"{n_BQ_missing} positions had missing BQ")

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

    
    