import os
import sys
import pandas as pd
import pysam
import random

def open_bam_w_index(input_bam): 
    bai_check1 = f"{input_bam}.bai"
    bai_check2 = input_bam[:-1] + 'i'
    if not os.path.exists(bai_check1) and not os.path.exists(bai_check2):
        pysam.index(input_bam)
    return pysam.AlignmentFile(input_bam, "rb")


class Read: # modified from BBQ github
    def __init__(self, pileup_read):
        # set attributes
        self.pos = pileup_read.query_position
        self.secondary = pileup_read.alignment.is_secondary
        self.supplementary = pileup_read.alignment.is_supplementary
        self.allel = pileup_read.alignment.query_sequence[self.pos]
        self.full_sequence = pileup_read.alignment.query_sequence
        self.base_qual = pileup_read.alignment.query_qualities[self.pos]
        self.query_name = pileup_read.alignment.query_name
        self.isR1 = pileup_read.alignment.is_read1
        
        # Process cigar stats
        cigar_stats = pileup_read.alignment.get_cigar_stats()[0]
        self.has_indel = sum(cigar_stats[1:4]) != 0


def change_alt_to_ref(change_dict, snv_id, ref, alt, read, mem_read = None):            
    if read.query_name not in change_dict: 
        change_dict[read.query_name] = []
    change_dict[read.query_name].append((snv_id, ref, alt, read.pos, read.base_qual, read.allel, read.isR1, len(read.full_sequence)))

    if mem_read != None: 
        change_dict[mem_read.query_name].append((snv_id, ref, alt, mem_read.pos, mem_read.base_qual, mem_read.allel, mem_read.isR1, len(mem_read.full_sequence)))
        n_exp_ch = 2
    else: 
        n_exp_ch = 1

    return(change_dict, n_exp_ch)

# log and printing
os.makedirs(os.path.dirname(snakemake.params.log), exist_ok=True)
sys.stderr = open(snakemake.params.log, "a+")
sys.stdout = sys.stderr

# input parameters
alt_keep_prob = float(snakemake.params.alt_keep_prob)
print(f'Probability of keeping the alternative allele: {alt_keep_prob}')
chrom = str(snakemake.params.contig)
print(f'Contig name from snakefile: {chrom}')

# input files
bam_in = snakemake.input.bam_file
bam = open_bam_w_index(bam_in)
snvs_in = snakemake.input.snvs_file
snvs = pd.read_csv(snvs_in, sep = "\t", header = None)
snvs.columns = ["chrom", "pos", "ref", "alt", "n_ref", "n_alt", "n", "n_total", "median_alt_BQ", "AF"]
snvs = snvs.loc[snvs['chrom'].isin([chrom])] 
print(f'SNV chroms, after filtering: {snvs.chrom.unique()}')
print(f'n SNVs after filtering: {snvs.shape[0]}')


change_dict = {}
n_exp_single = 0
n_exp_match_pairs = 0
n_exp_mismatch_pairs = 0
n_alt_not_obs = 0
change_fractions = []
n_snvs = 0
print_prog = False
for index, row in snvs.iterrows():
    if index % 10000 == 0:
        print_prog = True
    else: 
        print_prog = False
    chr_i = row["chrom"]
    pos_i = row["pos"]
    ref_i = row["ref"]
    alt_i = row["alt"]
    snv_id = f'{chr_i}:{str(pos_i)}_{ref_i}>{alt_i}'
    # check that alternative has length 1
    if len(alt_i) > 1:
        continue

    n_snvs += 1
    for pileupcolumn in bam.pileup(str(chr_i), int(pos_i-1), int(pos_i), truncate = True, ignore_overlaps = False, stepper = "nofilter", min_base_quality = 0, compute_baq = False):
        if print_prog:
            print(f'SNP nr {n_snvs} --> Chr: {chr_i}, Pos: {pos_i}, SNV: {ref_i}->{alt_i}')
        reads_mem = {}

        n_pairs_total = 0
        n_pairs_ch = 0
        n_pairs_alt_match = 0
        n_pairs_ref_match = 0
        n_pairs_other_match = 0
        n_pairs_alt_ref_mismatch = 0
        n_pairs_other_mismatch = 0
        
        for pileup_read in pileupcolumn.pileups:
            # fetch read information
            if pileup_read.is_del or pileup_read.is_refskip:
                continue
            read = Read(pileup_read)
            
            # skip AltToRef change if read is del, refskip, has an indel at the next position or is a secondary or supplementary alignment
            if read.has_indel or read.secondary == True or read.supplementary == True:
                continue

            # pick out all the read-pairs that overlap the snv position
            if read.query_name in reads_mem:
                n_pairs_total += 1
                # found partner -> process read pair
                mem_read = reads_mem.pop(read.query_name)
                if read.allel == mem_read.allel:  # overlap matches
                    if read.allel == alt_i:       # ALT match in overlap
                        n_pairs_alt_match += 1
                        change_status = random.random()
                        if change_status >= alt_keep_prob:
                            n_pairs_ch += 2  
                            change_dict, n_ch_i = change_alt_to_ref(change_dict, snv_id, ref_i, alt_i, read, mem_read) # change both ALT (in ALT/ALT overlap) to REF
                            n_exp_match_pairs += n_ch_i
                    elif read.allel == ref_i:       # REF match in overlap
                        n_pairs_ref_match += 1
                    else: 
                        n_pairs_other_match += 1  # OTHER match in overlap
                else:                             # overlap mismatch     
                    if read.allel == alt_i: 
                        if mem_read.allel == ref_i: # ALT -> REF mismatch
                            n_pairs_alt_ref_mismatch += 1
                            change_status = random.random()
                            if change_status >= alt_keep_prob:
                                n_pairs_ch += 1
                                change_dict, n_ch_i = change_alt_to_ref(change_dict, snv_id, ref_i, alt_i, read, None) # change ALT (in ALT/REF overlap) to REF
                                n_exp_mismatch_pairs += n_ch_i
                    elif read.allel == ref_i:
                        if mem_read.allel == alt_i: # REF -> ALT mismatch
                            n_pairs_alt_ref_mismatch += 1
                            change_status = random.random()
                            if change_status >= alt_keep_prob:
                                n_pairs_ch += 1
                                change_dict, n_ch_i = change_alt_to_ref(change_dict, snv_id, ref_i, alt_i, mem_read, None) # change ALT (in REF/ALT overlap) to REF
                                n_exp_mismatch_pairs += n_ch_i
                    else:
                        n_pairs_other_mismatch += 1         
            else:
                reads_mem[read.query_name] = read

        n_single_total = 0 # n_s_read_total
        n_single_ch = 0 # n_s_alt_changed
        n_single_alt = 0
        n_single_ref = 0
        n_single_other = 0
        
        # Handle single reads ie. without read-pair overlap
        for read in reads_mem.values():
            n_single_total += 1
            if read.allel == alt_i:
                n_single_alt += 1
                change_status = random.random()
                if change_status >= alt_keep_prob:
                    n_single_ch += 1
                    change_dict, n_ch_i = change_alt_to_ref(change_dict, snv_id, ref_i, alt_i, read, None) # change ALT to REF
                    n_exp_single += n_ch_i
            elif read.allel == ref_i: 
                n_single_ref += 1
            else: 
                n_single_other += 1

        n_alt_sum = n_single_alt + (n_pairs_alt_match*2) + n_pairs_alt_ref_mismatch
        n_ch = n_single_ch + n_pairs_ch
        if n_alt_sum == 0:
            n_alt_not_obs += 1
        else:
            ch_frac_i = round(n_ch/n_alt_sum, 4)
            change_fractions.append(ch_frac_i)
        
        if print_prog:
            print(f"Observed {n_pairs_total} pairs (AA: {n_pairs_alt_match}, RR: {n_pairs_ref_match}, RA: {n_pairs_alt_ref_mismatch}, BB: {n_pairs_other_match}, BX: {n_pairs_other_mismatch}) and {n_single_total} singles (A: {n_single_alt}, R: {n_single_ref}, B: {n_single_other}); Modification made to {n_ch} reads out of {n_alt_sum} ALT reads.")
            print(f'\n')
        sys.stdout.flush()

# close input file
bam.close()

##########################################################################################################################################
##########################################################################################################################################

print(f'Expect to see {n_exp_match_pairs + n_exp_single + n_exp_mismatch_pairs} changes: {n_exp_match_pairs} changes in ALT match pairs, {n_exp_mismatch_pairs} changes in REF/ALT mismatch pairs and {n_exp_single} in single reads')
print(f'Out of the total number of SNVs {n_snvs}, in {n_alt_not_obs} no alternative alleles was observed (and no changes will be made)')
print(f'\n')
sys.stdout.flush()

mod_bam = open_bam_w_index(bam_in)
changes = {} # actual_changes
n_query_qual_missing = 0
n_filtered_wrong_seq_length = 0
n_filtered_wrong_qual_length = 0

# Open outputfile for writing
os.makedirs(os.path.dirname(snakemake.output.bam_file), exist_ok=True)
out_file = pysam.AlignmentFile(snakemake.output.bam_file, "wb", template=mod_bam)


for read in mod_bam.fetch(chrom): 
    read_id = read.query_name

    if read_id in change_dict:
        for data in change_dict[read_id]: # (snv_id, ref_i, alt_i, read_i.pos, read_i.base_qual, read_i.allel, read_i.isR1, read.full_sequence)
            snv_id = data[0]
            ref_new = data[1]
            alt_old = data[2]
            pos = data[3]
            basequal = data[4]
            allele = data[5]
            isR1 = data[6]
            len_seq = data[7]

            if read.is_read1 == isR1:
                if read.is_secondary or read.is_supplementary:
                    continue
                if len(read.query_sequence) != len_seq:
                    n_filtered_wrong_seq_length += 1
                    continue
                if len(read.query_qualities) != len_seq:
                    n_filtered_wrong_qual_length += 1
                    continue
             
                assert(read.query_sequence[pos] == alt_old)
                assert(read.query_sequence[pos] == allele)
                if read.query_qualities == None:
                    n_query_qual_missing += 1
                    continue

                q = read.query_qualities
                list_seq = list(read.query_sequence)
                list_seq[pos] = ref_new
                read.query_sequence = ''.join(list_seq)
                read.query_qualities = q

                if snv_id in changes:
                    if read_id in changes[snv_id]:
                        changes[snv_id][read_id] += 1
                    else: 
                        changes[snv_id][read_id] = 1
                else: 
                    changes[snv_id] = {}
                    changes[snv_id][read_id] = 1
    out_file.write(read)

pair_changes = 0
single_changes = 0
total_changes = 0
other_n_changes = 0
for snp_key, reads_val in changes.items():
    for read, n_changes in reads_val.items():
        if n_changes == 2:
            pair_changes += 2
        elif n_changes == 1:
            single_changes += 1
        else:
            other_n_changes += n_changes
        total_changes += n_changes


print(f'{n_filtered_wrong_seq_length} reads with wrong sequence length, update skipped')
print(f'{n_filtered_wrong_qual_length} reads with wrong qual length, update skipped')
print(f'Query qualities missing in {n_query_qual_missing} reads, update skipped')
print(f'In total {pair_changes + single_changes + other_n_changes} = {total_changes} changes made:  {pair_changes} pair changes and {single_changes} single changes, (other n changes {other_n_changes})')

if len(change_fractions) == 0:
    frac_alt_changed_final = 0
else: 
    frac_alt_changed_final = round(sum(change_fractions)/len(change_fractions), 4)
    
print(f'Average fraction of alt changed to ref: {frac_alt_changed_final}')
sys.stdout.flush()

mod_bam.close()
out_file.close()





    