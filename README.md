# Synthetic low-frequency variant set for Illumina data

This repo contains scripts to create a set of low-frequency variants in Illumina sequencing data for benchmarking variant calling tools. 


The low-frequency variant set is created based on two samples obtained from the same biopsy and sequenced with Illumina (NovaSeq) and PacBio HiFi.  


The `call_snvs_from_HiFi.py` and `call_snvs_from_HiFi_all_qualities_no_filters.py` scripts call heterozygous SNVs from the PacBio HiFi data.

Subsequently, the `change_VAF_in_Illumina.py` script changes the VAF of the called SNVs in the Illumina data to create low-frequency variants. 

The recommended snakemake workflows for running these scripts can be found in the `snakefiles` folder.
