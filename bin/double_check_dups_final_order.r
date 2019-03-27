# Install packages
pacman::p_load('ape', 'rentrez', 'stringr', 'readxl', 'stringr', 'tidyverse', 'stringr', 'Biostrings', 'DECIPHER', 'phangorn', 'ggplot2', 'seqinr', 'bgafun')

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/jgi_genes_round_2/")

# single
sing <- read_excel("final_jgI_order/single_gene_vecs_for_jgi_20190322.xlsx")
duplicated(sing$primary_amino_acid_sequence)

duet <- read_excel("final_jgI_order/all_duet_vecs_for_jgi_20190322.xlsx")
duplicated(duet$mcs1_insert_primary_amino_acid_seq)
duplicated(duet$mcs2_insert_primary_amino_acid_seq)
intersect(duet$mcs1_insert_primary_amino_acid_seq, duet$mcs2_insert_primary_amino_acid_seq)
