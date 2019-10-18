# Install packages
pacman::p_load('readxl', 'janitor', 'dplyr', 'stringr', 'tidyverse', 'taxize', 'DECIPHER', 'Biostrings')

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/jgi_genes_round_2/")

# Pull the corresponding sequences from the old JGI set
old_dat <- read_excel("data/Proposal_503776_DOE_JGI_genes_final.xlsx") %>%
  janitor::clean_names() %>%
  dplyr::filter(role == "OleC") %>%
  dplyr::mutate(organism = gsub("\\]", "", word(gene, sep = "\\[", 2))) %>%
  dplyr::mutate(genus = word(organism, sep = "_", 1)) %>%
  dplyr::select(gene, seq, organism) %>%
  dplyr::rename(aa_seq = seq) 
old_dat$organism <- gsub("_", " ", old_dat$organism)


# Pull the corresponding sequences from the new JGI set
new_dat <- read_excel("data/single_gene_vecs_for_jgi_20190322_round2.xlsx") %>%
  janitor::clean_names() %>%
  dplyr::filter(role == "OleC" | role == "OleBC") %>%
  dplyr::rename(aa_seq = primary_amino_acid_sequence) %>%
  dplyr::mutate(genus = word(organism, sep = " ", 1)) %>%
  dplyr::select(gene, aa_seq, organism)
new_dat
new_dat$organism
