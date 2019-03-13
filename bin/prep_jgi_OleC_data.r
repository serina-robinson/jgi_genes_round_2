# Install packages
pacman::p_load('ape', 'rentrez', 'stringr', 'readxl', 'stringr', 'tidyverse', 'stringr', 'Biostrings', 'DECIPHER', 'phangorn', 'ggplot2', 'seqinr', 'bgafun')

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/jgi_genes_round_2/")

# Read in the data
rd <- read_excel("data/20190304JGIgeneRequest.xlsx") %>%
  janitor::clean_names()
rd$cut_sites

# write_delim(data.frame(rd$Gene), "data/Gene_accession_numbers.txt", col_names = F)

# Search for flanking OleA, B, and D genes
# "PF08545" OleA
# "PF01073" OleD
# "PF00561" OleB

# Read in the RODEO results
rodeo_raw <- read_csv("data/RODEO_output_JGI_OleC/output/main_co_occur.csv") %>%
  janitor::clean_names() %>%
  mutate(ole_gene_id = case_when(pfam_id1 == "PF00561" ~ "OleB",
                                  pfam_id1 == "PF08545" ~ "OleA",
                                  pfam_id1 == "PF08541" ~ "OleA",
                                  pfam_id1 == "PF01073" ~ "OleD",
                                  pfam_id1 == "PF01370" ~ "OleD",
                                  pfam_id1 == "PF00501" ~ "OleC",
                                  TRUE ~ NA_character_)) %>%
  dplyr::filter(grepl("Ole", ole_gene_id))

# Get the Granulosicoccus and Mobilicoccus OleBC & D
old_oleAs <- read_csv("data/920_OleA/main_co_occur.csv") %>%
  janitor::clean_names() %>%
  dplyr::filter(grepl("Granulosicoccus|Mobilicoccus", genus_species)) %>%
  mutate(ole_gene_id = case_when(pfam_id1 == "PF00561" ~ "OleB",
                                 pfam_id1 == "PF01073" ~ "OleD",
                                 pfam_id1 == "PF01370" ~ "OleD",
                                 pfam_id1 == "PF00501" ~ "OleC",
                                 TRUE ~ NA_character_)) %>%
  dplyr::filter(grepl("Ole", ole_gene_id)) 

rodeo <- bind_rows(old_oleAs, rodeo_raw)
head(rodeo) 
dim(rodeo)

# writeXStringSet(aa, "data/jgi_olec_protein_seqs.fasta")

# How long are the rest of the Ole sequences
# other_oles <- rodeo$protein_acc[rodeo$ole_gene_id != "OleC"]
# unique(other_oles) # 41 other sequences
all_oles <- unique(rodeo$protein_acc)
all_oles
#ncbi_sqs <- entrez_fetch(all_oles, db = "protein", rettype = "fasta",
#             api_key = "826357f5ff17c7ec62e583909071e94f9d08")
# write(ncbi_sqs, "data/all_Ole_sqs_from_NCBI.fasta")
ncbi <- readAAStringSet("data/all_Ole_sqs_from_NCBI.fasta")
other_oles_width <- sum(width(ncbi)) * 3
other_oles_width
# (other_oles_width + olec_width)/1000 # 136 kbp
orgs <- gsub("\\]", "", word(names(ncbi), 2, sep = "\\["))

names(ncbi) <- NULL
table(rodeo$ole_gene_id)


final_df <- data.frame(role = na.omit(rodeo$ole_gene_id),
                       gene = all_oles,
                       organism = orgs,
                       aa_seq_length = width(ncbi),
                       nucleotide_seq_length = width(ncbi) * 3,
           x6x_his_tag_terminus = ifelse(na.omit(rodeo$ole_gene_id)=="OleC", "C-term 6XHIS", "N-term 6XHIS"),
           cut_sites = ifelse(na.omit(rodeo$ole_gene_id)=="OleC", "Nco1/Xho1", "Nde1/Xho1"),
           add_stop_codon = ifelse(na.omit(rodeo$ole_gene_id)=="OleC", "after C_terminal_His_tag", "before Xho1"),
           seq_name = names(readAAStringSet("data/all_Ole_sqs_from_NCBI.fasta")),
           primary_amino_acid_sequence = as.character(ncbi))


ord_df <- final_df[order(final_df$organism),]  %>%
  dplyr::filter(!grepl("pelag", organism)) %>%
  mutate(role_updated = case_when(aa_seq_length > 700) ~ "OleBC",
                          TRUE ~ role))

# write_csv(ord_df, "output/Morgan_JGI_genes_to_order.csv")  
# approximately 300 kbp

# Read in the current JGI order
jgi_old <- read_excel("data/Proposal_503776_DOE_JGI_genes_final.xlsx")
jgi_old_orgs <- gsub("\\]", "", word(jgi_old$gene, 2, sep = "\\["))
jgi_old_orgs
jgi_old_orgs %in% orgs
