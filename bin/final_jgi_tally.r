# Install packages
pacman::p_load('ape', 'rentrez', 'stringr', 'readxl', 'stringr', 'tidyverse', 'stringr', 'Biostrings', 'DECIPHER', 'phangorn', 'ggplot2', 'seqinr', 'bgafun')

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/jgi_genes_round_2/")

# Read in the jgi genes
morgan1 <- read_csv("data/JGI_DuetPlasmids_with_Nocardia.csv") %>%
  janitor::clean_names()

morgan2 <- read_excel("data/JGI_genes_to_order_BCs_listed_20190314_NoVecOleB.xlsx") %>%
  janitor::clean_names()

megan <- read_excel("data/Smith_JGI_OleA_genes.xlsx") %>%
  janitor::clean_names()

cmbnd <- megan %>%
  mutate(vector = "pET28b+") %>%
  bind_rows(morgan2) %>%
  dplyr::filter(!duplicated(primary_amino_acid_sequence)) %>%
  dplyr::select(-seq_name) %>%
  dplyr::filter(!is.na(role))

# write_csv(cmbnd, "output/total_jgi_order.csv")

# Get the two structural ANL superfamily proteins
struct <- read_excel("data/anl_training_set_updated_20190215.xlsx")
checkC <- struct[grep("Brucella|Dyadobacter", struct$organism),] %>%
  dplyr::select(entry_name, organism, length, protein_names, aa_seq) %>%
  dplyr::mutate(gene = paste0(entry_name, "_", protein_names)) %>%
  dplyr::select(-entry_name, -protein_names) %>%
  dplyr::rename(aa_seq_length = length) %>%
  dplyr::rename(primary_amino_acid_sequence = aa_seq) %>%
  dplyr::mutate(role = "OleC",
                x6x_his_tag_terminus = "C-term 6XHIS",
                cut_sites = "Nde1/Xho1",
                vector = "pET30b+",
                add_stop_codon = "after C_terminal_His_tag") %>%
  dplyr::mutate(nucleotide_seq_length = as.numeric(aa_seq_length) * 3)
head(checkC)


cmbndp <- cmbnd %>%
  bind_rows(checkC)
write_csv(cmbndp, "output/total_jgi_order.csv")

# colnames(nltC)
# # Get the entire nocardiolactone cluster in Nocardia brasiliensis
# nltC <- read_csv("data/NltC_trimmed_co_occur.csv") %>%
#   janitor::clean_names() %>%
#   dplyr::slice(1:12) %>%
#   mutate(role = case_when(#pfam_id1 == "PF00561" ~ "OleB",
#                                  #pfam_id1 == "PF08545" ~ "OleA",
#                                  #pfam_id1 == "PF08541" ~ "OleA",
#                                  #pfam_id1 == "PF01073" ~ "OleD",
#                                  #pfam_id1 == "PF01370" ~ "OleD",
#                                  pfam_id1 == "PF00501" ~ "OleC",
#                                  TRUE ~ NA_character_)) %>%
#   dplyr::filter(grepl("Ole", role)) %>%
#   dplyr::rename(organism = genus_species) %>%
#   dplyr::mutate(aa_seq_length = abs(end - start)) %>%
#   dplyr::mutate(#role = "OleC",
#                 x6x_his_tag_terminus = "C-term 6XHIS",
#                 cut_sites = "Nde1/Xho1",
#                 vector = "pET28b+",
#                 add_stop_codon = "after C_terminal_His_tag") %>%
#                 #nucleotide_seq_length = aa_seq_length * 3) %>%
#   dplyr::select(protein_acc, organism, role, x6x_his_tag_terminus, cut_sites, vector, add_stop_codon)
# 
# nltABD <- read_csv("data/NltC_trimmed_co_occur.csv") %>%
#   janitor::clean_names() %>%
#   dplyr::slice(1:12) %>%
#   mutate(role = case_when(pfam_id1 == "PF00561" ~ "OleB",
#     pfam_id1 == "PF08545" ~ "OleA",
#     pfam_id1 == "PF08541" ~ "OleA",
#     pfam_id1 == "PF01073" ~ "OleD",
#     pfam_id1 == "PF01370" ~ "OleD",
#     # pfam_id1 == "PF00501" ~ "OleC",
#     TRUE ~ NA_character_)) %>%
#   dplyr::filter(grepl("Ole", role)) %>%
#   dplyr::rename(organism = genus_species) %>%
#   dplyr::mutate(x6x_his_tag_terminus = "N-term 6XHIS",
#                 cut_sites = "Nde1/Xho1",
#                 vector = "pET28b+",
#                 add_stop_codon = "before Xho1") %>%
#   dplyr::select(protein_acc, organism, role, x6x_his_tag_terminus, cut_sites, vector, add_stop_codon)
# 
# #nltabcd <- entrez_fetch(id = c("WP_042260942.1", "WP_042260944.1", "WP_042260945.1", "WP_042260949.1"), 
# #                        db = "protein", rettype = "fasta", ap_key="826357f5ff17c7ec62e583909071e94f9d08")
#                          
# #write(nltabcd, "output/NltABCD.fasta")
# fast <- readAAStringSet("output/NltABCD.fasta")
# head(fast)
# 
# # WP_042260949.1 NltD
# # WP_042260945.1 NltC
# # WP_042260944.1 NltB
# # WP_042260942.1 NltA
# head(fast)
# 
# 
# nltABD
# # Combine 
# noct <- nltC %>%
#   bind_rows(nltABD) 
# 
# 
# df2 <- data.frame(names(fast), fast)[c(3, 4, 2, 1),]
# names(fast) <- NULL
# # df2 <- data.frame(fast)
# df3 <- data.frame(fast)[c(3, 4, 2, 1),]
# df4 <- cbind(noct, df2)
# df4
# write.csv(df4, "output/NltABCD_output.csv")


