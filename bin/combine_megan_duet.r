# Install packages
pacman::p_load('ape', 'rentrez', 'stringr', 'readxl', 'stringr', 'tidyverse', 'stringr', 'Biostrings', 'DECIPHER', 'phangorn', 'ggplot2', 'seqinr', 'bgafun')

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/jgi_genes_round_2/")

# Read in the mobilicoccus & granulosicoccus
gran <- read_excel("data/")
grep("Granulosicoccus|Mobilicoccus", gran$)

# Read in the jgi genes
morgan1 <- read_csv("data/JGI_DuetPlasmids_with_Nocardia.csv") %>%
  janitor::clean_names()

megan <- read_excel("data/Smith_JGI_OleA_genes.xlsx") %>%
  janitor::clean_names() %>%
  dplyr::filter(!duplicated(primary_amino_acid_sequence)) %>%
  dplyr::filter(grepl("Lst", role)) 

megan_mutated <- megan %>%
  dplyr::select(role, gene, organism, x6x_his_tag_terminus, cut_sites, primary_amino_acid_sequence)

lsta <- megan_mutated %>%
  dplyr::filter(role == "LstA")

lstb <- megan_mutated %>%
  dplyr::filter(role == "LstB")
dim(lstb)

megan_comb <- data.frame(organism = lsta$organism, vector = "Pcola-Duet1 (kan resist, T7 IPTG inducible promoter)",
           mcs1_insert = lsta$gene, mcs2_insert = lstb$gene, combination_type = "LstA + LstB",
           mcs1_tags = "None", mcs2_tags = "None", mcs1_cut_sites = "NcoI/HindIII", mcs2_cut_sites = "NdeI/XhoI",
           mcs1_insert_primary_amino_acid_seq = lsta$primary_amino_acid_sequence,
           mcs2_insert_primary_amino_acid_seq = lstb$primary_amino_acid_sequence)

# Read in the jgi genes
final_comb <- read_csv("data/JGI_DuetPlasmids_with_Nocardia.csv") %>%
  janitor::clean_names() %>%
  bind_rows(megan_comb)

wid1 <- width(final_comb$mcs1_insert_primary_amino_acid_seq)
wid2 <- width(final_comb$mcs2_insert_primary_amino_acid_seq)
sum(wid1, wid2) * 3 # 100 kbp

head(final_comb)
write_csv(final_comb, "output/all_duet_vecs_for_jgi.csv")

# Read in the jgi genes
morgan2 <- read_excel("data/JGI_genes_to_order_BCs_listed_20190314_NoVecOleB.xlsx") %>%
  janitor::clean_names()

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
                vector = "pET30b",
                add_stop_codon = "after C_terminal_His_tag") %>%
  dplyr::mutate(nucleotide_seq_length = as.numeric(aa_seq_length) * 3)
head(checkC)

megan_ole <- read_excel("data/Smith_JGI_OleA_genes.xlsx") %>%
  janitor::clean_names() %>%
  dplyr::filter(!duplicated(primary_amino_acid_sequence)) %>%
  dplyr::filter(grepl("Ole", role)) %>%
  dplyr::mutate(vector = "pET28b") 

cmbndp <- morgan2 %>%
  bind_rows(checkC) %>% 
  bind_rows(megan_ole) %>%
  dplyr::filter(!is.na(role))

head(cmbndp)
cmbndp$nucleotide_seq_length
sum(cmbndp$nucleotide_seq_length) # 184323
 
184323 +  100233
cmbndp
write_csv(cmbndp, "output/single_gene_vecs_jgi.csv")
