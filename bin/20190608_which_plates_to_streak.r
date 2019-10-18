# Install packages
pacman::p_load('readxl', 'scales', 'janitor', 'dplyr', 'stringr', 'tidyverse', 'taxize', 'DECIPHER', 'Biostrings', 'ggtree', 'RColorBrewer')

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/jgi_genes_round_2/")

# Read in the taxonomy
tax <- read_excel("data/20190805_OleC_taxonomic_classification.xlsx")

# Read in the genes to screen
# Pull the corresponding sequences from the old JGI set
old_dat <- read_excel("data/Proposal_503776_DOE_JGI_genes_final.xlsx") %>%
  janitor::clean_names() %>%
  dplyr::filter(role == "OleC") %>%
  dplyr::mutate(organism = gsub("\\]", "", word(gene, sep = "\\[", 2))) %>%
  dplyr::mutate(genus = word(organism, sep = "_", 1)) %>%
  dplyr::select(gene, seq, genus, organism) %>%
  dplyr::rename(aa_seq = seq) 
old_dat$organism <- gsub("_", " ", old_dat$organism)

# Pull the corresponding sequences from the new JGI set
new_dat <- read_excel("data/single_gene_vecs_for_jgi_20190322_round2.xlsx") %>%
  janitor::clean_names() %>%
  dplyr::filter(role == "OleC" | role == "OleBC") %>%
  dplyr::rename(aa_seq = primary_amino_acid_sequence) %>%
  dplyr::mutate(genus = word(organism, sep = " ", 1)) %>%
  dplyr::select(gene, aa_seq, genus, organism)

merg_dat <- old_dat %>%
  bind_rows(new_dat) %>%
  mutate(acc = word(gene, sep = "\\.1", 1))
merg_dat$acc[merg_dat$acc == "C6W5A4_DYAFD_Acyl-CoA synthetase (AMP-forming)/AMP-acid ligase II-like protein"] <- "C6W5A4_DYAFD"

# Exclude the ones we don't want
# already_tested <- c("Xylella", "Saccharopolyspora", "Pseudonocardia", "Shewanella", "Chlamydiales")
# 
#                     #"Actinoplanes", "Saccharothrix")
# tr_dat <- merg_dat %>%
#   dplyr::filter(!grepl(paste0(already_tested, collapse = "|"), organism)) %>%
#   dplyr::filter(!grepl("WP_065302598|WP_053141833|WP_015749355", gene))
tr_dat <- merg_dat

to_test <- AAStringSet(tr_dat$aa_seq)
names(to_test) <- paste0(tr_dat$acc, "_", tr_dat$genus)

# writeXStringSet(to_test, "data/45_OleCs_to_screen.fasta")


# Find them in the JGI key
key <- read_excel("data/JGI_freezer_stock_KEY.xlsx", sheet = 1) %>%
  janitor::clean_names()
colnames(key)

dat2 <- read_excel("data/JGI_freezer_stock_KEY.xlsx", sheet = 2)
colnames(dat2) <- tolower(colnames(dat2))
colnames(dat2)[1:2] <- c("name", "users_name")
colnames(dat2) <- gsub(" ", "_", colnames(dat2))

dat3 <- read_excel("data/JGI_freezer_stock_KEY.xlsx", sheet = 3)
colnames(dat3) <- tolower(colnames(dat3))
colnames(dat3)[1:2] <- c("name", "users_name")
colnames(dat3) <- gsub(" ", "_", colnames(dat3))

dat4 <- key %>%
  bind_rows(dat2) %>%
  bind_rows(dat3) %>%
  mutate(acc = word(users_name, sep = "\\.1", 1))
dat4$users_name
dat4$acc

# Match up the desired genes with the well positions
notfound <- tr_dat$genus[!grepl(paste0(dat4$acc, collapse = "|"), tr_dat$acc)]
dat4$users_name[grepl(paste0(notfound, collapse = "|"), dat4$users_name)]

notfound
dat4$users_name
grep("Brucella", dat4$users_name)
tosearch <- c(tr_dat$acc, "DYAFD", "BRUC", "Isosphaera")

datun <- data.frame(unique(dat4$users_name[grep(paste0(tosearch, collapse = "|"), dat4$users_name)]))
# datun2 <- data.frame(unique(dat4$users_name[grep(paste0(tr_dat$acc, collapse = "|"), dat4$users_name)]))
# datun2
colnames(datun) <- "V1"


# Pull the key locations
findat <- dat4[dat4$users_name %in% datun$V1,]
write_csv(findat, "output/JGI_OleC_genes_to_plate.csv")


# datfin <- datun %>%
#   dplyr::filter(!grepl("NAD|oxoacyl|Oxoacyl|ketoacyl", V1))
# datfin
# 40 seqs

# dat4$users_name[grep("DYAFD", dat4$users_name)]

tr_dat$organism
# dat4$users_name[!grepl(paste0(tr_dat$acc, collapse = "|"), dat4$acc)]# 40 seqs

# 
