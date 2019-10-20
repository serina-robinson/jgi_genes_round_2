# Install packages
pacman::p_load('readxl', 'janitor', 'dplyr', 'stringr', 'tidyverse', 'taxize')

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/jgi_oleA/")

# First pull the original JGI genes
# perc_id <- read_csv("data/percent_identity_OleA_characterized.csv", col_names = F) %>%
#   arrange(desc(X12)) %>%
#   group_by(X2) %>%
#   dplyr::select(X2) %>%
#   pull()

old_dat <- read_excel("data/Proposal_503776_DOE_JGI_genes_final.xlsx") %>%
  janitor::clean_names() %>%
  dplyr::filter(role == "OleA") %>%
  dplyr::mutate(orgs = gsub("\\]", "", word(gene, sep = "\\[", 2))) %>%
  dplyr::mutate(genus = word(orgs, sep = "_", 1)) %>%
  dplyr::select(gene, orgs, genus) 
  #dplyr::rename(aa_seq = seq)

# Read in the data
dat1 <- read_excel("data/JGI_freezer_stock_KEY.xlsx", sheet = 1) %>%
  janitor::clean_names() %>%
  dplyr::filter(users_name %in% old_dat$gene) %>%
  dplyr::rename(user_id = users_name) %>%
  dplyr::mutate(master_well = paste0("1-", master_well))


colnames(dat1) <- tolower(colnames(dat1))
dim(dat1)
# colnames(dat1)[1:2] <-  c("jgi name", "user ID")
# colnames(dat1) # Actually dat1 is the same as JGI round 1

dat2 <- read_excel("data/JGI_freezer_stock_KEY.xlsx", sheet = 2)
colnames(dat2) <- tolower(colnames(dat2))
colnames(dat2)

dat3 <- read_excel("data/JGI_freezer_stock_KEY.xlsx", sheet = 3)
colnames(dat3) <- tolower(colnames(dat3))
colnames(dat3)

# Select the sequences that were screened
seq1 <- paste0("A", 6:10)
seq2 <- paste0(rep(LETTERS[2:8], 6), rep(5:10, 6))
seq3 <- c(seq1, seq2)
seq4 <- seq3[-grep("D9", seq3)] # remove the alpha/beta hydrolase
length(seq4) # 46
seq4

# Pull these from the original order
colnames(dat2)

# Combine all the data
dat <- dat2 %>%
  bind_rows(., dat3) %>%
  janitor::clean_names() %>%
  dplyr::filter(grepl(paste0(seq4, collapse = "|"), master_well)) %>%
  dplyr::bind_rows(dat1) %>%
  dplyr::mutate(orgs = gsub("\\]", "", word(user_id, sep = "\\[", 2))) %>%
  dplyr::mutate(genus = word(orgs, sep = " ", 1)) %>%
  dplyr::select(user_id, orgs, genus, master_well) %>%
  dplyr::rename(gene = user_id) %>%
  dplyr::filter(!grepl("H291A", gene)) 


write_csv(dat, "data/72_OleA_masterwell_org_key.csv")
write_csv(dat, "~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/data/72_OleA_masterwell_org_key.csv")


