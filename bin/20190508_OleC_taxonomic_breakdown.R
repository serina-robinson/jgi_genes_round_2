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
  bind_rows(new_dat)

# Read in the OleA taxonomy
tax <- read_excel("data/OleA_taxonomic_classification.xlsx")
head(tax)

# Pull the taxonomy
# phy1 <- taxize::classification(merg_dat$genus[1:11], db = "ncbi", "826357f5ff17c7ec62e583909071e94f9d08")
# phy2 <- taxize::classification(merg_dat$genus[12:20], db = "ncbi", "826357f5ff17c7ec62e583909071e94f9d08")
# phy3 <- taxize::classification(merg_dat$genus[21:26], db = "ncbi")
# phy4 <- taxize::classification(merg_dat$genus[27:33], db = "ncbi")
# phy5 <- taxize::classification(merg_dat$genus[34:41], db = "ncbi")
# phy6 <- taxize::classification(merg_dat$genus[42:47], db = "ncbi")
# phy7 <- taxize::classification(merg_dat$genus[48:53], db = "ncbi")
phy <- c(phy1, phy2, phy3, phy4, phy5, phy6, phy7)


# Check all results are the appropriate length
tax_trim <- phy[unlist(lapply(phy, length)) == 3]

# Create a data frame of results
class <-  lapply(1:length(tax_trim), function(x) { tax_trim[[x]][grep("^class$", tax_trim[[x]]$rank),1] })
class_l <- unlist(lapply(1:length(class), function(x) ifelse(length(class[[x]]) == 0, "Unknown", class[[x]])))
order2 <-  lapply(1:length(tax_trim), function(x) { tax_trim[[x]][grep("^order$", tax_trim[[x]]$rank),1] })
order_l <- unlist(lapply(1:length(order2), function(x) ifelse(length(order2[[x]]) == 0, "Unknown", order2[[x]])))
phylum <-  lapply(1:length(tax_trim), function(x) { tax_trim[[x]][grep("^phylum$", tax_trim[[x]]$rank),1] })
phylum_l <-  unlist(lapply(1:length(phylum), function(x) ifelse(length(phylum[[x]]) == 0, "Unknown", phylum[[x]])))
genus <- lapply(1:length(tax_trim), function(x) { tax_trim[[x]][grep("^genus$", tax_trim[[x]]$rank),1] })
genus_l <- unlist(lapply(1:length(genus), function(x) ifelse(length(genus[[x]])== 0, "Unknown", genus[[x]])))
family <- lapply(1:length(tax_trim), function(x) { tax_trim[[x]][grep("^family$", tax_trim[[x]]$rank),1] })
family_l <- unlist(lapply(1:length(family), function(x) ifelse(length(family[[x]]) == 0, "Unknown", family[[x]])))

dtf <- data.frame(merg_dat$organism, merg_dat$genus, merg_dat$gene, merg_dat$aa_seq, phylum_l, class_l, order_l, family_l)
colnames(dtf)[1:4] <- c("organism", "genus", "gene", "aa_seq")

# Write to file
# write_csv(dtf, "output/OleC_taxonomic_classification.csv")

# Make alignment
exel <- read_excel("data/20190805_OleC_taxonomic_classification.xlsx")
exel$acc <- word(exel$gene, sep = "\\.1", 1) 
exel$acc
olecs <- AAStringSet(exel$aa_seq)
names(olecs) <- paste0(exel$acc, "_", exel$organism)
names(olecs) <- gsub(" ", "_", names(olecs))
names(olecs)

# Align sequences locally
aa.al <- AlignSeqs(olecs)
BrowseSeqs(aa.al)
aa.tr <- AAStringSet(substr(aa.al, 590, 1390))
BrowseSeqs(aa.tr)
length(aa.tr)
writeXStringSet(aa.tr, "data/57_JGI_OleC_seqs_aligned.fasta")

# FastTree
# FastTree -gamma <57_JGI_OleC_seqs_aligned.fasta> 57_JGI_OleC_seqs.nwk