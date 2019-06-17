# Install packages
pacman::p_load('Biostrings', 'rentrez', 'tidyverse', 'stringr')

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/jgi_genes_round_2/")

# Read in the glycosyltransferases

dat <- read_csv("data/920_OleA/main_co_occur.csv") %>%
  janitor::clean_names()
  
gcsy <- dat %>%
  dplyr::filter(pfam_id1 == "PF00535") %>%
  distinct(query)

olecs <- dat %>%
  dplyr::filter(pfam_id1 == "PF00501") %>%
  distinct(query)

oleds <- dat %>%
  dplyr::filter(pfam_id1 %in% c("PF01370", "PF01073")) %>%
  distinct(query)

cyto <- dat %>%
  dplyr::filter(pfam_id1 == "PF00067") %>%
  distinct(query)

both <- intersect(pull(gcsy), intersect(pull(cyto), intersect(pull(olecs), pull(oleds))))

gcsy_nbs <- dat %>%
  dplyr::filter(query %in% both) 

unique(gcsy_nbs$query)
orgs2 <- gcsy_nbs %>%
  mutate(orgs = word(genus_species, 1, sep = " ")) %>%
  distinct(orgs)
 
write_delim(orgs2, "output/20_genera_with_glycosyltransferase_cytochrome.csv")
write_csv(gcsy_nbs, "output/glycosyltransferase_cytochrome_neighborhoods.csv")

abs <- gcsy_nbs %>%
  dplyr::filter(pfam_id1 == "PF00561" | pfam_id2 == "PF00561") %>%
  distinct(protein_acc)
abs
