# Install packages
pacman::p_load('Biostrings', 'rentrez', 'tidyverse')

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/jgi_genes_round_2/")

# Pull accession numbers
accs <- c("WP_068691872.1", "WP_068691874.1", "WP_084012646.1", "WP_068691877.1", "WP_068691880.1", "WP_068691907.1")

pull <- entrez_fetch(id = c(accs),
                        db = "protein", rettype = "fasta", ap_key="826357f5ff17c7ec62e583909071e94f9d08")

write(pull, "output/Thermobifida_LstABCDEF.fasta")
fast <- readAAStringSet("output/Thermobifida_LstABCDEF.fasta")
names(fast)


dtf <- data.frame(gene = word(names(fast), 1, sep = " "),
                  organism = gsub("\\]", "", word(names(fast), 2, sep = "\\[")),
                  aa_seq_length = width(fast),
                  nucleotide_seq_length = width(fast) * 3,
                  seq_name = names(fast),
                  primary_amino_acid_sequence = as.character(fast))
head(dtf)
write.csv(dtf, "output/Thermobifida_LstABCDEF_seqs.csv")


lstc <- readAAStringSet("data/Thermobifida_seqs_for_jgi.fa")
fast <- lstc
testdtf <- data.frame(gene = word(names(fast), 1, sep = " "),
                      organism = gsub("\\]", "", word(names(fast), 2, sep = "\\[")),
                      aa_seq_length = width(fast),
                      nucleotide_seq_length = width(fast) * 3,
                      seq_name = names(fast),
                      primary_amino_acid_sequence = as.character(fast))

testdtf
write_csv(testdtf, "output/Thermobifida_df.csv")
