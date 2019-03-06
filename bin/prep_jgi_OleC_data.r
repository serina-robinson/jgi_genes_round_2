## Install packages
pacman::p_load('ape', 'readxl', 'stringr', 'Biostrings', 'DECIPHER', 'phangorn', 'ggplot2', 'seqinr', 'bgafun')

# Set seed 
set.seed(123091)

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/jgi_genes_round_2/")

# Read in the data
rd <- read_excel("data/")