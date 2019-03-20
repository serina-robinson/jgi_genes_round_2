# Install packages
pacman::p_load('ape', 'rentrez', 'stringr', 'readxl', 'stringr', 'tidyverse', 'stringr', 'Biostrings', 'DECIPHER', 'phangorn', 'ggplot2', 'seqinr', 'bgafun')

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/jgi_genes_round_2/")

# Read in the data
rd <- read_excel("data/20190304JGIgeneRequest.xlsx") %>%
  janitor::clean_names()

# 