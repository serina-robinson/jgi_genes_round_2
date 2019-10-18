# Install packages
pacman::p_load("data.table", "ggrepel", "phangorn", "readxl", "ranger", "muscle", "tidyverse", "Biostrings","DECIPHER","igraph","RColorBrewer", "officer", "rvg")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/jgi_genes_round_2/")

# Old JGI genes
old_jgi <- read_excel("data/Proposal_503776_DOE_JGI_genes_final.xlsx")
old_oled_df <- old_jgi %>%
  dplyr::filter(Role %in% c("OleD", "LstD")) 

old_oled <- old_oled_df %>%
  pull(seq) %>%
  AAStringSet()
names(old_oled) <- old_oled_df$gene

# New JGI genes
fin <- read_csv("final_jgI_order/single_gene_vecs_for_jgi_20190322.csv")
new_oled_df <- fin %>%
  dplyr::filter(role %in% c("OleD", "LstD"))

new_oled <- new_oled_df %>%
  pull(primary_amino_acid_sequence) %>%
  AAStringSet()
names(new_oled) <- new_oled_df$seq_name

# Known NltD
nltd <- readAAStringSet("~/Documents/University_of_Minnesota/Wackett_Lab/gblocks/NltD_OleD/NltD_raw.fasta")

# Known OleD
oled <- readAAStringSet("~/Documents/University_of_Minnesota/Wackett_Lab/github/jgi_genes_round_2/data/CDhit_OleD.fasta")

# Align all sequences
comb <- AAStringSet(c(nltd, oled, old_oled, new_oled))
dedup <- comb[!duplicated(comb)]
width(dedup)
aa.al <- AlignSeqs(comb)
writeXStringSet(comb, "data/NltD_OleD.fasta")
BrowseSeqs(aa.al)
writeXStringSet(aa.al, "data/NltD_OleD_aligned.fasta")

# Read the sequence
phy<-read.phyDat("data/NltD_OleD_aligned.fasta", format="fasta", type="AA")

# Make a distance matrix
dm<-dist.ml(phy)
mat <- as.matrix(dm)
mdist <- dist(dm)

# Multidimensional scaling using cmdscale
mds <- cmdscale(mdist, eig=TRUE, k=3)
x <- mds$points[,1]
y <- mds$points[,2]
z <- mds$points[,3]

# Calculate percent explained by each principal component
pc1 <- mds$GOF[1]
pc2 <- (mds$GOF[2] - mds$GOF[1])
pc1
pc2

# Make x- and y-axis labels
xlab<-paste0("PC 1 (", round( pc1 * 100, 2), "% tot. explained var.)")
ylab<-paste0("PC 2 (", round( pc2 * 100, 2), "% tot. explained var.)")

cdat <- data.frame(cbind(x, y))

# Make a colored 2D plot
numseqs <- nrow(cdat)

pdf(paste0("output/", numseqs,"_NltD_OleD_PCoA_colored.pdf"))
par(mar=c(0.01, 0.01, 0.01, 0.01))
ggplot(data = cdat, aes(x=x, y=y)) + 
  # geom_text(label = cdat$labl) +
  # geom_text(x = cdat$x, y = cdat$y, label = cdat$labl, check_overlap = T) +
  geom_point() +
  geom_text_repel(label = rownames(cdat), size = 1) +
  #geom_point(fill = cdat$colrs, size = cdat$sz, 
  #           alpha = cdat$alph, shape = 21) +
  # scale_color_manual(pal) +
  scale_shape(solid = TRUE) +
  labs(x=xlab,y=ylab) +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        legend.key = element_rect(fill = "white"),
        legend.background = element_rect(fill = "white"),
        legend.position = c(0.14, 0.80),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title=element_text(size=20, face="bold", hjust=0)
  )
dev.off()


# Which have OleBs
fin <- read_csv("final_jgI_order/single_gene_vecs_for_jgi_20190322.csv")
new_olebcd_df <- fin %>%
  dplyr::filter(role %in% c("OleD", "OleB", "OleBC", "LstD")) %>%
  arrange(role) %>%
  distinct(organism, .keep_all = T)
new_olebcd_df

writeXStringSet()