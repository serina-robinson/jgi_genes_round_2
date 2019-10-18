# Install packages
pacman::p_load('readxl', 'scales', 'janitor', 'dplyr', 'stringr', 'tidyverse', 'taxize', 'DECIPHER', 'Biostrings', 'ggtree', 'RColorBrewer')

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/jgi_genes_round_2/")

# Run RaxML # Still running
# module load raxml/8.2.9_sse3_pthread 
# raxmlHPC-PTHREADS-SSE3 -f a -T 20 -p 1234 -x 1234 -m PROTGAMMAJTT -s 72_OleA_sequences_aligned.fasta -n 72_OleA_sequences_tree.nwk -# 100

# Alternativly run FastTree
# FastTree -gamma 72_OleA_sequences_aligned.fasta > 72_OleA_sequences_tree.nwk

# Read in the RaxML results
# nwk <- raxml2nwk("output/RAxML_bipartitionsBranchLabels.T15")
# raxml <- read.tree("output/230_FastTree.nwk")
#raxml <- read.raxml("output/RAxML_bipartitionsBranchLabels.T15")

# Read Fasttree
tr <- "data/57_JGI_OleC_seqs.nwk"
fastree <- treeio::read.tree(tr)
fastree$tip.label
fastree2 <- ape::drop.tip(fastree, grep("Bacillus", fastree$tip.label))

# Plot the tree using ggtree package
p <- ggtree(fastree2, layout = "circular") #color = "black",
# branch.length = 'none')
p$data$label

# Read in the taxonomy
tax <- read_excel("data/20190805_OleC_taxonomic_classification.xlsx")
tax$organism <- gsub(" ", "_", tax$organism)

ptax <- p$data$label
tax$acc <- NA
tax$acc
for(i in 1:length(tax$organism)) {
  ind <- grep(tax$organism[i], p$data$label)
  tax$acc[i] <- paste0(word(p$data$label[ind], sep = "_", 1), "_", word(p$data$label[ind], sep = "_", 2), "_", word(p$data$label[ind], sep = "_", 3))
  ptax[ind] <- paste0(p$data$label[ind], "_", tax$class_l[i], "_", tax$phylum_l[i])
}
tax[54:55,]
tax$acc[54] <- "C6W5A4_DYAFD_Acyl-CoA"
tax$acc[55] <- "A9MB96_BRUC2_AMP-dependent"
# C6W5A4_DYAFD_Acyl-CoA
# A9MB96_BRUC2_AMP-dependent

# Reorder the data frame to match
dt2 <- data.frame(p$data$label[p$data$isTip], word(p$data$label[p$data$isTip], sep = "\\.1", 1))
colnames(dt2) <- c("label", "acc")
dt2$acc <- paste0(word(dt2$acc, sep = "_", 1), "_", word(dt2$acc, sep = "_", 2), "_", word(dt2$acc, sep = "_", 3))
dt2$acc

# Merge with the tax df
merg <- dt2 %>%
  inner_join(tax, by = "acc")
merg$sz <- 0

dim(merg)
merg$label
active <- c("Xanthomonas", "Stenotrophomonas", "Nocardia", "Lysobacter")
inactive <- c("Chlamydiales_bacterium", "Bacteroides_luti", "Actinoplanes", "Saccharothrix_espanaensis", "Allokutzneria_albata", "Bacillus_sp")

p2$data$label
merg$is_active_colr <- "white"
merg <- merg %>%  
  dplyr::mutate(., is_active_colr = case_when(grepl(paste0(inactive, collapse = "|"), p2$data$label[p2$data$isTip]) ~ "maroon",
                                              grepl(paste0(active, collapse = "|"), p2$data$label[p2$data$isTip]) ~ "forestgreen",
                                              grepl("Photorhabdus_luminescens|Azospirillum_brasilense|C6W5A4_DYAFD", p2$data$label[p2$data$isTip]) ~ "white",
                                              TRUE ~ "white"))
p2$data$label

merg <- merg %>%
  dplyr::mutate(., is_active_sz = case_when(grepl(paste0(active, collapse = "|"), p2$data$label[p2$data$isTip]) ~ 1,
                      grepl(paste0(inactive, collapse = "|"), p2$data$label[p2$data$isTip]) ~ 1,
                      TRUE ~ 0))


data.frame(merg$organism, merg$is_active_colr)


# table(word(ptax, -2, sep = "_"))

# Annotate the phylogenetic tree
# mat <- data.frame(matrix(data = rep("black", 70*8), nrow = length(p2$data$label) - dim(merg)[1], ncol = dim(merg)[2]))
# mat
# colnames(mat) <- colnames(merg)
# merg_long <- merg %>%
#   bind_rows(mat)
# dim(merg_long)
# merg_long$colr <- "black"

p2 <- p %<+% merg

p2$data$is_active_colr
p2$data$class[!p2$data$isTip] <- "unknown"
p2$data$order[!p2$data$isTip] <- "unknown"
p2$data$sz[as.numeric(p2$data$label) > 0.9] <- 0.2

# Set the color palette
#pal<-colorRampPalette(brewer.pal(12,"Paired"))(length(unique(p2$data$class)))
pal <- brewer.pal(12,"Paired")
pal2 <- c(pal) #rep("black", 9))

# pal2[pal2 == '#6A3D9A'] <- "gray"
# pal2[pal2 == '#FDBF6F'] <- "#6A3D9A"
pal2[pal2 == "#B15928"] <- "gray"
pal2[pal2 == '#FFFF99'] <- "#B15928"
pal2[pal2 == '#FF7F00'] <- "goldenrod"


# Active

show_col(pal2, labels = TRUE, borders = NULL, cex_label = 1)

p2$data$is_active_colr[grepl("Photorhabdus_luminescens|Azospirillum_brasilense|C6W5A4_DYAFD", p2$data$label[p2$data$isTip])][1:3] <- "white"
p2$data$is_active_colr[grepl("Actinoplanes|Bacteroides|Bacillus", p2$data$label[p2$data$isTip])][1] <- "maroon"
p2$data$is_active_sz[grepl("Actinoplanes|Bacteroides|Bacillus", p2$data$label[p2$data$isTip])][1] <- 1

pdf("output/57_FastTree_JGI_OleC_circular_20191009.pdf", width = 6, height = 6)
par(mar=c(0.001,0.001,0.001,0.001))
ptree <- p2 +
  aes(color = class) +
  scale_color_manual(values = pal2) +
  ggplot2::xlim(-1, NA) +
  # geom_tippoint(color = p2$data$is_active_colr[p2$data$isTip], size = p2$data$is_active_sz[p2$data$isTip]) +
  geom_tiplab2(aes(label=genus), size = 2, hjust = -0.2) +
  geom_nodepoint(size = p2$data$sz[!p2$data$isTip], color = "black")
#   geom_tiplab2(aes(label=label), size=12, hjust=-.2)
ptree
#scale_color_manual(pal)
# geom_text2(aes(subset=!isTip,label=label),size=14,hjust=-.2)
dev.off()

lvls <- levels(as.factor(p2$data$class))[c(1:11)]
lvls
# lvls[10] <- "Green nonsulfur bacteria"
colrs <- pal2[c(1:11)]
colrs
pdf("56_OleC_legend.pdf", width = 5, height = 5)
plot.new()
legend("bottom", ncol = 1, legend =  lvls, fill = colrs,
       border=FALSE, bty="n")
dev.off()

