# Install packages
pacman::p_load(DECIPHER, viridis, ggthemes, treeio, tidyverse, seqinr, bgafun,
               RColorBrewer, ape, phangorn, Biostrings, data.table, purrr, dplyr, ggtree, bgafun)

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/beta_lactone_stereochemistry/")

# Read in the sequences
cis <- readAAStringSet("data/cis_beta_lactone_enzymes.fasta")
names(cis) <- paste0(names(cis), "_cis")
trans <- readAAStringSet("data/trans_beta_lactone_enzymes.fasta")
names(trans) <- paste0(names(trans), "_trans")
names(trans)
# Remove Saccharothrix 
trans <- trans[-grep("Saccharothrix", names(trans))]
# Combine and align with HMM
comb <- c(cis, trans)
names(comb) <- gsub(" ", "_", names(comb))
grep("Nocardia_brasiliensis", names(comb))
vmatchPattern("PRAI", comb[10]) # 172-173

# writeXStringSet(comb, "data/combined_cis_trans_D_enzymes.fasta")

# hmmalign --trim -o data/combined_cis_trans_D_enzymes.sto ~/Downloads/Epimerase.hmm data/combined_cis_trans_D_enzymes.fasta

# Read in the alignment
# rdaln <- read.alignment("data/HMMAligned_OleD.fasta", format = "fasta")
rdaln <- read.alignment("data/HMMAligned_no_Saccharothrix.fasta", format = "fasta")
rdaln$nam
aa.al <- readAAStringSet("data/HMMAligned_OleD.fasta")
aa.al <- aa.al[-grep("Saccharothrix", names(aa.al))]
writeXStringSet(aa.al, "data/HMMAligned_no_Saccharothrix.fasta", format = "fasta")
aa.al <- AAStringSet(toupper(aa.al))



# Convert alignemnt into a binary presence-abscence matrix 
amino <- convert_aln_amino(rdaln)

# Define groups: "positive" and "negative"
grps <- rownames(amino)
grps[grep("_trans", grps)] <- "Trans"
grps[grep("_cis", grps)] <- "Cis"
grps <- as.factor(grps)
table(grps)

# Remove gaps
amino.gapless <- remove_gaps_groups(x = amino, z = grps) 

# Run correspondance analysis (CA)
ca.aap <- bga(t(amino.gapless + 1), grps)

# Identify important residues
top_res <- top_residues_2_groups(ca.aap)
names(top_res) <- gsub("X", "", names(top_res))
top_res

# Conserved amino acid profiles
profiles <- create_profile_strings(amino.gapless, grps)
keep <- profiles[, colnames(profiles) %in% names(top_res)]
keep
table(grps)
inds <- which(duplicated(substr(colnames(keep), 1, nchar(colnames(keep))-1)))
inds
before.inds <- inds-1
before.inds
comb <- sort(c(inds, before.inds))
comb
keep[,comb]
write_csv(data.frame(keep[,comb]), "output/table_conserved_positions.csv")

# 8 candidate positions for mutagenesis
# 127, 212, 220, 221, 234, 262, 263, 273
aas3 <- colnames(keep[,comb])
aas3
aas.ind <- substr(aas3, 1, nchar(aas3)-1)

gp <- list()
aas.ind <-as.numeric(aas.ind)
aas.ind

rleg.pos <- aa.al[grep("Nocardia_brasiliensis", names(aa.al))]
rleg.pos
for(i in 1:length(aas.ind)) gp[[i]] <- substr(rleg.pos, start = aas.ind[i] - 2, stop = aas.ind[i] + 2)

unlist(gp)

patterns<-unlist(gp)
patterns
patt <- gsub("-","", patterns)

leg <- trans[grep("brasiliensis", names(trans))]
leg

tmp<-list()
for(i in 1:length(patt)){
  tmp[[i]]<-vmatchPattern(patt[i],leg) 
}

pattmatches <- unlist(tmp, use.names = T)

m <- t(data.frame(pattmatches))
m <- m[complete.cases(m)]
m
m2<-as.numeric(m[as.numeric(m)!=1])
m2

strpos1 <- m2[m2 > 10]
strpos1

strposA <- strpos1[seq(1, length(strpos1), 2)]
strposB <- strpos1[seq(2, length(strpos1), 2)]
length(strposA)
#mean(strposA,strposB)
mns <- list()
strposA
strposB


for(i in 1:length(strposA)) mns[[i]] <- (strposA[i] + strposB[i])/2

mnsl <- unlist(mns)
mnsl
dtinds <- data.frame(mnsl)
# dtinds$pos <- aas
mnsl
paste0(unique(mnsl), collapse = '+')
mnsl

write.csv(data.frame(mnsl), "data/Nocardia_brasiliensis_positions.csv")

rleg.pos <- aa.al[grep("Xanthomonas_campestris", names(aa.al))]
rleg.pos
for(i in 1:length(aas.ind)) gp[[i]] <- substr(rleg.pos, start = aas.ind[i] - 2, stop = aas.ind[i] + 2)

unlist(gp)

patterns<-unlist(gp)
patterns
patt <- gsub("-","", patterns)

leg <- cis[grep("campestris", names(cis))]
leg

tmp<-list()
for(i in 1:length(patt)){
  tmp[[i]]<-vmatchPattern(patt[i],leg) 
}

pattmatches <- unlist(tmp, use.names = T)

m <- t(data.frame(pattmatches))
m <- m[complete.cases(m)]
m
m2<-as.numeric(m[as.numeric(m)!=1])
m2

strpos1 <- m2[m2 > 10]
strpos1

strposA <- strpos1[seq(1, length(strpos1), 2)]
strposB <- strpos1[seq(2, length(strpos1), 2)]
length(strposA)
#mean(strposA,strposB)
mns <- list()
strposA
strposB


for(i in 1:length(strposA)) mns[[i]] <- (strposA[i] + strposB[i])/2

mnsl <- unlist(mns)
mnsl
dtinds <- data.frame(mnsl)
# dtinds$pos <- aas
mnsl
unique(mnsl)
paste0(unique(mnsl), collapse = '+')

# colnames(dtinds) <- c("res_pos", "identity")
write.csv(data.frame(mnsl), "data/Xanthomonas_campestris_positions.csv")
