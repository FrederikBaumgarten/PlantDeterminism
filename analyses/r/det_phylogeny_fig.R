# =============================================================================
# Phylogeny of tree-species growth determinacy
# -----------------------------------------------------------------------------
# Generates a phylogeny of the study species with the determinacy trait
# (indeterminate / intermediate / determinate) mapped onto the tree, and
# writes the pruned tree and the figure (phyloIntColor.pdf).
#
# Associated publication:
#   Baumgarten et al. (2026). (In)determinacy in woody plants: limits and opportunities for timing growth in a changing climate. Ecology Letters. [DOI]
#
# Data/code author: Frederik Martin Baumgarten
#   ORCID 0000-0002-8284-8384
#   University of Basel, Department of Environmental Sciences
#   frederik.baumgarten@unibas.ch
# 
# Reproducibility notes:
#   * Paths are RELATIVE to this script's location, so it runs without editing
#     as long as all files are in the same folder (see "Required files" below).
#   * Please run in a CLEAN R session. This script does NOT call
#     rm(list = ls()); just start a fresh session before sourcing it.
#   * R and package versions used in the original analysis are listed in
#     read_me.txt.
#
# Required files in the SAME folder as this script:
#   Determinism_data.csv   - species / determinacy data (included in this archive)
#   ALLMB.tre              - Smith & Brown (2018) seed-plant megatree.
#                            NOT included here (large file, separate licence).
#                            Download "ALLMB.tre v1.0" from
#                            https://doi.org/10.6084/m9.figshare.9747638
#                            and place it next to this script.
#                            Cite: Smith, S.A. & Brown, J.W. (2018) Constructing
#                            a broadly inclusive seed plant phylogeny. American
#                            Journal of Botany 105, 302-314.
#
# Outputs (written to the same folder):
#   determinacyPhylogeny.tre
#   phyloIntColor.pdf
# =============================================================================
options(stringsAsFactors = FALSE)  # read words as characters, not factors

#library(tidyverse)
library(stringr)
library(ape)
library(phytools)
library(geiger)
library(pez)
library(caper)
library(phangorn)

# read in the data

det <- read.csv("Determinism_data.csv")
det$sp.name <- paste(det$genus, det$species, sep = "_")

### combine with trait of determinacy
# 1 is indeterminate
# 2 is mixed
# 3 is determinate


#make a new column with this trait
det$trait_det <- NA

#if there is a "x" in column "Determinate_all_preformed"
det[!is.na(det$Determinate_all_preformed),"trait_det"] <- 3 #determinate assiged with 3
det[!is.na(det$Intermediate_preformed_neogrown),"trait_det"] <- 2 #intermediate assiged with 2
det[!is.na(det$Indeterminate_all_neogrown),"trait_det"] <- 1 #indeterminate assiged with 1

unique(det$trait_det)

## species to keep:
#remove species with NA in trait_det
det<-det[!is.na(det$trait_det),]

#remove all rows containing "spp" in "species"
det <- det[!grepl("spp|ssp", det$species),]

#calculate the mean of the trait_det for each species
df_trait <- aggregate(det["trait_det"], det[c("sp.name")], function(x) {
  mean(x, na.rm=TRUE)
})


unique(df_trait$sp.name)
sps.list <- sort(unique(df_trait$sp.name))
genus.list=sort(unique(df_trait$genus))


## load phylo (from Smith and Brown 2019)
phy.plants<-read.tree("ALLMB.tre")

## getting a list of genera in S&B's phylo
phy.genera<-unlist(
  lapply(strsplit(phy.plants$tip.label, "_"),function(x){return(x[1])})
)
phy.genera.uniq<-sort(unique(phy.genera))

## how many genera of det.csv are in the phylogeny?
phenosp.genus.inphylo<-genus.list[which(genus.list%in%phy.genera.uniq)] #182 out of our 185

## first prune the phylogeny to include$ only these genera
# phy.genera.egret<-drop.tip(phy.plants,
#                              which(!phy.genera %in% phenosp.genus.inphylo)) #34940 tips
# length(phy.genera.egret$tip.label)
tree <- drop.tip(phy.plants, which(!phy.plants$tip.label %in% sps.list))

length(tree$tip.label)
sort(tree$tip.label)

#basic plotting
plot(tree, type="fan")

write.tree(tree,"determinacyPhylogeny.tre")




namesphy <- tree$tip.label
tree$root.edge <- 0

is.rooted(tree)
tree$node.label<-NULL

dataPhy = comparative.data(tree, df_trait, names.col = "sp.name", na.omit = T,
                           vcv = T, warn.dropped = T)

phyloplot = dataPhy$phy
x = dataPhy$data$trait_det
names(x)=dataPhy$phy$tip.label

study <- contMap(tree, x, plot = T)

slopeCol <- setMap(study, colors=c("blue","yellow","red"))
h<-max(nodeHeights(slopeCol$tree))

pdf("phyloIntColor.pdf", height = 13, width = 10)
plot(slopeCol,legend = T, lwd=3, ylim=c(1-0.09*(Ntip(slopeCol$tree)),Ntip(slopeCol$tree)))

dev.off()

pdf("phyloIntColor.pdf", height = 10.5, width = 8.3)

# Plot without the legend
#plot(slopeCol, legend = FALSE, lwd = 3, ylim = c(1 - 0.09 * Ntip(slopeCol$tree), Ntip(slopeCol$tree)))
plot(slopeCol, legend = FALSE, lwd = 3, ylim = c(0, Ntip(slopeCol$tree)))

# Add the color legend to the right
add.color.bar(
  150,                          # Length of the color bar (you can adjust this)
  slopeCol$cols,              # Colors used
  title = "Trait value",      # Legend title (adjust to your variable)
  lims = slopeCol$lims,       # Value limits
  digits = 2,                 # Number of digits to show
  prompt = FALSE,
  x = 10, y = 0,              # Position – tweak y value as needed
  subtitle = "",              # Optional subtitle
  lwd = 3,                    # Match the tree line width
  fsize = 1                   # Font size for legend text
)
dev.off()
