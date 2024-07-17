# Started Nov 27, 2023 by Deirdre L.
# Aim of this code is to get a phylogeny for all of the spp we have germination data scraped

if(length(grep("deirdreloughnan", getwd()) > 0)) {
  setwd("~/Documents/github/egret")
} else if(length(grep("Lizzie", getwd()) > 0)) {
  setwd("~/Documents/git/projects/others/deirdre/egret")
} else if(length(grep("frederik", getwd())) > 0){
  setwd("/Users/frederik/github/PlantDeterminism/analyses")
}

#library(tidyverse)
library(stringr)
library(ape)
library(phytools)
library(geiger)
library(pez)
library(caper)
library(phangorn)
library(readxl)

rm(list = ls()) # Clear whatever is already in R's memory
options(stringsAsFactors=FALSE)# Make sure words are read in as characters rather than factors

# read in the data

det <- read_excel("input/Determinism.xlsx")
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
phy.plants<-read.tree("./input/ALLMB.tre")

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

write.tree(tree,"./output/determinacyPhylogeny.tre")




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

pdf("./output/phyloIntColor.pdf", height = 20, width = 7)
plot(slopeCol,legend = F, lwd=3, ylim=c(1-0.09*(Ntip(slopeCol$tree)),Ntip(slopeCol$tree)))

dev.off()




