#### Setup ####
library(plyr)
library(tidyverse)
library(fgsea)
library(venn)

#### Data ####
myGO <- fgsea::gmtPathways("data_clean/Broad_gmt/h.all.v7.4.symbols.gmt") %>% 
  #Convert to data frame
  plyr::ldply(., data.frame) %>% 
  dplyr::rename(pathway = `.id`, gene = `X..i..`)

#### Venn ####
# Get list of genes in pathways

venn.ls <- list()

for(term in c("HALLMARK_INTERFERON_ALPHA_RESPONSE",
              "HALLMARK_INTERFERON_GAMMA_RESPONSE",
              "HALLMARK_INFLAMMATORY_RESPONSE",
              "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")){
  
  #Subset to term of interest
  myGO.temp <- myGO %>% 
    filter(pathway == term) %>% 
    distinct(gene) %>% unlist(use.names = FALSE)
  
  #Save to list
  name <- gsub("HALLMARK_", "", term)
  name <- gsub("_", "\n", name)
  
  venn.ls[[name]] <- myGO.temp
}

# Plot venn and save
pdf("publication/fig/FigE4.GSEA.term.venn.pdf",
    width=5.5, height=5)
venn(ilab=FALSE, ilcs=1, sncs=1, box = FALSE,
     x=venn.ls)
dev.off()