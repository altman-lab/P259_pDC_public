library(tidyverse)
library(biomaRt)
source("https://raw.githubusercontent.com/kdillmcfarland/R_bioinformatic_scripts/master/STRING_network_fxn.R")
set.seed(8434)

#### DEG ####
fdr.cut <- 0.1
#Define change with virus
DEG1.v <- read_csv("results/gene_level/P259.1_gene_pval.csv") %>% 
  filter(group %in% c("none_HRV - none_none", "EOS.supp_HRV - EOS.supp_none") &
           adj.P.Val <= fdr.cut)
#Define change with treament
DEG1.t <- read_csv("results/gene_level/P259.1_gene_pval.csv") %>% 
  filter(group %in% c("EOS.supp_none - none_none", "EOS.supp_HRV - none_HRV") &
           adj.P.Val <= fdr.cut)

#Define change with virus
DEG2.v <- read_csv("results/gene_level/P259.2_gene_pval.csv") %>% 
  filter(group %in% c("none_HRV - none_none", "AntiIL5_HRV - AntiIL5_none") &
           adj.P.Val <= fdr.cut)
#Define change with treament
DEG2.t <- read_csv("results/gene_level/P259.2_gene_pval.csv") %>% 
  filter(group %in% c("AntiIL5_none - none_none", "AntiIL5_HRV - none_HRV") &
           adj.P.Val <= fdr.cut)

#Format in list
DEG <- c(intersect(DEG1.v$hgnc_symbol, DEG1.t$hgnc_symbol),
         intersect(DEG2.v$hgnc_symbol, DEG2.t$hgnc_symbol))
##replace HGNC symbol with one in STRING
DEG <- gsub("MARCHF2","MARCH2",DEG)
DEG <- gsub("H2BC21","HIST2H2BE",DEG)
DEG <- gsub("H4C14","HIST2H4A",DEG)

#### Enrichment ####
H <- read_csv("results/enrichment/enrich_DEG_H.csv") %>% 
  filter(p.adjust <= 0.05 & size.overlap.term >= 2)
C2 <- read_csv("results/enrichment/enrich_DEG_C2_CP.csv") %>% 
  filter(p.adjust <= 0.05 & size.overlap.term >= 2)

enrich <- bind_rows(H,C2) %>% 
  mutate(Description = recode_factor(factor(Description), 
           "HALLMARK_INTERFERON_ALPHA_RESPONSE"="IFNA response", 
           "HALLMARK_INTERFERON_GAMMA_RESPONSE"="IFNG response", 
           "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING"="IFNA/B signaling", 
           "REACTOME_INTERFERON_SIGNALING"="IFN signaling", 
           "REACTOME_GPCR_LIGAND_BINDING"="GPCR ligand binding", 
           "REACTOME_CLASS_A_1_RHODOPSIN_LIKE_RECEPTORS"="Class A1 rhodopsin-like receptors", 
           "REACTOME_PEPTIDE_LIGAND_BINDING_RECEPTORS"="Peptide ligand binding receptors", 
           "REACTOME_CHEMOKINE_RECEPTORS_BIND_CHEMOKINES"="Chemokine receptors bind chemokines", 
           "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION"="Cytokine-cytokine reception interaction")) %>% 
  arrange(Description) %>% 
  #remove duplicate terms
  filter(Description != "IFNA/B signaling")

col.vec <- c('#bdd7e7','#6baed6','#3182bd',#'#08519c', #IFN blue
             '#a1d99b','#31a354', #GPCR green
             '#fdbe85','#fd8d3c','#d94701', #Cytokine/chemokine receptors orange
             'grey')        

# DEG2 <- c("CXCL10","CXCL11","GMPR","IFIT1","IFIT2","MX2","PARP9","RSAD2","IFNA2","IFNG")
# DEG2 <- c(DEG,"IFNA2","IFNG")

#### Network ####
string.plot(genes=DEG, version="11", score_threshold=400,
            layout='fr',
            enrichment=enrich, size.overlap.term=2, p.adjust=0.05,
            ID="SYMBOLs",# discard = "edge.keep.enrich",
            colors=col.vec, 
            outdir="publication/fig/", basename="FigSX.DEG.STRING.",
            width=12, height=8)

# string.plot(genes=DEG2, version="11", score_threshold=400,
#             layout='fr',
#             enrichment=enrich, size.overlap.term=2, p.adjust=0.05,
#             ID="SYMBOLs",
#             colors=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99',
#                      '#e31a1c','#fdbf6f','#ff7f00','#cab2d6',#'#6a3d9a',
#                      'grey'), 
#             outdir="publication/fig/", basename="FigSX.DEG.STRING.",
#             width=12, height=8)
