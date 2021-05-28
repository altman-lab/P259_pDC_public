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

#### GSEA terms ####
H <- read_csv("results/enrichment/enrich_DEG_H.csv") %>% 
  filter(Description %in% c("HALLMARK_INTERFERON_ALPHA_RESPONSE",
                            "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                            "HALLMARK_INFLAMMATORY_RESPONSE",
                            "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"))  %>% 
  mutate(Description = recode_factor(factor(Description), 
           "HALLMARK_INTERFERON_ALPHA_RESPONSE"="IFNA response", 
           "HALLMARK_INTERFERON_GAMMA_RESPONSE"="IFNG response",
           "HALLMARK_INFLAMMATORY_RESPONSE"="Inflammatory response",
           "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"=
             "Epithelial mesenchymal transition")) %>% 
  arrange(Description) 

col.vec <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','grey')        

#### Network ####
string.plot(genes=DEG, version="11", score_threshold=400,
            layout='kk',
            enrichment=H, size.overlap.term=0, p.adjust=1,
            ID="SYMBOLs",
            colors=col.vec, 
            outdir="publication/fig/", basename="FigE3.",
            width=12, height=8)
