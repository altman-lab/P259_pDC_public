#Create directories of outputs
dir.create("publication/table/", showWarnings = FALSE, recursive = TRUE)

#Load packages used across all data
library(tidyverse)
library(limma)
#Set seed
set.seed(4389)


# Tables
## RNA-seq library and patient metadata
attach("data_clean/P259_pDC_clean.RData")

dat.pDC.voom$targets %>% 
  select(libID, everything()) %>% 
  arrange(libID) %>% 
  write_csv(., "publication/table/TableS.metadata.csv")

### Patient summary demographics
#Age
dat.pDC.voom$targets %>% 
  distinct(experiment, donorID, Age) %>% 
  group_by(experiment) %>% 
  summarise(.groups = "keep",
            n = n(),
            mean.age = mean(Age, na.rm=TRUE),
            sd.age = sd(Age, na.rm=TRUE)/sqrt(n))
#Sex
dat.pDC.voom$targets %>% 
  distinct(experiment, donorID, Sex) %>% 
  count(experiment, Sex) %>% 
  group_by(experiment) %>% 
  mutate(pct = n/sum(n))

## RNA-seq normalized log2 counts
as.data.frame(dat.pDC.voom$E) %>% 
  rownames_to_column("geneName") %>% 
  #Save
  write_csv(., "publication/table/TableS.norm.log2.counts.csv")

## All contrast model results
#Get all model results
read_csv("results/gene_level/P259.1_gene_pval.csv") %>% 
  #Add experiment variable
  mutate(experiment = "P259_1") %>% 
  bind_rows(read_csv("results/gene_level/P259.2_gene_pval.csv")) %>% 
  #Fill in experiment variable
  mutate(experiment = ifelse(is.na(experiment), "P259_2", experiment)) %>% 
  #Keep only contrasts model 
  filter(model == "contrasts") %>% 
  #Reorder and rename variables
  select(experiment, model, geneName, hgnc_symbol, group,
         AveExpr:group) %>% 
  rename(contrast = group) %>% 
  arrange(experiment, geneName, contrast) %>% 
  #Save
  write_csv(., "publication/table/TableS.models.csv")

## GSEA results
read_csv("results/GSEA_FoldChange/h_GSEA.result.csv") %>% 
  rename(contrast = group) %>% 
  #Add experiment variable
  mutate(experiment = ifelse(grepl(".1", contrast), "P259_1", "P259_2")) %>%
  #Set gene set variable
  mutate(gene.set = "HALLMARK",
         pathway = gsub("HALLMARK_", "", pathway)) %>% 
  #Format contrast to match model results
  mutate(contrast = gsub(".1", "", contrast),
         contrast = gsub(".2", "", contrast),
         contrast = gsub("[.]n"," - n", contrast),
         contrast = gsub("[.]A"," - A", contrast),
         contrast = gsub("[.]E"," - E", contrast)) %>% 
  #Reorder variables
  select(experiment, gene.set, pathway, contrast, everything()) %>% 
  arrange(experiment, pathway, contrast) %>% 
  #Save
  write_csv(., "publication/table/TableS.GSEA.csv") 