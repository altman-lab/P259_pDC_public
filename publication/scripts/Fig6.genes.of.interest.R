#### Setup ####
library(tidyverse)
library(cowplot)
library(limma)
set.seed(4389)

#### Data ####
load("data_clean/P259_pDC_clean.RData")

#### List all genes to plot ####
#Genes to label from Fig7
load("publication/fig/to.label.DEG.RData")
load("publication/fig/to.label.LE.RData")

LE <- unique(unlist(to.label[1:3], use.names=FALSE))
GOI <- intersect(c(DEG1,DEG2),LE)
#Add IFN genes
GOI2 <- c("IFNA2","IFNG")

#### Extract expression data ####
dat.GOI <- as.data.frame(dat.pDC.voom$E) %>% 
  rownames_to_column("geneName") %>% 
  pivot_longer(-geneName, names_to="libID", values_to="expression") %>% 
  left_join(dat.pDC.voom$targets)

#### Add hgnc symbol ####
#Get ref genome
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", 
                   host = "useast.ensembl.org")
#Check genome version
#searchDatasets(mart = ensembl, pattern = "hsapiens")

#Get HGNC symbols
dat.GOI.anno <- getBM(attributes = c('hgnc_symbol', 'ensembl_gene_id'),
                      mart = ensembl) %>%
  #Keep genes in RNAseq data
  filter(ensembl_gene_id %in% unique(dat.GOI$geneName)) %>% 
  inner_join(dat.GOI, by=c("ensembl_gene_id"="geneName")) %>% 
  #Beautify labels
  mutate(x.lab=paste(experiment, virus.detail,IL5, sep="_"),
         x.lab=gsub("oldH|newH", "", x.lab),
         x.lab=recode_factor(factor(x.lab),
                             "P259_1_none_none"='"RV            -\nEOS sup  -"',
                             "P259_2_none_none"='"-\n-"',
                             "P259_1_none_EOS.supp"='"-\n+"',
                             "P259_2_none_AntiIL5"='"-\n+"',
                             "P259_1_RV_none"='"+\n-"',
                             "P259_2_RV_none"='"+\n-"',
                             "P259_1_RV_EOS.supp"='"+\n+"',
                             "P259_2_RV_AntiIL5"='"+  RV\n+  Anti-IL-5/5R"*alpha')) %>% 
  mutate(facet.lab = experiment,
         facet.lab = recode_factor(factor(facet.lab),
                            "P259_1" = '"EOS sup"',
                            "P259_2" = '"Anti-IL-5/5R"*alpha')) %>% 
  arrange(facet.lab) 

#### Format pval data ####
library(ggpubr)
#Pvals for eos experiment
GOI.p1 <- read_csv("results/gene_level/P259.1_gene_pval.csv") %>% 
  filter(group != "(Intercept)" & model == "contrasts") %>% 
  dplyr::select(hgnc_symbol, group, P.Value, AveExpr) %>% 
  mutate(dataset="P259.1")

#Pvals for aniti-IL5 experiment
GOI.p <- read_csv("results/gene_level/P259.2_gene_pval.csv") %>% 
  filter(group != "(Intercept)" & model == "contrasts") %>% 
  dplyr::select(hgnc_symbol, group, P.Value,AveExpr) %>% 
  mutate(dataset="P259.2") %>% 
  #Combine with eos data
  bind_rows(GOI.p1) %>% 
  #Keep genes of interest
  filter(hgnc_symbol %in% c(GOI,GOI2)) %>% 
  #Make symbols for plots
  mutate(symbol = ifelse(P.Value <= 0.001,"***",
                         ifelse(P.Value <= 0.01, "**",
                                ifelse(P.Value <= 0.05, "*", NA)))) %>% 
  filter(!is.na(symbol)) %>% 
  
  #Match group labels to those in plots
  mutate(group.long = paste(dataset, group, sep="_"),
         group1 = recode_factor(factor(group.long),
                        "P259.2_none_HRV - none_none"='"-\n-"',
                        "P259.2_AntiIL5_HRV - AntiIL5_none"='"-\n+"',
                        "P259.2_AntiIL5_none - none_none"='"-\n-"',
                        "P259.2_AntiIL5_HRV - none_HRV"='"+\n-"',
                        "P259.1_none_HRV - none_none"='"RV            -\nEOS sup  -"',
                        "P259.1_EOS.supp_HRV - EOS.supp_none"='"-\n+"',
                        "P259.1_EOS.supp_none - none_none"='"RV            -\nEOS sup  -"',
                        "P259.1_EOS.supp_HRV - none_HRV"='"+\n-"'),
         group2 = recode_factor(factor(group.long),
                          "P259.2_none_HRV - none_none"='"+\n-"',
                          "P259.2_AntiIL5_HRV - AntiIL5_none"='"+  RV\n+  Anti-IL-5/5R"*alpha',
                          "P259.2_AntiIL5_none - none_none"='"-\n+"',
                          "P259.2_AntiIL5_HRV - none_HRV"='"+  RV\n+  Anti-IL-5/5R"*alpha',
                          "P259.1_none_HRV - none_none"='"+\n-"',
                          "P259.1_EOS.supp_HRV - EOS.supp_none"='"+\n+"',
                          "P259.1_EOS.supp_none - none_none"='"-\n+"',
                          "P259.1_EOS.supp_HRV - none_HRV"='"+\n+"'),
         facet.lab = recode_factor(factor(dataset),
                            "P259.1" = '"EOS sup"',
                            "P259.2" = '"Anti-IL-5/5R"*alpha')) %>% 
  dplyr::select(facet.lab, hgnc_symbol, group1, group2, symbol)

#Add y location for pval based on max expression in plot
GOI.pe <- dat.GOI.anno %>% 
  #Max expression per gene and experiment
  group_by(hgnc_symbol, experiment) %>% 
  summarise(max.e = max(expression, na.rm=TRUE)) %>% 
  ungroup() %>% 
  mutate(facet.lab = experiment,
         facet.lab = recode_factor(factor(facet.lab),
                            "P259_1" = '"EOS sup"',
                            "P259_2" = '"Anti-IL-5/5R"*alpha')) %>% 
  #Add to pval data
  right_join(GOI.p) %>% 
  arrange(hgnc_symbol, experiment, group1, group2)

#first entry per gene
#Set y position to 1
first <- GOI.pe %>% 
  group_by(hgnc_symbol, experiment) %>% 
  slice(1) %>% 
  mutate(y.position1 = 1)

#Add first position data back and fill in remaining
GOI.pey <- GOI.pe %>% 
  full_join(first) %>% 
  
  #Fill in positions 2 - N
  group_by(hgnc_symbol, experiment, facet.lab) %>% 
  mutate(y.position2 = lag(y.position1)+1,
         y.position3 = lag(y.position2)+1,
         y.position4 = lag(y.position3)+1) %>% 
  #Collapse positions into 1 column
  mutate(y.position = ifelse(!is.na(y.position1),y.position1,
                             ifelse(!is.na(y.position2),y.position2,
                                    ifelse(!is.na(y.position3),y.position3,
                                           ifelse(!is.na(y.position4),y.position4,NA))))) %>% 
  #Scale to max expression
  mutate(y.position = (2^max.e)*((y.position+2)/11+0.85)) %>% 
  ungroup()

#### Plot GSEA genes ####
plot.GOI1a <- dat.GOI.anno %>% 
  filter(hgnc_symbol %in% sort(GOI)[1:4]) %>% 
  droplevels() %>% 
  
  ggplot(aes(x=x.lab, y=2^expression, color=donorID)) +
  geom_jitter(width=0.1, height=0) +
  stat_summary(fun.data=mean_sdl, 
               fun.args = list(mult=1), 
               geom="errorbar", color="black", width=0.25) +
  stat_summary(fun=mean, geom="errorbar", 
               aes(ymax=..y.., ymin=..y..),
               color="black", width=0.5) +
  facet_grid(hgnc_symbol~facet.lab, scales="free", 
             labeller = labeller(facet.lab=label_parsed)) +
  #Add pval
  stat_pvalue_manual(data=filter(GOI.pey, hgnc_symbol %in% sort(GOI)[1:4]), 
                     label="symbol", xmin="group1", xmax="group2") +
  # #Beautify
  theme_bw() +
  labs(x="", y="Normalized expression") +
  theme(legend.position = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background =element_rect(fill="white"),
        axis.text.x=element_text(hjust=0.95,vjust=-0.5),
        plot.margin = margin(0.2,0.2,0.2,0.2,"cm")) +
  scale_x_discrete(labels = ggplot2:::parse_safe)

plot.GOI1b <- dat.GOI.anno %>% 
  filter(hgnc_symbol %in% sort(GOI)[5:8]) %>% 
  droplevels() %>% 
  
  ggplot(aes(x=x.lab, y=2^expression, color=donorID)) +
  geom_jitter(width=0.1, height=0) +
  stat_summary(fun.data=mean_sdl, 
               fun.args = list(mult=1), 
               geom="errorbar", color="black", width=0.25) +
  stat_summary(fun=mean, geom="errorbar", 
               aes(ymax=..y.., ymin=..y..),
               color="black", width=0.5) +
  facet_grid(hgnc_symbol~facet.lab, scales="free", 
             labeller = labeller(facet.lab=label_parsed)) +
  #Add pval
  stat_pvalue_manual(data=filter(GOI.pey, hgnc_symbol %in% sort(GOI)[5:8]), 
                     label="symbol", xmin="group1", xmax="group2") +
  # #Beautify
  theme_bw() +
  labs(x="", y="Normalized expression") +
  theme(legend.position = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background =element_rect(fill="white"),
        axis.text.x=element_text(hjust=0.95,vjust=-0.5),
        plot.margin = margin(0.2,0.2,0.2,0.2,"cm")) +
  scale_x_discrete(labels = ggplot2:::parse_safe)

#plot.GOI1a
#plot.GOI1b

#### Plot addtl IFN genes
plot.GOI2 <- dat.GOI.anno %>% 
  filter(hgnc_symbol %in% GOI2) %>% 
  droplevels() %>% 
  
  ggplot(aes(x=x.lab, y=2^expression, color=donorID)) +
  geom_jitter(width=0.1, height=0) +
  stat_summary(fun.data=mean_sdl, 
               fun.args = list(mult=1), 
               geom="errorbar", color="black", width=0.25) +
  stat_summary(fun=mean, geom="errorbar", 
               aes(ymax=..y.., ymin=..y..),
               color="black", width=0.5) +
  facet_grid(hgnc_symbol~facet.lab, scales="free", 
             labeller = labeller(facet.lab=label_parsed)) +
  #Add pval
  stat_pvalue_manual(data=filter(GOI.pey, hgnc_symbol %in% GOI2), 
                     label="symbol", xmin="group1", xmax="group2") +
  #Beautify
  theme_bw() +
  labs(x="", y="Normalized expression") +
  theme(legend.position = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background =element_rect(fill="white"),
        axis.text.x=element_text(hjust=1,vjust=-0.5),
        plot.margin = margin(0.2,2,0.2,0.2,"cm")) +
  scale_x_discrete(labels = ggplot2:::parse_safe)
#plot.GOI2

#### Save ####
plot.GOIB <- plot_grid(plot.GOI2, NULL, rel_heights=c(1.3,1),nrow=2)
plot.GOI.all <- plot_grid(plot.GOI1a,NULL,plot.GOI1b,NULL,plot.GOIB, 
                          rel_widths = c(1,0.15,1,0.15,1.25), nrow=1,
                          labels=c("A","","","","B"))

ggsave("publication/fig/Fig6.genes.of.interest.pdf", plot.GOI.all,
       height=8, width=10)
ggsave("publication/fig/Fig6.genes.of.interest.tiff", plot.GOI.all,
       height=8, width=10)
