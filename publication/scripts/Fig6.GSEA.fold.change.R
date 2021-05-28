#### Setup ####
library(tidyverse)
library(cowplot)
library(limma)
library(ggrepel)
set.seed(4389)

#### Data ####
#List genes in each Hallmark term of interest
myGO <- fgsea::gmtPathways("data_clean/Broad_gmt/h.all.v7.4.symbols.gmt")
GO.df <- plyr::ldply(myGO, rbind) %>% 
  dplyr::rename(term = `.id`) %>% 
  pivot_longer(-term, values_to="hgnc_symbol")

#Expression data
load("data_clean/P259_pDC_clean.RData")

#### Format data ####
#Select genes of interest and format
## P259.2
h.combo.plot_2 <- as.data.frame(dat.pDC.voom_2$E) %>% 
  rownames_to_column("geneName") %>% 
  #Get HGNC symbols
  left_join(dat.pDC.voom_2$genes) %>% 
  dplyr::select(hgnc_symbol, everything()) %>% 
  #Filter to genes in terms
  inner_join(GO.df) %>% 
  dplyr::select(term, hgnc_symbol, geneName, everything(), 
                -c(`Previous symbols`:gene_biotype), -name) %>% 
  #longer format
  pivot_longer(-c(term:geneName), names_to = "libID", 
               values_to = "expression") %>% 
  #Add metadata
  mutate(experiment = "P259.2") %>% 
  left_join(dplyr::select(dat.pDC.voom_2$targets, donorID, libID, 
                          IL5, virus.detail))

## P259.1
h.combo.plot_1 <- as.data.frame(dat.pDC.voom_1$E) %>% 
  rownames_to_column("geneName") %>% 
  #Get HGNC symbols
  left_join(dat.pDC.voom_1$genes) %>% 
  dplyr::select(hgnc_symbol, everything()) %>% 
  #Filter to genes in terms
  inner_join(GO.df) %>% 
  dplyr::select(term, hgnc_symbol, geneName, everything(), 
                -c(`Previous symbols`:gene_biotype), -name) %>% 
  #longer format
  pivot_longer(-c(term:geneName), names_to = "libID", 
               values_to = "expression") %>% 
  #Add metadata
  mutate(experiment = "P259.1") %>% 
  left_join(dplyr::select(dat.pDC.voom_1$targets, donorID,
                          libID, IL5, virus.detail))

#Combine data and calculate fold change
h.combo.plot <- bind_rows(h.combo.plot_1, h.combo.plot_2) %>% 
  mutate(group = paste(virus.detail, IL5, sep=".")) %>% 
  dplyr::select(-IL5, -virus.detail, -libID) %>% 
  #Averages
  group_by(term, hgnc_symbol, experiment, group) %>% 
  dplyr::summarise(meanE = mean(expression, na.rm=TRUE)) %>% 
  #Fold change
  pivot_wider(names_from = group, values_from = meanE) %>% 
  rowwise() %>% 
  mutate(
    #Virus in untreated
    virus.untreat_old = oldHRV.none-none.none,
    virus.untreat_new = newHRV.none-none.none,
    #Virus in treated
    virus.treat_old1 = oldHRV.EOS.supp-none.EOS.supp,
    virus.treat_old2 = oldHRV.AntiIL5-none.AntiIL5,
    virus.treat_new = newHRV.AntiIL5-none.AntiIL5,
    #Treatment in media
    treat.media_1 = none.EOS.supp-none.none,
    treat.media_2 = none.AntiIL5-none.none,
    #Treatment in virus
    treat.virus_old1 = oldHRV.EOS.supp-oldHRV.none,
    treat.virus_old2 = oldHRV.AntiIL5-oldHRV.none,
    treat.virus_new = newHRV.AntiIL5-newHRV.none) %>% 
  ungroup() %>% 
  
  dplyr::select(term:experiment, virus.untreat_old:treat.virus_new) %>% 
  pivot_longer(virus.untreat_old:treat.virus_new) %>% 
  drop_na(value) %>% 
  
  #Average old and new HRV values
  separate(name, into=c("group"), sep="_") %>% 
  group_by(term, hgnc_symbol, experiment, group) %>% 
  dplyr::summarise(value = mean(value, na.rm=TRUE)) %>% 
  ungroup() %>% 
  #Beautify for plotting
  mutate(term = gsub("HALLMARK_","", term)) %>% 
  mutate(col.group = ifelse(value <0, "down", "up")) %>% 
  mutate(group = factor(group, levels=c("treat.media","treat.virus",
                                        "virus.untreat","virus.treat")))

#### Plot parameters ####
#Set jitter
pos <- position_jitter(width = 0.3, seed = 589, height = 0)

#### List leading edge genes to label ####
gsea.LE <- read_csv("results/GSEA/h_GSEA.result.csv") %>% 
  filter(pathway %in% c("HALLMARK_INTERFERON_ALPHA_RESPONSE",
                        "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                        "HALLMARK_INFLAMMATORY_RESPONSE",
                        "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")) %>% 
  #Add experi group
  mutate(experi = c(rep("AntiIL5",nrow(.)/2),rep("EOSsup",nrow(.)/2))) %>% 
  dplyr::select(experi, group, pathway, fgsea.leadingEdge) %>% 
  #Add FC group
  mutate(group2 = ifelse(experi=="AntiIL5" & 
                           group %in% c("none_HRVnone_none", "AntiIL5_HRVAntiIL5_none"), 
                         "AntiIL5.RV",
                         ifelse(experi=="EOSsup" & 
                                  group %in% c("none_HRVnone_none", "EOSsupp_HRVEOSsupp_none"), 
                                "EOSsup.RV",
                                ifelse(group %in% c("AntiIL5_nonenone_none", "AntiIL5_HRVnone_HRV"),
                                       "AntiIL5",
                                       ifelse(group %in% c("EOSsupp_nonenone_none",
                                                           "EOSsupp_HRVnone_HRV"),
                                              "EOSsup", NA))))) %>% 
  group_by(group2, pathway) %>% 
  #group_by(experi, group, pathway) %>% 
  dplyr::summarise(genes = paste(fgsea.leadingEdge, collapse=";")) %>% 
  separate(genes, into = as.character(c(1:500)), sep=";") %>% 
  pivot_longer(as.character(c(1:500)), values_to = "gene") %>% 
  drop_na(gene) %>% 
  mutate(pathway = gsub("HALLMARK_","", pathway))

#Significant Hallmark terms (in Fig 6)
term.OI.ls <- c("INTERFERON_ALPHA_RESPONSE",
                "INTERFERON_GAMMA_RESPONSE",
                "INFLAMMATORY_RESPONSE",
                "EPITHELIAL_MESENCHYMAL_TRANSITION")

to.label <- list()

for(term.OI in term.OI.ls){
  
  temp <- gsea.LE %>% 
    filter(pathway == term.OI)
  
  overlap.vec1 <- intersect(unique(temp[temp$group2=="EOSsup",]$gene),
                            unique(temp[temp$group2=="EOSsup.RV",]$gene))
  overlap.vec2 <- intersect(unique(temp[temp$group2=="AntiIL5",]$gene),
                           unique(temp[temp$group2=="AntiIL5.RV",]$gene))
  
  #Save to list object 
  to.label[[paste(term.OI, "P259.1", sep="_")]] <- overlap.vec1
  to.label[[paste(term.OI, "P259.2", sep="_")]] <- overlap.vec2
}

## Save for use in Fig8
save(to.label,  file = "publication/fig/to.label.LE.RData")
           
#Convert genes to list to df
to.label.df <- plyr::ldply(to.label, data.frame, .id="term") %>% 
  mutate(to.label = "y") %>% 
  dplyr::rename('hgnc_symbol'=`X..i..`) %>% 
  separate(term, into=c("term","experiment"), sep="_(?=[^_]+$)")

#### DEGs to label ####
fdr.cut <- 0.1
DEG1.v <- read_csv("results/gene_level/P259.1_gene_pval.csv") %>% 
  filter(group %in% c("none_HRV - none_none", "EOS.supp_HRV - EOS.supp_none") &
           adj.P.Val <= fdr.cut)
DEG1.t <- read_csv("results/gene_level/P259.1_gene_pval.csv") %>% 
  filter(group %in% c("EOS.supp_none - none_none", "EOS.supp_HRV - none_HRV") &
           adj.P.Val <= fdr.cut)
DEG1 <- intersect(DEG1.v$hgnc_symbol, DEG1.t$hgnc_symbol)

DEG2.v <- read_csv("results/gene_level/P259.2_gene_pval.csv") %>% 
  filter(group %in% c("none_HRV - none_none", "AntiIL5_HRV - AntiIL5_none") &
           adj.P.Val <= fdr.cut)
DEG2.t <- read_csv("results/gene_level/P259.2_gene_pval.csv") %>% 
  filter(group %in% c("AntiIL5_none - none_none", "AntiIL5_HRV - none_HRV") &
           adj.P.Val <= fdr.cut)
DEG2 <- intersect(DEG2.v$hgnc_symbol, DEG2.t$hgnc_symbol)

## Save for use in Fig8
save(DEG1, DEG2,  file = "publication/fig/to.label.DEG.RData")

#add to data
h.combo.plot.lab <- h.combo.plot %>% 
 mutate(to.label = ifelse(experiment == "P259.1" & hgnc_symbol %in% c(DEG1,DEG2) &
                          hgnc_symbol %in% to.label.df[to.label.df$experiment=="P259.1",]$hgnc_symbol,
                          "y",
                          ifelse(experiment == "P259.2" & hgnc_symbol %in% c(DEG1,DEG2) &
                          hgnc_symbol %in% to.label.df[to.label.df$experiment=="P259.2",]$hgnc_symbol,
                                 "y", NA)))

#### All plot data ####
FC.plot.dat <- h.combo.plot.lab %>% 
  mutate(group = paste(experiment, group, sep="_")) %>% 
  mutate(group1 = recode_factor(factor(group),
                                "P259.2_virus.treat"='"+ Anti-IL-5/5R"*alpha',
                                "P259.2_virus.untreat"='"- Anti-IL-5/5R"*alpha',
                                "P259.1_virus.treat"='"+ EOS sup-"',
                                "P259.1_virus.untreat"='"- EOS sup"',
                                
                                "P259.1_treat.virus"="+ RV",
                                "P259.2_treat.virus"="+ RV",
                                "P259.1_treat.media"="- RV",
                                "P259.2_treat.media"="- RV"),
         group2 = recode_factor(factor(group),
            "P259.2_virus.treat"='italic("Ex vivo")~"RV-infected vs media"',
            "P259.2_virus.untreat"='italic("Ex vivo")~"RV-infected vs media"',
            "P259.1_virus.treat"='italic("Ex vivo")~"RV-infected vs media"',
            "P259.1_virus.untreat"='italic("Ex vivo")~"RV-infected vs media"',
                                
            "P259.1_treat.virus"='italic("Ex vivo")~"EOS supernatant vs none"',
            "P259.2_treat.virus"='italic("In vivo")~"Anti-IL-5/5R"*alpha~"vs none"',
            "P259.1_treat.media"='italic("Ex vivo")~"EOS supernatant vs none"',
            "P259.2_treat.media"='italic("In vivo")~"Anti-IL-5/5R"*alpha~"vs none"')) %>% 
  #Reorder terms
  mutate(term = gsub("HALLMARK_", "", term)) %>% 
  #Reorder terms
  mutate(term.ord = recode_factor(factor(term),
                                     "INTERFERON_ALPHA_RESPONSE"="IFNA response",
                                     "INTERFERON_GAMMA_RESPONSE"="IFNG response",
                                     "INFLAMMATORY_RESPONSE"="Inflammatory response",
                                     "EPITHELIAL_MESENCHYMAL_TRANSITION"="Epithelial mesenchymal transition"))

#### Plot A: VIRUS ####
dat.v <- FC.plot.dat %>% 
  filter(group2 == 'italic("Ex vivo")~"RV-infected vs media"' &
           term %in% c("INTERFERON_ALPHA_RESPONSE",
                           "INTERFERON_GAMMA_RESPONSE",
                           "INFLAMMATORY_RESPONSE",
                           "EPITHELIAL_MESENCHYMAL_TRANSITION")) %>% 
  droplevels()

plot.v <- dat.v %>% 
  ggplot(aes(x = group1, y=value)) +
  geom_violin() +
  #Add non-labeled points
  geom_jitter(data=filter(dat.v, is.na(to.label)),
              aes(color=col.group), position = pos) +
  #Add labeled points
  geom_point(data=filter(dat.v, !is.na(to.label)),
             aes(color=col.group), ) +
  #Add labels
  geom_text_repel(data=filter(dat.v, !is.na(to.label)),
                  aes(label=hgnc_symbol), direction="both",
                  nudge_x=-0.4, 
                  show.legend = FALSE, size=3, max.overlaps=100) +
 #Add mean lines
  stat_summary(fun="mean", geom="crossbar")+
  facet_grid(term.ord~group2, scales="free", 
             labeller = labeller(group2=label_parsed)) +
  coord_flip() +
  #Beautify
  theme_bw() +
  labs(x="", y="Gene mean log2 fold change") +
  scale_color_manual(values=c("down"="#6baed6","up"="#ef3b2c")) +
  #theme(strip.text.y = element_blank())+
  theme(legend.position = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background =element_rect(fill="white")) +
  scale_x_discrete(labels = ggplot2:::parse_safe)
#plot.v

#### Plot B/C: EOS supernatant or Anti-IL5 therapy ####
dat.t <- FC.plot.dat %>% 
  filter(group2 != 'italic("Ex vivo")~"RV-infected vs media"' &
           term %in% c("INTERFERON_ALPHA_RESPONSE",
                           "INTERFERON_GAMMA_RESPONSE",
                           "INFLAMMATORY_RESPONSE",
                           "EPITHELIAL_MESENCHYMAL_TRANSITION")) %>% 
  droplevels()

plot.t <- dat.t %>% 
  ggplot(aes(x = group1, y=value)) +
  geom_violin() +
  #Add non-labeled points
  geom_jitter(data=filter(dat.t, is.na(to.label)),
              aes(color=col.group), position = pos) +
  #Add labeled points
  geom_point(data=filter(dat.t, !is.na(to.label)),
             aes(color=col.group), ) +
  #Add labels left
  geom_text_repel(data=filter(dat.t, !is.na(to.label)),
                  aes(label=hgnc_symbol), direction="both",
                  nudge_x=-0.4,
                  show.legend = FALSE, size=3, max.overlaps = 100) +
  #Add mean lines
  stat_summary(fun="mean", geom="crossbar")+
  facet_grid(term.ord~group2, scales="free", 
             labeller = labeller(group2=label_parsed)) +
  coord_flip() +
  #Beautify
  theme_bw() +
  labs(x="", y="Gene mean log2 fold change") +
  scale_color_manual(values=c("down"="#6baed6","up"="#ef3b2c")) +
  theme(legend.position = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background =element_rect(fill="white"))
#plot.t

#### Save ####
fc.plot.all <- plot_grid(plot.v, plot.t, ncol=2, rel_widths = c(0.6,1),
                         labels=c("A","B"))
#fc.plot.all
ggsave("publication/fig/Fig6.GSEA.fold.change.pdf", fc.plot.all, 
       width = 20, height = 15)
