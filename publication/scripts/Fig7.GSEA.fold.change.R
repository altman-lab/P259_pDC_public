#### Setup ####
library(tidyverse)
library(cowplot)
library(limma)
library(ggrepel)
set.seed(4389)

#### Data ####
#List genes in each Hallmark term of interest
myGO <- fgsea::gmtPathways("data_clean/Broad_gmt/h.all.v7.1.symbols.gmt")
GO.df <- plyr::ldply(myGO, rbind) %>% 
  rename(term = `.id`) %>% 
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
  summarise(meanE = mean(expression, na.rm=TRUE)) %>% 
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
  summarise(value = mean(value, na.rm=TRUE)) %>% 
  ungroup() %>% 
  #Beautify for plotting
  mutate(term = gsub("HALLMARK_","", term),
         term = gsub("_"," ", term)) %>% 
  mutate(col.group = ifelse(value <0, "down", "up")) %>% 
  mutate(group = factor(group, levels=c("treat.media","treat.virus",
                                        "virus.untreat","virus.treat")))

#### Plot parameters ####
#Set jitter
pos <- position_jitter(width = 0.3, seed = 589, height = 0)

#### List genes to label ####
#Top N genes within panel
no.genes <- 2
#Significant Hallmark terms (in Fig 6)
term.OI.ls <- c("INTERFERON ALPHA RESPONSE", "INTERFERON GAMMA RESPONSE",
                "INFLAMMATORY RESPONSE",
                "EPITHELIAL MESENCHYMAL TRANSITION")

to.label <- list()

for(term.OI in term.OI.ls){
  #List genes in congruent direction
  #Virus
  if(term.OI == "EPITHELIAL MESENCHYMAL TRANSITION"){
    genes.v <- h.combo.plot %>% 
      #virus fold changes
      filter(term ==term.OI &
               group %in% c("virus.treat","virus.untreat")) %>% 
      #All values negative (fold change down)
      group_by(term, hgnc_symbol) %>% 
      filter(max(value, na.rm=TRUE)<0) %>% 
      ungroup() %>% 
      distinct(hgnc_symbol) %>% unlist(use.names = FALSE)
  } else{
    genes.v <- h.combo.plot %>% 
      #virus fold changes
      filter(term ==term.OI &
               group %in% c("virus.treat","virus.untreat")) %>% 
      #All values position (fold change up)
      group_by(term, hgnc_symbol) %>% 
      filter(min(value, na.rm=TRUE)>0) %>% 
      ungroup() %>% 
      distinct(hgnc_symbol) %>% unlist(use.names = FALSE)
  }
  
  #EOS: down fold change
  genes.eos <- h.combo.plot %>% 
    filter(term == term.OI & group %in% c("treat.media","treat.virus") &
             experiment == "P259.1") %>% 
    group_by(term, hgnc_symbol) %>% 
    filter(max(value, na.rm=TRUE)<0) %>% 
    ungroup() %>% 
    distinct(hgnc_symbol) %>% unlist(use.names = FALSE)
  
  #AntiIL5: up fold change
  genes.aIL5 <- h.combo.plot %>% 
    filter(term == term.OI & group %in% c("treat.media","treat.virus") &
             experiment == "P259.2") %>% 
    group_by(term, hgnc_symbol) %>% 
    filter(min(value, na.rm=TRUE)>0) %>% 
    ungroup() %>% 
    distinct(hgnc_symbol) %>% unlist(use.names = FALSE)
  
  #List congruent for all genes
  genes.all <- intersect(genes.v, genes.eos)
  genes.all <- intersect(genes.all, genes.aIL5)
  
  #Save to list object if not 0
  if(length(genes.all)>0){
    #Find max FC in each panel
    max.FC <- h.combo.plot %>% 
      filter(term ==term.OI & hgnc_symbol %in% genes.all) %>% 
      mutate(fc.group = ifelse(group %in% c("virus.untreat", "virus.treat"), "virus",
                               ifelse(experiment == "P259.1", "eos","aIL5"))) %>% 
      group_by(term, fc.group, hgnc_symbol) %>% 
      summarise(max.FC = max(abs(value), na.rm = TRUE)) %>% 
      ungroup() %>% 
      #keep top X per group
      group_by(term, fc.group) %>% 
      slice_max(max.FC,n=no.genes) %>%
      ungroup() %>% 
      arrange(hgnc_symbol) 
    
    to.label[[term.OI]] <- max.FC %>% distinct(hgnc_symbol)
  } else{
    to.label[[term.OI]] <- data.frame(hgnc_symbol="none")
  }
}

## Save for use in Fig8
save(to.label,  file = "publication/fig/to.label.RData")

#Convert genes to list to df
to.label.df <- data.table::rbindlist(to.label, idcol = "term") %>% 
  mutate(to.label = "y")

#add to data
h.combo.plot.lab <- h.combo.plot %>% 
  left_join(to.label.df, by = c("term", "hgnc_symbol"))

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
  mutate(term.ord = recode_factor(factor(term), 
            "INTERFERON ALPHA RESPONSE"="INTERFERON\nALPHA\nRESPONSE",
            "INTERFERON GAMMA RESPONSE"="INTERFERON\nGAMMA\nRESPONSE",
            "INFLAMMATORY RESPONSE"="INFLAMMATORY\nRESPONSE",
            "EPITHELIAL MESENCHYMAL TRANSITION"="EPITHELIAL\nMESENCHYMAL\nTRANSITION"))

#### Plot A: VIRUS ####
dat.v <- FC.plot.dat %>% 
  filter(group2 == 'italic("Ex vivo")~"RV-infected vs media"' &
           term.ord %in% c("INTERFERON\nALPHA\nRESPONSE",
                           "INTERFERON\nGAMMA\nRESPONSE",
                           "INFLAMMATORY\nRESPONSE",
                           "EPITHELIAL\nMESENCHYMAL\nTRANSITION")) %>% 
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
                  show.legend = FALSE, size=3) +
 #Add mean lines
  stat_summary(fun="mean", geom="crossbar")+
  facet_grid(term.ord~group2, scales="free", 
             labeller = labeller(group2=label_parsed)) +
  coord_flip() +
  #Beautify
  theme_bw() +
  labs(x="", y="Gene mean log2 fold change") +
  scale_color_manual(values=c("down"="#6baed6","up"="#ef3b2c")) +
  theme(strip.text.y = element_blank())+
  theme(legend.position = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background =element_rect(fill="white")) +
  scale_x_discrete(labels = ggplot2:::parse_safe)
#plot.v

#### Plot B/C: EOS supernatant or Anti-IL5 therapy ####
dat.t <- FC.plot.dat %>% 
  filter(group2 != 'italic("Ex vivo")~"RV-infected vs media"' &
           term.ord %in% c("INTERFERON\nALPHA\nRESPONSE",
                           "INTERFERON\nGAMMA\nRESPONSE",
                           "INFLAMMATORY\nRESPONSE",
                           "EPITHELIAL\nMESENCHYMAL\nTRANSITION")) %>% 
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
                  show.legend = FALSE, size=3) +
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
ggsave("publication/fig/Fig7.GSEA.fold.change.pdf", fc.plot.all, 
       width = 20, height = 15)
