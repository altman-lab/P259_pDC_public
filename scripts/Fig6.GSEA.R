#### Setup ####
library(tidyverse)
library(cowplot)
library(limma)
set.seed(4389)

#### Data ####
H <- read_csv("results/GSEA_FoldChange/h_GSEA.result.csv")
fdr.cut <- 0.1

#### Subset significant ####
#list all terms with at least 1 significant treatment contrast
term.OI.t <- H %>% 
  filter(fgsea.FDR <= fdr.cut & gage.FDR <= fdr.cut) %>% 
  filter(group %in% c("AntiIL5_none.none_none.2", "AntiIL5_HRV.none_HRV.2",
                      "EOS.supp_none.none_none.1", "EOS.supp_HRV.none_HRV.1")) %>% 
  distinct(pathway) %>% unlist(use.names=FALSE)

#list all terms with at least 1 significant virus contrast 
term.OI.v <- H %>% 
  filter(fgsea.FDR <= fdr.cut & gage.FDR <= fdr.cut) %>% 
  filter(group %in% c("none_HRV.none_none.2","AntiIL5_HRV.AntiIL5_none.2",
                      "none_HRV.none_none.1","EOS.supp_HRV.EOS.supp_none.1")) %>% 
  distinct(pathway) %>% unlist(use.names=FALSE)

#Find terms with both t and v cutoffs
term.OI <- intersect(term.OI.t, term.OI.v)

#### Combine data and format labels ####
h.plot.dat <- H %>% 
  filter(pathway %in% term.OI) %>% 
  mutate(group1 = recode_factor(factor(group),
                                "AntiIL5_HRV.AntiIL5_none.2"="+ Anti-IL-5",
                                "none_HRV.none_none.2"="- Anti-IL-5",
                                "EOS.supp_HRV.EOS.supp_none.1"="+ EOS sup",
                                "none_HRV.none_none.1"="- EOS sup",
                                
                                "EOS.supp_HRV.none_HRV.1"="+ RV",
                                "AntiIL5_HRV.none_HRV.2"="+ RV",
                                "EOS.supp_none.none_none.1"="- RV",
                                "AntiIL5_none.none_none.2"="- RV"),
         group2 = recode_factor(factor(group),
            "none_HRV.none_none.1"='italic("Ex vivo")~"RV-infected vs media"',
            "none_HRV.none_none.2"='italic("Ex vivo")~"RV-infected vs media"',
            "EOS.supp_HRV.EOS.supp_none.1"='italic("Ex vivo")~"RV-infected vs media"',
            "AntiIL5_HRV.AntiIL5_none.2"='italic("Ex vivo")~"RV-infected vs media"',
                                
            "EOS.supp_HRV.none_HRV.1"='italic("Ex vivo")~"EOS supernatant vs none"',
            "AntiIL5_HRV.none_HRV.2"='italic("In vivo")~"Anti-IL-5 vs none"',
            "EOS.supp_none.none_none.1"='italic("Ex vivo")~"EOS supernatant vs none"',
            "AntiIL5_none.none_none.2"='italic("In vivo")~"Anti-IL-5 vs none"')) %>% 
  mutate(Significance = ifelse(fgsea.FDR <= 0.05 & 
                                 gage.FDR <= 0.05 & 
                                 fgsea.FC == gage.FC, 
                               "FDR < 0.05",
                               ifelse(fgsea.FDR <= 0.1 & 
                                        gage.FDR <= 0.1 & 
                                        fgsea.FC == gage.FC, 
                                      "FDR < 0.1", "NS"))) %>% 
  mutate(pathway = gsub("HALLMARK_", "", pathway),
         pathway2 = gsub("_", "\n", pathway)) %>% 
  #Reorder terms
  mutate(pathway.ord = factor(pathway2,
                              levels=c("INTERFERON\nALPHA\nRESPONSE",
                              "INTERFERON\nGAMMA\nRESPONSE",
                              "INFLAMMATORY\nRESPONSE",
                              "EPITHELIAL\nMESENCHYMAL\nTRANSITION")))

#### Plot A: virus ####
h.plot1 <- h.plot.dat %>% 
  filter(group2=='italic("Ex vivo")~"RV-infected vs media"') %>% 
  
  ggplot(aes(x=group1, y=fgsea.NES)) +
  geom_segment(aes(x=group1, xend=group1, y=0, yend=fgsea.NES)) +
  geom_point(size=2, aes(fill = Significance),
             shape=21, stroke=1) +
  geom_hline(yintercept = 0) +
  
  scale_fill_manual(values=c("firebrick", "#fdbb84", "white")) +
  lims(y=c(-3.6,3.6)) +
  coord_flip() +
  labs(x="", y="Normalized Enrichment Score (NES)") + 
  facet_grid(pathway.ord~group2, scales="free_y", 
             labeller = labeller(group2=label_parsed)) +
  theme_bw() +
  theme(legend.position = "none") +
  theme(strip.text.y = element_blank())+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background =element_rect(fill="white"))

#h.plot1

#### Plot B/C: EOS supernatant and anti-IL% therapy ####
h.plot2 <- h.plot.dat %>% 
  filter(group2!='italic("Ex vivo")~"RV-infected vs media"') %>% 
  
  ggplot(aes(x=group1, y=fgsea.NES)) +
  geom_segment(aes(x=group1, xend=group1, y=0, yend=fgsea.NES)) +
  geom_point(size=2, aes(fill = Significance),
             shape=21, stroke=1) +
  geom_hline(yintercept = 0) +
  
  scale_fill_manual(values=c("firebrick", "#fdbb84", "white")) +
  lims(y=c(-3.6,3.6)) +
  coord_flip() +
  labs(x="", y="Normalized Enrichment Score (NES)") + 
  facet_grid(pathway.ord~group2, scales="free_y",
             labeller = labeller(group2=label_parsed)) +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background =element_rect(fill="white"))

#h.plot2

#### Save ####
h.plot.all <- plot_grid(h.plot1, h.plot2, ncol=2, rel_widths = c(0.6,1),
                        labels=c("A","B")) 
#h.plot.all
ggsave("publication/fig/Fig6.GSEA.pdf", h.plot.all, 
       width = 11, height = 5)
