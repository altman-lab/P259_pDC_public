#### Setup ####
library(tidyverse)
library(cowplot)
library(limma)
set.seed(4389)

#### Data ####
H <- read_csv("results/GSEA/h_GSEA.result.csv")
fdr.cut <- 0.1

#### Subset significant ####
#list all terms with at least 1 significant treatment contrast
term.OI.t <- H %>% 
  filter(fgsea.FDR <= fdr.cut & gage.FDR <= fdr.cut) %>% 
  filter(group %in% c("AntiIL5_nonenone_none", "AntiIL5_HRVnone_HRV",
                      "EOSsupp_nonenone_none", "EOSsupp_HRVnone_HRV")) %>% 
  distinct(pathway) %>% unlist(use.names=FALSE)

#list all terms with at least 1 significant virus contrast 
term.OI.v <- H %>% 
  filter(fgsea.FDR <= fdr.cut & gage.FDR <= fdr.cut) %>% 
  filter(group %in% c("none_HRVnone_none","AntiIL5_HRVAntiIL5_none",
                      "EOSsupp_HRVEOSsupp_none")) %>% 
  distinct(pathway) %>% unlist(use.names=FALSE)

#Find terms with both t and v cutoffs
term.OI <- intersect(term.OI.t, term.OI.v)


#### Combine data and format labels ####
h.plot.dat <- H %>% 
  filter(pathway %in% term.OI) %>% 
  #Add experi group
  mutate(experi = c(rep("AntiIL5",nrow(.)/2),rep("EOSsup",nrow(.)/2))) %>% 
  mutate(group1 = recode(group,
                                "AntiIL5_HRVAntiIL5_none"='"+ Anti-IL-5/5R"*alpha',
                                "EOSsupp_HRVEOSsupp_none"='"+ EOS sup"',

                                "EOSsupp_HRVnone_HRV"="+ RV",
                                "AntiIL5_HRVnone_HRV"="+ RV",
                                "EOSsupp_nonenone_none"="- RV",
                                "AntiIL5_nonenone_none"="- RV"),
         group1 = ifelse(group1 == "none_HRVnone_none" & experi == "AntiIL5",
                         '"- Anti-IL-5/5R"*alpha',
                         ifelse(group1 == "none_HRVnone_none" & experi == "EOSsup",
                                '"- EOS sup"', group1)),
         group2 = recode_factor(factor(group),
            "none_HRVnone_none"='italic("Ex vivo")~"RV-infected vs media"',
            "EOSsupp_HRVEOSsupp_none"='italic("Ex vivo")~"RV-infected vs media"',
            "AntiIL5_HRVAntiIL5_none"='italic("Ex vivo")~"RV-infected vs media"',
                                
            "EOSsupp_HRVnone_HRV"='italic("Ex vivo")~"EOS supernatant vs none"',
            "AntiIL5_HRVnone_HRV"='italic("In vivo")~"Anti-IL-5/5R"*alpha~"vs none"',
            "EOSsupp_nonenone_none"='italic("Ex vivo")~"EOS supernatant vs none"',
            "AntiIL5_nonenone_none"='italic("In vivo")~"Anti-IL-5/5R"*alpha~"vs none"')) %>% 
  #order factors
  mutate(group1 = factor(group1, levels=c('"+ Anti-IL-5/5R"*alpha','"- Anti-IL-5/5R"*alpha',
                                          '"+ EOS sup"','"- EOS sup"', 
                                          '+ RV', '- RV'))) %>% 
  mutate(Significance = ifelse(fgsea.FDR <= 0.05 & 
                                 gage.FDR <= 0.05 & 
                                 fgsea.FC == gage.FC, 
                               "FDR < 0.05",
                               ifelse(fgsea.FDR <= 0.1 & 
                                        gage.FDR <= 0.1 & 
                                        fgsea.FC == gage.FC, 
                                      "FDR < 0.1", "NS"))) %>% 
  mutate(pathway = gsub("HALLMARK_", "", pathway)) %>% 
  #Reorder terms
  mutate(pathway.ord = recode_factor(factor(pathway),
                                     "INTERFERON_ALPHA_RESPONSE"="IFNA response",
                                     "INTERFERON_GAMMA_RESPONSE"="IFNG response",
                                     "INFLAMMATORY_RESPONSE"="Inflammatory\nresponse",
                                     "EPITHELIAL_MESENCHYMAL_TRANSITION"="Epithelial\nmesenchymal\ntransition"))

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
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background =element_rect(fill="white")) +
  scale_x_discrete(labels = ggplot2:::parse_safe)

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
h.plot.all <- plot_grid(h.plot1, h.plot2, ncol=2, rel_widths = c(0.55,1),
                        labels=c("A","B")) 
#h.plot.all
ggsave("publication/fig/Fig6.GSEA.pdf", h.plot.all, 
       width = 11, height = 5)
