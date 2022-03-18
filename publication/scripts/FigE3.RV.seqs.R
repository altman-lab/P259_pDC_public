library(tidyverse)

#### Data ####
attach("data_clean/P259_pDC_clean.RData")
meta <- read_csv("data_clean/P259_pDC_metadata.csv")

dat2 <- read_csv("data_raw/211102_P259-1_P259-2_RhinoViruses_normlized_read_counts.csv") %>% 
  rename(libID=libid) %>% 
  #filter libraries in final analysis
  filter(libID %in% dat.pDC.voom$targets$libID) %>% 
  #add metadata
  left_join(select(meta, libID, experiment, donorID, IL5, virus, virus.detail)) %>% 
  mutate(contrast=paste(IL5,virus, sep="_")) %>% 
  mutate(contrast=factor(contrast, 
                         levels = c("none_none","none_HRV",
                                    "AntiIL5_none","AntiIL5_HRV",
                                    "EOS.supp_none","EOS.supp_HRV")))

# Format for limma
count1 <- dat2 %>% 
  filter(experiment == "P259_1") %>% 
  select(libID, RhinoVirusA_normCount) %>% 
  pivot_longer(-libID) %>% 
  arrange(libID) %>% 
  pivot_wider(names_from = libID) %>% 
  column_to_rownames("name")

meta1 <- dat2 %>% 
  filter(experiment == "P259_1") %>% 
  arrange(libID) %>% 
  droplevels()

count2 <- dat2 %>% 
  filter(experiment == "P259_2") %>% 
  select(libID, RhinoVirusA_normCount) %>% 
  pivot_longer(-libID) %>% 
  arrange(libID) %>% 
  pivot_wider(names_from = libID) %>% 
  column_to_rownames("name")

meta2 <- dat2 %>% 
  filter(experiment == "P259_2") %>% 
  arrange(libID) %>% 
  droplevels()

#### limma ####
#Run in 3.P259.pDC_response.to.review
pval_1.contrast <- read_csv('review_response/RV.EOS.model.csv')
pval_2.contrast <- read_csv('review_response/RV.antiIL5.model.csv')

#### Plot ####
#Pvals for eos experiment
GOI.p1 <- pval_1.contrast %>% 
  mutate(dataset="P259.1")

#Pvals for aniti-IL5 experiment
GOI.p <- pval_2.contrast%>% 
  mutate(dataset="P259.2") %>% 
  #Combine with eos data
  bind_rows(GOI.p1) %>% 
  #Make symbols for plots
  mutate(symbol = ifelse(P.Value > 0.05, "n.s.", NA)) %>% 
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
  dplyr::select(facet.lab, group1, group2, symbol)

#Add y location for pval based on max expression in plot
GOI.pe <- dat2 %>% 
  #Max expression per gene and experiment
  group_by(experiment) %>% 
  summarise(max.e = max(RhinoVirusA_normCount, na.rm=TRUE)) %>% 
  ungroup() %>% 
  mutate(facet.lab = experiment,
         facet.lab = recode_factor(factor(facet.lab),
                                   "P259_1" = '"EOS sup"',
                                   "P259_2" = '"Anti-IL-5/5R"*alpha')) %>% 
  #Add to pval data
  right_join(GOI.p) %>% 
  arrange(experiment, group1, group2) %>% 
  mutate(y.position = max.e+2E3) %>% 
  ungroup()

plot <- dat2 %>% 
  mutate(x.lab=paste(experiment, virus.detail, IL5, sep="_"),
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
  arrange(facet.lab) %>% 
  ggplot(aes(x=x.lab, y=RhinoVirusA_normCount, color=donorID)) +
  geom_jitter(width=0.1, height=0) +
  stat_summary(fun.data=mean_sdl, 
               fun.args = list(mult=1), 
               geom="errorbar", color="black", width=0.25) +
  stat_summary(fun=mean, geom="errorbar", 
               aes(ymax=..y.., ymin=..y..),
               color="black", width=0.5) +
  facet_grid(~facet.lab, scales="free", 
             labeller = labeller(facet.lab=label_parsed)) +
  #Add pval
  stat_pvalue_manual(data=GOI.pe,
                     label="symbol", xmin="group1", xmax="group2") +
  # #Beautify
  theme_bw() +
  labs(x="", y="Normalized RV expression (x 1E6)") +
  theme(legend.position = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background =element_rect(fill="white"),
        axis.text.x=element_text(hjust=0.95,vjust=-0.5),
        plot.margin = margin(0.2,0.2,0.2,0.2,"cm")) +
  scale_x_discrete(labels = ggplot2:::parse_safe)
plot

ggsave("publication/fig/FigEX.RV.seqs.pdf", width=4, height=4)
