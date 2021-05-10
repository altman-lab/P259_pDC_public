#### Setup ####
library(tidyverse)
library(cowplot)
library(limma)
set.seed(4389)

#Define ggplot colors
group.cols <- c("none:none"="#dadaeb",
                "none:AntiIL5"="#9e9ac8",
                "none:EOS.supp"="#54278f",
                "HRV:none"="#c7e9c0",
                "HRV:AntiIL5"="#74c476",
                "HRV:EOS.supp"="#006d2c",
                "flu:none"="#fdae6b",
                "flu:AntiIL5"="#e6550d")

#### Data ####
load("data_clean/P259_pDC_clean.RData")

#### Calculate PCA ####
## EOS supernatant experiment
PCA1 <- as.data.frame(dat.pDC.voom_1$E) %>% 
  t() %>% 
  prcomp(scale. = TRUE)

##Axes labels
PC1.label1 <- paste("PC1 (", 
                    round(summary(PCA1)$importance[2,1]*100, digits=1), 
                    "%)", sep="")
PC2.label1 <-paste("PC2 (", 
                   round(summary(PCA1)$importance[2,2]*100, digits=1), 
                   "%)", sep="")

## Anti-IL5 therapy experiment
PCA2 <- as.data.frame(dat.pDC.voom_2$E) %>% 
  t() %>% 
  prcomp(scale. = TRUE)

##Axes labels
PC1.label2 <- paste("PC1 (", 
                    round(summary(PCA2)$importance[2,1]*100, digits=1), 
                    "%)", sep="")
PC2.label2 <-paste("PC2 (", 
                   round(summary(PCA2)$importance[2,2]*100, digits=1), 
                   "%)", sep="")

# Extract PC values and merge
PCA.dat <- as.data.frame(PCA1$x) %>% 
  rownames_to_column("libID") %>%
  bind_rows(rownames_to_column(as.data.frame(PCA2$x),"libID")) %>% 
  # Select PCs
  dplyr::select(libID, PC1, PC2) %>% 
  # Merge with metadata
  left_join(dat.pDC.voom$targets, by="libID") %>% 
  #Add facet labels
  mutate(facet.lab = experiment,
         facet.lab = recode(facet.lab,
                            "P259_1" = 'italic("Ex vivo")~"EOS supernatant"',
                            "P259_2" = 'italic("In vivo")~"Anti-IL-5/5R"*alpha'))

#### Plots ####
PCA1 <- PCA.dat %>% 
  filter(experiment=="P259_1") %>% 
  
  ggplot(aes(PC1, PC2, color=virus:IL5)) +
  geom_point(size=5) +
  #Beautify
  theme_classic(base_size = 16) +
  labs(x=PC1.label1, y=PC2.label1) +
  coord_fixed(ratio=1)  +
  scale_color_manual(values = group.cols, 
                     labels = c("media", "+ EOS sup", 
                                "+ RV", "+ RV + EOS sup"),
                     name="") +
  facet_wrap(~facet.lab, labeller = label_parsed) +
  theme(legend.position = "bottom") +
  guides(color=guide_legend(nrow=2))

PCA2 <- PCA.dat %>% 
  filter(experiment=="P259_2") %>% 
  
  ggplot(aes(PC1, PC2, color=virus:IL5)) +
  geom_point(size=5) +
  #Beautify
  theme_classic(base_size = 16) +
  labs(x=PC1.label2, y=PC2.label2) +
  coord_fixed(ratio=1) +
  scale_color_manual(values = group.cols, 
                     labels = c('media', bquote(+ "Anti-IL-5/5R"*alpha), 
                                '+ RV', bquote("+ RV + Anti-IL-5/5R"*alpha)),
                     name="") +
  facet_wrap(~facet.lab, labeller = label_parsed) +
  theme(legend.position = "bottom") +
  guides(color=guide_legend(nrow=2))

#### Save ####
ggsave("publication/fig/Fig5.PCA.pdf", 
       plot_grid(PCA1,PCA2, align = "hv", labels = c("A","B")),
       width=10, height=5)
