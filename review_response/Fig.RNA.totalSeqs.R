library(tidyverse)
library(readxl)

#### Data ####
attach("data_clean/P259_pDC_clean.RData")

dat <- read_excel("data_raw/P259-2 Final Annotation.xlsx") %>% 
  select("library sampleId", "RNA conc (ng/ul)")

dat <- read_csv("data_clean/P259_pDC_metadata.csv") %>% 
  select(libID, total_sequences) %>% 
  inner_join(dat, by=c("libID"="library sampleId")) %>% 
  mutate(`Quality filter` = ifelse(libID %in% dat.pDC.voom$targets$libID, "Pass", "Fail"))

#### Plot ####

dat %>% 
  ggplot(aes(x=`RNA conc (ng/ul)`, y=total_sequences)) +
  geom_point(aes(color=`Quality filter`), size=2) +
  theme_classic() +
  labs(y = "Raw sequences", x="RNA concentration (ng/ul)") +
  ggsave(filename = "review_response/Fig.RNA.totalSeqs.png",
         height=4, width=5)
