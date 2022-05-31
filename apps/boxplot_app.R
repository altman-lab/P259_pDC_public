#P259
library(BiocManager)
options(repos = BiocManager::repositories())
#data
library(readr)
library(dplyr)
library(tibble)
library(tidyr)
library(limma)
#plots
library(ggplot2)
#shiny
library(shiny)

#### Data prep ####
# Expression data
# attach("data_clean/P259_pDC_clean.RData")
# dat1 <- as.data.frame(dat.pDC.voom_1$E) %>%
#   rownames_to_column("geneName") %>%
#   pivot_longer(-geneName, names_to = "libID") %>%
#   left_join(dat.pDC.voom_1$targets) %>%
#   left_join(dat.pDC.voom_1$genes) %>%
#   mutate(experiment="EOS supernatant") %>%
#   select(experiment, cellType, libID, geneName,
#          hgnc_symbol, IL5, virus, value) %>%
#   rename(treatment=IL5)
# dat2 <- as.data.frame(dat.pDC.voom_2$E) %>%
#   rownames_to_column("geneName") %>%
#   pivot_longer(-geneName, names_to = "libID") %>%
#   left_join(dat.pDC.voom_2$targets) %>%
#   left_join(dat.pDC.voom_2$genes) %>%
#   mutate(experiment = "anti-IL5/5Ra") %>%
#   select(experiment, cellType, libID, geneName,
#          hgnc_symbol, IL5, virus, value) %>%
#   rename(treatment=IL5)
# 
# attach("data_clean/P259_EOS_clean.RData")
# dat3 <- as.data.frame(dat.EOS.voom$E) %>%
#   rownames_to_column("geneName") %>%
#   pivot_longer(-geneName, names_to = "libID") %>%
#   left_join(dat.EOS.voom$targets) %>%
#   left_join(dat.EOS.voom$genes) %>%
#   mutate(experiment = "EOS supernatant") %>%
#   select(experiment, cellType, libID, geneName,
#          hgnc_symbol, IFNa, virus, value) %>%
#   rename(treatment=IFNa)
# 
# dat.all <- bind_rows(dat1,dat2,dat3) %>%
#   mutate(gene = paste0(hgnc_symbol, " (", geneName, ")"),
#          virus = recode(virus, "HRV"="RV"),
#          x=paste(virus, treatment, sep="\n"),
#          x = factor(x, c("none\nnone",
#                          "none\nAntiIL5",
#                          "none\nEOS.supp",
#                          "none\nIFNa",
#                          "RV\nnone",
#                          "RV\nAntiIL5",
#                          "RV\nEOS.supp")),
#          cellType=factor(cellType, c("pDC","EOS")))
# 
# #FDR data
# fdr1 <- read_csv("results/gene_level/P259.1_gene_pval.csv") %>%
#   mutate(gene = paste0(hgnc_symbol, " (", geneName, ")")) %>%
# 
#   filter(group %in% c("none_HRV - none_none","EOS.supp_HRV - EOS.supp_none",
#                       "EOS.supp_none - none_none","EOS.supp_HRV - none_HRV")) %>%
#   mutate(between = recode(group,
#                            "none_HRV - none_none"="RV vs none",
#                            "EOS.supp_HRV - EOS.supp_none"="RV vs none",
#                            "EOS.supp_none - none_none"="EOS sup vs none",
#                            "EOS.supp_HRV - none_HRV"="EOS sup vs none"),
#          within = recode(group,
#                             "none_HRV - none_none"="untreated",
#                             "EOS.supp_HRV - EOS.supp_none"="EOS sup treated",
#                             "EOS.supp_none - none_none"="uninfected",
#                             "EOS.supp_HRV - none_HRV"="RV infected"),
#          experiment="EOS supernatant",
#          cellType="pDC") %>%
#   rename(P=P.Value, FDR=adj.P.Val) %>%
#   select(experiment, cellType, gene, between, within, logFC, P, FDR)
# 
# fdr2 <- read_csv("results/gene_level/P259.2_gene_pval.csv") %>%
#   mutate(gene = paste0(hgnc_symbol, " (", geneName, ")")) %>%
# 
#   filter(group %in% c("none_HRV - none_none","AntiIL5_HRV - AntiIL5_none",
#                       "AntiIL5_none - none_none","AntiIL5_HRV - none_HRV")) %>%
#   mutate(between = recode(group,
#                           "none_HRV - none_none"="RV vs none",
#                           "AntiIL5_HRV - AntiIL5_none"="RV vs none",
#                           "AntiIL5_none - none_none"="AntiIL5 vs none",
#                           "AntiIL5_HRV - none_HRV"="AntiIL5 vs none"),
#          within = recode(group,
#                          "none_HRV - none_none"="untreated",
#                          "AntiIL5_HRV - AntiIL5_none"="antiIL5 treated",
#                          "AntiIL5_none - none_none"="uninfected",
#                          "AntiIL5_HRV - none_HRV"="RV infected"),
#          experiment="anti-IL5/5Ra",
#          cellType="pDC") %>%
#   rename(P=P.Value, FDR=adj.P.Val) %>%
#   select(experiment, cellType, gene, between, within, logFC, P, FDR)
# 
# fdr.all <- bind_rows(fdr1,fdr2) %>%
#   mutate(` ` = ifelse(FDR < 0.001, "***",
#                       ifelse(FDR < 0.01, "**",
#                              ifelse(FDR < 0.1, "*"," ")))) %>%
#   mutate(FDR = ifelse(FDR < 0.001, formatC(FDR, digits=2, format = "E"),
#                                           round(FDR, digits=3))) %>%
#   mutate(between = factor(between, c("AntiIL5 vs none", "EOS sup vs none","RV vs none")),
#          within = factor(within, c("uninfected","RV infected",
#                             "untreated", "antiIL5 treated", "EOS sup treated"))) %>%
#   arrange(experiment, between, within)
# 
# #Genes
# genes <- sort(unique(dat.all$gene))
# save(dat.all, fdr.all, genes, file="apps/boxplot_app_data.RData")

#### App ####
load("boxplot_app_data.RData")
#UI
shinyApp(
  ui = fluidPage(
    tags$head(
      tags$style(HTML(".selectize-input {
        height: 30px;
        width: 400px;
        font-size: 11pt;
        padding-top: 5px;
      }"
      ))),
    fluidRow(inputPanel(selectizeInput(inputId="choose_gene",
                                         label="Gene HGNC (ENSEMBL)",
                                         choices=NULL,
                                         selected=FALSE)),
               column(5, "pDC contrast model", tableOutput("tab1"), tableOutput("tab2")),
               column(7, "If one or more panels does not appear, this gene was not sufficiently expressed for analysis in the missing data set.", plotOutput("p1",height="400px")))
    ),
  
  server = function(input, output, session) {
    updateSelectizeInput(session = session, 
                         inputId = 'choose_gene', 
                         choices = genes, 
                         server = TRUE,
                         selected = FALSE,
                         options = list(placeholder = 'gene'))
    
    output$tab1 = renderTable({
      fdr_tab1 <- fdr.all %>% 
        filter(gene == input$choose_gene) %>% 
        select(-gene,-cellType) %>% 
        filter(experiment == "anti-IL5/5Ra")
      
      fdr_tab1
    })
    output$tab2 = renderTable({
      fdr_tab2 <- fdr.all %>% 
        filter(gene == input$choose_gene) %>% 
        select(-gene,-cellType) %>% 
        filter(experiment == "EOS supernatant")
      
      fdr_tab2
    })
    output$p1 = renderPlot({
      p1 <- dat.all %>%
        filter(gene == input$choose_gene) %>%
        # filter(gene=="TSPAN6 (ENSG00000000003)") %>% 
        ggplot(aes(y=value, x=x)) +
        geom_boxplot(outlier.shape=NA) +
        geom_jitter(width=0.2, height=0) +
        labs(y="Normalized log2 expression", x="",
             title = input$choose_gene) +
             # title="TSPAN6 (ENSG00000000003)") +
        facet_wrap(~experiment+cellType, scales="free") +
        theme_classic(base_size = 16) +
        theme(panel.border = element_rect(colour = "black", fill=NA))
      
      p1
    })
  },
  
  options = list(height = 500)
)
