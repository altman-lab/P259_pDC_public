# P259 pDC, U. of Washington

Analyses associated with 

> Dill-McFarland KA, Schwartz JT, Fulkerson PC, Zhao H, Shao B, Altman MC, Gill MA. 2021. Eosinophil-mediated suppression and Anti-IL-5 enhancement of plasmacytoid dendritic cell interferon responses in asthma. J Allergy Clin Immunol. In revision

## Abstract
Background: Virus-induced IFNα secretion by plasmacytoid dendritic cells (pDCs) is negatively impacted by IgE and has been linked to asthma exacerbations. Eosinophils, another contributor to Type 2 (T2) inflammation, are also associated with asthma severity.

Objective: To investigate the impact of eosinophils on pDC antiviral IFN responses and determine whether anti-IL-5/5Rα therapy enhances pDC antiviral function.

Methods: Blood pDCs purified from anonymous donors were stimulated ex vivo with rhinovirus-16 (RV) in the presence or absence of eosinophils/eosinophil supernatants. IFNα was measured in supernatants and RNA collected for RNA-sequencing and differential gene expression analysis. Next, purified pDCs from 8 individuals with moderate-severe asthma, treated or not treated with anti-IL-5/5Rα therapy, were cultured with or without RV; IFNα secretion and differential gene expression analysis was compared between groups.

Results: Exposure to either eosinophils or eosinophil supernatants inhibited rhinovirus (RV)-induced pDC IFNα secretion in a dose-dependent manner and did not impact pDC viability. Eosinophil-derived neurotoxin (EDN) partially recapitulated this inhibitory effect on pDC IFNα secretion. Transcriptome analysis revealed global repression of pDC IFN response patterns by eosinophils, most notably in basal expression of interferon stimulated genes (ISGs). Increased RV-induced IFNα secretion and transcription as well as increased basal ISG expression was detected in pDCs from participants treated with anti-IL-5/5Rα therapy.

Conclusion: Our findings highlight a novel mechanism through which T2 inflammation regulates pDC IFNα responses relevant to RV respiratory infections in the context of eosinophilic airway disease and suggest one potential mechanism through which eosinophil66 depleting therapies may reduce the severity of RV illnesses 

PMID: TBD

## data_clean
Normalized log2 counts per million (CPM) in genes, patient metadata, and Broad reference gene sets.

## data_raw
Unnormalized RNA-seq counts in genes and HGNC gene symbol key.

## publication
Figures and tables in publication. Reproducible scripts for figure and table generation.

## results

* **gene_level**:  linear models of gene expression
* **GSEA_FoldChange**:  Gene set enrichment analysis

## Rmarkdowns
All analyses run from the environment created in P259_pDC.Rproj

1. **data.cleaning**:  Data import and cleaning
		+ Final environment saved as data_clean/P259_pDC_clean.RData
2. **model.selection**:  Linear model exploration including interaction, pairwise contrasts, and post-pre (delta)
		+ Results saved to results/gene_level
3. **GSEA**:  Gene set enrichment analysis of Broad Hallmark, curated, and GO gene sets
		+ Results saved to results/GSEA_FoldChange
