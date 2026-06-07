# MDD_ACC_sex_specific

**RNA-seq analysis of sex-specific transcriptomic changes in anterior cingulate cortex (BA24) in major depressive disorder**

## Описание проекта
This repository contains code and results from an analysis of sex-specific differences in the transcriptome profile of the anterior cingulate cortex (Brodmann area 24) in patients with major depressive disorder (MDD).

## Data
- **Source:** GEO GSE80655
- **Tissue:** Передняя поясная кора (BA24), постмортальные образцы
- **Sample:** 20 образцов
  - Women with MDD : 5
  - Healthy women: 5
  - Men with MDD: 5
  - Healthy men: 6

## Analysis Methods
- **DESeq2** — Differential Expression Analysis
  - Model: `~ SV1 + SV2 + tin_scaled + age_scaled + sex + group + sex:group`
- **Visualization:** PCA, volcano plots, heatmaps
- **Functional Analysis:** GO (Gene Ontology), GSEA (MSigDB Hallmark)
