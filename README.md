# OCD antisense lncRNA overlap WGS analysis

This repository contains the analysis code for a whole-genome sequencing (WGS) study of rare, evolutionarily conserved variants in antisense lncRNA–protein-coding gene overlap regions and their contribution to obsessive–compulsive disorder (OCD) risk.

---

## Contents

- `scripts/` – Analysis scripts to reproduce each major result.
  - `01_define_antisense_overlaps.R`
  - `02_variant_annotation.sh`
  - `03_overlap_level_burden_and_SKATO.R`
  - `04_LOEUF_stratified_burden_and_bootstrap.ipynb`
  - `05_OCD_GWAS_enrichment_and_permutation.R`
  - `06_expression_and_coexpression_analysis_GTEx.R`
  
---

## Key analyses

The repository is organized around the main analyses described in the manuscript:

1. **Definition of antisense lncRNA–protein-coding overlap regions**
2. **Selection of rare, conserved variants**
3. **Overlap-level association tests (Fisher’s exact tests and SKAT-O)**
4. **LOEUF-stratified burden and fixed-K bootstrap**
5. **Enrichment of OCD GWAS risk genes in co-expressed networks**
6. **Expression and co-expression analyses for KNCN / MKNK1-AS1 (GTEx)**

Each of these corresponds to one or more scripts under `scripts/`.

---

### External resources (not included)

You will need to obtain the following from their original sources:

- **All of Us WGS and phenotype data** (controlled tier; accessed via the All of Us workbench)
- **Gene annotations** (e.g., GENCODE; antisense lncRNA and protein-coding gene models)
- **LOEUF scores** (e.g., gnomAD constraint metrics)
- **GERP++ conservation scores**
- **GTEx v8 expression data** for brain tissues


Paths to these resources are configured in `config/paths_example.yaml`.

---
