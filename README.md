***This repository contains code to reproduce the analyses presented in***
# APOE4 impairs myelination via cholesterol dysregulation in oligodendrocytes

Find the corresponding paper [here]().
### Data availability
- If you would like to *process the un-qc'ed counts matrix* and associated metadata, these files can be downloaded [here](link to synapse pre-qc data)
- If you would like to *access the fully-processed, annotated, and qc-ed data*, that data can be found [here](link to synapse qc_counts_data).
- All other datasets needed for this analysis are available through [the open science framework](https://osf.io/uyczk/).

### Reproduce analyses and plots
Follow these instructions to reproduce analyses and plots, as shown in the paper.

##### 1. Download the repository by running:

```bash
git clone https://github.com/djunamay/APOE4_impairs_myelination_via_cholesterol_dysregulation_in_oligodendrocytes.git
```

##### 2. Create the conda environment by running:

```bash
conda env create -f ./environment/apoe_env.yml
```

##### 3a. If you'd like to perform your own QC and celltype annotation from scratch
1. Download the un-qc'ed counts matrix and associated metadata from Synapse [here](link to pre-qc counts matrix and associated metadata) 
Please note, a data-use agreement must be submitted to access these data. Follow instructions on Synapse [here](https://www.synapse.org/#!RegisterAccount:0).

##### 3b. If you'd like to recapitulate our QC and celltype annotation
1. create a custom conda environment by running:
```bash
conda create -n actionet_legacy_env  -c conda-forge r-devtools
```
2. install the ACTIONet package (legacy version) from github by running the following:
```bash
conda activate actionet_legacy_env
R
```
In your R console, run this:
```R
library(devtools)
devtools::install_github("shmohammadi86/ACTIONet_legacy")
install.packages("harmony")
install.packages("R.utils")
install.packages("BiocManager")
BiocManager::install("BiocParallel")
install.packages("ggpubr")

```
3. Follow instructions from point 3a to access the pre-qc data.
4. Now run:
```bash
conda activate actionet_legacy_env
Rscript ./scripts/generate_raw_sce_object.r 
Rscript ./scripts/qc_and_annotation.r 
```

##### 3c. If you'd like to recapitulate any of the analyses presented in the paper
1. follow instructions [here](https://github.com/lhe17/nebula) to download the nebula package
2. follow instructions [here](https://github.com/immunogenomics/presto) to download the immunogenomics/presto package
4. Create the following /plots directory within this git repo
```
APOE4_impairs_myelination_via_cholesterol_dysregulation_in_oligodendrocytes
└───plots
    └───Extended_1
    └───Extended_2
    └───Extended_3
    └───Extended_4
    └───Extended_5
    └───Extended_6
    └───Figure_1
    └───Figure_2
    └───Figure_4
    └──qc_annotation
```
3. download necessary [data](https://osf.io/uyczk/) from OSF into a local directory named /data. This directory includes the following files:

| Data File                                                     | Description / Origin                                                                                                                                                                                                         |       
|---------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| pathway_databases/GO_Biological_Process_2018.txt              | from mayaan lab  [here](https://maayanlab.cloud/Enrichr/#libraries)                                                                                                                                                          |
| pathway_databases/HumanCyc_2016.txt                           | from mayaan lab [here](https://maayanlab.cloud/Enrichr/#libraries)                                                                                                                                                           |
| pathway_databases/KEGG_2019_Human.txt                         | from mayaan lab [here](https://maayanlab.cloud/Enrichr/#libraries)                                                                                                                                                           |
| pathway_databases/Reactome_2016.txt                           | from mayaan lab [here](https://maayanlab.cloud/Enrichr/#libraries)                                                                                                                                                           |
| iPSC_data/FPKM_table_AST.txt                                  | FPKM normalized counts from GEO accession number: [GSE102956](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102956)                                                                                                  |
| iPSC_data/FPKM_table_MIC.txt                                  | FPKM normalized counts from GEO accession number: [GSE102956](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102956)                                                                                                  |
| iPSC_data/FPKM_table_NEU.txt                                  | FPKM normalized counts from GEO accession number: [GSE102956](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102956)                                                                                                  |
| iPSC_data/FPKM_table_OPC.txt                                  | see methods section "Bulk RNA-sequencing from isogenic iPSC-derived oligodendroglia" in our paper; run ../scripts/get_opc_ipsc_counts_table.r                                                                                |
| iPSC_data/opc_ipsc_bulk_rnaseq_count_files/                   | see methods section "Bulk RNA-sequencing from isogenic iPSC-derived oligodendroglia" in our paper                                                                                                                            |                                                                                                                           |
| single_cell_data/Cell_group_colors.rds                        | colors assigned manually                                                                                                                                                                                                     |
| single_cell_data/expressed_genes_per_celltype.rds             | run ../scripts/get_expressed_genes_per_celltype.r                                                                                                                                                                            |
| single_cell_data/Mapping.rds                                  | colors assigned manually                                                                                                                                                                                                     |
| single_cell_data/RefCellTypeMarkers.adultBrain.rds            | Reference cell type marker genes were obtained from PsychENCODE, reported in Wang, D. et al. Comprehensive functional genomic resource and integrative model for the human brain. Science 362, (2018)                        |
| single_cell_data/PanglaoDB.by.organ.by.celltype.rds           | downloaded from Panglao database (https://panglaodb.se/); Franzén, O., Gan, L.-M. & Björkegren, J. L. M. PanglaoDB: a web server for exploration of mouse and human single-cell RNA sequencing data. Database  2019, (2019). |                                                                                                                                                                                                       |
| differentially_expressed_genes/E4_nebula_associations_Oli.rds | run ../scripts/nebula_degs.r                                                                                                                                                                                                 |
| differentially_expressed_genes/oli_wilcox_results.rds         | run ../scripts/get_wilcox_degs.r                                                                                                                                                                                             |
| differentially_expressed_genes/oli_wilcox_results_AD.rds      | run ../scripts/get_wilcox_degs.r                                                                                                                                                                                             |
| differentially_expressed_genes/oli_wilcox_results_noAD.rds    | run ../scripts/get_wilcox_degs.r                                                                                                                                                                                             |
| differentially_expressed_genes/OPC_deg_statistics.txt         | see methods section "Bulk RNA-sequencing from isogenic iPSC-derived oligodendroglia" in our paper                                                                                                                            |
| single_cell_data/ensembl.GRCh38p12.genes.complete.annot.rd    | |

4. Download the single-cell- and lipidomic-related data from Synapse [here] and add these data to the ./data directory according to the directories given in the table below. This includes the following files:

| Data File                                                       | Description / Origin                                                                                                                                  |       
|-----------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------|
| single_cell_data/individual_level_averages_per_celltype/Ast.csv | run ../scripts/get_individual_level_averages_object.r                                                                                                 |
| single_cell_data/individual_level_averages_per_celltype/Ex.csv  | run ../scripts/get_individual_level_averages_object.r                                                                                                 |
| single_cell_data/individual_level_averages_per_celltype/In.csv  | run ../scripts/get_individual_level_averages_object.r                                                                                                 |
| single_cell_data/individual_level_averages_per_celltype/Mic.csv | run ../scripts/get_individual_level_averages_object.r                                                                                                 |
| single_cell_data/individual_level_averages_per_celltype/Oli.csv | run ../scripts/get_individual_level_averages_object.r                                                                                                 |
| single_cell_data/individual_level_averages_per_celltype/Opc.csv | run ../scripts/get_individual_level_averages_object.r                                                                                                 |
| single_cell_data/qc_counts_data/qc_column_metadata.csv | collected and shared by ROSMAP (Dr David Bennett and colleagues)                                                                                      |
| single_cell_data/qc_counts_data/qc_counts.mtx | see methods sections "Quality control for cell inclusion" and "Clustering analysis and QC filtering" in our paper, and ../scripts/qc_and_annotation.r |
| single_cell_data/qc_counts_data/qc_gene_names.txt | see methods section "snRNA-seq data preprocessing" in our paper                                                                                       |
| single_cell_data/raw_counts_data/column_metadata.csv | collected and shared by ROSMAP (Dr David Bennett and colleagues)                                                                                      |
| single_cell_data/raw_counts_data/gene_names.csv | see methods section "snRNA-seq data preprocessing" in our paper                                                                                       |
| single_cell_data/raw_counts_data/raw_counts.mtx | see methods section "snRNA-seq data preprocessing" in our paper                                                                                       |
| cc_lipidomics/Lipidomics_RawData.csv | see methods section "Untargeted lipidomics of post-mortem corpus callosum" in our paper                                                               |
| cc_lipidomics/Lipidomics_RawData_2.csv | see methods section "Untargeted lipidomics of post-mortem corpus callosum" in our paper                                                               |
| pfc_lipidomics/ChE_summary_cyc_05312022_all_samples.csv | see methods section "Untargeted Lipidomics on post-mortem prefrontal cortex" in our paper                                                             | 
| pfc_lipidomics/ROSMAP_Lipidomics_Emory_biospecimen_metadata.csv | accession through Synapse [here](https://www.synapse.org/#!Synapse:syn26475187)                                                                       |

5. Download the iPSC RNA-seq counts tables for astrocytes, microglia, and neurons from from GEO accession number: [GSE102956](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102956)
6. create an empty directory in ./data titled "other_analyses_outputs"

7. Now run the following code snippets to recapitulate the analysis:

Run this first:
```bash
conda activate apoe_env
Rscript ./scripts/generate_qc_sce_object.r 
Rscript ./scripts/get_pathways.r
Rscript ./scripts/nebula_degs.r # check that this produces the deg table I was using
Rscript ./scripts/get_expressed_genes_per_celltype.r 
Rscript ./scripts/get_individual_level_averages_object.r 
Rscript ./scripts/get_opc_ipsc_counts_table.r 
Rscript ./scripts/make_metadata_file.r 
```

Figure 1
```bash
conda activate apoe_env
Rscript ./scripts/pathway_analyses.r 
Rscript ./scripts/get_figure_1_plots.r
```

Figure 2
```bash
conda activate apoe_env
Rscript ./scripts/dissecting_cholesterol_dysregulation.r
Rscript ./scripts/lipidomic_analysis_cc.r 
Rscript ./scripts/get_figure_2_plots.r  
```
Figure 4
```bash
conda activate apoe_env
Rscript ./scripts/get_wilcox_degs.r
Rscript ./scripts/get_figure_4_plots.r
```

Extended Data Figure 1
```bash
conda activate apoe_env
Rscript ./scripts/plots_for_extended_data_figure_1.r 
```

Extended Data Figure 2
```bash
conda activate apoe_env
Rscript ./scripts/fgsea_analysis.r
Rscript ./scripts/pseudo_bulk.r
Rscript ./scripts/plots_for_extended_data_figure2.r 
Rscript ./scripts/e4_effects_stratification_by_AD.r
Rscript ./scripts/e4_stratification_plots.r
```

Extended Data Figure 3
```bash
conda activate apoe_env
Rscript ./scripts/lipidomic_analysis_pfc_get_data.r  # check 
Rscript ./scripts/lipidomic_analysis_pfc_make_plots.r # check
```

Extended Data Figure 5
```bash
conda activate apoe_env
Rscript ./scripts/comparison_of_ipsc_and_brain_get_scaled_matrices.r
Rscript ./scripts/comparison_of_ipsc_and_brain_plots.r
Rscript ./scripts/APOE_expression_oligodendrocytes.r
Rscript ./scripts/apoe_expression_ipsc.r  
```

Extended Data Figure 7
```bash
conda activate apoe_env
Rscript ./scripts/get_postmortem_er_stress_pathways.r
Rscript ./scripts/er_postmortem_plots.r  
Rscript ./scripts/ipsc_gene_perturbations.r  
```

Extended Data Figure 8
```bash
conda activate apoe_env
Rscript ./scripts/get_wilcox_myelin_plots.r
```

Or, run the full analysis pipeline:
```bash
conda activate apoe_env
chmod +x run_all_analyses.sh
./scripts/run_all_analyses.sh
```
If you have any questions, notice any inconsistencies, need help, or would like to brainstorm future collaborations and ideas, please don't hesitate to reach out: djuna@mit.edu
