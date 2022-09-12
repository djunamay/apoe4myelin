***This repository contains code to reproduce the analyses presented in***
# APOE4 impairs myelination via cholesterol dysregulation in oligodendrocytes

Find the corresponding paper [here]().
### Data availability
- If you would like to *process the raw Fastq* files and associated metadata, these files can be downloaded [here](link to synapse fastq files)
- If you would like to *access the fully-processed, annotated, and qc-ed data*, that data can be found [here](link to single cell experiment object).
- All other datasets needed for this analysis are available through [OSF](https://osf.io/uyczk/).

### Reproduce analyses and plots
Follow these instructions to reproduce analyses and plots, as shown in the paper.

##### 1. Download the repository by running:

```bash
git clone https://github.com/djunamay/APOE4_impairs_myelination_via_cholesterol_dysregulation_in_oligodendrocytes.git
```

##### 2. Create the conda environment by running:

```bash
conda env create -f ../environment/apoe_env.yml
```

##### 3a. If you'd like to perform your own QC and celltype annotation from scratch
1. Download the FASTQ files from Synapse [here](...) 
2. Download the metadata files from Synapse [here](...)

##### 3b. If you'd like to recapitulate our QC and celltype annotation
1. Follow instructions [here](https://github.com/shmohammadi86/ACTIONet/tree/R-release) to install the ACTIONet package.
2. Download the single cell data counts matrix from Synapse [here](...)
3. Download the metadata files from Synapse [here](...)
3. Now run:
```bash
conda activate apoe_env
Rscript ../scripts/qc_and_annotation.r
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
    └───Figure_3
```
3. download necessary [data](https://osf.io/uyczk/) from OSF into a local directory named /data. This directory includes the following files:

| Data File                                                             | Description / Origin                                                |       
|-----------------------------------------------------------------------|---------------------------------------------------------------------|
| pathway_databases/GO_Biological_Process_2018.txt                      | from mayaan lab  [here](https://maayanlab.cloud/Enrichr/#libraries) |
| pathway_databases/HumanCyc_2016.txt                                   | from mayaan lab [here](https://maayanlab.cloud/Enrichr/#libraries)  |
| pathway_databases/KEGG_2019_Human.txt                                 | from mayaan lab [here](https://maayanlab.cloud/Enrichr/#libraries)  |
| pathway_databases/Reactome_2016.txt                                   | from mayaan lab [here](https://maayanlab.cloud/Enrichr/#libraries)  |
| iPSC_data/FPKM_table_AST.txt                                          |                                                                     |
| iPSC_data/FPKM_table_MIC.txt                                          |                                                                     |
| iPSC_data/FPKM_table_NEU.txt                                          |                                                                     |
| iPSC_data/FPKM_table_OPC.txt                                          |                                                                     |
| single_cell_data/Cell_group_colors.rds                                | NA                                                                  |
| single_cell_data/expressed_genes_per_celltype.rds                     |                                                                     |
| single_cell_data/Mapping.rds                                          | NA                                                                  |
| single_cell_data/RefCellTypeMarkers.adultBrain.rds                    |                                                                     |
| single_cell_data/PanglaoDB.by.organ.by.celltype.rds                   |                                                                     |
| differentially_expressed_genes/E4_nebula_associations_by_celltype.rds |                                                                     |
| differentially_expressed_genes/oli_wilcox_results.rds                 | run ../scripts/get_nebula_degs.r                                    |
| differentially_expressed_genes/oli_wilcox_results_AD.rds              | run ../scripts/get_nebula_degs.r                                    |
| differentially_expressed_genes/oli_wilcox_results_noAD.rds            | run ../scripts/get_nebula_degs.r                                    |
| differentially_expressed_genes/OPC_deg_statistics.txt                 |                                                                     |


6. Download the single-cell-related data from Synapse [here] and add these data to the ./data directory into their respective sub-folders. This includes the following files:

| Data File                                                                          | Description / Origin |       
|------------------------------------------------------------------------------------|----------------------|
| single_cell_data/single_cell_experiment_object.rds                                 | 
| single_cell_data/individual_level_averages_per_celltype.rds                        |
| single_cell_data/metadata_by_individual.csv                                        |
| single_cell_data/metadata_PFC_all_individuals_092520.tsv                           |
| lipidomic_datasets/cc_lipidomics/Lipidomics_RawData_2.csv                          |
| lipidomic_datasets/cc_lipidomics/Lipidomics_RawData.csv                            |
| lipidomic_datasets/pfc_lipidomics/ChE_summary_cyc_05342022_all_samples.csv         |
| lipidomic_datasets/pfc_lipidomics/ROSMAP_Lipidomics_Emory_biospecimen_metadata.csv |
| lipidomic_datasets/pfc_lipidomics/single_cell_counts.csv                           | pre-qc               
| lipidomic_datasets/pfc_lipidomics/single_cell_rowdata.csv                          | pre-qc               
| lipidomic_datasets/pfc_lipidomics/single_cell_coldata.csv                          | pre-qc               

3. create an empty directory in ./data titled other_analyses_outputs

3. Now run the following code snippets to recapitulate the analysis:

Figure 1
```bash
conda activate apoe_env
Rscript ../scripts/get_pathways.r
Rscript ../scripts/get_nebula_degs.r 
Rscript ../scripts/pathway_analyses.r
Rscript ../scripts/get_figure_1_plots.r
```

Figure 2
```bash
conda activate apoe_env
Rscript ../scripts/dissecting_cholesterol_dysregulation.r
Rscript ../scripts/lipidomic_analysis_cc.r
Rscript ../scripts/get_figure_2_plots.r  
```
Figure 3
```bash
conda activate apoe_env
Rscript ../scripts/get_wilcox_degs.r
Rscript ../scripts/get_figure_3_plots.r
```

Extended Data Figure 1
```bash
conda activate apoe_env
Rscript ../scripts/plots_for_extended_data_figure_1.r 
```

Extended Data Figure 2
```bash
conda activate apoe_env
Rscript ../scripts/fgsea_analysis.r
Rscript ../scripts/get_postmortem_er_stress_pathways.r
Rscript ../scripts/pseudo_bulk.r
Rscript ../scripts/plots_for_extended_data_figure2.r 
```

Extended Data Figure 3
```bash
conda activate apoe_env
Rscript ../scripts/e4_effects_stratification_by_AD.r
Rscript ../scripts/e4_stratification_plots.r
```

Extended Data Figure 4
```bash
conda activate apoe_env
Rscript ../scripts/lipidomic_analysis_pfc_get_data.r
```

Extended Data Figure 6
```bash
conda activate apoe_env
Rscript ../scripts/comparison_of_ipsc_and_brain_get_scaled_matrices.r
Rscript ../scripts/comparison_of_ipsc_and_brain_plots.r
Rscript ../scripts/APOE_expression_oligodendrocytes.r
Rscript ../scripts/apoe_expression_ipsc.r  
```

Extended Data Figure 8
```bash
conda activate apoe_env
Rscript ../scripts/er_postmortem_plots.r  
Rscript ../scripts/ipsc_gene_perturbations.r  
```

Extended Data Figure 9
```bash
conda activate apoe_env
Rscript ../scripts/get_wilcox_myelin_plots.r
```

