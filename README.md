***This repository contains code to reproduce the analyses presented in***
# APOE4 impairs myelination via cholesterol dysregulation in oligodendrocytes

**Joel W. Blanchard<sup>1,2,7</sup>, Leyla Anne Akay<sup>1,2,3</sup>, Jose Davila-Velderrain<sup>1,2,3,8</sup>, Djuna von Maydel<sup>1,2,3</sup>**, Hansruedi Mathys<sup>1,2,9, Shawns M. Davidson<sup>6</sup>, Audrey Effenberger<sup>1,2</sup>, Michael Bula1<sup>2</sup>, Martin Kahn<sup>1,2</sup>, Cristina Blanco<sup>1,2</sup>, Emre Agbas<sup>1,2</sup>, Ayesha Ng<sup>1,2</sup>, Xueqiao Jiang<sup>1,2</sup>, Yuan-Ta Lin<sup>1,2</sup>, Liwang Liu<sup>1,2</sup>, Tak Ko<sup>1</sup>, William T. Ralvenius<sup>1,2</sup>, David A. Bennett<sup>5</sup>, Hugh P. Cam<sup>1,2</sup>, *Manolis Kellis<sup>3,4</sup>, Li-Huei Tsai<sup>1,2,4</sup>*

1. Picower Institute for Learning and Memory, Massachusetts Institute of Technology, Massachusetts Institute of Technology, Cambridge, MA 02139, USA,
2. Department of Brain and Cognitive Sciences, Massachusetts Institute of Technology, Massachusetts Institute of Technology, Cambridge, MA 02139, USA.
3. MIT Computer Science and Artificial Intelligence Laboratory, Cambridge, MA 02139, USA.
4. Broad Institute of Harvard and MIT, Cambridge, MA 02139, USA.
5. Rush Alzheimer’s Disease Center, Rush University Medical Center, Chicago, IL USA.
6. Lewis-Sigler Institute for Integrative Genomics, Princeton University, Princeton New Jersey.
7. Current address: Department of Neuroscience, Black Family Stem Cell Institute, Ronald M. Loeb Center for Alzheimer’s Disease, Icahn School of Medicine at Mt. Sinai, New York, NY 10029
8. Current address: Human Technopole, Viale Rita Levi-Montalcini 1, 20157, Milan, Italy
9. Current address: Department of Neurobiology, University of Pittsburgh, Pittsburgh, PA 15261

**These authors contributed equally**\
*To whom correspondence should be addressed: lhtsai@mit.edu, manoli@mit.edu*B

### Data availability
- If you would like to *process the raw Fastq* files and associated metadata, these files can be downloaded [here](link to synapse). [coming soon]
- If you would like to *perform your own quality control and celltype annotation* on the aggregated counts matrix and associated metadata, that data can be found [here](link to synapse).
- If you would like to *access the fully-processed, annotated, and qc-ed data*, that data can be found [here](https://www.dropbox.com/sh/8i8hhvsyzoqdpzu/AAB8uV6pHN56OvFoRJF3Keqea?dl=0).


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
1. Download the FASTQ files [here](...) [coming soon]

##### 3b. If you'd like to recapitulate our QC and celltype annotation
1. Follow instructions [here](https://github.com/shmohammadi86/ACTIONet/tree/R-release) to install the ACTIONet package.
2. Download the cellranger aggregation outputs [here]
3. Now run:
```bash
conda activate apoe_env
Rscript ../scripts/qc_and_annotation.r #TODO: check/edit this
```

##### 3c. If you'd like to recapitulate any of the analyses presented in the paper
1. follow instructions [here](https://github.com/lhe17/nebula) to download the nebula package
2. follow instructions [here](https://github.com/immunogenomics/presto) to download the immunogenomics/presto package
3. Create the following /plots directory within this git repo

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

4. download the [data](https://www.dropbox.com/sh/8i8hhvsyzoqdpzu/AAB8uV6pHN56OvFoRJF3Keqea?dl=0) directory into a local directory named /data. This directory contains the following data:

```
data
└───single_cell_data
    └───Cell_group_colors.rds
    └───ensembl.GRCh38p12.genes.complete.annot.rds
    └───single_cell_experiment_object.rds
    └───expressed_genes_per_celltype.rds
    └───individual_level_averages_per_celltype.rds
    └───Mapping
    └───metadata_by_individual.csv
    └───metadata_PFC_all_individuals_092520.tsv
    └───Metadata.APOE.project.rds
    └───RefCellTypeMarkers.adultBrain.rds
    └───Summary.data.celltype.rds
    └───PanglaoDB.by.organ.by.celltype.rds
└───differentially_expressed_genes_data
    └───nebula_oli_degs.rds
    └───wilcoxon_degs_all.rds
    └───wilcoxon_degs_AD.rds
    └───wilcoxon_degs_noAD.rds
    └───pseudo_bulk_degs_single_cell_all_celltypes.csv
└───iPSC_data
    └───ipsc_bulk_rnaseq_count_files
        └───...
    └───OPC_DEG_statistics.txt
    └───FPKM_table_AST.txt
    └───FPKM_table_MIC.txt
    └───FPKM_table_OPC.txt
    └───FPKM_table_NEU.txt # modify this name
└───pathway_databases
    └───GO_Biological_Process_2018.txt
    └───HumanCyc_2016.txt
    └───KEGG_2019_Human.txt
    └───Reactome_2016.txt
└───lipidomics_dataset
    └───cc_lipidomics
        └───Lipidomics_RawData_2.csv
        └───Lipidomics_RawData.csv
    └───pfc_lipidomics
        └───ChE_summary_cyc_05342022_all_samples.csv
        └───metadata_PFC_all_individuals_092520.tsv
        └───ROSMAP_clinical.csv
        └───ROSMAP_Lipidomics_Emory_biospecimen_metadata.csv
└───other_analyses_outputs
    └───pathways.rds
    └───pathway_scores.rds    
└───supplementary_tables
    └───...
```
2. here is some info regarding the origins of each piece of data:

| Data File                                                      | Origin                                                              |
|-|-
| <td colspan=2>single_cell_data                                                                                                       |         
|----------------------------------------------------------------|---------------------------------------------------------------------|
| Cell_group_colors.rds                                          | NA                                                                  |
| ensembl.GRCh38p12.genes.complete.annot.rds                     |                                                                     |
| single_cell_experiment_object.rds                              | run ../scripts/run_qc_and_annotation.r                              |
| expressed_genes_per_celltype.rds                               | run ../scripts/run_qc_and_annotation.r                              |
| individual_level_averages_per_celltype.rds                     | run ../scripts/run_qc_and_annotation.r                              |
| Mapping                                                        |                                                                     |
| metadata_by_individual.csv                                     | provided by ROSMAP                                                  |
| metadata_PFC_all_individuals_092520.tsv                        | provided by ROSMAP                                                  |
| Metadata.APOE.project.rds                                      | provided by ROSMAP                                                  |
| RefCellTypeMarkers.adultBrain.rds                              | provided by ROSMAP                                                  |
| Summary.data.celltype.rds                                      | provided by ROSMAP                                                  |
| PanglaoDB.by.organ.by.celltype.rds                             |                                                                     |
| nebula_oli_degs.rds                                            | run ../scripts/get_nebula_degs.r                                    |
| wilcoxon_degs_all.rds                                          | run ../scripts/get_wilcox_degs.r                                    |
| wilcoxon_degs_AD.rds                                           | run ../scripts/get_wilcox_degs.r                                    |
| wilcoxon_degs_noAD.rds                                         | run ../scripts/get_wilcox_degs.r                                    |
| pseudo_bulk_degs_single_cell_all_celltypes                     | run ../scripts/get_pseudobulk_degs.r                                |
| FPKM_tables_.txt                                               | provided by MIT biomicro center core facility                       |
| ipsc_bulk_rnaseq_count_files/                                  | provided by MIT biomicro center core facility                       |
| ipsc_metadata.rds                                              | NA                                                                  |
| GO_Biological_Process_2018.txt                                 | from mayaan lab  [here](https://maayanlab.cloud/Enrichr/#libraries) |
| HumanCyc_2016.txt                                              | from mayaan lab [here](https://maayanlab.cloud/Enrichr/#libraries)  |
| KEGG_2019_Human.txt                                            | from mayaan lab [here](https://maayanlab.cloud/Enrichr/#libraries)  |
| Reactome_2016.txt                                              | from mayaan lab [here](https://maayanlab.cloud/Enrichr/#libraries)  |
| cc_lipidomics_values.rds                                       | provided by _                                                       |
| pfc_lipidomics_values/ChE_summary_cyc_05342022_all_samples.csv | provided by Emory University                                        |
| pfc_lipidomics_values/ROSMAP_clinical.csv                      | on synapse [here](https://www.synapse.org/#!Synapse:syn21682218)    |
| pfc_lipidomics_values/ROSMAP_Lipidomics_Emory_.._metadata.csv  | on Synapse [here](https://www.synapse.org/#!Synapse:syn26452615)    |
| processed_pathways.rds                                         | run ../scripts/get_pathways.r                                       |
| processed_pathway_fits.rds                                     | run ../scripts/pathway_analyses.r                                   |
| cholesterol_and_er_gene_level_results.rds                      | run ../scripts/dissecting_cholesterol_er_dysregulation.r            |
| lipidomics_pfc_wilcoxon_results.rds                            | run ../scripts/lipidomic_analysis_pfc.r                             |
| lipidomics_cc_wilcoxon_results.rds                             | run ../scripts/lipidomic_analysis_cc.r                              |
| pathways.rds                                                   | run ../scripts/get_pathways.r                                       |
| pathway_scores.rds                                             | run ../scripts/pathway_analyses.r                                   |  
| cholesterol_analysis.rds                                       | run ../scripts/dissecting_cholesterol_dysregulation.r               |
| cc_lipidomics_data.rds                                         | run ../scripts/lipidomic_analysis_cc.r                              |
| stratified_stats.csv                                           | run ../scripts/e4_effects_stratification_by_AD.r                    |
| stratified_anaylsis.rds                                        | run ../scripts/e4_effects_stratification_by_AD.r                    |
| pfc_lipidomics_data.rds                                        | run ../scripts/lipidomic_analysis_pfc.r                             |
| pfc_lipidomics_data_metadata.csv                               | run ../scripts/lipidomic_analysis_pfc.r                             |
| ipsc_post_mortem_pca_var_explained.csv                         | run ../scripts/comparison_of_ipsc_and_brain.r                       |
| scaled_expression_ipsc_post_mortem.csv                         | run ../scripts/comparison_of_ipsc_and_brain.r                       |
| APOE_expression_post_mortem_oligos.csv                         | run ../scripts/APOE_expression_oligodendrocytes.r                   |
| APOE_expression_post_mortem_oligos.csv                         | run ../scripts/APOE_expression_oligodendrocytes.r                   |
| er_stress_results.readRDS                                      | run ../scripts/get_postmortem_er_stress_pathways.r                  |

3. Now follow these steps to recapitulate the analysis:

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
Rscript ../scripts/plots_for_extended_data_figure_1.r #TODO: check / edit this
```

Extended Data Figure 2
```bash
conda activate apoe_env
Rscript ../scripts/fgsea_analysis.r # edit which input data using
Rscript ../scripts/pseudo_bulk.r
Rscript ../scripts/plots_for_extended_data_figure2.r #TODO: check / edit this
```

Extended Data Figure 3
```bash
conda activate apoe_env
Rscript ../scripts/e4_effects_stratification_by_AD.r
```

Extended Data Figure 4
```bash
conda activate apoe_env
Rscript ../scripts/lipidomic_analysis_pfc.r
```

Extended Data Figure 6
```bash
conda activate apoe_env
Rscript ../scripts/comparison_of_ipsc_and_brain.r
Rscript ../scripts/APOE_expression_oligodendrocytes.r
Rscript ../scripts/apoe_expression_ipsc.r
```

Extended Data Figure 8
```bash
conda activate apoe_env
Rscript ../scripts/get_postmortem_er_stress_pathways.r
Rscript ../scripts/er_postmortem_plots.r
Rscript ../scripts/ipsc_pathway_perturbations.r
```

Extended Data Figure 9
```bash
conda activate apoe_env
Rscript ../scripts/get_wilcox_myelin_plots.r
```
