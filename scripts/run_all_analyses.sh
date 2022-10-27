#!/bin/bash

conda activate apoe4myelin_env
echo '*** Running pre - analysis scripts...'
echo '***** generating SCE object...'
Rscript ../scripts/generate_qc_sce_object.r
conda activate actionet_legacy_env
echo '***** getting individual level averages...'
Rscript ../scripts/get_individual_level_averages_object.r
echo '***** curating pathways...'
conda activate apoe4myelin_env
Rscript ../scripts/get_pathways.r
echo '***** computing degs with nebula...'
Rscript ../scripts/nebula_degs.r
echo '***** getting expressed genes per celltype...'
Rscript ../scripts/get_expressed_genes_per_celltype.r
echo '***** getting ipsc opc counts table...'
Rscript ../scripts/get_opc_ipsc_counts_table.r
echo '***** making metadata file...'
Rscript ../scripts/make_metadata_file.r

# Figure 1
echo '*** Running scripts for figure 1...'
echo '***** running pathway analyses...'
Rscript ../scripts/pathway_analyses.r
echo '***** making plots for figure 1...'
Rscript ../scripts/get_figure_1_plots.r

# Figure 2
echo '*** Running scripts for figure 2...'
echo '***** dissecting cholesterol dysregulation...'
Rscript ../scripts/dissecting_cholesterol_dysregulation.r
echo '***** performing lipidomic analyses...'
Rscript ../scripts/lipidomic_analysis_cc.r
echo '***** making plots for figure 2...'
Rscript ../scripts/get_figure_2_plots.r

# Figure 4
echo '*** Running scripts for figure 4...'
echo '***** computing degs using wilcoxon...'
Rscript ../scripts/get_wilcox_degs.r
echo '***** making plots for figure 4...'
Rscript ../scripts/get_figure_4_plots.r

# Extended 1
echo '*** Running scripts for extended 1...'
Rscript ../scripts/plots_for_extended_data_figure_1.r

# Extended 2
echo '*** Running scripts for extended 2...'
Rscript ../scripts/fgsea_analysis.r
Rscript ../scripts/pseudo_bulk.r
Rscript ../scripts/plots_for_extended_data_figure2.r
Rscript ../scripts/e4_effects_stratification_by_AD.r
Rscript ../scripts/e4_stratification_plots.r

# Extended 3
echo '*** Running scripts for extended 3...'
Rscript ../scripts/lipidomic_analysis_pfc_get_data.r
Rscript ../scripts/lipidomic_analysis_pfc_make_plots.r

# Extended 5
echo '*** Running scripts for extended 5..'
Rscript ../scripts/comparison_of_ipsc_and_brain_get_scaled_matrices.r
Rscript ../scripts/comparison_of_ipsc_and_brain_plots.r
Rscript ../scripts/APOE_expression_oligodendrocytes.r
Rscript ../scripts/apoe_expression_ipsc.r

# Extended 7
echo '*** Running scripts for extended 7..'
Rscript ../scripts/get_postmortem_er_stress_pathways.r
Rscript ../scripts/er_postmortem_plots.r
Rscript ../scripts/ipsc_gene_perturbations.r

# Extended 8
echo '*** Running scripts for extended 8..'
Rscript ../scripts/get_wilcox_myelin_plots.r

ehco 'done.'
echo '*** See /plots directory for output plots..'

