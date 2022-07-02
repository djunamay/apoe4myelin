#!/bin/bash

conda activate apoe_env

# Figure 1
echo '*** Running Figure 1 analysis...'
Rscript ../scripts/get_pathways.r
Rscript ../scripts/get_nebula_degs.r
Rscript ../scripts/pathway_analyses.r
Rscript ../scripts/get_figure_1_plots.r
#Figure 2
echo '*** Running Figure 2 analysis...'
Rscript ../scripts/dissecting_cholesterol_dysregulation.r
Rscript ../scripts/lipidomic_analysis_cc.r
Rscript ../scripts/get_figure_2_plots.r
#Figure 3
echo '*** Running Figure 3 analysis...'
Rscript ../scripts/get_wilcox_degs.r
Rscript ../scripts/get_figure_3_plots.r
# Extended Data Figure 1
echo '*** Making Extended Data Figure 1 plots...'
Rscript ../scripts/plots_for_extended_data_figure_1.r
#Extended Data Figure 2
echo '*** Making Extended Data Figure 2 plots...'
Rscript ../scripts/fgsea_analysis.r
Rscript ../scripts/pseudo_bulk.r
Rscript ../scripts/plots_for_extended_data_figure2.r
#Extended Data Figure 3
echo '*** Making Extended Data Figure 3 plots...'
Rscript ../scripts/e4_effects_stratification_by_AD.r
Rscript ../scripts/e4_stratification_plots.r
#Extended Data Figure 4
echo '*** Making Extended Data Figure 4 plots...'
Rscript ../scripts/lipidomic_analysis_pfc.r
#Extended Data Figure 6
echo '*** Making Extended Data Figure 6 plots...'
Rscript ../scripts/comparison_of_ipsc_and_brain_get_scaled_matrices.r
Rscript ../scripts/comparison_of_ipsc_and_brain_plots.r
Rscript ../scripts/APOE_expression_oligodendrocytes.r
Rscript ../scripts/apoe_expression_ipsc.r
#Extended Data Figure 8
echo '*** Making Extended Data Figure 8 plots...'
Rscript ../scripts/get_postmortem_er_stress_pathways.r
Rscript ../scripts/er_postmortem_plots.r
Rscript ../scripts/ipsc_gene_perturbations.r
#Extended Data Figure 9
echo '*** Making Extended Data Figure 9 plots...'
Rscript ../scripts/get_wilcox_myelin_plots.r

ehco 'done.'
