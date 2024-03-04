# Structural connectome in depression 
<p>Scripts used for the following paper: "A comprehensive hierarchical comparison of structural connectomes in Major Depressive Disorder cases versus controls in two large population samples"</p>

## prep directory
* prep_conn_matrices.R: extract FA matrices for each subject and compile them into a 3D matrix (85 nodes x 85 nodes x no_of_subj) for cases and controls
* prep_subj_phenotypes_gs.R: prepare MDD phenotype and other necessary covariates in GS
* prep_subj_phenotypes_ukb.R: prepare MDD phenotype and other necessary covariates in UKB

## analysis directory 
* analy_richclub_calculation_and_permutation.R: derive group-averaged network for cases and controls, conduct proportional thresholding to remove spurious edges, calculate rich club coefficient for each group, and conduct permutation testing to determine significance of (1) rich club coefficient for each group, and (2) case-control differences
* analy_network_measures_calculation_and_permutation.R: conduct proportional thresholding for each individual, derive global, tier and nodal network measures, conduct permutation testing to determine significance of network measures, and conduct correction for multiple comparisons using the hierarchical FDR approach
* analy_compare_ukb_gs.R: to compare the case-control Cohen's d derived in UKB versus those derived in GS for all network measures

## plotting directory 
* plot_global_tier_cohend_fig2.R: to plot bargraphs to reflect case-control Cohen's d for the global and tier network measures (see Fig 2 in paper)
* plot_richclub_curves_fig3.R: to plot rich club curves with two axes to represent original and normalised rich club coefficient using plotly (see Fig 3 in paper)
* plot_nodal_ggseg_pplot_fig4.R: to plot case-control Cohen's d for nodal measures onto brain segmentation map, and plot p-values for significant measures (see Fig 4 in paper)

## functions directory 
(for large functions; smaller functions are within the scripts) 
* richclub_w_lists_FUN.R: to determine rich club coefficient for MDD cases and healthy controls
* richclub_permutation_lists_FUN.R: to determine significance of rich club coefficients
* nodal_measures_indiv_threshprop_FUN.R: used to determine network measures for each individual 
* nodal_tier_FUN.R: used to determine tiers and tier-based network measures
* nodal_cohend_wperm_FUN.R: to derive case-control Cohen's d for the network measures and to determine significance through permutation testing
* hfdr_FUN.R: to conduct hierarchical FDR correction in the order of global, tier and nodal level

