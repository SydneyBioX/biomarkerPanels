# The metric functions are internal (reached via registry names in user code,
# e.g. define_objectives(metrics = "auc")), but the unit tests exercise them
# directly. Alias them here so tests work both under devtools::load_all() and
# against the installed namespace during R CMD check.
metric_sensitivity <- biomarkerPanels:::metric_sensitivity
metric_specificity <- biomarkerPanels:::metric_specificity
metric_precision <- biomarkerPanels:::metric_precision
metric_npv <- biomarkerPanels:::metric_npv
metric_f1 <- biomarkerPanels:::metric_f1
metric_balanced_accuracy <- biomarkerPanels:::metric_balanced_accuracy
metric_auc <- biomarkerPanels:::metric_auc
metric_pauc <- biomarkerPanels:::metric_pauc
metric_num_features <- biomarkerPanels:::metric_num_features
metric_sensitivity_at_specificity <- biomarkerPanels:::metric_sensitivity_at_specificity
metric_min_cohort_auc <- biomarkerPanels:::metric_min_cohort_auc
metric_cohort_auc_gap <- biomarkerPanels:::metric_cohort_auc_gap
metric_cohort_auc_var <- biomarkerPanels:::metric_cohort_auc_var
metric_max_cohort_brier <- biomarkerPanels:::metric_max_cohort_brier
metric_cohort_leakage <- biomarkerPanels:::metric_cohort_leakage
