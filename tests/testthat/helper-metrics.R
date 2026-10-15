# The metric functions are internal (reached via registry names in user code,
# e.g. define_objectives(metrics = "auc")), but the unit tests exercise them
# directly. Alias them here so tests work both under devtools::load_all() and
# against the installed namespace during R CMD check.
metric_sensitivity <- MOBID:::metric_sensitivity
metric_specificity <- MOBID:::metric_specificity
metric_precision <- MOBID:::metric_precision
metric_npv <- MOBID:::metric_npv
metric_f1 <- MOBID:::metric_f1
metric_balanced_accuracy <- MOBID:::metric_balanced_accuracy
metric_auc <- MOBID:::metric_auc
metric_pauc <- MOBID:::metric_pauc
metric_num_features <- MOBID:::metric_num_features
metric_sensitivity_at_specificity <- MOBID:::metric_sensitivity_at_specificity
metric_min_cohort_auc <- MOBID:::metric_min_cohort_auc
metric_cohort_auc_gap <- MOBID:::metric_cohort_auc_gap
metric_cohort_auc_var <- MOBID:::metric_cohort_auc_var
metric_max_cohort_brier <- MOBID:::metric_max_cohort_brier
metric_cohort_leakage <- MOBID:::metric_cohort_leakage
