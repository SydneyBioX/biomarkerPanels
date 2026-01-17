#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <cmath>
using namespace Rcpp;

//' @title Compute stable gene scores (C++)
//' @param abs_t Matrix of absolute t-statistics (genes x cohorts)
//' @param se_matrix Matrix of standard errors (genes x cohorts)
//' @param method Scoring method: "precision_weighted", "cv_t_stats", or "inverse_t_se"
//' @return Numeric vector of scores (one per gene)
//' @keywords internal
// [[Rcpp::export(.select_stable_genes_cpp)]]
NumericVector select_stable_genes_cpp(NumericMatrix abs_t,
                                       NumericMatrix se_matrix,
                                       std::string method) {
  int n_genes = abs_t.nrow();
  int n_cohorts = abs_t.ncol();
  NumericVector scores(n_genes);

  for (int i = 0; i < n_genes; i++) {
    if (method == "precision_weighted") {
      // Weighted mean of 1/(t+0.01) with weights 1/(se+1e-8)
      double sum_weighted = 0.0;
      double sum_weights = 0.0;
      int n_valid = 0;

      for (int j = 0; j < n_cohorts; j++) {
        double t_val = abs_t(i, j);
        double se_val = se_matrix(i, j);

        if (!R_IsNA(t_val) && !R_IsNA(se_val) && R_finite(t_val) && R_finite(se_val)) {
          double precision = 1.0 / (se_val + 1e-8);
          if (precision > 0) {
            double val = 1.0 / (t_val + 0.01);
            sum_weighted += val * precision;
            sum_weights += precision;
            n_valid++;
          }
        }
      }

      if (n_valid == 0) {
        scores[i] = NA_REAL;
      } else {
        scores[i] = sum_weighted / sum_weights;
      }

    } else if (method == "cv_t_stats") {
      // 1 / ((sd + 1e-6) * (mean + 0.01))
      double sum = 0.0;
      double sum_sq = 0.0;
      int n_valid = 0;

      for (int j = 0; j < n_cohorts; j++) {
        double t_val = abs_t(i, j);
        if (!R_IsNA(t_val) && R_finite(t_val)) {
          sum += t_val;
          sum_sq += t_val * t_val;
          n_valid++;
        }
      }

      if (n_valid == 0) {
        scores[i] = NA_REAL;
      } else {
        double mean_abs = sum / n_valid;
        double sd_abs = 0.0;
        if (n_valid > 1) {
          double variance = (sum_sq - sum * sum / n_valid) / (n_valid - 1);
          sd_abs = variance > 0 ? std::sqrt(variance) : 0.0;
        }
        scores[i] = 1.0 / ((sd_abs + 1e-6) * (mean_abs + 0.01));
      }

    } else if (method == "inverse_t_se") {
      // mean((t+0.01) * (se+1e-8)), then invert at end
      double sum = 0.0;
      int n_valid = 0;

      for (int j = 0; j < n_cohorts; j++) {
        double t_val = abs_t(i, j);
        double se_val = se_matrix(i, j);
        if (!R_IsNA(t_val) && !R_IsNA(se_val) && R_finite(t_val) && R_finite(se_val)) {
          sum += (t_val + 0.01) * (se_val + 1e-8);
          n_valid++;
        }
      }

      if (n_valid == 0) {
        scores[i] = NA_REAL;
      } else {
        // Store mean first, will invert later
        scores[i] = sum / n_valid;
      }
    } else {
      scores[i] = NA_REAL;
    }
  }

  // For inverse_t_se, invert all scores at the end
  if (method == "inverse_t_se") {
    for (int i = 0; i < n_genes; i++) {
      if (!R_IsNA(scores[i])) {
        scores[i] = 1.0 / (scores[i] + 1e-6);
      }
    }
  }

  // Set non-finite values to NA
  for (int i = 0; i < n_genes; i++) {
    if (!R_finite(scores[i])) {
      scores[i] = NA_REAL;
    }
  }

  return scores;
}

//' @title Score transferable features (C++)
//' @param coefficient_matrix Matrix of coefficients (features x cohorts)
//' @return List with mean_abs, sd, min_abs, sign_agreement vectors
//' @keywords internal
// [[Rcpp::export(.score_transferable_features_cpp)]]
List score_transferable_features_cpp(NumericMatrix coefficient_matrix) {
  int n_features = coefficient_matrix.nrow();
  int n_cohorts = coefficient_matrix.ncol();

  NumericVector mean_abs(n_features);
  NumericVector sd_vals(n_features);
  NumericVector min_abs(n_features);
  NumericVector sign_agreement(n_features);

  for (int i = 0; i < n_features; i++) {
    // Compute mean of absolute values
    double sum_abs = 0.0;
    double sum = 0.0;
    double sum_sq = 0.0;
    double min_val = R_PosInf;
    int n_valid = 0;

    // Also track non-zero values for sign agreement
    std::vector<double> nz_vals;
    nz_vals.reserve(n_cohorts);

    for (int j = 0; j < n_cohorts; j++) {
      double val = coefficient_matrix(i, j);
      if (!R_IsNA(val) && R_finite(val)) {
        double abs_val = std::abs(val);
        sum_abs += abs_val;
        sum += val;
        sum_sq += val * val;
        if (abs_val < min_val) {
          min_val = abs_val;
        }
        n_valid++;

        if (abs_val > 1e-6) {
          nz_vals.push_back(val);
        }
      }
    }

    if (n_valid == 0) {
      mean_abs[i] = NA_REAL;
      sd_vals[i] = NA_REAL;
      min_abs[i] = NA_REAL;
      sign_agreement[i] = NA_REAL;
    } else {
      mean_abs[i] = sum_abs / n_valid;
      min_abs[i] = min_val;

      // Compute standard deviation of original values (not abs)
      if (n_valid > 1) {
        double variance = (sum_sq - sum * sum / n_valid) / (n_valid - 1);
        sd_vals[i] = variance > 0 ? std::sqrt(variance) : 0.0;
      } else {
        sd_vals[i] = NA_REAL;  // R's sd() returns NA for single value
      }

      // Compute sign agreement
      if (nz_vals.empty()) {
        sign_agreement[i] = 1.0;  // All zeros -> consistent
      } else {
        // Find median sign
        std::vector<double> sorted_nz = nz_vals;
        std::sort(sorted_nz.begin(), sorted_nz.end());
        double median_val;
        int n_nz = sorted_nz.size();
        if (n_nz % 2 == 0) {
          median_val = (sorted_nz[n_nz/2 - 1] + sorted_nz[n_nz/2]) / 2.0;
        } else {
          median_val = sorted_nz[n_nz/2];
        }
        int median_sign = (median_val > 0) ? 1 : ((median_val < 0) ? -1 : 0);

        // Count agreements
        int n_agree = 0;
        for (size_t k = 0; k < nz_vals.size(); k++) {
          int val_sign = (nz_vals[k] > 0) ? 1 : ((nz_vals[k] < 0) ? -1 : 0);
          if (val_sign == median_sign) {
            n_agree++;
          }
        }
        sign_agreement[i] = (double)n_agree / nz_vals.size();
      }
    }
  }

  return List::create(
    Named("mean_abs") = mean_abs,
    Named("sd") = sd_vals,
    Named("min_abs") = min_abs,
    Named("sign_agreement") = sign_agreement
  );
}

// Helper function for ranking (used in aggregate_de_pvalues)
std::vector<double> rank_vector(const std::vector<double>& vals) {
  int n = vals.size();
  std::vector<double> ranks(n);

  // Create index array
  std::vector<int> order(n);
  for (int i = 0; i < n; i++) order[i] = i;

  // Sort indices by values
  std::sort(order.begin(), order.end(), [&vals](int a, int b) {
    return vals[a] < vals[b];
  });

  // Assign ranks with tie handling (average)
  int i = 0;
  while (i < n) {
    int j = i;
    // Find all tied values
    while (j < n - 1 && vals[order[j]] == vals[order[j+1]]) {
      j++;
    }
    // Average rank for tied values
    double avg_rank = (i + j) / 2.0 + 1.0;  // +1 for 1-based ranks
    for (int k = i; k <= j; k++) {
      ranks[order[k]] = avg_rank;
    }
    i = j + 1;
  }

  return ranks;
}

//' @title Aggregate DE p-values across cohorts (C++)
//' @param t_matrix Matrix of t-statistics (genes x cohorts)
//' @param method Combination method: "Stouffer", "Fisher", "OSP", or "maxP"
//' @return Numeric vector of combined two-sided p-values (sorted ascending)
//' @keywords internal
// [[Rcpp::export(.aggregate_de_pvalues_cpp)]]
NumericVector aggregate_de_pvalues_cpp(NumericMatrix t_matrix, std::string method) {
  int n_genes = t_matrix.nrow();
  int n_cohorts = t_matrix.ncol();

  if (n_genes == 0 || n_cohorts == 0) {
    return NumericVector(0);
  }

  // Step 1: Convert to z-scores column by column
  NumericMatrix z_scores(n_genes, n_cohorts);

  for (int j = 0; j < n_cohorts; j++) {
    // Extract valid values for this column
    std::vector<double> valid_vals;
    std::vector<int> valid_idx;

    for (int i = 0; i < n_genes; i++) {
      double val = t_matrix(i, j);
      if (!R_IsNA(val) && R_finite(val)) {
        valid_vals.push_back(val);
        valid_idx.push_back(i);
      }
    }

    // Initialize all to NA
    for (int i = 0; i < n_genes; i++) {
      z_scores(i, j) = NA_REAL;
    }

    if (!valid_vals.empty()) {
      // Rank the valid values
      std::vector<double> ranks = rank_vector(valid_vals);

      // Convert ranks to probabilities and then to z-scores
      int n_valid = valid_vals.size();
      for (size_t k = 0; k < valid_idx.size(); k++) {
        double prob = ranks[k] / (n_valid + 1.0);
        z_scores(valid_idx[k], j) = R::qnorm(prob, 0.0, 1.0, 1, 0);
      }
    }
  }

  // Step 2: Combine p-values row by row
  NumericVector combined_p(n_genes);

  for (int i = 0; i < n_genes; i++) {
    // Collect valid z-scores for this row
    std::vector<double> valid_z;
    for (int j = 0; j < n_cohorts; j++) {
      double z = z_scores(i, j);
      if (!R_IsNA(z) && R_finite(z)) {
        valid_z.push_back(z);
      }
    }

    if (valid_z.empty()) {
      combined_p[i] = NA_REAL;
      continue;
    }

    double p;
    int n_z = valid_z.size();

    if (method == "Stouffer") {
      double sum_z = 0.0;
      for (int k = 0; k < n_z; k++) sum_z += valid_z[k];
      p = R::pnorm(sum_z, 0.0, std::sqrt((double)n_z), 0, 0);  // lower.tail = FALSE
    } else if (method == "OSP") {
      double sum_log_p = 0.0;
      for (int k = 0; k < n_z; k++) {
        double p_val = R::pnorm(valid_z[k], 0.0, 1.0, 1, 0);  // lower.tail = TRUE
        sum_log_p += std::log(p_val);
      }
      p = R::pchisq(-2.0 * sum_log_p, 2.0 * n_z, 1, 0);  // lower.tail = TRUE
    } else if (method == "Fisher") {
      double sum_log_p = 0.0;
      for (int k = 0; k < n_z; k++) {
        double p_val = R::pnorm(valid_z[k], 0.0, 1.0, 0, 0);  // lower.tail = FALSE
        sum_log_p += std::log(p_val);
      }
      p = R::pchisq(-2.0 * sum_log_p, 2.0 * n_z, 0, 0);  // lower.tail = FALSE
    } else if (method == "maxP") {
      double max_z = valid_z[0];
      for (int k = 1; k < n_z; k++) {
        if (valid_z[k] > max_z) max_z = valid_z[k];
      }
      p = R::pnorm(max_z, 0.0, 1.0, 0, 0);  // lower.tail = FALSE
    } else {
      p = NA_REAL;
    }

    // Clamp to minimum double value
    if (p < 2.225074e-308) {
      p = 2.225074e-308;
    }

    combined_p[i] = p;
  }

  // Step 3: Convert to two-sided p-values
  NumericVector two_sided(n_genes);
  for (int i = 0; i < n_genes; i++) {
    if (R_IsNA(combined_p[i])) {
      two_sided[i] = NA_REAL;
    } else {
      double combined_z = R::qnorm(combined_p[i], 0.0, 1.0, 0, 0);
      two_sided[i] = 2.0 * R::pnorm(-std::abs(combined_z), 0.0, 1.0, 1, 0);
    }
  }

  // Filter out NA values and sort
  std::vector<std::pair<double, int>> valid_results;
  CharacterVector row_names = rownames(t_matrix);

  for (int i = 0; i < n_genes; i++) {
    if (!R_IsNA(two_sided[i])) {
      valid_results.push_back(std::make_pair(two_sided[i], i));
    }
  }

  // Sort by p-value ascending
  std::sort(valid_results.begin(), valid_results.end());

  // Create output vector with names
  int n_valid = valid_results.size();
  NumericVector result(n_valid);
  CharacterVector result_names(n_valid);

  for (int i = 0; i < n_valid; i++) {
    result[i] = valid_results[i].first;
    if (row_names.size() > 0) {
      result_names[i] = row_names[valid_results[i].second];
    }
  }

  if (row_names.size() > 0) {
    result.names() = result_names;
  }

  return result;
}
