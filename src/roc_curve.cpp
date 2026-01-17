#include <Rcpp.h>
#include <algorithm>
#include <vector>
using namespace Rcpp;

//' @title Compute ROC curve efficiently (C++)
//' @param scores Numeric vector of predicted scores
//' @param is_positive Logical vector indicating positive class membership
//' @return DataFrame with columns: threshold, tpr, fpr
//' @keywords internal
// [[Rcpp::export(.compute_roc_curve_cpp)]]
DataFrame compute_roc_curve_cpp(NumericVector scores, LogicalVector is_positive) {
  int n = scores.size();
  if (n != is_positive.size()) {
    stop("scores and is_positive must have the same length");
  }

  // Count total positives and negatives
  int pos_total = 0;
  int neg_total = 0;
  for (int i = 0; i < n; i++) {
    if (is_positive[i]) {
      pos_total++;
    } else {
      neg_total++;
    }
  }

  // Create index array and sort by scores descending
  std::vector<int> order(n);
  for (int i = 0; i < n; i++) {
    order[i] = i;
  }
  std::sort(order.begin(), order.end(), [&scores](int a, int b) {
    return scores[a] > scores[b];
  });

  // Get unique thresholds (sorted descending)
  std::vector<double> unique_thresholds;
  unique_thresholds.reserve(n + 2);
  unique_thresholds.push_back(R_PosInf);  // Threshold above all scores

  double last_score = R_PosInf;
  for (int i = 0; i < n; i++) {
    double s = scores[order[i]];
    if (s != last_score) {
      unique_thresholds.push_back(s);
      last_score = s;
    }
  }
  unique_thresholds.push_back(R_NegInf);  // Threshold below all scores

  int n_thresholds = unique_thresholds.size();

  // Pre-allocate output vectors
  NumericVector out_threshold(n_thresholds);
  NumericVector out_tpr(n_thresholds);
  NumericVector out_fpr(n_thresholds);

  // Walk through thresholds, accumulating counts
  // At threshold t, we predict positive for all scores >= t
  int tp = 0;
  int fp = 0;
  int score_idx = 0;

  for (int t = 0; t < n_thresholds; t++) {
    double thresh = unique_thresholds[t];
    out_threshold[t] = thresh;

    // Accumulate all samples with score >= thresh
    while (score_idx < n && scores[order[score_idx]] >= thresh) {
      if (is_positive[order[score_idx]]) {
        tp++;
      } else {
        fp++;
      }
      score_idx++;
    }

    // Compute rates
    if (pos_total == 0) {
      out_tpr[t] = NA_REAL;
    } else {
      out_tpr[t] = (double)tp / pos_total;
    }

    if (neg_total == 0) {
      out_fpr[t] = NA_REAL;
    } else {
      out_fpr[t] = (double)fp / neg_total;
    }
  }

  // Create output DataFrame
  DataFrame result = DataFrame::create(
    Named("threshold") = out_threshold,
    Named("tpr") = out_tpr,
    Named("fpr") = out_fpr,
    Named("sensitivity") = out_tpr,
    Named("specificity") = 1.0 - out_fpr
  );

  return result;
}
