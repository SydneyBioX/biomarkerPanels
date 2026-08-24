#include <Rcpp.h>
using namespace Rcpp;

//' @title Compute pairwise difference between matrix columns (C++)
//' @param x A numeric matrix
//' @param col_names Character vector of column names
//' @return A matrix of pairwise differences
//' @noRd
// [[Rcpp::export(.pairwise_col_diff_cpp)]]
 NumericMatrix pairwise_col_diff_cpp(NumericMatrix x, CharacterVector col_names) {
   int n = x.nrow();
   int p = x.ncol();
   // Use long long to avoid integer overflow for large p
   // (p * (p-1) / 2 can exceed INT_MAX when p > ~46,340)
   long long n_pairs_ll = (static_cast<long long>(p) * (p - 1)) / 2;
   if (n_pairs_ll > INT_MAX) {
     Rcpp::stop("Too many feature pairs for pairwise aggregation");
   }
   int n_pairs = static_cast<int>(n_pairs_ll);

   // Pre-allocate result
   NumericMatrix result(n, n_pairs);
   CharacterVector result_names(n_pairs);

   int idx = 0;
   for (int i = 0; i < p - 1; i++) {
     for (int j = i + 1; j < p; j++) {
       // Compute difference for this pair
       for (int row = 0; row < n; row++) {
         result(row, idx) = x(row, i) - x(row, j);
       }
       // Create column name using std::string for concatenation
       std::string name1 = Rcpp::as<std::string>(col_names[i]);
       std::string name2 = Rcpp::as<std::string>(col_names[j]);
       result_names[idx] = name1 + "--" + name2;
       idx++;
     }
   }

   colnames(result) = result_names;
   return result;
 }

//' @title Compute differences relative to a reference column (C++)
//' @param x A numeric matrix (excluding reference column)
//' @param ref_col The reference column vector
//' @param ref_name Name of the reference feature
//' @param other_names Names of other columns
//' @return A matrix of differences (other - ref)
//' @noRd
// [[Rcpp::export(.reference_diff_cpp)]]
NumericMatrix reference_diff_cpp(NumericMatrix x, NumericVector ref_col,
                                 String ref_name, CharacterVector other_names) {
  int n = x.nrow();
  int p = x.ncol();

  NumericMatrix result(n, p);
  CharacterVector result_names(p);

  std::string ref_str = ref_name;

  for (int j = 0; j < p; j++) {
    for (int i = 0; i < n; i++) {
      result(i, j) = x(i, j) - ref_col[i];
    }
    std::string other_str = Rcpp::as<std::string>(other_names[j]);
    result_names[j] = other_str + "--" + ref_str;
  }

  colnames(result) = result_names;
  return result;
}
