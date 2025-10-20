#include <Rcpp.h>
using namespace Rcpp;

//' @title Compute pairwise difference between matrix columns (C++)
 //' @param x A numeric matrix
 //' @param col_names Character vector of column names
 //' @return A matrix of pairwise differences
 //' @keywords internal
 // [[Rcpp::export(.pairwise_col_diff_cpp)]]
 NumericMatrix pairwise_col_diff_cpp(NumericMatrix x, CharacterVector col_names) {
   int n = x.nrow();
   int p = x.ncol();
   int n_pairs = (p * (p - 1)) / 2;
   
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

//' @title Compute pairwise ratios to a specific feature (C++)
 //' @param x A numeric matrix
 //' @param feature_col The feature column vector
 //' @param feature_name Name of the feature
 //' @param other_names Names of other columns
 //' @return A matrix of pairwise ratios
 //' @keywords internal
 // [[Rcpp::export(.pairwise_ratios_cpp)]]
 NumericMatrix pairwise_ratios_cpp(NumericMatrix x, NumericVector feature_col, 
                                   String feature_name, CharacterVector other_names) {
   int n = x.nrow();
   int p = x.ncol();
   
   NumericMatrix result(n, p);
   CharacterVector result_names(p);
   
   // Convert feature_name to std::string once
   std::string feat_str = feature_name;
   
   for (int j = 0; j < p; j++) {
     for (int i = 0; i < n; i++) {
       result(i, j) = feature_col[i] / x(i, j);
     }
     // Create column name using std::string for concatenation
     std::string other_str = Rcpp::as<std::string>(other_names[j]);
     result_names[j] = feat_str + "/" + other_str;
   }
   
   colnames(result) = result_names;
   return result;
 }