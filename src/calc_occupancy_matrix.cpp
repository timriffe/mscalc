#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix calc_occupancy_matrix_cpp(NumericMatrix p_mat,
                                        CharacterVector transient_states,
                                        CharacterVector all_states,
                                        NumericVector init_probs,
                                        int n,
                                        std::string delim = "->") {

  int n_transient = transient_states.size();
  NumericMatrix state_occup(n, n_transient);
  CharacterVector colnames = transient_states;
  state_occup(0, _) = init_probs;

  // Grab colnames of matrix (should be transitions like "P->W", etc.)
  List dimnames = p_mat.attr("dimnames");
  CharacterVector mat_colnames = dimnames[1];
  std::vector<std::string> trans_names = as<std::vector<std::string>>(mat_colnames);
  std::vector<std::string> from_vec = as<std::vector<std::string>>(transient_states);
  std::vector<std::string> to_vec = as<std::vector<std::string>>(all_states);

  for (int i = 0; i < n - 1; i++) {
    for (int tf = 0; tf < n_transient; tf++) {
      std::string from = from_vec[tf];
      double occup_from = state_occup(i, tf);
      double outflow = 0.0;

      for (const std::string& to : to_vec) {
        std::string trans_name = from + delim + to;

        // Find column in p_mat
        auto it = std::find(trans_names.begin(), trans_names.end(), trans_name);
        double p = 0.0;

        if (it != trans_names.end()) {
          int col_idx = std::distance(trans_names.begin(), it);
          p = p_mat(i, col_idx);
        }

        // Update occupancy if destination is transient
        auto to_it = std::find(from_vec.begin(), from_vec.end(), to);
        if (to_it != from_vec.end()) {
          int to_col = std::distance(from_vec.begin(), to_it);
          state_occup(i + 1, to_col) += occup_from * p;
        }

        outflow += p;
      }

      // Residual self-transition
      double stay_prob = std::max(0.0, 1.0 - outflow);
      state_occup(i + 1, tf) += occup_from * stay_prob;
    }
  }

  state_occup.attr("dimnames") = List::create(R_NilValue, transient_states);
  return state_occup;
}
