#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix calc_occupancy_cpp(List p_list,
                                 CharacterVector transient_states,
                                 CharacterVector all_states,
                                 NumericVector init_probs,
                                 int n,
                                 std::string delim = "->") {

  int n_transient = transient_states.size();
  NumericMatrix state_occup(n, n_transient);
  CharacterVector colnames = transient_states;
  state_occup(0, _) = init_probs;

  // Caching names for speed
  std::vector<std::string> from_vec = as<std::vector<std::string>>(transient_states);
  std::vector<std::string> to_vec = as<std::vector<std::string>>(all_states);

  for (int i = 0; i < n - 1; i++) {
    for (int tf = 0; tf < n_transient; tf++) {
      std::string from = from_vec[tf];
      double occup_from = state_occup(i, tf);
      double outflow = 0.0;

      for (std::string to : to_vec) {
        std::string trans_name = from + delim + to;
        NumericVector trans_p = p_list.containsElementNamed(trans_name.c_str()) ?
        as<NumericVector>(p_list[trans_name]) :
          NumericVector(n, 0.0);
        double p = trans_p[i];

        if (std::find(from_vec.begin(), from_vec.end(), to) != from_vec.end()) {
          int to_col = std::distance(from_vec.begin(), std::find(from_vec.begin(), from_vec.end(), to));
          state_occup(i + 1, to_col) += occup_from * p;
        }
        outflow += p;
      }

      double stay_prob = std::max(0.0, 1.0 - outflow);
      state_occup(i + 1, tf) += occup_from * stay_prob;
    }
  }

  state_occup.attr("dimnames") = List::create(R_NilValue, transient_states);

  return state_occup;
}

