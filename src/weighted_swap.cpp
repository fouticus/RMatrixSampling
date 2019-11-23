#include <utility>
#include <random>
#include <unordered_set>
#include <Rcpp.h>
using namespace Rcpp;

// Structs for creating a hash set of a pair
struct pair_hash{
  std::size_t operator()(std::pair<int,int> const &v) const {
    return std::hash<int>()(v.first) ^ std::hash<int>()(v.second);
  }
};
struct pair_equal {
  bool operator()(const std::pair<int,int> & p1, const std::pair<int,int> & p2) const {
    return (p1.first == p2.first) && (p1.second == p2.second);
  }
};

// [[Rcpp::plugins("cpp11")]]

/*
Perform repeated checkerboard swaps in an adjacency matrix of a simple directed graph

This swapping produces a draws from the uniform distribution of all directed graphs
with the same row and column sums.

@Re_from: list of originating vertices for each edge
@Re_to: list of terminating vertices for each edge
@n: number of swaps to perform
@swap_p: base probability of performing a swap (need swap_p < 1 for aperiodicity).
@seed: optional seed for random number generators
*/
// [[Rcpp::export]]
List simple_swap_n(SEXP Re_from, SEXP Re_to, int n, double swap_p, double seed=-1) {
  // Create cpp versions of R objects
  Rcpp::NumericVector e_from(Re_from);
  Rcpp::NumericVector e_to(Re_to);
  int m = e_from.size();  // number of edges
  // Make random number generators
  if(seed==-1){
      std::random_device rd;  // used to create a seed
      seed = rd();
  }
  std::mt19937 gen(seed);  // random number generator
  std::uniform_int_distribution<> randint(0, m-1);  // for selecting an edge
  std::uniform_real_distribution<> randu(0.0, 1.0);  // for deciding whether to swap
  double u;
  // Build hashmap for edges
  std::unordered_set< std::pair<int,int>, pair_hash, pair_equal> edgeSet;
  for (int i=0; i<m; i++){
    edgeSet.insert(std::make_pair(e_from[i], e_to[i]));
  }
  // initialize stuff before swapping
  int i;
  int j;
  int e_from_i;
  int e_to_i;
  int e_from_j;
  int e_to_j;
  // perform swaps
  for (int k=0; k<n; k++){
    // select two edges
    i = randint(gen);
    j = randint(gen);
    // make sure i and j are not the same
    if (i == j){
      continue;
    }
    e_from_i = e_from[i];
    e_to_i = e_to[i];
    e_from_j = e_from[j];
    e_to_j = e_to[j];
    // check that it's a checkerboard (0's in anti-diagonal)
    if (edgeSet.find(std::make_pair(e_from_i, e_to_j)) == edgeSet.end() &&
        edgeSet.find(std::make_pair(e_from_j, e_to_i)) == edgeSet.end() ){
      // it is a checkerboard, generate random uniform number
      u = randu(gen);
      if(u < swap_p){
        // make the swap in the hash set
        edgeSet.insert(std::make_pair(e_from_i, e_to_j));
        edgeSet.insert(std::make_pair(e_from_j, e_to_i));
        edgeSet.erase(std::make_pair(e_from_i, e_to_i));
        edgeSet.erase(std::make_pair(e_from_j, e_to_j));
        // make the swap in the edge lists
        e_from[i] = e_from_j;
        e_from[j] = e_from_i;
      }
    }
  }
  List edges;
  edges["from"] = e_from;
  edges["to"] = e_to;
  return edges;
}


/*
Perform repeated checkerboard swaps in an adjacency matrix of a simple directed graph,
respecting structural zeros and swapping according to weights.

This swapping produces a draws from the non-uniform distribution of all directed graphs
with the same row and column sums, with probabilities determined by the edge weights.

@Re_from: list of originating vertices for each edge
@Re_to: list of terminating vertices for each edge
@n: number of swaps to perform
@Rz_from: list of originating vertices for each structural zero
@Rz_to: list of terminating vertices for each structural zero
@swap_p: base probability of performing a swap (need swap_p < 1 for aperiodicity).
@seed: optional seed for random number generators
*/
// [[Rcpp::export]]
List swap_n(SEXP Re_from, SEXP Re_to, int n, SEXP Rw, SEXP Rz_from, SEXP Rz_to, double seed=-1) {
  // Create cpp versions of R objects
  Rcpp::NumericVector e_from(Re_from); // edge tails
  Rcpp::NumericVector e_to(Re_to);  // edge heads
  Rcpp::NumericVector z_from(Rz_from); // forbidden edge tails
  Rcpp::NumericVector z_to(Rz_to);  // forbidden edge heads
  Rcpp::NumericMatrix w(Rw);  // matrix of edge weights
  int m = e_from.size();  // number of edges
  // Make random number generators
  if(seed==-1){
      std::random_device rd;  // used to create a seed
      seed = rd();
  }
  std::mt19937 gen(seed);  // random number generator
  std::uniform_int_distribution<> randint(0, m-1);  // for selecting an edge
  std::uniform_real_distribution<> randu(0.0, 1.0);  // for deciding whether to swap
  double u;
  //Rcpp::rcout << "Building edge hashmap\n";
  // Build hashmap for edges
  std::unordered_set< std::pair<int,int>, pair_hash, pair_equal> edgeSet;
  for (int i=0; i<m; i++){
    edgeSet.insert(std::make_pair(e_from[i], e_to[i]));
  }
  //Rcpp::Rcout << "Building zeros hashmap\n";
  // Build hashmap for struct zeros
  std::unordered_set< std::pair<int,int>, pair_hash, pair_equal> zeroSet;
  for (int i=0; i<m; i++){
    zeroSet.insert(std::make_pair(z_from[i], z_to[i]));
  }
  // perform swaps
  //Rcpp::Rcout << "Allocating variables\n";
  int i;
  int j;
  int e_from_i;
  int e_to_i;
  int e_from_j;
  int e_to_j;
  double w_pre;
  double w_post;
  double swap_p;
  Rcpp::LogicalVector same_edge (n);
  Rcpp::LogicalVector is_checkerboard (n);
  Rcpp::LogicalVector is_not_struct_zeros (n);
  Rcpp::LogicalVector can_swap (n);
  Rcpp::LogicalVector did_swap (n);
  Rcpp::NumericVector swap_ps (n);
  //Rcpp::Rcout << "Swapping\n";
  for (int k=0; k<n; k++){
    // select two edges
    i = randint(gen);
    j = randint(gen);
    // make sure i and j are not the same
    if (i == j){
      same_edge[k] = true;
      continue;
    }
    e_from_i = e_from[i];
    e_to_i = e_to[i];
    e_from_j = e_from[j];
    e_to_j = e_to[j];
    // check that it's a checkerboard (diagonal edges don't exist)
    if (edgeSet.find(std::make_pair(e_from_i, e_to_j)) != edgeSet.end() ||
        edgeSet.find(std::make_pair(e_from_j, e_to_i)) != edgeSet.end() ){
      is_checkerboard[k] = false;
      continue;
    }
    // check that there are no structral zeros in the diagonals
    if (zeroSet.find(std::make_pair(e_from_i, e_to_j)) != edgeSet.end() ||
        edgeSet.find(std::make_pair(e_from_j, e_to_i)) != edgeSet.end() ){
      is_not_struct_zeros[k] = false;
      continue;
    }
    can_swap[k] = true;
    // it is a checkerboard and no struct zeros on opposite diag,
    // compute transition prob
    w_pre = w(e_from_i, e_to_i) * w(e_from_j, e_to_j);
    w_post = w(e_from_i, e_to_j) * w(e_from_j, e_to_i);
    swap_p = w_post/w_pre;
    swap_ps[k] = swap_p;
    u = randu(gen);
    if(u < swap_p){
      did_swap[k] = true;
      // make the swap in the hash set
      edgeSet.insert(std::make_pair(e_from_i, e_to_j));
      edgeSet.insert(std::make_pair(e_from_j, e_to_i));
      edgeSet.erase(std::make_pair(e_from_i, e_to_i));
      edgeSet.erase(std::make_pair(e_from_j, e_to_j));
      // make the swap in the edge lists
      e_from[i] = e_from_j;
      e_from[j] = e_from_i;
    }
  }
  //Rcpp::Rcout << "Creating list\n";
  List edges;
  edges["from"] = e_from;
  edges["to"] = e_to;
  edges["same_edge"] = same_edge;
  edges["is_checkerboard"] = is_checkerboard;
  edges["is_not_struct_zeros"] = is_not_struct_zeros;
  edges["can_swap"] = can_swap;
  edges["did_swap"] = did_swap;
  edges["swap_p"] = swap_ps;
  return edges;
}
