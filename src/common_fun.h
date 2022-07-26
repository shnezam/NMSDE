// Define a header file containing function signatures for functions in
// the A/ subdirectory of src/

// Protect signatures using an inclusion guard.
#ifndef common_fun_H
#define common_fun_H

arma::cube permute(arma::cube a,int order);
arma::mat inChol2spec_c(arma::mat inchol,Rcpp::List init);

#endif
