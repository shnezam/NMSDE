#include <RcppArmadillo.h>
#include "common_fun.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]

using namespace Rcpp;
using namespace arma;
using namespace std;

// in a cube, permute {1,3,2} (permute slices by columns) OR permute {2,3,1}

// [[Rcpp::export]]
arma::cube permute(arma::cube a,int order){
  arma::uword D1=a.n_rows;
  arma::uword D2=a.n_cols;
  arma::uword D3=a.n_slices;
  if(order==132){
    arma::cube output(D1,D3,D2);
    for(arma::uword s = 0; s < D3; s++){
      for(arma::uword c = 0; c < D2; c++){
        for(arma::uword r = 0; r < D1; r++){
          output.at(r, s, c) = a.at(r, c, s);
        }
      }
    }  
    return output;
  }
  else{
    arma::cube output(D2,D3,D1);
    for(arma::uword s = 0; s < D3; s++){
      for(arma::uword c = 0; c < D2; c++){
        for(arma::uword r = 0; r < D1; r++){
          output.at(c,s,r ) = a.at(r, c, s);
        }
      }
    } 
    return output;
  }
}

// [[Rcpp::export]]
arma::mat inChol2spec_c(arma::mat inchol,List init){
  arma::uword P2=inchol.n_rows;
  arma::uword P=sqrt(P2);
  arma::uword Nreal = P * (P + 1) / 2;
  arma::uword K=inchol.n_cols;
  arma::mat indexer=init["indexer"];
  
  arma::mat spec=zeros(P2,K);
  for (arma::uword k =0; k <K ; k++ ){
    
    // inv cholesky at freq k
    arma::cx_mat cholMat= cx_mat(zeros(P,P),zeros(P,P));
    for (arma::uword p=0; p <Nreal ; p++ ){
      cholMat(indexer(p,2)-1,indexer(p,3)-1)=cx_double(inchol(p,k),0);
    }
    for(arma::uword p=Nreal; p<P2; p++){
      cholMat(indexer(p,2)-1,indexer(p,3)-1)=cholMat(indexer(p,2)-1,indexer(p,3)-1)+cx_double(0,inchol(p,k));
    }
    
    //spec matrix at freq k
    arma::cx_mat F_k = inv(cholMat*cholMat.t());
    arma::mat F_k_re = arma::real(F_k);
    arma::mat F_k_im = arma::imag(F_k);
    for (arma::uword p=0; p<Nreal ; p++ ){
      spec(p,k) = F_k_re(indexer(p,2)-1,indexer(p,3)-1);
    }
    for(arma::uword p=Nreal; p<P2; p++){
      spec(p,k) = F_k_im(indexer(p,2)-1,indexer(p,3)-1);
    }
  }
  
  return spec;
}



