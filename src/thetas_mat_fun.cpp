#include <RcppArmadillo.h>

// Load header files for common functions
#include "common_fun.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]

using namespace Rcpp;
using namespace arma;
using namespace std;

/////////////////////////// tA2G_c_mat ////////////////////////////////
// [[Rcpp::export]]
arma::cube tA2G_c_mat(arma::mat thetas, arma::cube As, Rcpp::List init) {

  // Takes thetas and As and reconstructs G
  arma::uword P = init["P"];
  arma::uword K = init["K"];
  arma::mat B = init["B"];
  arma::uword P2 = pow(P, 2);

  arma::cube hatG_m = zeros(P2, K, As.n_rows);
  arma::cube theta_tA = zeros(B.n_cols, As.n_rows, P2);
  for (arma::uword j = 0; j < P2; j++) {
    arma::mat As_j = As.slice(j);
    theta_tA.slice(j) = thetas * As_j.t();
    hatG_m.row(j) = B * theta_tA.slice(j);
  }

  return hatG_m;
}


/////////////////////////// loglike_c_mat ////////////////////////////////
////' @export
// // [[Rcpp::export]]
// Rcpp::List loglike_c_mat(arma::mat thetas, arma::cube As, arma::vec lambda, Rcpp::List init, Rcpp::List tildeX) {
//   Rcout << "theta1 " << thetas << std::endl;
//   //Whittle negative loglikelihood
//   arma::uword P = init["P"];
//   arma::uword P2 = pow(P, 2);
//   arma::uword Nreal = P * (P + 1) / 2;
//   arma::uword K = init["K"];
//   arma::mat indexer = init["indexer"];
//   //Rcpp::List freqdom = init["freqdom"];
//   arma::mat D = init["D"];
// 
//   arma::cube hatG_m = tA2G_c_mat(thetas, As, init);
//   arma::cx_mat l = cx_mat(zeros(K, As.n_rows), zeros(K, As.n_rows));
//   for (arma::uword i = 0; i < As.n_rows; i++) {
// 
//     //Rcpp::List freqdom_i = freqdom[i];
//     //arma::cx_mat tildeX_i = freqdom_i[0];
//     arma::cx_mat tildeX_i = tildeX[i];
//     for (arma::uword k = 0; k < K; k++) {
//       //get the inverse cholesky matrix
//       arma::cx_mat Cholmat_i = cx_mat(zeros(P, P), zeros(P, P));
//       for (arma::uword p = 0; p < Nreal; p++) {
//         Cholmat_i(indexer(p, 2) - 1, indexer(p, 3) - 1) = cx_double(hatG_m(p, k, i), 0);
//       }
//       for (arma::uword p = Nreal; p < P2; p++) {
//         Cholmat_i(indexer(p, 2) - 1, indexer(p, 3) - 1) = Cholmat_i(indexer(p, 2) - 1, indexer(p, 3) - 1) + cx_double(0, hatG_m(p, k, i));
//       }
//       
//       //get the loglikelihood
//       arma::cx_mat tildeX_k_i = tildeX_i(span::all, span(k));
//       //l(k, i) = -log_det_sympd(Cholmat_i * Cholmat_i.t()) + det(tildeX_k_i.t() * Cholmat_i * Cholmat_i.t() * tildeX_k_i);
//       l(k,i) = -log(prod(arma::eig_sym(Cholmat_i*Cholmat_i.t())))+det(tildeX_k_i.t()*Cholmat_i*Cholmat_i.t()*tildeX_k_i);
//     }
//   }
//   double loglike = accu(arma::real(l));
//   //Rcout << "loglike " << loglike << std::endl;
//   //Penalized loglikelihood
// 
//   if (lambda.n_elem == 1) {
//     arma::vec lambda = repelem(lambda, thetas.n_cols,1);
//   }
//   //Rcout << "lampen " << lambda << std::endl;
//   //Rcout << "diag " << diagmat(lambda)  << std::endl;
//   double lam_pen = sum(diagmat(lambda) * diagvec(thetas.t() * D * thetas));//lambda:scalar
//   
//   double ploglike = loglike +lam_pen;
// 
//   Rcpp::List L = List::create(Named("ploglike") = ploglike, Named("loglike") = loglike);
//   return L;
// }
////////////////////////////////

// [[Rcpp::export]]
Rcpp::List loglike_c_mat(arma::mat thetas,arma::cube As,arma::vec lambda,List init, Rcpp::List tildeX){
  NumericVector l=wrap(lambda);
  l.attr("dim")=R_NilValue;//vec into NumericVector (2dim into 1dim)
  Function h("loglike");
  return h(thetas,As,l,init,tildeX);
}


/////////////////////////// Workingfish_mat ////////////////////////////////
// [[Rcpp::export]]
Rcpp::List Workingfish_mat(arma::mat thetas,arma::cube As,arma::vec lambda,Rcpp::List init, Rcpp::List tildeX) { 
  arma::uword P=init["P"]; arma::uword K=init["K"]; arma::uword Nreal=init["Nreal"];
  arma::vec omegas=init["omegas"]; arma::vec knots=init["knots"];
  arma::mat B=init["B"]; arma::mat D=init["D"]; arma::mat indexer=init["indexer"];
  arma::cube thetastA=init["theta.tA"];
  Rcpp::List drivG=init["drivG"]; 
  //Rcpp::List freqdom=init["freqdom"];
  arma::uword nbasis=B.n_cols;
  arma::uword norder=nbasis-knots.n_elem+2;
  arma::uword P2=pow(P,2);
  arma::uword R=thetas.n_cols;
  if(lambda.n_elem==1){
    arma::vec lambda_R=repelem(lambda,R,1);
    lambda.set_size(R);
    lambda=lambda_R;
  }	
  arma::cube hatG_m=tA2G_c_mat(thetas,As,init);
  arma::cube spec(size(hatG_m));
  for(arma::uword i=0; i<hatG_m.n_slices; i++){
    spec.slice(i)=inChol2spec_c(hatG_m.slice(i),init);
  }
  arma::mat U_theta=zeros(nbasis,R);
  arma::cube F_theta=zeros(nbasis,nbasis,R);
  arma::cube tmpli=zeros(nbasis,nbasis,R);
  arma::mat theta_chng=zeros(nbasis,R);
  arma::mat curr_thetas=zeros(nbasis,R); 
  arma::mat A_chng=zeros(R*P2,As.n_rows);
  arma::cube curr_As=zeros(As.n_rows,R,P2);
  arma::cube As_perm=permute(As,231); //[m*R*P2] cube into [R*P2*m] cube
  As_perm.reshape((As_perm.n_rows)*(As_perm.n_cols),As_perm.n_slices,1);
  arma::mat A_RP2_m=As_perm.slice(0);//[R*P2*m] cube into [RP2*m] matrix
  Rcpp::List dG=drivG[0];
  cx_cube HGGvec_c=drivG[2];
  field<cx_mat>HGG(P2,P2);
  cx_mat blankm2(P,P);
  HGG.fill(blankm2);
  vec DF(R);
  mat temp1_i(P2,K);
  mat temp2_k_i(P2,P2);
  mat temp8_k_i(P2,P2);
  for (uword i =0; i <As.n_rows ; i++ ){
    mat UA_i=zeros(R*P2,1);
    field<cube>FA_k_i(R,P2); //array(R,P2,K,R,P2)
    cube blankc1=zeros(R,P2,K);
    FA_k_i.fill(blankc1);
    cx_mat G_i(P,P);
    cx_mat F_i(P,P);
    temp1_i=zeros(P2,K);
    field<cube>temp4_i(R,1);//array(L,P2,K,R)---field(R,1)(L,P2,K)
    cube blank=zeros(nbasis,P2,K);
    temp4_i.fill(blank);
    field<cube>temp9_i(nbasis,nbasis,R); //array(P2,P2,K,L,L,R)
    cube blankc2=zeros(P2,P2,K);
    temp9_i.fill(blankc2);
    for(uword k=0; k<K; k++){
      //Rcpp::List freqdom_i=freqdom[i];
      arma::cx_mat tildeX_i=tildeX[i];
      //cx_mat tildeX_i=freqdom_i[0];
      cx_mat tildeX_k_i=tildeX_i(span::all,span(k));
      cx_mat CjtildeX_k_i=tildeX_k_i.t();
      for(uword p=0; p<Nreal; p++){
        G_i(indexer(p,2)-1,indexer(p,3)-1)=hatG_m(p,k,i);
        F_i(indexer(p,2)-1,indexer(p,3)-1)=spec(p,k,i);
        F_i(indexer(p,3)-1,indexer(p,2)-1)=spec(p,k,i);
      }
      for(uword p=Nreal; p<P2; p++){
        G_i(indexer(p,2)-1,indexer(p,3)-1)=G_i(indexer(p,2)-1,indexer(p,3)-1)+cx_double(0,hatG_m(p,k,i));
        F_i(indexer(p,2)-1,indexer(p,3)-1)=F_i(indexer(p,2)-1,indexer(p,3)-1)+cx_double(0,spec(p,k,i));
        F_i(indexer(p,3)-1,indexer(p,2)-1)=F_i(indexer(p,3)-1,indexer(p,2)-1)-cx_double(0,spec(p,k,i));
      }
      mat UA_k_i=zeros(R,P2);
      temp2_k_i=zeros(P2,P2);
      temp8_k_i=zeros(P2,P2);
      uword b=max(find(knots<=omegas[k]));
      vec Btheta_k(R);
      for(uword r=0; r<R; r++){
        Btheta_k[r]=sum(B(k,span(b,norder+b-1))*thetas(span(b,norder+b-1),r));
      }
      //calculate temp-s
      for(uword j1=0;j1<P2;j1++){
        cx_mat dG_j=dG[j1];
        //numeric_limits<double>::epsilon()
        if(abs(hatG_m(j1,k,i))<pow(10,-153)){
          hatG_m(j1,k,i)=pow(10,-153);
        }
        temp1_i(j1,k)=as_scalar(real(CjtildeX_k_i*(dG_j*G_i.t()+ G_i*dG_j.t())*tildeX_k_i))+indexer(j1,4)*((-2)/hatG_m(j1,k,i));
        temp2_k_i.diag()[j1]=indexer(j1,4)*(2/(pow(hatG_m(j1,k,i),2)));
        uvec nzero_l=find(B.row(k)!=0); //OR find(B(k,span::all)!=0)
        for(uword l=0; l<nzero_l.size(); l++){
          for(uword r=0;r<R;r++){
            temp4_i(r,0)(nzero_l(l),j1,k)=temp1_i(j1,k)*As(i,r,j1)*B(k,nzero_l(l));
          }
        }	
        for(uword j2=0; j2<P2; j2++){
          HGG(j1,j2)=reshape(HGGvec_c.slice(j2).col(j1),P,P);
          temp8_k_i(j1,j2)=trace(real(HGG(j1,j2)*F_i))+temp2_k_i(j1,j2);
          for(uword r1=0; r1<R; r1++){
            UA_k_i(r1,j1)=temp1_i(j1,k)*Btheta_k(r1);
            for(uword r2=0; r2<R; r2++){
              FA_k_i(r2,j2)(r1,j1,k)=temp8_k_i(j1,j2)*Btheta_k(r1)*Btheta_k(r2);
            }
          }
          uvec nzero_basis=find(B(k,span::all)!=0);
          for(uword l1=0; l1<nzero_basis.size(); l1++){
            for(uword l2=0; l2<nzero_basis.size(); l2++){
              for(uword r=0; r<R; r++){
                temp9_i(nzero_basis(l1),nzero_basis(l2),r)(j1,j2,k)=temp8_k_i(j1,j2)*B(k,nzero_basis(l1))*B(k,nzero_basis(l2))*As(i,r,j1)*As(i,r,j2);
              }
            }
          }	
        }
      }
      UA_i+= UA_k_i.as_col();
    }//k loop is closed
    
    field<mat>FA_sumk(R,P2);//[R*P2*K*R*P2] into [R*P2*R*P2]
    mat blankm3(R,P2);
    FA_sumk.fill(blankm3);
    cube FA_vec12=zeros(R*P2,R,P2);//[R*P2*R*P2] into [RP2*R*P2]
    mat FA_i=zeros(R*P2,R*P2);//[RP2*R*P2] into [RP2*RP2]
    cube temp4_i_sumk=zeros(nbasis,P2,R);//[L*P2*K*R] into [L*P2*R]
    cube temp10_i=zeros(nbasis,nbasis,R);//[P2*P2*K*L*L*R] into [L*L*R]
    for(uword r2=0; r2<R; r2++){
      for(uword j2=0; j2<P2; j2++){
        FA_sumk(r2,j2)=sum(FA_k_i(r2,j2),2);//step1----sum on k
        FA_vec12.slice(j2).col(r2)=FA_sumk(r2,j2).as_col();//step2-----vectorise on dim (1,2)
        for(uword rp=0; rp<R*P2; rp++){
          FA_i(j2*R+r2,rp)=FA_vec12(rp,r2,j2);//step3----permute(FA_vec12,c(2,3,1))<<(array(R,P2,RP2)) and vectorise dim(1,2)<<(matrix(RP2,RP2))
        }
      }
      temp4_i_sumk.slice(r2)=sum(temp4_i(r2,0),2);//[L*P2*K*R] into [L*P2*R]----som on k
      for(uword l1=0; l1<nbasis; l1++){
        for(uword l2=0; l2<nbasis; l2++){
          temp10_i(l1,l2,r2)=accu(temp9_i(l1,l2,r2));//[P2,P2,K,L,L,R] into [L,L,R]
        }
      }	 
    }
    mat temp5_i=sum(temp4_i_sumk,1);
    U_theta+=temp5_i;
    F_theta=F_theta+temp10_i; 
    A_chng.col(i)=inv_sympd(FA_i)*UA_i;//OR A_chng.col(i)=solve(FA_i,UA_i);
    if(i==As.n_rows-1){
      for(uword r=0;r<R;r++){
        tmpli.slice(r)=inv(F_theta.slice(r)+lambda(r)*D);
        theta_chng.col(r)=tmpli.slice(r)*(U_theta.col(r)+lambda(r)*D*thetas.col(r));
        mat tmpli_F=tmpli.slice(r)*F_theta.slice(r);
        DF(r)=trace(tmpli.slice(r)*F_theta.slice(r));
      }
    }
  }//i loop is closed
  //step-halving algorithm
  bool repit=true;
  double counter=0;
  while(repit){
    double tau=pow(0.5,counter);
    for(uword i=0; i<As.n_rows; i++){
      vec curr_As_vec=A_RP2_m.col(i)-tau*A_chng.col(i);
      curr_As.tube(span(i),span::all)=reshape(curr_As_vec,R,P2);//OR curr_As(span(i),span::all,span::all)=reshape(curr_As_vec,R,P2); 
    }
    for(uword r=0; r<R; r++){
      curr_thetas.col(r)=thetas.col(r)-tau*theta_chng.col(r);//curr_thetas(span::all,span(r),span::all)=reshape(curr_thetas_vec,nbasis,P2);
    }
    double ploglike_curr=loglike_c_mat(curr_thetas,curr_As,lambda,init,tildeX)["ploglike"];
    double ploglike_old=loglike_c_mat(thetas,As,lambda,init,tildeX)["ploglike"];
    if(ploglike_curr<=ploglike_old){
      repit=false;
    }
    else{
      counter+=1;
    }
  }
  Rcpp::List L=List::create(Named("thetas")=curr_thetas, Named("As")=curr_As, _["DF"]=real(DF), _["stepH"]=counter);
  return L;
}


