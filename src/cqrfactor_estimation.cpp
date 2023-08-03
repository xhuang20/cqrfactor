// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <RcppEigen.h>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/Eigenvalues>
#include <Eigen/SVD>
#include <cmath>
#include <algorithm>         // for the sort function in Quantile function.
#include <vector>            // for the quantile function.
#include <stdio.h>           // for the quantile function
#include <string>            // for the quantile funciton
#include "Quantile.h"        // include the Quantile function in a separate file.

// [[Rcpp::depends(RcppEigen)]]
using Eigen::Map;                       
using Eigen::MatrixXd;                  
using Eigen::VectorXd;

//' @title Function to compute the composite quantile factors and loadings.
//' 
//' @description This function uses the MM algorithm to compute the composite 
//'  quantile factors and factor loadings for a given number of factors.
//' 
//' @param y A t by n matrix with n variables and t observations each.
//' @param tau A vector of quantile positions, and \code{tau = (1:q)/(q+1)}, 
//'  where q is the numbrer of quantiles.
//' @param q The number of quantiles.
//' @param r The number of factors.
//' @param tol Tolerance for convergence of the MM algorithm. Defaults to 1e-6
//'  when the minconver option is 1. When the minconver option is 0, the 
//'  algorithm's convergence depends on the average parameter value change between
//'  iteration steps, and it can take long time to converge. We recommend to use 
//'  a bigger number such as 1e-3 when minconver = 1.
//' @param maxit The maximum iteration number between factor estimation and factor 
//'  loading estimation. Defaults to 500.
//' @param maxit_factor The maximum iteration number in the MM algorithm in factor
//'  estimation. Defaults to 500.
//' @param maxit_loading The maximum iteration number in the MM algorithm in factor
//'  loading estimation. Defaults to 500.
//' @param convergence If the value is 0, the algorithm converges when the minimum
//'  change among all estimated parameters between two iteration steps is less
//'  than or equal to tol. If the value is 1, the algorithm converges when the 
//'  average change of all estimated parameters is less than or equal to tol.
//'  This can take a long computation time. If minconver = 1, set tol to a bigger
//'  number such as 1e-3. Defaults to 0. 
//' @param seed A non-negative integer for random seed used in the C++ function \code{srand} to
//'  initialize factor values. If seed is -1, the initial value of the
//'  factors will be set to the eigenvectors of the outer product of the data
//'  matrix y. Defaults to 1.
//' 
//' @return \code{mmfactor} returns a list with the following components:
//'  \item{fmat}{The estimated factors, a t by r matrix, where r is the number
//'  of factors.} \item{lmat}{The estimated factor loading matrix, an n by r 
//'  matrix.} \item{quantiles}{The q estimated quantiles of the error term.} 
//'  \item{iteration_mat}{A \code{maxit} by 3 matrix. Each row corresponds to an
//'  iteration between factor estimation and factor loading estimation. The 
//'  second element is the number of iterations the MM uses for estimating the 
//'  factors while the third element is the number of iteration for the factor
//'  loadings.}
//' 
//' @examples
//' \dontrun{
//' 
//' }
//' 
// [[Rcpp::export]]
Rcpp::List mmfactor(const Eigen::MatrixXd & y,       
                    const Eigen::VectorXd & tau,     
                    int q = 5,                       
                    int r = 3,                       
                    double tol = 1e-6,               
                    int maxit = 500,
                    int maxit_factor = 500,
                    int maxit_loading = 500,
                    int convergence = 0,
                    int seed = 1) {
  
  int i, j, s;
  int t = y.rows();   
  int n = y.cols();   
  
  // Get the initial estimates for factors and factor loadings.
  MatrixXd outpy = y * y.transpose();
  MatrixXd innpy = y.transpose() * y;
  MatrixXd fmat0,lmat0;
  if (seed == -1) {
    Eigen::SelfAdjointEigenSolver<MatrixXd> es(outpy);
    fmat0 = es.eigenvectors().rowwise().reverse().leftCols(r);
    lmat0 = ((fmat0.transpose() * fmat0).llt().solve(fmat0.transpose() * y)).transpose();
  } else {
    srand(seed);
    fmat0 = MatrixXd::Random(t,r);
    lmat0 = ((fmat0.transpose() * fmat0).llt().solve(fmat0.transpose() * y)).transpose();
  }
  
  MatrixXd rmat = y - fmat0 * lmat0.transpose();
  
  Rcpp::NumericVector tau_vec(tau.data(),  tau.data()  + tau.size());
  Rcpp::NumericVector res_vec(rmat.data(), rmat.data() + rmat.size());
  Rcpp::NumericVector quantiles = Quantile(res_vec, tau_vec);
  Rcpp::NumericVector quantiles0;
  
  MatrixXd fmat(t,r), lmat(n,r);
  MatrixXd rmat0(t,n);
  
  int iter = 0, iter_loading= 0, iter_factor = 0;
  double delta_fmat = 10000.0, delta_lmat = 10000.0, delta = 10000.0;
  double tn, e0, epsilon;
  tn = tol / (t * n);
  e0 = -tn / log(tn);
  epsilon = (e0 - tn) / (1 + log(e0));
  
  MatrixXd A(q*t,n);
  MatrixXd c  = (2 * tau).array() - 1;
  MatrixXd c1 = MatrixXd::Zero(r,q);
  MatrixXd c2 = MatrixXd::Zero(t,q);    
  MatrixXd c3 = MatrixXd::Zero(r,q);    
  MatrixXd c4 = MatrixXd::Zero(t,q);    
  MatrixXd c5 = MatrixXd::Zero(1,r);
  MatrixXd c6 = MatrixXd::Zero(r,r);    
  MatrixXd tc = MatrixXd::Zero(r,n);    
  MatrixXd tc2 = MatrixXd::Zero(r,n);   
  
  MatrixXd cdenor = MatrixXd::Zero(r,r);
  MatrixXd cnumer = MatrixXd::Zero(r,1);
  
  MatrixXd d1 = MatrixXd::Zero(r,q);
  MatrixXd d2 = MatrixXd::Zero(n,q);    
  MatrixXd d3 = MatrixXd::Zero(r,q);    
  MatrixXd d4 = MatrixXd::Zero(n,q);    
  MatrixXd d5 = MatrixXd::Zero(1,r);
  MatrixXd d6 = MatrixXd::Zero(r,r);    
  MatrixXd td = MatrixXd::Zero(r,t);    
  MatrixXd td2 = MatrixXd::Zero(r,t);  
  
  MatrixXd ddenor = MatrixXd::Zero(r,r);
  MatrixXd dnumer = MatrixXd::Zero(r,1);
  MatrixXd imat   = MatrixXd::Identity(r,r);
  MatrixXd itermat = MatrixXd::Zero(maxit,3);
  
  // Iterate between factor and factor loading estimation.
  while (iter < maxit && delta > tol) {
    iter++;
    
    // Given the factors, update the factor loadings and error quantiles.
    d5 = fmat0.colwise().sum();    
    iter_loading = 0;       
    delta_lmat   = 10000.0; 
    while (delta_lmat > tol && iter_loading < maxit_loading) {
      iter_loading++;
      
      quantiles0 = clone(quantiles);
      
      for (s = 0; s < q; s++) {
        rmat0 = (y - fmat0 * lmat0.transpose()).array() - quantiles[s]; 
        A.block(s*t,0,t,n) = (rmat0.cwiseAbs().array() + epsilon).matrix().cwiseInverse();
        d2.col(s) = y.cwiseProduct(A.block(s*t,0,t,n)).colwise().sum().transpose();
        d4.col(s) = A.block(s*t,0,t,n).colwise().sum().transpose(); 
      }
      
      for (i = 0; i < n; i++) {
        
        ddenor.setZero();
        dnumer.setZero();
        
        for (s = 0; s < q; s++) {
          td        = A.block(s*t,0,t,n).col(i).transpose().replicate(r,1).cwiseProduct(y.col(i).transpose().replicate(r,1));  
          d1.col(s) = fmat0.transpose().cwiseProduct(td).rowwise().sum();     
          td2       = fmat0.transpose().cwiseProduct(A.block(s*t,0,t,n).col(i).transpose().replicate(r,1).cwiseSqrt());       
          d6        = td2 * td2.transpose();           
          d3.col(s) = fmat0.transpose().cwiseProduct(A.block(s*t,0,t,n).col(i).transpose().replicate(r,1)).rowwise().sum();
          
          ddenor += d6 - ((d3.col(s) * d3.col(s).transpose()).array() / d4.coeff(i,s)).matrix(); 
          dnumer += d1.col(s) + c(s) * d5.transpose() - d2.col(s).coeff(i) * d3.col(s) / d4.col(s).coeff(i) -
            c(s) * t * d3.col(s) / d4.col(s).coeff(i);                                          
        }
        
        lmat.row(i) = (ddenor.completeOrthogonalDecomposition().solve(dnumer)).transpose();  
      }
      
      for (s = 0; s < q; s++) {
        quantiles[s] = ((y - fmat0 * lmat.transpose()).cwiseProduct(A.block(s*t,0,t,n)).sum() + c(s) * n * t) / A.block(s*t,0,t,n).sum();
      }
      
      if (convergence == 0) {
        delta_lmat = Rcpp::min(Rcpp::abs(quantiles - quantiles0)) +
          (lmat - lmat0).cwiseAbs().minCoeff();
      } else if (convergence == 1) {
        delta_lmat = Rcpp::mean(Rcpp::abs(quantiles - quantiles0)) +
          (lmat - lmat0).cwiseAbs().mean();
      }
      
      lmat0 = lmat;
    }
    
    // Given the factor loadings, update the factors and error quantiles.
    c5 = lmat0.colwise().sum();    
    iter_factor = 0;      
    delta_fmat = 10000.0; 
    while (delta_fmat > tol && iter_factor < maxit_factor) {
      iter_factor++;
      quantiles0 = clone(quantiles);
      
      for (s = 0; s < q; s++) {
        rmat0 = (y - fmat0 * lmat0.transpose()).array() - quantiles[s]; 
        A.block(s*t,0,t,n) = (rmat0.cwiseAbs().array() + epsilon).matrix().cwiseInverse();
        c2.col(s) = y.cwiseProduct(A.block(s*t,0,t,n)).rowwise().sum();
        c4.col(s) = A.block(s*t,0,t,n).rowwise().sum(); 
      }
      
      for (j = 0; j < t; j++) {
        cdenor.setZero();
        cnumer.setZero();
        for (s = 0; s < q; s++) {
          tc        = A.block(s*t,0,t,n).row(j).replicate(r,1).cwiseProduct(y.row(j).replicate(r,1));  
          c1.col(s) = lmat0.transpose().cwiseProduct(tc).rowwise().sum();                   
          tc2       = lmat0.transpose().cwiseProduct(A.block(s*t,0,t,n).row(j).replicate(r,1).cwiseSqrt());   
          c6        = tc2 * tc2.transpose();   
          c3.col(s) = lmat0.transpose().cwiseProduct(A.block(s*t,0,t,n).row(j).replicate(r,1)).rowwise().sum();
          
          cdenor += c6 - ((c3.col(s) * c3.col(s).transpose()).array() / c4.coeff(j,s)).matrix(); 
          cnumer += c1.col(s) + c(s) * c5.transpose() - c2.col(s).coeff(j) * c3.col(s) / c4.col(s).coeff(j) -
            c(s) * n * c3.col(s) / c4.col(s).coeff(j);                                         
        }
        
        fmat.row(j) = (cdenor.completeOrthogonalDecomposition().solve(cnumer)).transpose();
      }
      
      for (s = 0; s < q; s++) {
        quantiles[s] = ((y - fmat * lmat0.transpose()).cwiseProduct(A.block(s*t,0,t,n)).sum() + c(s) * n * t) / A.block(s*t,0,t,n).sum();
      }
      
      if (convergence == 0) {
        delta_fmat = Rcpp::min(Rcpp::abs(quantiles - quantiles0)) +
          (fmat - fmat0).cwiseAbs().minCoeff();
      } else if (convergence == 1) {
        delta_fmat = Rcpp::mean(Rcpp::abs(quantiles - quantiles0)) +
          (fmat - fmat0).cwiseAbs().mean();
      }
      
      fmat0 = fmat;
    }
 
    delta = std::max(delta_fmat,delta_lmat);
    itermat(iter-1, 0) = iter;
    itermat(iter-1, 1) = iter_factor;
    itermat(iter-1, 2) = iter_loading;
  }   
  
  MatrixXd innpfmat = fmat.transpose() * fmat / t;
  MatrixXd innplmat = lmat.transpose() * lmat / n;
  
  Eigen::SelfAdjointEigenSolver<MatrixXd> es_innpfmat(innpfmat);
  MatrixXd innpfmat_sqrt = es_innpfmat.operatorSqrt();
  MatrixXd svdmat = innpfmat_sqrt * innplmat * innpfmat_sqrt;
  Eigen::JacobiSVD<MatrixXd> svd(svdmat, Eigen::ComputeFullU | Eigen::ComputeFullV);
  MatrixXd rotamat = innpfmat_sqrt.inverse() * svd.matrixU();
  fmat = fmat * rotamat;
  lmat = lmat * rotamat.inverse().transpose();
  
  return Rcpp::List::create(
    Rcpp::Named("fmat") = fmat,
    Rcpp::Named("lmat") = lmat,
    Rcpp::Named("quantiles") = quantiles,
    Rcpp::Named("iteration_mat") = itermat);
}
