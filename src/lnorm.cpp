#include <RcppArmadillo.h>
/* Two Parameter (location, scale) Lognormal Distribution Added by Abiola ADEGBOYEGA
Part of graduate work conducted at the Electrical Engineering Dept, University of Calgary 
Area of Study: Right-tailed distributions in Cloud Traffic */

using namespace Rcpp;
using namespace arma;

arma::vec lnorm_Score(double dY, arma::vec vTheta){

  double dMu=vTheta(0);
  double dSigma2=vTheta(1);
  
  if(dY <= 0){dY = 1;}

  double dMu_s     =(log(dY) - dMu)/dSigma2;
  
  double dSigma2_s =-0.5*(1.0 - pow(log(dY) - dMu,2)/dSigma2 )/dSigma2;

  arma::vec vScore(2);

  vScore(0)=dMu_s;
  vScore(1)=dSigma2_s;

  return vScore;

}
arma::mat lnorm_IM(arma::vec vTheta){

  
  double dSigma2 = vTheta(1);
  arma::mat mIM=zeros(2,2);

  mIM(0,0) = 1.0/dSigma2;
  mIM(1,1) = 1.0/(2.0*pow(dSigma2,2.0));
  return mIM;

}

double dLNORM(double dY, double dMu, double dSigma2, bool bLog = false){

  if(dY <= 0){dY = 1;}
  double dLPDF = -0.5 * log(2.0 * M_PI * dSigma2) - log(dY) - 0.5*pow(log(dY) - dMu,2.0)/dSigma2;
  if(!bLog){
   dLPDF = exp(dLPDF);
  }

  return dLPDF;

}


//
double pLNORM(double dY, double dMu, double dSigma2, bool bLog = false, bool lTail = false) {

  double dP = Rf_plnorm(dY, dMu, pow(dSigma2,0.5),1,0);

  return dP;

}
double qLNORM(double dP, double dMu, double dSigma2, bool bLog = false, bool lTail = false){

  double dQ = Rf_qlnorm(dP, dMu, pow(dSigma2,0.5), 1, 0);

  return dQ;

}
double rLNORM(double dMu, double dSigma2){
  double dY = Rf_rlnorm(dMu, pow(dSigma2,0.5));

  return dY;
}




arma::vec mLNORM(double dMu, double dSigma2){
  arma::vec vMoments(4);
  vMoments(0) = exp(dMu + (dSigma2/2.0));
  vMoments(1) = exp(2*dMu + 2*dSigma2)*(exp(dSigma2) - 1);
  vMoments(2) = 0.0;
  vMoments(3) = 0.0;
  return vMoments;
}
