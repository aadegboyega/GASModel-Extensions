#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

double dWEIBULL(double dY, double dAlpha, double dBeta, bool bLog=false) {
  //if(dY <= 0){dY = 1.0;}
  //double dLPDF = (dAlpha - 1.0)*log(dY) + (dBeta - 1.0)*log(1.0 - dY) + Rf_lgammafn(dAlpha + dBeta) - Rf_lgammafn(dAlpha) - Rf_lgammafn(dBeta);//
  //double dLPDF = log (dY) + (1/dBeta) - (1/dAlpha)* (pow(dY,dAlpha)* log (dY));
  double dLPDF = log (dAlpha) -  (dAlpha * log (dBeta)) - (pow(dY/dBeta,dAlpha)) + ((dAlpha - 1.0) * log (dY));
  if(!bLog) dLPDF=exp(dLPDF);

  return dLPDF;

}
double pWEIBULL(double dY, double dAlpha, double dBeta) {

  double dP = Rf_pweibull(dY, dAlpha, dBeta, 1, 0);

  return dP;

}
double qWEIBULL(double dP, double dAlpha, double dBeta) {

  double dQ = Rf_qweibull(dP,dAlpha, dBeta,1,0);

  return dQ;

}
//
double rWEIBULL(double dAlpha, double dBeta){
  double dY = Rf_rweibull(dAlpha, dBeta);

  return dY;
}
//
arma::vec mWEIBULL(double dAlpha, double dBeta){
  arma::vec vMoments(4);
  
  //double dMu = mean(dY);
  double dSigmaa = (1.0 + (2.0/dAlpha));
  double dSigmab = (1.0 + (1.0/dAlpha));
  vMoments(0) = dBeta * Rf_lgammafn(dSigmab);
  vMoments(1) = pow(dBeta,2.0) * (Rf_lgammafn(dSigmaa) - pow(Rf_lgammafn(dSigmab),2.0));
  //vMoments(1) = pow(dBeta,2.0) * (Rf_lgammafn(dSigmaa) - pow(dMu,2.0));
  vMoments(2) = 0;
  vMoments(3) = 0;
  return vMoments;
}
//
arma::vec weibull_Score(double dY, arma::vec vTheta){

  double dAlpha = vTheta(0);
  double dBeta  = vTheta(1);

  arma::vec vScore(2);
  //if(dY <= 0){dY = 1.0;}
  double dAlpha_s = dBeta * (pow(dY, dAlpha)* log(dY)) - log(dY) - (1.0/dAlpha);
  double dBeta_s  =  pow(dY, dAlpha) - (1.0/dBeta);
  //double dAlpha_s = (Rf_digamma(dAlpha)/dBeta) - log(dY) - (1/dAlpha);
  //double dBeta_s  = dAlpha/dBeta*(pow((dY * dBeta), dAlpha) - 1.0);

  vScore(0) = dAlpha_s;
  vScore(1) = dBeta_s;

  return vScore;

}
arma::mat weibull_IM( arma::vec vTheta){

  double dAlpha = vTheta(0);
  double dBeta  = vTheta(1);

  arma::mat mIM = zeros(2,2);

  //double dTrig =  Rf_trigamma(dAlpha + dBeta);

  mIM(0,0) = 1.822/pow(dAlpha,2.0);
  mIM(1,1) = pow(dAlpha,2.0)/pow(dBeta,2.0);
  mIM(1,0) = 0.4228/dBeta;
  mIM(0,1) = 0.4228/dBeta;

  return mIM;
}
//
//
