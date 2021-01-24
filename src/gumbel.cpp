#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


double dGUMBEL(double dY, double dAlpha, double dBeta, bool bLog=false) {

  double dLPDF = - log(dBeta) - ((dY - dAlpha)/dBeta) - exp(-(dY - dAlpha)/dBeta);

  if(!bLog) dLPDF=exp(dLPDF);

  return dLPDF;

}

double pGUMBEL(double dY, double dAlpha, double dBeta) {

  double dP = Rf_pgamma(dY, dAlpha, dBeta, 1, 0);

  return dP;

}
double qGUMBEL(double dP, double dAlpha, double dBeta){

  double dQ = Rf_qgamma(dP, dAlpha, dBeta, 1 ,0 );

  return dQ;

}
double rGUMBEL(double dAlpha, double dBeta){
  double dY = Rf_rgamma(dAlpha, dBeta);

  return dY;
}

arma::vec mGUMBEL(double dAlpha, double dBeta){
  arma::vec vMoments(4);
  vMoments(0) = dAlpha + (0.5772 * dBeta);
  vMoments(1) = 0.1666* pow(3.142*dBeta,2.0);
  vMoments(2) = 0;
  vMoments(3) = 0;
  return vMoments;
}

arma::vec gumbel_Score(double dY, arma::vec vTheta){

  double dAlpha = vTheta(0);
  double dBeta  = vTheta(1);

  arma::vec vScore(2);

  double z = (dY - dAlpha)/dBeta;
  double dAlpha_s = 1.0/dBeta - ((1.0/dBeta)* exp(-z));
  double dBeta_s  = 1.0/dBeta + (dY - dAlpha)/pow(dBeta,2.0) - ((dY - dAlpha)/pow(dBeta,2.0)* exp(-z));

  vScore(0) = dAlpha_s;
  vScore(1) = dBeta_s;

  return vScore;

}
arma::vec gumbel_Scaler(double dY, arma::vec vTheta){
 
  double dAlpha = vTheta(0);
  double dBeta  = vTheta(1);
  double zed = (dY - dAlpha)/dBeta;

  //arma::mat gumbel_X(arma::vec vTheta);
  //arma::mat mIM(2,2);
  //mIM(0,0) = (-1.0/pow(dBeta, 2.0)) * exp(-zed);
  //mIM(0,1) = (1.0/pow(dBeta, 2.0)) + (1.0/pow(dBeta, 2.0))* exp(-zed) - ((1.0/pow(dBeta, 3.0))* (dY - dAlpha) * exp(-zed));
  //mIM(1,0) = (1.0/pow(dBeta, 2.0)) + (1.0/pow(dBeta, 2.0))* exp(-zed) - ((1.0/pow(dBeta, 3.0))* (dY - dAlpha) * exp(-zed));
  //mIM(1,1) = (1.0/pow(dBeta, 2.0)) - ((2.0/pow(dBeta, 3.0))* (dY - dAlpha)) + ((2.0/pow(dBeta, 3.0))* (dY - dAlpha) * exp(-zed)) - ((1.0/pow(dBeta, 4.0))* pow((dY - dAlpha),2.0) * exp(-zed));
  //return mIM;
  
  arma::mat vtheta2(2,2);
  vtheta2(0,0) = (-1.0/pow(dBeta, 2.0)) * exp(-zed);
  vtheta2(0,1) = (-1.0/pow(dBeta, 2.0)) + (1.0/pow(dBeta, 2.0))* exp(-zed) - ((1.0/pow(dBeta, 3.0))* (dY - dAlpha) * exp(-zed));
  vtheta2(1,0) = (-1.0/pow(dBeta, 2.0)) + (1.0/pow(dBeta, 2.0))* exp(-zed) - ((1.0/pow(dBeta, 3.0))* (dY - dAlpha) * exp(-zed));
  vtheta2(1,1) = (1.0/pow(dBeta, 2.0)) - ((2.0/pow(dBeta, 3.0))* (dY - dAlpha)) + ((2.0/pow(dBeta, 3.0))* (dY - dAlpha) * exp(-zed)) - ((1.0/pow(dBeta, 4.0))* pow((dY - dAlpha),2.0) * exp(-zed));
  return vtheta2;

}

arma::mat gumbel_IM(arma::vec vTheta2){

  //double dAlpha = vTheta(0);
  //double dBeta  = vTheta(1);

  arma::mat mIM(2,2);

  mIM(0,0) = vTheta2(0,0);
  mIM(0,1) = vTheta2(0,1);
  mIM(1,0) = vTheta2(1,0);
  mIM(1,1) = vTheta2(1,1);

  return mIM;
}


