#ifndef WEIBULL_H
#define WEIBULL_H

double dWEIBULL(double dY, double dAlpha, double dBeta, bool bLog=false);
double pWEIBULL(double dY, double dAlpha, double dBeta);
double qWEIBULL(double dP, double dAlpha, double dBeta);
double rWEIBULL(double dAlpha, double dBeta);
arma::vec mWEIBULL(double dAlpha, double dBeta);
arma::vec weibull_Score(double dY, arma::vec vTheta);
arma::mat weibull_IM(arma::vec vTheta);

#endif
