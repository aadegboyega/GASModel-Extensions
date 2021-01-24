#ifndef GUMBEL_H
#define GUMBEL_H

double dGUMBEL(double dY, double dAlpha, double dBeta, bool bLog=false);
double pGUMBEL(double dY, double dAlpha, double dBeta);
double qGUMBEL(double dP, double dAlpha, double dBeta);
double rGUMBEL(double dAlpha, double dBeta);
arma::vec mGUMBEL(double dAlpha, double dBeta);
arma::vec gumbel_Score(double dY, arma::vec vTheta);
arma::vec gumbel_Scaler(double dY, arma::vec vTheta);
arma::mat gumbel_IM(arma::vec vTheta2);

#endif
