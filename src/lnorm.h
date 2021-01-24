#ifndef LNORM_H
#define LNORM_H

arma::vec lnorm_Score(double dY, arma::vec vTheta);
arma::mat lnorm_IM(arma::vec vTheta);
double dLNORM(double dY, double dMu, double dSigma2, bool bLog = false);
double pLNORM(double dY, double dMu, double dSigma2, bool bLog = false, bool lTail = false);
double qLNORM(double dP, double dMu, double dSigma2, bool bLog = false, bool lTail = false);
double rLNORM(double dMu, double dSigma2);
arma::vec mLNORM(double dMu, double dSigma2);



#endif
