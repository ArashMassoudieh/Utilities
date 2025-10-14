#ifndef NORMALDIST_H
#define NORMALDIST_H

#include "Vector.h"
#include "Matrix.h"

class CNormalDist
{
public:
    CNormalDist(void);
    ~CNormalDist(void);

    // Static methods for random number generation
    static double unitrandom();
    static double getstdnormalrand();
    static double getnormalrand(double mu, double std);
    static double getlognormalrand(double mu, double std);

    // Likelihood calculation
    static double likelihood_mixed(double x_mod, double x_obs, double std_ln, double std_n);

    // Free functions for vector/matrix operations
    CVector getnormal(CVector &mu, CMatrix &sigma);
    CMatrix getnormal(int m, int n, double mu, double std);
    CVector getnormal(int m, double mu, double std);

    CVector getlognormal(CVector &mu, CMatrix &sigma);
    CMatrix getlognormal(int m, int n, double mu, double std);
    CVector getlognormal(int m, double mu, double std);

    // CDF and inverse CDF functions
    static double stdnormal_cdf(double u);
    static double getnormalcdf(double x, double mu = 0, double std = 1);
    static double stdnormal_inv(double p);

// PDF functions (optional - enable with arma flag)
#ifdef _arma
    static double getpdfnormal(CVector &X, CVector &mu, CMatrix &std);
    static double getpdflognormal(CVector &X, CVector &mu, CMatrix &std);
#endif


};



#endif // NORMALDIST_H
