#pragma once
#include "Vector.h"
#include "Matrix.h"
#include "NormalDist.h"
#include <random>

class CNormalDist
{
public:
	CNormalDist(void);
	~CNormalDist(void);
	double unitrandom();
	double getstdnormalrand();
	double getnormalrand(double mu, double std);
private:
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen; //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis;

};
CVector getnormal(CVector &mu, CMatrix &sigma);
CMatrix getnormal(int m, int n, double mu, double std);
CVector getnormal(int m, double mu, double std);
double getlognormalrand(double mu, double std);
CVector getlognormal(CVector &mu, CMatrix &sigma);
CMatrix getlognormal(int m, int n, double mu, double std);
CVector getlognormal(int m, double mu, double std);
double getnormalrand(double mu, double std);
double getnormalcdf(double x, double mu = 0, double std = 1);
double getstdnormalrand();
double stdnormal_inv(double p);
double unitrandom();

