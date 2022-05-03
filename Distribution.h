
#pragma once
#include<vector>
#include<string>
#include "Vector.h"
#include "Matrix.h"
#include "BTC.h"
#include "interface.h"

using namespace std;

enum class distribution_type {nonparameteric, normal, lognormal, gamma, levy, exponential, inverse_gaussian};

class CDistribution: public Interface
{
public:
	CDistribution(void);
	CDistribution(string name);
	~CDistribution(void);
    static int NumberOfCoreParameters(string _name);
    bool CreateDistribution(const map<string,string> &Arguments);
	vector<double> params;
	string name;
	double evaluate(double x);
	double evaluate_CDF(double x, bool flux_w = false);
    static bool HasCommand(const string &cmd);
    vector<string> commands();
    static vector<string> Commands();
    static vector<string> list_of_commands;
    double pi;
	int n;
	vector<double> s;
	vector<double> e;
    CTimeSeries<double> inverse_cumulative;
    CTimeSeries<double> density;
    distribution_type DistributionType;
    CDistribution(int nn);
	CDistribution(const CDistribution &C);
    bool SetInverseCumulative(const map<string,string> &arguments);
    bool WriteInverseCumulativeToFile(const map<string,string> &arguments);
    double InverseCumulativeValue(double u);
    double CumulativeValue(double x);
	CDistribution operator = (const CDistribution &C);
	int GetRand();
	double inverseCDF(double u, bool flux_weight=false);
	double map_normal_to(double z);
    bool readfromfile(const string &filename);
    bool WriteToFile(const map<string,string> Arguments);
    vector<double> SetRangeBasedOnMeanStd(const double &stdcoeff=3);
    bool Execute(const string &cmd, const map<string,string> &arguments);
};

//double erf(double x);
//double erfc(double x);
double Gammapdf(double x, double k, double theta);
double gamma(double x);
double unitrandom();
double getstdnormalrand();
double getnormalrand(double mu, double std);
CVector getnormal(CVector &mu, CMatrix &sigma);
CMatrix getnormal(int m, int n, double mu, double std);
CVector getnormal(int m, double mu, double std);
double getlognormalrand(double mu, double std);
CVector getlognormal(CVector &mu, CMatrix &sigma);
CMatrix getlognormal(int m, int n, double mu, double std);
CVector getlognormal(int m, double mu, double std);
double stdnormal_cdf(double u);
double GetRndUniF(double xmin, double xmax);
double std_normal_phi_inv(double u);
