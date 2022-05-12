#include "Distribution.h"
#include "gsl/gsl_cdf.h"
#include "Utilities.h"
#include "interface.h"

vector<string> CDistribution::list_of_commands = vector<string>({"CreateDistribution","WriteToFile", "SetInverseCumulative", "WriteInverseCumulativeToFile"});

CDistribution::CDistribution(void):Interface()
{
	pi = 4 * atan(1.0);
}

bool CDistribution::CreateDistribution(const map<string,string> &Arguments)
{
    if (Arguments.count("type")==0)
        return false;
    name = Arguments.at("type");
    params.clear();
    for (map<string,string>::const_iterator it=Arguments.begin(); it!=Arguments.end(); it++)
    {
        if (aquiutils::left(it->first,1) == "p" && aquiutils::isnumber(aquiutils::right(it->first,it->first.size()-1)))
            params.push_back(atof(it->second.c_str()));
    }
    return true;
}

CDistribution::~CDistribution(void)
{
}

CDistribution::CDistribution(string _name)
{
	pi = 4 * atan(1.0);
	name = _name;
    if (name == "normal") {DistributionType=distribution_type::normal;}
    if (name == "lognormal") {DistributionType=distribution_type::lognormal;}
    if (name == "levy") {DistributionType=distribution_type::levy;}
    if (name == "exp") {DistributionType=distribution_type::exponential;}
    if (name == "invgaussian") {DistributionType=distribution_type::inverse_gaussian;}
    if (name == "gamma") {DistributionType=distribution_type::gamma;}
    if (name == "nonparametric") {DistributionType=distribution_type::nonparameteric;}
    params.resize(NumberOfCoreParameters(name)+1);

}

int CDistribution::NumberOfCoreParameters(string _name)
{
    if (_name == "normal") {return 2;}
    if (_name == "lognormal") {return 2;}
    if (_name == "levy") {return 1;}
    if (_name == "exp") {return 1;}
    if (_name == "invgaussian") {return 2;}
    if (_name == "gamma") {return 2;}
    if (_name == "nonparametric") {return 0;}
    return -1;
}


double CDistribution::evaluate(double x)
{
    double out = 0;
    if (name == "nonparameteric")
        return density.interpol(x);
    int m = NumberOfCoreParameters(name)+1;
    int n=params.size()/m;
    for (unsigned int i=0; i<n; i++)
    {
        if (name == "normal")
            out+= params[i*m] / (sqrt(2*pi)*params[i*m+2])*exp(-pow(x - params[i*m+1], 2) / (2 * params[i*m+2] * params[i*m+2]));
        if (name == "lognormal")
            out+= params[i*m] / (sqrt(2*pi)*params[i*m+2] * x)*exp(-pow(log(x) - params[i*m+1], 2) / (2 * params[i*m+2] * params[i*m+2]));
        if (name == "levy")
            out += params[i*m]*sqrt(params[i*m+1] / (2 * pi))*exp(-params[i*m+1] / (2 * x)) / pow(x, 1.5);
        if (name == "exp")
            out +=  params[i*m] / params[i*m+1] * exp(x / params[i*m+1]);
        if (name == "invgaussian")
            out +=  params[i*m]*sqrt(params[i*m+2] / (2 * pi*pow(x, 3)))*exp(-params[i*m+2] * pow(x - params[i*m+1], 2) / (2 * params[i*m+1] * params[i*m+1] * x));
        if (name == "gamma")
            out+= params[i*m]*Gammapdf(x,params[i*m+1],params[i*m+2]);

    }
    return out;
}

double Gammapdf(double x, double k, double theta)
{
	return 1 / pow(theta, k)*pow(x, k - 1)*exp(-x / theta) / gamma(k);
}

double gamma(double x)
{
	int i, k, m;
	double ga, gr, r, z;
	double m_PI = atan(1.0)*4.0;
	static double g[] = {
		1.0,
		0.5772156649015329,
		-0.6558780715202538,
		-0.420026350340952e-1,
		0.1665386113822915,
		-0.421977345555443e-1,
		-0.9621971527877e-2,
		0.7218943246663e-2,
		-0.11651675918591e-2,
		-0.2152416741149e-3,
		0.1280502823882e-3,
		-0.201348547807e-4,
		-0.12504934821e-5,
		0.1133027232e-5,
		-0.2056338417e-6,
		0.6116095e-8,
		0.50020075e-8,
		-0.11812746e-8,
		0.1043427e-9,
		0.77823e-11,
		-0.36968e-11,
		0.51e-12,
		-0.206e-13,
		-0.54e-14,
		0.14e-14 };

	if (x > 171.0) return 1e308;    // This value is an overflow flag.
	if (x == (int)x) {
		if (x > 0.0) {
			ga = 1.0;               // use factorial
			for (i = 2; i<x; i++) {
				ga *= i;
			}
		}
		else
			ga = 1e308;
	}
	else {
		if (fabs(x) > 1.0) {
			z = fabs(x);
			m = (int)z;
			r = 1.0;
			for (k = 1; k <= m; k++) {
				r *= (z - k);
			}
			z -= m;
		}
		else
			z = x;
		gr = g[24];
		for (k = 23; k >= 0; k--) {
			gr = gr*z + g[k];
		}
		ga = 1.0 / (gr*z);
		if (fabs(x) > 1.0) {
			ga *= r;
			if (x < 0.0) {
				ga = -m_PI / (x*ga*sin(m_PI*x));
			}
		}
	}
	return ga;
}


CDistribution::CDistribution(int nn)
{
	pi = 4 * atan(1.0);
	n = nn;
	e.resize(n);
	s.resize(n);
}

CDistribution::CDistribution(const CDistribution &C)
{
	pi = 4 * atan(1.0);
	n = C.n;
	e.resize(n);
	s.resize(n);
	for (int i = 0; i<n; i++)
	{
		e[i] = C.e[i];
		s[i] = C.s[i];
	}
    inverse_cumulative = C.inverse_cumulative;
    density = C.density;

}

CDistribution CDistribution::operator = (const CDistribution &C)
{
	pi = 4 * atan(1.0);
	n = C.n;
	e.resize(n);
	s.resize(n);
	for (int i = 0; i<n; i++)
	{
		e[i] = C.e[i];
		s[i] = C.s[i];
	}
    inverse_cumulative = C.inverse_cumulative;
    density = C.density;
	return *this;

}

int CDistribution::GetRand()
{
	double x = GetRndUniF(0, 1);
	int ii = 0;
	for (int i = 0; i<n - 1; i++)
	{
		if (x<e[i] && x>s[i])
			ii = i;
	}
	return ii;

}



double CDistribution::inverseCDF(double u,bool flux_w)
{
	if (aquiutils::tolower(name) == "exp" || aquiutils::tolower(name) == "exponential")
	{
		return -params[0] * log(1.0 - u);
	}
	if (aquiutils::tolower(name) == "gamma")
	{
		return gsl_cdf_gamma_Pinv(u, params[0], params[1]);
	}
	if (aquiutils::tolower(name) == "lognormal" || aquiutils::tolower(name) == "log-normal")
	{
		if (flux_w)
            return gsl_cdf_lognormal_Pinv(u, params[0]+pow(params[1],2), params[1]);
        else
            return gsl_cdf_lognormal_Pinv(u, params[0], params[1]);
	}
	if (aquiutils::tolower(name) == "normal" || aquiutils::tolower(name) == "gaussian")
	{
		return gsl_cdf_gaussian_Pinv(u, params[1])+ params[0];
	}
	if (aquiutils::tolower(name) == "power" || aquiutils::tolower(name) == "power_law")
	{
		return pow(2 * u, 1 / (params[0] - 1))*params[1];
	}
	if (aquiutils::tolower(name) == "nonparameteric")
	{
        return inverse_cumulative.interpol(u);
	}
	return 0;

}

double GetRndUniF(double xmin, double xmax)
{
	double a = double(rand());
	double k = double(RAND_MAX);
	return a / k*(xmax - xmin) + xmin;
}

double std_normal_phi_inv(double u)
{
	double u1 = gsl_cdf_ugaussian_Pinv(u);
	double out = 1.0 / sqrt(2.0 * 4.0 * atan(1.0))*exp(-0.5*pow(u1, 2.0));
	return out;
}

double CDistribution::evaluate_CDF(double x, bool flux_w)
{
    if (aquiutils::tolower(name) == "lognormal" || aquiutils::tolower(name) == "log-normal")
	{
		if (flux_w)
            return gsl_cdf_lognormal_P(x, params[0]+pow(params[1],2), params[1]);
        else
            return gsl_cdf_lognormal_P(x, params[0], params[1]);
	}
}

bool CDistribution::readfromfile(const string &filename)
{
    if (!density.readfile(filename))
    {
        cout<< "File [" << filename << "] was not found!"<<endl;
        return false;
    }
    CTimeSeries<double> cumulative = density.getcummulative();
    cumulative = cumulative/cumulative.GetC(cumulative.n-1);
    inverse_cumulative = cumulative.inverse_cumulative_uniform();
    return true;
}

vector<string> CDistribution::commands()
{
    return Commands();
}

vector<string> CDistribution::Commands()
{
    //return vector<string>();
    return list_of_commands;
}

bool CDistribution::HasCommand(const string &cmd)
{
    if (aquiutils::lookup(Commands(),cmd)!=-1)
        return true;
    else
        return false;
}

bool CDistribution::WriteToFile(const map<string,string> Arguments)
{
    int nbins = 100;
    if (Arguments.count("filename")==0)
        return false;
    if (params.size()==0)
        return false;
    if (Arguments.count("nbins")>0)
        nbins = aquiutils::atof(Arguments.at("nbins"));

    vector<double> range = SetRangeBasedOnMeanStd(3);
    CTimeSeries<double> timeseries;
    for (double x = range[0]; x<=range[1]; x+=(range[1]-range[0])/nbins)
    {
        timeseries.append(x, evaluate(x));
    }
    timeseries.writefile(Arguments.at("filename"));
    return true;
}

vector<double> CDistribution::SetRangeBasedOnMeanStd(const double &stdcoeff)
{
    vector<double> range(2);
    int n=NumberOfCoreParameters(name);
    range[0] = 1e12;
    range[1] = -1e12;
    for (int i=0; i<params.size()/NumberOfCoreParameters(name); i++)
    {   if (name=="normal")
        {
            range[0]=min(params[i*n+1]-stdcoeff*params[i*n+2],range[0]);
            range[1]=max(params[i*n+1]+stdcoeff*params[i*n+2],range[1]);
        }
        else if (name=="lognormal")
        {
            range[0]=min(exp(params[i*n+1])*exp(-stdcoeff*params[i*n+2]),range[0]);
            range[1]=max(exp(params[i*n+1])*exp(stdcoeff*params[i*n+2]),range[1]);
        }
        else if (name=="exponential")
        {
            range[0] = 0;
            range[1] = max(stdcoeff/params[i*n+1],range[1]);

        }
    }
    return range;
}

FunctionOutPut CDistribution::Execute(const string &cmd, const map<string,string> &arguments)
{
    FunctionOutPut output;
    if (cmd=="CreateDistribution")
    {   output.success = CreateDistribution(arguments);
        output.output = this;

    }
    if (cmd=="WriteToFile")
    {   output.success = WriteToFile(arguments);

    }
    if (cmd=="SetInverseCumulative")
    {   output.success = SetInverseCumulative(arguments);

    }
    if (cmd=="WriteInverseCumulativeToFile")
    {   output.success = WriteInverseCumulativeToFile(arguments);

    }
    return output;
}

bool CDistribution::SetInverseCumulative(const map<string,string> &arguments)
{
    if (arguments.count("ninc")==0) return false;
    int ninc = aquiutils::atoi(arguments.at("ninc"));
    inverse_cumulative.clear();
    double epsilon = 1e-8;

    for (int i = 0; i < ninc+1; i++)
    {
        double u = double(i) / double(ninc)*(1 - 2 * epsilon) + epsilon;
        inverse_cumulative.append(u, InverseCumulativeValue(u));
    }
    return true;
}

bool CDistribution::WriteInverseCumulativeToFile(const map<string,string> &arguments)
{
    if (inverse_cumulative.n==0)
        return false;
    if (arguments.count("filename")==0) return false;
    return inverse_cumulative.writefile(arguments.at("filename"));
}

double CDistribution::InverseCumulativeValue(double u)
{
    double x2 = 1000;
    double x1 = 0.001;
    double tol = 1e-8;
    double err1 = log(CumulativeValue(x1)) - log(u);
    double err2 = log(CumulativeValue(x2)) - log(u);
    while (err1 > 0)
    {
        x1 /= 2;
        err1 = log(CumulativeValue(x1)) - log(u);
    }

    while (err2 < 0)
    {
        x2 *= 2;
        err2 = log(CumulativeValue(x2)) - log(u);
    }

    while (min(fabs(err1),fabs(err2)) > tol && fabs(x1-x2)>tol)
    {
        double slope = (err2 - err1) / (log(x2) - log(x1));
        double x_p = exp(log(x1) - err1 / slope);
        double err_p = log(CumulativeValue(x_p)) - log(u);
        if (err_p > 0)
        {
            x2 = x_p;
            err2 = log(CumulativeValue(x2)) - log(u);
        }
        else
        {
            x1 = x_p;
            err1 = log(CumulativeValue(x1)) - log(u);
        }

    }
    if (fabs(err1) > fabs(err2))
        return x2;
    else
        return x1;
}

double CDistribution::CumulativeValue(double x)
{
    double out = 0;
    int m = NumberOfCoreParameters(name)+1;
    int n = params.size() / (NumberOfCoreParameters(name)+1);
    if (this->name == "lognormal")
    {
        for (int i = 0; i < n; i++)
        {
            double _erf = erf((log(x) - params[m * i + 1]) / (sqrt(2.0)*params[m * i + 2]));
            out += params[m * i] * 0.5*(1.0 + _erf);

        }
    }
    else
    {

    }
    return out;

}
