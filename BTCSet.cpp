#include "BTCSet.h"
#include "string.h"
#include <fstream>
#include "Distribution.h"
#include <iostream>
#include "Vector.h"
#include "Matrix.h"
#include "omp.h"
#if Qt_version
#include <qdebug.h>
#endif

//#define GNUplot

using namespace std;

CTimeSeriesSet::CTimeSeriesSet(void)
{
	nvars = 0;
	BTC.resize(nvars);
	unif = true;
}

CTimeSeriesSet::~CTimeSeriesSet(void)
{

}

CTimeSeriesSet::CTimeSeriesSet(int n)
{
	nvars = n;
	BTC.resize(nvars);
	names.resize(nvars);
	for (int i=0; i<nvars; i++) BTC[i] = CBTC();
	unif = true;
}
CTimeSeriesSet::CTimeSeriesSet(vector<vector<double> > &data, int writeInterval) :CTimeSeriesSet(data[0].size())
{
	for (int i = 0; i<int(data.size()); i++)
		if (i%writeInterval == 0) append(i, data[i]);
}
CTimeSeriesSet::CTimeSeriesSet(int numberofBTCs, int sizeofBTCvector)
{
	nvars = numberofBTCs;
	BTC.resize(nvars);
	names.resize(nvars);
	for (int i = 0; i<nvars; i++)
		BTC[i] = CBTC(sizeofBTCvector);
	unif = true;
}

CTimeSeriesSet merge(CTimeSeriesSet A, CTimeSeriesSet &B)
{
	CTimeSeriesSet C = A;
	for (int i=0; i<B.nvars; i++)
	{	if (int(B.names.size())>i) C.names.push_back(B.names[i]);
		C.BTC.push_back(B.BTC[i]);
		C.nvars++;
	}
	return C;
}

CTimeSeriesSet merge(vector<CTimeSeriesSet> &A)
{
	CTimeSeriesSet C = A[0];
	for (int i=1; i<int(A.size()); i++)
		C = merge(C,A[i]);
	return C;
}

void CTimeSeriesSet::writetofile(char outputfile[])
{
	FILE *Fil;
	Fil = fopen(outputfile, "w");
#ifndef GNUplot
        for (int i=0; i<int(names.size()); i++)
		fprintf(Fil , "t, %s, ", names[i].c_str());
	fprintf(Fil, "\n");
	for (int j=0; j<maxnumpoints(); j++)
	{
		for (int i=0; i<nvars; i++)
		{
			//MM
			//if (i<BTC[i].n)
			if (i<BTC[i].n+1)
				fprintf(Fil, "%lf, %le,", BTC[i].t[j], BTC[i].C[j]);
			else
				fprintf(Fil, ", ,");

		}
		fprintf(Fil, "\n");
	}
#else
        fprintf(Fil , "#");
        for (int i=0; i<int(names.size()); i++)
		fprintf(Fil , " t \t %s ", names[i].c_str());
	fprintf(Fil, "\n");
	for (int j=0; j<maxnumpoints(); j++)
	{
		for (int i=0; i<nvars; i++)
		{
			//MM
			//if (i<BTC[i].n)
			if (i<BTC[i].n+1)
				fprintf(Fil, "%lf\t %le\t", BTC[i].t[j], BTC[i].C[j]);
			else
				fprintf(Fil, "\t \t");

		}
		fprintf(Fil, "\n");
	}
#endif
	fclose(Fil);

}

void CTimeSeriesSet::writetofile(string outputfile, bool writeColumnNameHeaders)
{
	FILE *Fil;
	Fil = fopen(outputfile.c_str() , "w");
#ifndef GNUplot
        if (writeColumnNameHeaders)
	{
		fprintf(Fil, "names, ");
		for (int i = 0; i < int(names.size()); i++)
			fprintf(Fil, "%s, ", names[i].c_str());
		fprintf(Fil, "\n");
	}
	fprintf(Fil , "//");
	for (int i=0; i<int(names.size()); i++)
            fprintf(Fil , "t, %s, ", names[i].c_str());
	fprintf(Fil, "\n");
	for (int j = 0; j<maxnumpoints(); j++)
	{
		for (int i=0; i<nvars; i++)
		{
			if (j<BTC[i].n)
				fprintf(Fil, "%le, %le,", BTC[i].t[j], BTC[i].C[j]);
			else
				fprintf(Fil, ", ,");

		}
		fprintf(Fil, "\n");
	}
#else
        if (writeColumnNameHeaders)
	{
		fprintf(Fil, "# names\t ");
		for (int i = 0; i < int(names.size()); i++)
			fprintf(Fil, "# %s\t ", names[i].c_str());
		fprintf(Fil, "\n");
	}
	fprintf(Fil , "#");
	for (int i=0; i<int(names.size()); i++)
            fprintf(Fil , " t \t %s\t ", names[i].c_str());
	fprintf(Fil, "\n");
	for (int j = 0; j<maxnumpoints(); j++)
	{
		for (int i=0; i<nvars; i++)
		{
			if (j<BTC[i].n)
				fprintf(Fil, "%lf\t %le\t", BTC[i].t[j], BTC[i].C[j]);
			else
				fprintf(Fil, "\t \t");

		}
		fprintf(Fil, "\n");
	}
#endif
	fclose(Fil);

}

void CTimeSeriesSet::writetofile(string outputfile, int outputwriteinterval, bool onlyonetimecolumn)
{
	FILE *Fil;
	Fil = fopen(outputfile.c_str() , "w");

#ifndef GNUplot
        if (onlyonetimecolumn)
                fprintf(Fil , "t, ");
        for (int i=0; i<int(names.size()); i++)
        {   if (onlyonetimecolumn)
                fprintf(Fil , "%s, ", names[i].c_str());
            else
                fprintf(Fil , "t, %s, ", names[i].c_str());
        }
	fprintf(Fil, "\n");
	for (int j=0; j<maxnumpoints(); j++)
	{
		int column_to_be_written = 0;
		bool time_written=false;
		for (int i=0; i<nvars; i++)
		{
			if (onlyonetimecolumn)
            {
                if (j<BTC[i].n & !time_written)
                {
                    fprintf(Fil, "%lf, ", BTC[column_to_be_written].t[j]);
                    time_written=true;
                }
            }
            column_to_be_written++;
		}
        for (int i=0; i<nvars; i++)
		{
			if (j%outputwriteinterval==0)
			{	if (j<BTC[i].n)
				{
					if (onlyonetimecolumn)
                        fprintf(Fil, "%lf, ", BTC[i].C[j]);
					else
                        fprintf(Fil, "%lf, %le,", BTC[i].t[j], BTC[i].C[j]);
				}
				else
					if (!onlyonetimecolumn)
                        fprintf(Fil, ", ,");
                    else
                        fprintf(Fil, ", ");
            }

		}
		if (j%outputwriteinterval==0)
			fprintf(Fil, "\n");
	}
#else
        for (int i=0; i<int(names.size()); i++)
		fprintf(Fil , "# t \t %s\t ", names[i].c_str());
	fprintf(Fil, "\n");
	for (int j=0; j<maxnumpoints(); j++)
	{
		for (int i=0; i<nvars; i++)
		{
			if (j%outputwriteinterval==0)
			{	if (j<BTC[i].n)
					fprintf(Fil, "%lf\t %le\t", BTC[i].t[j], BTC[i].C[j]);
				else
					fprintf(Fil, "\t \t");
                        }

		}
		if (j%outputwriteinterval==0)
			fprintf(Fil, "\n");
	}
#endif
	fclose(Fil);

}

int CTimeSeriesSet::maxnumpoints()
{
	int m = 0;
	for (int i=0; i<nvars; i++)
		if (BTC[i].n>m) m = BTC[i].n;

	return m;
}

CTimeSeriesSet::CTimeSeriesSet(const CTimeSeriesSet &B)
{
	nvars = B.nvars;
	BTC.resize(nvars);
	names = B.names;
	BTC = B.BTC;
	unif = B.unif;

}

CTimeSeriesSet::CTimeSeriesSet(const CBTC &B)
{
	nvars = 1;
	BTC.resize(1);

	BTC[0] = B;
}

CTimeSeriesSet& CTimeSeriesSet::operator = (const CTimeSeriesSet &B)
{
	nvars = B.nvars;
	BTC.resize(nvars);
	names = B.names;
	for (int i=0; i<nvars; i++)
		BTC[i] = B.BTC[i];
	unif = B.unif;
	return *this;

}

vector<double> CTimeSeriesSet::interpolate(double t)
{
	vector<double> out;
	out.resize(nvars);
	for (int i=0; i<nvars; i++)
		out[i] = BTC[i].interpol(t);

	return out;
}

vector<double> CTimeSeriesSet::interpolate(double t, int fnvars)
{
	vector<double> out;
	out.resize(fnvars);
	for (int i=0; i<min(nvars,fnvars); i++)
		out[i] = BTC[i].interpol(t);

	return out;
}

CTimeSeriesSet::CTimeSeriesSet(string filename, bool varytime)
{
	unif = false;
	vector<string> units;
	ifstream file(filename.c_str());
	vector<string> s;
	nvars = 0;
	if (file.good() == false)
	{
		file_not_found = true;
		return;
	}

        if (varytime == false)
            while (file.eof() == false)
            {
                s = aquiutils::getline(file);
                if (s.size()>0)
                {
                    if (s[0] == "names")
                            for (int i = 1; i < int(s.size()); i++) names.push_back(s[i]);
                    if (s[0] == "units")
                            for (int i = 1; i < int(s.size()); i++) units.push_back(s[i]);
                    if ((s[0].substr(0, 2) != "//") && (s[0] != "names") && (s[0] != "units"))
                    {
                        if (nvars == 0) { nvars = int(s.size()) - 1; BTC.resize(nvars); }
                        if (int(s.size()) == nvars + 1)
                            for (int i = 0; i < nvars; i++)
                            {
                                BTC[i].t.push_back(atof(s[0].c_str()));
                                BTC[i].C.push_back(atof(s[i + 1].c_str()));
                                BTC[i].n++;
                                if (BTC[i].t.size()>2)
                                        if ((BTC[i].t[BTC[i].t.size() - 1] - BTC[i].t[BTC[i].t.size() - 2]) != (BTC[i].t[BTC[i].t.size() - 2] - BTC[i].t[BTC[i].t.size() - 3]))
                                                BTC[i].structured = false;

                            }

                    }
                }
            }
	else
            while (file.eof() == false)
            {
                s = aquiutils::getline(file);
                if (s.size() > 0)
                {
                    if (s[0] == "names")
                            for (int i = 1; i < int(s.size()); i++) names.push_back(s[i]);
                    if (s[0] == "units")
                            for (int i = 1; i < int(s.size()); i++) units.push_back(s[i]);
                    if ((s[0].substr(0, 2) != "//") && (s[0] != "names") && (s[0] != "units"))
                    {
                        if (nvars == 0) { nvars = s.size() / 2; BTC.resize(nvars); }

                        for (int i = 0; i < nvars; i++)
                        {
                            if (int(s.size()) >= 2 * (i + 1))
                                if ((aquiutils::trim(s[2 * i]) != "") && (aquiutils::trim(s[2 * i + 1]) != ""))
                                {
                                    BTC[i].t.push_back(atof(s[2 * i].c_str()));
                                    BTC[i].C.push_back(atof(s[2 * i + 1].c_str()));
                                    BTC[i].n++;
                                    if (BTC[i].t.size()>2)
                                        if ((BTC[i].t[BTC[i].t.size() - 1] - BTC[i].t[BTC[i].t.size() - 2]) != (BTC[i].t[BTC[i].t.size() - 2] - BTC[i].t[BTC[i].t.size() - 3]))
                                            BTC[i].structured = false;
                                }
                        }
                    }
                }
            }



        file.close();
	for (int i = 0; i < min(int(names.size()), nvars); i++)
		BTC[i].name = names[i];

	for (int i = 0; i < min(int(units.size()), nvars); i++)
		BTC[i].unit = units[i];

	//for (int i=0; i<nvars; i++)
	//	BTC[i].assign_D();

	if (int(names.size()) < nvars)
	{
		names.resize(nvars);
	}
	if (nvars == 1 && names[0] == "")
		names[0] = "Data";
	if (nvars > 1)
		for (int i = 0; i < nvars; i++)
			if (names[i] == "")
				names[i] = "Data (" + i + ')';
}




void CTimeSeriesSet::getfromfile(string filename, bool varytime)
{
	unif = false;
	vector<string> units;
	ifstream file(filename.c_str());
	vector<string> s;
	nvars = 0;
	if (varytime==false)
		while (file.eof()== false)
		{
			s = aquiutils::getline(file,',');
			if (s.size()>0)
			{
                            if (s[0] == "names")
                                    for (int i = 1; i < int(s.size()); i++) names.push_back(s[i]);
                            if (s[0] == "units")
                                    for (int i = 1; i < int(s.size()); i++) units.push_back(s[i]);
                            if ((s[0].substr(0, 2) != "//") && (s[0] != "names") && (s[0] != "units"))
                            {
                                if (nvars==0) {nvars = s.size()-1; BTC.resize(nvars); for (int i=0; i<nvars; i++) BTC[i].structured = true;}
                                if (int(s.size())==nvars+1)
                                    for (int i=0; i<nvars; i++)
                                    {
                                        BTC[i].t.push_back(atof(s[0].c_str()));
                                        BTC[i].C.push_back(atof(s[i+1].c_str()));
                                        BTC[i].n++;
                                        if (BTC[i].n>2)
                                                        if ((BTC[i].t[BTC[i].n-1]-BTC[i].t[BTC[i].n-2]) != (BTC[i].t[BTC[i].n-2]-BTC[i].t[BTC[i].n-3]))
                                                                BTC[i].structured = false;
                                    }

				}
			}
		}
	else
		while (file.eof()== false)
		{
			s = aquiutils::getline(file,',');
			if (s.size() > 0)
			{
				if (s[0] == "names")
					for (int i = 1; i < int(s.size()); i++) names.push_back(s[i]);
                                if (aquiutils::trim(s[0]).substr(0,1)=="#")
                                    for (int i=1; i<int(s.size()); i+=2) names.push_back(aquiutils::trim(s[i]));
				if (s[0] == "units")
					for (int i = 1; i < int(s.size()); i++) units.push_back(s[i]);
				if ((s[0].substr(0, 2) != "//") && (s[0] != "names") && (s[0] != "units") && (aquiutils::trim(s[0]).substr(0,1)!="#"))			{
					if (nvars == 0) { nvars = s.size() / 2; BTC.resize(nvars); for (int i = 0; i < nvars; i++) BTC[i].structured = true; }

					int n_line = s.size() / 2;
					for (int i = 0; i < n_line; i++)
					{
						if ((aquiutils::trim(s[2 * i]) != "") && (aquiutils::trim(s[2 * i + 1]) != ""))
						{
							BTC[i].t.push_back(atof(s[2 * i].c_str()));
							BTC[i].C.push_back(atof(s[2 * i + 1].c_str()));
							BTC[i].n++;
							if (BTC[i].n>2)
								if ((BTC[i].t[BTC[i].n - 1] - BTC[i].t[BTC[i].n - 2]) != (BTC[i].t[BTC[i].n - 2] - BTC[i].t[BTC[i].n - 3]))
									BTC[i].structured = false;
						}
					}
				}
			}
		}

	for (int i = 0; i < min(int(names.size()), nvars); i++)
		BTC[i].name = names[i];

	for (int i = 0; i < min(int(units.size()), nvars); i++)
		BTC[i].unit = units[i];

	file.close();
}

double CTimeSeriesSet::maxtime()
{
	return BTC[0].t[BTC[0].n-1];

}

double CTimeSeriesSet::mintime()
{
	return BTC[0].t[0];

}

double diff(CTimeSeriesSet B1, CTimeSeriesSet B2)
{
	double sum = 0;
	for (int i=0; i<B1.nvars; i++)
		sum += diff(B1.BTC[i],B2.BTC[i]);

	return sum;

}

CTimeSeriesSet operator * (const CTimeSeriesSet &BTC, const double &C)
{
	CTimeSeriesSet A = BTC;
	A.BTC[0] = A.BTC[0]*C;
	return A;
}

vector<double> CTimeSeriesSet::getrandom()
{
	int a = int(GetRndUniF(0,BTC[0].n));
	vector<double> res(nvars);
	for (int i=0; i<nvars; i++)
		res[i] = BTC[i].C[a];

	return res;
}

vector<double> CTimeSeriesSet::getrandom(int burnin)
{
	int a = int(GetRndUniF(0,BTC[0].n-burnin));
	vector<double> res(nvars);
	for (int i=0; i<nvars; i++)
		res[i] = BTC[i].C[a+burnin];

	return res;
}


vector<double> CTimeSeriesSet::getrow(int a)
{

	vector<double> res(nvars);
	for (int i = 0; i<nvars; i++)
		res[i] = BTC[i].C[a];

	return res;
}

vector<double> CTimeSeriesSet::percentile(double x)
{
	vector<double> v;
	for (int i=0; i<nvars; i++)
		v.push_back(BTC[i].percentile(x));

	return v;
}

vector<double> CTimeSeriesSet::mean(int limit)
{
	vector<double> v;
	for (int i=0; i<nvars; i++)
		v.push_back(BTC[i].mean(limit));
	return v;

}

vector<double> CTimeSeriesSet::mean(int limit, vector<int> index)
{
	vector<double> v;
	for (int i = 0; i<int(index.size()); i++)
		v.push_back(BTC[index[i]].mean(limit));
	return v;

}

vector<double> CTimeSeriesSet::std(int limit)
{
	vector<double> v;
	for (int i=0; i<nvars; i++)
		v.push_back(BTC[i].std(limit));
	return v;

}

vector<double> CTimeSeriesSet::std(int limit, vector<int> index)
{
	vector<double> v;
	for (int i = 0; i<int(index.size()); i++)
		v.push_back(BTC[index[i]].std(limit));
	return v;

}

CMatrix CTimeSeriesSet::correlation(int limit, int n)
{
	CMatrix r_xy(n);

	for (int i=0; i<n; i++)
		for (int j=0; j<=i; j++)
			r_xy[i][j] = R(BTC[i], BTC[j], limit);

	return r_xy;

}

vector<double> CTimeSeriesSet::average()
{
	vector<double> v;
	for (int i=0; i<nvars; i++)
		v.push_back(BTC[i].average());
	return v;

}

vector<double> CTimeSeriesSet::integrate()
{
	vector<double> v;
	for (int i=0; i<nvars; i++)
		v.push_back(BTC[i].integrate());
	return v;

}

vector<double> CTimeSeriesSet::percentile(double x, int limit)
{
	vector<double> v;
	for (int i=0; i<nvars; i++)
		v.push_back(BTC[i].percentile(x,limit));

	return v;
}

vector<double> CTimeSeriesSet::percentile(double x, int limit, vector<int> index)
{
	vector<double> v;
	for (int i = 0; i<int(index.size()); i++)
		v.push_back(BTC[index[i]].percentile(x, limit));

	return v;
}

CTimeSeriesSet CTimeSeriesSet::sort(int burnOut, bool reverse_cumulative)
{
	CTimeSeriesSet r(nvars);
	if (burnOut < 0)
		burnOut = 0;
	vector<vector<double> > temp;
	temp.resize(nvars);
	vector<double> tempVec;

	int counter = 0;
	clock_t tt0 = clock();

#pragma omp parallel

	#pragma omp for
	for (int i = 0; i < nvars; i++)
	{
		counter++;

		clock_t t0 = clock();
//		r.BTC[i].C.resize(BTC[i].n - burnOut);
		tempVec.resize(BTC[i].n - burnOut);
		for (int j = 0; j < int(tempVec.size()); j++)
			tempVec[j] = BTC[i].C[j + burnOut];

		temp[i] = QSort(tempVec);
		//r.BTC[i].C = QSort(temp);
		clock_t t1 = clock() - t0;
		float run_time = ((float)t1) / CLOCKS_PER_SEC;

	}
	for (int i = 0; i < nvars; i++)
	{
		r.BTC[i].C.resize(BTC[i].n - burnOut);
		r.BTC[i].C = temp[i];
		r.BTC[i].n = temp[i].size();
		r.BTC[i].t.resize(r.BTC[i].n);
		for (int j=0; j<r.BTC[i].n; j++)
        {
            if (!reverse_cumulative)
                r.BTC[i].t[j] = double(j+1)/double(r.BTC[i].n);
            else
                r.BTC[i].t[j] = 1-double(j)/double(r.BTC[i].n);
        }
	}
	clock_t tt1 = clock() - tt0;
	float run_time = ((float)tt1) / CLOCKS_PER_SEC;


	return r;
}



CTimeSeriesSet CTimeSeriesSet::detivative()
{
	CTimeSeriesSet out(nvars);
	out.names = names;
	for (int i = 0; i < int(BTC.size()); i++)
		out.BTC[i] = BTC[i].derivative();

	return out;
}
CTimeSeriesSet CTimeSeriesSet::distribution(int n_bins, int n_columns, int limit, double smoothing_factor)
{

    CTimeSeriesSet A(n_columns);
    for (int i = 0; i < n_columns; i++)
    {
        A.BTC[i] = BTC[i].distribution(n_bins, (BTC[i].maxC()-BTC[i].minC())*smoothing_factor, limit);
        A.names[i] = names[i];

    }

    return A;
}

CVector norm2dif(CTimeSeriesSet &A, CTimeSeriesSet &B)
{
	CVector res;
	for (int i=0; i<min(A.nvars,B.nvars); i++)
	{   CBTC AA = A.BTC[i].Log(1e-5);
            CBTC BB = B.BTC[i].Log(1e-5);
            double x = diff_abs(AA,BB)/B.BTC[i].n;
            res.append(x);

        }

	return res;

}

void CTimeSeriesSet::append(double t, vector<double> c, double weight)
{
    for (int i=0; i<min(int(c.size()), nvars); i++)
    {   BTC[i].structured = true;
        if (weight!=1.0) BTC[i].weighted = true;
        BTC[i].append(t,c[i],weight);
        if (BTC[i].n>2)
            if ((BTC[i].t[BTC[i].n-1]-BTC[i].t[BTC[i].n-2]) != (BTC[i].t[BTC[i].n-2]-BTC[i].t[BTC[i].n-3]))
                BTC[i].structured = false;
    }
}

CBTC CTimeSeriesSet::add(vector<int> ii)
{
	CBTC A = BTC[ii[0]];
	A.structured = BTC[ii[0]].structured;
	for (int i=1; i<int(ii.size()); i++)
	if (unif==false)
	{	A+=BTC[ii[i]];
		A.structured = (A.structured && BTC[ii[i]].structured);
	}
	else
	{	A%=BTC[ii[i]];
		A.structured = (A.structured && BTC[ii[i]].structured);
	}

	return A;
}

CBTC CTimeSeriesSet::add_mult(vector<int> ii, vector<double> mult)
{
	CBTC A;
	if (ii.size()>0)
	{	A = mult[0]*BTC[ii[0]];
		A.structured = BTC[ii[0]].structured;
		for (int i=1; i<int(ii.size()); i++)
		if (unif==false)
		{   CBTC AA = mult[i]*BTC[ii[i]];
                    A+=AA;
                    A.structured = (A.structured && BTC[ii[i]].structured);
		}
		else
		{
                        CBTC AA = mult[i]*BTC[ii[i]];
                        A%=AA;
			A.structured = (A.structured && BTC[ii[i]].structured);
		}
	}
	else
	{
		A.setnumpoints(2);
		A.t[0] = mintime();
		A.t[1] = maxtime();
	}
	return A;
}

CBTC CTimeSeriesSet::add_mult(vector<int> ii, CTimeSeriesSet &mult)
{
	CBTC A;
	if (ii.size()>0)
	{	A = mult.BTC[0]*BTC[ii[0]];
		A.structured = BTC[ii[0]].structured;
		for (int i=1; i<int(ii.size()); i++)
		if (unif==false)
		{   CBTC AA = BTC[ii[i]]*mult.BTC[i];
                    A+=AA;
                    A.structured = (A.structured && BTC[ii[i]].structured && mult.BTC[i].structured);
		}
		else
		{   CBTC AA = BTC[ii[i]]*mult.BTC[i];
                    A%=AA;
                    A.structured = (A.structured && BTC[ii[i]].structured && mult.BTC[i].structured);
		}
	}
	else
	{
		A.setnumpoints(2);
		A.t[0] = mintime();
		A.t[1] = maxtime();
	}
	return A;
}


CBTC CTimeSeriesSet::divide(int ii, int jj)
{
	CBTC A;
	A.structured = (BTC[ii].structured && BTC[jj].structured);
	if (unif==false)
		A=BTC[ii]/BTC[jj];
	else
		A=BTC[ii]%BTC[jj];

	return A;


}

CTimeSeriesSet CTimeSeriesSet::make_uniform(double increment)
{
	if (nvars==0) return CTimeSeriesSet();
	CTimeSeriesSet out(nvars);
	out.names = names;

	if (unif == true)
	{
		for (int i = 0; i < nvars; i++)
		{
			out.BTC[i].append(BTC[i].t[0], BTC[i].C[0]);
			if (BTC[i].D.size() == 0) BTC[i].assign_D();
		}
		for (int i=0; i<BTC[0].n-1; i++)
		{
			int i1 = int((BTC[0].t[i]-BTC[0].t[0])/increment);
			int i2 = int((BTC[0].t[i+1]-BTC[0].t[0])/increment);
			for (int j=i1+1; j<=i2; j++)
			{
				double x = j*increment+BTC[0].t[0];
				for (int k=0; k<nvars; k++)
				{
					double CC = (x-BTC[k].t[i])/(BTC[k].t[i+1]-BTC[k].t[i])*(BTC[k].C[i+1]-BTC[k].C[i])+BTC[k].C[i];
					double DD = (x - BTC[k].t[i]) / (BTC[k].t[i + 1] - BTC[k].t[i])*(BTC[k].D[i + 1] - BTC[k].D[i]) + BTC[k].D[i];
					out.BTC[k].append(x,CC);
					out.BTC[k].D.push_back(DD);
				}
			}
		}
	}
	else
	{
		for (int k=0; k<nvars; k++)
			out.BTC[k] = BTC[k].make_uniform(increment);
	}
	for (int k=0; k<nvars; k++)
		out.BTC[k].structured = true;

	out.unif=true;

	return out;

}

CTimeSeriesSet CTimeSeriesSet::getpercentiles(vector<double> percents)
{
	CTimeSeriesSet X(1+percents.size());

	X.names.clear();
	X.setname(0, "Mean");
	for (int j=0; j<int(percents.size()); j++)
	{
		string Xname = aquiutils::numbertostring(percents[j]*100,false)+" %";
		X.setname(j + 1, Xname);
	}

	vector<double> XX(1+percents.size());
	vector<double> XX_prc(percents.size());

	double meanX;
	for (int i=0; i<BTC[0].n; i++)
	{
		vector<double> x;
		int count = 0;
		for (int j=0; j<nvars; j++)
			if (i<BTC[j].n)
			{	x.push_back(BTC[j].C[i]);
				count++;
			}

		meanX = CVector(x).sum()/count;

		XX[0] = meanX;
		XX_prc = prcntl(x,percents);
		for (int j=0; j<int(percents.size()); j++)
			XX[j+1] = XX_prc[j];

		X.append(BTC[0].t[i],XX);
	}

	return X;
}

CVector CTimeSeriesSet::out_of_limit(double limit)
{
	CVector v(nvars);
	for (int i=0; i<nvars; i++)
	{
		double n_tot = BTC[i].n;
		double n_exceed = 0;
		for (int j=0; j<BTC[i].n; j++)
		{		if (BTC[i].C[j] > limit)
				n_exceed++;		}

		v[i] = n_exceed/n_tot;
	}

	return v;
}

CTimeSeriesSet CTimeSeriesSet::add_noise(vector<double> std, bool logd)
{
	CTimeSeriesSet X(nvars);
	for (int i=0; i<nvars; i++)
		X.BTC[i] = BTC[i].add_noise(std[i],logd);

	return X;
}

CVector sum_interpolate(vector<CTimeSeriesSet> &BTC, double t)
{
	if (BTC.size()==0) return CVector(1);
	CVector sum(max(max_n_vars(BTC),2)); //Be chacked later?
	for (int i=0; i<int(BTC.size()); i++)   //We can have several BTCs (Precipitation, Evaporation,...) and each one can have several variables (flow, concentration, ...)
	{
		for (int j=0; j<BTC[i].nvars; j++)
			sum[j]+=BTC[i].BTC[j].interpol(t);
	}
	return sum;
}

double sum_interpolate(vector<CTimeSeriesSet> &BTC, double t, string name)
{
	if (BTC.size() == 0) return 0;
	double sum=0;
	for (int i = 0; i<int(BTC.size()); i++)
	{
		int ii = BTC[i].lookup(name);
		if (ii!=-1)
			sum += BTC[i].BTC[ii].interpol(t);
	}
	return sum;
}

void CTimeSeriesSet::clear()
{
	BTC.clear();
        names.clear();
	nvars = 0;
}

vector<double> CTimeSeriesSet::max_wiggle()
{
	double max_wig=0;
	int wiggle_id=-1;
	for (int i=0; i<nvars; i++)
	{	double a = BTC[i].wiggle();
		if (a>max_wig) wiggle_id = i;
		max_wig = max(max_wig,a);

	}
	vector<double> out;
	out.push_back(max_wig);
	out.push_back(wiggle_id);
	return out;
}

vector<double> CTimeSeriesSet::max_wiggle_corr(int _n)
{
	double max_wig = 0;
	int wiggle_id = -1;
	for (int i = 0; i<nvars; i++)
	{
		double a = exp(-5*(1+BTC[i].wiggle_corr(_n)));
		if (a>max_wig) wiggle_id = i;
		max_wig = max(max_wig, a);

	}
	vector<double> out;
	out.push_back(max_wig);
	out.push_back(wiggle_id);
	return out;
}

vector<int> CTimeSeriesSet::max_wiggle_sl(int ii, double tol)
{
	double max_wig = 0;
	int wiggle_id = -1;
	for (int i = 0; i<min(ii,nvars); i++)
	{
		int a = int(BTC[i].wiggle_sl(tol));
		if (a==1) wiggle_id = i;
		max_wig = max_wig || a;

	}
	vector<int> out;
	out.push_back(max_wig);
	out.push_back(wiggle_id);
	return out;
}

int max_n_vars(vector<CTimeSeriesSet> &BTC)
{
	int k = 0;
	for (int i=0; i<int(BTC.size()); i++)
	{
		if (BTC[i].nvars>k) k = BTC[i].nvars;
	}
	return k;
}

void CTimeSeriesSet::knockout(double t)
{
	for (int i=0; i<nvars; i++) BTC[i].knock_out(t);
}

int CTimeSeriesSet::lookup(string S)
{
	int out = -1;
	for (int i = 0; i < int(names.size()); i++)
		if (aquiutils::tolower(S) == aquiutils::tolower(names[i]))
		{
			out = i;
			return out;
		}

	for (int i = 0; i < int(BTC.size()); i++)
		if (aquiutils::tolower(S) == aquiutils::tolower(BTC[i].name))
		{
			out = i;
			return out;
		}


	return out;
}
CBTC &CTimeSeriesSet::operator[](int index)
{
	while (int(names.size()) < index+1)
		names.push_back("");

	if (BTC[index].name == "")
		if (names[index] != "")
			BTC[index].name = names[index];
	if (names[index] == "")
		if (BTC[index].name != "")
			names[index] = BTC[index].name;
	return BTC[index];
}
CBTC &CTimeSeriesSet::operator[](string BTCName)
{
	if (lookup(BTCName) != -1)
            return BTC[lookup(BTCName)];
	else
	{   CBTC AA;
            return AA;
        }
}

void CTimeSeriesSet::pushBackName(string name)
{
	names.push_back(name);
	BTC[names.size() - 1].name = name;
}
void CTimeSeriesSet::append(CBTC &_BTC, string name)
{
	if (name != "")
            _BTC.name = name;
	BTC.push_back(_BTC);
	pushBackName(name);
	nvars = BTC.size();
}


void CTimeSeriesSet::append(const CBTC &_BTC, string name)
{

	BTC.push_back(_BTC);
	pushBackName(name);
	nvars = BTC.size();
	if (name != "")
        BTC[nvars-1].name = name;
}


void CTimeSeriesSet::setname(int index, string name)
{
	while (names.size() < BTC.size())
		names.push_back("");

	names[index] = name;
	BTC[index].name = name;

}


CTimeSeries CTimeSeriesSet::GetGradientDistribution(double dx, int nbins)
{
    CTimeSeries gradients;
    for (int i=0; i<BTC[0].n; i++)
    {
        gradients.append(i,(BTC[0].C[i]-BTC[1].C[i])/dx);
    }
    return gradients.distribution(nbins);

}

double CTimeSeriesSet::get_correlation()
{
    CVector X(nvars);
	double sum_prod = 0;
	double sum_2 = 0;
    cout<<"     Calculating Correlation..."<<endl;
	for (int i=0; i<BTC[0].n; i++)
	{
            for (int ii=0; ii<nvars-1; ii++)
            {   sum_prod += BTC[ii].C[i] * BTC[ii+1].C[i]*BTC[0].weight[i];
                sum_2 += (pow(BTC[ii].C[i], 2) + pow(BTC[ii+1].C[i], 2))*0.5*BTC[0].weight[i];
            }
	}
	return sum_prod / sum_2;;
}

CVector CTimeSeriesSet::get_kappa_gamma(double delta_x)
{
	CVector X(nvars);
	double sum_prod = 0;
	double sum_2 = 0;
    cout<<"     Calculating Gamma..."<<endl;
	for (int i=0; i<BTC[0].n; i++)
	{
            for (int ii=0; ii<nvars-1; ii++)
            {   sum_prod += BTC[ii].C[i] * BTC[ii+1].C[i]*BTC[0].weight[i];
                sum_2 += (pow(BTC[ii].C[i], 2) + pow(BTC[ii+1].C[i], 2))*0.5*BTC[0].weight[i];
            }
	}
	X[0] = (1.0-sum_prod / sum_2)/delta_x;
	double err = 0;
    double sum_weight = 0;
    cout<<"     Calculating Kappa..."<<endl;
    for (int i = 0; i<BTC[0].n; i++)
    {
        for (int ii=0; ii<nvars-1; ii++)
        {   err += (pow(BTC[ii].C[i] - (sum_prod / sum_2) * BTC[ii+1].C[i], 2)+pow(BTC[ii+1].C[i] - (sum_prod / sum_2) * BTC[ii].C[i], 2))*BTC[0].weight[i]*0.5;
            sum_weight += BTC[ii].weight[i];
        }
    }
    X[1] = err /(2.0*delta_x)/sum_weight;

    if (nvars>2)
    {
        cout<<"     Calculating Error Correlation..."<<endl;
        CBTCSet residuals(nvars-1);
        for (int i = 0; i<BTC[0].n; i++)
        {
            {
                vector<double> res;
                for (int j=0; j<nvars-1;j++)
                    res.push_back(BTC[j].C[i] - (sum_prod / sum_2) * BTC[j+1].C[i]);

                residuals.append(i,res);

            }
        }
        cout<<"     Calculating Correlation..."<<endl;
        if (nvars>2)
            X[2] = TS::correlation(residuals.BTC[0], residuals.BTC[1]);


    }



	return X;
}

void CTimeSeriesSet::append(string name, double t, double c, double weight)
{
    if (lookup(name)!=-1)
    {   if (weight!=1.0) BTC[lookup(name)].weighted = true;
        BTC[lookup(name)].append(t,c,weight);

    }
    else
    {
        CBTC BTC;
        BTC.name = name;
        if (weight!=1.0) BTC.weighted = true;
        BTC.append(t,c,weight);
        append(BTC,name);
    }
}

CTimeSeriesSet CTimeSeriesSet::getcummulative()
{
    CTimeSeriesSet out;
    for (int i=0; i<nvars; i++)
    {
        CTimeSeries B = BTC[i].getcummulative();
        out.append(B);
        out.setname(i,names[i]);
    }
    return out;
}


CTimeSeriesSet CTimeSeriesSet::Transpose(const double &dt, const string &column_name)
{
    if (nvars == 0) return CTimeSeries();
    CTimeSeriesSet transposed(BTC[0].n);
    for (int i=0; i<nvars; i++)
    {
        for (int j=0; j<BTC[i].n; j++)
        {
            transposed.BTC[j].append(dt*i,BTC[i].C[j]);
        }
    }
    for (int i=0; i<transposed.nvars; i++)
        transposed.setname(i, aquiutils::numbertostring(BTC[0].t[i],false));

    return transposed;

}

double CTimeSeriesSet::FrankCopulaLogLikelihood(const double &alpha)
{
    double sum=0;
    double n = BTC[0].n;
    sum += n*(log(alpha)+2*alpha)+n*log(1-exp(-alpha))-alpha*(BTC[0].sum()+BTC[1].sum());
    for (int i=0; i<BTC[0].n; i++)
    {
        sum+= -2*log(exp(-alpha*(BTC[0].C[i]-1))+exp(-alpha*(BTC[1].C[i]-1))-exp(-alpha*(BTC[0].C[i]+BTC[1].C[i]-1))-1);
    }
    return sum/n;

}

double CTimeSeriesSet::FrankCopulaLogLikelihood_deriv(const double &alpha)
{
    double sum=0;
    double n = BTC[0].n;
    sum += (1.0/alpha+2.0)+ exp(-alpha)/(1-exp(-alpha))-(BTC[0].sum()+BTC[1].sum())/n;
    for (int i=0; i<BTC[0].n; i++)
    {
        double numerator = -(BTC[0].C[i]-1)*exp(-alpha*(BTC[0].C[i]-1))-(BTC[1].C[i]-1)*exp(-alpha*(BTC[1].C[i]-1))+(BTC[0].C[i]+BTC[1].C[i]-1)*exp(-alpha*(BTC[0].C[i]+BTC[1].C[i]-1));
        double denumerator = exp(-alpha*(BTC[0].C[i]-1))+exp(-alpha*(BTC[1].C[i]-1))-exp(-alpha*(BTC[0].C[i]+BTC[1].C[i]-1))-1;
        sum+= -2.0/n*numerator/denumerator;
    }

    double numerical_derivative = (FrankCopulaLogLikelihood(alpha*1.000001)-FrankCopulaLogLikelihood(alpha*1))/(alpha*0.000001);
    return sum;
}

double addup(const vector<double> x)
{
    double sum=0;
    for (int i=0; i<x.size(); i++)
        sum+=x[i];

}

double CTimeSeriesSet::FrankCopulaLogLikelihood_deriv_derive(const double &alpha)
{
    double epsilon = 1e-6;
    return (FrankCopulaLogLikelihood_deriv(alpha+epsilon) - FrankCopulaLogLikelihood_deriv(alpha))/epsilon;
}

double CTimeSeriesSet::Estimate_Frank_Alpha(double initial_guess, double error_tolerance)
{
    double err = 10000;
    double alpha = initial_guess;

    while (err>error_tolerance)
    {
        err = FrankCopulaLogLikelihood_deriv(alpha);
        double err_prime = FrankCopulaLogLikelihood_deriv_derive(alpha);
        alpha -= err/err_prime;

    }
    return alpha;
}




