#pragma once
#ifndef UTILITIES_H
#define UTILITIES_H

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>

#define SMALLNUMBER 1e-23
#define PI 3.14159265359

using namespace std;

namespace aquiutils
{

    int lookup(const vector<string> &s, const string &s1);
    int lookup(const vector<int> &s, const int &s1);
    int lookup(const vector<vector<int> > &s, const vector<int> &s1);
    int corresponding_parenthesis(string S, int i);
    int count(const string &s, const string &s1);
    bool parantheses_balance(string S);
    bool contains(const string &s, const string &s1);
    string left(const string &s, int i);
    string right(const string &s, int i);
    void remove(string &s,unsigned int i);
    void replace(string &s,unsigned int i,string p);
    void remove(string &s,unsigned int pos, unsigned int len);
    bool remove(string &s, string s1);
    void insert(string &s, unsigned int pos, string s1);
    bool replace(string &s,string s1, string s2);
    void replace(string &s,unsigned int i, unsigned int j, string p);
    bool isnumber(char S);
    bool isnumber(string S);
    bool isintegernumber(string S);
    double atof(const string &S);
    double atoi(const string &S);
    string trim(const string &s);
    void push_back(vector<string> &s, const vector<string> &s1);
    bool isnegative(const double &x);
    bool ispositive(const double &x);
    bool iszero(const double &x);
    vector<string> split(const string &s, char del);
    vector<string> getline(ifstream& file);
    vector<string> getline(ifstream& file, char del1);
    vector<vector<string>> getline_op(ifstream& file,char del1);
    vector<string> split(const string &s, const vector<char> &del);
    vector<vector<string>> getline_op(ifstream& file,vector<char> del1);
    vector<vector<string>> getline_op_eqplus(ifstream& file);
    vector<string> split_curly_semicolon(string s);
    vector<int> look_up(string s, char del);  //Returns a vector with indices of "del"
    vector<int> ATOI(vector<string> ii);
    vector<double> ATOF(vector<string> ii);
    string tolower(const string &S);
    vector<string> tolower(const vector<string> &S);
    void writeline(ofstream& f, vector<string> s, string del=",");
    void writeline(ofstream& f, vector<vector<string>> s, string del=",", string del2="&");
    void writestring(ofstream& f, string s);
    void writestring(string filename, string s);
    void writenumber(ofstream& f, double s);
    void writeendl(ofstream& f);
    double Heavyside(double x);
    double Pos(double x);
    double Neg(double x);
    string numbertostring(const double &x, bool scientific=false);
    string numbertostring(vector<double> &x, bool scientific=false);
    string numbertostring(int x);
    string numbertostring(unsigned int x);
    string numbertostring(unsigned int x, int number_of_digits);
    string numbertostring(vector<int> &x, bool scientific=false);
    string tail(std::string const& source, size_t const length);
    string tabs(int i);
    bool And(vector<bool> x);
    //double max(vector<double> x);
    //int max(vector<int> x);
    double Max(const vector<double>&);
    double Min(const vector<double>&);
    int MaxElement(const vector<double> &x);
    int MinElement(const vector<double> &x);

    int Max(vector<int>);
    string remove_backslash_r(const string &ss);
    string GetOnlyFileName(const string &fullfilename);
    double avg(double x, double y, string type="arithmetic");
    template<class T>
    T randompick(const vector<T> &vec)
    {
        unsigned int i=rand()%vec.size();
        return vec[i];
    };
    template<typename T> bool isfinite(T arg);
    template<typename T> T sum(vector<T> vec)
    {
        T out=0;
        for (int i=0; i<vec.size(); i++)
        {
            out+=vec[i];
        }
        return out;
    }
    template<typename T> bool isfinite(T arg)
    {
        if (arg==arg)
            return true;
        else
            return false;
    }
    vector<unsigned int> Rank(const vector<double> &v);
    template<typename T, typename T1, typename T2>
    int CountLessThan(const vector<T> &v, const T1 &x, T2 boundary=0, bool fromend=false)
    {
        int out = 1;
        if (boundary==0 && !fromend)
        {
            boundary = v.size();
        }
        if (!fromend)
            for (unsigned int i=0; i<boundary; i++)
            {
                if (v[i]<x)
                    out++;
            }
        else
            for (unsigned int i=boundary; i<v.size(); i++)
            {
                if (v[i]<x)
                    out++;
            }
        return out;
    }
    template<typename T, typename T1, typename T2>
    int CountGreaterThan(const vector<T> &v, const T1 &x, T2 boundary=0, bool fromend=false)
    {
        int out = 0;
        if (boundary==0 && !fromend)
        {
            boundary = v.size();
        }
        if (!fromend)
            for (unsigned int i=0; i<boundary; i++)
            {
                if (v[i]>x)
                    out++;
            }
        else
            for (unsigned int i=boundary; i<v.size(); i++)
            {
                if (v[i]>x)
                    out++;
            }
        return out;
    }



}


#endif // UTILITIES_H
