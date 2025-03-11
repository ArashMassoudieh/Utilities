#pragma once

#include <map>
#include <string>
#include <vector>

// External Libraries
#include "QuickSort.h"

#ifdef _arma
#include "armadillo"
#endif

#ifdef QT_version
#include <QList>
#include <QMap>
#include <QVariant>
#endif

// Forward Declaration
class CDistribution;

// Regression Parameters Structure
struct RegressionParameters {
    std::vector<double> parameters;
    double MSE = 0;
    double R2 = 0;
    enum class _regress_type { linear, power, exponential } regress_type;
};

// CTimeSeries Template Class
template<class T>
class CTimeSeries {
public:
    bool structured = true;
    bool error = false;
    bool file_not_found = false;
    int n = 0;
    T max_fabs = 0;

    std::string filename;
    std::string name = "";
    std::string unit = "";
    std::string defaultUnit = "";
    std::vector<std::string> unitsList;

    // Constructors & Destructor
    CTimeSeries();
    explicit CTimeSeries(int n);
    explicit CTimeSeries(const std::string& filename);
    CTimeSeries(const CTimeSeries& C);
    CTimeSeries(std::vector<T> &data, int writeInterval);
    CTimeSeries(T a, T b, const std::vector<T> &x);
    CTimeSeries(T a, T b, const CTimeSeries<T> &btc);
    CTimeSeries(const std::vector<T> &t, const std::vector<T> &C);
    ~CTimeSeries();

#ifdef _arma
    CTimeSeries(arma::mat& x, arma::mat& y);
#endif

    // File I/O
    bool readfile(const std::string&);
    bool writefile(const std::string&);

    // Interpolation Functions
    T interpol(const T& x) const;
    T interpol_D(const T& x);
    CTimeSeries interpol(const std::vector<T>& x);
    CTimeSeries interpol(const CTimeSeries& x) const;
    CTimeSeries interpol(CTimeSeries* x) const;

    // Data Manipulation
    bool SetRow(int i, const double& _t, const double& value);
    void setnumpoints(int);
    bool append(T x);
    bool append(T tt, T xx);
    void append(CTimeSeries& CC, bool continuous_time = false);
    void ResizeIfNeeded(int _increment);
    void clear();
    void adjust_size();

    // Statistical Functions
    T maxC() const;
    T minC() const;
    T maxt() const;
    T mint() const;
    T variance() const;
    T sum() const;
    T sum_squared() const;
    T percentile(T x, int limit = 0) const;
    T mean_log(int limit = 0) const;
    T mean(int limit = 0) const;
    T std(int limit = 0) const;
    T slope() const;
    std::vector<T> trend() const;

    // Integration & Aggregation
    T integrate() const;
    T integrate(T t) const;
    T integrate(T t1, T t2) const;
    T average() const;
    T average(T t) const;
    int lookupt(T t) const;
    T mean_t();

    // Signal Processing
    CTimeSeries MA_smooth(int span);
    CTimeSeries Log();
    CTimeSeries Log(T min);
    CTimeSeries Exp();
    CTimeSeries fabs();
    CTimeSeries derivative();
    CTimeSeries distribution(int n_bins, double smoothing_span=0, int limit=0);
    CTimeSeries make_uniform(T increment, T t0 = -999, bool asn_D = false);
    CTimeSeries extract(T t1, T t2);
    CTimeSeries<T> inverse_cumulative_uniform(int n_intervals = 100);
    CTimeSeries<T> LogTransformX();
    CTimeSeries<T> getcummulative();
    CTimeSeries<T> KernelSmooth(CDistribution* dist, int span = 100);
    CTimeSeries<T> KernelSmooth(CDistribution* dist, const double& span = 1);

    // Noise Handling
    CTimeSeries add_noise(T std, bool);
    void assign_D();
    T wiggle() const;
    T wiggle_corr(int _n = 10) const;
    bool wiggle_sl(T tol) const;
    T maxfabs() const;
    void knock_out(T t);
    T AutoCor1(int i = 0);

    // Regression
    RegressionParameters LinearRegress(const CTimeSeries<T>& other_series);
    RegressionParameters PowerRegress(const CTimeSeries<T>& other_series);
    CTimeSeries<T> Predict(const RegressionParameters& regression_parameters);

    // Operators
    CTimeSeries& operator=(const CTimeSeries& C);
    CTimeSeries& operator=(const double& value);
    CTimeSeries& operator+=(CTimeSeries& v);
    CTimeSeries& operator%=(CTimeSeries& v);

    // Memory Management
    bool resize(unsigned int _size);
    unsigned int Capacity() const;
    unsigned int CSize() const { return C.size(); }
    unsigned int tSize() const { return t.size(); }
    unsigned int DSize() const { return D.size(); }

    // Accessors
    T GetLastItemValue() const;
    T GetLastItemTime() const;
    T& lastD();
    T& lastC();
    T& lastt();
    bool SetT(int i, const T& value);
    bool SetC(int i, const T& value);
    bool SetD(int i, const T& value);
    T GetT(int i) const;
    T GetC(int i) const;
    T GetD(int i) const;
    void AppendD(const T& value) { D.push_back(value); }

    // Conversion
    std::vector<double> tToStdVector() const { return t; }
    std::vector<double> ValuesToStdVector() const { return C; }

    // Time-Series Processing
    void CreatePeriodicStepFunction(const T& t_start, const T& t_end, const T& duration, const T& gap, const T& magnitude);

#ifdef QT_version
    CTimeSeries(QList<QMap<QVariant, QVariant>> data);
    void compact(QDataStream& data) const;
    static CTimeSeries unCompact(QDataStream& data);
#endif

private:
    std::vector<T> t;
    std::vector<T> C;
    std::vector<T> D;
};

// Free Functions for Time Series Operations
template<class T> T diff(const CTimeSeries<T>& BTC_p, const CTimeSeries<T>& BTC_d);
template<class T> T diff_abs(const CTimeSeries<T>& BTC_p, const CTimeSeries<T>& BTC_d);
template<class T> T diff_log(const CTimeSeries<T>& BTC_p, const CTimeSeries<T>& BTC_d, double lowlim);
template<class T> T diff_norm(const CTimeSeries<T>& BTC_p, const CTimeSeries<T>& BTC_d);
template<class T> T ADD(const CTimeSeries<T>& BTC_p, const CTimeSeries<T>& BTC_d);
template<class T> T R2(const CTimeSeries<T>& BTC_p, const CTimeSeries<T>& BTC_d);
template<class T> CTimeSeries<T> operator*(T, const CTimeSeries<T>&);
template<class T> CTimeSeries<T> operator*(const CTimeSeries<T>&, double);
template<class T> CTimeSeries<T> operator/(const CTimeSeries<T>&, const CTimeSeries<T>&);
template<class T> CTimeSeries<T> operator+(const CTimeSeries<T>&, const CTimeSeries<T>&);
template<class T> CTimeSeries<T> operator-(const CTimeSeries<T>&, const CTimeSeries<T>&);
template<class T> CTimeSeries<T> operator%(const CTimeSeries<T>&, const CTimeSeries<T>&);

#include "BTC.hpp"
