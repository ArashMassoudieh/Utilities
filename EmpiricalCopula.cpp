#include "EmpiricalCopula.h"
#include "Matrix.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <numeric>

namespace
{
    struct RankedValue
    {
        double value;
        int index;
    };

    inline bool isFiniteNumber(double x)
    {
        return std::isfinite(x);
    }
}

EmpiricalCopula::EmpiricalCopula()
    : valid_(false)
{
}

void EmpiricalCopula::clear()
{
    x_.clear();
    y_.clear();
    u_.clear();
    v_.clear();
    valid_ = false;
}

bool EmpiricalCopula::fit(const std::vector<double>& x,
                          const std::vector<double>& y,
                          TieMethod tieMethod)
{
    clear();

    if (x.size() != y.size()) return false;
    if (x.empty()) return false;

    // Keep only finite paired values
    for (size_t i = 0; i < x.size(); ++i)
    {
        if (isFiniteNumber(x[i]) && isFiniteNumber(y[i]))
        {
            x_.push_back(x[i]);
            y_.push_back(y[i]);
        }
    }

    if (x_.empty()) return false;
    if (x_.size() != y_.size()) return false;

    u_ = makePseudoObservations(x_, tieMethod);
    v_ = makePseudoObservations(y_, tieMethod);

    if (u_.size() != x_.size() || v_.size() != y_.size())
    {
        clear();
        return false;
    }

    valid_ = true;
    return true;
}

bool EmpiricalCopula::isValid() const
{
    return valid_;
}

int EmpiricalCopula::size() const
{
    return static_cast<int>(u_.size());
}

const std::vector<double>& EmpiricalCopula::x() const
{
    return x_;
}

const std::vector<double>& EmpiricalCopula::y() const
{
    return y_;
}

const std::vector<double>& EmpiricalCopula::u() const
{
    return u_;
}

const std::vector<double>& EmpiricalCopula::v() const
{
    return v_;
}

double EmpiricalCopula::clampUnit(double z)
{
    if (z < 0.0) return 0.0;
    if (z > 1.0) return 1.0;
    return z;
}

std::vector<double> EmpiricalCopula::computeRanks(const std::vector<double>& values,
                                                  TieMethod tieMethod)
{
    const int n = static_cast<int>(values.size());
    std::vector<double> ranks(n, 0.0);
    if (n == 0) return ranks;

    std::vector<RankedValue> arr(n);
    for (int i = 0; i < n; ++i)
    {
        arr[i].value = values[i];
        arr[i].index = i;
    }

    std::sort(arr.begin(), arr.end(),
              [](const RankedValue& a, const RankedValue& b)
              {
                  if (a.value < b.value) return true;
                  if (a.value > b.value) return false;
                  return a.index < b.index;
              });

    int i = 0;
    while (i < n)
    {
        int j = i;
        while (j + 1 < n && arr[j + 1].value == arr[i].value)
            ++j;

        // ranks in 1..n
        if (tieMethod == AverageRank)
        {
            const double avgRank = 0.5 * ((i + 1) + (j + 1));
            for (int k = i; k <= j; ++k)
                ranks[arr[k].index] = avgRank;
        }
        else // MinRank
        {
            const double minRank = static_cast<double>(i + 1);
            for (int k = i; k <= j; ++k)
                ranks[arr[k].index] = minRank;
        }

        i = j + 1;
    }

    return ranks;
}

std::vector<double> EmpiricalCopula::makePseudoObservations(const std::vector<double>& values,
                                                            TieMethod tieMethod)
{
    const int n = static_cast<int>(values.size());
    std::vector<double> out(n, 0.0);
    if (n == 0) return out;

    const std::vector<double> ranks = computeRanks(values, tieMethod);

    // Use (r - 0.5) / n to avoid exact 0 and 1
    for (int i = 0; i < n; ++i)
        out[i] = (ranks[i] - 0.5) / static_cast<double>(n);

    return out;
}

double EmpiricalCopula::cdf(double uVal, double vVal) const
{
    if (!valid_ || u_.empty()) return 0.0;

    uVal = clampUnit(uVal);
    vVal = clampUnit(vVal);

    double count = 0.0;
    for (size_t i = 0; i < u_.size(); ++i)
    {
        if (u_[i] <= uVal && v_[i] <= vVal)
            count += 1.0;
    }

    return count / static_cast<double>(u_.size());
}

double EmpiricalCopula::upperTailCDF(double uVal, double vVal) const
{
    if (!valid_ || u_.empty()) return 0.0;

    uVal = clampUnit(uVal);
    vVal = clampUnit(vVal);

    double count = 0.0;
    for (size_t i = 0; i < u_.size(); ++i)
    {
        if (u_[i] > uVal && v_[i] > vVal)
            count += 1.0;
    }

    return count / static_cast<double>(u_.size());
}

CopulaBinnedMatrix EmpiricalCopula::makeBinnedMatrix(int nBins, bool normalize) const
{
    CopulaBinnedMatrix M(nBins);

    if (!valid_ || nBins <= 0) return M;

    M.accumulateUnitPairs(u_, v_);
    if (normalize)
        M.normalizeToUnitMass();

    return M;
}

double EmpiricalCopula::spearmanRho() const
{
    if (!valid_ || u_.size() < 2) return 0.0;

    const double mean = 0.5;
    double sxx = 0.0, syy = 0.0, sxy = 0.0;

    for (size_t i = 0; i < u_.size(); ++i)
    {
        const double du = u_[i] - mean;
        const double dv = v_[i] - mean;
        sxx += du * du;
        syy += dv * dv;
        sxy += du * dv;
    }

    if (sxx <= 0.0 || syy <= 0.0) return 0.0;
    return sxy / std::sqrt(sxx * syy);
}

double EmpiricalCopula::kendallTau() const
{
    if (!valid_ || u_.size() < 2) return 0.0;

    long long concordant = 0;
    long long discordant = 0;
    long long tiedPairs = 0;

    for (size_t i = 0; i < u_.size(); ++i)
    {
        for (size_t j = i + 1; j < u_.size(); ++j)
        {
            const double du = u_[i] - u_[j];
            const double dv = v_[i] - v_[j];
            const double prod = du * dv;

            if (prod > 0.0) ++concordant;
            else if (prod < 0.0) ++discordant;
            else ++tiedPairs;
        }
    }

    const long long denom = concordant + discordant + tiedPairs;
    if (denom == 0) return 0.0;

    // Simple tau-a style form with ties counted in denominator
    return static_cast<double>(concordant - discordant) /
           static_cast<double>(denom);
}

std::vector<std::vector<double>> EmpiricalCopula::pseudoObservationPairs() const
{
    std::vector<std::vector<double>> out;
    if (!valid_) return out;

    out.resize(u_.size(), std::vector<double>(2, 0.0));
    for (size_t i = 0; i < u_.size(); ++i)
    {
        out[i][0] = u_[i];
        out[i][1] = v_[i];
    }
    return out;
}

bool EmpiricalCopula::writePseudoObservations(const std::string& filename) const
{
    if (!valid_) return false;

    std::ofstream f(filename.c_str());
    if (!f.good()) return false;

    f << "x,y,u,v\n";
    for (size_t i = 0; i < u_.size(); ++i)
        f << x_[i] << "," << y_[i] << "," << u_[i] << "," << v_[i] << "\n";

    return true;
}

bool EmpiricalCopula::writeCDFGrid(const std::string& filename, int nBins) const
{
    if (!valid_) return false;
    if (nBins <= 0) return false;

    std::ofstream f(filename.c_str());
    if (!f.good()) return false;

    f << "u,v,C\n";

    for (int i = 0; i < nBins; ++i)
    {
        const double uVal = (i + 0.5) / static_cast<double>(nBins);
        for (int j = 0; j < nBins; ++j)
        {
            const double vVal = (j + 0.5) / static_cast<double>(nBins);
            f << uVal << "," << vVal << "," << cdf(uVal, vVal) << "\n";
        }
    }

    return true;
}
