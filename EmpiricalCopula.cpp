#include "EmpiricalCopula.h"
#include "Matrix.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <numeric>

// ============================================================================
// Internal helpers (not exposed)
// ============================================================================
namespace
{
    // Helper struct used for sorting values while keeping original indices
    struct RankedValue
    {
        double value;
        int index;
    };

    // Check if value is finite (not NaN/Inf)
    inline bool isFiniteNumber(double x)
    {
        return std::isfinite(x);
    }
}

// ============================================================================
// Constructor / lifecycle
// ============================================================================

EmpiricalCopula::EmpiricalCopula()
    : valid_(false) // Initially not fitted
{
}

// Clear all stored data
void EmpiricalCopula::clear()
{
    x_.clear();
    y_.clear();
    u_.clear();
    v_.clear();
    valid_ = false;
}

// ============================================================================
// Fit copula from raw paired data
// ============================================================================

bool EmpiricalCopula::fit(const std::vector<double>& x,
                          const std::vector<double>& y,
                          TieMethod tieMethod)
{
    clear();

    // Must be same size and non-empty
    if (x.size() != y.size()) return false;
    if (x.empty()) return false;

    // ------------------------------------------------------------------------
    // Filter only finite paired values
    // ------------------------------------------------------------------------
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

    // ------------------------------------------------------------------------
    // Convert to pseudo-observations in (0,1)
    // ------------------------------------------------------------------------
    u_ = makePseudoObservations(x_, tieMethod);
    v_ = makePseudoObservations(y_, tieMethod);

    // Sanity check
    if (u_.size() != x_.size() || v_.size() != y_.size())
    {
        clear();
        return false;
    }

    valid_ = true;
    return true;
}

// ============================================================================
// Basic accessors
// ============================================================================

bool EmpiricalCopula::isValid() const
{
    return valid_;
}

int EmpiricalCopula::size() const
{
    return static_cast<int>(u_.size());
}

const std::vector<double>& EmpiricalCopula::x() const { return x_; }
const std::vector<double>& EmpiricalCopula::y() const { return y_; }
const std::vector<double>& EmpiricalCopula::u() const { return u_; }
const std::vector<double>& EmpiricalCopula::v() const { return v_; }

// ============================================================================
// Utility: clamp value to [0,1]
// ============================================================================

double EmpiricalCopula::clampUnit(double z)
{
    if (z < 0.0) return 0.0;
    if (z > 1.0) return 1.0;
    return z;
}

// ============================================================================
// Rank computation (core step of empirical copula)
// ============================================================================

std::vector<double> EmpiricalCopula::computeRanks(const std::vector<double>& values,
                                                  TieMethod tieMethod)
{
    const int n = static_cast<int>(values.size());
    std::vector<double> ranks(n, 0.0);
    if (n == 0) return ranks;

    // Pair values with indices
    std::vector<RankedValue> arr(n);
    for (int i = 0; i < n; ++i)
    {
        arr[i].value = values[i];
        arr[i].index = i;
    }

    // Sort by value (stable via index)
    std::sort(arr.begin(), arr.end(),
              [](const RankedValue& a, const RankedValue& b)
              {
                  if (a.value < b.value) return true;
                  if (a.value > b.value) return false;
                  return a.index < b.index;
              });

    // Assign ranks (handle ties)
    int i = 0;
    while (i < n)
    {
        int j = i;

        // Find tie block
        while (j + 1 < n && arr[j + 1].value == arr[i].value)
            ++j;

        if (tieMethod == AverageRank)
        {
            // Average rank across tied values
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

// ============================================================================
// Convert ranks → pseudo-observations in (0,1)
// ============================================================================

std::vector<double> EmpiricalCopula::makePseudoObservations(
    const std::vector<double>& values,
    TieMethod tieMethod)
{
    const int n = static_cast<int>(values.size());
    std::vector<double> out(n, 0.0);
    if (n == 0) return out;

    const std::vector<double> ranks = computeRanks(values, tieMethod);

    // Use (r - 0.5) / n to avoid exactly 0 and 1
    for (int i = 0; i < n; ++i)
        out[i] = (ranks[i] - 0.5) / static_cast<double>(n);

    return out;
}

// ============================================================================
// Empirical copula CDF: C(u,v) = P(U ≤ u, V ≤ v)
// ============================================================================

double EmpiricalCopula::cdf(double uVal, double vVal) const
{
    if (!valid_ || u_.empty()) return 0.0;

    uVal = clampUnit(uVal);
    vVal = clampUnit(vVal);

    double count = 0.0;

    // Count samples in lower-left quadrant
    for (size_t i = 0; i < u_.size(); ++i)
    {
        if (u_[i] <= uVal && v_[i] <= vVal)
            count += 1.0;
    }

    return count / static_cast<double>(u_.size());
}

// ============================================================================
// Upper-tail copula: P(U > u, V > v)
// ============================================================================

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

// ============================================================================
// Build binned copula (2D histogram in unit square)
// ============================================================================

CopulaBinnedMatrix EmpiricalCopula::makeBinnedMatrix(int nBins, bool normalize) const
{
    CopulaBinnedMatrix M(nBins);

    if (!valid_ || nBins <= 0) return M;

    // Accumulate all (u,v) pairs
    M.accumulateUnitPairs(u_, v_);

    if (normalize)
        M.normalizeToUnitMass(); // Convert counts → probabilities

    return M;
}

// ============================================================================
// Rank correlation (Spearman rho)
// ============================================================================

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

// ============================================================================
// Kendall tau (O(N^2))
// ============================================================================

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

    return static_cast<double>(concordant - discordant) /
           static_cast<double>(denom);
}

// ============================================================================
// Return all (u,v) pairs
// ============================================================================

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

// ============================================================================
// Write raw + copula data to CSV
// ============================================================================

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

// ============================================================================
// Write empirical CDF grid
// ============================================================================

bool EmpiricalCopula::writeCDFGrid(const std::string& filename, int nBins) const
{
    if (!valid_) return false;
    if (nBins <= 0) return false;

    std::ofstream f(filename.c_str());
    if (!f.good()) return false;

    f << "u,v,C\n";

    // Evaluate C(u,v) on grid centers
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
