#pragma once

#include <vector>
#include <string>

class CopulaBinnedMatrix;

class EmpiricalCopula
{
public:
    enum TieMethod
    {
        AverageRank,
        MinRank
    };

    EmpiricalCopula();

    bool fit(const std::vector<double>& x,
             const std::vector<double>& y,
             TieMethod tieMethod = AverageRank);

    void clear();

    bool isValid() const;
    int size() const;

    const std::vector<double>& x() const;
    const std::vector<double>& y() const;
    const std::vector<double>& u() const;
    const std::vector<double>& v() const;

    // Empirical copula C_n(u,v) = (1/n) sum I(U_i <= u, V_i <= v)
    double cdf(double uVal, double vVal) const;

    // Survival copula-like empirical tail count:
    // (1/n) sum I(U_i > u, V_i > v)
    double upperTailCDF(double uVal, double vVal) const;

    // Build a binned representation on the unit square.
    CopulaBinnedMatrix makeBinnedMatrix(int nBins, bool normalize = true) const;

    // Rank-based dependence diagnostics from pseudo-observations.
    double spearmanRho() const;
    double kendallTau() const;

    // Useful diagnostics / exports
    std::vector<std::vector<double>> pseudoObservationPairs() const;
    bool writePseudoObservations(const std::string& filename) const;
    bool writeCDFGrid(const std::string& filename, int nBins) const;

    static std::vector<double> makePseudoObservations(const std::vector<double>& values,
                                                      TieMethod tieMethod = AverageRank);

private:
    std::vector<double> x_;
    std::vector<double> y_;
    std::vector<double> u_;
    std::vector<double> v_;
    bool valid_;

    static double clampUnit(double z);
    static std::vector<double> computeRanks(const std::vector<double>& values,
                                            TieMethod tieMethod);
};
