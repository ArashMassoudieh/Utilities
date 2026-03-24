/*
 * OpenHydroQual - Environmental Modeling Platform
 * Copyright (C) 2025 Arash Massoudieh
 *
 * This file is part of OpenHydroQual.
 *
 * OpenHydroQual is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * If you use this file in a commercial product, you must purchase a
 * commercial license. Contact arash.massoudieh@enviroinformatics.co for details.
 */

#pragma once

#include <vector>
#include <string>
#include <optional>
#include <fstream>

#include "TimeSeries.h"

/**
 * @brief Mean grid mode for mean_ts_grid()
 * - Intersection   : only where ALL selected series cover the time (common overlap)
 * - UnionAvailable : from min start to max end; averages only available series at each time
 */
enum class MeanGridMode
{
    Intersection,
    UnionAvailable
};

/**
 * @brief A class that represents a collection of TimeSeries objects.
 */
template<typename T>
class TimeSeriesSet : public std::vector<TimeSeries<T>>
{
public:

    // ============================================================
    // Constructors
    // ============================================================

    TimeSeriesSet();
    TimeSeriesSet(const TimeSeriesSet<T>& other);
    TimeSeriesSet(const TimeSeries<T>& ts);
    TimeSeriesSet(int num_series);
    TimeSeriesSet(const std::vector<std::string>& filenames);
    TimeSeriesSet(const std::string& filename, bool has_header = true);

    TimeSeriesSet<T>& operator=(const TimeSeriesSet<T>& other);
    TimeSeriesSet<T>& operator=(TimeSeriesSet<T>&& other) noexcept;

    // ============================================================
    // File I/O
    // ============================================================

    bool read(const std::string& filename, bool has_header = true);
    void write(const std::string& filename, const std::string& delimiter = ",") const;
    void appendtofile(const std::string& filename, bool include_time = true) const;

    // ============================================================
    // Accessors
    // ============================================================

    TimeSeries<T>& operator[](int index);
    TimeSeries<T>& operator[](const std::string& name);
    const TimeSeries<T>& operator[](int index) const;
    const TimeSeries<T>& operator[](const std::string& name) const;

    // ============================================================
    // Metadata
    // ============================================================

    void setSeriesName(int index, const std::string& name);
    std::string getSeriesName(int index) const;
    std::vector<std::string> getSeriesNames() const;

    void SetAllSeriesSize(int n);
    bool SetRow(int i, const T& t, const std::vector<T>& c);
    void setname(int index, const std::string& name);

    // ============================================================
    // Ranges
    // ============================================================

    T maxval() const;
    T minval() const;
    T mintime() const;
    T maxtime() const;

    // ============================================================
    // Query
    // ============================================================

    bool Contains(const std::string& name) const;
    int indexOf(const std::string& name) const;
    int lookup(const std::string& name) const;
    int maxnumpoints() const;
    size_t seriesCount() const;
    size_t pointCount(int series_index) const;

    // ============================================================
    // Append & Modify
    // ============================================================

    void append(const std::vector<T>& values);
    void append(double t, const std::vector<T>& values);
    bool append(const TimeSeries<T>& ts, const std::string& name = "");
    void append(const std::string& name);

    void clear();
    void clearContents();
    void clearContentsExceptLastRow();
    void knockout(T t);
    void resize(size_t num_series);
    void removeNaNs();
    TimeSeriesSet<T> removeNaNs() const;

    // ============================================================
    // Interpolation & Extraction
    // ============================================================

    std::vector<T> interpolate(T t) const;
    std::vector<T> interpolate(T t, int n) const;

    TimeSeries<T> extract(int index, T t1, T t2) const;
    TimeSeriesSet<T> extract(T t1, T t2) const;

    TimeSeriesSet<T> derivative() const;

    // ============================================================
    // Series Statistics
    // ============================================================

    std::vector<T> getrandom() const;
    std::vector<T> getrandom(int start_item) const;
    std::vector<T> getrow(int index) const;

    std::vector<T> percentile(T p) const;
    std::vector<T> percentile(T p, int start_item) const;
    std::vector<T> percentile(T p, int start_item, const std::vector<int>& indices) const;

    std::vector<T> mean(int start_item = 0) const;
    std::vector<T> mean(int start_item, const std::vector<int>& indices) const;

    std::vector<T> standardDeviation(int start_item = 0) const;
    std::vector<T> standardDeviation(int start_item, const std::vector<int>& indices) const;

    std::vector<T> min(int start_item = 0) const;
    std::vector<T> max(int start_item = 0) const;

    std::vector<T> integrate() const;
    std::vector<T> average() const;

    // ============================================================
    // Correlation & Linear Combinations
    // ============================================================

    CMatrix correlation(int start_item, int end_item) const;

    TimeSeries<T> add(const std::vector<int>& indices) const;
    TimeSeries<T> add_mult(const std::vector<int>& indices, const std::vector<T>& weights) const;
    TimeSeries<T> add_mult(const std::vector<int>& indices, const TimeSeriesSet<T>& other) const;
    TimeSeries<T> divide(int numerator_index, int denominator_index) const;

    // ============================================================
    // Transformations
    // ============================================================

    TimeSeriesSet<T> make_uniform(T increment, bool assign_d = true) const;
    TimeSeriesSet<T> getpercentiles(const std::vector<T>& fractions) const;
    TimeSeriesSet<T> distribution(int n_bins, int start_index = 0, int end_index = -1) const;
    TimeSeriesSet<T> distributionLog(int n_bins, int start_index, int end_index = -1) const;
    TimeSeriesSet<T> add_noise(const std::vector<T>& stddevs, bool log_noise = false) const;
    TimeSeriesSet<T> sort(int column_index = 0) const;
    TimeSeriesSet<T> ConverttoNormalScore() const;
    TimeSeriesSet<T> AutoCorrelation(const double& span, const double& increment) const;
    TimeSeriesSet<T> cummulative() const;
    TimeSeriesSet<T> GetCummulativeDistribution() const;
    TimeSeriesSet<T> Log() const;

    // ============================================================
    // Wiggle Analysis
    // ============================================================

    std::vector<T> max_wiggle() const;
    std::vector<T> max_wiggle_corr(int back_steps) const;
    std::vector<int> max_wiggle_sl(int back_steps, T tolerance) const;

#ifdef Q_JSON_SUPPORT
    QJsonObject toJson() const;
    void fromJson(const QJsonObject& json);
#endif

    // ============================================================
    // TimeSeriesSet mean as a single TimeSeries
    //
    // mean_ts(...)        : DEFAULT (ROBUST)
    //   - uses LONGEST series time grid
    //   - averages available series at each time
    //   - NO extrapolation
    //
    // mean_ts_union(...)  : ROBUST
    //   - UNION of timestamps (>= start_item)
    //   - averages interpolated values from series that cover time
    //   - NO extrapolation
    //
    // mean_ts_grid(...)   : FAST + ROBUST
    //   - uniform grid with dt
    //   - Intersection / UnionAvailable modes
    //   - streaming interpolation
    //   - NO extrapolation
    //
    // mean_ts_longest(...) : ROBUST + SIMPLE
    //   - explicit longest-series-grid version
    // ============================================================

    TimeSeries<T> mean_ts(int start_item = 0) const;
    TimeSeries<T> mean_ts(int start_item, const std::vector<int>& indices) const;

    TimeSeries<T> mean_ts_union(int start_item = 0,
                                T time_merge_tol = (T)1e-12) const;

    TimeSeries<T> mean_ts_union(int start_item,
                                const std::vector<int>& indices,
                                T time_merge_tol = (T)1e-12) const;

    TimeSeries<T> mean_ts_grid(T dt,
                               int start_item = 0,
                               MeanGridMode mode = MeanGridMode::Intersection,
                               T time_eps = (T)1e-12) const;

    TimeSeries<T> mean_ts_grid(T dt,
                               int start_item,
                               const std::vector<int>& indices,
                               MeanGridMode mode = MeanGridMode::Intersection,
                               T time_eps = (T)1e-12) const;

    TimeSeries<T> mean_ts_longest(int start_item = 0,
                                  T time_eps = (T)1e-10) const;

    TimeSeries<T> mean_ts_longest(int start_item,
                                  const std::vector<int>& indices,
                                  T time_eps = (T)1e-10) const;

#ifdef TORCH_SUPPORT
    static TimeSeriesSet<T> fromTensor(const torch::Tensor& tensor,
                                       T start_time,
                                       T end_time,
                                       const std::vector<std::string>& series_names = {});

    torch::Tensor toTensor(bool include_time = false,
                           torch::Device device = torch::kCPU) const;

    torch::Tensor toTensorAtIntervals(double t_start, double t_end, double dt,
                                      bool include_time = false,
                                      torch::Device device = torch::kCPU) const;
#endif

    // ============================================================
    // Properties
    // ============================================================

    std::string filename;
    std::string name;
    bool file_not_found = false;
    bool unif = false;

private:
    static size_t countRows(std::ifstream& file, bool has_header);
};

// ============================================================
// Helper functions
// ============================================================

template<class T>
T diff(TimeSeriesSet<T> A, TimeSeriesSet<T> B);

template<class T>
TimeSeriesSet<T> merge(TimeSeriesSet<T> A, const TimeSeriesSet<T>& B);

template<class T>
TimeSeriesSet<T> merge(std::vector<TimeSeriesSet<T>>& sets);

template<class T>
TimeSeriesSet<T> operator*(const TimeSeriesSet<T>& set, const T& scalar);

template<class T>
CVector norm2dif(TimeSeriesSet<T>& A, TimeSeriesSet<T>& B);

template<class T>
CVector sum_interpolate(std::vector<TimeSeriesSet<T>>& sets, double t);

template<class T>
T sum_interpolate(std::vector<TimeSeriesSet<T>>& sets, T t, std::string name);

template<class T>
int max_n_vars(std::vector<TimeSeriesSet<T>>& sets);

// ✅ NEW: mean over ensemble of multi-column TimeSeriesSet<T>
// Each column is averaged independently using longest-grid logic.
// Each output column may have its own t_end.
//
// IMPORTANT: avoid default-arg redeclaration errors by placing defaults
// in ONLY ONE overload (see below).

// Full overload (NO defaults here)
template<typename T>
TimeSeriesSet<T> mean_ts_longest_cols(const std::vector<TimeSeriesSet<T>>& sets,
                                      int start_item,
                                      T time_eps);

// Convenience overload (defaults live ONLY here)
template<typename T>
TimeSeriesSet<T> mean_ts_longest_cols(const std::vector<TimeSeriesSet<T>>& sets,
                                      int start_item = 0);

#include "TimeSeriesSet.hpp"
