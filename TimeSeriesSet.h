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

 // TimeSeriesSet.h

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
class TimeSeriesSet : public std::vector<TimeSeries<T>> {
public:
    // Constructors
    TimeSeriesSet(); ///< Default constructor
    TimeSeriesSet(const TimeSeriesSet<T>& other); ///< Copy constructor
    TimeSeriesSet(const TimeSeries<T>& ts); ///< Construct with a single TimeSeries
    TimeSeriesSet(int num_series); ///< Construct with specified number of empty series
    TimeSeriesSet(const std::vector<std::string>& filenames); ///< Construct from multiple file inputs
    TimeSeriesSet(const std::string& filename, bool has_header = true); ///< Construct from CSV file

    TimeSeriesSet<T>& operator=(const TimeSeriesSet<T>& other); ///< Copy assignment
    TimeSeriesSet<T>& operator=(TimeSeriesSet<T>&& other) noexcept; ///< Move assignment

    // File I/O
    bool read(const std::string& filename, bool has_header = true); ///< Read from CSV
    void write(const std::string& filename, const std::string& delimiter = ",") const; ///< Write to CSV
    void appendtofile(const std::string& filename, bool include_time = true) const; ///< Append to existing CSV

    // Accessors
    TimeSeries<T>& operator[](int index); ///< Access by index
    TimeSeries<T>& operator[](const std::string& name); ///< Access by name (mutable)
    const TimeSeries<T>& operator[](int index) const; ///< Access by index (const)
    const TimeSeries<T>& operator[](const std::string& name) const; ///< Access by name (const)

    // Metadata
    void setSeriesName(int index, const std::string& name); ///< Set series name by index
    std::string getSeriesName(int index) const; ///< Get series name by index
    std::vector<std::string> getSeriesNames() const; ///< Get all series names
    void SetAllSeriesSize(int n); ///<< Resizes all the timeseries
    bool SetRow(int i, const T& t, const std::vector<T>& c); /// Set values of an entire row
    void setname(int index, const std::string& name); ///< Set name on internal series object

    //Ranges
    T maxval() const;
    T minval() const;
    T mintime() const;
    T maxtime() const;

    // Query
    bool Contains(const std::string& name) const; ///< Check if series exists by name
    int indexOf(const std::string& name) const; ///< Get index of named series
    int lookup(const std::string& name) const; ///< Alias for indexOf
    int maxnumpoints() const; ///< Get max number of points across all series
    size_t seriesCount() const; ///< Number of series in the set
    size_t pointCount(int series_index) const; ///< Number of points in a specific series

    // Append & Modify
    void append(const std::vector<T>& values); ///< Append a row of values
    void append(double t, const std::vector<T>& values); ///< Append a row with time
    bool append(const TimeSeries<T>& ts, const std::string& name = ""); ///< Append a named TimeSeries
    void append(const std::string& name); ///< Add empty TimeSeries with given name
    void clear(); ///< Clear all series
    void clearContents(); ///< Clear content of each series, retain structure
    void clearContentsExceptLastRow(); ///< Keep last row only
    void knockout(T t); ///< Knock out all points beyond time t
    void resize(size_t num_series); ///< Resize to given number of series
    //void ResizeIfNeeded(size_t new_size); ///< Expand size only if needed
    void removeNaNs(); ///< Remove NaN values from all series
    TimeSeriesSet<T> removeNaNs() const; ///< Return a new set with NaN values removed

    // Interpolation & Extraction
    std::vector<T> interpolate(T t) const; ///< Interpolate all series at given time
    std::vector<T> interpolate(T t, int n) const; ///< Interpolate for first n series
    TimeSeries<T> extract(int index, T t1, T t2) const; ///< Extract a subset of one series

    TimeSeriesSet<T> derivative() const;

    TimeSeriesSet<T> extract(T t1, T t2) const;

    // Series Statistics
    std::vector<T> getrandom() const; ///< Random row of values
    std::vector<T> getrandom(int start_item) const; ///< Random row after index
    std::vector<T> getrow(int index) const; ///< Get values from specific row
    std::vector<T> percentile(T p) const; ///< Percentiles from each series
    std::vector<T> percentile(T p, int start_item) const;
    std::vector<T> percentile(T p, int start_item, const std::vector<int>& indices) const;
    std::vector<T> mean(int start_item = 0) const; ///< Means from each series
    std::vector<T> mean(int start_item, const std::vector<int>& indices) const;
    std::vector<T> standardDeviation(int start_item = 0) const; ///< Stddev from each series
    std::vector<T> standardDeviation(int start_item, const std::vector<int>& indices) const;
    std::vector<T> min(int start_item = 0) const; ///< Minimums from each series
    std::vector<T> max(int start_item = 0) const; ///< Maximums from each series
    std::vector<T> integrate() const; ///< Trapezoidal integration
    std::vector<T> average() const; ///< Time-averaged values

    // Correlation & Linear Combinations
    CMatrix correlation(int start_item, int end_item) const; ///< Correlation matrix
    TimeSeries<T> add(const std::vector<int>& indices) const; ///< Sum multiple series
    TimeSeries<T> add_mult(const std::vector<int>& indices, const std::vector<T>& weights) const; ///< Weighted sum
    TimeSeries<T> add_mult(const std::vector<int>& indices, const TimeSeriesSet<T>& other) const; ///< Elementwise multiplication with another set
    TimeSeries<T> divide(int numerator_index, int denominator_index) const; ///< Pointwise division

    // Transformations
    TimeSeriesSet<T> make_uniform(T increment, bool assign_d = true) const; ///< Make all series uniform
    TimeSeriesSet<T> getpercentiles(const std::vector<T>& fractions) const; ///< Get specified percentiles at all times
    TimeSeriesSet<T> distribution(int n_bins, int start_index=0, int end_index=-1) const; ///< Density histograms
    TimeSeriesSet<T> distributionLog(int n_bins, int start_index, int end_index=-1) const; ///< Density histograms with logarithmic intervals
    TimeSeriesSet<T> add_noise(const std::vector<T>& stddevs, bool log_noise = false) const; ///< Add Gaussian noise
    TimeSeriesSet<T> sort(int column_index = 0) const; ///< Sort based on a column
    TimeSeriesSet<T> ConverttoNormalScore() const; ///< Convert each to normal score
    TimeSeriesSet<T> AutoCorrelation(const double& span, const double& increment) const; ///< Autocorrelation function
    TimeSeriesSet<T> cummulative() const; ///< convert to cummulative distribution
    TimeSeriesSet<T> GetCummulativeDistribution() const; ///< Extract CDF for each series
    TimeSeriesSet<T> Log() const; ///< Log-transform values

    // Wiggle Analysis
    std::vector<T> max_wiggle() const; ///< Max wiggle metric across series
    std::vector<T> max_wiggle_corr(int back_steps) const; ///< Max wiggle correlation
    std::vector<int> max_wiggle_sl(int back_steps, T tolerance) const; ///< Wiggle slope test

    // File Serialization (Qt)
#ifdef Q_JSON_SUPPORT
    QJsonObject toJson() const; ///< Convert to QJsonObject
    void fromJson(const QJsonObject& json); ///< Load from QJsonObject
#endif // Q_JSON_SUPPORT

    // ============================================================
    // TimeSeriesSet mean as a single TimeSeries
    //
    // mean_ts(...)       : FAST / legacy behavior
    //   - assumes aligned time grids (same timestamps by index)
    //
    // mean_ts_union(...) : ROBUST (but can be slower)
    //   - uses UNION of timestamps from series (>= start_item)
    //   - at each union time, averages interpolated values from series
    //     that cover that time (available-only behavior)
    //
    // mean_ts_grid(...)  : FAST + ROBUST (recommended)
    //   - builds a UNIFORM time grid with dt you choose (e.g., 0.01)
    //   - supports:
    //       Intersection   -> only where ALL series cover the time
    //       UnionAvailable -> extend to longest end; average available series
    //   - streaming interpolation (no ts.interpol() inside inner loops)
    // ============================================================

    TimeSeries<T> mean_ts(int start_item = 0) const;
    TimeSeries<T> mean_ts(int start_item, const std::vector<int>& indices) const;

    TimeSeries<T> mean_ts_union(int start_item = 0, T time_merge_tol = (T)1e-12) const;
    TimeSeries<T> mean_ts_union(int start_item,
                                const std::vector<int>& indices,
                                T time_merge_tol = (T)1e-12) const;

    // ✅ NEW: uniform-grid mean (fast)
    TimeSeries<T> mean_ts_grid(T dt,
                              int start_item = 0,
                              MeanGridMode mode = MeanGridMode::Intersection,
                              T time_eps = (T)1e-12) const;

    TimeSeries<T> mean_ts_grid(T dt,
                              int start_item,
                              const std::vector<int>& indices,
                              MeanGridMode mode = MeanGridMode::Intersection,
                              T time_eps = (T)1e-12) const;

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

    // Properties
    std::string filename; ///< File associated with the set
    std::string name; ///< Name of the set
    bool file_not_found = false; ///< Flag for file read failure
    bool unif = false; ///< Whether time steps are uniform

private:
    static size_t countRows(std::ifstream& file, bool has_header); // estimate the number of rows in a csv file
};

// Helper functions

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

#include "TimeSeriesSet.hpp"
