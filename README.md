# Utilities

This repository contains reusable C++ utility classes used by environmental
and numerical modeling components (vector/matrix math, probability
distributions, and time-series data structures).

## Repository map

- `Vector.h` / `Vector.cpp`: Dynamic numeric vector class (`CVector`) and
  arithmetic helpers.
- `Matrix.h` / `Matrix.cpp`: Matrix class (`CMatrix`) and linear algebra
  helpers.
- `Distribution.h` / `Distribution.cpp`: Probability distribution wrapper and
  random sampling utilities.
- `TimeSeries.h` / `TimeSeries.hpp`: Generic `TimeSeries<T>` implementation.
- `TimeSeriesSet.h` / `TimeSeriesSet.hpp`: Collections of time series and
  ensemble statistics.
- `BTC.h` / `BTC.hpp`: Legacy-compatible `CTimeSeries<T>` class template.
- `Utilities.h` / `Utilities.cpp`: Generic string/file/math helper functions.

## Commenting approach used in this codebase

To make maintenance easier, comments should explain:

1. **Intent** (why a type/function exists),
2. **Contracts** (input assumptions + output guarantees),
3. **Edge cases** (empty vectors, out-of-range indices, NaNs),
4. **Interoperability** (optional Armadillo/Qt/GSL build paths).

The headers now include section-level comments so contributors can quickly
navigate API groups without scanning implementation files.
