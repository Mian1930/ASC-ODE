#pragma once
// Lightweight wrapper typedefs that expose nanoblas matrix types in ASC_ode.
#include "../nanoblas/src/matrix.hpp"

namespace ASC_ode {

template<typename T = double>
using Matrix = nanoblas::Matrix<T>;

template<typename T = double>
using MatrixView = nanoblas::MatrixView<T>;

} // namespace ASC_ode
