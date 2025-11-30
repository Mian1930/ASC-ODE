#pragma once
// Single canonical global aliases for the project's linear algebra types.
#include "../nanoblas/src/vector.hpp"
#include "../nanoblas/src/matrix.hpp"

template<typename T = double>
using Vector = nanoblas::Vector<T>;

template<typename T = double>
using VectorView = nanoblas::VectorView<T>;

template<typename T = double>
using Matrix = nanoblas::Matrix<T>;

template<typename T = double>
using MatrixView = nanoblas::MatrixView<T>;
