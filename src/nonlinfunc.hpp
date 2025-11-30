#pragma once
// NonlinearFunction interface using the canonical global Vector/VectorView.
#include "vector.hpp"   // provides Vector, VectorView, MatrixView

class NonlinearFunction
{
public:
    virtual ~NonlinearFunction() = default;

    // dimension of input
    virtual std::size_t dimX() const = 0;
    // dimension of output
    virtual std::size_t dimF() const = 0;

    // evaluate f(x)
    virtual void evaluate(VectorView<double> x, VectorView<double> f) const = 0;

    // evaluate derivative df/dx
    virtual void evaluateDeriv(VectorView<double> x, MatrixView<double> df) const = 0;
};
