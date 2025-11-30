// include/timestepper.hpp
#pragma once
#include <memory>
#include <vector>
#include <iostream>
#include <cmath>

// include the project's canonical NonlinearFunction and Vector types
#include "../src/nonlinfunc.hpp"

// Base TimeStepper used by the steppers in include/
class TimeStepper
{
protected:
  std::shared_ptr<NonlinearFunction> m_rhs;
public:
  TimeStepper(std::shared_ptr<NonlinearFunction> rhs) : m_rhs(rhs) {}
  virtual ~TimeStepper() = default;
  virtual void doStep(double tau, VectorView<double> y) = 0;
};
