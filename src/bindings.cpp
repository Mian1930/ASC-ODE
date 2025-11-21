#include <pybind11/pybind11.h>

namespace py = pybind11;

PYBIND11_MODULE(asc_ode, m) {
    m.doc() = "ASC ODE library – winter semester 2025/26";
    m.attr("__version__") = "0.1.0";

   
    m.def("hello", []() { return "Hello from ASC ODE!"; });
}