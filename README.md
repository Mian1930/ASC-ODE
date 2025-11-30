# ASC-ODE — Quick Build & Run Guide

This repository contains ODE solvers and demos for the ASC course, including explicit/implicit Euler, Runge–Kutta, and automatic differentiation examples.

## Requirements
- g++ with C++20 support  
- No external libraries needed (nanoblas included)

## Build & Run

### 1. Mass–Spring ODE (Explicit Euler)
```bash
g++ -std=c++20 -O2 -I. -Isrc -Iinclude -Inanoblas/src \
    demos/test_ode.cpp -o build/demos/test_ode
./build/demos/test_ode
2. Runge–Kutta Demo
bash
Copy code

g++ -std=c++20 -O2 -I. -Isrc -Iinclude -Inanoblas/src \
    demos/demo_rk.cpp -o build/demos/demo_rk
./build/demos/demo_rk
3. Automatic Differentiation Demo
bash
Copy code

g++ -std=c++20 -O2 -I. -Isrc -Iinclude -Inanoblas/src \
    demos/demo_autodiff.cpp -o build/demos/demo_autodiff
./build/demos/demo_autodiff
Notes
Some demos write output files such as .txt or .csv.

All demos share the same headers in src/ and include/.

Compilation with -O2 is recommended.
