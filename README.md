# Multithreaded C++ Polynomial Library

A performant C++ library for polynomial arithmetic. It supports addition, subtraction, multiplication (standard & FFT-accelerated), and modulo operations, all with optional multithreading for large inputs.

## Features

- Addition, subtraction, multiplication, modulo
- Thread-safe data structures with `std::mutex`
- FFT acceleration for dense polynomials (degree > 1000)
- Canonical form output and simple printing
- Sparse-friendly standard multiplication

## Usage

Include the header and link the implementation:

```cpp
#include "poly.h"
// ... build your polynomials ...
polynomial p1{{1,1},{0,1}}; // x + 1
polynomial p2{{1,2},{0,3}}; // 2x + 3
auto sum = p1 + p2;         // 3x + 4
auto prod = p1 * p2;        // 2x^2 + 5x + 3
```

Use `canonical_form()` to get a sorted vector of terms, or `print()` for quick debug.

## Testing

- `simple_poly.txt` and `result.txt` drive the sample tests in `main.cpp`.  
- Modify those files to add new test cases; each line is `coeffx^power` separated by newlines, groups split by `;`.
