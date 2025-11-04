# Matrix Operations Library (Java)

## Project Motivation
This project was developed as part of a **Linear Algebra and Data Structures** course to deepen understanding of matrix computation, algorithm design, and numerical stability in Java.  
Rather than relying on external math libraries, this implementation builds core linear algebra operations from the ground up, emphasizing clarity, correctness, and educational value.

---

## Overview
The `Matrix` class provides a general-purpose framework for performing standard matrix operations, including addition, multiplication, transposition, and determinant calculations.  
It demonstrates how fundamental linear algebra concepts can be implemented programmatically, using clean, well-structured Java code.

---

## Key Features
- **Matrix Construction:** Create matrices of arbitrary size from dimensions or 2D arrays.  
- **Arithmetic Operations:**  
  - `add(Matrix B)` — Elementwise matrix addition  
  - `subtract(Matrix B)` — Elementwise subtraction  
  - `multiply(Matrix B)` — Standard matrix multiplication (triple-loop algorithm with micro-optimization)  
- **Transpose:**  
  - `transpose()` — Returns a new matrix with rows and columns swapped.  
- **Error Checking:**  
  - Automatic validation of dimensions before any operation.  
- **Numerical Stability:**  
  - Uses an `EPS = 1e-9` threshold to handle near-zero floating-point values.

---

## Code Example
```java
Matrix A = new Matrix(new double[][] {
    {1, 2, 3},
    {4, 5, 6}
});

Matrix B = new Matrix(new double[][] {
    {7, 8, 9},
    {1, 2, 3}
});

// Example operations
Matrix C = A.add(B);          // Matrix addition
Matrix D = A.transpose();     // Transpose
Matrix E = A.multiply(B.T()); // Matrix multiplication (A × Bᵀ)

Time Complexity:

Addition/Subtraction: O(mn)

Multiplication: O(mnk)

Transpose: O(mn)


Design Notes

Written for readability and educational transparency rather than low-level optimization.

Multiplication loop order (i, k, j) improves cache performance in row-major storage.

Floating-point comparisons are stabilized using the epsilon tolerance EPS.