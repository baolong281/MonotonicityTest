#ifndef MATRIX_H
#define MATRIX_H

#include <Rcpp.h>

namespace matrix_utils {

// Rcpp doesnt have a built in operator for matrix multiplication, so we have to
// make our own
Rcpp::NumericMatrix operator*(const Rcpp::NumericMatrix &a,
                              const Rcpp::NumericMatrix &b) {
  int nrow = a.nrow();
  int ncol = b.ncol();
  int ninner = a.ncol();

  if (ninner != b.nrow()) {

    Rcpp::stop("Incompatible matrix dimensions for multiplication.");
  }

  Rcpp::NumericMatrix result(nrow, ncol);
  // Normal naive matrix multiplication
  // Our matrices are 2x2 max so this is fine
  for (int i = 0; i < nrow; ++i) {
    for (int j = 0; j < ncol; ++j) {
      double sum = 0.0;
      for (int k = 0; k < ninner; ++k) {
        sum += a(i, k) * b(k, j);
      }
      result(i, j) = sum;
    }
  }
  return result;
}

// overload the plus operator for matrix addition
Rcpp::NumericMatrix operator+(const Rcpp::NumericMatrix &a,
                              const Rcpp::NumericMatrix &b) {
  if (a.nrow() != b.nrow() || a.ncol() != b.ncol()) {
    Rcpp::stop("Incompatible matrix dimensions for addition.");
  }
  Rcpp::NumericMatrix result(clone(a));
  int nrow = result.nrow();
  int ncol = result.ncol();
  for (int i = 0; i < nrow; ++i) {
    for (int j = 0; j < ncol; ++j) {
      result(i, j) += b(i, j);
    }
  }
  result.attr("dim") = Rcpp::Dimension(a.nrow(), a.ncol());
  return result;
}

// overload the minus operator for matrix subtraction
Rcpp::NumericMatrix operator-(const Rcpp::NumericMatrix &a,
                              const Rcpp::NumericMatrix &b) {
  return a + (-1.0 * b);
}

// Overload the * operator for scalar multiplication
Rcpp::NumericMatrix operator*(const Rcpp::NumericMatrix &a,
                              const double &scalar) {
  Rcpp::NumericMatrix result(clone(a));
  int nrow = result.nrow();
  int ncol = result.ncol();
  for (int i = 0; i < nrow; ++i) {
    for (int j = 0; j < ncol; ++j) {
      result(i, j) *= scalar;
    }
  }
  return result;
}

// Define it the other way around
Rcpp::NumericMatrix operator*(const double &scalar,
                              const Rcpp::NumericMatrix &a) {
  return a * scalar;
}

Rcpp::NumericMatrix operator/(const Rcpp::NumericMatrix &a,
                              const double &scalar) {
  return a * (1.0 / scalar);
}

Rcpp::NumericMatrix create_colvec(Rcpp::NumericVector v) {
  if (v.length() != 2) {
    Rcpp::stop("Input vector must be a single element");
  }

  v.attr("dim") = Rcpp::Dimension(2, 1);
  return Rcpp::as<Rcpp::NumericMatrix>(v);
}

} // namespace matrix_utils

#endif
