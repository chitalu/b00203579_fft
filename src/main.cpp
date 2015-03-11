#define _USE_MATH_DEFINES
#include <cmath>

#include <cstdio>
#include <cstdlib>
#include <vector>

// fft in the west!
#include "fftw3.h"

#define REAL_PART 0
#define IMAGINARY_PART 1

#define Imag(x_) x_[IMAGINARY_PART]
#define Real(x_) x_[REAL_PART]

#define N (1024)

const double two_pi = M_PI * 2.0;
const int freq = 440;
const double t = 0.0001;

const double w = two_pi * freq;

int main(int argc, char const *argv[]) {
  printf("%s\n", "hello fft!");

  // time domain decomposed signal impulses
  fftw_complex *x = (fftw_complex *)malloc(sizeof(fftw_complex) * N);
  // frequency domain amplitudes
  fftw_complex *X = (fftw_complex *)malloc(sizeof(fftw_complex) * N);
  ;

  for (int n = 0; n < N; ++n) {
    Real(x[n]) = sin(w * (double)n * t);
    Imag(x[n]) = 0;

    Real(X[n]) = 0;
    Imag(X[n]) = 0;
  }

  fftw_plan forward_fft_plan, inverse_fft_plan;

  forward_fft_plan = fftw_plan_dft_1d(N, x, X, FFTW_FORWARD, FFTW_ESTIMATE);

  fftw_execute(forward_fft_plan);

  inverse_fft_plan = fftw_plan_dft_1d(N, X, x, FFTW_BACKWARD, FFTW_ESTIMATE);

  fftw_execute(forward_fft_plan);

  fftw_destroy_plan(forward_fft_plan);
  fftw_destroy_plan(inverse_fft_plan);

  return 0;
}