#define _USE_MATH_DEFINES
#include <cmath>

#include <cstdio>
#include <cstdlib>
#include <vector>

// fft in the west!
#include "fftw3.h"

#include "base.h"

#define N (1024)

const int freq = 440;
const double t = 0.0001;
const double w = two_pi * freq;

int main(int argc, char const *argv[]) {
  printf("%s\n", "hello fft!");

  fftw_complex *x = NULL, *X = NULL;

  START_PROFILING(allocate_mem) {
    x = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);
    X = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);
  }
  STOP_PROFILING()

  START_PROFILING(initialise_mem) {
    for (int n = 0; n < N; ++n) {
	  // time domain decomposed signal impulses
      Re(x[n]) = sin(w * (double)n * t);
      Im(x[n]) = 0;

	  // frequency domain amplitudes
      Re(X[n]) = 0;
      Im(X[n]) = 0;
    }
  }
  STOP_PROFILING();

  fftw_plan forward_fft_plan, inverse_fft_plan;

  START_PROFILING(complex_fft) {
    forward_fft_plan = fftw_plan_dft_1d(N, x, X, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(forward_fft_plan);
  }
  STOP_PROFILING();

  START_PROFILING(complex_ifft) {
    inverse_fft_plan = fftw_plan_dft_1d(N, X, x, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(forward_fft_plan);
  }
  STOP_PROFILING();

  START_PROFILING(deallocate_mem) {
    fftw_destroy_plan(forward_fft_plan);
    fftw_destroy_plan(inverse_fft_plan);
    fftw_free(x);
    x = NULL;
    fftw_free(X);
    X = NULL;
  }
  STOP_PROFILING();

  for (ptime_iter_t i = g_tstamps.cbegin(); i != g_tstamps.cend(); ++i) {
    printf("task: %s\ntime: %f microseconds\n\n", i->first.c_str(), i->second);
  }
  return 0;
}