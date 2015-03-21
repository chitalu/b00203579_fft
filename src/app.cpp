#include "base.h"

// fft in the west!
#include "fftw3.h"

std::map<std::string, double> g_tstamps;

void rfft(const std::size_t &N)
{
	const int freq = 440;
	const double t = 0.0001;
	const double w = M_2_PI * freq;

	double *x = NULL;
	fftw_complex *X = NULL;

	START_PROFILING("alloc_mem") {
		x = (double *)fftw_malloc(sizeof(double) * N);
		X = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)* N);
	}
	STOP_PROFILING();

	START_PROFILING("init_mem") {
		for (std::size_t n = 0u; n < N; ++n) {
			// time domain decomposed signal impulses
			x[n] = sin(w * (double)n * t);
			x[n] = 0;

			// frequency domain amplitudes
			Re(X[n]) = 0;
			Im(X[n]) = 0;
		}
	}
	STOP_PROFILING();

	fftw_plan forward, inverse;

	START_PROFILING("forward") {
		//fftw_plan_dft_1d
		forward = fftw_plan_dft_r2c_1d(N, x, X, FFTW_ESTIMATE);
		fftw_execute(forward);
	}
	STOP_PROFILING();

	START_PROFILING("inverse") {
		inverse = fftw_plan_dft_c2r_1d(N, X, x, FFTW_ESTIMATE);
		fftw_execute(forward);
	}
	STOP_PROFILING();

	START_PROFILING("free_mem") {
		fftw_destroy_plan(forward);
		fftw_destroy_plan(inverse);
		fftw_free(x);
		x = NULL;
		fftw_free(X);
		X = NULL;
	}
	STOP_PROFILING();
}

void cfft(const std::size_t &N)
{

}

DEF_FUNCS_(1024);

//DEF_FUNCS_(1023)
//{
//
//}
//
//DEF_FUNCS_(65536)
//{
//
//}
//
//DEF_FUNCS_(65535)
//{
//
//}
//
//DEF_FUNCS_(4294967296)
//{
//
//}
//
//DEF_FUNCS_(4294967295)
//{
//
//}