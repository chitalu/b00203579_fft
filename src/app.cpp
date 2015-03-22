#include "base.h"

// fft in the west!
#include "fftw3.h"

std::map<std::string, std::list<double>> g_tstamps;

std::string gen_name(const std::string &task, const std::string &dtype, const std::size_t &N)
{
	char buf[64];
	return std::string(_itoa(N, buf, 10)) + "::" + dtype + "::" + task;
}

#define PNAME(task_name_ )\
	gen_name(task_name_, __FUNCTION__, N)

void rfft(const std::size_t &N)
{
	const int freq = 440;
	const double t = 0.0001;
	const double w = M_2_PI * freq;

	double *x = NULL;
	fftw_complex *X = NULL;

	START_PROFILING(PNAME("alloc_mem")) {
		x = (double *)fftw_malloc(sizeof(double) * N);
		X = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)* N);
	}
	STOP_PROFILING();  

	START_PROFILING(PNAME("init_mem")) {
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

	START_PROFILING(PNAME("forward")) {
		//fftw_plan_dft_1d
		forward = fftw_plan_dft_r2c_1d(N, x, X, FFTW_ESTIMATE);
		fftw_execute(forward);
	}
	STOP_PROFILING();

	START_PROFILING(PNAME("inverse")) {
		inverse = fftw_plan_dft_c2r_1d(N, X, x, FFTW_ESTIMATE);
		fftw_execute(inverse);
	}
	STOP_PROFILING();

	START_PROFILING(PNAME("free_mem")) {
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
	const int freq = 440;
	const double t = 0.0001;
	const double w = M_2_PI * freq;
	// fftw_complex *x = NULL, *X = NULL;

	// START_PROFILING("allocate_mem") {
	//  x = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);
	//  X = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);
	//}
	// STOP_PROFILING()

	// START_PROFILING("initialise_mem") {
	//  for (int n = 0; n < N; ++n) {
	// // time domain decomposed signal impulses
	//    Re(x[n]) = sin(w * (double)n * t);
	//    Im(x[n]) = 0;

	// // frequency domain amplitudes
	//    Re(X[n]) = 0;
	//    Im(X[n]) = 0;
	//  }
	//}
	// STOP_PROFILING();

	// fftw_plan forward_fft_plan, inverse_fft_plan;

	// START_PROFILING("complex_fft") {
	//  forward_fft_plan = fftw_plan_dft_1d(N, x, X, FFTW_FORWARD, FFTW_ESTIMATE);
	//  fftw_execute(forward_fft_plan);
	//}
	// STOP_PROFILING();

	// START_PROFILING("complex_ifft") {
	//  inverse_fft_plan = fftw_plan_dft_1d(N, X, x, FFTW_BACKWARD,
	//  FFTW_ESTIMATE);
	//  fftw_execute(forward_fft_plan);
	//}
	// STOP_PROFILING();

	// START_PROFILING("deallocate_mem") {
	//  fftw_destroy_plan(forward_fft_plan);
	//  fftw_destroy_plan(inverse_fft_plan);
	//  fftw_free(x);
	//  x = NULL;
	//  fftw_free(X);
	//  X = NULL;
	//}
	// STOP_PROFILING();
}

DEF_FUNCS_(1023);
DEF_FUNCS_(1024);

//DEF_FUNCS_(1023)
//DEF_FUNCS_(65536)
//DEF_FUNCS_(65535)
//DEF_FUNCS_(4294967296)
//DEF_FUNCS_(4294967295)