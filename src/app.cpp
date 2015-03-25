#define _CRT_SECURE_NO_WARNINGS

#include <cassert>

#include "base.h"

// fft in the west!
#include "fftw3.h"

std::map<std::string, std::list<double>> g_tstamps;

std::string gen_name(const std::string &task, const std::string &dtype, const std::size_t &N)
{
	char buf[64];
	return std::string(_itoa(N, buf, 10)) + "-" + dtype + "-" + task;
}

#define PNAME(task_name_ )\
	gen_name(task_name_, __FUNCTION__, N)

void rfft(const std::size_t &N)
{
	double *x = NULL;
	fftw_complex *X = NULL;
	 
	x = (double *)fftw_malloc(sizeof(double) * N);
	X = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)* N);

	for (std::size_t n = 0u; n < N; ++n) {
		// time domain decomposed signal impulses
		x[n] = rand_norm();

		// frequency domain amplitudes
		Re(X[n]) = 0;
		Im(X[n]) = 0;
	}

	//forward...
	fftw_plan forward = fftw_plan_dft_r2c_1d(N, x, X, FFTW_ESTIMATE);
	START_PROFILING(PNAME("forward")) {
		fftw_execute(forward);
	}
	STOP_PROFILING();
	fftw_destroy_plan(forward);

	//inverse...
	fftw_plan inverse = fftw_plan_dft_c2r_1d(N, X, x, FFTW_ESTIMATE);
	START_PROFILING(PNAME("inverse")) {
		fftw_execute(inverse);
	}
	STOP_PROFILING();
	fftw_destroy_plan(inverse);

	//free resources
	fftw_free(x);
	x = NULL;
	fftw_free(X);
	X = NULL;
}

void cfft(const std::size_t &N)
{
	struct F{
	private:
		fftw_complex	*x, //input sequence
						*X;//output sequence
		//number of samples in sequences
		const std::size_t N;

		void alloc_data_mem(void)
		{
			x = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);
			assert(x != NULL && "failed to allocate memory for time domain sequence");
			X = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);
			assert(X != NULL && "failed to allocate memory for frequency domain sequence");
		}

		void dealloc_data_mem(void)
		{
			fftw_free(x);
			x = NULL;
			fftw_free(X);
			X = NULL;
		}

		void fill_data(bool real_valued = false)
		{
			for (unsigned n = 0; n < N; ++n) {
			// time domain decomposed signal impulses
			// fill with random values between 0.0 and 1.0
			Re(x[n]) = rand_norm();
			//if we are using the complex form of the dft but 
			//are using real valued data then set the imaginary part to 
			//zero
			Im(x[n]) = real_valued ? 0 : rand_norm();

			// frequency domain amplitudes
			// fill with zeroes1
			Re(X[n]) = 0;
			Im(X[n]) = 0;
			}
		}

		void wipe_data(bool real_valued = false)
		{
			for (unsigned n = 0; n < N; ++n) {
				// time domain decomposed signal impulses
				Re(x[n]) = 0;
				if(!real_valued) Im(x[n]) = 0;

				// frequency domain amplitudes
				Re(X[n]) = 0;
				Im(X[n]) = 0;
			}
		}

	public:
		F(const std::size_t &N_): x(NULL), X(NULL), N(N_){
			alloc_data_mem();
		}

		~F(void){
			dealloc_data_mem();
		}

		//complex valued complex fft i.e effectively two input signals
		//one for real values and the other for the imaginary
		void f0(void)
		{
			//forward...
			fftw_plan forward_fft_plan = fftw_plan_dft_1d(N, x, X, FFTW_FORWARD, FFTW_ESTIMATE);
			START_PROFILING(PNAME("complex_forward")) {
			  fftw_execute(forward_fft_plan);
			}
			STOP_PROFILING();
			fftw_destroy_plan(forward_fft_plan);
			
			//inverse...
			fftw_plan inverse_fft_plan = fftw_plan_dft_1d(N, X, x, FFTW_BACKWARD, FFTW_ESTIMATE);
			START_PROFILING(PNAME("complex_inverse")) {
			  fftw_execute(inverse_fft_plan);
			}
			STOP_PROFILING();
			fftw_destroy_plan(inverse_fft_plan);
		}

		//real valued complex fft i.e one signal stored in the real component of 
		//every element in the input sequence
		void f1(void)
		{
			//forward...
			fftw_plan forward_fft_plan = fftw_plan_dft_1d(N, x, X, FFTW_FORWARD, FFTW_ESTIMATE);
			START_PROFILING(PNAME("rcomplex_forward")) {
			  fftw_execute(forward_fft_plan);
			}
			STOP_PROFILING();
			fftw_destroy_plan(forward_fft_plan);
			
			//inverse...
			fftw_plan inverse_fft_plan = fftw_plan_dft_1d(N, X, x, FFTW_BACKWARD, FFTW_ESTIMATE);
			 START_PROFILING(PNAME("rcomplex_inverse")) {
			 
			  fftw_execute(inverse_fft_plan);
			}
			 STOP_PROFILING();
			 
			fftw_destroy_plan(inverse_fft_plan);
		}

		//use function call operator to kill two birds with one stone
		void operator()(void)
		{
			//fill data for complex valued fft
			fill_data();
			//execute forward and inverse operation
			f0();
			//wipe generated temporary results
			wipe_data();

			//---------------------------------------

			//fill data for real valued complex fft
			fill_data(true);
			//execute forward and inverse operation
			f1();
			//wipe generated temporary results
			wipe_data();
		}
	}cfft_func(N);

	//do analysis...
	cfft_func();
}

DEF_FUNCS_(1023);
DEF_FUNCS_(1024);

//DEF_FUNCS_(1023)
//DEF_FUNCS_(65536)
//DEF_FUNCS_(65535)
//DEF_FUNCS_(4294967296)
//DEF_FUNCS_(4294967295)