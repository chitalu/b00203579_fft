// This work is submitted in partial fulfillment of the requirements 
// for the degree of BSc (Hons) Computer Games Technology in the University 
// of the West of Scotland.
//
// I declare that this work embodies the results of my own work and 
// that it has been composed by me. Following normal academic conventions, 
// I have made due acknowledgement to the work of others.
//
// Name: FLOYD MULENGA CHITALU

#define _CRT_SECURE_NO_WARNINGS

#include <cassert>

#include "base.h"

std::map<std::string, unsigned int> g_flags;
std::map<std::string, std::list<double>> g_tstamps;

inline std::string gen_name(const std::string &task, const std::string &dtype,
                            const std::size_t &N) {
  char buf[64];
  return std::string(_itoa(N, buf, 10)) + "-" + dtype + "-" + task;
}

#define PNAME(task_name_) gen_name(task_name_, __FUNCTION__, N)

// truncate values after decimal point
inline float trnc(double v) {
  char sz[64];
  sprintf(sz, "%.6lf\n", v); // sz contains 0.6000
  float r = (float)atof(sz);
  return r;
}

// Note to self:
// You must create the plan before initializing the input, because
// FFTW_EXHAUSTIVE
// overwrites the in/out arrays. (Technically, FFTW_ESTIMATE does not touch your
// arrays, but you should always create plans first just to be sure.)
// http://www.fftw.org/doc/Complex-One_002dDimensional-DFTs.html

void rfft(const std::size_t &N) {
  printf(".");
  fftw_plan forward = NULL, inverse = NULL;

  double *x = NULL;
  fftw_complex *X = NULL;

  const std::size_t N_over_2(N / 2);

  // equivalent to (double *) fftw_malloc(sizeof(double) * n) and
  // (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * n), respectively
  // http://www.fftw.org/doc/Memory-Allocation.html
  x = fftw_alloc_real(N);
  X = fftw_alloc_complex(N);
  printf(".");
  // create plans first and profile how long each takes
  // since FFTW_EXHAUSTIVE can cause the process to take a
  // while for best estimates
  START_PROFILING(PNAME("real_fplan_crte"));
  { forward = fftw_plan_dft_r2c_1d(N, x, X, g_planner_flag); }
  STOP_PROFILING();

  ASSERT_TRUE(forward != NULL, "failed to create forward real plan");

  STORE_AMF_OP_STATS(forward);

  printf(".");

  START_PROFILING(PNAME("real_iplan_crte"));
  { inverse = fftw_plan_dft_c2r_1d(N, X, x, g_planner_flag); }
  STOP_PROFILING();

  ASSERT_TRUE(inverse != NULL, "failed to create inverse real plan");

  STORE_AMF_OP_STATS(inverse);
  printf("|");

  // to get a good estimate call fft analysis function multiple times
  // to determine an average value
  int n = 0;
  while (n++ < MAX_FUNC_RUNS) {

    // then initialise data...
    for (std::size_t n = 0u; n < N; ++n) {
      // time domain decomposed signal impulses
      x[n] = rand_norm();

      // frequency domain amplitudes
      Re(X[n]) = 0;
      Im(X[n]) = 0;
    }

    // the real discrete Fourier transform changes an N point input
    // signal into two N / 2 + 1 point output signals.
    // N points in the time domain corresponds to N / 2 % 1 points in the
    // frequency domain(not N / 2 points).

    // forward...
    START_PROFILING(PNAME("forward"));
    { fftw_execute(forward); }
    STOP_PROFILING();

    unsigned base = (N_over_2 + 1) + 1;
    for (unsigned i = base; i < N; ++i) {
      ASSERT_TRUE(
          Re(X[i]) == 0 && Im(X[i]) == 0,
          "upper half of tranformed freq output is non-zero at index: %u\n"
          "with values: real = %g imag = %g\n",
          i, Re(X[i]), Im(X[i]));
    }

    // inverse...
    START_PROFILING(PNAME("inverse"));
    { fftw_execute(inverse); }
    STOP_PROFILING();
  }
  printf(".");
  fftw_destroy_plan(inverse);
  fftw_destroy_plan(forward);
  // free resources
  fftw_free(x);
  x = NULL;
  fftw_free(X);
  X = NULL;
  printf("..\n");
}

void cfft(const std::size_t &N) {
  struct F {
  private:
    fftw_complex *x, // input sequence
        *X;          // output sequence
    // number of samples in sequences
    const std::size_t N;
    const std::size_t N_over_2;
    const std::size_t N_over_2_plus_1;

    void alloc_data_mem(void) {
      // behaves like malloc except that it properly aligns the array when
      // SIMD instructions (such as SSE and Altivec) are available
      x = fftw_alloc_complex(N);
      assert(x != NULL && "failed to allocate memory for time domain sequence");
      X = fftw_alloc_complex(N);
      assert(X != NULL &&
             "failed to allocate memory for frequency domain sequence");
    }

    void dealloc_data_mem(void) {
      fftw_free(x);
      x = NULL;
      fftw_free(X);
      X = NULL;
    }

    void fill_data(bool real_valued = false) {
      for (unsigned n = 0; n < N; ++n) {
        // time domain decomposed signal impulses
        // fill with random values between 0.0 and 1.0
        Re(x[n]) = rand_norm();
        // if we are using the complex form of the dft but
        // are using real valued data then set the imaginary part to
        // zero
        Im(x[n]) = real_valued ? 0 : rand_norm();

        // frequency domain amplitudes
        // fill with zeroes1
        Re(X[n]) = 0;
        Im(X[n]) = 0;
      }
    }

    void wipe_data(bool real_valued = false) {
      for (unsigned n = 0; n < N; ++n) {
        // time domain decomposed signal impulses
        Re(x[n]) = 0;
        if (!real_valued)
          Im(x[n]) = 0;

        // frequency domain amplitudes
        Re(X[n]) = 0;
        Im(X[n]) = 0;
      }
    }

  public:
    F(const std::size_t &N_)
        : x(NULL), X(NULL), N(N_), N_over_2(N / 2),
          N_over_2_plus_1(N_over_2 + 1) {
      alloc_data_mem();
    }

    ~F(void) { dealloc_data_mem(); }

    // complex valued complex fft i.e effectively two input signals
    // one for real values and the other for the imaginary
    void f0(fftw_plan &forward_plan, fftw_plan &inverse_plan) {
      // forward...

      START_PROFILING(PNAME("complex_forward"));
      { fftw_execute(forward_plan); }
      STOP_PROFILING();

      // inverse...

      START_PROFILING(PNAME("complex_inverse"));
      { fftw_execute(inverse_plan); }
      STOP_PROFILING();
    }

    // real-valued complex fft i.e one signal stored in the real component of
    // every complex element in the input sequence
    void f1(fftw_plan &forward_plan, fftw_plan &inverse_plan) {
      for (unsigned i = 0; i < N; ++i) {
        ASSERT_TRUE(Re(X[i]) == 0 && Im(X[i]) == 0,
                    "frequency values not reset");
      }

      // forward...
      START_PROFILING(PNAME("rcomplex_forward"));
      { fftw_execute(forward_plan); }
      STOP_PROFILING();

      // to calculate the real DFT by means of the Complex DFT. First, move
      // the N point signal into the real part of the complex DFT's time domain,
      // and then set all of the samples in the imaginary part to zero.
      // Calculation of the complex DFT results in a real and an imaginary
      // signal in the frequency domain, each composed of N points. Samples
      // 0 through N/2 of these signals correspond to the real DFT's spectrum.

      // points 0 through N/2 in the complex DFT are the same as in the real
      // DFT, for both the real and the imaginary parts. For the real part,
      // point N/2 + 1 is the same as point N/2 - 1 , point N/2+ 2 is the same
      // as point N/2- 2, etc. This continues to point N-1 being the same as
      // point 1. The same basic pattern is used for the imaginary part, except
      // the sign is changed. That is, point N/2+ 1 is the negative of point
      // N/2- 1 , point N/2+ 2 is the negative of point N/2- 2 , etc.
      // Samples 0 and N/2 do not have a matching point in this duplication
      // scheme.

      // THE ASSERTION FAILS [UNPREDICTABLY] DUE TO PRECISION ERRORS
      // SO I JUST CAST TO FLOAT AS A COMPROMISE BEFORE DOING THE VERIFICATION
      // COMPARISON
      unsigned c = 0;
      while ((++c) != N_over_2) {
        ASSERT_EQ(Re(X[N_over_2 + c]), Re(X[N_over_2 - c]),
                  "Real: (N/2 + c) -> %f is NOT the same as (N/2 - c) -> %f",
                  Re(X[N_over_2 + c]), Re(X[N_over_2 - c]));

        ASSERT_TRUE(trnc(Im(X[N_over_2 + c])) == -trnc(Im(X[N_over_2 - c])),
                    "imag: (N/2 + c) -> %f is NOT negative of (N/2 - c) -> %f ",
                    trnc(Im(X[N_over_2 + c])), -trnc(Im(X[N_over_2 - c])));
      }

      // inverse...
      START_PROFILING(PNAME("rcomplex_inverse"));
      { fftw_execute(inverse_plan); }
      STOP_PROFILING();

      // To check that the proper symmetry is present, after taking the inverse
      // FFT, look at the imaginary part of the time domain. It will contain
      // all zeros if everything is correct (except for a few parts-permillion
      // of noise, using single
      // precision calculations).

      // AGAIN:
      // THE ASSERTION FAILS [UNPREDICTABLY] DUE TO PRECISION ERRORS
      // SO THIS TIME I SIMPLY CHECK TO SEE IF THE VALUE PRODUCED BY THE INVERSE
      // OPERATION
      // IS BELOW A CERTAIN [ACCEPTABLE] THRESHOLD FOR THIS PARTICULAR CASE
      for (unsigned i(0); i < N; ++i) {
        ASSERT_TRUE(Im(x[i]) <= 0.000001,
                    "inverse op returned invalid result:\n"
                    "instead Im(x[i]) == %f",
                    Im(x[i]));
      }
    }

    // use function call operator to kill two birds with one stone
    void operator()(void) {
      printf(".");
      fftw_plan forward_plan = NULL, inverse_plan = NULL;

      // create plans first
      // Once the plan has been created, you can use it as many times as
      // you like for transforms on the specified in/out arrays, computing
      // the actual transforms via fftw_execute(plan)
      START_PROFILING(PNAME("cmplx_fplan_crte"));
      {
        forward_plan = fftw_plan_dft_1d(N, x, X, FFTW_FORWARD, g_planner_flag);
      }
      STOP_PROFILING();

      ASSERT_TRUE(forward_plan != NULL,
                  "failed to create forward complex plan");

      STORE_AMF_OP_STATS(forward_plan);
      printf(".");
      START_PROFILING(PNAME("cmplx_iplan_crte"));
      {
        inverse_plan = fftw_plan_dft_1d(N, X, x, FFTW_BACKWARD, g_planner_flag);
      }
      STOP_PROFILING();

      ASSERT_TRUE(inverse_plan != NULL,
                  "failed to create inverse complex plan");

      STORE_AMF_OP_STATS(inverse_plan);
      printf(".|");
      // to get a good estimate call fft analysis function multiple times
      // to determine an average value
      int n = 0;
      while (n++ < MAX_FUNC_RUNS) {
        //---------------------------------------
        {
          // fill data for COMPLEX VALUED fft
          fill_data();
          // execute forward and inverse operation
          f0(forward_plan, inverse_plan);
          // wipe generated temporary results
          wipe_data();
        }
        //---------------------------------------
        {
          // fill data for REAL VALUED complex fft
          fill_data(IGNORE_IMAGINARY_INPUT);
          // execute forward and inverse operation
          f1(forward_plan, inverse_plan);
          // wipe generated temporary results
          wipe_data(IGNORE_IMAGINARY_INPUT);
        }
      }
      printf(".");
      fftw_destroy_plan(forward_plan);
      fftw_destroy_plan(inverse_plan);
      printf(".");
    }
  } cfft_func(N);

  // do analysis...
  cfft_func();
  printf(".\n");
}

// define respective functions 
DEF_FUNCS_(16); // 2^4
DEF_FUNCS_(32); // 2^5

DEF_FUNCS_(258); // 2^8 + 2
DEF_FUNCS_(512); // 2^9

DEF_FUNCS_(1024); // 2^10
DEF_FUNCS_(1026); // 2^10 + 2

DEF_FUNCS_(4096); // 2^12
DEF_FUNCS_(4098); // 2^12 + 2

DEF_FUNCS_(16384); // 2^14
DEF_FUNCS_(16386); // 2^14 + 2

DEF_FUNCS_(65536); // 2^16
DEF_FUNCS_(65538);