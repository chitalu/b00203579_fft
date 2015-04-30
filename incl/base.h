// This work is submitted in partial fulfillment of the requirements 
// for the degree of BSc (Hons) Computer Games Technology in the University 
// of the West of Scotland.
//
// I declare that this work embodies the results of my own work and 
// that it has been composed by me. Following normal academic conventions, 
// I have made due acknowledgement to the work of others.
//
// Name: FLOYD MULENGA CHITALU

#ifndef __BASE_H__
#define __BASE_H__

#include <Windows.h>

#define _USE_MATH_DEFINES
#include <math.h>

#include <string>
#include <map>
#include <list>
#include <vector>

#include <cstdarg>

#include <cassert>

// fft in the west!
#include "fftw3.h"

// please note: function declarations and definitions are 
// defined using macro expansions. The reason for this is 
// to allow for automatic generation of the names used to 
// label all profiling statistics produced by the program 
// upon execution

// number of runs per number of samples. useful for 
// calculating a decent set of results.
#define MAX_FUNC_RUNS (128)

// debugging helper macro
#define ASSERT_TRUE(cond_, msg_, ...)\
do{\
	if(!(cond_)){\
	fprintf(stderr, "runtime error: \n@  -> " #cond_ "\nmsg: " msg_ "\n", ##__VA_ARGS__);\
	__debugbreak();\
	}\
}while(false);

// used for comparing floating point values. not as pedantic 
// as a regular comparison due to precesion error when
// calculating FFT operations
#define ASSERT_EQ(arg0, arg1, msg_, ...)\
	ASSERT_TRUE(fabs(arg0 - arg1) < 0.000001, msg_, ##__VA_ARGS__);

#define IGNORE_IMAGINARY_INPUT (true)

// type used to signify function type
typedef void (*fft_func_t)(void);
// storage for all analysis functions
extern std::map<std::string, fft_func_t> g_fft_funcs;

//profiling time stamps
extern std::map<std::string, std::list<double>> g_tstamps;
// interator for time stamps container
typedef std::map<std::string, std::list<double>>::const_iterator ptime_iter_t;
extern std::string g_last_profiled_task;

// stores the details of the determined number of floating-point additions, 
// multiplications, and fused multiply-add operations involved in the
// execution of a "named" plan i.e. one corresponding the number of samples
// analysed. For example, this will store the information of a plan used to execute
// an [inverse] fft with [N] samples of [complex] data. The key in the std::map
// is a string corresponding to a unique name of the task and the value is an
// array of three doubles holding the aforementioned info on the task.
// This is written to a file at program teardown.
// ...
//
// Taken from: http://www.fftw.org/doc/Using-Plans.html#Using-Plans
// void fftw_flops(const fftw_plan plan, double *add, double *mul, double *fma);
// 
// Given a plan, set add, mul, and fma to an exact count of the number of floating-point 
// additions, multiplications, and fused multiply-add operations involved in the plan's 
// execution. The total number of floating-point operations (flops) is add + mul + 2*fma, 
// or add + mul + fma if the hardware supports fused multiply-add instructions (although 
// the number of FMA operations is only approximate because of compiler voodoo). (The 
// number of operations should be an integer, but we use double to avoid overflowing int 
// for large transforms; the arguments are of type double even for single and long-double 
// precision versions of FFTW.) 
extern std::map<std::string, std::vector<double>> g_op_stats;
typedef std::map<std::string, std::vector<double>>::const_iterator op_stats_iter_t;

// just some integral constants to aid readability when storing 
// FLOPS associated with plan execution
#define ADD_OPERATIONS (0)
#define MUL_OPERATIONS (1)
#define FMA_OPERATIONS (2)

// read [add] [mul] and [fma] op stats from fftw API
// and store them into g_op_stats
#define STORE_AMF_OP_STATS(plan_)\
	do{\
		std::vector<double> ops; ops.reserve(3); ops.resize(3);\
		fftw_flops(plan_, &ops[ADD_OPERATIONS], &ops[MUL_OPERATIONS], &ops[FMA_OPERATIONS]);\
		g_op_stats.insert(std::make_pair(g_last_profiled_task, ops));\
	}while(false);

// flag determing how much we want FFTw to care about how fast we want the 
// application to be.
extern unsigned int g_planner_flag;

// system clock frequency value used for calculating timing values.
extern LARGE_INTEGER g_system_clock_freq;

// profiling struct, an instance begins the timer
// the deconstructor stores timing values
struct cprofile_t {
  cprofile_t(const std::string &desc_in)
      : desc(desc_in) {
    start_counter();
  }
  ~cprofile_t(void) {
	  g_last_profiled_task = desc;
    // store time results when object leaves macrro scope
    g_tstamps[desc].push_back(get_counter());

    if (g_tstamps[desc].size() >= MAX_FUNC_RUNS) {
      g_tstamps[desc].pop_front();
    }
  }

private:
  std::string desc;
  LARGE_INTEGER m_start_time;

  // records the number of ticks the performance counter has
  // in the start_time variable
  void start_counter(void) {
    QueryPerformanceCounter(&m_start_time);
  }

  // returns the number of milliseconds since start_counter() was
  // last called as a double, so if get_counter() returns 0.001 then
  // it has been about 1 microsecond since start_counter() was called
  double get_counter(void) {
    LARGE_INTEGER end_time, elapsed;
    // On a multiprocessor computer, it should not matter which processor
    // is called. However, you can get different results on different processors
    // due to bugs in the basic input/output system (BIOS) or the hardware
    // abstraction layer (HAL).
    QueryPerformanceCounter(&end_time);

	elapsed.QuadPart = end_time.QuadPart - m_start_time.QuadPart;

	// We now have the elapsed number of ticks, along with the
	// number of ticks-per-second. We use these values
	// to convert to the number of elapsed microseconds.
	// To guard against loss-of-precision, we convert
	// to microseconds *before* dividing by ticks-per-second.
	elapsed.QuadPart *= 1000000;

	//return micro seconds
	return ((double)elapsed.QuadPart / (double)g_system_clock_freq.QuadPart);
  }
};

// start profiling codeblock
#define START_PROFILING(name_)                                                 \
  {                                                                            \
    cprofile_t(name_);

// end profiling codeblock
#define STOP_PROFILING() }

// complex number auxilliary macros, helpful for reading and modifying an
// FFTw complex value as it would look in equation form 
#define REAL_PART 0
#define IMAGINARY_PART 1
#define Im(arg_) arg_[IMAGINARY_PART]
#define Re(arg_) arg_[REAL_PART]


// analysis function declaration macros. 

#define DECL_FUNC_(dtype, dsize) extern "C" void dtype##_fft_op_##dsize(void);

#define DECL_FUNCS_(dsize)                                                     \
  DECL_FUNC_(real, dsize);                                                     \
  DECL_FUNC_(complex, dsize)

// fft analysis function declarations
//
// each with a specifies different number of samples to be tested. 
// sizes that are products of small factors are transformed most 
// efficiently (although prime sizes still use an O(n log n) algorithm)
// http://www.fftw.org/doc/Complex-One_002dDimensional-DFTs.html

DECL_FUNCS_(16); // 2^4
DECL_FUNCS_(32); // 2^5

DECL_FUNCS_(258); // 2^8 + 2
DECL_FUNCS_(512); // 2^9

DECL_FUNCS_(1024); // 2^10
DECL_FUNCS_(1026); // 2^10 + 2

DECL_FUNCS_(4096); // 2^12
DECL_FUNCS_(4098); // 2^12 + 2

DECL_FUNCS_(16384); // 2^14
DECL_FUNCS_(16386); // 2^14 + 2

DECL_FUNCS_(65536); // 2^16
DECL_FUNCS_(65538); // 2^16 + 2

#define DEF_FUNC_(dtype, dsize) void dtype##_fft_op_##dsize(void)

#define DEF_FUNCS_(dsize)                                                      \
  DEF_FUNC_(real, dsize) { rfft(dsize); }                                      \
  DEF_FUNC_(complex, dsize) { cfft(dsize); }

//generate value between 0.0f and 1.0f
extern "C" float rand_norm(void);

//generate value between 0.0f and "hi"
extern "C" float rand_1(float hi);

//generate value between "lo" and "hi"
extern "C" float rand_2(float lo, float hi);

#endif