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

// fft in the west!
#include "fftw3.h"

#define MAX_FUNC_RUNS (128)

#define ASSERT_TRUE(cond_, msg_, ...)\
do{\
	if(!(cond_)){\
	fprintf(stderr, "runtime error: \n@  -> " #cond_ "\nmsg: " msg_ "\n", ##__VA_ARGS__);\
	__debugbreak();\
	}\
}while(false);

#define ASSERT_EQ(arg0, arg1, msg_, ...)\
	ASSERT_TRUE(fabs(arg0 - arg1) < 0.000001, msg_, ##__VA_ARGS__);

#define IGNORE_IMAGINARY_INPUT (true)

typedef void (*fft_func_t)(void);
extern std::map<std::string, fft_func_t> g_fft_funcs;

extern std::map<std::string, std::list<double>> g_tstamps;
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

extern unsigned int g_planner_flag;

struct cprofile_t {
  cprofile_t(const std::string &desc_in)
      : desc(desc_in), m_PC_freq(0.0), m_start_time(0) {
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
  double m_PC_freq;
  __int64 m_start_time;

  // records the number of ticks the performance counter has
  // in the start_time variable
  void start_counter(void) {
    LARGE_INTEGER li;
    // Retrieves the frequency of the performance counter. The frequency
    // of the performance counter is fixed at system boot and is consistent
    // across all processors. Therefore, the frequency need only be queried
    // upon application initialization, and the result can be cached.
    //**though the value is fixed i still keep querying, it works the same!
    if (!QueryPerformanceFrequency(&li))
      printf("QueryPerformanceFrequency failed!\n");

    m_PC_freq = double(li.QuadPart) / 1000000.0;

    QueryPerformanceCounter(&li);
    m_start_time = li.QuadPart;
  }

  // returns the number of milliseconds since start_counter() was
  // last called as a double, so if get_counter() returns 0.001 then
  // it has been about 1 microsecond since start_counter() was called
  double get_counter(void) {
    LARGE_INTEGER li;
    // On a multiprocessor computer, it should not matter which processor
    // is called. However, you can get different results on different processors
    // due to bugs in the basic input/output system (BIOS) or the hardware
    // abstraction layer (HAL).
    QueryPerformanceCounter(&li);
    return double(li.QuadPart - m_start_time) / m_PC_freq;
  }
};

// start profiling codeblock
#define START_PROFILING(name_)                                                 \
  {                                                                            \
    cprofile_t(name_);

// end profiling codeblock
#define STOP_PROFILING() }

#define REAL_PART 0
#define IMAGINARY_PART 1
#define Im(arg_) arg_[IMAGINARY_PART]
#define Re(arg_) arg_[REAL_PART]

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