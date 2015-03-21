#ifndef __BASE_H__
#define __BASE_H__

#include <string>
#include <map>
#include <cstdint>
#include <Windows.h>

std::map<std::string, double> g_tstamps;
typedef std::map<std::string, double>::const_iterator ptime_iter_t;

struct cprofile_t {
  cprofile_t(const char *desc_in) : desc(desc_in), m_PC_freq(0.0), 
	  m_start_time(0) { start_counter(); }
  ~cprofile_t() {
    // store time results when object leaves macrro scope
    g_tstamps.insert(std::make_pair(desc, get_counter()));
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
    cprofile_t(#name_);

// end profiling codeblock
#define STOP_PROFILING() }

#define REAL_PART 0
#define IMAGINARY_PART 1

#define Imag(arg_) arg_[IMAGINARY_PART]
#define Real(arg_) arg_[REAL_PART]

const double two_pi = M_PI * 2.0;

#endif