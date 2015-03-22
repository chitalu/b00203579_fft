#include <cstdio>

#include "base.h"

#define REG_FUNC(arg)                                                          \
  g_fft_funcs["real_fft_op_" #arg] = real_fft_op_##arg;                        \
  g_fft_funcs["complex_fft_op_" #arg] = complex_fft_op_##arg;

std::map<std::string, fft_func_t> g_fft_funcs;

void init(void) {
  REG_FUNC(1023);
  REG_FUNC(1024);

  /*REG_FUNC(65536);
  REG_FUNC(65535);

  REG_FUNC(4294967296);
  REG_FUNC(4294967295);*/
}

int main(int argc, char const *argv[]) {
  printf("%s\n\n", "hello fftw!");

  printf("runs per-analysis func: %d\n\n", MAX_TIME_SAMPLES);

  // initialise function pointer vars
  init();

  // loop for the analysis function
  for (std::map<std::string, fft_func_t>::const_iterator f_iter =
           g_fft_funcs.begin();
       f_iter != g_fft_funcs.cend(); ++f_iter) {

    // print the name of the fft function about to run
    printf("run: %s\n", f_iter->first.c_str());

    // to get a good estimate call fft analysis function multiple times
    // to determine an average value
    int n = 0;
    while (n++ < MAX_TIME_SAMPLES) {
      // call the fft function we want to analyse
      f_iter->second();
    }
  }

  struct {
    double operator()(const std::list<double> &arg) {
      double out = 0;
      for (std::list<double>::const_iterator i = arg.begin(); i != arg.end();
           ++i) {
        out += *i;
      }
      out /= arg.size();
      return out;
    }
  } mean_reduce;

  printf("\nelapsed times:\n");
  for (ptime_iter_t i = g_tstamps.cbegin(); i != g_tstamps.cend(); ++i) {
    printf("task: %s\ntime: %f microseconds\n\n", i->first.c_str(),
           mean_reduce(i->second));
  }
  return 0;
}