// This work is submitted in partial fulfillment of the requirements 
// for the degree of BSc (Hons) Computer Games Technology in the University 
// of the West of Scotland.
//
// I declare that this work embodies the results of my own work and 
// that it has been composed by me. Following normal academic conventions, 
// I have made due acknowledgement to the work of others.
//
// Name: FLOYD MULENGA CHITALU

#include <cstdio>
#include <ctime>
#include <sstream>
#include <fstream>

#include "base.h"

#define REG_FUNC(arg)                                                          \
  g_fft_funcs["real_fft_op_" #arg] = real_fft_op_##arg;                        \
  g_fft_funcs["complex_fft_op_" #arg] = complex_fft_op_##arg;

// extern...
std::map<std::string, fft_func_t> g_fft_funcs;
std::map<std::string, std::vector<double>> g_op_stats;
std::string g_last_profiled_task;

unsigned int g_planner_flag;

unsigned program_seed = 0;

LARGE_INTEGER g_system_clock_freq;

// generate value between 0.0f and 1.0f
float rand_norm(void) {
  return static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
}

// generate value between 0.0f and "hi"
float rand_1(float hi) {
  return static_cast<float>(rand()) / (static_cast<float>(RAND_MAX / hi));
}

// generate value between "lo" and "hi"
float rand_2(float lo, float hi) {
  return lo +
         static_cast<float>(rand()) /
             (static_cast<float>(RAND_MAX / (hi - lo)));
}

// this function simply saves the profiling data that has been
// generated while the program was running
void save_pdata(void) {
  using namespace std;

  struct {
    // reduce list elements by summation
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

  // anonymous struct used to parse strings for particular values
  struct {
    string num_samples(const string &s) {
      return s.substr(0, s.find_first_of('-'));
    }

    string op_stats(const string &s) {
      return s.substr(s.find_first_of(',') + 1);
    }

    string transform_dir(const string &s) {
      return (s.find("forward") != string::npos ? "forward" : "inverse");
    }

    string input_data_type(const string &s) {
      std::size_t off = s.find_first_of('-') + 1;
      string subs = s.substr(off, 4);
      if (subs == "rfft") {
        subs = "real";
      } else {
        subs = "complex";
        if (s.find("::f1") != string::npos) {
          subs.append("-ro");
        }
      }

      return subs;
    }
  } get_;

  // map containing per file data
  std::map<string, stringstream> file_data;

  // loop for the time stamps per task
  for (ptime_iter_t i = g_tstamps.cbegin(); i != g_tstamps.cend(); ++i) {
    string gname = i->first;
    if (gname.find("crte") != string::npos) {
      continue;
    }

    double time_taken = mean_reduce(i->second);

    string input_data_type = get_.input_data_type(gname);

    string fname(input_data_type + "-");

    string transform_dir = get_.transform_dir(gname);
    fname.append(transform_dir + ".csv");

    string samples = get_.num_samples(gname);

    if (file_data.find(fname) == file_data.end()) {
      printf("\npreparing file: %s\n", fname.c_str());
      file_data[fname]
          << ("samples, adds, multiplies, fused-madds, time (microsecs)\n");
    }

    file_data[fname] << samples << ", ";

    bool found = false;
	// now we loop for the FLOPS stats report by every FFTw plan used
	// during execution
    for (op_stats_iter_t i_ = g_op_stats.cbegin(); i_ != g_op_stats.cend();
         ++i_) {
      string gname_ = i_->first;

      if (samples == get_.num_samples(gname_) &&
          input_data_type.find(get_.input_data_type(gname_).substr(0, 4)) !=
              string::npos) {
        if ((transform_dir == "forward" &&
             gname_.find("fplan") != string::npos) ||
            (transform_dir == "inverse" &&
             gname_.find("iplan") != string::npos)) {
		  // type "long" is only 4 bytes on windows!!
          file_data[fname] << (long long)i_->second[ADD_OPERATIONS] << ", "
                           << (long long)i_->second[MUL_OPERATIONS] << ", "
                           << (long long)i_->second[FMA_OPERATIONS] << ", ";
          found = true;
          break;
        } else {
          continue;
        }
      }
    }

    assert(found);

    file_data[fname] << time_taken << "\n";
  }

  // write data associated with every file 
  for (std::map<string, stringstream>::iterator i = file_data.begin();
       i != file_data.end(); ++i) {
    ofstream f;

    f.open(i->first);

    if (f.is_open()) {
      f << i->second.str();
    }

    f.close();
  }
}

void cmd_args(int argc, char const *argv[]) {
  static std::map<std::string, unsigned int> flags;
  flags.insert(std::make_pair("-lo", FFTW_ESTIMATE));
  flags.insert(std::make_pair("-mid", FFTW_MEASURE));
  flags.insert(std::make_pair("-hi", FFTW_EXHAUSTIVE));

  for (int i(1); i < argc; ++i) {
    if (std::string(argv[i]) == "-h") {
      printf("USAGE\n"
             "'-hi' :: optimize for high performance during "
             "execution\n\t(takes a long time to initialise)\n"
             "'-mid' :: similar to \"-hi\" but takes slightly less time\n"
             "'-lo' :: fastest initialisation times but results in "
             "slow\n\tperfomance during transformation\n"
             "'<arg>' :: specify seed value used to initialise "
             "srand\n\texample: 'B00203579_fft.exe 219'\n"
             "'-h' :: display this menu\n\n");

      exit(0);
    }

    bool found = false;
    for (std::map<std::string, unsigned int>::iterator it = flags.begin();
         it != flags.end(); ++it) {
      if (std::string(argv[i]) == it->first) {
        found = true;
        g_planner_flag = it->second;
        break;
      }
    }
    if (found)
      continue;

    int seed = atoi(argv[i]);
    if (seed > 0) {
      program_seed = seed;
      continue;
    }

    fprintf(stderr, "invalid arg: %s", argv[i]);
    exit(1);
  }
}

int main(int argc, char const *argv[]) {
  printf("%s\n\n", "hello fftw!");

  printf("runs per-analysis func: %d\n\n", MAX_FUNC_RUNS);

  //did the user specify any arguments
  if (argc)
    cmd_args(argc, argv);

  if (program_seed == 0) {
    program_seed = static_cast<unsigned>(time(NULL));
    srand(program_seed);
  }

  printf("using seed: %u\n\n", program_seed);

  // Retrieves the frequency of the performance counter. The frequency
  // of the performance counter is fixed at system boot and is consistent
  // across all processors. Therefore, the frequency need only be queried
  // upon application initialization, and the result can be cached.
  //**though the value is fixed i still keep querying, it works the same!

  if (!QueryPerformanceFrequency(&g_system_clock_freq)) {
    fprintf(stderr, "QueryPerformanceFrequency failed!\n");
    exit(1);
  }
  assert(g_system_clock_freq.QuadPart != 0);

  printf("clock tick frequency: %lld\n\n", g_system_clock_freq.QuadPart);

  // initialise function pointer vars denoting the sequence lengths too!

  REG_FUNC(16); // 2^4
  REG_FUNC(32); // 2^5

  REG_FUNC(258); // 2^8 + 2
  REG_FUNC(512); // 2^9

  REG_FUNC(1024); // 2^10
  REG_FUNC(1026); // 2^10 + 2

  REG_FUNC(4096); // 2^12
  REG_FUNC(4098); // 2^12 + 2

  REG_FUNC(16384); // 2^14
  REG_FUNC(16386); // 2^14 + 2

  REG_FUNC(65536); // 2^16
  REG_FUNC(65538); // 2^16 + 2

  // run every analysis function
  for (std::map<std::string, fft_func_t>::const_iterator f_iter =
           g_fft_funcs.begin();
       f_iter != g_fft_funcs.cend(); ++f_iter) {

    // print the name of the fft function about to run
    printf("run: %s\n", f_iter->first.c_str());
    // call the fft function we want to analyse
    f_iter->second();
  }

  save_pdata();
  return 0;
}