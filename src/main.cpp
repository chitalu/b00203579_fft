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

  stringstream strm;

  string fname0("fftw-time-stats.csv");
  string pdata0("Task, time (ms)\n");

  for (ptime_iter_t i = g_tstamps.cbegin(); i != g_tstamps.cend(); ++i) {
    strm << i->first.c_str() << ", " << mean_reduce(i->second) << "\n";
  }
  pdata0.append(strm.str());

  strm.clear(); // clear any bits set
  strm.str(std::string());

  ofstream fhndl0, fhndl1;

  fhndl0.open(fname0);

  if (fhndl0.is_open()) {
    fhndl0 << pdata0;
  }
  fhndl0.close();

  string fname1("fftw-op-stats.csv");
  string pdata1("fftw-plan, add ops, mul ops, f-madd ops\n");

  for (op_stats_iter_t i = g_op_stats.cbegin(); i != g_op_stats.cend(); ++i) {
    strm << i->first.c_str() << ", " << i->second[ADD_OPERATIONS] << ", "
         << i->second[MUL_OPERATIONS] << ", " << i->second[FMA_OPERATIONS]
         << "\n";
  }
  pdata1.append(strm.str());

  fhndl1.open(fname1);
  if (fhndl1.is_open()) {
    fhndl1 << pdata1;
  }
  fhndl1.close();
}

void cmd_args(int argc, char const *argv[])
{
	static std::map<std::string, unsigned int> flags;
	flags.insert(std::make_pair("-lo", FFTW_ESTIMATE));
	flags.insert(std::make_pair("-mid", FFTW_MEASURE));
	flags.insert(std::make_pair("-hi", FFTW_EXHAUSTIVE));

	for (int i(1); i < argc; ++i)
	{
		if (std::string(argv[i]) == "-h")
		{
			printf(
				"USAGE\n"
				"'-hi' :: optimize for high performance during execution\n\t(takes a long time to initialise)\n"
				"'-mid' :: similar to \"-hi\" but takes slightly less time\n"
				"'-lo' :: fastest initialisation times but results in slow\n\tperfomance during transformation\n"
				"'<arg>' :: specify seed value used to initialise srand\n\texample: 'B00203579_fft.exe 219'\n"
				"'-h' :: display this menu\n\n"
				);

			exit(0);
		}

		bool found = false;
		for (std::map<std::string, unsigned int>::iterator it = flags.begin();
			it != flags.end(); ++it)
		{
			if (std::string(argv[i]) == it->first)
			{
				found = true;
				g_planner_flag = it->second;
				break;
			}
		}
		if (found)	continue;

		int seed = atoi(argv[i]);
		if (seed > 0)
		{
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

  if (argc)
	  cmd_args(argc, argv);

  if (program_seed == 0)
  {
	  program_seed = static_cast<unsigned>(time(NULL));
	  srand(program_seed);
  }

  printf("using seed: %u\n\n", program_seed);

  // initialise function pointer vars
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