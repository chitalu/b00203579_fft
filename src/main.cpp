#include <cstdio>
//complex valued fft
#include "fftw3.h"

#define N (128)

int main(int argc, char const *argv[])
{
	printf("%s\n", "hello fft!");

	/*fftw_complex in[N], out[N];
    fftw_plan p;
  
    p = fftw_create_plan(N, FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_one(p, in, out);

    fftw_destroy_plan(p);*/

	return 0;
}