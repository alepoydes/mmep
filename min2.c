#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <getopt.h>

#include "vector.h"
#include "skyrmion.h"
#include "optim.h"
#include "plot.h"
#include "parse.h"

float epsilon=1e-6;
int max_iter=1000;
float mode_param=2;
int mode=SDM_CONSTANT;
int debug_plot=0;
int debug_every=1;

void skyrmion_display(int iter, float* restrict a, float* restrict grad_f, float f, float res, float alpha) {
  static float prev_f=NAN;
  static float prev_iter=-1;
  if(iter>=prev_iter+debug_every) {
  	fprintf(stderr, "%d: E %g%+g R %g A %g\n", iter, f, f-prev_f, res, alpha);
  	if(debug_plot) plot_field3(stdout,a);
  };
  prev_f=f; prev_iter=iter;
};

int skyrmion_lagrange_conjugate(float* restrict x, int mode, float mode_param, 
	float epsilon, int max_iter) 
{
	return lagrange_conjugate(
		sizeu*sizex*sizey*sizez*3,  // int N 
		sizeu*sizex*sizey*sizez, // int M 
		x, // float* x0, 
        hamiltonian_hessian, // void (*Q)(const float* x, float* y)
        subtract_field, // void (*L)(float* y),
        mode, // int mode
        mode_param, // float mode_param
        epsilon, // float epsilon
        max_iter, // int max_iter
        skyrmion_display, // void (*display)(int iter, float* a, float* grad_f, float f, float res, float alpha),
        skyrmion_constrain, // void (*C)(const float* x, float* r),
        skyrmion_constrain_gradient, // void (*D)(const float* x, const float* u, float* r),
        skyrmion_constrain_adjucent // void (*P)(const float* x, const float* y, float* r)
    );
};

void showUsage(const char* program) {
  fprintf(stderr, "Compute stable states.\
\nUsage:\
\n    %s [options] [lattice description file]\
\nOptions:\
\n   -h|--help              Show this message and exit\
\n   -p|--plot              Show graphics\
\n   -e|--epsilon REAL      Desired residual\
\n   -i           INT       Set maximum number of iterations\
\n   -r           INT       Progress will be shown every given iteration\
\n   -m|--mode    INT       Optimization method\
\n   -a           REAL      A parameter for optimization methods\
\n", program);
};

int parseCommandLine(int argc, char** argv) {
  int c;
  while(1) {
    static struct option long_options[] = {      
      {"help", no_argument, 0, 'h'},
      {"plot", no_argument, 0, 'p'},
      {"epsilon", required_argument, 0, 'e'},
      {"mode", required_argument, 0, 'm'},
      {0, 0, 0, 0}
    };
    int option_index = 0;
    c = getopt_long(argc, argv, "hpe:i:r:m:a:", long_options, &option_index);
    if (c==-1) break;
    switch(c) {
      case 'p': debug_plot=1; break;
      case 'h': showUsage(argv[0]); exit(0);
      case 'e': epsilon=atof(optarg); break;
      case 'i': max_iter=atoi(optarg); break;
      case 'r': debug_every=atoi(optarg); break;
      case 'm': mode=atoi(optarg); break;
      case 'a': mode_param=atof(optarg); break;
      case '?': break;
      default: fprintf(stderr,"Unprocessed option '%c'\n", c); exit(1);
    };
  };
  return optind;
};

int main(int argc, char** argv) {
  srand(time(NULL));
  int i=parseCommandLine(argc,argv);
  if(i<argc) {
  	FILE* file=fopen(argv[i],"r");
  	if(!file) { fprintf(stderr, "Can not open file '%s'\n", argv[i]); exit(1); };
  	parse_lattice(file);
  	fclose(file);
  } else parse_lattice(stdin);
  int size=sizeu*sizex*sizey*sizez*3;
  float* spins=(float*)malloc(sizeof(float)*size); assert(spins);
  random_vector(size, spins); normalize(spins);
  int status=skyrmion_lagrange_conjugate(spins, mode, mode_param, epsilon, max_iter);
  fprintf(stderr, "Status: %d\n", status);
  free(spins);
  return 0;
};