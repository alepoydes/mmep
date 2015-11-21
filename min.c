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

real epsilon=1e-6;
int max_iter=1000;
real mode_param=0.05;
int mode=SDM_CONSTANT;
int debug_plot=0;
int debug_every=1;

void skyrmion_display(int iter, real* restrict a, real* restrict grad_f, real f, real res, real alpha) {
  static real prev_f=NAN;
  if(iter%debug_every==0) {
  	fprintf(stderr, "%d: E %"RF"g%+"RF"g R %"RF"g A %"RF"g\n", iter, f, f-prev_f, res, alpha);
  	if(debug_plot) plot_field3(stdout,a);
  };
  prev_f=f;
};

int skyrmion_steepest_descent(real* restrict x, int mode, real mode_param, 
	real epsilon, int max_iter) 
{
	return steepest_descend(
		sizeu*sizex*sizey*sizez*3, (real*)x, 
		hamiltonian_hessian,	subtract_field,
		mode, mode_param, epsilon, max_iter,
		skyrmion_display, 
		normalize,project_to_tangent
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
  real* spins=(real*)malloc(sizeof(real)*size); assert(spins);
  if(initial_state) {
    copy_vector(size, initial_state, spins);
    free(initial_state);
  } else random_vector(size, spins); 
  skyrmion_steepest_descent(spins, mode, mode_param, epsilon, max_iter);
  free(spins);
  return 0;
};