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
#include "debug.h"

int size=0; // Dimenion of vector containing skyrmionic solutions
int sizep=10; // Number of nodes on path
real epsilon=1e-6;
int max_iter=1000;
real mode_param=0.2;
int mode=SDM_CONSTANT;
int debug_plot=0;
int debug_every=1;

void skyrmion_display(int iter, real* restrict a, real* restrict grad_f, real f, real res, real alpha) {
  static real prev_f=NAN;
  if(iter%debug_every==0 || iter<0) {
  	fprintf(stderr, "%d: E %"RF"g%+"RF"g R %"RF"g A %"RF"g\n", iter, f, f-prev_f, res, alpha);
  	if(debug_plot) plot_field3(stdout,a);
  };
  prev_f=f;
};

int skyrmion_steepest_descent(real* restrict x, int mode, real mode_param, 
	real epsilon, int max_iter) 
{
	return steepest_descend(
		size, (real*)x, 
		hamiltonian_hessian,	subtract_field,
		mode, mode_param, epsilon, max_iter,
		skyrmion_display, 
		normalize,project_to_tangent
	);
};

void path_hessian(const real* restrict arg, real* restrict out) {
  for(int p=0; p<sizep; p++) {
    hamiltonian_hessian(arg+size*p, out+size*p);
  };
};

void path_subtract_field(real* restrict inout) {
  for(int p=0; p<sizep; p++) {
    subtract_field(inout+size*p);
  };
};

void path_display(int iter, real* restrict mep, real* restrict grad_f, real f, real res, real alpha) {
  static real prev_f=NAN;
  if(iter%debug_every==0 || iter<0) {
    fprintf(stderr, "%d: E %"RF"g%+"RF"g R %"RF"g A %"RF"g\n", iter, f, f-prev_f, res, alpha);
    if(debug_plot) plot_path(stdout, sizep, mep);
  };
  prev_f=f;
};

void path_normalize(real* mep) {
  for(int p=0; p<sizep; p++) {
    normalize(mep+size*p);
  };
};

void path_tangent(const real* restrict mep, real* restrict grad) {
  for(int p=0; p<sizep; p++) {
    project_to_tangent(mep+size*p,grad+size*p);
  };
}

int path_steepest_descent(real* restrict path, int mode, 
  real mode_param, real epsilon, int max_iter) 
{
  return steepest_descend(
    size*sizep, (real*)path, 
    path_hessian, path_subtract_field,
    mode, mode_param, epsilon, max_iter,
    path_display, 
    path_normalize, path_tangent
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
\n   -n           INT       Nodes along path\
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
    c = getopt_long(argc, argv, "hpe:n:i:r:m:a:", long_options, &option_index);
    if (c==-1) break;
    switch(c) {
      case 'p': debug_plot=1; break;
      case 'h': showUsage(argv[0]); exit(0);
      case 'e': epsilon=atof(optarg); break;
      case 'n': sizep=atoi(optarg); break;
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
  // Read parameters
  int i=parseCommandLine(argc,argv);
  if(i<argc) {
  	FILE* file=fopen(argv[i],"r");
  	if(!file) { fprintf(stderr, "Can not open file '%s'\n", argv[i]); exit(1); };
  	parse_lattice(file);
  	fclose(file);
  } else parse_lattice(stdin);
  // Initializaton
  srand(time(NULL));
  size=sizeu*sizex*sizey*sizez*3;
  real* path=(real*)malloc(sizeof(real)*size*sizep); assert(path);
  // find two minima
  fprintf(stderr, COLOR_GREEN"Calculating initial state\n"COLOR_RESET);
  random_vector(size, path); 
  skyrmion_steepest_descent(path, mode, mode_param, epsilon, max_iter);
  fprintf(stderr, COLOR_GREEN"Calculating final state\n"COLOR_RESET);
  random_vector(size, path+size*(sizep-1)); 
  skyrmion_steepest_descent(path+size*(sizep-1), mode, mode_param, epsilon, max_iter);
  // Set initial path as geodesic approximation
  fprintf(stderr, COLOR_GREEN"Calculating MEP\n"COLOR_RESET);
  skyrmion_geodesic(sizep, path);
  // MEP calculation
  path_steepest_descent(path, mode, mode_param, epsilon, max_iter);
  // Ouput result
  const char* filename="../fields/mep.gnuplot";
  FILE* file=fopen(filename,"w");
  if(file) {
    animate_path(file, sizep, path, "");
    fclose(file);
  } else {
    fprintf(stderr, COLOR_RED"Can not open '%s' for writing\n"COLOR_RESET,filename);
  };
  // Deinitialization
  free(path);
  return 0;
};