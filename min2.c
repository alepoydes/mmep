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

real epsilon=1e-6;
int max_iter=1000;
real mode_param=0.1;
int mode=SDM_CONSTANT;
int debug_plot=0;
int debug_every=1;

void skyrmion_display(int iter, real* restrict a, real* restrict grad_f, real f, real res, real alpha) {
  static real prev_f=NAN; static real prev_res=NAN; 
  static real prev_alpha=NAN; 
  static int prev_iter=-1;
  if(iter>=prev_iter+debug_every || iter<0) {
    // Compute projection on tangent space
    int size=sizeu*sizex*sizey*sizez*3;
    real* tmp=(real*)malloc(sizeof(real)*size); assert(tmp);
    hamiltonian_hessian(a,tmp); 
    subtract_field(tmp); 
    project_to_tangent(a, tmp); 
    real tanres=rsqrt(normsq(size, tmp)/size);
    free(tmp);
    fprintf(stderr, "%d: E %"RF"g", abs(iter), f);
    if(f<prev_f) fprintf(stderr, COLOR_GREEN"%+"RF"g "COLOR_RESET, f-prev_f);
    else fprintf(stderr, COLOR_RED"%+"RF"g "COLOR_RESET, f-prev_f);
    if(res<prev_res) fprintf(stderr, "R "COLOR_GREEN"%"RF"g"COLOR_RESET, res);
    else fprintf(stderr, "R "COLOR_RED"%"RF"g"COLOR_RESET, res);
    fprintf(stderr, "/%"RF"g ", tanres);
    if(alpha<=prev_alpha) fprintf(stderr, "A "COLOR_YELLOW"%"RF"g\n"COLOR_RESET, alpha);
    else fprintf(stderr, "A %"RF"g\n"COLOR_RESET, alpha);
  	if(debug_plot) plot_field3(stdout,a);
    prev_iter=iter; prev_res=res; prev_alpha=alpha;
  };
  //fprintf(stderr, "%d: prev iter %d\n", iter, prev_iter);  
  prev_f=f; 
};

int skyrmion_lagrange_conjugate(real* restrict x, int mode, real mode_param, 
	real epsilon, int max_iter) 
{
  if(mode==0)
	return lagrange_conjugate_quad(
		sizeu*sizex*sizey*sizez*3,  // int N 
		sizeu*sizex*sizey*sizez, // int M 
		x, // real* x0, 
        hamiltonian_hessian, // void (*Q)(const real* x, real* y)
        subtract_field, // void (*L)(real* y),
        mode, // int mode
        mode_param, // real mode_param
        epsilon, // real epsilon
        max_iter, // int max_iter
        skyrmion_display, // void (*display)(int iter, real* a, real* grad_f, real f, real res, real alpha),
        skyrmion_constrain, // void (*C)(const real* x, real* r),
        skyrmion_constrain_gradient, // void (*D)(const real* x, const real* u, real* r),
        skyrmion_constrain_adjucent, // void (*P)(const real* x, const real* y, real* r)
        -skyrmion_minimum_energy()
    );
  else if(mode==1)
  return lagrange_conjugate(
    sizeu*sizex*sizey*sizez*3,  // int N 
    sizeu*sizex*sizey*sizez, // int M 
    x, // real* x0, 
        hamiltonian_hessian, // void (*Q)(const real* x, real* y)
        subtract_field, // void (*L)(real* y),
        mode, // int mode
        mode_param, // real mode_param
        epsilon, // real epsilon
        max_iter, // int max_iter
        skyrmion_display, // void (*display)(int iter, real* a, real* grad_f, real f, real res, real alpha),
        skyrmion_constrain, // void (*C)(const real* x, real* r),
        skyrmion_constrain_gradient, // void (*D)(const real* x, const real* u, real* r),
        skyrmion_constrain_adjucent, // void (*P)(const real* x, const real* y, real* r)
        -skyrmion_minimum_energy()        
    );    
  else {
    fprintf(stderr, "Unknown solver mode: %d\n",mode);
    exit(1);
  };
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
  int status=skyrmion_lagrange_conjugate(spins, mode, mode_param, epsilon, max_iter);
  fprintf(stderr, "Status: %d\n", status);
  // Saving result
  const char* filename="fields/state.gnuplot";
  FILE* file=fopen(filename,"w");
  if(file) {
    fprintf(file,"set terminal pngcairo\n");
    fprintf(file,"set output 'fields/state.png'\n");
    plot_field3(file, spins);
    fclose(file);
  } else {
    fprintf(stderr, COLOR_RED"Can not open '%s' for writing\n"COLOR_RESET,filename);
  };
  free(spins);
  return 0;
};