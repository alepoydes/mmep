#include "vector.h"
#include "skyrmion.h"
#include "optim.h"
#include "plot.h"
#include "parse.h"
#include "debug.h"
#include "octave.h"

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <getopt.h>

#define OUTDIR "fields"

real epsilon=1e-6;
int max_iter=1000;
real mode_param=0.2;
real param2=1;
int mode=2;
int debug_plot=0;
int debug_every=1000;
int save_octave=0;
real dipole_negligible=0.001;

void skyrmion_display(int iter, real* restrict a, real* restrict grad_f, real f, real res, real constres, real alpha) {
  static real prev_f=NAN; static real prev_res=NAN;
  if(iter%debug_every==0 || iter<0) {
  	fprintf(stderr, "%d: E ", abs(iter));
    watch_number(f,prev_f,16);
    fprintf(stderr, " R ");
    watch_number(res,prev_res,16);
    fprintf(stderr, " %+"RF"g", constres);
    fprintf(stderr, " A %"RF"g\n", alpha);
  	if(debug_plot) plot_field3(stdout,a);
    prev_f=f; prev_res=res;
  };
};

real quasynorm(real* restrict x) {
  return rsqrt(seminormalize(param2, x)/sizeu/sizex/sizey/sizez);
}

int skyrmion_steepest_descent(real* restrict x, int mode, real mode_param, 
	real epsilon, int max_iter) 
{
	return steepest_descend(
		SIZE*3, (real*)x, 
		hamiltonian_hessian,	subtract_field,
		mode, mode_param, epsilon, max_iter,
		skyrmion_display, 
		quasynorm,project_to_tangent
	);
};

void showUsage(const char* program) {
  fprintf(stderr, "Compute stable states.\
\nUsage:\
\n    %s [options] [lattice description file]\
\nOptions:\
\n   -h|--help              Show this message and exit\
\n   -p|--plot              Show graphics\
\n   -o                     Save result in Octave format\
\n   -e|--epsilon REAL      Desired residual\
\n   -E           REAL      Neglible value of dipole interaction\
\n   -i           INT       Set maximum number of iterations\
\n   -r           INT       Progress will be shown every given iteration\
\n   -m|--mode    INT       Optimization method\
\n   -a           REAL      First parameter for optimization methods\
\n   -b           REAL      Second parameter for optimization methods\
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
    c = getopt_long(argc, argv, "hpE:e:i:r:m:a:b:o", long_options, &option_index);
    if (c==-1) break;
    switch(c) {
      case 'p': debug_plot=1; break;
      case 'h': showUsage(argv[0]); exit(0);
      case 'e': epsilon=atof(optarg); break;
      case 'E': dipole_negligible=atof(optarg); break;      
      case 'i': max_iter=atoi(optarg); break;
      case 'r': debug_every=atoi(optarg); break;
      case 'm': mode=atoi(optarg); break;
      case 'a': mode_param=atof(optarg); break;
      case 'b': param2=atof(optarg); break;
      case 'o': save_octave=1; break;
      case '?': break;
      default: fprintf(stderr,"Unprocessed option '%c'\n", c); exit(1);
    };
  };
  return optind;
};

int main(int argc, char** argv) {
  srand(time(NULL));
  init_signal();
  int i=parseCommandLine(argc,argv);
  if(i<argc) {
  	FILE* file=fopen(argv[i],"r");
  	if(!file) { fprintf(stderr, "Can not open file '%s'\n", argv[i]); exit(1); };
  	parse_lattice(file);
  	fclose(file);
  } else parse_lattice(stdin);

  prepare_dipole_table(dipole_negligible);

  int size=SIZE*3;
  real* spins=(real*)malloc(sizeof(real)*size); assert(spins);
  if(initial_state) {
    copy_vector(size, initial_state, spins);
    free(initial_state);
  } else {
    random_vector(size, spins); 
    /*forall(u,x,y,z) {
      int i=3*INDEX(u,x,y,z);
      real s,c; rsincos(2*M_PI*(x-y)/(real)sizex,&s,&c);
      spins[i]=-s; spins[i+1]=s; spins[i+2]=c;
    };*/
  };
  skyrmion_steepest_descent(spins, mode, mode_param, epsilon, max_iter);
  // Saving result
  real energy[5]; skyrmion_energy(spins, energy);
  fprintf(stderr,COLOR_YELLOW"Zeeman energy:"COLOR_RESET" %.10"RF"f\n",energy[1]);
  fprintf(stderr,COLOR_YELLOW"Exchange energy:"COLOR_RESET" %.10"RF"f\n",energy[2]);
  fprintf(stderr,COLOR_YELLOW"Anisotropy energy:"COLOR_RESET" %.10"RF"f\n",energy[0]);
  fprintf(stderr,COLOR_YELLOW"Dzyaloshinskii-Moriya energy:"COLOR_RESET" %.10"RF"f\n",energy[3]);
  fprintf(stderr,COLOR_YELLOW"Dipole interaction energy:"COLOR_RESET" %.10"RF"f\n",energy[4]);
  fprintf(stderr,COLOR_YELLOW"TOTAL energy:"COLOR_RESET" %.10"RF"f\n",energy[0]+energy[1]+energy[2]+energy[3]+energy[4]);  

  fprintf(stderr, COLOR_YELLOW"Saving gnuplot\n"COLOR_RESET);
  const char* filename=OUTDIR"/state.gnuplot";
  FILE* file=fopen(filename,"w");
  if(file) {
    fprintf(file,"set terminal pngcairo\n");
    fprintf(file,"set output '"OUTDIR"/state.png'\n");
    plot_field3(file, spins);
    fclose(file);
  } else {
    fprintf(stderr, COLOR_RED"Can not open '%s' for writing\n"COLOR_RESET,filename);
  };
  // Octave output
  if(save_octave) {
    const char* octfile=OUTDIR"/state.oct";
    fprintf(stderr, COLOR_YELLOW"Saving %s\n"COLOR_RESET,octfile);
    file=fopen(octfile,"w");
    if(file) {
      oct_save_init(file);
      oct_save_hessian(file);
      oct_save_linear(file);
      oct_save_state(file,"initial",spins);
      oct_save_vertices(file);
      oct_save_finish(file);
      fclose(file);
    } else 
      fprintf(stderr, COLOR_RED"Can not open '%s' for writing\n"COLOR_RESET,octfile);
  };
  // deinitializing
  free(spins);
  return 0;
};