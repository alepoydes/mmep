#include "vector.h"
#include "skyrmion.h"
#include "optim.h"
#include "plot.h"
#include "parse.h"
#include "debug.h"
#include "octave.h"
#include "integra.h"

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <getopt.h>

#define OUTDIR "fields"

real epsilon=1e-6;
int max_iter=1000;
real mode_param=0.2;
real param2=NAN;
int mode=2;
int integrator=0;
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
  normalize(x);
  return 0;
}

void integrate(int N, void (*F)(const real* x, real* g, real* E), real T, real* X, real* E, int* iter) {
  real err=0;
  switch(integrator) {
    case 0: euler(N, F, T, X, E, iter); break;
    case 1: runge_kutta(N, F, T, X, E, iter); break;
    case 2: err=gauss_integrator(N, F, T, X, 1e-7, 10, E, iter); break;
    case 3: err=radau_integrator(N, F, T, X, 1e-7, 10, E, iter); break;
    default: fprintf(stderr, COLOR_RED"Unknown integrator:"COLOR_RESET" %d\n", integrator); exit(1);
  };
  if(err>1e-5) 
    fprintf(stderr, COLOR_RED"Warning:"COLOR_RESET" Runge-Kutta convergence error: %"RF"g\n", err);
}

void projected_gradient(const real* restrict arg, real* restrict grad, real* restrict energy) {
    skyrmion_gradient(arg, grad, energy);
    project_to_tangent(arg, grad);
};

int skyrmion_steepest_descent(real* restrict x) {
  return flow_descend(
    SIZE*3, (real*)x, 
    projected_gradient,
    mode, mode_param, epsilon, max_iter,
    skyrmion_display, 
    quasynorm,
    integrate
  );
/*
	return steepest_descend(
		SIZE*3, (real*)x, 
		skyrmion_gradient,
		mode, mode_param, epsilon, max_iter,
		skyrmion_display, 
		quasynorm,project_to_tangent
	);
*/
};

void showUsage(const char* program) {
  fprintf(stderr, "Compute stable states.\
\nUsage:\
\n    %s [options] [lattice description file]\
\nOptions:\
\n   -h|--help              Show this message and exit\
\n   -p|--plot              Show graphics\
\n   -O                     Save result in Octave format\
\n   -o                     Save result in Octave format but Hessian matrix\
\n   -e|--epsilon REAL      Desired residual\
\n   -E           REAL      Neglible value of dipole interaction\
\n   -i           INT       Set maximum number of iterations\
\n   -r           INT       Progress will be shown every given iteration\
\n   -m|--mode    INT       Optimization method\
\n   -I           INT       Integrator\
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
    c = getopt_long(argc, argv, "hpE:e:i:I:r:m:a:b:oO", long_options, &option_index);
    if (c==-1) break;
    switch(c) {
      case 'p': debug_plot=1; break;
      case 'h': showUsage(argv[0]); exit(0);
      case 'e': epsilon=atof(optarg); break;
      case 'E': dipole_negligible=atof(optarg); break;      
      case 'i': max_iter=atoi(optarg); break;
      case 'I': integrator=atoi(optarg); break;
      case 'r': debug_every=atoi(optarg); break;
      case 'm': mode=atoi(optarg); break;
      case 'a': mode_param=atof(optarg); break;
      case 'b': param2=atof(optarg); break;
      case 'o': save_octave=1; break;      
      case 'O': save_octave=2; break;      
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
  if(active) fprintf(stderr, "Active spins: %d / %d\n", number_of_active, SIZE);
    else fprintf(stderr, "Active spins: all / %d\n", SIZE);

  prepare_dipole_table(dipole_negligible);

  int size=SIZE*3;
  real* spins=(real*)malloc(sizeof(real)*size); assert(spins);
  if(initial_state) {
    copy_vector(size, initial_state, spins);
    free(initial_state);
  } else {
    random_vector(size, spins); 
  };
  skyrmion_steepest_descent(spins);
  // Saving result
  real energy[6]; skyrmion_energy(spins, energy);
  fprintf(stderr,COLOR_YELLOW"Zeeman energy:"COLOR_RESET" %.10"RF"f\n",energy[1]);
  fprintf(stderr,COLOR_YELLOW"Exchange energy:"COLOR_RESET" %.10"RF"f\n",energy[2]);
  fprintf(stderr,COLOR_YELLOW"Anisotropy energy:"COLOR_RESET" %.10"RF"f\n",energy[0]);
  fprintf(stderr,COLOR_YELLOW"Dzyaloshinskii-Moriya energy:"COLOR_RESET" %.10"RF"f\n",energy[3]);
  fprintf(stderr,COLOR_YELLOW"Dipole interaction energy:"COLOR_RESET" %.10"RF"f\n",energy[4]);
  real ergy=energy[5];
  fprintf(stderr,COLOR_YELLOW"TOTAL energy:"COLOR_RESET" %.10"RF"f\n",ergy);  

  // Checking gradient
  real* grad=(real*)malloc(sizeof(real)*size); assert(grad);
  projected_gradient(spins, grad, NULL);
  real nrm=rsqrt(normsq(size, grad));
  fprintf(stderr,COLOR_YELLOW"Gradient norm:"COLOR_RESET" %.10"RF"f\n",nrm/size);    
  free(grad);

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
      /*oct_save_linear(file);
      oct_save_state(file,"initial",path);
      oct_save_state(file,"final",path+3*SIZE*(sizep-1));
      real mx=energy[0]; int mi=0; for(int p=1; p<sizep; p++) if(energy[p]>mx) { mx=energy[p]; mi=p; }; 
      if(mi>0 && mi<sizep-1) oct_save_state(file,"saddle",path+3*SIZE*mi);
      oct_save_vertices(file);*/
      if(nonuniform_field) oct_save_field(file,"H",nonuniform_field);
      else oct_save_vector(file,"H",magnetic_field,3);
      int size[3]={sizex,sizey,sizez};
      oct_save_vector_int(file,"SZ",size,3);
      oct_save_vector_int(file,"BC",boundary_conditions,3);
      oct_save_matrix(file,"TRANSLATIONS",(real*)translation_vectors,3,3);
      oct_save_matrix(file,"CELL",(real*)atom_positions,sizeu,3);
      oct_save_matrix_int(file,"BONDS",(int*)neighbours,sizen,5);
      oct_save_real(file,"K",magnetic_anisotropy_norm);
      oct_save_real(file,"mu",dipole);
      oct_save_vector(file,"K0",magnetic_anisotropy_unit,3);
      oct_save_vector(file,"J",exchange_constant,sizen);
      oct_save_matrix(file,"D",dzyaloshinskii_moriya_vector,sizen,3);
      oct_save_path(file,"PATH",spins,1);
      oct_save_vector(file,"ENERGY",&ergy,1);
      if(save_octave>1) oct_save_hessian(file);
      oct_save_finish(file);
      fclose(file);
    } else 
      fprintf(stderr, COLOR_RED"Can not open '%s' for writing\n"COLOR_RESET,octfile);
  };
  // deinitializing
  free(spins);
  return 0;
};