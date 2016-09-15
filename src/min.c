#include "vector.h"
#include "skyrmion.h"
#include "optim.h"
#include "plot.h"

#include "debug.h"
#include "octave.h"
#include "integra.h"
#include "cmd.h"

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <getopt.h>

void skyrmion_display(int iter, real* restrict a, real* restrict grad_f, realp f, real res, 
real constres, real alpha, realp last_f, real last_res) {
  static realp prev_f=NAN; static real prev_res=NAN; 
  if(iter==1) { prev_f=f; prev_res=res; };
  if((iter%debug_every==0 && iter>0) || (debug_every==1) || iter<0) {
  	if(iter!=0) fprintf(stderr, "%6d", abs(iter));
    else fprintf(stderr, "%6s", "");
    fprintf(stderr, " "COLOR_YELLOW"E"COLOR_RESET);
    watch_number(f,debug_every==1?last_f:prev_f,16);
    fprintf(stderr, " "COLOR_YELLOW"R"COLOR_RESET);
    watch_number(res,debug_every==1?last_res:prev_res,16);
    fprintf(stderr, "%+"RF"g", constres);
    fprintf(stderr, " "COLOR_YELLOW"A"COLOR_RESET"%"RF"g", alpha);
    //fprintf(stderr, " "COLOR_YELLOW"d"COLOR_RESET"%.2g", pow(last_res/res,1./alpha));    
    fprintf(stderr, "\n");
  	if(debug_plot && iter!=0) 
      plot_field3(stdout,a);
      //plot_field3(stdout,grad_f);
    prev_f=f; prev_res=res;
  };
};

real quasynorm(real* restrict x) {
  return rsqrt(normalize(x));
};

void integrate(int N, void (*F)(const real* x, real* g, realp* E), real T, real* X, realp* E, int* iter) {
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
};

int skyrmion_steepest_descent(real* restrict x) {
/*
  return flow_descend(
    SIZE*3, (real*)x, 
    projected_gradient,
    mode, mode_param, epsilon, max_iter,
    skyrmion_display, 
    quasynorm,
    integrate
  );
*/
	return steepest_descend(
		SIZE*3, (real*)x, 
    projected_gradient,
		mode, mode_param, epsilon, max_iter,
		skyrmion_display, 
		quasynorm
	);
};

int main(int argc, char** argv) {
  int i=init_program(argc,argv,
    "Compute stable states.", "",
    NULL,NULL);
  if(i<argc) {
    fprintf(stderr, COLOR_RED"There are unused parameters:"COLOR_RESET"\n");
    while(i<argc) fprintf(stderr, "  %s\n", argv[i++]);
  };

  print_settings();

  int size=SIZE*3;
  real* spins=(real*)malloc(sizeof(real)*size); assert(spins);
  if(initial_state) {
    copy_vector(size, initial_state, spins);
    free(initial_state);
  } else {
    skyrmion_random(spins); 
  };

  skyrmion_steepest_descent(spins);

  // Saving result
  realp energy[6]; skyrmion_energy(spins, energy);
  fprintf(stderr,COLOR_YELLOW"Zeeman energy:"COLOR_RESET" %.10"RPF"f\n",energy[1]);
  fprintf(stderr,COLOR_YELLOW"Exchange energy:"COLOR_RESET" %.10"RPF"f\n",energy[2]);
  fprintf(stderr,COLOR_YELLOW"Anisotropy energy:"COLOR_RESET" %.10"RPF"f\n",energy[0]);
  fprintf(stderr,COLOR_YELLOW"Dzyaloshinskii-Moriya energy:"COLOR_RESET" %.10"RPF"f\n",energy[3]);
  fprintf(stderr,COLOR_YELLOW"Dipole interaction energy:"COLOR_RESET" %.10"RPF"f\n",energy[4]);
  realp ergy=energy[5];
  fprintf(stderr,COLOR_YELLOW"TOTAL energy:"COLOR_RESET" %.10"RPF"f\n",ergy);  

  // Checking gradient
  real* grad=(real*)malloc(sizeof(real)*size); assert(grad);
  projected_gradient(spins, grad, NULL);
  real nrm=rpsqrt(normsq(size, grad));
  fprintf(stderr,COLOR_YELLOW"Gradient norm:"COLOR_RESET" %.10"RF"f\n",nrm/size);    

  fprintf(stderr, COLOR_YELLOW"Saving gnuplot\n"COLOR_RESET);
  FILE* file=open_file(outdir,"/state.gnuplot",TRUE);
  if(file) {
    fprintf(file,"set terminal pngcairo\n");
    fprintf(file,"set output '%s/state.png'\n", outdir);
    plot_field3(file, spins);
    fclose(file);
  };
  // Octave output
  if(save_octave) {
    file=open_file(outdir, "/state.oct", TRUE);
    if(file) {
      oct_save_init(file);
      oct_save_lattice(file);
      oct_save_path(file,"PATH",spins,1);
      oct_save_vectorp(file,"ENERGY",&ergy,1);
      if(save_octave>1) oct_save_hessian(file);
      oct_save_path(file,"GRAD",grad,1);
      oct_save_finish(file);
      fclose(file);
    };
  };
  // deinitializing
  free(grad);
  free(spins);
  return 0;
};