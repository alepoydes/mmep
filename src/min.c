#include "vector.h"
#include "skyrmion.h"
#include "optim.h"
#include "plot.h"
#include "debug.h"
#include "octave.h"
#include "integra.h"
#include "cmd.h"
#include "magopt.h"

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <getopt.h>

int main(int argc, char** argv) {
  int i=init_program(argc,argv,
    "Compute stable states.", "",
    NULL,NULL);
  if(i<argc) {
    fprintf(stderr, COLOR_RED "There are unused parameters:" COLOR_RESET "\n");
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

  skyrmion_steepest_descent(spins, mode, mode_param, epsilon, max_iter);

  // Saving result
  realp energy[6]; skyrmion_energy(spins, energy);
  fprintf(stderr,COLOR_YELLOW "Zeeman energy:" COLOR_RESET " %.10" RPF "f\n",RT(energy[1]));
  fprintf(stderr,COLOR_YELLOW "Exchange energy:" COLOR_RESET " %.10" RPF "f\n",RT(energy[2]));
  fprintf(stderr,COLOR_YELLOW "Anisotropy energy:" COLOR_RESET " %.10" RPF "f\n",RT(energy[0]));
  fprintf(stderr,COLOR_YELLOW "Dzyaloshinskii-Moriya energy:" COLOR_RESET " %.10" RPF "f\n",RT(energy[3]));
  fprintf(stderr,COLOR_YELLOW "Dipole interaction energy:" COLOR_RESET " %.10" RPF "f\n",RT(energy[4]));
  realp ergy=energy[5];
  fprintf(stderr,COLOR_YELLOW "TOTAL energy:" COLOR_RESET " %.10" RPF "f\n",RT(ergy));  

  // Checking gradient
  real* grad=(real*)malloc(sizeof(real)*size); assert(grad);
  projected_gradient(spins, grad, NULL);
  real nrm=rpsqrt(normsq(size, grad));
  fprintf(stderr,COLOR_YELLOW "Gradient norm:" COLOR_RESET " %.10" RF "f\n",RT(nrm/size));    

  fprintf(stderr, COLOR_YELLOW "Saving gnuplot\n" COLOR_RESET);
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