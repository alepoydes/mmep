#include "vector.h"
#include "skyrmion.h"
#include "optim.h"
#include "plot.h"
#include "debug.h"
#include "octave.h"
#include "cmd.h"
#include "magopt.h"

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <getopt.h>

int min_sizep=7; // Number of nodes on path
int max_sizep=25; // Number of nodes on path
real random_noise=0;
int do_not_relax_ends=0;

const char options_desc[]="\
\n   -P      Enable GNUPlot output of MEP\
\n   -n INT  Minimum number of images on MEP\
\n   -N INT  Maximum number of images on MEP\
\n   -z      Disable translations preserving energy\
\n   -R REAL Noise amplitude for initial path\
\n   -q      Skip ends relaxation\
\n";

char handle_option(char opt, const char* arg) {
  switch(opt) {
    case 'P': debug_plot=1; debug_plot_path=1; break;
    case 'n': min_sizep=atoi(optarg); break;
    case 'N': max_sizep=atoi(optarg); break;
    case 'z': remove_zero_modes=1; break;
    case 'R': random_noise=atof(optarg); break;
    case 'q': do_not_relax_ends=1; break;
    default: return FALSE;
  };
  return TRUE;
};

int main(int argc, char** argv) {
  mode=SDM_CONSTANT;
  int i=init_program(argc,argv,
    "Calculate MEP for magnetic systems.", options_desc,
    "zPN:n:R:q", handle_option);
  if(i<argc) {
    fprintf(stderr, COLOR_RED "There are unused parameters:" COLOR_RESET "\n");
    while(i<argc) fprintf(stderr, "  %s\n", argv[i++]);
  };

  // Initializaton
  print_settings();
  if(remove_zero_modes)
    fprintf(stderr, "Zero modes (translations) are removed\n");
  if(random_noise>0) 
    fprintf(stderr, "Initial path noise amplitude: %" RF "g\n", RT(random_noise));

  int size=3*SIZE; // Dimension of vector containing skyrmionic solutions
  real* path=ralloc(size*max_sizep); 
  path_steepest_descent_init(max_sizep);

  // Set initial path size
  sizep=min_sizep; 
  if(sizep<initial_states_count) sizep=initial_states_count;
  if(max_sizep<sizep) sizep=max_sizep;

  fprintf(stderr, "Images on path: %d in [%d, %d]\n", sizep, min_sizep, max_sizep);  
  // find two minima
  fprintf(stderr, COLOR_YELLOW COLOR_BOLD "Initializing path\n" COLOR_RESET);
  if(initial_states_count<2) {
    fprintf(stderr, COLOR_RED COLOR_BOLD "There should be at least two images in the path\n" COLOR_RESET);
    exit(1);
  };
  copy_vector(size, initial_state, path);
  int last_image=0;
  for (int n=1; n<initial_states_count; n++) {
    int next_image=n*sizep/(initial_states_count-1)-1;
    while(last_image>=next_image) next_image=last_image+1;
    assert(next_image<sizep); 
    copy_vector(size, initial_state+size*n, path+size*next_image);
    // Set initial path as geodesic approximation between given states
    skyrmion_geodesic(random_noise/(next_image-last_image), next_image-last_image+1, path+size*last_image);
    last_image=next_image;
  };
  assert(last_image==sizep-1);

  // minimize initial and final states
  
  if(!do_not_relax_ends) {
    fprintf(stderr, COLOR_YELLOW COLOR_BOLD "Relaxing initial state\n" COLOR_RESET);
    skyrmion_steepest_descent(path, SDM_PROGR, 0.1, epsilon, max_iter*sizep);
    fprintf(stderr, COLOR_YELLOW COLOR_BOLD "Relaxing final state\n" COLOR_RESET);
    skyrmion_steepest_descent(path+size*(sizep-1), SDM_PROGR, 0.1, epsilon, max_iter*sizep);
  };

  fprintf(stderr, COLOR_YELLOW COLOR_BOLD "Calculating MEP\n" COLOR_RESET);
  // MEP calculation
  post_optimization=0;
  path_steepest_descent(path, mode, mode_param, epsilon, max_iter);
  while(2*(sizep-1)+1<=max_sizep) {
    // interpolating path
    for(int p=sizep-1; p>0; p--) copy_vector(size, path+p*size, path+2*p*size);
    for(int p=1; p<sizep; p++) {
      if(p==1) skyrmion_middle_third_order(path+2*size*(p-1), path+2*size*p, path+2*size*(p+1), path+size*(2*p-1));
      else if(p==sizep-1) skyrmion_middle_third_order(path+2*size*p, path+2*size*(p-1), path+2*size*(p-2), path+size*(2*p-1));
      else skyrmion_middle_fourth_order(path+2*size*(p-2), path+2*size*(p-1), path+2*size*p, path+2*size*(p+1), path+size*(2*p-1));
      //skyrmion_middle(path+2*size*(p-1), path+2*size*p, path+size*(2*p-1));
    };
    sizep=2*(sizep-1)+1;
    fprintf(stderr, COLOR_YELLOW COLOR_BOLD "Increasing number of images to %d\n" COLOR_RESET, sizep);
    //post_optimization=0;    
    path_steepest_descent(path, mode, mode_param, epsilon, max_iter);
  };
  // Ouput result
  flat_distance=0;
  energy_evaluate(path);
  realp max_energy=energy[0];
  realp min_energy=energy[0];  
  for(int p=1; p<sizep; p++) {
    if(max_energy<energy[p]) max_energy=energy[p];
    if(min_energy>energy[p]) min_energy=energy[p];    
  };
  fprintf(stderr, COLOR_BLUE "Energy:" COLOR_RESET " initial %.8" RPF "f maximum %.8" RPF "f minimum %.8" RPF "f final %.8" RPF "f\n", RT(energy[0]),RT(max_energy),RT(min_energy),RT(energy[sizep-1]));
  fprintf(stderr, COLOR_BLUE "Barriers:" COLOR_RESET " forward %" RPF "g backward %" RPF "g\n", RT(max_energy-energy[0]), RT(max_energy-energy[sizep-1]));  
  if(!debug_plot) {
    for(int p=0;p<sizep;p++) printf("%.8" RPF "g ", RT(energy[p]));
    printf("\n");
    for(int p=0;p<sizep;p++) printf("%.8" RPF "g ", RT(distance[p]));
    printf("\n");
  };
  // Compare distances between images
  real mean=distance[sizep-1]/(sizep-1); 
  real var=0; for(int p=1; p<sizep; p++) {
    real d=distance[p]-distance[p-1]-mean;
    var+=d*d;
  };  
  fprintf(stderr, COLOR_BLUE COLOR_BOLD "Std of distances between images:" COLOR_RESET " %" RF "g\n", RT(rsqrt(var)));

  // save energy
  fprintf(stderr, COLOR_YELLOW COLOR_BOLD "Saving result\n" COLOR_RESET);
  FILE* file=open_file(outdir,"/energy.gnuplot",TRUE);
  if(file) {
    fprintf(file,"set terminal png\nset output '%s/energy.png'\n", outdir);
    energy_display(file);
    fclose(file);
  }; 
  // save field
  file=open_file(outdir, "/mep.gnuplot", TRUE);
  if(file) {
    fprintf(file,"set terminal gif animate delay 10\n");
    fprintf(file,"set output '%s/mep.gif'\n",outdir);
    animate_path(file, sizep, path);
    fclose(file);
  };
  // Octave output
  if(save_octave>0) {
    fprintf(stderr, COLOR_YELLOW "Computing energy contributions\n" COLOR_RESET);
    realp* contr=(realp*)malloc(sizeof(realp)*sizep*6);
    for(int p=0; p<sizep; p++) skyrmion_energy(path+p*SIZE*3, contr+6*p);
    file=open_file(outdir, "/mep.oct", TRUE);    
    if(file) {
      oct_save_init(file);
      oct_save_lattice(file);
      oct_save_path(file,"PATH",path,sizep);
      oct_save_vectorp(file,"ENERGY",energy,sizep);
      oct_save_vectorp(file,"DISTANCE",distance,sizep);
      oct_save_vectorp(file,"DIFF",diff,sizep);
      oct_save_vectorp(file,"TDIFF",tdiff,sizep);
      oct_save_matrixp(file,"CONTRIBUTIONS",contr,sizep,6);
      if(save_octave>1) oct_save_hessian(file);
      oct_save_finish(file);
      fclose(file);
    };
  };
  // Deinitialization
  path_steepest_descent_deinit();
  free(path); 
  return 0;
};