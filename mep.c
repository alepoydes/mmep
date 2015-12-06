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
int sizep=0; // Number of nodes on path
int max_sizep=65; // Number of nodes on path
real epsilon=1e-6;
int max_iter=10000;
real mode_param=0.2;
int mode=2;
int debug_plot=0;
int debug_plot_path=0;
int debug_every=100;

real *distance=NULL, *energy=NULL, *diff=NULL, *tdiff=NULL;
int post_optimization=0;

void energy_evaluate(real* path) {
  real* q=(real*)malloc(sizeof(real)*size); assert(q);
  real* u=(real*)malloc(sizeof(real)*size); assert(u);
  for(int p=0; p<sizep; p++) {
    distance[p]=(p<=0)?0:distance[p-1]+rsqrt(distsq(size,path+size*(p-1),path+size*p));
    hamiltonian_hessian(path+size*p, q);
    energy[p]=-dot(size,path+size*p, q)/2;
    subtract_field(q);
    energy[p]+=dot(size,path+size*p, q);
    project_to_tangent(path+size*p,q);
    diff[p]=rsqrt(normsq(size,q)/size);
    tdiff[p]=diff[p];
    if(p>0 && p<sizep-1) {
      three_point_tangent(path+size*(p-1),path+size*(p+1),path+size*p,u);
      project_to_tangent(u,q);
      tdiff[p]=rsqrt(normsq(size,q)/size);
    };
  };
  free(u); free(q);
};

void skyrmion_display(int iter, real* restrict a, real* restrict grad_f, real f, real res, real constres, real alpha) {
  static real prev_f=NAN; 
  static real prev_res=NAN; 
  if(iter%(debug_every*sizep)==0 || iter<0) {
  	fprintf(stderr, "%d: E %"RF"g", abs(iter), f);
    if(f<prev_f) fprintf(stderr, COLOR_GREEN"%+"RF"g "COLOR_RESET, f-prev_f);
    else fprintf(stderr, COLOR_RED"%+"RF"g "COLOR_RESET, f-prev_f);
    if(res<prev_res) fprintf(stderr, "R "COLOR_GREEN"%"RF"g "COLOR_RESET, res);
    else fprintf(stderr, "R "COLOR_RED"%"RF"g "COLOR_RESET, res);
    fprintf(stderr, "A %"RF"g\n", alpha);

  	if(debug_plot && size<=1024) plot_field3(stdout,a);
    prev_f=f; prev_res=res;
  };
};

real quasynorm(real* restrict x) {
  normalize(x); return 0;
}


int skyrmion_steepest_descent(real* restrict x, int mode, real mode_param, 
	real epsilon, int max_iter) 
{
	return steepest_descend(
		size, (real*)x, 
		hamiltonian_hessian,	subtract_field,
		mode, mode_param, epsilon, max_iter,
		skyrmion_display, 
		quasynorm, project_to_tangent
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

void energy_display(FILE* file) {
  fprintf(file,"set ytics nomirror\nset y2tics nomirror\nset log y2\n");
  fprintf(file,"set autoscale xy\n");
  fprintf(file,"plot '-' using 1:2 with linespoints axes x1y1 title 'energy', '' using 1:3 with lines axes x1y2 title 'grad.', '' using 1:4 with lines axes x1y2 title 'orth. grad.'\n");
  for(int k=0; k<3; k++) {
    for(int p=0; p<sizep; p++)
      fprintf(file,"%.8"RF"e %.8"RF"e %.8"RF"e %.8"RF"e\n",distance[p],energy[p],diff[p],tdiff[p]);
    fprintf(file,"EOF\n\n");
  };
  fflush(file);
};

void path_display(int iter, real* restrict mep, real* restrict grad_f
, real f, real res, real restres, real alpha) {
  static real prev_f=NAN;
  static real prev_res=NAN;
  if(iter%debug_every==0 || iter<0) {
    fprintf(stderr, "%d: E %"RF"g", abs(iter), f);
    if(f<prev_f) fprintf(stderr, COLOR_GREEN"%+"RF"g "COLOR_RESET, f-prev_f);
    else fprintf(stderr, COLOR_RED"%+"RF"g "COLOR_RESET, f-prev_f);
    if(res<prev_res) fprintf(stderr, "R "COLOR_GREEN"%"RF"g "COLOR_RESET, res);
    else fprintf(stderr, "R "COLOR_RED"%"RF"g "COLOR_RESET, res);
    fprintf(stderr, "A %"RF"g\n", alpha);    
    if(debug_plot) { 
      if(debug_plot_path) plot_path(stdout, sizep, mep);
      else { /*energy_evaluate(mep);*/ energy_display(stdout); };
    };
    prev_f=f; prev_res=res;
  };
};

real path_normalize(real* mep) {
  for(int p=0; p<sizep; p++) normalize(mep+size*p);
  
  /*real* shifted=(real*)malloc(sizeof(real)*size*sizep); assert(shifted);
  for(int p=1; p<sizep-1; p++) {
    //if(energy[p]>=energy[p-1] && energy[p]>=energy[p+1]) {
      // local maxima
      // moving to the minimum along gradiend
    //} else if(energy[p]<energy[p-1] && energy[p]<energy[p+1]) {
      // local minimum
      // moving to the maxima
    //} else {
      // general point move orthogonal to path and equalizing distance between nodes
      three_point_equalizer(mep+size*(p-1), mep+size*p, mep+size*(p+1), shifted+size*p);
    //};
  };
  copy_vector(size*(sizep-2), shifted+size, mep+size);
  free(shifted);*/
 
  for(int p=1; p<sizep-1; p++) {
    if(post_optimization && energy[p]>=energy[p-1] && energy[p]>=energy[p+1]) {
    } else if(post_optimization && energy[p]<energy[p-1] && energy[p]<energy[p+1]) { 
    } else {
      three_point_equalize(mep+size*(p-1), mep+size*(p+1), mep+size*p);
    };
  };
  
  energy_evaluate(mep);
  /*real res=0; for(int p=0; p<sizep; p++) res+=tdiff[p]*tdiff[p]; res=rsqrt(res/sizep);
  if(res<10*epsilon) {
    post_optimization=1;
    fprintf(stderr, COLOR_YELLOW COLOR_BOLD"Post-optimization\n"COLOR_RESET);
  };*/
  return 0;
};

void path_tangent(const real* restrict mep, real* restrict grad) {
  for(int p=0; p<sizep; p++) {
    project_to_tangent(mep+size*p,grad+size*p);
    if(p>0 && p<sizep-1) {
      if(post_optimization && energy[p]>=energy[p-1] && energy[p]>=energy[p+1]) {
        three_point_reverse(mep+size*(p-1), mep+size*(p+1), grad+p*size);
      } else if(post_optimization && energy[p]<energy[p-1] && energy[p]<energy[p+1]) { 
      } else {
        three_point_project(mep+size*(p-1), mep+size*(p+1), grad+p*size);
      };
    };
  };
}

int path_steepest_descent(real* restrict path, int mode, 
  real mode_param, real epsilon, int max_iter) 
{
  real updated_param=mode_param;
  if(mode==2) updated_param=mode_param/sizep;
  return steepest_descend(
    size*sizep, (real*)path, 
    path_hessian, path_subtract_field,
    mode, updated_param, epsilon, max_iter,
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
\n   -p|--plot              Enable GNUPlot output of energy\
\n   -P|--plot-path         Enable GNUPlot output of MEP\
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
      {"plot-path", no_argument, 0, 'P'},
      {"epsilon", required_argument, 0, 'e'},
      {"mode", required_argument, 0, 'm'},
      {0, 0, 0, 0}
    };
    int option_index = 0;
    c = getopt_long(argc, argv, "hpPe:n:i:r:m:a:", long_options, &option_index);
    if (c==-1) break;
    switch(c) {
      case 'p': debug_plot=1; break;
      case 'P': debug_plot=1; debug_plot_path=1; break;
      case 'h': showUsage(argv[0]); exit(0);
      case 'e': epsilon=atof(optarg); break;
      case 'n': sizep=atoi(optarg); 
        max_sizep=3; while(max_sizep<sizep) max_sizep=2*(max_sizep-1)+1;
        break;
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
  if(max_sizep<2) {
    fprintf(stderr, "Number of nodes is too small: %d < 2\n", sizep);
    exit(1);
  };
  // Initializaton
  fprintf(stderr, "Size of real: %zd\n", sizeof(real));
  fprintf(stderr, "Nodes on path: %d\n", max_sizep);
  srand(time(NULL));
  size=sizeu*sizex*sizey*sizez*3;
  real* path=(real*)malloc(sizeof(real)*size*max_sizep); assert(path);
  distance=(real*)malloc(sizeof(real)*max_sizep); assert(distance);
  energy=(real*)malloc(sizeof(real)*max_sizep); assert(energy);
  diff=(real*)malloc(sizeof(real)*max_sizep); assert(diff);
  tdiff=(real*)malloc(sizeof(real)*max_sizep); assert(tdiff);
  // Set initla path size
  sizep=9; if(max_sizep<sizep) sizep=max_sizep;
  // find two minima
  fprintf(stderr, COLOR_YELLOW COLOR_BOLD"Calculating initial state\n"COLOR_RESET);
  if(initial_state) {
    copy_vector(size, initial_state, path);
    free(initial_state);
  } else random_vector(size, path); 
  skyrmion_steepest_descent(path, mode, mode_param, epsilon, max_iter);
  fprintf(stderr, COLOR_YELLOW COLOR_BOLD"Calculating final state\n"COLOR_RESET);
  if(final_state) {
    copy_vector(size, final_state, path+size*(sizep-1));
    free(final_state);
  } else set_to_field(path+size*(sizep-1));
  skyrmion_steepest_descent(path+size*(sizep-1), mode, mode_param, epsilon, max_iter);
  // Set initial path as geodesic approximation
  fprintf(stderr, COLOR_YELLOW COLOR_BOLD"Calculating MEP\n"COLOR_RESET);
  skyrmion_geodesic(sizep, path);
  // MEP calculation
  path_steepest_descent(path, mode, mode_param, epsilon, max_iter);
  while(2*(sizep-1)+1<=max_sizep) {
    // interpolating path
    for(int p=sizep-1; p>0; p--) copy_vector(size, path+p*size, path+2*p*size);
    for(int p=1; p<sizep; p++) 
      skyrmion_middle(path+2*size*(p-1), path+2*size*p, path+size*(2*p-1));
    sizep=2*(sizep-1)+1;
    fprintf(stderr, COLOR_YELLOW COLOR_BOLD"Increasing number of nodes: %d\n"COLOR_RESET, sizep);
    path_steepest_descent(path, mode, mode_param, epsilon, max_iter);
  };
  // Ouput result
  energy_evaluate(path);
  real max_energy=energy[0];
  for(int p=1; p<sizep; p++) if(max_energy<energy[p]) max_energy=energy[p];
  fprintf(stderr, COLOR_BLUE"Energy:"COLOR_RESET" initial %.8"RF"f maximum %.8"RF"f final %.8"RF"f\n", energy[0],max_energy,energy[sizep-1]);
  if(!debug_plot) {
    for(int p=0;p<sizep;p++) printf("%.8"RF"e ", energy[p]);
    printf("\n");
  };

  // save energy
  fprintf(stderr, COLOR_YELLOW COLOR_BOLD"Saving result\n"COLOR_RESET);
  const char* energyname="fields/energy.gnuplot";
  FILE* file=fopen(energyname,"w");
  if(file) {
    fprintf(file,"set terminal png\nset output 'fields/energy.png'\n");
    energy_display(file);
    fclose(file);
   } else {
    fprintf(stderr, COLOR_RED"Can not open '%s' for writing\n"COLOR_RESET,energyname);
  }; 
  // save field
  const char* filename="fields/mep.gnuplot";
  file=fopen(filename,"w");
  if(file) {
    fprintf(file,"set terminal gif animate delay 10\n");
    fprintf(file,"set output 'fields/mep.gif'\n");
    animate_path(file, sizep, path);
    fclose(file);
  } else {
    fprintf(stderr, COLOR_RED"Can not open '%s' for writing\n"COLOR_RESET,filename);
  };

  // Deinitialization
  free(path); free(tdiff); free(diff); free(energy); free(distance);
  return 0;
};