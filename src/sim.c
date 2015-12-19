#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <getopt.h>

#include <GL/freeglut.h>

#include "vector.h"
#include "skyrmion.h"
#include "optim.h"
#include "plot.h"
#include "parse.h"
#include "debug.h"
#include "display.h"
#include "integra.h"

#define OUTDIR "../fields"

real epsilon=1e-6;
int max_iter=1000;
real mode_param=0.2;
real param2=1;
int mode=2;
int debug_plot=0;
int debug_every=1;
real time_step=0;
real sim_time=0.0;
real damping=0;
int powered=0;
real power=-1;

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
    c = getopt_long(argc, argv, "hpe:i:r:m:a:b:", long_options, &option_index);
    if (c==-1) break;
    switch(c) {
      case 'p': debug_plot=1; break;
      case 'h': showUsage(argv[0]); exit(0);
      case 'e': epsilon=atof(optarg); break;
      case 'i': max_iter=atoi(optarg); break;
      case 'r': debug_every=atoi(optarg); break;
      case 'm': mode=atoi(optarg); break;
      case 'a': mode_param=atof(optarg); break;
      case 'b': param2=atof(optarg); break;
      case '?': break;
      default: fprintf(stderr,"Unprocessed option '%c'\n", c); exit(1);
    };
  };
  return optind;
};

void dspins(real* X, real* G) {
  hamiltonian_hessian(X, G);
  subtract_field(G);
  forall(u,x,y,z) {
    int i=INDEX(u,x,y,z)*3;
    real vec[3]; cross3(X+i,G+i,vec);
    real vec2[3]; cross3(vec,X+i,vec2);
    for3(c) G[i+c]=vec[c]-damping*vec2[c];
  };
};

real doStep(real* spins) {
  // computing next state
  real delta=time_step;
  runge_kutta(SIZE*3, dspins, delta, spins);
  normalize(spins);
  sim_time+=delta;
  // check if energy conserves
  real* g=malloc(sizeof(real)*SIZE*3); assert(g);
  hamiltonian_hessian(spins, g);
  real E=-dot(3*SIZE,spins,g)/2;
  subtract_field(g);
  E+=dot(3*SIZE,spins,g);
  free(g);
  return E;
};

void keyboard_function(unsigned char key) {
  switch(key) {
    case '[': damping/=1.1; break; // SYNCHRONIZE !!!
    case ']': if(damping>0) damping*=1.1; else damping=0.01; break;
    case '\\': damping=-damping; break;
    case '-': time_step/=1.1; break;
    case '=': if(time_step>0) time_step*=1.1; else time_step=0.001; break;
    case '\'': power+=0.1; break;
    case ';': power-=0.1; break;
  }
};

void mouse_function(int button, int state, real p[3]) {
  if(state==GLUT_DOWN) {
    powered=1;
  } else {
    powered=0;
    set_to_field(nonuniform_field);
  };
};

void motion_function(real p[3]) {
  if(powered) {
    //fprintf(stderr, "%"RF"g %"RF"g %"RF"g \n", p[0], p[1], p[2]);
    forall(u,x,y,z) {
      int i=INDEX(u,x,y,z)*3;
      real vec[3]; COORDS(u,x,y,z,vec);
      vec[0]-=p[0]; vec[1]-=p[1]; 
      real dist=rsqrt(vec[0]*vec[0]+vec[1]*vec[1]);
      dist=1+power/(1+dist);
      for3(c) nonuniform_field[i+c]=magnetic_field[c]*dist;
    };
  };
};

int main(int argc, char** argv) {
  srand(time(NULL));
  // Initialize lattice
  int i=parseCommandLine(argc,argv);
  if(i<argc) {
  	FILE* file=fopen(argv[i],"r");
  	if(!file) { fprintf(stderr, "Can not open file '%s'\n", argv[i]); exit(1); };
  	parse_lattice(file);
  	fclose(file);
  } else parse_lattice(stdin);
  // Initialize states
  real* spins=(real*)malloc(sizeof(real)*SIZE*3); assert(spins);
  if(initial_state) {
    copy_vector(SIZE*3, initial_state, spins);
    free(initial_state);
  } else random_vector(SIZE*3, spins); 
  // Initialize graphics
  is_aborting=initDisplay(&argc, argv);
  fprintf(stderr,"Mouse actions:\n");
  fprintf(stderr,"  wheel - change scale\n");
  fprintf(stderr,"  right drag - translate scene\n");
  fprintf(stderr,"  middle drag - rotate scene\n");
  fprintf(stderr,"  hold left - emit magnetic field\n");
  fprintf(stderr,"Key bindings:\n");
  fprintf(stderr,"  b - show/hide boundung box\n");
  fprintf(stderr,"  r - reset camers\n");
  fprintf(stderr,"  v - cycle arrow styles\n");
  fprintf(stderr,"  t - turn transparency on/off\n");
  fprintf(stderr,"  m - select spins/external field\n");
  fprintf(stderr,"  n - choose background\n");
  fprintf(stderr,"  -/= - decrease/increase time speed\n");
  fprintf(stderr,"  [/] - decrease/increase damping\n");
  fprintf(stderr,"  ;/' - decrease/increase emited field strength\n");

  // initialize magnetic field
  nonuniform_field=malloc(sizeof(real)*SIZE*3); assert(nonuniform_field);
  set_to_field(nonuniform_field);
  // Main loop
  int iter=0;
  real E=NAN;
  while(!is_aborting) { 
    if(is_new_frame) {
      fprintf(stderr, "%d: T=%"RF"g E=%"RF"g dT=%"RF"g dE=%"RF"g P=%"RF"g        \r", iter, sim_time, E, time_step, damping,power);    
      lockDisplay(); 
      copy_vector(3*SIZE,spins,display_buffer); 
      releaseDisplay();
      displayRedraw();
    };
    E=doStep(spins);
    iter++;
  };
  // Deinitialization
  deinitDisplay();
  free(spins); free(nonuniform_field);
  return 0;
};