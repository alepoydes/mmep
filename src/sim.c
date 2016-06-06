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
#include "bitmap.h"

#define OUTDIR "fields"

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
int integrator=1; // 0 - RK, 1 - simplectic RK

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

// state
int iter=0;
real E=NAN;
real L=NAN;

void screen() {
  static real E0=NAN;
  static real sim_time0=0;
  static real L0=NAN;
  fprintf(stderr,RESET_CURSOR);
  fprintf(stderr,"%d: ", iter); 
  fprintf(stderr,"Time ");
  watch_number(sim_time,sim_time0,4);
  fprintf(stderr," Time step %"RF"g "COLOR_YELLOW"-/="COLOR_RESET" inc./dec.\n",time_step);

  fprintf(stderr,"Energy "); 
  watch_number(E,E0,6);
  fprintf(stderr," Damping %"RF"g ", damping);   
  fprintf(stderr,COLOR_YELLOW"[/]"COLOR_RESET" dec./inc. ");
  fprintf(stderr,COLOR_YELLOW"\\"COLOR_RESET" neg.    \n");

  fprintf(stderr,"Total length "); 
  watch_number(L,L0,10);
  fprintf(stderr,"           \n");  

  fprintf(stderr,"Field %s ",powered?"on":"off");
  fprintf(stderr,"Strength %"RF"g ", power);
  fprintf(stderr,COLOR_YELLOW";/'"COLOR_RESET" dec./inc.      \n");

  fprintf(stderr,"Integrator: ");
  switch(integrator) {
    case 0: fprintf(stderr,"Runge-Kutta"); break;
    case 1: fprintf(stderr,"simplectic Runge-Kutta"); break;
    default: fprintf(stderr,"unknown"); break;
  };
  fprintf(stderr," "COLOR_YELLOW"i"COLOR_RESET" chng.                      \n");  

  fprintf(stderr,"Key bindings:\n");
  fprintf(stderr,"  "COLOR_YELLOW"ENTER"COLOR_RESET" - make screenshot\n");  
  fprintf(stderr,"  "COLOR_YELLOW"SPACE"COLOR_RESET" - pause\n");
  fprintf(stderr,"  "COLOR_YELLOW"m"COLOR_RESET" - select spins/external field\n");
  fprintf(stderr,"  "COLOR_YELLOW"n"COLOR_RESET" - choose background\n");
  fprintf(stderr,"  "COLOR_YELLOW"b"COLOR_RESET" - show/hide boundung box\n");
  fprintf(stderr,"  "COLOR_YELLOW"r"COLOR_RESET" - reset camers\n");
  fprintf(stderr,"  "COLOR_YELLOW"v"COLOR_RESET" - cycle arrow styles\n");
  fprintf(stderr,"  "COLOR_YELLOW"t"COLOR_RESET" - turn transparency on/off\n");
  fprintf(stderr,"Mouse actions:\n");
  fprintf(stderr,"  "COLOR_YELLOW"wheel"COLOR_RESET" - change scale\n");
  fprintf(stderr,"  "COLOR_YELLOW"right drag"COLOR_RESET" - translate scene\n");
  fprintf(stderr,"  "COLOR_YELLOW"middle drag"COLOR_RESET" - rotate scene\n");
  fprintf(stderr,"  "COLOR_YELLOW"hold left"COLOR_RESET" - emit magnetic field\n");

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
  if(time_step>0) {
    // computing next state
    real delta=time_step;
    if(integrator==0) {
      runge_kutta(SIZE*3, dspins, delta, spins);
      normalize(spins);
    } else {
      //radau_integrator(SIZE*3, dspins, delta, spins, 12*SIZE*EPSILON, 10);
      radau_integrator(SIZE*3, dspins, delta, spins, 1e-7, 10);      
    };
    sim_time+=delta;
    iter++;
  };
  // check if energy conserves
  real* g=malloc(sizeof(real)*SIZE*3); assert(g);
  hamiltonian_hessian(spins, g);
  real E=-dot(3*SIZE,spins,g)/2;
  subtract_field(g);
  E+=dot(3*SIZE,spins,g);
  free(g);
  L=dot(3*SIZE,spins,spins);
  return E;
};

void do_print(int width, int height, void* buffer) {
  print_screen=NULL;
  fprintf(stderr, "Saving screenshot\n");
  FILE* file=fopen(OUTDIR"/screen.png","wb");
  write_png(file, width, height, (unsigned char*) buffer);
  fclose(file);
  fprintf(stderr, "Saved\n");
};

void keyboard_function(unsigned char key) {
  switch(key) {
    case '[': damping/=1.1; break; // SYNCHRONIZE !!!
    case ']': if(damping!=0) damping*=1.1; else damping=0.01; break;
    case '\\': damping=-damping; break;
    case '-': time_step/=1.1; break;
    case '=': if(time_step!=0) time_step*=1.1; else time_step=0.001; break;
    case '\'': power+=0.1; break;
    case ';': power-=0.1; break;
    case 13: print_screen=do_print; break;
    case ' ': time_step=0; break;
    case 'i': integrator=(integrator+1)%2; break;
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

  // initialize magnetic field
  nonuniform_field=malloc(sizeof(real)*SIZE*3); assert(nonuniform_field);
  set_to_field(nonuniform_field);
  // Main loop
  fprintf(stderr, CLEAR_SCREEN);
  display_buffer=malloc(sizeof(real)*SIZE*3); assert(display_buffer);
  while(!is_aborting) { 
    if(is_new_frame) {
      screen();
      lockDisplay(); 
      copy_vector(3*SIZE,spins,display_buffer); 
      releaseDisplay();
      displayRedraw();
    };
    E=doStep(spins);
  };
  // Deinitialization
  deinitDisplay();
  free(spins); free(nonuniform_field); free(display_buffer);
  return 0;
};