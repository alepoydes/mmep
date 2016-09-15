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
#include "cmd.h"

real time_step=0;
real sim_time=0.0;
real damping=0;
int powered=0;
real power=-1;

// state
int iter=0;
realp E=NAN;
real L=NAN;

void screen() {
  static real E0=NAN;
  static real sim_time0=0;
  static real L0=NAN;
  fprintf(stderr,RESET_CURSOR);
  fprintf(stderr,"%d: ", iter); 
  fprintf(stderr,"Time ");
  watch_number(sim_time,sim_time0,4); sim_time0=sim_time;
  fprintf(stderr," Time step %"RF"g ",time_step);
  fprintf(stderr,COLOR_YELLOW"-/="COLOR_RESET" inc./dec. ");
  fprintf(stderr,COLOR_YELLOW"`123"COLOR_RESET" speed \n");

  fprintf(stderr,"Energy "); 
  watch_number(E,E0,6); E0=E;
  fprintf(stderr," Damping %"RF"g ", damping);   
  fprintf(stderr,COLOR_YELLOW"[/]"COLOR_RESET" dec./inc. ");
  fprintf(stderr,COLOR_YELLOW"\\"COLOR_RESET" neg.    \n");

  fprintf(stderr,"Total length "); 
  watch_number(L,L0,6); L0=L;
  fprintf(stderr,"           \n");  

  fprintf(stderr,"Field %s ",powered?"on":"off");
  fprintf(stderr,"Strength %"RF"g ", power);
  fprintf(stderr,COLOR_YELLOW";/'"COLOR_RESET" dec./inc.      \n");

  fprintf(stderr,"Integrator: ");
  switch(integrator) {
    case 0: fprintf(stderr,"Runge-Kutta"); break;
    case 1: fprintf(stderr,"Gauss-Legendre Runge-Kutta"); break;
    case 2: fprintf(stderr,"Radau"); break;
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
  fprintf(stderr,"  "COLOR_YELLOW"middle drag"COLOR_RESET" - translate scene\n");
  fprintf(stderr,"  "COLOR_YELLOW"right drag"COLOR_RESET" - rotate scene\n");
  fprintf(stderr,"  "COLOR_YELLOW"hold left"COLOR_RESET" - emit magnetic field\n");

  fprintf(stderr,"\nDipole int. neigh. %d\n", dipole_count);

};

void dspins(const real* X, real* G, realp* E) {
  skyrmion_gradient(X, G, E);
  forall(u,x,y,z) {
    int i=INDEX(u,x,y,z)*3;
    real vec[3]; cross3(X+i,G+i,vec);
    real vec2[3]; cross3(vec,X+i,vec2);
    for3(c) G[i+c]=vec[c]-damping*vec2[c];
  };
};

void doStep(real* spins) {
  if(time_step>0) {
    // computing next state
    real delta=time_step;
    switch(integrator) {
      case 0:
        runge_kutta(SIZE*3, dspins, delta, spins, &E, NULL);
        normalize(spins);
        break;
      case 1:
        gauss_integrator(SIZE*3, dspins, delta, spins, 1e-7, 10, &E, NULL);
        break;
      case 2: 
        radau_integrator(SIZE*3, dspins, delta, spins, 1e-7, 10, &E, NULL);
        break;
    };
    sim_time+=delta;
    iter++;
  // check if energy conserves
    /*real* g=malloc(sizeof(real)*SIZE*3); assert(g);
    hamiltonian_hessian(spins, g);
    E=-dot(3*SIZE,spins,g)/2;
    subtract_field(g);
    E+=dot(3*SIZE,spins,g);
    free(g);*/
    L=dot(3*SIZE,spins,spins);
  };
};

void do_print(int width, int height, void* buffer) {
  print_screen=NULL;
  fprintf(stderr, "Saving screenshot\n");
  FILE* file=open_file(outdir, "/screen.png", TRUE);
  if(file) {
    write_png(file, width, height, (unsigned char*) buffer);
    fclose(file);
    fprintf(stderr, "Saved\n");
  };
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
    case 'i': integrator=(integrator+1)%3; break;
    case '`': time_step=0; break;
    case '1': time_step=0.001; break;
    case '2': time_step=0.01; break;
    case '3': time_step=0.1; break;
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

const char options_desc[]="\
\n";

char handle_option(char opt, const char* arg) {
  switch(opt) {
    default: return FALSE;
  };
  return TRUE;
};

int main(int argc, char** argv) {
  epsilon=1e-6;
  integrator=1; // 0 - RK, 1 - simplectic RK

  int i=init_program(argc,argv,
    "Integrate LLG equation and visializate dynamics.", options_desc,
    "", handle_option);
  if(i<argc) {
    fprintf(stderr, COLOR_RED"There are unused parameters:"COLOR_RESET"\n");
    while(i<argc) fprintf(stderr, "  %s\n", argv[i++]);
  };

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
    doStep(spins);
  };
  // Deinitialization
  deinitDisplay();
  free(spins); free(nonuniform_field); free(display_buffer);
  return 0;
};