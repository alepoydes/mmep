#include "vector.h"
#include "skyrmion.h"
#include "optim.h"
#include "plot.h"
#include "debug.h"
#include "display.h"
#include "bitmap.h"
#include "cmd.h"

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <getopt.h>

#include <GL/freeglut.h>

int frame=0; // current frame
int sizef=0; // Number of loaded frames
int sizea=0; // Number of frames allocated in memory
real* mep=NULL;
real* ergy=NULL;

int play_mode=0;
int speed=1;
int step=0;

real* allocate_image_on_mep() {
  if(sizef>=sizea) { // Allocating buffer
      sizea=sizea*2+1;
      mep=(real*)realloc(mep,sizeof(real)*SIZE*3*sizea);
  };
  return mep+3*SIZE*(sizef++);
};

void parseMEP(FILE* file) {
  int count=load_path_from_gnuplot(file, allocate_image_on_mep);
  assert(count==sizef); 
};

void evaluateEnergy() {
  real* g=ralloc(SIZE*3); 
  for(int f=0; f<sizef; f++) {
    real* spins=mep+3*SIZE*f;
    hamiltonian_hessian(spins, g);
    ergy[f]=-dot(3*SIZE,spins,g)/2;
    subtract_field(g);
    ergy[f]+=dot(3*SIZE,spins,g);
    //ergy[f]=NAN;
    fprintf(stderr, "Frame %d Energy %" RF "g\n",f,RT(ergy[f]));
  };
  free(g);
};

void switchFrame() {
  assert(frame>=0); assert(frame<sizef); 
  lockDisplay(); display_buffer=mep+frame*SIZE*3; releaseDisplay();
  displayRedraw();
};

void do_print(int width, int height, void* buffer) {
  print_screen=NULL;
  fprintf(stderr, "Saving screenshot                 \r");
  FILE* file=open_file(outdir,"/screen.png",TRUE);
  if(!file) exit(1); 
  write_png(file, width, height, (unsigned char*) buffer);
  fclose(file);
  fprintf(stderr, "Saved\r");
};

void do_animation(int width, int height, void* buffer) {
  char name[128]; sprintf(name,"/%04d.png",frame);
  FILE* file=open_file(outdir,name,TRUE);
  if(!file) exit(1);
  write_png(file, width, height, (unsigned char*) buffer);
  fclose(file);
  frame++; 
  if(frame>=sizef) {
    print_screen=NULL;
    fprintf(stderr, "Saved                         \r");
    play_mode=-play_mode; frame=0;
    return;
  };
  switchFrame();
};

void keyboard_function(unsigned char key) {
  switch(key) {
    case '-': frame=(sizef+frame-1)%sizef; switchFrame(); break;
    case '=': frame=(frame+1)%sizef; switchFrame(); break;
    case ']': if(abs(speed)>1) speed=speed>0?speed-1:speed+1; break;
    case '[': speed=speed>0?speed+1:speed-1; break;
    case ' ': play_mode=(play_mode+1)%3; break;
    case 13: lockDisplay(); print_screen=do_print; releaseDisplay(); break;
    case 's': play_mode=-play_mode; frame=0; switchFrame();
      lockDisplay(); print_screen=do_animation; releaseDisplay();
      break;
  }
};

void mouse_function(int button, int state, real p[3]) {
  //if(state==GLUT_DOWN) { } else { };
};

void motion_function(real p[3]) {
};


int main(int argc, char** argv) {
  int i=init_program(argc,argv,
    "Render MEP interactively.", "",
    NULL, NULL);
  // Read MEP
  if(i<argc) {
    FILE* file=open_file(NULL,argv[i],FALSE);
    if(!file) exit(1);
    parseMEP(file);
    fclose(file); 
    i++;
  } else { 
    fprintf(stderr, "No MEP file is provided\n"); 
    exit(1); 
  };
  assert(mep);
  if(i<argc) {
    fprintf(stderr, COLOR_RED "There are unused parameters:" COLOR_RESET "\n");
    while(i<argc) fprintf(stderr, "  %s\n", argv[i++]);
  };

  ergy=ralloc(sizef);
  evaluateEnergy();
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
  fprintf(stderr,"  -/= - previous/next frame\n");
  fprintf(stderr,"  Space - autoplay mode\n");
  fprintf(stderr,"  Enter - save screen\n");
  fprintf(stderr,"  s - save all frames\n");

  // Main loop
  display_buffer=mep;
  while(!is_aborting) { 
    usleep(100);
    if(is_new_frame) {
      fprintf(stderr, "Frame %d/%d Energy %" RF "g       \r", frame, sizef,RT(ergy[frame]));
      if(play_mode==1) {
        if(step%abs(speed)==0) {
          frame=(frame+1)%sizef;
          switchFrame();
        };
      } else if(play_mode==2) {
        if(step%abs(speed)==0) {
          if(speed>0) {
            frame++; if(frame>=sizef) { frame=sizef-1; speed=-speed; };
          } else {
            frame--; if(frame<0) { frame=0; speed=-speed; };
          };
          switchFrame();
        };
      };
    };
    step++;
  };
  // Deinitialization
  deinitDisplay();
  free(mep); free(ergy);
  return 0;
};