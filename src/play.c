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
real* energy=NULL;

int play_mode=0;
int speed=1;
int step=0;

void parseMEP(FILE* file) {
  int pos=0; int line=0;
  char buf[128]; 
  while(fgets(buf,sizeof(buf),file)) {
    //buf[127]=0;
    line++;
    if(sizef>=sizea) { // Allocating buffer
      sizea=sizea*2+1;
      mep=(real*)realloc(mep,sizeof(real)*SIZE*3*sizea);
    };
    long double pd[3],vd[3];
    int c=sscanf(buf,"%Lg %Lg %Lg %Lg %Lg %Lg",pd,pd+1,pd+2,vd,vd+1,vd+2);
    //real p[3]={pd[0],pd[1],pd[2]};
    real v[3]={vd[0],vd[1],vd[2]};
    if(c!=6) {
      if(pos==0) { // skipping head
        continue;
      };
      // end of frame
      if(pos!=SIZE) {
        fprintf(stderr,"Line %d: Frame %d is too short\n",line,sizef);
        break;
      };
      sizef++; pos=0; 
      fprintf(stderr, "line %d frame %d\n",line,sizef);
    } else {
      // next point of field
      for3(c) mep[3*(SIZE*sizef+pos)+c]=v[c]; pos++;
      if(pos>SIZE) {
        fprintf(stderr,"Line %d: Frame %d is too large\n",line,sizef);
        break;
      }
    };
  };
  if(pos!=SIZE) {
    fprintf(stderr,"Incomplete last frame\n");
  } else sizef++;
  if(sizef<1) {
    fprintf(stderr,"No frames read\n");
  };
};

void evaluateEnergy() {
  real* g=ralloc(SIZE*3); 
  for(int f=0; f<sizef; f++) {
    real* spins=mep+3*SIZE*f;
    hamiltonian_hessian(spins, g);
    energy[f]=-dot(3*SIZE,spins,g)/2;
    subtract_field(g);
    energy[f]+=dot(3*SIZE,spins,g);
    //energy[f]=NAN;
    fprintf(stderr, "Frame %d Energy %" RF "g\n",f,RT(energy[f]));
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

  energy=ralloc(sizef);
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
      fprintf(stderr, "Frame %d/%d Energy %" RF "g       \r", frame, sizef,RT(energy[frame]));
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
  free(mep); free(energy);
  return 0;
};