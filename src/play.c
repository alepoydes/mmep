#include <unistd.h>
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

#define OUTDIR "fields"

int frame=0; // current frame
int sizef=0; // Number of loaded frames
int sizea=0; // Number of frames allocated in memory
real* mep=NULL;

int play_mode=0;
int speed=1;
int step=0;

void showUsage(const char* program) {
  fprintf(stderr, "Render MEP interactively.\
\nUsage:\
\n    %s [options] <lattice description file> <MEP file>\
\nMEP file is normally saved by mepx to "OUTDIR"/mep.hnuplot\
\nOptions:\
\n   -h|--help              Show this message and exit\
\n", program);
};

int parseCommandLine(int argc, char** argv) {
  int c;
  while(1) {
    static struct option long_options[] = {      
      {"help", no_argument, 0, 'h'},
      {0, 0, 0, 0}
    };
    int option_index = 0;
    c = getopt_long(argc, argv, "h", long_options, &option_index);
    if (c==-1) break;
    switch(c) {
      case 'h': showUsage(argv[0]); exit(0);
      case '?': break;
      default: fprintf(stderr,"Unprocessed option '%c'\n", c); exit(1);
    };
  };
  return optind;
};

void parseMEP(FILE* file) {
  int pos=0; int line=0;
  char buf[128]; 
  while(fgets(buf,sizeof(buf),file)) {
    //buf[127]=0;
    line++;
    if(sizef>=sizea) { // Allocating buffer
      sizea=sizea*2+1;
      mep=realloc(mep,sizeof(real)*SIZE*3*sizea);
    };
    real p[3],v[3];
    int c=sscanf(buf,"%"RF"g %"RF"g %"RF"g %"RF"g %"RF"g %"RF"g",p,p+1,p+2,v,v+1,v+2);
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

void switchFrame() {
  display_buffer=mep+frame*SIZE*3;
  displayRedraw();
};

void keyboard_function(unsigned char key) {
  switch(key) {
    case '-': frame=(sizef+frame-1)%sizef; switchFrame(); break;
    case '=': frame=(frame+1)%sizef; switchFrame(); break;
    case ']': if(abs(speed)>1) speed=speed>0?speed-1:speed+1; break;
    case '[': speed=speed>0?speed+1:speed-1; break;
    case ' ': play_mode=(play_mode+1)%3; break;
  }
};

void mouse_function(int button, int state, real p[3]) {
  //if(state==GLUT_DOWN) { } else { };
};

void motion_function(real p[3]) {
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
  } else { fprintf(stderr, "No lattice file is provided\n"); exit(1); };
  i++;
  // Read MEP
  if(i<argc) {
    FILE* file=fopen(argv[i],"r");
    if(!file) { fprintf(stderr, "Can not open file '%s'\n", argv[i]); exit(1); };
    parseMEP(file);
    fclose(file);
  } else { fprintf(stderr, "No MEP file is provided\n"); exit(1); };
  assert(mep);
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

  // Main loop
  while(!is_aborting) { 
    step++;
    usleep(100);
    if(is_new_frame) {
      fprintf(stderr, "Frame %d/%d        \r", frame, sizef);
      if(play_mode==1) {
        if(step%abs(speed)==0) {
          frame=(frame+1)%sizef;
          switchFrame();
        };
      } else if(play_mode==2) {
        if(step%abs(speed)==0) {
          if(speed>0) {
            frame++; if(frame>=sizen) { frame=sizen-1; speed=-speed; };
          } else {
            frame--; if(frame<0) { frame=0; speed=-speed; };
          };
          frame=(frame+1)%sizef;
          switchFrame();
        };
      };
    };
  };
  // Deinitialization
  deinitDisplay();
  free(mep); 
  return 0;
};