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
#include "display.h"

real epsilon=1e-6;
int max_iter=1000;
real mode_param=0.2;
real param2=1;
int mode=2;
int debug_plot=0;
int debug_every=1;

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
  lockDisplay(); copy_vector(3*SIZE,spins,display_buffer); releaseDisplay();
  // Main loop
  while(!is_aborting) { 
    
  };
  // Deinitialization
  deinitDisplay();
  free(spins);
  return 0;
};