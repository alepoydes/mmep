#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <getopt.h>
#include <string.h>

#include "parser.h"
#include "debug.h"
#include "skyrmion.h"
#include "octave.h"
#include "optim.h"

real epsilon=1e-9;
int long max_iter=10000;
real mode_param=0.2;
real param2=NAN;
int mode=SDM_PROGR;
int single_mode=0;
int integrator=0;
int debug_plot=0;
int debug_every=1000;
int save_octave=0;
real dipole_negligible=0.001;
real active_threshold=0;
int active_iterations=0;
int do_energy_shift=1;

const char* outdir=NULL;
const char default_outdir[]="tmp";

void showUsage(const char* program, const char* info, const char* options_desc) {
  fprintf(stderr, "%s\
\nUsage:\
\n    %s [options] [lattice description file]\
\nOptions:\
\n   -h      Show this message and exit\
\n   -p      Show graphics\
\n   -O      Save result in Octave format\
\n   -o      Save result in Octave format but Hessian matrix\
\n   -e REAL Desired gradient norm per atom\
\n   -E REAL Neglible value of dipole interaction\
\n   -i INT  Set maximum number of iterations\
\n   -r INT  Progress will be shown every given iteration\
\n   -I INT  Integrator\
\n   -m INT  Optimization method for MEP\
\n   -M INT  Optimization method for single image\
\n   -a REAL First parameter for optimization methods\
\n   -b REAL Second parameter for optimization methods\
\n   -D PATH Directory to output results\
\n   -t REAL Threshold of inactivity\
\n   -u INT  Update active spins every given interation\
\n   -S      Disable energy shift\
%s\
\n", info, program, options_desc);
};

int parseCommandLine(int argc, char** argv, 
const char* info, const char* options_desc,
const char* optstr, char (*handle)(char opt, const char* arg)
) {
  char buf[1024];
  const char common_opts[]="hpE:e:i:I:r:m:a:b:oOD:t:u:SM:";
  strncpy(buf, common_opts, sizeof(buf));
  if(optstr) strncat(buf, optstr, sizeof(buf)-strlen(buf));
  int c;
  while(1) {
    static struct option long_options[] = {      
      {"help", no_argument, 0, 'h'},
      {"debug", no_argument, 0, 1000},
      {0, 0, 0, 0}
    };
    int option_index = 0;
    c = getopt_long(argc, argv, buf, long_options, &option_index);
    if (c==-1) break;
    switch(c) {
      case 'p': debug_plot=1; break;
      case 'h': showUsage(argv[0], info, options_desc); exit(0);
      case 'e': epsilon=atof(optarg); break;
      case 'E': dipole_negligible=atof(optarg); break;      
      case 'i': max_iter=atoi(optarg); break;
      case 'I': integrator=atoi(optarg); break;
      case 'r': debug_every=atoi(optarg); break;
      case 'm': mode=atoi(optarg); break;
      case 'M': single_mode=atoi(optarg); break;
      case 'a': mode_param=atof(optarg); break;
      case 'b': param2=atof(optarg); break;
      case 'o': save_octave=1; break;      
      case 'O': save_octave=2; break;      
      case 'D': outdir=strdup(optarg); break;
      case 't': active_threshold=atof(optarg); break;
      case 'u': active_iterations=atoi(optarg); break;
      case 'S': do_energy_shift=0; break;
      case 1000: yydebug=1; break;
      case '?': break;
      default: 
        if(!handle || !handle(c, optarg)) {
      	  fprintf(stderr,COLOR_RED "Unprocessed option" COLOR_RESET " '%c'\n", c); 
      	  exit(1);
      	};
    };
  };
  return optind;
};

int init_program(int argc, char** argv,
const char* info, const char* options_desc,
const char* optstr, char (*handle)(char opt, const char* arg)
) {
  srand(time(NULL));
  init_signal();
  int i=parseCommandLine(argc,argv,info,options_desc,optstr,handle);
  if(i<argc) {
  	FILE* file=fopen(argv[i],"r");
  	if(!file) { fprintf(stderr, "Can not open file '%s'\n", argv[i]); exit(1); };
  	parse_lattice(file);
  	fclose(file);
  	i++;
  } else parse_lattice(stdin);

  prepare_energy_shift(do_energy_shift);
  prepare_dipole_table(dipole_negligible);

  if(!outdir) outdir=default_outdir;

  return i;
};

void print_settings() {
  fprintf(stderr, COLOR_BOLD "Machine epsilon:" COLOR_RESET " %"RF"g (%zd bytes per real)\n", RT(EPSILON), sizeof(real));
  fprintf(stderr,COLOR_BOLD "Using lattice %dx%dx%dx%d with %d bonds\n",sizeu,sizex,sizey,sizez,sizen);
  if(all_active) fprintf(stderr, COLOR_BOLD "Active spins:" COLOR_RESET " %d / %d\n", number_of_active, SIZE);
  else fprintf(stderr, COLOR_BOLD "Active spins:" COLOR_RESET " all / %d\n", SIZE);
  fprintf(stderr,COLOR_BOLD "Active spins gradient thr.:" COLOR_RESET " %"RF"g\n", RT(active_threshold));
  fprintf(stderr,COLOR_BOLD "Update active every" COLOR_RESET " %d-th iteration\n", active_iterations);
  fprintf(stderr,COLOR_BOLD "Energy shift" COLOR_RESET);
  for(int u=0; u<sizeu; u++) fprintf(stderr, " %"RF"g", RT(energy_shift_per_atom[u]));
  fprintf(stderr, " per atom\n");
  fprintf(stderr,COLOR_BOLD "Image optimizer" COLOR_RESET " %s\n", single_mode==0?"old":"new");  
};

void oct_save_lattice(FILE* file) {
	if(nonuniform_field) oct_save_field(file,"H",nonuniform_field);
  else oct_save_vector(file,"H",magnetic_field,3);
  int size[3]={sizex,sizey,sizez};
  oct_save_vector_int(file,"SZ",size,3);
  oct_save_vector_int(file,"BC",boundary_conditions,3);
  oct_save_matrix(file,"TRANSLATIONS",(real*)translation_vectors,3,3);
  oct_save_matrix(file,"CELL",(real*)atom_positions,sizeu,3);
  oct_save_matrix_int(file,"BONDS",(int*)neighbours,sizen,5);
  real K[magnetic_anisotropy_count];
  real K0[magnetic_anisotropy_count][3];
  for(int n=0; n<magnetic_anisotropy_count; n++) {
    K[n]=magnetic_anisotropy[n].norm;
    for3(j) K0[n][j]=magnetic_anisotropy[n].unit[j];
  };
  oct_save_vector(file,"K",K,magnetic_anisotropy_count);
  oct_save_matrix(file,"K0",(real*)K0,magnetic_anisotropy_count,3);
  oct_save_real(file,"mu",dipole);
  oct_save_vector(file,"J",exchange_constant,sizen);
  oct_save_matrix(file,"D",dzyaloshinskii_moriya_vector,sizen,3);
};