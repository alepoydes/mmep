#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <getopt.h>
#include <unistd.h>
#include <alloca.h>
#include <string.h>

#include "vector.h"
#include "skyrmion.h"
#include "parse.h"
#include "debug.h"
#include "integra.h"

#define OUTDIR "fields"

long int max_iter=5000000;
real step=0;
real dipole_negligible=0.001;
real debug_every_sec=10;
int update_by_one=1;
int rejection_rule=0;
int random_number_generator=1;
char filename[1024]="";

#define WRITE_ANISOTROPY 1
#define WRITE_ZEEMAN 2
#define WRITE_EXCHANGE 4
#define WRITE_DMI 8
#define WRITE_DIPOLAR 16
#define WRITE_TOTAL 32
int what_to_write=WRITE_TOTAL | WRITE_DMI;

void showUsage(const char* program) {
  fprintf(stderr, "Compute stable states.\
\nUsage:\
\n    %s [options] [lattice description file]\
\nOptions:\
\n   -h|--help              Show this message and exit\
\n   -s           REAL      Monte-Carlo step size\
\n   -T           REAL      kT=Temperature by Bolzman factor\
\n   -E           REAL      Neglible value of dipole interaction\
\n   -R           0..1      Rejection rule\
\n   -N           0..1      Noise type\
\n   -i           INT       Set maximum number of steps\
\n   -r           INT       Delay between reports\
\n   -u                     Switch simulaneous/one by one updates\
\n   -o           FILE      Where to write result\
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
    c = getopt_long(argc, argv, "uhs:i:E:T:r:R:N:o:", long_options, &option_index);
    if (c==-1) break;
    switch(c) {
      case 'u': update_by_one=!update_by_one; break;
      case 'h': showUsage(argv[0]); exit(0);
      case 's': step=atof(optarg); break;
      case 'T': temperature=atof(optarg); break;            
      case 'E': dipole_negligible=atof(optarg); break;
      case 'N': random_number_generator=atoi(optarg); break;
      case 'i': max_iter=atof(optarg); break;
      case 'r': debug_every_sec=atof(optarg); break;
      case 'R': rejection_rule=atoi(optarg); break;
      case 'o': strncpy(filename,optarg,sizeof(filename)); break;
      case '?': break;
      default: fprintf(stderr,"Unprocessed option '%c'\n", c); exit(1);
    };
  };
  return optind;
};

int is_to_reject(real denergy) {
    // exp(denergy) = P_new/P_current
    // P_new>P_current <=> exp(denergy)>1 <=> denergy>0
    // check if to reject state
    // state is accepted if
    // 1. definitely if P_new>P_current
    // 2. with probability P_new/P_current
    // If alpha is uniformly distributed in [0,1].
    // then accept if alpha<P_new/P_old.  
    denergy/=temperature;
    switch(rejection_rule) {
      case 0: return rexp(denergy)<random_real();
      case 1: return  1./(1.+rexp(denergy))>random_real();
      default: fprintf(stderr, "Unknown rejection rule: %d\n", rejection_rule); exit(1);
    };
};

int main(int argc, char** argv) {
  init_signal();
  int i=parseCommandLine(argc,argv);
  if(i<argc) {
  	FILE* file=fopen(argv[i],"r");
  	if(!file) { fprintf(stderr, "Can not open file '%s'\n", argv[i]); exit(1); };
  	parse_lattice(file);
  	fclose(file);
  } else parse_lattice(stdin);
  if(active) fprintf(stderr, COLOR_YELLOW"Active spins:"COLOR_RESET" %d / %d\n", number_of_active, SIZE);
    else fprintf(stderr, COLOR_YELLOW"Active spins:"COLOR_RESET" all / %d\n", SIZE);
  if(temperature==0) {
    fprintf(stderr, COLOR_RED"Temperature is not set or zero\n");
    exit(1);
  };
  fprintf(stderr, COLOR_YELLOW"Temperature:"COLOR_RESET" %"RF"g\n", temperature); 
  fprintf(stderr, COLOR_YELLOW"Step size:"COLOR_RESET" %"RF"g\n", step); 
  fprintf(stderr, COLOR_YELLOW"Number of steps:"COLOR_RESET" %ld\n", max_iter); 
  fprintf(stderr, COLOR_YELLOW"Update type:"COLOR_RESET" %s\n", update_by_one?"one by one":"simultaneous"); 
  fprintf(stderr, COLOR_YELLOW"Rejection rule:"COLOR_RESET" %s\n", rejection_rule==0?"Metropolis-Hastings":"Nowak");   
  fprintf(stderr, COLOR_YELLOW"Noise type:"COLOR_RESET" %s\n", random_number_generator==0?"fast noise":"symmetric noise"); 
  // precompute some data
  prepare_dipole_table(dipole_negligible);
  // Initialize timer
  time_t last_time_reported;
  srand(time(&last_time_reported));
  time_t start_time=last_time_reported;   
  // open file for writing
  if(filename[0]==0) {
    struct tm* timeinfo=localtime(&start_time);
    char buf[100]; strftime(buf, sizeof(buf), "%H:%M:%S", timeinfo);
    pid_t pid=getpid();
    sprintf(filename, "fields/mc.%s.%d.bin", buf, pid);
  };
  FILE* file=fopen(filename,"wb");
  if(!file) {
    fprintf(stderr, COLOR_RED"Can not open '%s' for writing\n"COLOR_RESET,filename);
    exit(1);
  };
  fprintf(stderr, "Writing to file '"COLOR_BLUE"%s"COLOR_RESET"'\n", filename);
  // compute parameters
  if(step<=0) step=0.1;
  // initialize MC state
  int size=SIZE*3;
  real* spins=(real*)malloc(sizeof(real)*size); assert(spins);
  real* spins2=(real*)malloc(sizeof(real)*size); assert(spins2);
  if(initial_state) {
    copy_vector(size, initial_state, spins);
    free(initial_state);
  } else random_vector(size, spins); 
  normalize(spins);
  real energy[6]; skyrmion_energy(spins, energy);
  long int mcs=1; long int iter=0;
  // start timer
  int per_iter=(update_by_one)?SIZE:1; // updates per MCS
  max_iter*=per_iter;
  const int initial_iter_rate=per_iter;
  long int iter_rate=initial_iter_rate;
  long int iter_last=0;
  while(mcs<max_iter) {
    if(stop_signal) { 
      fprintf(stderr, COLOR_RED"Aborting\n"COLOR_RESET); 
      break; 
    };
    iter++;
    // show progress
    if(iter>iter_rate+iter_last) {
      time_t now; time(&now); 
      time_t ltime=difftime(now,start_time);
      if(ltime>0) {
        iter_rate=mcs/ltime; if(iter_rate<=0) iter_rate=initial_iter_rate;
        iter_last=iter;
      };
      if(difftime(now,last_time_reported)>=debug_every_sec) {
        long int acceptence=100*mcs/iter;
        long int iter_per_sec=iter/ltime;
        last_time_reported=now;
        struct tm* timeinfo=localtime(&now);
        char buf[100]; strftime(buf, sizeof(buf), "%H:%M:%S", timeinfo);
        fprintf(stderr,COLOR_BLUE"%s"COLOR_RESET" ",buf);
        fprintf(stderr,COLOR_YELLOW"#"COLOR_RESET"%ld ", mcs/per_iter);
        fprintf(stderr,COLOR_YELLOW"Energy"COLOR_RESET" %.1"RF"g(%.1"RF"g) ",energy[5], energy[3]);
        fprintf(stderr,COLOR_YELLOW"Acc."COLOR_RESET" %ld%% ",acceptence);
        fprintf(stderr,COLOR_YELLOW"Iter./sec"COLOR_RESET" %ld ",iter_per_sec/per_iter);
        fprintf(stderr,COLOR_YELLOW"ETA"COLOR_RESET" ");
        fprint_timediff(stderr, (max_iter-mcs)*ltime/mcs);
        fprintf(stderr,"\n");

        skyrmion_energy(spins, energy);
        //fprintf(stderr,COLOR_YELLOW"Energy"COLOR_RESET" %"RF"g(%"RF"g)\n",energy[5], energy[3]);
      };    
    };
    // compute new state
    if(update_by_one) {
      // choose node
      int x=rand()%sizex, y=rand()%sizey, z=rand()%sizez, u=rand()%sizeu;
      int idx=INDEX(u,x,y,z)*3;
      real old[3]={spins[idx], spins[idx+1], spins[idx+2]}; 
      real oldenergy[6]; node_energy(u,x,y,z,spins,oldenergy);
      switch (random_number_generator) {
        case 0: add_random_vector3(2*step, spins+idx, spins+idx); break;
        case 1: add_random_cone3(step, spins+idx, spins+idx); break;
        default: fprintf(stderr, "Unknown noise type: %d\n", random_number_generator); exit(1);
      };
      normalize3(spins+idx);
      real newenergy[6]; node_energy(u,x,y,z,spins,newenergy);
      for(int k=0; k<6; k++) newenergy[k]-=oldenergy[k];

      real denergy=-newenergy[5];
      if(is_to_reject(denergy)) {
        spins[idx]=old[0]; spins[idx+1]=old[1]; spins[idx+2]=old[2]; 
        continue;
      };
      for(int k=0; k<6; k++) energy[k]+=newenergy[k];
 
      //real e0[6]; skyrmion_energy(spins, e0);
      //fprintf(stderr, "dT %d %d %d %d %"RF"g %"RF"g %"RF"g %"RF"g %"RF"g %"RF"g\n", u,x,y,z,e0[0]-energy[0], e0[1]-energy[1], e0[2]-energy[2], e0[3]-energy[3], e0[4]-energy[4], e0[5]-energy[5]);
      //fprintf(stderr, "   %"RF"g %"RF"g %"RF"g %"RF"g %"RF"g %"RF"g\n", old[0],old[1],old[2], spins[idx],spins[idx+1],spins[idx+2]);      
      //for(int k=0; k<6; k++) energy[k]=e0[k];      
      
    } else {
      switch (random_number_generator) {
        case 0: add_random_vector(2*step, size, spins, spins2); break;
        case 1: add_random_cone(step, SIZE, spins, spins2); break;
        default: fprintf(stderr, "Unknown noise type: %d\n", random_number_generator); exit(1);
      };
      normalize(spins2);
      // compute energy
      real last_energy=energy[5];
      skyrmion_energy(spins2, energy);

      real denergy=(last_energy-energy[5]);
      if(is_to_reject(denergy)) continue;
      real* tmp; tmp=spins; spins=spins2; spins2=tmp;
    };
    mcs++;
    if(mcs%per_iter==0) {
      // save result
      int ne=0; float ergy[6];
      if(what_to_write & WRITE_ANISOTROPY) ergy[ne++]=energy[0];
      if(what_to_write & WRITE_ZEEMAN) ergy[ne++]=energy[1];
      if(what_to_write & WRITE_EXCHANGE) ergy[ne++]=energy[2];
      if(what_to_write & WRITE_DMI) ergy[ne++]=energy[3];
      if(what_to_write & WRITE_DIPOLAR) ergy[ne++]=energy[4];
      if(what_to_write & WRITE_TOTAL) ergy[ne++]=energy[5];
      fwrite(ergy, sizeof(float), ne, file);
    };
  }
  fprintf(stderr, "Finished after %ld iterations %ld steps\n", iter/SIZE, mcs/SIZE);
  // deinitialization
  fclose(file);
  free(spins); free(spins2);
  return 0;
};