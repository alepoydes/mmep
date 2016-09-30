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
#include "debug.h"
#include "integra.h"
#include "cmd.h"

real step=0;
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

char options_desc[]="\
\n   -s REAL Monte-Carlo step size\
\n   -T REAL kT=Temperature by Bolzman factor\
\n   -R 0..1 Rejection rule\
\n   -N 0..1 Noise type\
\n   -u      Switch simulaneous/one by one updates\
\n   -F FILE Where to write result\
\n";

char handle_option(char opt, const char* arg) {
  switch(opt) {
    case 'u': update_by_one=!update_by_one; break;
    case 's': step=atof(arg); break;
    case 'T': temperature=atof(arg); break;            
    case 'N': random_number_generator=atoi(arg); break;
    case 'R': rejection_rule=atoi(arg); break;
    case 'F': strncpy(filename,arg,sizeof(filename)); break;
    default: return FALSE;
  };
  return TRUE;
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
  max_iter=5000000;

  int i=init_program(argc,argv,
    "Run Metropolis-Hastings algorithm.", options_desc,
    "us:T:R:N:F:", handle_option);
  if(i<argc) {
    fprintf(stderr, COLOR_RED "There are unused parameters:" COLOR_RESET "\n");
    while(i<argc) fprintf(stderr, "  %s\n", argv[i++]);
  };

  real debug_every_sec=debug_every;

  if(temperature==0) {
    fprintf(stderr, COLOR_RED "Temperature is not set or zero\n");
    exit(1);
  };

  print_settings();
  fprintf(stderr, COLOR_YELLOW "Temperature:" COLOR_RESET " %" RF "g\n", RT(temperature)); 
  fprintf(stderr, COLOR_YELLOW "Step size:" COLOR_RESET " %" RF "g\n", RT(step)); 
  fprintf(stderr, COLOR_YELLOW "Number of steps:" COLOR_RESET " %ld\n", max_iter); 
  fprintf(stderr, COLOR_YELLOW "Update type:" COLOR_RESET " %s\n", update_by_one?"one by one":"simultaneous"); 
  fprintf(stderr, COLOR_YELLOW "Rejection rule:" COLOR_RESET " %s\n", rejection_rule==0?"Metropolis-Hastings":"Nowak");   
  fprintf(stderr, COLOR_YELLOW "Noise type:" COLOR_RESET " %s\n", random_number_generator==0?"fast noise":"symmetric noise"); 

  // Initialize timer
  time_t last_time_reported;
  srand(time(&last_time_reported));
  time_t start_time=last_time_reported;   
  // open file for writing
  if(filename[0]==0) {
    struct tm* timeinfo=localtime(&start_time);
    char buf[100]; strftime(buf, sizeof(buf), "%H:%M:%S", timeinfo);
    pid_t pid=getpid();
    sprintf(filename, "mc.%s.%d.bin", buf, pid);
  };
  FILE* file=open_file(outdir, filename, TRUE);
  if(!file) exit(1);

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
  realp energy[6]; skyrmion_energy(spins, energy);
  long int mcs=0;
  // start timer
  int per_iter=(update_by_one)?SIZE:1; // updates per MCS
  max_iter*=per_iter;
  const int initial_iter_rate=per_iter;
  long int iter_rate=initial_iter_rate;
  long int iter_last=0;
  for(long int iter=0; iter<max_iter; iter++) {
    if(stop_signal) { 
      fprintf(stderr, COLOR_RED "Aborting\n" COLOR_RESET); 
      break; 
    };
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
        fprintf(stderr,COLOR_BLUE "%s" COLOR_RESET " ",buf);
        fprintf(stderr,COLOR_YELLOW "#" COLOR_RESET "%ld ", mcs/per_iter);
        fprintf(stderr,COLOR_YELLOW "Energy" COLOR_RESET " %.1" RPF "g(%.1" RPF "g) ",RT(energy[5]), RT(energy[3]));
        fprintf(stderr,COLOR_YELLOW "Acc." COLOR_RESET " %ld%% ",acceptence);
        fprintf(stderr,COLOR_YELLOW "Iter./sec" COLOR_RESET " %ld ",iter_per_sec/per_iter);
        fprintf(stderr,COLOR_YELLOW "ETA" COLOR_RESET " ");
        fprint_timediff(stderr, (max_iter-iter)*ltime/iter);
        fprintf(stderr,"\n");
        fflush(stderr);

        skyrmion_energy(spins, energy);
        //fprintf(stderr,COLOR_YELLOW "Energy" COLOR_RESET " %" RF "g(%" RF "g)\n",energy[5], energy[3]);
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
        //continue;
      } else {
        mcs++;
        for(int k=0; k<6; k++) energy[k]+=newenergy[k];
      };
 
      //real e0[6]; skyrmion_energy(spins, e0);
      //fprintf(stderr, "dT %d %d %d %d %" RF "g %" RF "g %" RF "g %" RF "g %" RF "g %" RF "g\n", u,x,y,z,e0[0]-energy[0], e0[1]-energy[1], e0[2]-energy[2], e0[3]-energy[3], e0[4]-energy[4], e0[5]-energy[5]);
      //fprintf(stderr, "   %" RF "g %" RF "g %" RF "g %" RF "g %" RF "g %" RF "g\n", old[0],old[1],old[2], spins[idx],spins[idx+1],spins[idx+2]);      
      //for(int k=0; k<6; k++) energy[k]=e0[k];      
      
    } else {
      switch (random_number_generator) {
        case 0: add_random_vector(2*step, size, spins, spins2); break;
        case 1: add_random_cone(step, SIZE, spins, spins2); break;
        default: fprintf(stderr, "Unknown noise type: %d\n", random_number_generator); exit(1);
      };
      normalize(spins2);
      // compute energy
      realp newenergy[6]; 
      skyrmion_energy(spins2, newenergy);
      realp denergy=(energy[5]-newenergy[5]);
      if(is_to_reject(denergy)) {
        //continue;
      } else {
        real* tmp; tmp=spins; spins=spins2; spins2=tmp;
        for(int k=0; k<6; k++) energy[k]=newenergy[k];
        mcs++;
      };
    };
    if(iter%per_iter==0) {
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
  };
  fprintf(stderr, COLOR_YELLOW "Done\n" COLOR_RESET);
  // deinitialization
  fclose(file);
  free(spins); free(spins2);
  return 0;
};