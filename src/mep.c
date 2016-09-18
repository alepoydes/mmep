#include "vector.h"
#include "skyrmion.h"
#include "optim.h"
#include "plot.h"
#include "parse.h"
#include "debug.h"
#include "octave.h"
#include "cmd.h"

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <getopt.h>

int size=0; // Dimenion of vector containing skyrmionic solutions
int sizep=0; // Number of nodes on path
int min_sizep=7; // Number of nodes on path
int max_sizep=25; // Number of nodes on path
int debug_plot_path=0;
int use_ftt=0;
int use_first_order_repar=0;
int remove_zero_modes=0;
real random_noise=0;
int skip_projection=0;
int do_not_relax_ends=0;

static int flat_distance=1;

realp *distance=NULL, *energy=NULL, *diff=NULL, *tdiff=NULL, *inflation=NULL;
int post_optimization=0;
int single_maximum=1;

void energy_evaluate_neb(real* path) {
  real* q=(real*)malloc(sizeof(real)*size); assert(q);
  real* u=(real*)malloc(sizeof(real)*size); assert(u);
  for(int p=0; p<sizep; p++) {
    if(flat_distance)
      distance[p]=(p<=0)?0:distance[p-1]+rpsqrt(distsq(size,path+size*(p-1),path+size*p)/SIZE);
    else distance[p]=(p<=0)?0:distance[p-1]+rpsqrt(dist_sphere_sq(SIZE,path+size*(p-1),path+size*p)/SIZE);
    hamiltonian_hessian(path+size*p, q);
    energy[p]=-dot(size,path+size*p, q)/2;
    subtract_field(q);
    energy[p]+=dot(size,path+size*p, q);
    // gradient of energy is in q
    project_to_tangent(path+size*p,q);
    // q is orthogonal to normales to all spheres
    diff[p]=rpsqrt(normsq(size,q)/size);
    tdiff[p]=diff[p];
    if(p>0 && p<sizep-1) {
      three_point_tangent_stable(energy[p-1],energy[p],energy[p+1],path+size*(p-1),path+size*p,path+size*(p+1),u);
      //three_point_tangent(path+size*(p-1),path+size*p,path+size*(p+1),u);
      real proj=dot(size,u,q)/dot(size,u,u);
      mult_sub(size, proj, u, q);
      // q is orthogonal to tangent to MEP
      tdiff[p]=rpsqrt(normsq(size,q)/size);
    };
    inflation[p]=NAN;
  };
  free(u); free(q);
};

#define circperm(a,b,c) { real t=a; a=b; b=c; c=t; }
void energy_evaluate_ftt(real* path) {
  real* q=(real*)malloc(sizeof(real)*size*sizep); assert(q);
  real* u=(real*)malloc(sizeof(real)*size*sizep); assert(u);
  real* v=(real*)malloc(sizeof(real)*size); assert(v);
  for(int p=0; p<sizep; p++) {
    if(flat_distance)
      distance[p]=(p<=0)?0:distance[p-1]+rpsqrt(distsq(size,path+size*(p-1),path+size*p)/SIZE);
    else distance[p]=(p<=0)?0:distance[p-1]+rpsqrt(dist_sphere_sq(SIZE,path+size*(p-1),path+size*p)/SIZE);
    hamiltonian_hessian(path+size*p, q+size*p);
    energy[p]=-dot(size,path+size*p, q+size*p)/2;
    subtract_field(q+size*p);
    energy[p]+=dot(size,path+size*p, q+size*p);
    // gradient of energy is in q
    project_to_tangent(path+size*p,q+size*p);
    // q is orthogonal to normales to all spheres
    tdiff[p]=diff[p]=rpsqrt(normsq(size,q+size*p)/size);
  };
  for(int p=0; p<sizep; p++) {
    if(p>0 && p<sizep-1) {
      three_point_tangent_stable(energy[p-1],energy[p],energy[p+1],path+size*(p-1),path+size*p,path+size*(p+1),u+size*p);
      //three_point_tangent(path+size*(p-1),path+size*p,path+size*(p+1),u+size*p);
    } else if (p==0) {
      two_point_tangent0(path+size*p,path+size*(p+1),u+size*p);
    } else {
      two_point_tangent1(path+size*(p-1),path+size*p,u+size*p);
    };
    real unorm=rpsqrt(dot(size,u+size*p,u+size*p));
    real proj=dot(size,u+size*p,q+size*p)/unorm;
    copy_vector(size,q+size*p,v);
    mult_sub(size, proj/unorm, u+size*p, v);
    // v is orthogonal to the tangent to MEP at p
    tdiff[p]=rpsqrt(normsq(size,v)/size);
    if(p>0 && p<sizep-1) {
      sub(size,q+(p+1)*size,q+(p-1)*size,v);
      inflation[p]=-dot(size,v,u+p*size)/2/unorm;
    } else if(p==0) {
      sub(size,q+(p+1)*size,q+(p)*size,v);
      inflation[p]=-dot(size,v,u+p*size)/unorm;
    } else {
      sub(size,q+(p)*size,q+(p-1)*size,v);
      inflation[p]=-dot(size,v,u+p*size)/unorm;
    };
  };

  //for(int p=1; p<sizep; p++)
  //  fprintf(stderr,"%d: %" RF "g %" RF "g %" RF "g\n",p,distance[p]-distance[p-1], inflation[p],energy[p]);
  free(v); free(u); free(q); 
};

void energy_evaluate(real* path) {
  if(use_ftt) energy_evaluate_ftt(path);
  else energy_evaluate_neb(path);
}

void skyrmion_display(int iter, real* __restrict__ a, real* __restrict__ grad_f, realp f, 
real res, real constres, real alpha, realp last_f, real last_res) {
  static realp prev_f=NAN; static realp prev_res=NAN; 
  if(iter==1) { prev_f=f; prev_res=res; };
  if(iter%debug_every==0 || iter<0) {
    if(iter!=0) fprintf(stderr, "%6d", abs(iter));
    else fprintf(stderr, "%6s", "");
    fprintf(stderr, " " COLOR_YELLOW "E" COLOR_RESET);
    watch_number(f,debug_every==1?last_f:prev_f,16);
    fprintf(stderr, " " COLOR_YELLOW "R" COLOR_RESET);
    watch_number(res,debug_every==1?last_res:prev_res,16);
    fprintf(stderr, "%+" RF "g", constres);
    fprintf(stderr, " " COLOR_YELLOW "A" COLOR_RESET "%" RF "g", alpha);
    //fprintf(stderr, " " COLOR_YELLOW "d" COLOR_RESET "%.2g", pow(last_res/res,1./alpha));    
    fprintf(stderr, "\n");
  	if(debug_plot && iter!=0 && size<=1024) plot_field3(stdout,a);
    prev_f=f; prev_res=res;
  };
};

real quasynorm(real* __restrict__ x) {
  normalize(x); return 0;
};

int skyrmion_steepest_descent(real* __restrict__ x, int mode, real mode_param, 
	real epsilon, int max_iter) 
{
	return steepest_descend(
		size, (real*)x, 
    projected_gradient,
		mode, mode_param, epsilon, max_iter,
		skyrmion_display, 
		quasynorm
	);
};

void path_gradient(const real* __restrict__ arg, real* __restrict__ out, realp* __restrict__ E) {
  if(E) *E=0;
  for(int p=0; p<sizep; p++) 
    if(E) {
      realp locE; skyrmion_gradient(arg+size*p, out+size*p, &locE); (*E)+=locE;
    } else skyrmion_gradient(arg+size*p, out+size*p, NULL);
};

void energy_display(FILE* file) {
  fprintf(file,"set ytics nomirror\nset y2tics nomirror\nset log y2\n");
  fprintf(file,"set format y2 '10^{%%T}'\nset autoscale xy\n");
  fprintf(file,"plot '-' using 1:2 with linespoints axes x1y1 title 'energy', '' using 1:3 with lines axes x1y2 title 'grad.', '' using 1:4 with lines axes x1y2 title 'orth. grad.'\n");
  for(int k=0; k<3; k++) {
    for(int p=0; p<sizep; p++)
      fprintf(file,"%.8" RPF "e %.8" RPF "e %.8" RPF "e %.8" RPF "e\n",distance[p],energy[p],diff[p],tdiff[p]);
    fprintf(file,"EOF\n\n");
  };
  fflush(file);
};

void path_display(int iter, real* __restrict__ mep, real* __restrict__ grad_f
, realp f, real res, real constres, real alpha, realp last_f, real last_res) {
  if(!post_optimization && res<10*epsilon) {
      //fprintf(stderr,COLOR_YELLOW "Climbing is on." COLOR_RESET " Residue %" RF "g\n",res);
      post_optimization=1;
  }; 
  static realp prev_f=NAN; static realp prev_res=NAN;
  if(iter==1) { prev_f=f; prev_res=res; };
  if(iter%debug_every==0 || iter<0) {
    if(iter!=0) fprintf(stderr, "%6d", abs(iter));
    else fprintf(stderr, "%6s", "");
    fprintf(stderr, " " COLOR_YELLOW "E" COLOR_RESET);
    watch_number(f/(sizep-1),debug_every==1?last_f/(sizep-1):prev_f/(sizep-1),16);
    fprintf(stderr, " " COLOR_YELLOW "R" COLOR_RESET);
    watch_number(res,debug_every==1?last_res:prev_res,16);
    fprintf(stderr, "%+.2" RF "g", constres);
    fprintf(stderr, " " COLOR_YELLOW "A" COLOR_RESET "%" RF "g", alpha);
    fprintf(stderr, COLOR_FAINT " %s" COLOR_RESET, post_optimization?"Climbing":"");        
    fprintf(stderr, "\n");
    if(debug_plot && iter!=0) {
      if(debug_plot_path) plot_path(stdout, sizep, mep);
      else { /*energy_evaluate(mep);*/ energy_display(stdout); };
    };
    prev_f=f; prev_res=res;
  };
};


// Move nodes from+1..to-1 along path so that 
// all intervals between nodes in interval from..to are 
// of the same length
void path_equilize_rec(real* mep, int from, int to) {  
  //fprintf(stderr, "path_equilize_rec: %d %d\n",from, to);
  assert(from>=0); assert(to<sizep);
  if(to-from<2) return;
  int* idx=(int*)alloca(sizeof(int)*sizep);
  realp loc[sizep];
  realp dist[sizep];
  for(int p=0; p<sizep; p++) dist[p]=distance[p];  
  for(int p=to-1; p>from; p--) {
    realp d=distance[from]+(p-from)*(distance[to]-distance[from])/(to-from); // desired position
    int f=to; while(f>from) if(distance[--f]<=d) break; 
    assert(f>=from); assert(f<=to); 
    // move the image to interval [dist(f),dist(f+1)]
    idx[p]=f;
    loc[p]=(distance[f+1]>distance[f])?(d-distance[f])/(distance[f+1]-distance[f]):0; // local coordinate
    //fprintf(stderr,"[%d] %" RF "g [%d] %" RF "g + %" RF "g * %" RF "g\n",p,d,f,distance[f],loc,distance[f+1]-distance[f]);
    assert(loc[p]>=0); assert(loc[p]<=1);
    dist[p]=d;
  };
  idx[from]=from; loc[from]=0;
  idx[to]=to; loc[to]=0;

  //for(int p=from; p<=to; p++) fprintf(stderr,COLOR_GREEN "%d" COLOR_RESET ":%.3" RF "g(%.3" RF "g %d %.3" RF "g) ",p,distance[p],dist[p],idx[p], loc[p]); fprintf(stderr,"\n");      
  /*for(int p=from; p<=to; p++) {
     if(abs(idx[p])>p) fprintf(stderr, COLOR_RED); else if(abs(idx[p])<p) fprintf(stderr, COLOR_BLUE);
     if(idx[p]<0) fprintf(stderr, COLOR_BOLD);
     fprintf(stderr, "%d " COLOR_RESET, abs(idx[p])-from);
  }; fprintf(stderr, "\n");  
  */

  // check if the images can be interpolated from right to left
  int f; for(f=from+1; f<to; f++) if(idx[f]>=f) break;
  if(f>=to) { // all idx[f]<f, hence we can proceed interpolate inplace
    for(int j=to-1; j>from; j--) {
      int p=j;
      assert(p>from && p<to);
      assert(idx[p]>=from && idx[p+1]<=to);
      assert(idx[idx[p]]>=0 && idx[idx[p]+1]>=0); 
      linear_comb(size,1-loc[p],mep+idx[p]*size,loc[p],mep+(idx[p]+1)*size,mep+p*size);
      idx[p]=-idx[p];
    };
  } else {
    //fprintf(stderr, COLOR_FAINT "Failed to interpolate inplace\n" COLOR_RESET);
    // copying old images to a buffer
    real* buf=(real*)malloc(sizeof(real)*size*(to-from+1)); assert(buf);
    copy_vector(size*(to-from+1), mep+size*from, buf);
    for(int p=to-1; p>from; p--) {
      assert(p>from && p<to);
      assert(idx[p]>=from && idx[p+1]<=to);
      linear_comb(size,1-loc[p],buf+(idx[p]-from)*size,loc[p],buf+(idx[p]-from+1)*size,mep+p*size);
    };
    free(buf);
  };

  for(int p=0; p<sizep; p++) distance[p]=dist[p];
};

void path_equilize(real* mep) {
  if(sizep<2) return;
  if(post_optimization) {
    if(single_maximum) {
      real max=-INFINITY; int f=-1;
      for(int p=sizep-1; p>=0; p--) 
        if(energy[p]>max) { max=energy[p]; f=p; };
      if(f>0 && f<sizep-1) {
        path_equilize_rec(mep, 0, f);
        path_equilize_rec(mep, f, sizep-1);
      } else path_equilize_rec(mep, 0, sizep-1);
    } else {
      int end=sizep-1;
      for(int p=end-1; p>0; p--) 
        if((energy[p]<energy[p-1] && energy[p]<energy[p+1]) 
          || (energy[p]>energy[p-1] && energy[p]>energy[p+1]))
        { // Found a stationary point. Let it move to the extremum 
          path_equilize_rec(mep, p, end);
          end=p;
        };
      path_equilize_rec(mep, 0, end);
    };
  } else path_equilize_rec(mep, 0, sizep-1);
};

real path_normalize(real* mep) {
  if(!use_ftt || !use_first_order_repar) {
    energy_evaluate(mep);
    path_equilize(mep);
  };
  real sum=0;
  for(int p=0; p<sizep; p++) sum+=normalize(mep+size*p);  
  //energy_evaluate(mep);
  return sum;
};

void path_tangent_rec(const real* __restrict__ mep, real* __restrict__ grad, int from, int to) {
  real* u=ralloc(size);
  int C=to-from+1;  
  //fprintf(stderr, "[%d-%d] ", from,to);
  real* q=ralloc(C); 
  real* l=ralloc(C); 
  l[0]=0; q[0]=0;
  for(int p=from+1; p<=to; p++) {
    l[p-from]=distance[p]-distance[from];
    q[p-from]=q[p-from-1]+(inflation[p]+inflation[p-1])/2;
  };
  // compute tangent force for every non-extremal node
  for(int p=from; p<=to; p++) {
    if(p>from && p<to) {
      three_point_tangent_stable(energy[p-1],energy[p],energy[p+1],mep+size*(p-1), mep+p*size, mep+size*(p+1), u);
      //three_point_tangent(mep+size*(p-1), mep+p*size, mep+size*(p+1), u);
      real normu=rsqrt(normsq(size,u));
      real proj=(l[p-from]/l[C-1]*q[C-1]-q[p-from]);// inflation compensation
      //proj=dot(size,u,grad+p*size)/normu; // gradient tangential projection
      if(use_first_order_repar)
        proj+=(l[C-1]*(p-from)/(C-1)-l[p-from])*rsqrt(size);
      //fprintf(stderr,"%d$ %" RF "g %" RF "g %" RF "g %" RF "g %" RF "g\n",p,proj,l[p],q[p],distance[p],inflation[p]);
      mult_add(size, -proj/normu, u, grad+p*size);
    } else if(
        post_optimization 
        && ( (p==from && from>0 && energy[from]>energy[to]) 
           //||(p==to && to<sizep-1 && energy[from]<energy[to])
           )
      ) 
    {
      three_point_tangent_stable(energy[p-1],energy[p],energy[p+1],mep+size*(p-1), mep+p*size, mep+size*(p+1), u);
      //three_point_tangent(mep+size*(p-1), mep+p*size, mep+size*(p+1), u);
      real proj=dot(size,u,grad+p*size)/dot(size,u,u);
      mult_add(size, -2*proj, u, grad+p*size);
    };
  };
  free(u); free(q); free(l);
}

void path_tangent_ftt(const real* __restrict__ mep, real* __restrict__ grad) {
  // Project to the tangent subspaces
  for(int p=0; p<sizep; p++) 
    project_to_tangent(mep+size*p,grad+size*p);
  // Apply spring forces on every monotone segment of MEP
  if(sizep<2) return;
  if(post_optimization) {
    if(single_maximum) {
      realp max=-INFINITY; int f=-1;
      for(int p=sizep-1; p>=0; p--) 
        if(energy[p]>max) { max=energy[p]; f=p; };
      if(f>0 && f<sizep-1) {
        path_tangent_rec(mep, grad, 0, f);
        path_tangent_rec(mep, grad, f, sizep-1);
      } else path_tangent_rec(mep, grad, 0, sizep-1);
    } else {
      int end=sizep-1;
      for(int p=end-1; p>0; p--) 
        if((energy[p]<energy[p-1] && energy[p]<energy[p+1]) 
          || (energy[p]>energy[p-1] && energy[p]>energy[p+1]))
        { // Found a stationary point. Let it move to the extremum 
          path_tangent_rec(mep, grad, p, end);
          end=p;
        };
      path_tangent_rec(mep, grad, 0, end);
    };
  } else path_tangent_rec(mep, grad, 0, sizep-1);  
};

// INVALID
void path_tangent_rohart(const real* __restrict__ mep, real* __restrict__ grad) {
  real* u=ralloc(size);
  for(int p=0; p<sizep; p++) {
    //real l1=normsq(size, grad+size*p);
    project_to_tangent(mep+size*p,grad+size*p);
    //real l2=normsq(size, grad+size*p);
    //real l3=0, l4=0;
    if(p>0 && p<sizep-1) {
      three_point_tangent_stable(energy[p-1],energy[p],energy[p+1],mep+size*(p-1), mep+p*size, mep+size*(p+1), u);
      normalize(u);
      project_to_tangent(u,grad+size*p);
      //l3=normsq(size, grad+size*p);
      //project_to_tangent(mep+size*p,grad+size*p);
      //l4=normsq(size, grad+size*p);
    };
    //fprintf(stderr, "  %" RF "g  %" RF "g  %" RF "g  %" RF "g\n", l1, l2, l3, l4);
  };
  free(u); 
}

void path_tangent_neb(const real* __restrict__ mep, real* __restrict__ grad) {
  real* u=ralloc(size);
  real *g1=NULL, *g2=NULL, *g3=NULL; 
  if(remove_zero_modes) {
    g1=ralloc(size);
    g2=ralloc(size);
    g3=ralloc(size);
  };
  // Find index f of node with
  realp max=-INFINITY; int f=-1;
  for(int p=sizep-1; p>=0; p--) if(energy[p]>max) { max=energy[p]; f=p; };
  if(f==0 || f==sizep-1) f=-1;
  for(int p=0; p<sizep; p++) {
    project_to_tangent(mep+size*p,grad+size*p);
    if(p>0 && p<sizep-1) {
      //if(post_optimization && energy[p]<energy[p-1] && energy[p]<energy[p+1]) {
      //} else 
      if(post_optimization && (p==f || (energy[p]>=energy[p-1] && energy[p]>=energy[p+1] && !single_maximum))) {
        three_point_tangent_stable(energy[p-1],energy[p],energy[p+1],mep+size*(p-1), mep+p*size, mep+size*(p+1), u);
        //three_point_tangent(mep+size*(p-1), mep+p*size, mep+size*(p+1), u);
        real proj=dot(size,u,grad+p*size)/dot(size,u,u);
        mult_add(size, -2*proj, u, grad+p*size);
      } else {
        //three_point_tangent(mep+size*(p-1), mep+p*size, mep+size*(p+1), u);
        if(remove_zero_modes) {
          three_point_tangent_stable(energy[p-1],energy[p],energy[p+1],mep+size*(p-1), mep+p*size, mep+size*(p+1), u);
          group_generator(mep+size*p, 0, g1);
          group_generator(mep+size*p, 1, g2);
          group_generator(mep+size*p, 2, g3);
          real* basis[5]={g1,g2,g3,u,grad+p*size};
          gram_schmidt(size, 5, basis);
        } else if(!skip_projection) {
          three_point_tangent_stable(energy[p-1],energy[p],energy[p+1],mep+size*(p-1), mep+p*size, mep+size*(p+1), u);
          real* basis[2]={u,grad+p*size};
          gram_schmidt(size, 2, basis);
          //real normu=rsqrt(normsq(size,u));
          //real proj=dot(size,u,grad+p*size)/normu;
          //mult_add(size, -proj/normu, u, grad+p*size);
        };
      };
    };
  };
  free(u); if(g1)free(g1); if(g2)free(g2); if(g3)free(g3);
}

void path_tangent(const real* __restrict__ mep, real* __restrict__ grad) {
  // Rohart way
  //path_tangent_rohart(mep, grad); return;
  if(use_ftt) path_tangent_ftt(mep, grad);
  else path_tangent_neb(mep, grad);
}  

void projected_path_gradient(const real* __restrict__ arg, real* __restrict__ out, realp* __restrict__ E) {
  path_gradient(arg, out, E);
  path_tangent(arg, out);
};

int path_steepest_descent(real* __restrict__ path, int mode, 
  real mode_param, real epsilon, int max_iter) 
{
  real updated_param=mode_param;
  if(mode==2 || mode==0) updated_param=mode_param/sizep;
  return steepest_descend(
    size*sizep, (real*)path, 
    projected_path_gradient,
    mode, updated_param, epsilon, max_iter,
    path_display, 
    path_normalize
  );
};

const char options_desc[]="\
\n   -P      Enable GNUPlot output of MEP\
\n   -n INT  Minimum number of images on MEP\
\n   -N INT  Maximum number of images on MEP\
\n   -0      Use FTT instead of NEB method\
\n   -z      Disable translations preserving energy\
\n   -R REAL Noise amplitude for initial path\
\n   -q      Skip ends relaxation\
\n";

char handle_option(char opt, const char* arg) {
  switch(opt) {
    case 'P': debug_plot=1; debug_plot_path=1; break;
    case 'n': min_sizep=atoi(optarg); break;
    case 'N': max_sizep=atoi(optarg); break;
    case 'z': remove_zero_modes=1; break;
    case 'R': random_noise=atof(optarg); break;
    case 'q': do_not_relax_ends=1; break;
    case '0': use_ftt=1; break;
    default: return FALSE;
  };
  return TRUE;
};

int main(int argc, char** argv) {
  int i=init_program(argc,argv,
    "Calculate MEP for magnetic systems.", options_desc,
    "zPN:n:R:q0", handle_option);
  if(i<argc) {
    fprintf(stderr, COLOR_RED "There are unused parameters:" COLOR_RESET "\n");
    while(i<argc) fprintf(stderr, "  %s\n", argv[i++]);
  };

  // Initializaton
  print_settings();
  fprintf(stderr, "%s method in use\n", use_ftt?"FTT":"NEB");
  if(remove_zero_modes)
    fprintf(stderr, "Zero modes (translations) are removed\n");
  if(random_noise>0) 
    fprintf(stderr, "Initial path noise amplitude: %" RF "g\n", random_noise);

  size=SIZE*3;
  real* path=ralloc(size*max_sizep); 
  distance=(realp*)malloc(sizeof(realp)*max_sizep); assert(distance);
  energy=(realp*)calloc(sizeof(realp),max_sizep); assert(energy);
  diff=(realp*)malloc(sizeof(realp)*max_sizep); assert(diff);
  tdiff=(realp*)malloc(sizeof(realp)*max_sizep); assert(tdiff);
  inflation=(realp*)malloc(sizeof(realp)*max_sizep); assert(inflation);
  // Set initial path size
  sizep=min_sizep; 
  if(sizep<initial_states_count) sizep=initial_states_count;
  if(max_sizep<sizep) sizep=max_sizep;

  fprintf(stderr, "Images on path: %d in [%d, %d]\n", sizep, min_sizep, max_sizep);  
  // find two minima
  fprintf(stderr, COLOR_YELLOW COLOR_BOLD "Initializing path\n" COLOR_RESET);
  if(initial_states_count<2) {
    fprintf(stderr, COLOR_RED COLOR_BOLD "There should be at least two images in the path\n" COLOR_RESET);
    exit(1);
  };
  copy_vector(size, initial_state, path);
  int last_image=0;
  for (int n=1; n<initial_states_count; n++) {
    int next_image=n*sizep/(initial_states_count-1)-1;
    while(last_image>=next_image) next_image=last_image+1;
    assert(next_image<sizep); 
    copy_vector(size, initial_state+size*n, path+size*next_image);
    // Set initial path as geodesic approximation between given states
    skyrmion_geodesic(random_noise/(next_image-last_image), next_image-last_image+1, path+size*last_image);
    last_image=next_image;
  };
  assert(last_image==sizep-1);

  // minimize initial and final states
  
  if(!do_not_relax_ends) {
    fprintf(stderr, COLOR_YELLOW COLOR_BOLD "Relaxing initial state\n" COLOR_RESET);
    skyrmion_steepest_descent(path, mode, mode_param, epsilon, max_iter*max_sizep);
    fprintf(stderr, COLOR_YELLOW COLOR_BOLD "Relaxing final state\n" COLOR_RESET);
    skyrmion_steepest_descent(path+size*(sizep-1), mode, mode_param, epsilon, max_iter*max_sizep);
  };

  fprintf(stderr, COLOR_YELLOW COLOR_BOLD "Calculating MEP\n" COLOR_RESET);
  // MEP calculation
  post_optimization=0;
  path_steepest_descent(path, mode, mode_param, epsilon, max_iter);
  while(2*(sizep-1)+1<=max_sizep) {
    // interpolating path
    for(int p=sizep-1; p>0; p--) copy_vector(size, path+p*size, path+2*p*size);
    for(int p=1; p<sizep; p++) {
      if(p==1) skyrmion_middle_third_order(path+2*size*(p-1), path+2*size*p, path+2*size*(p+1), path+size*(2*p-1));
      else if(p==sizep-1) skyrmion_middle_third_order(path+2*size*p, path+2*size*(p-1), path+2*size*(p-2), path+size*(2*p-1));
      else skyrmion_middle_fourth_order(path+2*size*(p-2), path+2*size*(p-1), path+2*size*p, path+2*size*(p+1), path+size*(2*p-1));
      //skyrmion_middle(path+2*size*(p-1), path+2*size*p, path+size*(2*p-1));
    };
    sizep=2*(sizep-1)+1;
    fprintf(stderr, COLOR_YELLOW COLOR_BOLD "Increasing number of images to %d\n" COLOR_RESET, sizep);
    post_optimization=0;    
    path_steepest_descent(path, mode, mode_param, epsilon, max_iter);
  };
  // Ouput result
  energy_evaluate(path);
  realp max_energy=energy[0];
  realp min_energy=energy[0];  
  for(int p=1; p<sizep; p++) {
    if(max_energy<energy[p]) max_energy=energy[p];
    if(min_energy>energy[p]) min_energy=energy[p];    
  };
  fprintf(stderr, COLOR_BLUE "Energy:" COLOR_RESET " initial %.8" RPF "f maximum %.8" RPF "f minimum %.8" RPF "f final %.8" RPF "f\n", energy[0],max_energy,min_energy,energy[sizep-1]);
  fprintf(stderr, COLOR_BLUE "Barriers:" COLOR_RESET " forward %" RPF "g backward %" RPF "g\n", max_energy-energy[0], max_energy-energy[sizep-1]);  
  if(!debug_plot) {
    for(int p=0;p<sizep;p++) printf("%.8" RPF "g ", energy[p]);
    printf("\n");
    for(int p=0;p<sizep;p++) printf("%.8" RPF "g ", distance[p]);
    printf("\n");
  };
  // Compare distances between images
  real mean=distance[sizep-1]/(sizep-1); 
  real var=0; for(int p=1; p<sizep; p++) {
    real d=distance[p]-distance[p-1]-mean;
    var+=d*d;
  };  
  fprintf(stderr, COLOR_BLUE COLOR_BOLD "Std of distances between images:" COLOR_RESET " %" RF "g\n", rsqrt(var));

  // save energy
  fprintf(stderr, COLOR_YELLOW COLOR_BOLD "Saving result\n" COLOR_RESET);
  FILE* file=open_file(outdir,"/energy.gnuplot",TRUE);
  if(file) {
    fprintf(file,"set terminal png\nset output '%s/energy.png'\n", outdir);
    energy_display(file);
    fclose(file);
  }; 
  // save field
  file=open_file(outdir, "/mep.gnuplot", TRUE);
  if(file) {
    fprintf(file,"set terminal gif animate delay 10\n");
    fprintf(file,"set output '%s/mep.gif'\n",outdir);
    animate_path(file, sizep, path);
    fclose(file);
  };
  // Octave output
  if(save_octave>0) {
    fprintf(stderr, COLOR_YELLOW "Computing energy contributions\n" COLOR_RESET);
    realp* contr=(realp*)malloc(sizeof(realp)*sizep*6);
    for(int p=0; p<sizep; p++) skyrmion_energy(path+p*SIZE*3, contr+6*p);
    file=open_file(outdir, "/mep.oct", TRUE);    
    if(file) {
      oct_save_init(file);
      oct_save_lattice(file);
      oct_save_path(file,"PATH",path,sizep);
      oct_save_vectorp(file,"ENERGY",energy,sizep);
      oct_save_vectorp(file,"DISTANCE",distance,sizep);
      oct_save_vectorp(file,"DIFF",diff,sizep);
      oct_save_vectorp(file,"TDIFF",tdiff,sizep);
      oct_save_matrixp(file,"CONTRIBUTIONS",contr,sizep,6);
      if(save_octave>1) oct_save_hessian(file);
      oct_save_finish(file);
      fclose(file);
    };
  };
  // Deinitialization
  free(path); free(tdiff); free(diff); free(energy); 
  free(distance); free(inflation);
  return 0;
};