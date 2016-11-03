#include "magopt.h"
#include "cmd.h"
#include "debug.h"
#include "plot.h"
#include "skyrmion.h"
#include "integra.h"
#include "optim.h"

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

//////////////////////////////////////////////////////
// Single image optimization

void skyrmion_display(int iter, real* __restrict__ a, real* __restrict__ grad_f, realp f, real res, 
real constres, real alpha, realp last_f, real last_res) {
  static realp prev_f=NAN; static real prev_res=NAN; 
  if(iter==1) { prev_f=f; prev_res=res; };
  if(iter%debug_every==0 || iter<0) {
  	if(iter!=0) fprintf(stderr, "%6d", abs(iter));
    else fprintf(stderr, "%6s", "");
    fprintf(stderr, " " COLOR_YELLOW "E" COLOR_RESET);
    watch_number(f,debug_every==1?last_f:prev_f,16);
    fprintf(stderr, " " COLOR_YELLOW "R" COLOR_RESET);
    watch_number(res,debug_every==1?last_res:prev_res,16);
    fprintf(stderr, "%+.2" RF "g", RT(constres));
    fprintf(stderr, " " COLOR_YELLOW "A" COLOR_RESET "%" RF "g", RT(alpha));
    //fprintf(stderr, " " COLOR_YELLOW "d" COLOR_RESET "%.2g", pow(last_res/res,1./alpha));
    fprintf(stderr, "\n");
  	if(debug_plot && iter!=0) 
      plot_field3(stdout,a);
      //plot_field3(stdout,grad_f);
    prev_f=f; prev_res=res;
  };
};

real quasynorm(real* __restrict__ x) {
  return rsqrt(normalize(x));
};

void integrate(int N, void (*F)(const real* x, real* g, realp* E), real T, real* X, realp* E, int* iter) {
  real err=0;
  switch(integrator) {
    case 0: euler(N, F, T, X, E, iter); break;
    case 1: runge_kutta(N, F, T, X, E, iter); break;
    case 2: err=gauss_integrator(N, F, T, X, 1e-7, 10, E, iter); break;
    case 3: err=radau_integrator(N, F, T, X, 1e-7, 10, E, iter); break;
    default: fprintf(stderr, COLOR_RED "Unknown integrator:" COLOR_RESET " %d\n", integrator); exit(1);
  };
  if(err>1e-5) 
    fprintf(stderr, COLOR_RED "Warning:" COLOR_RESET " Runge-Kutta convergence error: %" RF "g\n", RT(err));
};

int skyrmion_steepest_descent(real* __restrict__ x, 
int mode, real mode_param, real epsilon, int max_iter) 
{
/*
  return flow_descend(
    SIZE*3, (real*)x, 
    projected_gradient,
    mode, mode_param, epsilon, max_iter,
    skyrmion_display, 
    quasynorm,
    integrate
  );
*/
	return steepest_descend(
		SIZE*3, (real*)x, 
    projected_gradient,
		mode, mode_param, epsilon, max_iter,
		skyrmion_display, 
		quasynorm
	);
};


//////////////////////////////////////////////////////
// EMP optimization

int sizep=0; // Number of nodes on path
int debug_plot_path=0;
int remove_zero_modes=0;
int flat_distance=1;
int post_optimization=0;
int single_maximum=1;

realp *distance=NULL, *energy=NULL, *diff=NULL, *tdiff=NULL, *inflation=NULL;

void path_steepest_descent_init(int max_sizep) {
  distance=(realp*)malloc(sizeof(realp)*max_sizep); assert(distance);
  energy=(realp*)calloc(sizeof(realp),max_sizep); assert(energy);
  diff=(realp*)malloc(sizeof(realp)*max_sizep); assert(diff);
  tdiff=(realp*)malloc(sizeof(realp)*max_sizep); assert(tdiff);
  inflation=(realp*)malloc(sizeof(realp)*max_sizep); assert(inflation);
};

void path_steepest_descent_deinit() {
  free(distance);
  free(energy);
  free(diff);
  free(tdiff);
  free(inflation);
};

void energy_evaluate_neb(real* path) {
  real* q=ralloc(3*SIZE);
  real* u=ralloc(3*SIZE); 
  for(int p=0; p<sizep; p++) {
    if(flat_distance)
      distance[p]=(p<=0)?0:distance[p-1]+rpsqrt(distsq(3*SIZE,path+3*SIZE*(p-1),path+3*SIZE*p));
    else distance[p]=(p<=0)?0:distance[p-1]+rpsqrt(dist_sphere_sq(SIZE,path+3*SIZE*(p-1),path+3*SIZE*p));
    hamiltonian_hessian(path+3*SIZE*p, q);
    energy[p]=-dot(3*SIZE,path+3*SIZE*p, q)/2;
    subtract_field(q);
    energy[p]+=dot(3*SIZE,path+3*SIZE*p, q);
    // gradient of energy is in q
    project_to_tangent(path+3*SIZE*p,q);
    // q is orthogonal to normales to all spheres
    diff[p]=rpsqrt(normsq(3*SIZE,q)/(SIZE));
    tdiff[p]=diff[p];
    if(p>0 && p<sizep-1) {
      three_point_tangent_stable(energy[p-1],energy[p],energy[p+1],path+3*SIZE*(p-1),path+3*SIZE*p,path+3*SIZE*(p+1),u);
      //three_point_tangent(path+3*SIZE*(p-1),path+3*SIZE*p,path+3*SIZE*(p+1),u);
      real proj=dot(3*SIZE,u,q)/dot(3*SIZE,u,u);
      mult_sub(3*SIZE, proj, u, q);
      // q is orthogonal to tangent to MEP
      tdiff[p]=rpsqrt(normsq(3*SIZE,q)/(SIZE));
    };
    inflation[p]=NAN;
  };
  free(u); free(q);
};

void energy_evaluate(real* path) {
  energy_evaluate_neb(path);
}

void path_gradient(const real* __restrict__ arg, real* __restrict__ out, realp* __restrict__ E) {
  if(E) *E=0;
  for(int p=0; p<sizep; p++) 
    if(E) {
      realp locE; skyrmion_gradient(arg+3*SIZE*p, out+3*SIZE*p, &locE); (*E)+=locE;
    } else skyrmion_gradient(arg+3*SIZE*p, out+3*SIZE*p, NULL);
};

void energy_display(FILE* file) {
  fprintf(file,"set ytics nomirror\nset y2tics nomirror\nset log y2\n");
  fprintf(file,"set format y2 '10^{%%T}'\nset autoscale xy\n");
  fprintf(file,"plot '-' using 1:2 with linespoints axes x1y1 title 'energy', '' using 1:3 with lines axes x1y2 title 'grad.', '' using 1:4 with lines axes x1y2 title 'orth. grad.'\n");
  for(int k=0; k<3; k++) {
    for(int p=0; p<sizep; p++)
      fprintf(file,"%.8" RPF "e %.8" RPF "e %.8" RPF "e %.8" RPF "e\n",RT(distance[p]),RT(energy[p]),RT(diff[p]),RT(tdiff[p]));
    fprintf(file,"EOF\n\n");
  };
  fflush(file);
};

void path_display(int iter, real* __restrict__ mep, real* __restrict__ grad_f
, realp f, real res, real constres, real alpha, realp last_f, real last_res) {
  if(!post_optimization && (res<100*epsilon || iter>max_iter/2)) {
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
    fprintf(stderr, "%+.2" RF "g", RT(constres));
    fprintf(stderr, " " COLOR_YELLOW "A" COLOR_RESET "%" RF "g", RT(alpha));
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
      linear_comb(3*SIZE,1-loc[p],mep+idx[p]*3*SIZE,loc[p],mep+(idx[p]+1)*3*SIZE,mep+p*3*SIZE);
      idx[p]=-idx[p];
    };
  } else {
    //fprintf(stderr, COLOR_FAINT "Failed to interpolate inplace\n" COLOR_RESET);
    // copying old images to a buffer
    real* buf=(real*)malloc(sizeof(real)*3*SIZE*(to-from+1)); assert(buf);
    copy_vector(3*SIZE*(to-from+1), mep+3*SIZE*from, buf);
    for(int p=to-1; p>from; p--) {
      assert(p>from && p<to);
      assert(idx[p]>=from && idx[p+1]<=to);
      linear_comb(3*SIZE,1-loc[p],buf+(idx[p]-from)*3*SIZE,loc[p],buf+(idx[p]-from+1)*3*SIZE,mep+p*3*SIZE);
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
  energy_evaluate(mep);
  path_equilize(mep);
  real sum=0;
  for(int p=0; p<sizep; p++) sum+=normalize(mep+3*SIZE*p);  
  //energy_evaluate(mep);
  return sum;
};

void path_tangent_rec(const real* __restrict__ mep, real* __restrict__ grad, int from, int to) {
  real* u=ralloc(3*SIZE);
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
      three_point_tangent_stable(energy[p-1],energy[p],energy[p+1],mep+3*SIZE*(p-1), mep+p*3*SIZE, mep+3*SIZE*(p+1), u);
      //three_point_tangent(mep+3*SIZE*(p-1), mep+p*3*SIZE, mep+3*SIZE*(p+1), u);
      real normu=rsqrt(normsq(3*SIZE,u));
      real proj=(l[p-from]/l[C-1]*q[C-1]-q[p-from]);// inflation compensation
      mult_add(3*SIZE, -proj/normu, u, grad+p*3*SIZE);
    } else if(
        post_optimization 
        && ( (p==from && from>0 && energy[from]>energy[to]) 
           //||(p==to && to<sizep-1 && energy[from]<energy[to])
           )
      ) 
    {
      three_point_tangent_stable(energy[p-1],energy[p],energy[p+1],mep+3*SIZE*(p-1), mep+p*3*SIZE, mep+3*SIZE*(p+1), u);
      //three_point_tangent(mep+3*SIZE*(p-1), mep+p*3*SIZE, mep+3*SIZE*(p+1), u);
      real proj=dot(3*SIZE,u,grad+p*3*SIZE)/dot(3*SIZE,u,u);
      mult_add(3*SIZE, -2*proj, u, grad+p*3*SIZE);
    };
  };
  free(u); free(q); free(l);
}

// INVALID
void path_tangent_rohart(const real* __restrict__ mep, real* __restrict__ grad) {
  real* u=ralloc(3*SIZE);
  for(int p=0; p<sizep; p++) {
    //real l1=normsq(3*SIZE, grad+3*SIZE*p);
    project_to_tangent(mep+3*SIZE*p,grad+3*SIZE*p);
    //real l2=normsq(3*SIZE, grad+3*SIZE*p);
    //real l3=0, l4=0;
    if(p>0 && p<sizep-1) {
      three_point_tangent_stable(energy[p-1],energy[p],energy[p+1],mep+3*SIZE*(p-1), mep+p*3*SIZE, mep+3*SIZE*(p+1), u);
      normalize(u);
      project_to_tangent(u,grad+3*SIZE*p);
      //l3=normsq(3*SIZE, grad+3*SIZE*p);
      //project_to_tangent(mep+3*SIZE*p,grad+3*SIZE*p);
      //l4=normsq(3*SIZE, grad+3*SIZE*p);
    };
    //fprintf(stderr, "  %" RF "g  %" RF "g  %" RF "g  %" RF "g\n", l1, l2, l3, l4);
  };
  free(u); 
}

void path_tangent_neb(const real* __restrict__ mep, real* __restrict__ grad) {
  real* u=ralloc(3*SIZE);
  real *g1=NULL, *g2=NULL, *g3=NULL; 
  if(remove_zero_modes) {
    g1=ralloc(3*SIZE);
    g2=ralloc(3*SIZE);
    g3=ralloc(3*SIZE);
  };
  // Find index f of node with
  realp max=-INFINITY; int f=-1;
  for(int p=sizep-1; p>=0; p--) if(energy[p]>max) { max=energy[p]; f=p; };
  if(f==0 || f==sizep-1) f=-1;
  for(int p=0; p<sizep; p++) {
    project_to_tangent(mep+3*SIZE*p,grad+3*SIZE*p);
    if(p>0 && p<sizep-1) {
      //if(post_optimization && energy[p]<energy[p-1] && energy[p]<energy[p+1]) {
      //} else 
      if(post_optimization && (p==f || (energy[p]>=energy[p-1] && energy[p]>=energy[p+1] && !single_maximum))) {
        three_point_tangent_stable(energy[p-1],energy[p],energy[p+1],mep+3*SIZE*(p-1), mep+p*3*SIZE, mep+3*SIZE*(p+1), u);
        //three_point_tangent(mep+3*SIZE*(p-1), mep+p*3*SIZE, mep+3*SIZE*(p+1), u);
        real proj=dot(3*SIZE,u,grad+p*3*SIZE)/dot(3*SIZE,u,u);
        mult_add(3*SIZE, -2*proj, u, grad+p*3*SIZE);
      } else {
        //three_point_tangent(mep+3*SIZE*(p-1), mep+p*3*SIZE, mep+3*SIZE*(p+1), u);
        if(remove_zero_modes) {
          three_point_tangent_stable(energy[p-1],energy[p],energy[p+1],mep+3*SIZE*(p-1), mep+p*3*SIZE, mep+3*SIZE*(p+1), u);
          group_generator(mep+3*SIZE*p, 0, g1);
          group_generator(mep+3*SIZE*p, 1, g2);
          group_generator(mep+3*SIZE*p, 2, g3);
          real* basis[5]={g1,g2,g3,u,grad+p*3*SIZE};
          gram_schmidt(3*SIZE, 5, basis);
        } else {
          three_point_tangent_stable(energy[p-1],energy[p],energy[p+1],mep+3*SIZE*(p-1), mep+p*3*SIZE, mep+3*SIZE*(p+1), u);
          real* basis[2]={u,grad+p*3*SIZE};
          gram_schmidt(3*SIZE, 2, basis);
          //real normu=rsqrt(normsq(3*SIZE,u));
          //real proj=dot(3*SIZE,u,grad+p*3*SIZE)/normu;
          //mult_add(3*SIZE, -proj/normu, u, grad+p*3*SIZE);
        };
      };
    };
  };
  free(u); if(g1)free(g1); if(g2)free(g2); if(g3)free(g3);
}

void path_tangent(const real* __restrict__ mep, real* __restrict__ grad) {
  // Rohart way
  //path_tangent_rohart(mep, grad); return;
  path_tangent_neb(mep, grad);
}  

void projected_path_gradient(const real* __restrict__ arg, real* __restrict__ out, realp* __restrict__ E) {
  path_gradient(arg, out, E);
  path_tangent(arg, out);
};

int path_steepest_descent(real* __restrict__ path, int mode, 
  real mode_param, real epsilon, int max_iter) 
{
  real updated_param=mode_param;
  //if(mode==SDM_CONSTANT) updated_param=mode_param/sizep;
  return steepest_descend(
    3*SIZE*sizep, (real*)path, 
    projected_path_gradient,
    mode, 
    updated_param, epsilon, max_iter,
    path_display, 
    path_normalize
  );
};
