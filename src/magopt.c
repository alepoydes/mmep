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

void show_grad(const real* arg, const real* grad) {
  forall(u,x,y,z) {
    int i=INDEX(u,x,y,z);
    if(ISACTIVE(active, i)) fprintf(stderr, COLOR_BOLD);
    fprintf(stderr, COLOR_BLUE "%d %d %d %d: " COLOR_RESET, u, x, y, z);
    const real* a=arg+3*i;
    fprintf(stderr, COLOR_GREEN "%" RF "g %" RF "g %" RF "g " COLOR_RESET, RT(a[0]), RT(a[1]), RT(a[2]));
    const real* g=grad+3*i;
    fprintf(stderr, "%" RF "g %" RF "g %" RF "g ", RT(g[0]), RT(g[1]), RT(g[2]));
    //fprintf(stderr, COLOR_RED "%" RF "g ", RT(dot3(g,a)));
    //fprintf(stderr, COLOR_YELLOW "%" RF "g ", RT(dot3(a,a)-1));
    fprintf(stderr, "\n");
  };
};

void skyrmion_display(int iter, real* __restrict__ a, real* __restrict__ grad_f, realp f, real res, 
real constres, real alpha, realp last_f, real last_res, real err) {
  // show progress
  static realp prev_f=NAN; static real prev_res=NAN; 
  if(iter==1) { prev_f=f; prev_res=res; };
  if(iter%debug_every==0 || iter<0) {
  	if(iter!=0) fprintf(stderr, "%6d", abs(iter));
    else fprintf(stderr, "%6s", "");
    fprintf(stderr, " " COLOR_YELLOW "E" COLOR_RESET);
    watch_number(f,isnan(last_f)?f:debug_every==1?last_f:prev_f,16);
    fprintf(stderr, " " COLOR_YELLOW "R" COLOR_RESET);
    watch_number(res,isnan(last_res)?res:debug_every==1?last_res:prev_res,16);
    fprintf(stderr, "%+.2" RF "g", RT(constres));
    fprintf(stderr, " " COLOR_YELLOW "A" COLOR_RESET "%" RF "g", RT(alpha));
    //fprintf(stderr, " " COLOR_YELLOW "d" COLOR_RESET "%.2g", pow(last_res/res,1./alpha));
    fprintf(stderr, " " COLOR_YELLOW "#" COLOR_RESET "%d", number_of_used);    
    if(!isnan(err)) fprintf(stderr, " " "%.0" RF "f" COLOR_YELLOW "%%" COLOR_RESET , RT(100*err));
    fprintf(stderr, "\n");
  	if(debug_plot && iter!=0) 
      plot_field3(stdout,a);
      //plot_field3(stdout,grad_f);
    prev_f=f; prev_res=res;
  };
  //show_grad(a,grad_f);
};

void update_active(const real* a, real* grad_f, real* f) {
  static int iter=0;
  if(active_iterations<1 || iter++%active_iterations) return;
  // activate all
  COPYMASK(all_active, active);
  // compute gradient in all points
  projected_gradient(a, grad_f, f);
  // make fast spins active
  number_of_used=deactivate_slow(grad_f, 1, active_threshold);
  // recompute energy and gradient
  projected_gradient(a, grad_f, f);
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

void projected_gradient2(const real* __restrict__ arg, real* __restrict__ grad, realp* __restrict__ energy) {
    skyrmion_gradient(arg, grad, energy);
    real t=rsqrt(normsq(3*SIZE, grad));
    project_to_tangent(arg, grad);
    real p=rsqrt(normsq(3*SIZE, grad));
    fprintf(stderr, "proj grad rat %" RF "g%%\n", RT(p/t*100));
};

/*
void projected_gradient_activate(const real* __restrict__ arg, real* __restrict__ grad, realp* __restrict__ energy) {
  projected_gradient(arg, grad, energy);
  activate_fast_and_adjacent(grad, 1, active_threshold);
};
*/

int skyrmion_steepest_descent(real* __restrict__ x, 
int mode, real mode_param, real epsilon, int max_iter) 
{
  COPYMASK(all_active, active);
  number_of_used=number_of_active;
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
    //projected_gradient_activate,
    projected_gradient,
		mode, mode_param, epsilon, max_iter,
		skyrmion_display, 
		quasynorm,
    update_active
	);

  COPYMASK(all_active, active);
  number_of_used=number_of_active;  
};

int skyrmion_better_steepest_descent(
real* __restrict__ x, real epsilon, int max_iter
) 
{
  // init mask
  COPYMASK(all_active, active);
  number_of_used=number_of_active;

  // allocate memory
  real* grad=ralloc(3*SIZE);
  real* dir=ralloc(3*SIZE);
  real* lambda=ralloc(SIZE);

  // initialize state
  realp last_f=-INFINITY, f=NAN;
  real res=NAN, last_res=0;
  real constres=NAN;
  int status=1;
  int count_constant_f=0;

  constres=normalize(x);
  // start iterations
  int iter; for(iter=1; iter<max_iter; iter++) {
    hamiltonian_hessian(x, grad);
    f=skyrmion_energy_given_hessian(x, grad);
    for(int n=0; n<SIZE; n++) {
      if(!ISACTIVE(active,n)) continue;
      real l=lambda[n]=dot3(grad+3*n, x+3*n);
      for3(j) grad[3*n+j]-=l*x[3*n+j];
    };
    real res2=normsq(3*SIZE, grad);

    // check exit condition
    res=rsqrt(res2/SIZE); assert(!isnan(res));
    if(res+constres<epsilon) { // If solution is found
      status=0; break; 
    };
    if(stop_signal>0) {
      fprintf(stderr, COLOR_YELLOW "Optimization aborted\n" COLOR_RESET);
      stop_signal=0;
      break;
    }; 
    if(last_f==f) {
      if(count_constant_f++>10) { 
        fprintf(stderr, "Maximum precision is reached\n");
        status=2; break; 
      };
    };

    // compute step size
    hamiltonian_hessian(grad, dir);
    real pAp=0;
    for(int n=0; n<SIZE; n++) {
      if(!ISACTIVE(active,n)) continue;
      tangent3(x+3*n,dir+3*n);
      for3(j) dir[3*n+j]-=lambda[n]*grad[3*n+j];
      pAp+=dot3(grad+3*n, dir+3*n);
      //fprintf(stderr, COLOR_BLUE "%d" COLOR_RESET " %"RF"g\n",n,RT(pAp));
    };
    real alpha=rabs(res2/pAp);
    // compute next approximatiob
    constres=0;
    for(int n=0; n<SIZE; n++) {
      if(!ISACTIVE(active,n)) continue;
      for3(j) x[3*n+j]-=alpha*grad[3*n+j];
      constres+=normalize3(x+3*n);
    };

    // update active spins
    update_active(x, grad, &f);
    // report progress
    skyrmion_display(iter,x,grad,f,res,constres,alpha,last_f,last_res,NAN);

    // save statistics
    last_f=f; last_res=res;
  };
  skyrmion_display(-iter,x,grad,f,res,constres,NAN,last_f,last_res,NAN);
  //free memory
  free(grad); free(dir); free(lambda);
  // reset mask
  COPYMASK(all_active, active);
  number_of_used=number_of_active; 

  return status; 
};


int skyrmion_minimize(real* __restrict__ x, real epsilon, int max_iter) 
{
  switch(single_mode) {
    case 0: 
      return skyrmion_steepest_descent(x, SDM_PROGR, 0.1, epsilon, max_iter);
    case 1:
      return skyrmion_better_steepest_descent(x, epsilon, max_iter);
    default:
      fprintf(stderr, COLOR_RED COLOR_BOLD "Error:" COLOR_RESET "Unknown image minimizer %d\n", single_mode);
      exit(1);
  };
};

//////////////////////////////////////////////////////
// MEP optimization

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
, realp f, real res, real constres, real alpha, realp last_f, real last_res, real err) {
  /*if(!post_optimization && (res<100*epsilon || iter>max_iter/2)) {
      //fprintf(stderr,COLOR_YELLOW "Climbing is on." COLOR_RESET " Residue %" RF "g\n",res);
      post_optimization=1;
  };*/ 
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
    fprintf(stderr, " %" RF "g" COLOR_YELLOW "%%" COLOR_RESET , RT(100*err));
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
  COPYMASK(all_active, active);
  number_of_used=number_of_active;
  real updated_param=mode_param;
  //if(mode==SDM_CONSTANT) updated_param=mode_param/sizep;
  int retcode=steepest_descend(
    3*SIZE*sizep, (real*)path, 
    projected_path_gradient,
    mode, 
    updated_param, epsilon, max_iter,
    path_display, 
    path_normalize,
    NULL
  );
  if(post_optimization) return retcode;
  post_optimization=1;
  retcode=steepest_descend(
    3*SIZE*sizep, (real*)path, 
    projected_path_gradient,
    mode, 
    updated_param, epsilon, max_iter,
    path_display, 
    path_normalize,
    NULL
  );
  return retcode;
};
