#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "vector.h"
#include "optim.h"
#include "debug.h"

// Compute minimum doing gradient descend.
// Minimizing function f(x) is of special form:
//   f(x)=<x|Hx>/2+<x|B>.
// Hessian part H is given by operator Q:x->H(x) with signature
//   void Q(const float* x, float* y)
// Linear part is given as an operator L:y->y+B
//   void L(float* y)
// The minimum is obtained as the limit of the sequence
// x(k+1)=x(k)-alpha(k)*grad_f(x(k))
// where grad_f(x) is the gradient of the minimizing function f(x):
//   grad_f(x)=Hx+B.
// Sequence alpha(k) is determined by mode:
//   SDM_CONSTANT : alpha(k)=mode_param
//   SDM_INERTIAL : alpha(k)=alpha(k-1)+mode_param*|grad_f(x(k-1))| if f(x(k))<f(x(k-1))
//                  alpha(k)=mode_param*|grad_f(x(k-1))| otherwise
// Initial approxiamtion is x(0)=a.
// Result is stored in the bufer 'a'.
// Length of the vector 'a' is 'n'.
// Progress can be displayed by function display:
//   void display(int d, int n, float* a, float* grad_f, float f, float res)
// where 
//   iter : iteration counter,
//   n : vector length,
//   a : current aproximation x(k),
//   grad_f : gradient of f(x) at a=x(k),
//   f : value of f(x(k)),
//   res : residual = |grad_f|.
//   alpha : last multiplier
// The method can be used for optimization with constrains,in the case
// routine 'P' must project approximation 'a' to manifold of solutions
// satisfying constrains, and 'T' must project vector field 't'
// to tangent space at the point 'a':
//   void P(int n, float* a)
//   void T(int n, const float* a, float* t)
// Return 0 on success (residual is less then epsilon)
//        1 if maximum number of iteration max_iter is reached
//        2 if machine precision is reached for 'f' evaluation
int steepest_descend(
	int n, float* a, 
	void (*Q)(const float* x, float* y),
	void (*L)(float* y),
	int mode, float mode_param, float epsilon, int max_iter,
	void (*display)(int iter, float* a, float* grad_f, float f, float res, float alpha),
  void (*P)(float* a),
  void (*T)(const float* a, float* t)
	)
{
  int iter;
  float alpha;
  assert(Q && L); assert((P && T) || (!P && !T));
  switch(mode) {
  	case SDM_CONSTANT: alpha=mode_param; break;
	  case SDM_INERTIAL: alpha=NAN; break;
    case SDM_PROGR: alpha=1e-4; break;
    case SDM_CAUCHY: alpha=0; break;
	  default: assert(1);
  };
  float* grad=(float*)malloc(sizeof(float)*n); assert(grad);
  float* hess=(float*)malloc(sizeof(float)*n); assert(hess);
  int status=1;
  float last_res=0;
  float last_f=-INFINITY;
  if(P) P(a);
  for(iter=0; iter<max_iter; iter++) {
    float next_alpha;
  	Q(a,grad);
  	float f=-dot(n,a,grad);
  	L(grad);
  	f+=dot(n,a,grad);
  	// grad contains Hx+B
  	// f value of f(a)
    //ill conditioned quadratic optimiation
    //fprintf(stderr,"res: %g\n",sqrtf(normsq(n,grad)));
    if(T) T(a,grad);
    float res2=normsq(n,grad); float res=sqrtf(res2);
    assert(!isnan(res));
    switch(mode) {
      case SDM_INERTIAL: 
        if(f<=last_f) alpha+=mode_param*last_res; else alpha=mode_param*res;
        break;
      case SDM_PROGR: 
        if(f<=last_f) alpha*=mode_param; else alpha/=2*mode_param;
        break;
      case SDM_CAUCHY: 
        Q(grad, hess); 
        next_alpha=fabs(res2/dot(n,grad,hess));
        if(next_alpha>mode_param) next_alpha=mode_param;
        break; 
    };
    if(display) display(iter,a,grad,f,res,alpha);
    if(res<epsilon) { status=0; break; }; // If solution is found
    //if(last_f==f) { status=2; break; }; // Iterations stop changing
    // Calculation next aproximation
    mult_sub(n,alpha,grad,a);
    if(P) P(a);
    last_res=res;
    last_f=f;
    if(mode==SDM_CAUCHY) alpha=next_alpha;
  };
  free(grad); free(hess);
  return status;
};

// Solve minimization problem
//   argmin <x|Hx>/2+<x|B> 
// subjected to constrains
//   <x|P_j x>-1=0
// Lagrange function is constructed
//   L(x,u)=<x|Hx>/2+mu*<x|x>/2+<x|B>+sum_l u_l*(<x|P_jx>-1)/2
// where mu - barrier multiplier, tends to 0 as x tends to minimum
//       u - Lagrange multipliers.
// Minimum is calculated as the limit of sequence x_k where
//   x_0=x0 - initial approximation
//   x_k - solution of auxilliary problem
//         grad_{x,u} L(x,u)=0
//         obtained by conjugated residual
// grad_x L(x,u)=Hx+mu*x+B+sum_l u_l*P_j x
// Hessian_x L(x,u)=H+mu+sum_l u_l*P_j
// If during solution of the auxilliary problem
//   <y|Hessian_x L|y> < 0 for some y
// then mu is increased mu+=<y|Hessian_x L|y>+correction
// otherwise mu/=ratio.
// dim(x)=N, dim(u)=M
// Hessian part H is given by operator Q:x->H(x) with signature
//   void Q(const float* x, float* y)
// Linear part is given as an operator L:y->y+B
//   void L(float* y)
// Constrains gradient is given by 
//   C:x->(<x|P_j x>/2-1/2)_j
//   D:x,u,r->r+sum_l u_l P_j x
//   P:x,y->(<x|P_j y>)_l


int lagrange_conjugate(
  int N, int M, float* x0, 
  void (*Q)(const float* x, float* y),
  void (*L)(float* y),
  int mode, float mode_param, float epsilon, int max_iter,
  void (*display)(int iter, float* a, float* grad_f, float f, float res, float alpha),
  void (*C)(const float* x, float* r),
  void (*D)(const float* x, const float* u, float* r),
  void (*P)(const float* x, const float* y, float* r)
  )
{
  float min_mu=0.1; // Minimum required bpttpm of spectrum of Q+mu
  float dot_precision=1.2e-07*(N+M)*10; // Precision of inner product
  // allocate mempry
  float* u=(float*)malloc(sizeof(float)*M); assert(u);
  float* gradx=(float*)malloc(sizeof(float)*N); assert(gradx);
  float* gradu=(float*)malloc(sizeof(float)*M); assert(gradu);
  float* conjx=(float*)malloc(sizeof(float)*N); assert(conjx);
  float* conju=(float*)malloc(sizeof(float)*M); assert(conju);
  float* xn=(float*)malloc(sizeof(float)*N); assert(xn);
  float* un=(float*)malloc(sizeof(float)*M); assert(un);
  float* hgradx=(float*)malloc(sizeof(float)*N); assert(hgradx);
  float* hgradu=(float*)malloc(sizeof(float)*M); assert(hgradu);
  float* hconjx=(float*)malloc(sizeof(float)*N); assert(hconjx);
  float* hconju=(float*)malloc(sizeof(float)*M); assert(hconju);
  // initialize state
  int status=1;  
  float mu=0;
  zero_vector(M, u);

  for(int iter=0; iter<max_iter;) {
    // Calculate gradient at the point
    restart: Q(x0,gradx); iter++;
    float f=-dot(N,x0,gradx);
    L(gradx);
    f+=dot(N,x0,gradx);
    D(x0,u,gradx);
    mult_add(N,mu,x0,gradx);
    C(x0,gradu); 
    // Calculate norm of residual
    float resx=normsq(N,gradx), resu=normsq(M,gradu);
    float res=sqrtf(resx+resu);
    float last_res=res;
    // Diplay progress
    #ifdef DEBUG
      fprintf(stderr,COLOR_RED"%d: res=%g %+g mu=%g\n"COLOR_RESET,iter,sqrtf(resx),sqrtf(resu),mu);
      fprintf(stderr, COLOR_YELLOW);
      for(int j=0;j<N;j++) fprintf(stderr,"%g ",x0[j]); fprintf(stderr,":");
      for(int j=0;j<M;j++) fprintf(stderr," %g",u[j]); 
      fprintf(stderr,COLOR_RESET"\n");
    #endif
    // Residula is negative gradient
    negate_inplace(N,gradx); negate_inplace(M,gradu);
    // Reporting result
    if(display) display(iter,x0,gradx,f,res,mu);
    if(res<epsilon) { 
      #ifdef DEBUG
        fprintf(stderr,COLOR_RED"%d: Converged %.3g < %.3g\n"COLOR_RESET,iter,res,epsilon);
      #endif    
      status=0; break; 
    }; // If solution is found
    // Solving auxilliary problem with respect to xn,un
    // [...]*[xn;un]=[gradx;gradu]
    // initial guess
    zero_vector(N, xn); zero_vector(M,un);
    // computing quadratic part
    //h=[A+2*lambda+shift,2*x;2*x',0];
    iter++;
    Q(gradx,hgradx); 
    D(gradx,u,hgradx);
    mult_add(N,mu,gradx,hgradx);
    float positive=dot(N,gradx,hgradx)/resx; // Must be positive if form is positive
    // If the quadratic form is not positive, update shift 'mu'
    if(positive<min_mu) { 
      #ifdef DEBUG      
        fprintf(stderr,COLOR_GREEN"  %d: Hessian is small: %.3e < %.3e\n"COLOR_RESET,iter,positive,min_mu);
      #endif      
      mu+=2*min_mu-positive;  
      goto restart; 
    };
    float low_boundary=positive; // Estimate of bottom of Hessian spectrum
    D(x0,gradu,hgradx); P(x0,gradx,hgradu);
    //zero_vector(M,hgradu);
    // choosing initial conjugate direction
    copy_vector(N, gradx, conjx); copy_vector(M, gradu, conju);
    copy_vector(N, hgradx, hconjx); copy_vector(M, hgradu, hconju);
    // Initialization of conjugate residual
    float q=dot(N,gradx,hgradx)+dot(M,gradu,hgradu); 
    while(iter<max_iter) {
      if(fabs(q)<dot_precision) { 
        #ifdef DEBUG        
          fprintf(stderr,COLOR_GREEN"  %d: Small q: %.3e < %.3e\n"COLOR_RESET,iter,q,dot_precision);
        #endif        
        mu+=min_mu; 
        goto restart; 
      };
      float hconj_normsq=normsq(N,hconjx)+normsq(M,hconju);
      #ifdef DEBUG 
        // Debug
        fprintf(stderr, COLOR_YELLOW"  X:"COLOR_RESET); for(int j=0;j<N;j++) fprintf(stderr," %g",xn[j]); fprintf(stderr,":"); for(int j=0;j<M;j++) fprintf(stderr," %g",un[j]); fprintf(stderr,"\n");
        fprintf(stderr, COLOR_YELLOW"  G:"COLOR_RESET); for(int j=0;j<N;j++) fprintf(stderr," %g",gradx[j]); fprintf(stderr,":"); for(int j=0;j<M;j++) fprintf(stderr," %g",gradu[j]); fprintf(stderr,"\n");
        fprintf(stderr, COLOR_YELLOW"  C:"COLOR_RESET); for(int j=0;j<N;j++) fprintf(stderr," %g",conjx[j]); fprintf(stderr,":"); for(int j=0;j<M;j++) fprintf(stderr," %g",conju[j]); fprintf(stderr,"\n");
        // <grad|A grad>=<grad|A conj>
        float hgradhconj=dot(N,hgradx,hconjx)+dot(M,hgradu,hconju);
        fprintf(stderr, "Biorthogonality %g\n",hgradhconj-hconj_normsq);
        // Self-adjointness test
        float gradhconj=dot(N,gradx,hconjx)+dot(M,gradu,hconju);
        float conjhgrad=dot(N,hgradx,conjx)+dot(M,hgradu,conju);
        fprintf(stderr, "Self-adjointness %g\n",gradhconj-conjhgrad);
      #endif

      float alpha=q/hconj_normsq;
      mult_add(N, alpha, conjx, xn); mult_add(M, alpha, conju, un);
      mult_sub(N, alpha, hconjx, gradx); mult_add(M, alpha, hconju, gradu);
      resx=normsq(N,gradx); resu=normsq(M,gradu); res=sqrtf(resx+resu);
      assert(!isnan(res));
      if(res<epsilon) { // Auxilliary problem is solved
        fprintf(stderr,COLOR_GREEN"  %d: Residual is small: R %.3e %+.3e \n"COLOR_RESET,iter,sqrtf(resx),sqrtf(resu));
        break; 
      };
      if(res>100*last_res) { // Unstability detected ???
        fprintf(stderr,COLOR_GREEN"  %d: Residual increases: R %.3e %+.3e \n"COLOR_RESET,iter,sqrtf(resx),sqrtf(resu));
        break; 
      };
      last_res=res;
      // Aplying quadratic form to residual
      iter++;
      Q(gradx,hgradx); 
      D(gradx,u,hgradx);
      mult_add(N,mu,gradx,hgradx);
      float positive=dot(N,gradx,hgradx)/resx; // Must be positive if form is positive
      // If the quadratic form is not positive, update shift 'mu'
      if(positive<min_mu) { 
        fprintf(stderr,COLOR_GREEN"  %d: Hessian is small: %.3e < %.3e\n"COLOR_RESET,iter,positive,min_mu);
        mu+=2*min_mu-positive; 
        goto restart; 
      };
      if(positive<low_boundary) low_boundary=positive;
      D(x0,gradu,hgradx); P(x0,gradx,hgradu);
      //zero_vector(M,hgradu);
      #ifdef DEBUG 
        // Self-adjointness test 2
        gradhconj=dot(N,gradx,hconjx)+dot(M,gradu,hconju);
        conjhgrad=dot(N,hgradx,conjx)+dot(M,hgradu,conju);
        fprintf(stderr, "Self-adjointness (2) %g\n",gradhconj-conjhgrad);
      #endif
      // Calculating projections
      float qn=dot(N,gradx,hgradx)+dot(M,gradu,hgradu);
      float beta=qn/q; q=qn;
      #ifdef DEBUG 
        // Checking orthogoaity <A conj_k|A conj_{k+1}>=0
        float hgradnexthconj=dot(N,hgradx,hconjx)+dot(M,hgradu,hconju);
        fprintf(stderr, "Orthogonality %g\n",hgradnexthconj+beta*hconj_normsq);
      #endif      
      // Update conjugate direction
      add_mult(N, gradx, beta, conjx); add_mult(M, gradu, beta, conju); 
      add_mult(N, hgradx, beta, hconjx); add_mult(M, hgradu, beta, hconju); 
      fprintf(stderr,COLOR_GREEN"  %d: R %.3e %+.3e A %.3e B %.3e P %.3e q %.3e\n"COLOR_RESET,iter,sqrtf(resx),sqrtf(resu),alpha,beta,positive,q);
    };
    #ifdef DEBUG 
      fprintf(stderr,COLOR_BLUE"  %d: Hessian bottom %.3e\n"COLOR_RESET,iter,mu-low_boundary);
    #endif    
    mu=(mu-low_boundary)/2+min_mu;
    sub_inplace(N,xn,x0); sub_inplace(M,un,u); 
  };
  //stop:{};
  free(u); free(xn); free(un); 
  free(gradx); free(gradu); free(conjx); free(conju);
  free(hgradx); free(hgradu); free(hconjx); free(hconju);
  return status; 
}