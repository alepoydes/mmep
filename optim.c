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
//   void Q(const real* x, real* y)
// Linear part is given as an operator L:y->y+B
//   void L(real* y)
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
//   void display(int d, int n, real* a, real* grad_f, real f, real res)
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
//   void P(int n, real* a)
//   void T(int n, const real* a, real* t)
// Return 0 on success (residual is less then epsilon)
//        1 if maximum number of iteration max_iter is reached
//        2 if machine precision is reached for 'f' evaluation
int steepest_descend(
	int n, real* a, 
	void (*Q)(const real* x, real* y),
	void (*L)(real* y),
	int mode, real mode_param, real epsilon, int max_iter,
	void (*display)(int iter, real* a, real* grad_f, real f, real res, real alpha),
  void (*P)(real* a),
  void (*T)(const real* a, real* t)
	)
{
  int iter;
  real alpha;
  real factor=rsqrt(n);
  assert(Q && L); assert((P && T) || (!P && !T));
  switch(mode) {
  	case SDM_CONSTANT: alpha=mode_param; break;
	  case SDM_INERTIAL: alpha=NAN; break;
    case SDM_PROGR: alpha=mode_param; break;
    case SDM_CAUCHY: alpha=0; break;
	  default: assert(1);
  };
  real* grad=(real*)malloc(sizeof(real)*n); assert(grad);
  real* hess=(real*)malloc(sizeof(real)*n); assert(hess);
  int status=1;
  real last_res=0;
  real last_f=-INFINITY;
  if(P) P(a);
  for(iter=0; iter<max_iter; iter++) {
    real next_alpha;
  	Q(a,grad);
  	real f=-dot(n,a,grad)/2;
  	L(grad);
  	f+=dot(n,a,grad);
  	// grad contains Hx+B
  	// f value of f(a)
    //ill conditioned quadratic optimiation
    //fprintf(stderr,"res: %g\n",rsqrt(normsq(n,grad)));
    assert(!isnan(normsq(n,grad)));
    if(T) T(a,grad);
    real res2=normsq(n,grad); real res=rsqrt(res2);
    assert(!isnan(res));
    switch(mode) {
      case SDM_INERTIAL: 
        if(f<=last_f) alpha+=mode_param*last_res/factor; else alpha=mode_param*res/factor;
        break;
      case SDM_PROGR: 
        //if(f<=last_f) alpha*=1.1; else alpha=mode_param;
        if(res<=last_res) alpha*=1.1; else alpha=mode_param;
        break;
      case SDM_CAUCHY: 
        Q(grad, hess); 
        next_alpha=rabs(res2/dot(n,grad,hess));
        if(next_alpha>mode_param) next_alpha=mode_param;
        break; 
    };
    if(display) display(iter,a,grad,f,res/factor,alpha);
    if(res<factor*epsilon) { status=0; break; }; // If solution is found
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
//   void Q(const real* x, real* y)
// Linear part is given as an operator L:y->y+B
//   void L(real* y)
// Constrains gradient is given by 
//   C:x->(<x|P_j x>/2-1/2)_j
//   D:x,u,r->r+sum_l u_l P_j x
//   P:x,y->(<x|P_j y>)_l


int lagrange_conjugate(
  int N, int M, real* x0, 
  void (*Q)(const real* x, real* y),
  void (*L)(real* y),
  int mode, real mode_param, real epsilon, int max_iter,
  void (*display)(int iter, real* a, real* grad_f, real f, real res, real alpha),
  void (*C)(const real* x, real* r),
  void (*D)(const real* x, const real* u, real* r),
  void (*P)(const real* x, const real* y, real* r)
  )
{
  #define MU_INC mu-positive+mode_param
  real min_mu=1e-4; // Minimum required bpttpm of spectrum of Q+mu
  real dot_precision=EPSILON*10/**(N+M)*/; // Precision of inner product
  // allocate mempry
  real* u=(real*)malloc(sizeof(real)*M); assert(u);
  real* gradx=(real*)malloc(sizeof(real)*N); assert(gradx);
  real* gradu=(real*)malloc(sizeof(real)*M); assert(gradu);
  real* conjx=(real*)malloc(sizeof(real)*N); assert(conjx);
  real* conju=(real*)malloc(sizeof(real)*M); assert(conju);
  real* xn=(real*)malloc(sizeof(real)*N); assert(xn);
  real* un=(real*)malloc(sizeof(real)*M); assert(un);
  real* hgradx=(real*)malloc(sizeof(real)*N); assert(hgradx);
  real* hgradu=(real*)malloc(sizeof(real)*M); assert(hgradu);
  real* hconjx=(real*)malloc(sizeof(real)*N); assert(hconjx);
  real* hconju=(real*)malloc(sizeof(real)*M); assert(hconju);
  // initialize state
  int status=1;  
  real mu=0;
  zero_vector(M, u);

  for(int iter=0; iter<max_iter;) {
    // Calculate gradient at the point
    restart: Q(x0,gradx); iter++;
    real f=-dot(N,x0,gradx)/2;
    L(gradx);
    f+=dot(N,x0,gradx);
    D(x0,u,gradx);
    //real rest=normsq(N,gradx);
    mult_add(N,mu,x0,gradx);
    C(x0,gradu); 
    // Calculate norm of residual
    real resx=normsq(N,gradx), resu=normsq(M,gradu);
    real res=rsqrt(resx+resu);
    real last_res=res;
    // Diplay progress
    #ifdef DEBUG
      fprintf(stderr,COLOR_RED"%d: res=%"RF"g %+"RF"g mu=%"RF"g\n"COLOR_RESET,iter,rsqrt(resx),rsqrt(resu),mu);
      if(N+M<8) {
        fprintf(stderr, COLOR_YELLOW);
        for(int j=0;j<N;j++) fprintf(stderr,"%"RF"g ",x0[j]); fprintf(stderr,":");
        for(int j=0;j<M;j++) fprintf(stderr," %"RF"g",u[j]); 
        fprintf(stderr,COLOR_RESET"\n");
      };
    #endif
    // Residula is negative gradient
    negate_inplace(N,gradx); negate_inplace(M,gradu);
    // Reporting result
    if(display) display(iter,x0,gradx,f,res,mu);
    if(res<epsilon) { 
      #ifdef DEBUG
        fprintf(stderr,COLOR_RED"%d: Converged %.3"RF"g < %.3"RF"g\n"COLOR_RESET,iter,res,epsilon);
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
    //vector_copy(N,gradx,hgradx);
    real positive=dot(N,gradx,hgradx)/resx; // Must be positive if form is positive
    // If the quadratic form is not positive, update shift 'mu'
    if(positive<min_mu) { 
      #ifdef DEBUG      
        fprintf(stderr,COLOR_GREEN"  %d: Hessian is small: %.3"RF"e < %.3"RF"e\n"COLOR_RESET,iter,positive,min_mu);
      #endif      
      mu=MU_INC;  
      goto restart; 
    };
    real low_boundary=positive; // Estimate of bottom of Hessian spectrum
    D(x0,gradu,hgradx); P(x0,gradx,hgradu);
    //vector_copy(M,gradu,hgradu);
    // choosing initial conjugate direction
    copy_vector(N, gradx, conjx); copy_vector(M, gradu, conju);
    copy_vector(N, hgradx, hconjx); copy_vector(M, hgradu, hconju);
    // Initialization of conjugate residual
    real q=dot(N,gradx,hgradx)+dot(M,gradu,hgradu); 
    while(iter<max_iter) {
      if(rabs(q)<dot_precision) { 
        #ifdef DEBUG        
          fprintf(stderr,COLOR_GREEN"  %d: Small q: %.3"RF"e < %.3"RF"e\n"COLOR_RESET,iter,q,dot_precision);
        #endif        
        //mu+=min_mu; goto restart; 
        break;
      };
      real hconj_normsq=normsq(N,hconjx)+normsq(M,hconju);
      #ifdef DEBUG 
        // <grad|A grad>=<grad|A conj>
        real hgradhconj=dot(N,hgradx,hconjx)+dot(M,hgradu,hconju);
        assert(rabs(hgradhconj-hconj_normsq)<1e-4);
        //fprintf(stderr, "Biorthogonality %g\n",hgradhconj-hconj_normsq);
        // Self-adjointness test
        real gradhconj=dot(N,gradx,hconjx)+dot(M,gradu,hconju);
        real conjhgrad=dot(N,hgradx,conjx)+dot(M,hgradu,conju);
        assert(rabs(gradhconj-conjhgrad)<1e-4);
        //fprintf(stderr, "Self-adjointness %g\n",gradhconj-conjhgrad);
      #endif

      real alpha=q/hconj_normsq;
      mult_add(N, alpha, conjx, xn); mult_add(M, alpha, conju, un);
      mult_sub(N, alpha, hconjx, gradx); mult_sub(M, alpha, hconju, gradu);

      #ifdef DEBUG 
        // Debug
        if(N+M<8) {
          fprintf(stderr, COLOR_YELLOW"  X:"COLOR_RESET); for(int j=0;j<N;j++) fprintf(stderr," %"RF"g",xn[j]); fprintf(stderr,":"); for(int j=0;j<M;j++) fprintf(stderr," %"RF"g",un[j]); fprintf(stderr,"\n");
          fprintf(stderr, COLOR_YELLOW"  G:"COLOR_RESET); for(int j=0;j<N;j++) fprintf(stderr," %"RF"g",gradx[j]); fprintf(stderr,":"); for(int j=0;j<M;j++) fprintf(stderr," %"RF"g",gradu[j]); fprintf(stderr,"\n");
          fprintf(stderr, COLOR_YELLOW"  C:"COLOR_RESET); for(int j=0;j<N;j++) fprintf(stderr," %"RF"g",conjx[j]); fprintf(stderr,":"); for(int j=0;j<M;j++) fprintf(stderr," %"RF"g",conju[j]); fprintf(stderr,"\n");
        };
      #endif      
      resx=normsq(N,gradx); resu=normsq(M,gradu); res=rsqrt(resx+resu);
      assert(!isnan(res));
      if(res<epsilon) { // Auxilliary problem is solved
        #ifdef DEBUG 
          fprintf(stderr,COLOR_GREEN"  %d: Residual is small: R %.3"RF"e %+.3"RF"e \n"COLOR_RESET,iter,rsqrt(resx),rsqrt(resu));
        #endif  
        break; 
      };
      if(res>100*last_res) { // Unstability detected ???
        #ifdef DEBUG 
          fprintf(stderr,COLOR_GREEN"  %d: Residual increases: R %.3"RF"e %+.3"RF"e \n"COLOR_RESET,iter,rsqrt(resx),rsqrt(resu));
        #endif  
        break; 
      };
      last_res=res;
      // Aplying quadratic form to residual
      iter++;
      Q(gradx,hgradx); 
      D(gradx,u,hgradx);
      mult_add(N,mu,gradx,hgradx);
      //vector_copy(N,gradx,hgradx);
      real positive=dot(N,gradx,hgradx)/resx; // Must be positive if form is positive
      // If the quadratic form is not positive, update shift 'mu'
      if(positive<min_mu) { 
        #ifdef DEBUG 
          fprintf(stderr,COLOR_GREEN"  %d: Hessian is small: %.3"RF"e < %.3"RF"e\n"COLOR_RESET,iter,positive,min_mu);
        #endif  
        mu=MU_INC; 
        goto restart; 
      };
      if(positive<low_boundary) low_boundary=positive;
      D(x0,gradu,hgradx); P(x0,gradx,hgradu);
      //vector_copy(M,gradu,hgradu);
      #ifdef DEBUG 
        // Self-adjointness test 2
        gradhconj=dot(N,gradx,hconjx)+dot(M,gradu,hconju);
        conjhgrad=dot(N,hgradx,conjx)+dot(M,hgradu,conju);
        assert(rabs(gradhconj-conjhgrad)<1e-4);
        //fprintf(stderr, "Self-adjointness (2) %g\n",gradhconj-conjhgrad);
      #endif
      // Calculating projections
      real qn=dot(N,gradx,hgradx)+dot(M,gradu,hgradu);
      real beta=qn/q; q=qn;
      #ifdef DEBUG 
        // Checking orthogoaity <A conj_k|A conj_{k+1}>=0
        real hgradnexthconj=dot(N,hgradx,hconjx)+dot(M,hgradu,hconju);
        assert(rabs(hgradnexthconj+beta*hconj_normsq)<1e0);
        //fprintf(stderr, "Orthogonality %g\n",hgradnexthconj+beta*hconj_normsq);
      #endif
      // Update conjugate direction
      add_mult(N, gradx, beta, conjx); add_mult(M, gradu, beta, conju); 
      add_mult(N, hgradx, beta, hconjx); add_mult(M, hgradu, beta, hconju); 
      #ifdef DEBUG 
        fprintf(stderr,COLOR_GREEN"  %d: R %.3"RF"e %+.3"RF"e A %.3"RF"e B %.3"RF"e P %.3"RF"e q %.3"RF"e\n"COLOR_RESET,iter,rsqrt(resx),rsqrt(resu),alpha,beta,positive,q);
      #endif   
    };
    #ifdef DEBUG 
      fprintf(stderr,COLOR_BLUE"  %d: Hessian bottom %.3"RF"e\n"COLOR_RESET,iter,low_boundary-mu);
    #endif    
    mu=mu*0.9+(mu-low_boundary+min_mu)*0.1;
    add_inplace(N,xn,x0); add_inplace(M,un,u); 
  };
  //stop:{};
  free(u); free(xn); free(un); 
  free(gradx); free(gradu); free(conjx); free(conju);
  free(hgradx); free(hgradu); free(hconjx); free(hconju);
  return status; 
}