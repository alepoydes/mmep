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
// Gradient of grad f=Hx+B is given by operator F:x->F(x) with signature
//   void F(const real* x, real* y, real* E)
// If E is not null, the value of f is stored in E.
// The minimum is obtained as the limit of the sequence
// x(k+1)=x(k)-alpha(k)*grad_f(x(k)).
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
// The method can be used for optimization with constrains, in the case
// routine 'P' must project approximation 'a' to manifold of solutions
// satisfying constrains, and 'T' must project vector field 't'
// to tangent space at the point 'a':
//   real P(int n, real* a)
//   void T(int n, const real* a, real* t)
// 'P' returns squared constrains residual divided by number of constrains.
// Return 0 on success (residual is less then epsilon)
//        1 if maximum number of iteration max_iter is reached
//        2 if machine precision is reached for 'f' evaluation
int steepest_descend(
	int n, real* a, 
  void (*TF)(const real* x, real* y, realp* E),
	int mode, real mode_param, real epsilon, int max_iter,
	void (*display)(int iter, real* a, real* grad_f, realp f, real res, real constres, real alpha, realp last_f, real last_grad),
  real (*P)(real* a)
	)
{
  real alpha;
  assert(TF);
  realp last_f=-INFINITY, f=NAN;
  real constres=NAN, res=NAN, last_res=0;

  fprintf(stderr, "Mode: "COLOR_BLUE"%d"COLOR_RESET"\n", mode);

  real init_alpha() {
    switch(mode) {
  	  case SDM_CONSTANT: return mode_param;
      case SDM_INERTIAL: return (res>1?1/res:1)*0.1; 
      case SDM_PROGR: return (res>1?1/res:1)*0.1;
    };
    fprintf(stderr,"Unknown optimization mode: " COLOR_RED "%d" COLOR_RESET "\n", mode);
    exit(1);
  };

  real update_alpha(real alpha) {
    switch(mode) {
      case SDM_CONSTANT: return alpha;
      case SDM_INERTIAL: return alpha*(1+mode_param); 
      case SDM_PROGR: return alpha*(1+mode_param); 
    };
    assert(0);
  };

  real decrease_alpha(real alpha) {
    switch(mode) {
      case SDM_CONSTANT: return alpha;
      case SDM_INERTIAL: return alpha/(1+mode_param);
      case SDM_PROGR: return alpha/2;
    };
    assert(0);
  };

  real* grad=(real*)malloc(sizeof(real)*n); assert(grad);
  real* bufa=(real*)malloc(sizeof(real)*n); assert(bufa);
  real* bufgrad=(real*)malloc(sizeof(real)*n); assert(bufgrad);
  real* toreturn=a; real* torelease=bufa; 
  int status=1;

  if(P) constres=P(a); else constres=0;
  TF(a,grad,&f);
  real res2=normsq(n,grad); res=rsqrt(res2/n); assert(!isnan(res));  
  last_res=res; last_f=f;
  alpha=init_alpha();

  int iter=1;
  while(iter<max_iter && status!=2) {
    if(res+constres<epsilon) { // If solution is found
      status=0; break; 
    };
    if(stop_signal>0) {
      fprintf(stderr, COLOR_YELLOW "Optimization aborted\n" COLOR_RESET);
      stop_signal=0;
      break;
    };      
    //if(display) display(iter,a,grad,f,res,constres,alpha,last_f,last_res); 
    //alpha=(res<=last_res)?update_alpha(alpha):decrease_alpha(alpha);
    //alpha=update_alpha(rabs(alpha));
    last_res=res;
    last_f=f;
    // update minimum
    mult_sub_ext(n,alpha,grad,a,bufa);
    if(P) constres=P(bufa);
    while(iter<max_iter) {
      TF(bufa,bufgrad,&f); iter++;
      //print_vector(n, bufa);
      //print_vector(n, bufgrad);
      res2=normsq(n,bufgrad); res=rsqrt(res2/n); assert(!isnan(res));
      if(last_f==f) { 
        fprintf(stderr, "Maximum precision is reached\n");
        status=2; break; 
      }; 
      if(display) display(iter,bufa,bufgrad,f,res,constres,alpha,last_f,last_res);
      if(stop_signal>0) break;  
      if(f<last_f || mode==SDM_CONSTANT) {
        alpha=update_alpha(rabs(alpha));  
        break;
      };
      /*if(alpha<0) { 
        fprintf(stderr, "Minimum with non-zero gradient\n");
        status=2; break; 
      };*/
      if(rabs(alpha)<1e-5) {
        fprintf(stderr, "Gradient error is too large\n");
        status=2; break;
      };
      if(mode==SDM_INERTIAL) {
        alpha=-decrease_alpha(alpha);
        //break;
      } else alpha=decrease_alpha(alpha);
      // step was too large
      mult_sub_ext(n,alpha,grad,a,bufa);
      if(P) constres=P(bufa);
    };
    real* tmp=a; a=bufa; bufa=tmp;
    tmp=bufgrad; bufgrad=grad; grad=tmp;
  };
  if(display) display(-iter,a,grad,f,res,constres,alpha,last_f,last_res);
  if(a!=toreturn) copy_vector(n, a, toreturn);
  free(bufgrad); free(torelease); free(grad); 
  return status;
};

// I is an integrator
int flow_descend(
  int n, real* a, 
  void (*TF)(const real* x, real* y, realp* E),
  int mode, real mode_param, real epsilon, int max_iter,
  void (*display)(int iter, real* a, real* grad_f, realp f, real res, real constres, real alpha),
  real (*P)(real* a),
  void (*I)(int N, void (*F)(const real* x, real* g, realp* E), real T, real* X, realp* E, int* iter)
  )
{
  int iter=0; real alpha; real res=NAN;
  assert(TF); 
  switch(mode) {
    case SDM_CONSTANT: alpha=mode_param; break;
    case SDM_INERTIAL: alpha=NAN; break;
    case SDM_PROGR: alpha=mode_param; break;
    default: assert(1);
  };
  int status=1;
  realp last_f=-INFINITY; realp f=NAN;
  P(a);
  while(iter<max_iter) {
    int iter_plus; I(n, TF, -alpha, a, &f, &iter_plus); iter+=iter_plus;
    assert(!isnan(f));
    P(a);
    res=fabs((last_f-f)/alpha);
    switch(mode) {
      case SDM_CONSTANT: break;
      case SDM_INERTIAL: 
        if(f<=last_f) alpha+=mode_param*res; else alpha=mode_param*res;
        break;
      case SDM_PROGR: 
        if(f<=last_f) alpha+=mode_param; else alpha=mode_param;
        break;
      default: assert(0);
    };
    if(res<epsilon) { // If solution is found
      status=0; break; 
    };
    if(stop_signal>0) {
      fprintf(stderr, COLOR_YELLOW "Optimization aborted\n" COLOR_RESET);
      stop_signal=0;
      break;
    };
    if(display) display(iter,a,NULL,f,res,NAN,alpha); 
    //if(last_f==f) { status=2; break; }; // Iterations stop changing
    last_f=f;
  };
  if(display) display(-iter,a,NULL,f,res,NAN,alpha);
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
int lagrange_conjugate_quad(
  int N, int M, real* x0, 
  void (*Q)(const real* x, real* y),
  void (*L)(real* y),
  int mode, real mode_param, real epsilon, int max_iter,
  void (*display)(int iter, real* a, real* grad_f, realp f, real res, real alpha),
  void (*C)(const real* x, real* r),
  void (*D)(const real* x, const real* u, real* r),
  void (*P)(const real* x, const real* y, real* r),
  real initial_mu
  )
{
  #define MU_INC mu-positive+mode_param
  real improvement=0.01;
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
  real* kktx=(real*)malloc(sizeof(real)*N); assert(kktx);
  real* kktu=(real*)malloc(sizeof(real)*M); assert(kktu);
  real* hkktx=(real*)malloc(sizeof(real)*N); assert(hkktx);
  // initialize state
  int iter=0;
  int status=1;  
  real mu=initial_mu;
  // Exact solution for orthonormal P_j x0 and exact minimizer x0.
  Q(x0,gradx); iter++; L(gradx); P(gradx,x0,u); negate_inplace(M,u);
  
  //zero_vector(M, u);
  real q=INFINITY;
  while(iter<max_iter) {
    // Calculate gradient at the point
    restart: Q(x0,gradx); iter++;
    realp f=-dot(N,x0,gradx)/2;
    L(gradx);
    f+=dot(N,x0,gradx);
    D(x0,u,gradx);
    mult_add(N,mu,x0,gradx);
    C(x0,gradu); 
    // Calculate norm of residual
    real resx=normsq(N,gradx), resu=normsq(M,gradu);
    real res=rsqrt(resx/N+resu/M); real top_res=res*improvement;
    real last_res=res;
    // Diplay progress
    #ifdef DEBUG
      fprintf(stderr,COLOR_RED "%d: res=%" RF "g %+" RF "g mu=%" RF "g\n" COLOR_RESET,iter,RT(rsqrt(resx/N)),RT(rsqrt(resu/M)),RT(mu));
      if(N+M<=8) {
        fprintf(stderr, COLOR_YELLOW);
        for(int j=0;j<N;j++) fprintf(stderr,"%" RF "g ",RT(x0[j])); fprintf(stderr,":");
        for(int j=0;j<M;j++) fprintf(stderr," %" RF "g",RT(u[j])); 
        fprintf(stderr,COLOR_RESET "\n");
      };
    #endif
    // Residula is negative gradient
    negate_inplace(N,gradx); negate_inplace(M,gradu);
    // Reporting result
    if(res<epsilon) { 
      #ifdef DEBUG
        fprintf(stderr,COLOR_RED "%d: Converged %.3" RF "g < %.3" RF "g\n" COLOR_RESET,iter,RT(res),RT(epsilon));
      #endif    
      if(display) display(-iter,x0,gradx,f,res,mu);
      status=0; break; 
    }; // If solution is found
    if(display) display(iter,x0,gradx,f,res,mu);
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

    // Chech suffucuent conditions for minimum
    // D:x,u,r->r+sum_l u_l P_j x
    // P:x,y->(<x|P_j y>)_l
    copy_vector(N,gradx,kktx);
    P(gradx,x0,kktu); P(x0,x0,hgradu); 
    negate_div(M,hgradu,kktu); D(x0,kktu,kktx);
    // kktx contains projection of gradx to tangent space
    Q(kktx,hkktx); D(kktx,u,hkktx); 
    //mult_add(N,mu,kktx,hkktx);
    real positive=mu+dot(N,kktx,hkktx)/normsq(N,kktx); // Must be positive if x0 is minimum
    //fprintf(stderr,"orth %" RF "e\n",dot(N,kktx,x0));
    // If the quadratic form is not positive, update shift 'mu'
    if(positive<min_mu) { 
      #ifdef DEBUG
        fprintf(stderr,"  %d: " COLOR_BLUE COLOR_BOLD "Hessian is small:" COLOR_RESET " %.3" RF "e < %.3" RF "e\n",iter,RT(positive),RT(min_mu));
      #endif
      mu=MU_INC;  
      continue; 
    };
    real low_boundary=positive; // Estimate of bottom of Hessian spectrum
    D(x0,gradu,hgradx); P(x0,gradx,hgradu);
    //vector_copy(M,gradu,hgradu);
    // choosing initial conjugate direction
    copy_vector(N, gradx, conjx); copy_vector(M, gradu, conju);
    copy_vector(N, hgradx, hconjx); copy_vector(M, hgradu, hconju);
    // Initialization of conjugate residual
    q=dot(N,gradx,hgradx)+dot(M,gradu,hgradu); 
    while(iter<max_iter) {
      if(rabs(q)<dot_precision) { 
        #ifdef DEBUG        
          fprintf(stderr,COLOR_BLUE COLOR_ITALIC "  %d: Small q: %.3" RF "e < %.3" RF "e\n" COLOR_RESET,iter,RT(q),RT(dot_precision));
        #endif        
        //mu+=min_mu; goto restart; 
        add_inplace(N,xn,x0); 
        if(display) display(-iter,x0,gradx,f,res,mu);
        goto stop;
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
        if(N+M<=8) {
          fprintf(stderr, COLOR_YELLOW "  X:" COLOR_RESET); for(int j=0;j<N;j++) fprintf(stderr," %" RF "g",RT(xn[j])); fprintf(stderr,":"); for(int j=0;j<M;j++) fprintf(stderr," %" RF "g",RT(un[j])); fprintf(stderr,"\n");
          fprintf(stderr, COLOR_YELLOW "  G:" COLOR_RESET); for(int j=0;j<N;j++) fprintf(stderr," %" RF "g",RT(gradx[j])); fprintf(stderr,":"); for(int j=0;j<M;j++) fprintf(stderr," %" RF "g",RT(gradu[j])); fprintf(stderr,"\n");
          fprintf(stderr, COLOR_YELLOW "  C:" COLOR_RESET); for(int j=0;j<N;j++) fprintf(stderr," %" RF "g",RT(conjx[j])); fprintf(stderr,":"); for(int j=0;j<M;j++) fprintf(stderr," %" RF "g",RT(conju[j])); fprintf(stderr,"\n");
        };
      #endif      
      resx=normsq(N,gradx); resu=normsq(M,gradu); res=rsqrt(resx/N+resu/N);
      assert(!isnan(res));
      //if(res<epsilon) { // Auxilliary problem is solved
      if(/*res<epsilon ||*/ res<top_res) { // Auxilliary problem is solved
        #ifdef DEBUG 
          fprintf(stderr,COLOR_GREEN "  %d: Residual is small: R %.3" RF "e %+.3" RF "e \n" COLOR_RESET,iter,RT(rsqrt(resx/N)),RT(rsqrt(resu/M)));
        #endif  
        break; 
      };
      if(res>100*last_res) { // Unstability detected ???
        #ifdef DEBUG 
          fprintf(stderr,COLOR_GREEN "  %d: Residual increases: R %.3" RF "e %+.3" RF "e \n" COLOR_RESET,iter,RT(rsqrt(resx/N)),RT(rsqrt(resu/M)));
        #endif  
        break; 
      };
      last_res=res;
      // Aplying quadratic form to residual
      iter++;
      Q(gradx,hgradx); 
      D(gradx,u,hgradx);
      mult_add(N,mu,gradx,hgradx);
      
      // Chech suffucuent conditions for minimum
      // D:x,u,r->r+sum_l u_l P_j x
      // P:x,y->(<x|P_j y>)_l
      copy_vector(N,gradx,kktx);
      P(gradx,x0,kktu); P(x0,x0,hgradu); 
      negate_div(M,hgradu,kktu); D(x0,kktu,kktx);
      // kktx contains projection of gradx to tangent space
      Q(kktx,hkktx); D(kktx,u,hkktx); D(kktx,un,hkktx); 
      real positive=mu+dot(N,kktx,hkktx)/normsq(N,kktx); // Must be positive if x0 is minimum
      if(positive<min_mu) { 
        #ifdef DEBUG 
          fprintf(stderr,"  %d: " COLOR_BLUE COLOR_BOLD "Hessian is small:" COLOR_RESET " %.3" RF "e < %.3" RF "e\n",iter,RT(positive),RT(min_mu));
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
        fprintf(stderr,COLOR_GREEN "  %d: R %.3" RF "e %+.3" RF "e A %.3" RF "e B %.3" RF "e P %.3" RF "e q %.3" RF "e\n" COLOR_RESET,iter,RT(rsqrt(resx/N)),RT(rsqrt(resu/M)),RT(alpha),RT(beta),RT(positive),RT(q));
      #endif   
    };
    #ifdef DEBUG 
      fprintf(stderr,COLOR_BLUE "  %d: Hessian bottom %.3" RF "e\n" COLOR_RESET,iter,RT(low_boundary-mu));
    #endif    
    real delta_mu=(-low_boundary+min_mu)*0.3;
    //delta_mu=0;
    if(mu+delta_mu<0) delta_mu=-mu;
    mu+=delta_mu;
    add_constant_inplace(M,-delta_mu,u);
    add_inplace(N,xn,x0); add_inplace(M,un,u); 
  };
  stop: {}; 
  free(u); free(xn); free(un); 
  free(gradx); free(gradu); free(conjx); free(conju);
  free(hgradx); free(hgradu); free(hconjx); free(hconju);
  free(kktx); free(kktu); free(hkktx);
  return status; 
}


int lagrange_conjugate(
  int N, int M, real* x0, 
  void (*Q)(const real* x, real* y),
  void (*L)(real* y),
  int mode, real mode_param, real epsilon, int max_iter,
  void (*display)(int iter, real* a, real* grad_f, realp f, real res, real alpha),
  void (*C)(const real* x, real* r),
  void (*D)(const real* x, const real* u, real* r),
  void (*P)(const real* x, const real* y, real* r),
  real initial_mu
  )
{
  #define MU_INC mu-positive+mode_param
  real improvement=0.001;
  real min_mu=1e-4; // Minimum required bpttpm of spectrum of Q+mu
  real dot_precision=EPSILON*10; // Precision of inner product
  // allocate mempry
  //real* u=(real*)malloc(sizeof(real)*M); assert(u);
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
  real* kktx=(real*)malloc(sizeof(real)*N); assert(kktx);
  real* kktu=(real*)malloc(sizeof(real)*M); assert(kktu);
  real* hkktx=(real*)malloc(sizeof(real)*N); assert(hkktx);
  // initialize state
  int iter=0;
  int status=1;  
  real mu=initial_mu;
  real low_boundary=0;
  //zero_vector(M, u);
  real q=INFINITY;
  while(iter<max_iter) {
    // Calculate gradient at the point
    restart: 
    Q(x0,gradx); iter++;
    real f=-dot(N,x0,gradx)/2;
    L(gradx);
    f+=dot(N,x0,gradx);
    //D(x0,u,gradx);
    mult_add(N,mu,x0,gradx);
    C(x0,gradu); 

    // initial guess
    zero_vector(N, xn); 
    P(gradx,x0,un); 
    negate_inplace(M,un); 
    D(x0,un,gradx);
    negate_inplace(N,gradx); negate_inplace(M,gradu); 

    // Calculate norm of residual
    real resx=normsq(N,gradx); real resu=normsq(M,gradu);
    real res=rsqrt(resx/N+resu/M); real top_res=res*improvement;
    real last_res=res;
    // Diplay progress
    #ifdef DEBUG
      fprintf(stderr,COLOR_RED "%d: res=%" RF "g %+" RF "g mu=%" RF "g\n" COLOR_RESET,iter,RT(rsqrt(resx/N)),RT(rsqrt(resu/M)),RT(mu));
      if(N+M<=8) {
        fprintf(stderr, COLOR_YELLOW);
        for(int j=0;j<N;j++) fprintf(stderr,"%" RF "g ",RT(x0[j])); 
        //fprintf(stderr,":"); for(int j=0;j<M;j++) fprintf(stderr," %" RF "g",u[j]); 
        fprintf(stderr,COLOR_RESET "\n");
      };
    #endif
    // Reporting result
    if(res<epsilon) { 
      #ifdef DEBUG
        fprintf(stderr,COLOR_RED "%d: Converged %.3" RF "g < %.3" RF "g\n" COLOR_RESET,iter,RT(res),RT(epsilon));
      #endif    
      if(display) display(-iter,x0,gradx,f,res,mu);
      status=0; break; 
    }; // If solution is found
    if(display) display(iter,x0,gradx,f,res,mu);
    // Solving auxilliary problem with respect to xn,un
    // [...]*[xn;un]=[gradx;gradu]

    // computing quadratic part
    //h=[A+2*lambda+shift,2*x;2*x',0];

    // Chech suffucuent conditions for minimum
    // D:x,u,r->r+sum_l u_l P_j x
    // P:x,y->(<x|P_j y>)_l
    copy_vector(N,gradx,kktx);
    P(gradx,x0,kktu); P(x0,x0,hgradu); 
    negate_div(M,hgradu,kktu); D(x0,kktu,kktx);
    // kktx contains projection of gradx to tangent space
    Q(kktx,hkktx); 
    real positive=mu+dot(N,kktx,hkktx)/normsq(N,kktx); // Must be positive if x0 is minimum
    //fprintf(stderr,"orth %" RF "e\n",dot(N,kktx,x0));
    // If the quadratic form is not positive, update shift 'mu'
    if(positive<min_mu) { 
      #ifdef DEBUG
        fprintf(stderr,"  %d: " COLOR_BLUE COLOR_BOLD "Hessian is small:" COLOR_RESET " %.3" RF "e < %.3" RF "e\n",iter,RT(positive),RT(min_mu));
      #endif
      mu=MU_INC;  
      continue; 
    };
    if(low_boundary>positive-mu) low_boundary=positive-mu; // Estimate of bottom of Hessian spectrum

    iter++;
    Q(gradx,hgradx); 
    mult_add(N,mu,gradx,hgradx);
    D(x0,gradu,hgradx); P(x0,gradx,hgradu);
    //vector_copy(M,gradu,hgradu);
    // choosing initial conjugate direction
    copy_vector(N, gradx, conjx); copy_vector(M, gradu, conju);
    copy_vector(N, hgradx, hconjx); copy_vector(M, hgradu, hconju);
    // Initialization of conjugate residual
    q=dot(N,gradx,hgradx)+dot(M,gradu,hgradu); 
    while(iter<max_iter) {
      if(rabs(q)<dot_precision) { 
        #ifdef DEBUG        
          fprintf(stderr,COLOR_BLUE COLOR_ITALIC "  %d: Small q: %.3" RF "e < %.3" RF "e\n" COLOR_RESET,iter,RT(q),RT(dot_precision));
        #endif        
        //mu+=min_mu; goto restart; 
        add_inplace(N,xn,x0); 
        if(display) display(-iter,x0,gradx,f,res,mu);
        goto stop;
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
        if(N+M<=8) {
          fprintf(stderr, COLOR_YELLOW "  X:" COLOR_RESET); for(int j=0;j<N;j++) fprintf(stderr," %" RF "g",RT(xn[j])); fprintf(stderr,":"); for(int j=0;j<M;j++) fprintf(stderr," %" RF "g",RT(un[j])); fprintf(stderr,"\n");
          fprintf(stderr, COLOR_YELLOW "  G:" COLOR_RESET); for(int j=0;j<N;j++) fprintf(stderr," %" RF "g",RT(gradx[j])); fprintf(stderr,":"); for(int j=0;j<M;j++) fprintf(stderr," %" RF "g",RT(gradu[j])); fprintf(stderr,"\n");
          fprintf(stderr, COLOR_YELLOW "  C:" COLOR_RESET); for(int j=0;j<N;j++) fprintf(stderr," %" RF "g",RT(conjx[j])); fprintf(stderr,":"); for(int j=0;j<M;j++) fprintf(stderr," %" RF "g",RT(conju[j])); fprintf(stderr,"\n");
        };
      #endif      
      resx=normsq(N,gradx); resu=normsq(M,gradu); res=rsqrt(resx/N+resu/N);
      assert(!isnan(res));
      //if(res<epsilon) { // Auxilliary problem is solved
      if(res<top_res) { // Auxilliary problem is solved
        #ifdef DEBUG 
          fprintf(stderr,COLOR_GREEN "  %d: Residual is small: R %.3" RF "e %+.3" RF "e \n" COLOR_RESET,iter,RT(rsqrt(resx/N)),RT(rsqrt(resu/M)));
        #endif  
        break; 
      };
      if(res>100*last_res) { // Unstability detected ???
        #ifdef DEBUG 
          fprintf(stderr,COLOR_GREEN "  %d: Residual increases: R %.3" RF "e %+.3" RF "e \n" COLOR_RESET,iter,RT(rsqrt(resx/N)),RT(rsqrt(resu/M)));
        #endif  
        break; 
      };
      last_res=res;
      // Aplying quadratic form to residual
      iter++;
      Q(gradx,hgradx); 
      mult_add(N,mu,gradx,hgradx);
      
      // Chech suffucuent conditions for minimum
      // D:x,u,r->r+sum_l u_l P_j x
      // P:x,y->(<x|P_j y>)_l
      copy_vector(N,gradx,kktx);
      P(gradx,x0,kktu); P(x0,x0,hgradu); 
      negate_div(M,hgradu,kktu); D(x0,kktu,kktx);
      // kktx contains projection of gradx to tangent space
      Q(kktx,hkktx); 
      real positive=mu+dot(N,kktx,hkktx)/normsq(N,kktx); // Must be positive if x0 is minimum
      if(positive<min_mu) { 
        #ifdef DEBUG 
          fprintf(stderr,"  %d: " COLOR_BLUE COLOR_BOLD "Hessian is small:" COLOR_RESET " %.3" RF "e < %.3" RF "e\n",iter,RT(positive),RT(min_mu));
        #endif  
        mu=MU_INC; 
        goto restart; 
      };
      if(positive-mu<low_boundary) low_boundary=positive-mu;
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
        fprintf(stderr,COLOR_GREEN "  %d: R %.3" RF "e %+.3" RF "e A %.3" RF "e B %.3" RF "e P %.3" RF "e q %.3" RF "e\n" COLOR_RESET,iter,RT(rsqrt(resx/N)),RT(rsqrt(resu/M)),RT(alpha),RT(beta),RT(positive),RT(q));
      #endif   
    };
    #ifdef DEBUG 
      fprintf(stderr,COLOR_BLUE "  %d: Hessian bottom %.3" RF "e\n" COLOR_RESET,iter,RT(low_boundary));
    #endif    
    real delta_mu=(-mu-2*low_boundary+min_mu)*0.1;
    delta_mu=0;
    if(mu+delta_mu<0) delta_mu=-mu;
    mu+=delta_mu;
    add_inplace(N,xn,x0); 
  };
  stop: {}; 
  //free(u); 
  free(xn); free(un); 
  free(gradx); free(gradu); free(conjx); free(conju);
  free(hgradx); free(hgradu); free(hconjx); free(hconju);
  free(kktx); free(kktu); free(hkktx);
  return status; 
}
