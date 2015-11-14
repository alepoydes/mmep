#define SDM_CONSTANT 0
#define SDM_INERTIAL 1
#define SDM_PROGR 2
#define SDM_CAUCHY 3
// Compute minimum doing gradient descend.
// Minimizing function f(x) is of special form:
//   f(x)=<x|Hx>/2+<x|B>.
// Hessian part H is given by operator Q:x->H(x) with signature
//   void Q(int n, const float* x, float* y)
// Linear part is given as an operator L:y->y+B
//   void L(int n, float* y)
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
// Return 0 on success (residual is less then epsilon)
//        1 if maximum number of iteration max_iter is reached
int steepest_descend(
	int n, float* a, 
	void (*Q)(const float* x, float* y),
	void (*L)(float* y),
	int mode, float mode_param, float epsilon, int max_iter,
	void (*display)(int iter, float* a, float* grad_f, float f, float res, float alpha),
    void (*P)(float* a),
    void (*T)(const float* a, float* t)
);

int lagrange_conjugate(
  int N, int M, float* x0, 
  void (*Q)(const float* x, float* y),
  void (*L)(float* y),
  int mode, float mode_param, float epsilon, int max_iter,
  void (*display)(int iter, float* a, float* grad_f, float f, float res, float alpha),
  void (*C)(const float* x, float* r),
  void (*D)(const float* x, const float* u, float* r),
  void (*P)(const float* x, const float* y, float* r)
 );