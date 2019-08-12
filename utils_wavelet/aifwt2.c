/* Fast 2D wavelet transform 
 * Wriiten by: Justin Romberg, Georgia Tech
 */

#include <stdlib.h>
/* #include "mex.h" */
#include "fwt_level.c"

/* these are in aifwt_level.c */
void adjmod_inverse_wavelet_one_level(double*, double*, int, double*, double*, int, int);
/* int checkPowerTwo(int); */

void adj_inverse_wavelet_multi_level2D(double *u, double *v, int n, double *g0, double *g1, int l, int J, int sym)
{
  int j, nj, m1, m2;
  double *tmp, *tmpo;
  
  
  /* tmp = (double*) mxCalloc(n, sizeof(double)); */
  /* tmpo = (double*) mxCalloc(n, sizeof(double)); */
  tmp = (double*) calloc(n, sizeof(double));
  tmpo = (double*) calloc(n, sizeof(double));
    
  nj = n;
  for (j=0; j < J; j++) {
    /* columns */
    for (m1=0; m1 < nj; m1++) {  
      if (j==0)
        /* adj_inverse_wavelet_one_level(&u[m1*n], &v[m1*n], nj, g0, g1, l, sym); */
		adjmod_inverse_wavelet_one_level(&u[m1*n], &v[m1*n],nj,g0,g1,l,sym);
      else {
       for (m2=0; m2 < nj; m2++)
         tmp[m2] = u[m1*n+m2];
       /* adj_inverse_wavelet_one_level(&u[m1*n], tmp, nj, g0, g1, l, sym); */
	   adjmod_inverse_wavelet_one_level(&u[m1*n], tmp, nj, g0, g1, l, sym);

      }
    }
    /* rows */
    for (m1=0; m1 < nj; m1++) {
      for (m2=0; m2 < nj; m2++)
        tmp[m2] = u[n*m2+m1];
      /* adj_inverse_wavelet_one_level(tmpo, tmp, nj, g0, g1, l, sym); */
	  adjmod_inverse_wavelet_one_level(tmpo, tmp, nj, g0, g1, l, sym);
      
      for (m2=0; m2 < nj; m2++)
        u[n*m2+m1] = tmpo[m2];
    }
    
    nj >>= 1;
  }
   
  /*mxFree(tmp); 
  mxFree(tmpo); */
  free(tmp);
  free(tmpo);
}



/* The gateway routine. *
/* w = fwt2(x, h0, h1, J, sym); *
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *v, *u, *g0, *g1;
  int n, n1, l, l1, sym, J, Jmax;
  
  /* Check for the proper number of arguments. *
  if (nrhs != 5) {
    mexErrMsgTxt("Exactly five inputs required");
  }
  if (nlhs > 1) {
    mexErrMsgTxt("Too many output arguments");
  }
  
  n = mxGetM(prhs[0]);
  n1 = mxGetN(prhs[0]);
  if (n != n1) {
    mexErrMsgTxt("Input x must be a square image.");
  }
  if (mxIsComplex(prhs[0])) {
    mexErrMsgTxt("Input x must be real");    
  }
  Jmax = checkPowerTwo(n);
  l = (mxGetM(prhs[1]) > mxGetN(prhs[1])) ? mxGetM(prhs[1]) : mxGetN(prhs[1]);
  l1 = (mxGetM(prhs[2]) > mxGetN(prhs[2])) ? mxGetM(prhs[2]) : mxGetN(prhs[2]);
  if (l != l1) {
    mexErrMsgTxt("Filters must be the same length");
  }
  J = mxGetScalar(prhs[3]);
  if ((J < 0) || (J > Jmax)) {
    mxErrMsgTxt("Input J must be an integer between 0 and log2(n)");
  }
  sym = mxGetScalar(prhs[4]);
  if (sym > 2) {
    mexErrMsgTxt("Symmetry flag must be 0, 1, or 2");
  }
  
  /* Create matrix for the return argument. *
  plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL);
  
  /* Assign pointers to each input and output. *
  v = mxGetPr(prhs[0]);
  g0 = mxGetPr(prhs[1]);
  g1 = mxGetPr(prhs[2]);
  u = mxGetPr(plhs[0]);
  
  /* Call the C subroutine. *
  adj_inverse_wavelet_multi_level2D(u, v, n, g0, g1, l, J, sym);
  return;

}*/