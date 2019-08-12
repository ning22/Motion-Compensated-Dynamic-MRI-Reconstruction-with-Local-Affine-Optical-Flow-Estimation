/* Fast 2D inverse wavelet transform using Complex Wavelets
 * Wriiten by: Salman Asif, Georgia Tech
 *
 * ------------------------------
 * Copyright (c) 2011 Salman Asif
 * ------------------------------
 */

#include <stdlib.h>
#include "mex.h"

/* these are in fwt_level.c */
void inv_wavelet_one_level_gen(double*, double*, int, double*, double*, int, int,int);
int checkPowerTwo(int,int);


void inv_wavelet_multi_level2D_CWT(double *x, double *w, int nR, int nC, double *g0_r, double *g1_r, double *g0_c, double *g1_c, int l, int J, int symr, int symc, int rosym_flip, int cosym_flip)
{
  int j, nRj, nCj, n, m1, m2;
  double *tmp, *tmpo;
  
  n = (nR > nC) ? nR: nC;
  /* Arrays for row/col processing */
  tmp = (double*) mxCalloc(n, sizeof(double));
  tmpo = (double*) mxCalloc(n, sizeof(double));

  /* find the right value of nj */
  nRj = nR;
  nCj = nC;
  for (j=0; j < J; j++){
    nRj >>= 1;
    nCj >>= 1;
  }
  
  for (j=J-1; j >= 0; j--) {
    nRj <<= 1;
    nCj <<= 1;
    
    /* columns */
    for (m2=0; m2 < nCj; m2++) {
      if (j==J-1)
        inv_wavelet_one_level_gen(&x[m2*nR], &w[m2*nR], nRj, g0_c, g1_c, l, symc,cosym_flip);
      else {
        for (m1=0; m1 < nRj/2; m1++) {
          tmp[m1] = (m2 < nCj/2) ? x[m2*nR+m1] : w[m2*nR+m1];
          tmp[nRj/2+m1] = w[m2*nR+nRj/2+m1];
        }
        inv_wavelet_one_level_gen(&x[m2*nR], tmp, nRj, g0_c, g1_c, l, symc,cosym_flip);
      }
    }
    /* rows */
    for (m1=0; m1 < nRj; m1++) {
      for (m2=0; m2 < nCj; m2++) {
        tmp[m2] = x[nR*m2+m1];
      }
      inv_wavelet_one_level_gen(tmpo, tmp, nCj, g0_r, g1_r, l, symr,rosym_flip);
      for (m2=0; m2 < nCj; m2++)
        x[nR*m2+m1] = tmpo[m2];
    }
      
  }
   
   mxFree(tmp);
   mxFree(tmpo);
}
  

/* The gateway routine. */
/* x = ifwt2_CWT(w, g0_r, g1_r, g0_c, g1_c, J, symr, symc, rosym_flip, cosym_flip); */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *x, *w, *g0_r, *g1_r, *g0_c, *g1_c;
  int nR, nC, l, l1, symr, symc, rosym_flip, cosym_flip, J, JmaxR, JmaxC;
  
  /* Check for the proper number of arguments. */
  if (nrhs != 10) {
    mexErrMsgTxt("Exactly ten inputs required");
  }
  if (nlhs > 1) {
    mexErrMsgTxt("Too many output arguments");
  }
  
  nR = mxGetM(prhs[0]);
  nC = mxGetN(prhs[0]);
  /*
   * if (n != n1) {
    mexErrMsgTxt("Input w must be a square image.");
  }
   */
  if (mxIsComplex(prhs[0])) {
    mexErrMsgTxt("Input w must be real");    
  }
  J = mxGetScalar(prhs[5]);
  JmaxR = checkPowerTwo(nR, J);
  JmaxC = checkPowerTwo(nC, J);
  if ((J < 0) || (J > JmaxR) || (J > JmaxC) ) {
        mexErrMsgTxt("Input J must be an integer between 0 and log2(n), and dyadic for ROW and COL---use smaller J.");
  }
  l = (mxGetM(prhs[1]) > mxGetN(prhs[1])) ? mxGetM(prhs[1]) : mxGetN(prhs[1]);
  l1 = (mxGetM(prhs[2]) > mxGetN(prhs[2])) ? mxGetM(prhs[2]) : mxGetN(prhs[2]);
  if (l != l1) {
    mexErrMsgTxt("Filters must be the same length");
  }
  
  symr = mxGetScalar(prhs[6]);
  symc = mxGetScalar(prhs[7]);
  rosym_flip = mxGetScalar(prhs[8]); /* tells if we need to flip odd symmetry of scaling and wavelet coeffs */
  cosym_flip = mxGetScalar(prhs[9]);
  
  /* mexPrintf("rosym_flip is %d, cosym_flip is %d, \n", rosym_flip, cosym_flip); */
  if ((symr > 2) || (symc > 2)) {
    mexErrMsgTxt("Symmetry flag must be 0, 1, or 2");
  }
  
  /* Create matrix for the return argument. */
  plhs[0] = mxCreateDoubleMatrix(nR, nC, mxREAL);
  
  /* Assign pointers to each input and output. */
  w = mxGetPr(prhs[0]);
  g0_r = mxGetPr(prhs[1]);
  g1_r = mxGetPr(prhs[2]);
  g0_c = mxGetPr(prhs[3]);
  g1_c = mxGetPr(prhs[4]);
  
  x = mxGetPr(plhs[0]);
  
  /* Call the C subroutine. */
  inv_wavelet_multi_level2D_CWT(x, w, nR, nC, g0_r, g1_r, g0_c, g1_c, l, J, symr, symc, rosym_flip, cosym_flip);
  return;
}