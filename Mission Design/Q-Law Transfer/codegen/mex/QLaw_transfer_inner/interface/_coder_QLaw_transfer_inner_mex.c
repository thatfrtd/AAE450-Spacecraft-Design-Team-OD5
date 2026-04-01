/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_QLaw_transfer_inner_mex.c
 *
 * Code generation for function '_coder_QLaw_transfer_inner_mex'
 *
 */

/* Include files */
#include "_coder_QLaw_transfer_inner_mex.h"
#include "QLaw_transfer_inner_data.h"
#include "QLaw_transfer_inner_initialize.h"
#include "QLaw_transfer_inner_terminate.h"
#include "_coder_QLaw_transfer_inner_api.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void QLaw_transfer_inner_mexFunction(int32_T nlhs, mxArray *plhs[8],
                                     int32_T nrhs, const mxArray *prhs[22])
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  const mxArray *outputs[8];
  int32_T i;
  st.tls = emlrtRootTLSGlobal;
  /* Check for proper number of arguments. */
  if (nrhs != 22) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 22, 4,
                        19, "QLaw_transfer_inner");
  }
  if (nlhs > 8) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 19,
                        "QLaw_transfer_inner");
  }
  /* Call the function. */
  QLaw_transfer_inner_api(prhs, nlhs, outputs);
  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    i = 1;
  } else {
    i = nlhs;
  }
  emlrtReturnArrays(i, &plhs[0], &outputs[0]);
}

void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs,
                 const mxArray *prhs[])
{
  mexAtExit(&QLaw_transfer_inner_atexit);
  QLaw_transfer_inner_initialize();
  QLaw_transfer_inner_mexFunction(nlhs, plhs, nrhs, prhs);
  QLaw_transfer_inner_terminate();
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLSR2022a(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1,
                           NULL, "windows-1252", true);
  return emlrtRootTLSGlobal;
}

/* End of code generation (_coder_QLaw_transfer_inner_mex.c) */
