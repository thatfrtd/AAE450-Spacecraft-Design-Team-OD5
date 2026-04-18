/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * QLaw_transfer_inner_initialize.c
 *
 * Code generation for function 'QLaw_transfer_inner_initialize'
 *
 */

/* Include files */
#include "QLaw_transfer_inner_initialize.h"
#include "QLaw_transfer_inner_data.h"
#include "_coder_QLaw_transfer_inner_mex.h"
#include "rt_nonfinite.h"

/* Function Declarations */
static void QLaw_transfer_inner_once(void);

/* Function Definitions */
static void QLaw_transfer_inner_once(void)
{
  mex_InitInfAndNan();
}

void QLaw_transfer_inner_initialize(void)
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtBreakCheckR2012bFlagVar = emlrtGetBreakCheckFlagAddressR2022b(&st);
  emlrtClearAllocCountR2012b(&st, false, 0U, NULL);
  emlrtEnterRtStackR2012b(&st);
  if (emlrtFirstTimeR2012b(emlrtRootTLSGlobal)) {
    QLaw_transfer_inner_once();
  }
}

/* End of code generation (QLaw_transfer_inner_initialize.c) */
