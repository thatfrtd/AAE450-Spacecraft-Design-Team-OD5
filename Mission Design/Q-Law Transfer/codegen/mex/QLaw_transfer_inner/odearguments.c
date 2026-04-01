/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * odearguments.c
 *
 * Code generation for function 'odearguments'
 *
 */

/* Include files */
#include "odearguments.h"
#include "QLaw_transfer_inner.h"
#include "QLaw_transfer_inner_data.h"
#include "rt_nonfinite.h"
#include "warning.h"
#include "mwmathutil.h"

/* Variable Definitions */
static emlrtRSInfo ed_emlrtRSI = {
    40,             /* lineNo */
    "odearguments", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2025b\\toolbox\\eml\\lib\\matlab\\funfun\\private\\odeargu"
    "ments.m" /* pathName */
};

static emlrtRSInfo fd_emlrtRSI = {
    22,             /* lineNo */
    "odearguments", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2025b\\toolbox\\eml\\lib\\matlab\\funfun\\private\\odeargu"
    "ments.m" /* pathName */
};

static emlrtRTEInfo h_emlrtRTEI = {
    93,             /* lineNo */
    5,              /* colNo */
    "odearguments", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2025b\\toolbox\\eml\\lib\\matlab\\funfun\\private\\odeargu"
    "ments.m" /* pName */
};

static emlrtRTEInfo i_emlrtRTEI = {
    66,             /* lineNo */
    27,             /* colNo */
    "odearguments", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2025b\\toolbox\\eml\\lib\\matlab\\funfun\\private\\odeargu"
    "ments.m" /* pName */
};

static emlrtRTEInfo j_emlrtRTEI = {
    44,             /* lineNo */
    23,             /* colNo */
    "odearguments", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2025b\\toolbox\\eml\\lib\\matlab\\funfun\\private\\odeargu"
    "ments.m" /* pName */
};

static emlrtRTEInfo k_emlrtRTEI = {
    36,             /* lineNo */
    23,             /* colNo */
    "odearguments", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2025b\\toolbox\\eml\\lib\\matlab\\funfun\\private\\odeargu"
    "ments.m" /* pName */
};

static emlrtRTEInfo l_emlrtRTEI = {
    19,             /* lineNo */
    23,             /* colNo */
    "odearguments", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2025b\\toolbox\\eml\\lib\\matlab\\funfun\\private\\odeargu"
    "ments.m" /* pName */
};

static emlrtRTEInfo m_emlrtRTEI = {
    17,             /* lineNo */
    1,              /* colNo */
    "odearguments", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2025b\\toolbox\\eml\\lib\\matlab\\funfun\\private\\odeargu"
    "ments.m" /* pName */
};

static emlrtRTEInfo n_emlrtRTEI = {
    16,             /* lineNo */
    1,              /* colNo */
    "odearguments", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2025b\\toolbox\\eml\\lib\\matlab\\funfun\\private\\odeargu"
    "ments.m" /* pName */
};

/* Function Definitions */
real_T odearguments(const emlrtStack *sp,
                    const real_T c_odeFun_workspace_ODEFunction_[3],
                    real_T d_odeFun_workspace_ODEFunction_,
                    const real_T tspan[2], const real_T b_y0[7],
                    real_T options_AbsTol,
                    const real_T options_InitialStep_data[],
                    const real_T options_MaxStep_data[], real_T options_RelTol,
                    real_T tspanOut[2], real_T *tfinal, real_T *tdir,
                    real_T f0[7], real_T *threshold, real_T *rtol, real_T *hmax,
                    real_T htry_data[], int32_T htry_size[2], real_T *htspan)
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack st;
  real_T htryIn_data;
  real_T t0;
  int32_T k;
  boolean_T tf;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  e_st.prev = &d_st;
  e_st.tls = d_st.tls;
  t0 = tspan[0];
  htryIn_data = tspan[1] - tspan[0];
  *htspan = muDoubleScalarAbs(htryIn_data);
  *tfinal = tspan[1];
  tf = false;
  for (k = 0; k < 2; k++) {
    real_T d;
    d = tspan[k];
    tspanOut[k] = d;
    if (tf || muDoubleScalarIsNaN(d)) {
      tf = true;
    }
  }
  if (tf) {
    emlrtErrorWithMessageIdR2018a(sp, &n_emlrtRTEI,
                                  "MATLAB:odearguments:TspanNaNValues",
                                  "MATLAB:odearguments:TspanNaNValues", 0);
  }
  if (tspan[0] == tspan[1]) {
    emlrtErrorWithMessageIdR2018a(
        sp, &m_emlrtRTEI, "MATLAB:odearguments:TspanEndpointsNotDistinct",
        "MATLAB:odearguments:TspanEndpointsNotDistinct", 0);
  }
  *tdir = muDoubleScalarSign(htryIn_data);
  tf = true;
  if (!(tspan[1] > tspan[0])) {
    tf = (tspan[0] > tspan[1]);
  }
  if (!tf) {
    emlrtErrorWithMessageIdR2018a(sp, &l_emlrtRTEI,
                                  "MATLAB:odearguments:TspanNotMonotonic",
                                  "MATLAB:odearguments:TspanNotMonotonic", 0);
  }
  st.site = &fd_emlrtRSI;
  b_st.site = &wb_emlrtRSI;
  c_st.site = &gd_emlrtRSI;
  d_st.site = &hd_emlrtRSI;
  e_st.site = &wb_emlrtRSI;
  QLaw_transfer_inner_anonFcn2(&e_st, c_odeFun_workspace_ODEFunction_,
                               d_odeFun_workspace_ODEFunction_, b_y0, f0);
  *rtol = options_RelTol;
  if (!(options_RelTol > 0.0)) {
    emlrtErrorWithMessageIdR2018a(sp, &k_emlrtRTEI,
                                  "MATLAB:odearguments:RelTolNotPosScalar",
                                  "MATLAB:odearguments:RelTolNotPosScalar", 0);
  }
  if (options_RelTol < 2.2204460492503131E-14) {
    *rtol = 2.2204460492503131E-14;
    st.site = &ed_emlrtRSI;
    warning(&st);
  }
  if (!(options_AbsTol > 0.0)) {
    emlrtErrorWithMessageIdR2018a(sp, &j_emlrtRTEI,
                                  "MATLAB:odearguments:AbsTolNotPos",
                                  "MATLAB:odearguments:AbsTolNotPos", 0);
  }
  *threshold = options_AbsTol / *rtol;
  if (!(options_MaxStep_data[0] > 0.0)) {
    emlrtErrorWithMessageIdR2018a(sp, &i_emlrtRTEI,
                                  "MATLAB:odearguments:MaxStepLEzero",
                                  "MATLAB:odearguments:MaxStepLEzero", 0);
  }
  *hmax = muDoubleScalarMin(*htspan, options_MaxStep_data[0]);
  htryIn_data = muDoubleScalarAbs(options_InitialStep_data[0]);
  htry_size[0] = 1;
  htry_size[1] = 1;
  htry_data[0] = htryIn_data;
  if (!(htryIn_data > 0.0)) {
    emlrtErrorWithMessageIdR2018a(sp, &h_emlrtRTEI,
                                  "MATLAB:odearguments:InitialStepLEzero",
                                  "MATLAB:odearguments:InitialStepLEzero", 0);
  }
  return t0;
}

/* End of code generation (odearguments.c) */
