/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_QLaw_transfer_inner_api.c
 *
 * Code generation for function '_coder_QLaw_transfer_inner_api'
 *
 */

/* Include files */
#include "_coder_QLaw_transfer_inner_api.h"
#include "QLaw_transfer_inner.h"
#include "QLaw_transfer_inner_data.h"
#include "QLaw_transfer_inner_emxutil.h"
#include "QLaw_transfer_inner_mexutil.h"
#include "QLaw_transfer_inner_types.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRTEInfo w_emlrtRTEI = {
    1,                                /* lineNo */
    1,                                /* colNo */
    "_coder_QLaw_transfer_inner_api", /* fName */
    ""                                /* pName */
};

/* Function Declarations */
static const mxArray *b_emlrt_marshallOut(emxArray_real_T *u);

static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[7];

static const mxArray *c_emlrt_marshallOut(emxArray_real_T *u);

static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[7];

static const mxArray *d_emlrt_marshallOut(real_T u[3]);

static real_T e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                 const char_T *identifier);

static const mxArray *e_emlrt_marshallOut(const boolean_T u);

static real_T f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                 const emlrtMsgIdentifier *parentId);

static real_T (*g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[5];

static real_T (*h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[5];

static boolean_T i_emlrt_marshallIn(const emlrtStack *sp,
                                    const mxArray *nullptr,
                                    const char_T *identifier);

static boolean_T j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                    const emlrtMsgIdentifier *parentId);

static real_T (*l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[7];

static real_T m_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                 const emlrtMsgIdentifier *msgId);

static real_T (*n_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[5];

static boolean_T o_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                    const emlrtMsgIdentifier *msgId);

/* Function Definitions */
static const mxArray *b_emlrt_marshallOut(emxArray_real_T *u)
{
  static const int32_T i = 0;
  const mxArray *m;
  const mxArray *y;
  real_T *u_data;
  u_data = u->data;
  y = NULL;
  m = emlrtCreateNumericArray(1, (const void *)&i, mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m, &u_data[0]);
  emlrtSetDimensions((mxArray *)m, &u->size[0], 1);
  u->canFreeData = false;
  emlrtAssign(&y, m);
  return y;
}

static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[7]
{
  emlrtMsgIdentifier thisId;
  real_T(*y)[7];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = d_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

static const mxArray *c_emlrt_marshallOut(emxArray_real_T *u)
{
  static const int32_T iv[2] = {0, 0};
  const mxArray *m;
  const mxArray *y;
  real_T *u_data;
  u_data = u->data;
  y = NULL;
  m = emlrtCreateNumericArray(2, (const void *)&iv[0], mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m, &u_data[0]);
  emlrtSetDimensions((mxArray *)m, &u->size[0], 2);
  u->canFreeData = false;
  emlrtAssign(&y, m);
  return y;
}

static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[7]
{
  real_T(*y)[7];
  y = l_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static const mxArray *d_emlrt_marshallOut(real_T u[3])
{
  static const int32_T i = 0;
  static const int32_T i1 = 3;
  const mxArray *m;
  const mxArray *y;
  y = NULL;
  m = emlrtCreateNumericArray(1, (const void *)&i, mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m, &u[0]);
  emlrtSetDimensions((mxArray *)m, &i1, 1);
  emlrtAssign(&y, m);
  return y;
}

static real_T e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                 const char_T *identifier)
{
  emlrtMsgIdentifier thisId;
  real_T y;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = f_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

static const mxArray *e_emlrt_marshallOut(const boolean_T u)
{
  const mxArray *m;
  const mxArray *y;
  y = NULL;
  m = emlrtCreateLogicalScalar(u);
  emlrtAssign(&y, m);
  return y;
}

static real_T f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                 const emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = m_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T (*g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[5]
{
  emlrtMsgIdentifier thisId;
  real_T(*y)[5];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = h_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

static real_T (*h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[5]
{
  real_T(*y)[5];
  y = n_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static boolean_T i_emlrt_marshallIn(const emlrtStack *sp,
                                    const mxArray *nullptr,
                                    const char_T *identifier)
{
  emlrtMsgIdentifier thisId;
  boolean_T y;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = j_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

static boolean_T j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                    const emlrtMsgIdentifier *parentId)
{
  boolean_T y;
  y = o_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T (*l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[7]
{
  static const int32_T dims = 7;
  real_T(*ret)[7];
  int32_T i;
  boolean_T b = false;
  emlrtCheckVsBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 1U,
                            (const void *)&dims, &b, &i);
  ret = (real_T(*)[7])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T m_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                 const emlrtMsgIdentifier *msgId)
{
  static const int32_T dims = 0;
  real_T ret;
  emlrtCheckBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 0U,
                          (const void *)&dims);
  ret = *(real_T *)emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T (*n_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[5]
{
  static const int32_T dims = 5;
  real_T(*ret)[5];
  int32_T i;
  boolean_T b = false;
  emlrtCheckVsBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 1U,
                            (const void *)&dims, &b, &i);
  ret = (real_T(*)[5])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static boolean_T o_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                    const emlrtMsgIdentifier *msgId)
{
  static const int32_T dims = 0;
  boolean_T ret;
  emlrtCheckBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "logical", false, 0U,
                          (const void *)&dims);
  ret = *emlrtMxGetLogicals(src);
  emlrtDestroyArray(&src);
  return ret;
}

void QLaw_transfer_inner_api(const mxArray *const prhs[22], int32_T nlhs,
                             const mxArray *plhs[8])
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  emxArray_real_T *t_step;
  emxArray_real_T *x_step;
  real_T(*nd_scalar)[7];
  real_T(*x_me_mass_nd)[7];
  real_T(*Q_params_W_oe)[5];
  real_T(*oe_t)[5];
  real_T(*u_nd)[3];
  real_T F_max_nd;
  real_T P;
  real_T Q;
  real_T Q_params_Theta_rot;
  real_T Q_params_eta_a_min;
  real_T Q_params_eta_r_min;
  real_T Q_params_m;
  real_T Q_params_n;
  real_T Q_params_r;
  real_T alpha;
  real_T angular_step;
  real_T beta;
  real_T char_star_t;
  real_T char_star_v;
  real_T integration_tolerance;
  real_T num_start_points;
  real_T penalty_params_W_p;
  real_T penalty_params_k;
  real_T r_p_min_nd;
  real_T spacecraft_params_Isp;
  real_T t_nd;
  boolean_T not_coast;
  boolean_T thrust_during_eclipse;
  st.tls = emlrtRootTLSGlobal;
  u_nd = (real_T(*)[3])mxMalloc(sizeof(real_T[3]));
  emlrtHeapReferenceStackEnterFcnR2012b(&st);
  /* Marshall function inputs */
  x_me_mass_nd = c_emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "x_me_mass_nd");
  t_nd = e_emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "t_nd");
  oe_t = g_emlrt_marshallIn(&st, emlrtAlias(prhs[2]), "oe_t");
  Q_params_W_oe = g_emlrt_marshallIn(&st, emlrtAlias(prhs[3]), "Q_params_W_oe");
  Q_params_m = e_emlrt_marshallIn(&st, emlrtAliasP(prhs[4]), "Q_params_m");
  Q_params_n = e_emlrt_marshallIn(&st, emlrtAliasP(prhs[5]), "Q_params_n");
  Q_params_r = e_emlrt_marshallIn(&st, emlrtAliasP(prhs[6]), "Q_params_r");
  Q_params_Theta_rot =
      e_emlrt_marshallIn(&st, emlrtAliasP(prhs[7]), "Q_params_Theta_rot");
  Q_params_eta_a_min =
      e_emlrt_marshallIn(&st, emlrtAliasP(prhs[8]), "Q_params_eta_a_min");
  Q_params_eta_r_min =
      e_emlrt_marshallIn(&st, emlrtAliasP(prhs[9]), "Q_params_eta_r_min");
  penalty_params_W_p =
      e_emlrt_marshallIn(&st, emlrtAliasP(prhs[10]), "penalty_params_W_p");
  penalty_params_k =
      e_emlrt_marshallIn(&st, emlrtAliasP(prhs[11]), "penalty_params_k");
  F_max_nd = e_emlrt_marshallIn(&st, emlrtAliasP(prhs[12]), "F_max_nd");
  thrust_during_eclipse =
      i_emlrt_marshallIn(&st, emlrtAliasP(prhs[13]), "thrust_during_eclipse");
  r_p_min_nd = e_emlrt_marshallIn(&st, emlrtAliasP(prhs[14]), "r_p_min_nd");
  nd_scalar = c_emlrt_marshallIn(&st, emlrtAlias(prhs[15]), "nd_scalar");
  char_star_t = e_emlrt_marshallIn(&st, emlrtAliasP(prhs[16]), "char_star_t");
  char_star_v = e_emlrt_marshallIn(&st, emlrtAliasP(prhs[17]), "char_star_v");
  spacecraft_params_Isp =
      e_emlrt_marshallIn(&st, emlrtAliasP(prhs[18]), "spacecraft_params_Isp");
  num_start_points =
      e_emlrt_marshallIn(&st, emlrtAliasP(prhs[19]), "num_start_points");
  integration_tolerance =
      e_emlrt_marshallIn(&st, emlrtAliasP(prhs[20]), "integration_tolerance");
  angular_step = e_emlrt_marshallIn(&st, emlrtAliasP(prhs[21]), "angular_step");
  /* Invoke the target function */
  emxInit_real_T(&st, &t_step, 1, &w_emlrtRTEI);
  emxInit_real_T(&st, &x_step, 2, &w_emlrtRTEI);
  QLaw_transfer_inner(
      &st, *x_me_mass_nd, t_nd, *oe_t, *Q_params_W_oe, Q_params_m, Q_params_n,
      Q_params_r, Q_params_Theta_rot, Q_params_eta_a_min, Q_params_eta_r_min,
      penalty_params_W_p, penalty_params_k, F_max_nd, thrust_during_eclipse,
      r_p_min_nd, *nd_scalar, char_star_t, char_star_v, spacecraft_params_Isp,
      num_start_points, integration_tolerance, angular_step, t_step, x_step, &Q,
      &P, &alpha, &beta, *u_nd, &not_coast);
  /* Marshall function outputs */
  t_step->canFreeData = false;
  plhs[0] = b_emlrt_marshallOut(t_step);
  emxFree_real_T(&st, &t_step);
  if (nlhs > 1) {
    x_step->canFreeData = false;
    plhs[1] = c_emlrt_marshallOut(x_step);
  }
  emxFree_real_T(&st, &x_step);
  if (nlhs > 2) {
    plhs[2] = emlrt_marshallOut(Q);
  }
  if (nlhs > 3) {
    plhs[3] = emlrt_marshallOut(P);
  }
  if (nlhs > 4) {
    plhs[4] = emlrt_marshallOut(alpha);
  }
  if (nlhs > 5) {
    plhs[5] = emlrt_marshallOut(beta);
  }
  if (nlhs > 6) {
    plhs[6] = d_emlrt_marshallOut(*u_nd);
  }
  if (nlhs > 7) {
    plhs[7] = e_emlrt_marshallOut(not_coast);
  }
  emlrtHeapReferenceStackLeaveFcnR2012b(&st);
}

/* End of code generation (_coder_QLaw_transfer_inner_api.c) */
