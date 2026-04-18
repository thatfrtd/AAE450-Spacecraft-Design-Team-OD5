/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * QLaw_transfer_inner.c
 *
 * Code generation for function 'QLaw_transfer_inner'
 *
 */

/* Include files */
#include "QLaw_transfer_inner.h"
#include "D_Q_Qdot_P_partial_Q_partial_oe_func.h"
#include "QLaw_transfer_inner_data.h"
#include "QLaw_transfer_inner_emxutil.h"
#include "QLaw_transfer_inner_internal_types.h"
#include "QLaw_transfer_inner_types.h"
#include "eml_int_forloop_overflow_check.h"
#include "linspace.h"
#include "ode45.h"
#include "odeset.h"
#include "rt_nonfinite.h"
#include "mwmathutil.h"
#include <emmintrin.h>
#include <string.h>

/* Type Definitions */
#ifndef struct_emxArray_real_T_1x1
#define struct_emxArray_real_T_1x1
struct emxArray_real_T_1x1 {
  real_T data[1];
};
#endif /* struct_emxArray_real_T_1x1 */
#ifndef typedef_emxArray_real_T_1x1
#define typedef_emxArray_real_T_1x1
typedef struct emxArray_real_T_1x1 emxArray_real_T_1x1;
#endif /* typedef_emxArray_real_T_1x1 */

#ifndef typedef_b_struct_T
#define typedef_b_struct_T
typedef struct {
  emxArray_real_T_1x1 InitialStep;
  emxArray_real_T_1x1 MaxStep;
} b_struct_T;
#endif /* typedef_b_struct_T */

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = {
    4,                     /* lineNo */
    "QLaw_transfer_inner", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QLaw_transfer_inner.m" /* pathName */
};

static emlrtRSInfo b_emlrtRSI = {
    5,                     /* lineNo */
    "QLaw_transfer_inner", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QLaw_transfer_inner.m" /* pathName */
};

static emlrtRSInfo c_emlrtRSI = {
    12,                    /* lineNo */
    "QLaw_transfer_inner", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QLaw_transfer_inner.m" /* pathName */
};

static emlrtRSInfo d_emlrtRSI = {
    13,                    /* lineNo */
    "QLaw_transfer_inner", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QLaw_transfer_inner.m" /* pathName */
};

static emlrtRSInfo e_emlrtRSI = {
    16,                    /* lineNo */
    "QLaw_transfer_inner", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QLaw_transfer_inner.m" /* pathName */
};

static emlrtRSInfo f_emlrtRSI = {
    44,                    /* lineNo */
    "QLaw_transfer_inner", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QLaw_transfer_inner.m" /* pathName */
};

static emlrtRSInfo g_emlrtRSI = {
    49,                    /* lineNo */
    "QLaw_transfer_inner", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QLaw_transfer_inner.m" /* pathName */
};

static emlrtRSInfo rb_emlrtRSI = {
    10,                    /* lineNo */
    "QLaw_thrust_mapping", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QLaw Functions\\QLaw_thrust_mapping.m" /* pathName */
};

static emlrtRSInfo sb_emlrtRSI = {
    19,                    /* lineNo */
    "Qdot_extremize_fast", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QLaw Functions\\Qdot_extremize_fast.m" /* pathName */
};

static emlrtRSInfo tb_emlrtRSI = {
    22,                    /* lineNo */
    "Qdot_extremize_fast", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QLaw Functions\\Qdot_extremize_fast.m" /* pathName */
};

static emlrtRSInfo ub_emlrtRSI = {
    28,                    /* lineNo */
    "Qdot_extremize_fast", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QLaw Functions\\Qdot_extremize_fast.m" /* pathName */
};

static emlrtRSInfo vb_emlrtRSI = {
    31,                    /* lineNo */
    "Qdot_extremize_fast", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QLaw Functions\\Qdot_extremize_fast.m" /* pathName */
};

static emlrtRSInfo xb_emlrtRSI = {
    16,                                                /* lineNo */
    "@(L)-norm(partial_Q_partial_oe*B_oe_func(oe,L))", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QLaw Functions\\Qdot_extremize_fast.m" /* pathName */
};

static emlrtRSInfo yb_emlrtRSI = {
    16,          /* lineNo */
    "B_oe_func", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QSymbolic\\B_oe_func.m" /* pathName */
};

static emlrtRSInfo ac_emlrtRSI = {
    17,          /* lineNo */
    "B_oe_func", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QSymbolic\\B_oe_func.m" /* pathName */
};

static emlrtRSInfo bc_emlrtRSI = {
    18,          /* lineNo */
    "B_oe_func", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QSymbolic\\B_oe_func.m" /* pathName */
};

static emlrtRSInfo cc_emlrtRSI = {
    19,          /* lineNo */
    "B_oe_func", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QSymbolic\\B_oe_func.m" /* pathName */
};

static emlrtRSInfo dc_emlrtRSI = {
    38,          /* lineNo */
    "B_oe_func", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QSymbolic\\B_oe_func.m" /* pathName */
};

static emlrtRSInfo ec_emlrtRSI = {
    39,          /* lineNo */
    "B_oe_func", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QSymbolic\\B_oe_func.m" /* pathName */
};

static emlrtRSInfo fc_emlrtRSI = {
    15,    /* lineNo */
    "min", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2025b\\toolbox\\eml\\lib\\matlab\\datafun\\min.m" /* pathName
                                                                        */
};

static emlrtRSInfo gc_emlrtRSI =
    {
        75,         /* lineNo */
        "minOrMax", /* fcnName */
        "C:\\Program "
        "Files\\MATLAB\\R2025b\\toolbox\\eml\\eml\\+coder\\+internal\\minOrMax."
        "m" /* pathName */
};

static emlrtRSInfo hc_emlrtRSI =
    {
        121,       /* lineNo */
        "minimum", /* fcnName */
        "C:\\Program "
        "Files\\MATLAB\\R2025b\\toolbox\\eml\\eml\\+coder\\+internal\\minOrMax."
        "m" /* pathName */
};

static emlrtRSInfo ic_emlrtRSI = {
    273,             /* lineNo */
    "unaryMinOrMax", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2025b\\toolbox\\eml\\eml\\+coder\\+"
    "internal\\unaryMinOrMax.m" /* pathName */
};

static emlrtRSInfo jc_emlrtRSI = {
    962,                    /* lineNo */
    "minRealVectorOmitNaN", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2025b\\toolbox\\eml\\eml\\+coder\\+"
    "internal\\unaryMinOrMax.m" /* pathName */
};

static emlrtRSInfo kc_emlrtRSI = {
    73,                      /* lineNo */
    "vectorMinOrMaxInPlace", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2025b\\toolbox\\eml\\eml\\+coder\\+"
    "internal\\vectorMinOrMaxInPlace.m" /* pathName */
};

static emlrtRSInfo lc_emlrtRSI = {
    65,                      /* lineNo */
    "vectorMinOrMaxInPlace", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2025b\\toolbox\\eml\\eml\\+coder\\+"
    "internal\\vectorMinOrMaxInPlace.m" /* pathName */
};

static emlrtRSInfo mc_emlrtRSI = {
    114,         /* lineNo */
    "findFirst", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2025b\\toolbox\\eml\\eml\\+coder\\+"
    "internal\\vectorMinOrMaxInPlace.m" /* pathName */
};

static emlrtRSInfo nc_emlrtRSI = {
    20,                               /* lineNo */
    "eml_int_forloop_overflow_check", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2025b\\toolbox\\eml\\lib\\matlab\\eml\\eml_int_forloop_"
    "overflow_check.m" /* pathName */
};

static emlrtRSInfo oc_emlrtRSI = {
    131,                        /* lineNo */
    "minOrMaxRealVectorKernel", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2025b\\toolbox\\eml\\eml\\+coder\\+"
    "internal\\vectorMinOrMaxInPlace.m" /* pathName */
};

static emlrtRSInfo pc_emlrtRSI = {
    15,    /* lineNo */
    "max", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2025b\\toolbox\\eml\\lib\\matlab\\datafun\\max.m" /* pathName
                                                                        */
};

static emlrtRSInfo qc_emlrtRSI =
    {
        73,         /* lineNo */
        "minOrMax", /* fcnName */
        "C:\\Program "
        "Files\\MATLAB\\R2025b\\toolbox\\eml\\eml\\+coder\\+internal\\minOrMax."
        "m" /* pathName */
};

static emlrtRSInfo rc_emlrtRSI =
    {
        108,       /* lineNo */
        "maximum", /* fcnName */
        "C:\\Program "
        "Files\\MATLAB\\R2025b\\toolbox\\eml\\eml\\+coder\\+internal\\minOrMax."
        "m" /* pathName */
};

static emlrtRSInfo sc_emlrtRSI = {
    255,             /* lineNo */
    "unaryMinOrMax", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2025b\\toolbox\\eml\\eml\\+coder\\+"
    "internal\\unaryMinOrMax.m" /* pathName */
};

static emlrtRSInfo tc_emlrtRSI = {
    966,                    /* lineNo */
    "maxRealVectorOmitNaN", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2025b\\toolbox\\eml\\eml\\+coder\\+"
    "internal\\unaryMinOrMax.m" /* pathName */
};

static emlrtRSInfo id_emlrtRSI = {
    49, /* lineNo */
    "@(t,x)[gauss_planetary_eqn(f0_modified_equinoctial(x,1),B_modified_"
    "equinoctial(x,1),a_control(x));mdot]", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QLaw_transfer_inner.m" /* pathName */
};

static emlrtRSInfo jd_emlrtRSI = {
    10,                        /* lineNo */
    "f0_modified_equinoctial", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Helper Functi"
    "ons\\Orbital Elements\\Gauss Planetary "
    "Equations\\f0_modified_equinoctial.m" /* pathName */
};

static emlrtRTEInfo emlrtRTEI = {
    21,                    /* lineNo */
    9,                     /* colNo */
    "Qdot_extremize_fast", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QLaw Functions\\Qdot_extremize_fast.m" /* pName */
};

static emlrtBCInfo emlrtBCI = {
    -1,                    /* iFirst */
    -1,                    /* iLast */
    22,                    /* lineNo */
    41,                    /* colNo */
    "L_start",             /* aName */
    "Qdot_extremize_fast", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QLaw Functions\\Qdot_extremize_fast.m", /* pName */
    0                                           /* checkKind */
};

static emlrtRTEInfo c_emlrtRTEI = {
    198,             /* lineNo */
    27,              /* colNo */
    "unaryMinOrMax", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2025b\\toolbox\\eml\\eml\\+coder\\+"
    "internal\\unaryMinOrMax.m" /* pName */
};

static emlrtDCInfo emlrtDCI = {
    20,                    /* lineNo */
    20,                    /* colNo */
    "Qdot_extremize_fast", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QLaw Functions\\Qdot_extremize_fast.m", /* pName */
    4                                           /* checkKind */
};

static emlrtDCInfo b_emlrtDCI = {
    20,                    /* lineNo */
    20,                    /* colNo */
    "Qdot_extremize_fast", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QLaw Functions\\Qdot_extremize_fast.m", /* pName */
    1                                           /* checkKind */
};

static emlrtBCInfo b_emlrtBCI = {
    -1,                    /* iFirst */
    -1,                    /* iLast */
    22,                    /* lineNo */
    16,                    /* colNo */
    "Qdot_start",          /* aName */
    "Qdot_extremize_fast", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QLaw Functions\\Qdot_extremize_fast.m", /* pName */
    0                                           /* checkKind */
};

static emlrtRTEInfo o_emlrtRTEI = {
    20,                    /* lineNo */
    1,                     /* colNo */
    "Qdot_extremize_fast", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QLaw Functions\\Qdot_extremize_fast.m" /* pName */
};

static emlrtRTEInfo p_emlrtRTEI = {
    19,                    /* lineNo */
    1,                     /* colNo */
    "Qdot_extremize_fast", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QLaw Functions\\Qdot_extremize_fast.m" /* pName */
};

/* Function Definitions */
void QLaw_transfer_inner(
    const emlrtStack *sp, const real_T x_me_mass_nd[7], real_T t_nd,
    const real_T oe_t[5], const real_T Q_params_W_oe[5], real_T Q_params_m,
    real_T Q_params_n, real_T Q_params_r, real_T Q_params_Theta_rot,
    real_T Q_params_eta_a_min, real_T Q_params_eta_r_min,
    real_T penalty_params_W_p, real_T penalty_params_k, real_T F_max_nd,
    boolean_T thrust_during_eclipse, real_T r_p_min_nd,
    const real_T nd_scalar[7], real_T char_star_t, real_T char_star_v,
    real_T spacecraft_params_Isp, real_T num_start_points,
    real_T integration_tolerance, real_T angular_step, emxArray_real_T *t_step,
    emxArray_real_T *x_step, real_T *Q, real_T *P, real_T *alpha, real_T *beta,
    real_T u_nd[3], boolean_T *not_coast)
{
  b_struct_T b_expl_temp;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack f_st;
  emlrtStack g_st;
  emlrtStack h_st;
  emlrtStack i_st;
  emlrtStack j_st;
  emlrtStack st;
  emxArray_real_T *L_start;
  emxArray_real_T *Qdot_start;
  struct_T expl_temp;
  real_T b_oe[15];
  real_T c_min_ab_Qdot_workspace_partial[5];
  real_T oe[5];
  real_T D[3];
  real_T b_t_nd[2];
  real_T Qdot;
  real_T b_Q;
  real_T b_alpha;
  real_T b_beta;
  real_T d = 0.0;
  real_T d1 = 0.0;
  real_T e;
  real_T t;
  real_T t10;
  real_T t15 = 0.0;
  real_T t20;
  real_T t21 = 0.0;
  real_T t26 = 0.0;
  real_T u_star_tmp;
  real_T w;
  real_T *L_start_data;
  real_T *Qdot_start_data;
  int32_T a;
  int32_T i;
  int32_T idx;
  int32_T k;
  int32_T last;
  int32_T s;
  boolean_T exitg1;
  boolean_T p = false;
  (void)Q_params_Theta_rot;
  (void)thrust_during_eclipse;
  (void)nd_scalar;
  (void)char_star_t;
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
  f_st.prev = &e_st;
  f_st.tls = e_st.tls;
  g_st.prev = &f_st;
  g_st.tls = f_st.tls;
  h_st.prev = &g_st;
  h_st.tls = g_st.tls;
  i_st.prev = &h_st;
  i_st.tls = h_st.tls;
  j_st.prev = &i_st;
  j_st.tls = i_st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  /* QLAW_TRANSFER Orbit transfer using Q-Law */
  /*  Create slow QLaw orbital elements (a, f, g, h, k) */
  st.site = &emlrtRSI;
  b_st.site = &h_emlrtRSI;
  c_st.site = &i_emlrtRSI;
  st.site = &emlrtRSI;
  b_st.site = &h_emlrtRSI;
  c_st.site = &i_emlrtRSI;
  st.site = &emlrtRSI;
  t10 = x_me_mass_nd[1] * x_me_mass_nd[1] + x_me_mass_nd[2] * x_me_mass_nd[2];
  e = muDoubleScalarSqrt(t10);
  st.site = &b_emlrtRSI;
  b_st.site = &h_emlrtRSI;
  c_st.site = &i_emlrtRSI;
  /*  a = p / (1 - e ^ 2) */
  w = (x_me_mass_nd[1] * muDoubleScalarCos(x_me_mass_nd[5]) + 1.0) +
      x_me_mass_nd[2] * muDoubleScalarSin(x_me_mass_nd[5]);
  /*  r = p / (1 + f .* cos(L) + g .* sin(L)) */
  oe[0] = x_me_mass_nd[0] / (1.0 - e * e);
  oe[1] = x_me_mass_nd[1];
  oe[2] = x_me_mass_nd[2];
  oe[3] = x_me_mass_nd[3];
  oe[4] = x_me_mass_nd[4];
  /*  Calculate control */
  st.site = &c_emlrtRSI;
  b_Q = c_D_Q_Qdot_P_partial_Q_partial_(
      &st, oe, x_me_mass_nd[5], oe_t, Q_params_W_oe, Q_params_m, Q_params_n,
      Q_params_r, F_max_nd, penalty_params_W_p, r_p_min_nd, penalty_params_k, D,
      &Qdot, P, c_min_ab_Qdot_workspace_partial);
  st.site = &d_emlrtRSI;
  /* QLAW_THRUST_MAPPING Calculate optimal thrusting angles to minimize
   * Q-function */
  /*    Detailed explanation goes here */
  /*  Q derivative w.r.t. RTN directions (D1, D2, D3) */
  /*  Force applied */
  b_alpha = muDoubleScalarAtan2(-D[0], -D[1]);
  b_st.site = &rb_emlrtRSI;
  c_st.site = &h_emlrtRSI;
  d_st.site = &i_emlrtRSI;
  b_st.site = &rb_emlrtRSI;
  c_st.site = &h_emlrtRSI;
  d_st.site = &i_emlrtRSI;
  b_st.site = &rb_emlrtRSI;
  b_beta =
      muDoubleScalarAtan(-D[2] / muDoubleScalarSqrt(D[0] * D[0] + D[1] * D[1]));
  u_star_tmp = muDoubleScalarCos(b_beta);
  /*  radial */
  /*  theta */
  /*  normal */
  /*   */
  /*  Qdot_ab = D * [cos(beta_star) * sin(alpha_star); ... */
  /*                 cos(beta_star) * cos(alpha_star); */
  /*                 sin(beta_star)]; */
  /*  Determine if thrusting should happen based on heuristic */
  st.site = &e_emlrtRSI;
  /* QDOT_EXTREMIZE Summary of this function goes here */
  /*    The minimum and maximum w.r.t. L of the minimum Qdot w.r.t. thrust  */
  /*  direction needs to be found so that the current L can be compared to the
   */
  /*  best and worst Ls to compute effeciences used to determine if coasting */
  /*  should be done or not. There is usually two minimums and maximums so the
   */
  /*  right ones should try to be chosen. Care must be taken to make this */
  /*  process fast. */
  /*  Know min w.r.t. alpha, beta of Qdot in terms of L */
  /*  Initialize search for L */
  emxInit_real_T(&st, &L_start, 2, &p_emlrtRTEI);
  b_st.site = &sb_emlrtRSI;
  linspace(&b_st, num_start_points + 1.0, L_start);
  L_start_data = L_start->data;
  if (!(num_start_points >= 0.0)) {
    emlrtNonNegativeCheckR2012b(num_start_points, &emlrtDCI, &st);
  }
  if (num_start_points != (int32_T)muDoubleScalarFloor(num_start_points)) {
    emlrtIntegerCheckR2012b(num_start_points, &b_emlrtDCI, &st);
  }
  emxInit_real_T(&st, &Qdot_start, 1, &o_emlrtRTEI);
  idx = (int32_T)num_start_points;
  a = Qdot_start->size[0];
  Qdot_start->size[0] = (int32_T)num_start_points;
  emxEnsureCapacity_real_T(&st, Qdot_start, a, &o_emlrtRTEI);
  Qdot_start_data = Qdot_start->data;
  for (s = 0; s < idx; s++) {
    Qdot_start_data[s] = 0.0;
  }
  emlrtForLoopVectorCheckR2021a(1.0, 1.0, num_start_points, mxDOUBLE_CLASS,
                                (int32_T)num_start_points, &emlrtRTEI, &st);
  if ((int32_T)num_start_points - 1 >= 0) {
    t15 = (x_me_mass_nd[3] * x_me_mass_nd[3] +
           x_me_mass_nd[4] * x_me_mass_nd[4]) +
          1.0;
    t20 = oe[0] * (t10 - 1.0);
    t21 = 1.0 / (t10 - 1.0);
    p = (-t20 < 0.0);
    t26 = muDoubleScalarSqrt(-t20);
    d = muDoubleScalarAtan(x_me_mass_nd[2] / x_me_mass_nd[1]);
    d1 = e;
    b_oe[3] = 0.0;
    b_oe[4] = 0.0;
    b_oe[8] = 0.0;
    b_oe[9] = 0.0;
    b_oe[10] = 0.0;
  }
  for (s = 0; s < idx; s++) {
    real_T t17;
    real_T t19;
    real_T t2;
    real_T t23;
    real_T t3;
    b_st.site = &tb_emlrtRSI;
    if ((s + 1 < 1) || (s + 1 > L_start->size[1])) {
      emlrtDynamicBoundsCheckR2012b(s + 1, 1, L_start->size[1], &emlrtBCI,
                                    &b_st);
    }
    c_st.site = &wb_emlrtRSI;
    d_st.site = &xb_emlrtRSI;
    /* B_oe_func */
    /*     B_oe = B_oe_func(IN1,L1) */
    /*     This function was generated by the Symbolic Math Toolbox
     * version 25.2. */
    /*     14-Feb-2026 17:40:22 */
    /* Inputs: [{oe}; {L}] */
    t2 = muDoubleScalarCos(L_start_data[s]);
    t3 = muDoubleScalarSin(L_start_data[s]);
    e_st.site = &yb_emlrtRSI;
    f_st.site = &i_emlrtRSI;
    e_st.site = &ac_emlrtRSI;
    f_st.site = &i_emlrtRSI;
    e_st.site = &bc_emlrtRSI;
    f_st.site = &i_emlrtRSI;
    e_st.site = &cc_emlrtRSI;
    f_st.site = &i_emlrtRSI;
    t20 = oe[1] * t2;
    t10 = oe[2] * t3;
    t17 = (t20 + t10) + 1.0;
    t19 = oe[4] * t2 - oe[3] * t3;
    t23 = 1.0 / t17;
    t20 = 1.0 / ((t20 * 2.0 + t10 * 2.0) + 2.0);
    e_st.site = &dc_emlrtRSI;
    if (p) {
      emlrtErrorWithMessageIdR2018a(
          &e_st, &d_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
          "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
    }
    e_st.site = &ec_emlrtRSI;
    b_oe[0] =
        oe[0] * t21 * t26 * muDoubleScalarSin(L_start_data[s] - d) * d1 * -2.0;
    b_oe[1] = t3 * t26;
    b_oe[2] = -t2 * t26;
    b_oe[5] = oe[0] * t17 * t21 * t26 * -2.0;
    t = t23 * t26;
    b_oe[6] = t * (oe[1] + t2 * (t17 + 1.0));
    b_oe[7] = t * (oe[2] + t3 * (t17 + 1.0));
    b_oe[11] = oe[2] * t19 * t23 * t26;
    b_oe[12] = -oe[1] * t19 * t23 * t26;
    b_oe[13] = t2 * t15 * t20 * t26;
    b_oe[14] = t3 * t15 * t20 * t26;
    memset(&D[0], 0, 3U * sizeof(real_T));
    t10 = 0.0;
    t2 = 3.3121686421112381E-170;
    for (k = 0; k < 3; k++) {
      t20 = D[k];
      for (i = 0; i < 5; i++) {
        t20 += c_min_ab_Qdot_workspace_partial[i] * b_oe[i + 5 * k];
      }
      D[k] = t20;
      t20 = muDoubleScalarAbs(t20);
      if (t20 > t2) {
        t = t2 / t20;
        t10 = t10 * t * t + 1.0;
        t2 = t20;
      } else {
        t = t20 / t2;
        t10 += t * t;
      }
    }
    t10 = t2 * muDoubleScalarSqrt(t10);
    if (((int32_T)((uint32_T)s + 1U) < 1) ||
        ((int32_T)((uint32_T)s + 1U) > Qdot_start->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)s + 1U), 1,
                                    Qdot_start->size[0], &b_emlrtBCI, &c_st);
    }
    Qdot_start_data[s] = -t10;
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(&st);
    }
  }
  emxFree_real_T(&st, &L_start);
  /*  Just look at the min and max of start point array, no iteration */
  /*  Min w.r.t. L (min w.r.t. alpha, beta (Qdot)) */
  b_st.site = &ub_emlrtRSI;
  c_st.site = &fc_emlrtRSI;
  d_st.site = &gc_emlrtRSI;
  e_st.site = &hc_emlrtRSI;
  if (Qdot_start->size[0] < 1) {
    emlrtErrorWithMessageIdR2018a(&e_st, &c_emlrtRTEI,
                                  "Coder:toolbox:eml_min_or_max_varDimZero",
                                  "Coder:toolbox:eml_min_or_max_varDimZero", 0);
  }
  f_st.site = &ic_emlrtRSI;
  g_st.site = &jc_emlrtRSI;
  last = Qdot_start->size[0];
  if (Qdot_start->size[0] <= 2) {
    if (Qdot_start->size[0] == 1) {
      t = Qdot_start_data[0];
    } else if ((Qdot_start_data[0] > Qdot_start_data[1]) ||
               (muDoubleScalarIsNaN(Qdot_start_data[0]) &&
                (!muDoubleScalarIsNaN(Qdot_start_data[1])))) {
      t = Qdot_start_data[1];
    } else {
      t = Qdot_start_data[0];
    }
  } else {
    h_st.site = &lc_emlrtRSI;
    if (!muDoubleScalarIsNaN(Qdot_start_data[0])) {
      idx = 1;
    } else {
      idx = 0;
      i_st.site = &mc_emlrtRSI;
      if (Qdot_start->size[0] > 2147483646) {
        j_st.site = &nc_emlrtRSI;
        check_forloop_overflow_error(&j_st);
      }
      a = 2;
      exitg1 = false;
      while ((!exitg1) && (a <= last)) {
        if (!muDoubleScalarIsNaN(Qdot_start_data[a - 1])) {
          idx = a;
          exitg1 = true;
        } else {
          a++;
        }
      }
    }
    if (idx == 0) {
      t = Qdot_start_data[0];
    } else {
      h_st.site = &kc_emlrtRSI;
      t = Qdot_start_data[idx - 1];
      a = idx + 1;
      i_st.site = &oc_emlrtRSI;
      if ((idx + 1 <= Qdot_start->size[0]) &&
          (Qdot_start->size[0] > 2147483646)) {
        j_st.site = &nc_emlrtRSI;
        check_forloop_overflow_error(&j_st);
      }
      for (s = a; s <= last; s++) {
        t20 = Qdot_start_data[s - 1];
        if (t > t20) {
          t = t20;
        }
      }
    }
  }
  /*  Max w.r.t. L (min w.r.t. alpha, beta (Qdot)) */
  b_st.site = &vb_emlrtRSI;
  c_st.site = &pc_emlrtRSI;
  d_st.site = &qc_emlrtRSI;
  e_st.site = &rc_emlrtRSI;
  f_st.site = &sc_emlrtRSI;
  g_st.site = &tc_emlrtRSI;
  if (Qdot_start->size[0] <= 2) {
    if (Qdot_start->size[0] == 1) {
      t20 = Qdot_start_data[0];
    } else if ((Qdot_start_data[0] < Qdot_start_data[1]) ||
               (muDoubleScalarIsNaN(Qdot_start_data[0]) &&
                (!muDoubleScalarIsNaN(Qdot_start_data[1])))) {
      t20 = Qdot_start_data[1];
    } else {
      t20 = Qdot_start_data[0];
    }
  } else {
    h_st.site = &lc_emlrtRSI;
    if (!muDoubleScalarIsNaN(Qdot_start_data[0])) {
      idx = 1;
    } else {
      idx = 0;
      i_st.site = &mc_emlrtRSI;
      if (Qdot_start->size[0] > 2147483646) {
        j_st.site = &nc_emlrtRSI;
        check_forloop_overflow_error(&j_st);
      }
      a = 2;
      exitg1 = false;
      while ((!exitg1) && (a <= last)) {
        if (!muDoubleScalarIsNaN(Qdot_start_data[a - 1])) {
          idx = a;
          exitg1 = true;
        } else {
          a++;
        }
      }
    }
    if (idx == 0) {
      t20 = Qdot_start_data[0];
    } else {
      h_st.site = &kc_emlrtRSI;
      t20 = Qdot_start_data[idx - 1];
      a = idx + 1;
      i_st.site = &oc_emlrtRSI;
      if ((idx + 1 <= Qdot_start->size[0]) &&
          (Qdot_start->size[0] > 2147483646)) {
        j_st.site = &nc_emlrtRSI;
        check_forloop_overflow_error(&j_st);
      }
      for (s = a; s <= last; s++) {
        t10 = Qdot_start_data[s - 1];
        if (t20 < t10) {
          t20 = t10;
        }
      }
    }
  }
  emxFree_real_T(&g_st, &Qdot_start);
  /* QLAW_EFFICIENCIES Summary of this function goes here */
  /*    Used as a heuristic to determine if thrusting should happen. If above a
   */
  /*    threshold then thrusting happens */
  /*  Absolute efficiency */
  /*  Relative efficiency */
  /*  if ~thrust_during_eclipse % Should also add state of charge based
   * thrusting */
  /*      % Eclipse */
  /*      R_E = 6378.1; % [km] Earth radius */
  /*      AU = 149597898; % [km] astronautical unit, Earth-Sun difference */
  /*      mu_sun = 132712440017.99; % [km3 / s2] Sun gravitational parameter */
  /*      n_E = sqrt(mu_sun / AU ^ 3); % [rad / s] Earth mean motion */
  /*   */
  /*      x_cartesian = x_me_nd_to_cartesian(x_me_mass_nd(1:6), nd_scalar,
   * Q_params_Theta_rot, mu); */
  /*      rvec = x_cartesian(1:3); */
  /*   */
  /*      [~, eclipsed] = check_eclipse(t_nd * char_star_t, rvec, R_E, AU, n_E);
   * % Need to properly get sun position and Earth tilt */
  /*   */
  /*      eclipse_no_thrust = eclipsed(iter); */
  /*  else % Thrust during eclispe */
  /*  end */
  if ((Qdot / t >= Q_params_eta_a_min) &&
      ((Qdot - t20) / (t - t20) >= Q_params_eta_r_min)) {
    p = true;
  } else {
    p = false;
  }
  u_nd[0] = F_max_nd * (u_star_tmp * muDoubleScalarSin(b_alpha)) * (real_T)p;
  u_nd[1] = F_max_nd * (u_star_tmp * muDoubleScalarCos(b_alpha)) * (real_T)p;
  u_nd[2] = F_max_nd * muDoubleScalarSin(b_beta) * (real_T)p;
  /* a_control = @(x) u_nd(:, iter) / x_me_mass_nd(7, end); % a = F / m */
  /*  [km / s2] */
  /*  Propagate orbit */
  st.site = &f_emlrtRSI;
  b_st.site = &h_emlrtRSI;
  c_st.site = &i_emlrtRSI;
  t20 = muDoubleScalarPower(x_me_mass_nd[0] / w, 3.0);
  st.site = &f_emlrtRSI;
  if (t20 < 0.0) {
    emlrtErrorWithMessageIdR2018a(
        &st, &d_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
        "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
  }
  t20 = muDoubleScalarSqrt(t20);
  st.site = &f_emlrtRSI;
  if (w < 0.0) {
    emlrtErrorWithMessageIdR2018a(
        &st, &d_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
        "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
  }
  t20 = angular_step * t20 * muDoubleScalarSqrt(w) / (e + 1.0);
  expl_temp = odeset(integration_tolerance, integration_tolerance);
  b_expl_temp.InitialStep.data[0] = t20;
  b_expl_temp.MaxStep.data[0] = t20;
  b_t_nd[0] = t_nd;
  b_t_nd[1] = t_nd + t20;
  st.site = &g_emlrtRSI;
  ode45(&st, u_nd,
        -F_max_nd / (spacecraft_params_Isp * 0.00981 / char_star_v) * (real_T)p,
        b_t_nd, x_me_mass_nd, expl_temp.AbsTol, b_expl_temp.InitialStep.data,
        b_expl_temp.MaxStep.data, expl_temp.RelTol, t_step, x_step);
  *Q = b_Q;
  *alpha = b_alpha;
  *beta = b_beta;
  *not_coast = p;
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

void QLaw_transfer_inner_anonFcn2(const emlrtStack *sp,
                                  const real_T a_control_workspace_u_nd[3],
                                  real_T mdot, const real_T x[7],
                                  real_T varargout_1[7])
{
  __m128d r;
  __m128d r1;
  emlrtStack b_st;
  emlrtStack st;
  real_T c_a[18];
  real_T d_a[6];
  real_T dv[6];
  real_T b_a_control_workspace_u_nd[3];
  real_T a;
  real_T a_tmp;
  real_T b_a;
  real_T b_a_tmp;
  real_T c_a_tmp;
  real_T cons_1;
  real_T d_a_tmp;
  int32_T i;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  st.site = &id_emlrtRSI;
  a_tmp = muDoubleScalarSin(x[5]);
  b_a_tmp = muDoubleScalarCos(x[5]);
  c_a_tmp = (x[1] * b_a_tmp + 1.0) + x[2] * a_tmp;
  a = c_a_tmp / x[0];
  b_st.site = &jd_emlrtRSI;
  if (x[0] < 0.0) {
    emlrtErrorWithMessageIdR2018a(
        &b_st, &d_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
        "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
  }
  st.site = &id_emlrtRSI;
  cons_1 = (x[3] * a_tmp - x[4] * b_a_tmp) / c_a_tmp;
  b_a = muDoubleScalarSqrt(x[0]);
  /*   */
  /*  % Make B convert disturbance into RTN frame from cartesian frame */
  /*  i = atan2(2 * sqrt(h ^ 2 + k ^ 2), 1 - h ^ 2 - k ^ 2); */
  /*  Omega = atan2(k, h); */
  /*  omega = atan2(g * h - f * k, f * h + g * k); */
  /*  nu = L - (Omega + omega); */
  /*  B = B * cartesian_to_RTN_DCM(i, Omega, omega, nu)'; */
  d_a_tmp = b_a * 0.0;
  c_a[0] = d_a_tmp;
  c_a[6] = b_a * (2.0 * x[0] / c_a_tmp);
  c_a[12] = d_a_tmp;
  c_a[1] = b_a * a_tmp;
  c_a[7] = b_a * (((c_a_tmp + 1.0) * b_a_tmp + x[1]) / c_a_tmp);
  c_a[13] = b_a * (-x[2] * cons_1);
  c_a[2] = b_a * -b_a_tmp;
  c_a[8] = b_a * (((c_a_tmp + 1.0) * a_tmp + x[2]) / c_a_tmp);
  c_a[14] = b_a * (x[1] * cons_1);
  c_a[3] = d_a_tmp;
  c_a[9] = d_a_tmp;
  c_a_tmp = ((x[3] * x[3] + 1.0) + x[4] * x[4]) / (2.0 * c_a_tmp);
  c_a[15] = b_a * (c_a_tmp * b_a_tmp);
  c_a[4] = d_a_tmp;
  c_a[10] = d_a_tmp;
  c_a[16] = b_a * (c_a_tmp * a_tmp);
  c_a[5] = d_a_tmp;
  c_a[11] = d_a_tmp;
  c_a[17] = b_a * cons_1;
  _mm_storeu_pd(&b_a_control_workspace_u_nd[0],
                _mm_div_pd(_mm_loadu_pd(&a_control_workspace_u_nd[0]),
                           _mm_set1_pd(x[6])));
  b_a_control_workspace_u_nd[2] = a_control_workspace_u_nd[2] / x[6];
  for (i = 0; i < 5; i++) {
    dv[i] = 0.0;
  }
  dv[5] = b_a * (a * a);
  memset(&d_a[0], 0, 6U * sizeof(real_T));
  for (i = 0; i < 3; i++) {
    __m128d r2;
    r = _mm_loadu_pd(&c_a[6 * i]);
    r1 = _mm_loadu_pd(&d_a[0]);
    r2 = _mm_set1_pd(b_a_control_workspace_u_nd[i]);
    _mm_storeu_pd(&d_a[0], _mm_add_pd(r1, _mm_mul_pd(r, r2)));
    r = _mm_loadu_pd(&c_a[6 * i + 2]);
    r1 = _mm_loadu_pd(&d_a[2]);
    _mm_storeu_pd(&d_a[2], _mm_add_pd(r1, _mm_mul_pd(r, r2)));
    r = _mm_loadu_pd(&c_a[6 * i + 4]);
    r1 = _mm_loadu_pd(&d_a[4]);
    _mm_storeu_pd(&d_a[4], _mm_add_pd(r1, _mm_mul_pd(r, r2)));
  }
  r = _mm_loadu_pd(&dv[0]);
  r1 = _mm_loadu_pd(&d_a[0]);
  _mm_storeu_pd(&varargout_1[0], _mm_add_pd(r, r1));
  r = _mm_loadu_pd(&dv[2]);
  r1 = _mm_loadu_pd(&d_a[2]);
  _mm_storeu_pd(&varargout_1[2], _mm_add_pd(r, r1));
  r = _mm_loadu_pd(&dv[4]);
  r1 = _mm_loadu_pd(&d_a[4]);
  _mm_storeu_pd(&varargout_1[4], _mm_add_pd(r, r1));
  varargout_1[6] = mdot;
}

/* End of code generation (QLaw_transfer_inner.c) */
