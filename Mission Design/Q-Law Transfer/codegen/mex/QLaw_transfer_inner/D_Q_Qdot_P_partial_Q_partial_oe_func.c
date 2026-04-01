/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * D_Q_Qdot_P_partial_Q_partial_oe_func.c
 *
 * Code generation for function 'D_Q_Qdot_P_partial_Q_partial_oe_func'
 *
 */

/* Include files */
#include "D_Q_Qdot_P_partial_Q_partial_oe_func.h"
#include "QLaw_transfer_inner_data.h"
#include "rt_nonfinite.h"
#include "mwmathutil.h"

/* Variable Definitions */
static emlrtRSInfo j_emlrtRSI = {
    38,                                     /* lineNo */
    "D_Q_Qdot_P_partial_Q_partial_oe_func", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QSymbolic\\D_Q_Qdot_P_partial_Q_partial_oe_func.m" /* pathName */
};

static emlrtRSInfo k_emlrtRSI = {
    39,                                     /* lineNo */
    "D_Q_Qdot_P_partial_Q_partial_oe_func", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QSymbolic\\D_Q_Qdot_P_partial_Q_partial_oe_func.m" /* pathName */
};

static emlrtRSInfo l_emlrtRSI = {
    40,                                     /* lineNo */
    "D_Q_Qdot_P_partial_Q_partial_oe_func", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QSymbolic\\D_Q_Qdot_P_partial_Q_partial_oe_func.m" /* pathName */
};

static emlrtRSInfo m_emlrtRSI = {
    41,                                     /* lineNo */
    "D_Q_Qdot_P_partial_Q_partial_oe_func", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QSymbolic\\D_Q_Qdot_P_partial_Q_partial_oe_func.m" /* pathName */
};

static emlrtRSInfo n_emlrtRSI = {
    42,                                     /* lineNo */
    "D_Q_Qdot_P_partial_Q_partial_oe_func", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QSymbolic\\D_Q_Qdot_P_partial_Q_partial_oe_func.m" /* pathName */
};

static emlrtRSInfo o_emlrtRSI = {
    64,                                     /* lineNo */
    "D_Q_Qdot_P_partial_Q_partial_oe_func", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QSymbolic\\D_Q_Qdot_P_partial_Q_partial_oe_func.m" /* pathName */
};

static emlrtRSInfo p_emlrtRSI = {
    65,                                     /* lineNo */
    "D_Q_Qdot_P_partial_Q_partial_oe_func", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QSymbolic\\D_Q_Qdot_P_partial_Q_partial_oe_func.m" /* pathName */
};

static emlrtRSInfo q_emlrtRSI = {
    80,                                     /* lineNo */
    "D_Q_Qdot_P_partial_Q_partial_oe_func", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QSymbolic\\D_Q_Qdot_P_partial_Q_partial_oe_func.m" /* pathName */
};

static emlrtRSInfo r_emlrtRSI = {
    92,                                     /* lineNo */
    "D_Q_Qdot_P_partial_Q_partial_oe_func", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QSymbolic\\D_Q_Qdot_P_partial_Q_partial_oe_func.m" /* pathName */
};

static emlrtRSInfo s_emlrtRSI = {
    93,                                     /* lineNo */
    "D_Q_Qdot_P_partial_Q_partial_oe_func", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QSymbolic\\D_Q_Qdot_P_partial_Q_partial_oe_func.m" /* pathName */
};

static emlrtRSInfo t_emlrtRSI = {
    94,                                     /* lineNo */
    "D_Q_Qdot_P_partial_Q_partial_oe_func", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QSymbolic\\D_Q_Qdot_P_partial_Q_partial_oe_func.m" /* pathName */
};

static emlrtRSInfo u_emlrtRSI = {
    95,                                     /* lineNo */
    "D_Q_Qdot_P_partial_Q_partial_oe_func", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QSymbolic\\D_Q_Qdot_P_partial_Q_partial_oe_func.m" /* pathName */
};

static emlrtRSInfo v_emlrtRSI = {
    96,                                     /* lineNo */
    "D_Q_Qdot_P_partial_Q_partial_oe_func", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QSymbolic\\D_Q_Qdot_P_partial_Q_partial_oe_func.m" /* pathName */
};

static emlrtRSInfo w_emlrtRSI = {
    99,                                     /* lineNo */
    "D_Q_Qdot_P_partial_Q_partial_oe_func", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QSymbolic\\D_Q_Qdot_P_partial_Q_partial_oe_func.m" /* pathName */
};

static emlrtRSInfo x_emlrtRSI = {
    100,                                    /* lineNo */
    "D_Q_Qdot_P_partial_Q_partial_oe_func", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QSymbolic\\D_Q_Qdot_P_partial_Q_partial_oe_func.m" /* pathName */
};

static emlrtRSInfo y_emlrtRSI = {
    101,                                    /* lineNo */
    "D_Q_Qdot_P_partial_Q_partial_oe_func", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QSymbolic\\D_Q_Qdot_P_partial_Q_partial_oe_func.m" /* pathName */
};

static emlrtRSInfo ab_emlrtRSI = {
    109,                                    /* lineNo */
    "D_Q_Qdot_P_partial_Q_partial_oe_func", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QSymbolic\\D_Q_Qdot_P_partial_Q_partial_oe_func.m" /* pathName */
};

static emlrtRSInfo bb_emlrtRSI = {
    110,                                    /* lineNo */
    "D_Q_Qdot_P_partial_Q_partial_oe_func", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QSymbolic\\D_Q_Qdot_P_partial_Q_partial_oe_func.m" /* pathName */
};

static emlrtRSInfo cb_emlrtRSI = {
    119,                                    /* lineNo */
    "D_Q_Qdot_P_partial_Q_partial_oe_func", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QSymbolic\\D_Q_Qdot_P_partial_Q_partial_oe_func.m" /* pathName */
};

static emlrtRSInfo db_emlrtRSI = {
    130,                                    /* lineNo */
    "D_Q_Qdot_P_partial_Q_partial_oe_func", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QSymbolic\\D_Q_Qdot_P_partial_Q_partial_oe_func.m" /* pathName */
};

static emlrtRSInfo eb_emlrtRSI = {
    132,                                    /* lineNo */
    "D_Q_Qdot_P_partial_Q_partial_oe_func", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QSymbolic\\D_Q_Qdot_P_partial_Q_partial_oe_func.m" /* pathName */
};

static emlrtRSInfo fb_emlrtRSI = {
    133,                                    /* lineNo */
    "D_Q_Qdot_P_partial_Q_partial_oe_func", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QSymbolic\\D_Q_Qdot_P_partial_Q_partial_oe_func.m" /* pathName */
};

static emlrtRSInfo gb_emlrtRSI = {
    143,                                    /* lineNo */
    "D_Q_Qdot_P_partial_Q_partial_oe_func", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QSymbolic\\D_Q_Qdot_P_partial_Q_partial_oe_func.m" /* pathName */
};

static emlrtRSInfo hb_emlrtRSI = {
    144,                                    /* lineNo */
    "D_Q_Qdot_P_partial_Q_partial_oe_func", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QSymbolic\\D_Q_Qdot_P_partial_Q_partial_oe_func.m" /* pathName */
};

static emlrtRSInfo ib_emlrtRSI = {
    145,                                    /* lineNo */
    "D_Q_Qdot_P_partial_Q_partial_oe_func", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QSymbolic\\D_Q_Qdot_P_partial_Q_partial_oe_func.m" /* pathName */
};

static emlrtRSInfo jb_emlrtRSI = {
    146,                                    /* lineNo */
    "D_Q_Qdot_P_partial_Q_partial_oe_func", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QSymbolic\\D_Q_Qdot_P_partial_Q_partial_oe_func.m" /* pathName */
};

static emlrtRSInfo kb_emlrtRSI = {
    160,                                    /* lineNo */
    "D_Q_Qdot_P_partial_Q_partial_oe_func", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QSymbolic\\D_Q_Qdot_P_partial_Q_partial_oe_func.m" /* pathName */
};

static emlrtRSInfo lb_emlrtRSI = {
    161,                                    /* lineNo */
    "D_Q_Qdot_P_partial_Q_partial_oe_func", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QSymbolic\\D_Q_Qdot_P_partial_Q_partial_oe_func.m" /* pathName */
};

static emlrtRSInfo mb_emlrtRSI = {
    245,                                    /* lineNo */
    "D_Q_Qdot_P_partial_Q_partial_oe_func", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QSymbolic\\D_Q_Qdot_P_partial_Q_partial_oe_func.m" /* pathName */
};

static emlrtRSInfo nb_emlrtRSI = {
    252,                                    /* lineNo */
    "D_Q_Qdot_P_partial_Q_partial_oe_func", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QSymbolic\\D_Q_Qdot_P_partial_Q_partial_oe_func.m" /* pathName */
};

static emlrtRSInfo ob_emlrtRSI = {
    255,                                    /* lineNo */
    "D_Q_Qdot_P_partial_Q_partial_oe_func", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QSymbolic\\D_Q_Qdot_P_partial_Q_partial_oe_func.m" /* pathName */
};

static emlrtRSInfo pb_emlrtRSI = {
    259,                                    /* lineNo */
    "D_Q_Qdot_P_partial_Q_partial_oe_func", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QSymbolic\\D_Q_Qdot_P_partial_Q_partial_oe_func.m" /* pathName */
};

static emlrtRSInfo qb_emlrtRSI = {
    264,                                    /* lineNo */
    "D_Q_Qdot_P_partial_Q_partial_oe_func", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "450\\AAE450-Spacecraft-Design-Team-OD5\\Mission Design\\Q-Law Transfe"
    "r\\QSymbolic\\D_Q_Qdot_P_partial_Q_partial_oe_func.m" /* pathName */
};

/* Function Definitions */
real_T c_D_Q_Qdot_P_partial_Q_partial_(
    const emlrtStack *sp, const real_T in1[5], real_T L1, const real_T in3[5],
    const real_T in4[5], real_T m1, real_T n1, real_T r1, real_T F_max1,
    real_T W_p1, real_T r_p_min1, real_T k1, real_T D[3], real_T *Qdot,
    real_T *P, real_T partial_Q_partial_oe[5])
{
  emlrtStack b_st;
  emlrtStack st;
  real_T Q;
  real_T t102;
  real_T t103;
  real_T t105;
  real_T t106;
  real_T t108;
  real_T t109;
  real_T t110;
  real_T t111;
  real_T t112;
  real_T t113;
  real_T t114;
  real_T t117;
  real_T t124;
  real_T t125;
  real_T t126;
  real_T t16;
  real_T t17;
  real_T t175;
  real_T t176;
  real_T t176_tmp_tmp;
  real_T t177;
  real_T t178;
  real_T t179;
  real_T t179_tmp;
  real_T t180;
  real_T t181;
  real_T t182;
  real_T t2;
  real_T t215;
  real_T t216;
  real_T t231;
  real_T t239;
  real_T t24;
  real_T t25;
  real_T t3;
  real_T t37;
  real_T t38;
  real_T t40;
  real_T t42;
  real_T t43;
  real_T t44;
  real_T t52;
  real_T t58;
  real_T t6;
  real_T t69;
  real_T t7;
  real_T t70;
  real_T t71;
  real_T t72;
  real_T t73;
  real_T t74;
  real_T t77;
  real_T t79;
  real_T t8;
  real_T t80;
  real_T t81;
  real_T t82;
  real_T t83;
  real_T t85;
  real_T t86;
  real_T t87;
  real_T t89;
  real_T t90;
  real_T t92;
  real_T t96;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  /* D_Q_Qdot_P_partial_Q_partial_oe_func */
  /*     [D,Q,Qdot,P,partial_Q_partial_oe] =
   * D_Q_Qdot_P_partial_Q_partial_oe_func(IN1,L1,IN3,IN4,M1,N1,R1,F_max1,W_p1,R_P_MIN1,K1)
   */
  /*     This function was generated by the Symbolic Math Toolbox version 25.2.
   */
  /*     14-Feb-2026 18:39:07 */
  /* Inputs: [{oe}; {L}; {oe_t}; {W_oe}; {m}; {n}; {r}; {F_max}; {W_p};
   * {r_p_min}; {k_p}] */
  t2 = muDoubleScalarCos(L1);
  t3 = muDoubleScalarSin(L1);
  t6 = in1[0] * 2.0;
  t7 = in1[1] * 2.0;
  t8 = in1[2] * 2.0;
  st.site = &j_emlrtRSI;
  t16 = in1[1] * in1[1];
  st.site = &k_emlrtRSI;
  t17 = in1[2] * in1[2];
  st.site = &l_emlrtRSI;
  st.site = &m_emlrtRSI;
  st.site = &n_emlrtRSI;
  t24 = 1.0 / (F_max1 * F_max1);
  t25 = 1.0 / m1;
  t37 = 1.0 / in1[0];
  t42 = 1.0 / in3[0];
  t43 = 1.0 / r1;
  t44 = 1.0 / r_p_min1;
  st.site = &o_emlrtRSI;
  t38 = t37 * t37;
  st.site = &p_emlrtRSI;
  t40 = muDoubleScalarPower(t37, 3.0);
  t52 = in1[0] - in3[0];
  t83 = in1[1] - in3[1];
  t239 = in1[2] - in3[2];
  t81 = in1[3] - in3[3];
  t108 = in1[4] - in3[4];
  t58 = t16 + t17;
  t74 = (in1[3] * in1[3] + in1[4] * in1[4]) + 1.0;
  st.site = &q_emlrtRSI;
  st.site = &r_emlrtRSI;
  t69 = t52 * t52;
  st.site = &s_emlrtRSI;
  t70 = t83 * t83;
  st.site = &t_emlrtRSI;
  t71 = t239 * t239;
  st.site = &u_emlrtRSI;
  t72 = t81 * t81;
  st.site = &v_emlrtRSI;
  t73 = t108 * t108;
  t77 = (in1[1] * t2 + in1[2] * t3) + 1.0;
  st.site = &w_emlrtRSI;
  t79 = muDoubleScalarSqrt(t58);
  st.site = &x_emlrtRSI;
  t80 = 1.0 / (t74 * t74);
  st.site = &y_emlrtRSI;
  t81 = 1.0 / muDoubleScalarPower(t74, 3.0);
  t82 = in1[4] * t2 - in1[3] * t3;
  t83 = in1[0] * (t58 - 1.0);
  t85 = 1.0 / (t58 - 1.0);
  t86 = 1.0 / t79;
  st.site = &ab_emlrtRSI;
  if (-t16 + 1.0 < 0.0) {
    emlrtErrorWithMessageIdR2018a(
        &st, &d_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
        "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
  }
  t89 = muDoubleScalarSqrt(-t16 + 1.0);
  st.site = &bb_emlrtRSI;
  if (-t17 + 1.0 < 0.0) {
    emlrtErrorWithMessageIdR2018a(
        &st, &d_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
        "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
  }
  t90 = muDoubleScalarSqrt(-t17 + 1.0);
  t96 = 1.0 / t77;
  t108 = 1.0 / ((t2 * t7 + t3 * t8) + 2.0);
  t113 = t25 * t42 * muDoubleScalarAbs(t52);
  st.site = &cb_emlrtRSI;
  t87 = t85 * t85;
  t102 = in1[1] + t90;
  t103 = in1[2] + t89;
  t111 = 1.0 / (t79 + 1.0);
  st.site = &db_emlrtRSI;
  b_st.site = &i_emlrtRSI;
  if ((t113 < 0.0) && (!muDoubleScalarIsNaN(n1)) &&
      (muDoubleScalarFloor(n1) != n1)) {
    emlrtErrorWithMessageIdR2018a(&b_st, &e_emlrtRTEI,
                                  "Coder:toolbox:power_domainError",
                                  "Coder:toolbox:power_domainError", 0);
  }
  t117 = muDoubleScalarPower(t113, n1);
  st.site = &eb_emlrtRSI;
  st.site = &fb_emlrtRSI;
  b_st.site = &i_emlrtRSI;
  if ((t113 < 0.0) && (!muDoubleScalarIsNaN(n1 - 1.0)) &&
      (muDoubleScalarFloor(n1 - 1.0) != n1 - 1.0)) {
    emlrtErrorWithMessageIdR2018a(&b_st, &e_emlrtRTEI,
                                  "Coder:toolbox:power_domainError",
                                  "Coder:toolbox:power_domainError", 0);
  }
  t92 = muDoubleScalarSin(L1 - muDoubleScalarAtan(in1[2] * (1.0 / in1[1])));
  t105 = in1[1] + t2 * (t77 + 1.0);
  t106 = in1[2] + t3 * (t77 + 1.0);
  st.site = &gb_emlrtRSI;
  t109 = t102 * t102;
  st.site = &hb_emlrtRSI;
  t110 = t103 * t103;
  st.site = &ib_emlrtRSI;
  t112 = t111 * t111;
  st.site = &jb_emlrtRSI;
  if (-t83 < 0.0) {
    emlrtErrorWithMessageIdR2018a(
        &st, &d_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
        "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
  }
  t114 = muDoubleScalarSqrt(-t83);
  st.site = &kb_emlrtRSI;
  b_st.site = &i_emlrtRSI;
  if ((t117 + 1.0 < 0.0) && (!muDoubleScalarIsNaN(t43)) &&
      (muDoubleScalarFloor(t43) != t43)) {
    emlrtErrorWithMessageIdR2018a(&b_st, &e_emlrtRTEI,
                                  "Coder:toolbox:power_domainError",
                                  "Coder:toolbox:power_domainError", 0);
  }
  t124 = muDoubleScalarPower(t117 + 1.0, t43);
  st.site = &lb_emlrtRSI;
  b_st.site = &i_emlrtRSI;
  if ((t117 + 1.0 < 0.0) && (!muDoubleScalarIsNaN(t43 - 1.0)) &&
      (muDoubleScalarFloor(t43 - 1.0) != t43 - 1.0)) {
    emlrtErrorWithMessageIdR2018a(&b_st, &e_emlrtRTEI,
                                  "Coder:toolbox:power_domainError",
                                  "Coder:toolbox:power_domainError", 0);
  }
  t125 = muDoubleScalarExp(k1 * (in1[0] * t44 * (t79 - 1.0) + 1.0));
  t126 = W_p1 * t125;
  t231 = in4[3] * t24;
  t58 = t231 * t37;
  t175 = (-(t58 * (in1[3] * 2.0 - in3[3] * 2.0) * t80 * t85 * t109 * 4.0) +
          in4[3] * in1[3] * t24 * t37 * t72 * t81 * t85 * t109 * 16.0) +
         in1[3] * in4[4] * t24 * t37 * t73 * t81 * t85 * t110 * 16.0;
  t176_tmp_tmp = in4[4] * t24;
  t16 = t176_tmp_tmp * t37;
  t176 = (-(t16 * (in1[4] * 2.0 - in3[4] * 2.0) * t80 * t85 * t110 * 4.0) +
          in4[3] * in1[4] * t24 * t37 * t72 * t81 * t85 * t109 * 16.0) +
         in4[4] * in1[4] * t24 * t37 * t73 * t81 * t85 * t110 * 16.0;
  t17 = in4[1] * t24;
  t239 = in4[2] * t24;
  t83 = in4[0] * t24;
  Q = t83 * t40;
  t216 = t17 * t37;
  t215 = t58 * t72 * t80 * t85;
  t179_tmp = t239 * t37;
  t81 = t16 * t73 * t80 * t85;
  t179 = (((t216 * t70 * t85 / 4.0 + t179_tmp * t71 * t85 / 4.0) +
           t215 * t109 * 4.0) +
          t81 * t110 * 4.0) +
         Q * t69 * (t79 - 1.0) * t111 * t124 / 4.0;
  t180 = k1 * t44 * (t79 - 1.0) * t126 * t179;
  t58 = k1 * in1[0];
  t181 = t58 * in1[1] * t44 * t86 * t126 * t179;
  t182 = t58 * in1[2] * t44 * t86 * t126 * t179;
  t177 = t2 * t74 * t108 * t114 * (t126 + 1.0) * t175;
  t178 = t3 * t74 * t108 * t114 * (t126 + 1.0) * t176;
  t42 = (t126 + 1.0) *
        ((((((t17 * t38 * t70 * t85 / 4.0 + t239 * t38 * t71 * t85 / 4.0) +
             t231 * t38 * t72 * t80 * t85 * t109 * 4.0) +
            t176_tmp_tmp * t38 * t73 * t80 * t85 * t110 * 4.0) +
           t83 * (t38 * t38) * t69 * (t79 - 1.0) * t111 * t124 * 0.75) -
          Q * (t6 - in3[0] * 2.0) * (t79 - 1.0) * t111 * t124 / 4.0) -
         in4[0] * n1 * t24 * t25 * t40 * t42 * t43 * t69 * (t52 + t52) *
             (t79 - 1.0) * t111 * (1.0 / muDoubleScalarSqrt(t69)) *
             muDoubleScalarPower(t113, n1 - 1.0) *
             muDoubleScalarPower(t117 + 1.0, t43 - 1.0) / 8.0);
  t58 = in1[1] * in4[4] * t24 * t37 * t73 * t80;
  t16 = in4[0] * in1[1] * t24 * t40 * t69 * t86;
  t25 = (t126 + 1.0) *
        ((((((((-(t216 * (t7 - in3[1] * 2.0) * t85 / 4.0) +
                in4[1] * in1[1] * t24 * t37 * t70 * t87 / 2.0) +
               in1[1] * in4[2] * t24 * t37 * t71 * t87 / 2.0) +
              in1[1] * in4[3] * t24 * t37 * t72 * t80 * t87 * t109 * 8.0) +
             t58 * t87 * t110 * 8.0) -
            t215 * (t7 + t90 * 2.0) * 4.0) +
           t58 * t85 * (1.0 / t89) * t103 * 8.0) -
          t16 * t111 * t124 / 4.0) +
         t16 * (t79 - 1.0) * t112 * t124 / 4.0);
  t58 = in1[2] * in4[3] * t24 * t37 * t72 * t80;
  t83 = in4[0] * in1[2] * t24 * t40 * t69 * t86;
  t102 = (t126 + 1.0) *
         ((((((((-(t179_tmp * (t8 - in3[2] * 2.0) * t85 / 4.0) +
                 in4[1] * in1[2] * t24 * t37 * t70 * t87 / 2.0) +
                in4[2] * in1[2] * t24 * t37 * t71 * t87 / 2.0) +
               t58 * t87 * t109 * 8.0) +
              in1[2] * in4[4] * t24 * t37 * t73 * t80 * t87 * t110 * 8.0) -
             t81 * (t8 + t89 * 2.0) * 4.0) +
            t58 * t85 * (1.0 / t90) * t102 * 8.0) -
           t83 * t111 * t124 / 4.0) +
          t83 * (t79 - 1.0) * t112 * t124 / 4.0);
  t111 = t180 - t42;
  t17 = in1[0] * t77 * t85 * t114 * t111 * 2.0;
  t109 = t182 - t102;
  t110 = -t2 * t114;
  t112 = t181 - t25;
  t103 = -t3 * t114 * t112;
  t16 = in1[0] * t79 * t85 * t92 * t114 * t111;
  t113 = t16 * 0.0;
  t117 = t16 * -2.0;
  t74 = -t96 * t105 * t114 * t112;
  t90 = -t96 * t106 * t114 * t109;
  t87 = t96 * t105 * t114;
  t231 = t87 * t112;
  t179_tmp = t96 * t106 * t114;
  t108 = t179_tmp * t109;
  t215 = t117;
  t216 = t113;
  t44 = in1[1] * t82 * t96 * t114;
  t124 = ((t177 + t178) + -in1[2] * t82 * t96 * t114 * t112) + t44 * t109;
  st.site = &mb_emlrtRSI;
  t81 = t2 * t114 * t109;
  t83 = t103 + t81;
  t176_tmp_tmp = t6 * t79 * t85 * t92 * t114;
  t58 = t83 + t176_tmp_tmp * t111;
  D[0] = t83 + t16 * 2.0;
  D[1] = (t17 + t74) + t90;
  D[2] = t124;
  Q = -(t126 + 1.0) * t179;
  st.site = &nb_emlrtRSI;
  t83 = (-t17 + t231) + t108;
  st.site = &ob_emlrtRSI;
  t239 = t58 * t58 + t83 * t83;
  t16 = t3 * t114;
  t83 = t16 * t112;
  t113 = ((((-t6 * t77 * t85 * t114 * t111 + t113) + t81 * 0.0) + t83 * 0.0) +
          t231) +
         t108;
  t117 = (t117 - t81) + t83;
  st.site = &pb_emlrtRSI;
  t58 = 1.0 / muDoubleScalarSqrt(t239);
  t108 = 1.0 / muDoubleScalarHypot(t113, t117);
  st.site = &qb_emlrtRSI;
  t81 = 1.0 / muDoubleScalarSqrt(t124 * t124 * (1.0 / t239) + 1.0);
  t239 = (t110 * t109 - t103) + t215;
  t83 = ((-in1[0] * t77 * t85 * t114 * t111 * 2.0 + t216) - t74) - t90;
  *Qdot = (((-t112 * ((t16 * t108 * t81 * t239 + t87 * t108 * t81 * t83) -
                      in1[2] * t82 * t96 * t114 * t124 * t58 * t81) -
             t109 * ((t110 * t108 * t81 * t239 + t179_tmp * t108 * t81 * t83) +
                     t44 * t124 * t58 * t81)) +
            (t6 * t77 * t85 * t114 * t108 * t81 * t83 +
             t176_tmp_tmp * t108 * t81 * t239) *
                t111) -
           t177 * t124 * t58 * t81) -
          t178 * t124 * t58 * t81;
  *P = t125;
  partial_Q_partial_oe[0] = -t180 + t42;
  partial_Q_partial_oe[1] = -t181 + t25;
  partial_Q_partial_oe[2] = -t182 + t102;
  partial_Q_partial_oe[3] = (t126 + 1.0) * t175;
  partial_Q_partial_oe[4] = (t126 + 1.0) * t176;
  return Q;
}

/* End of code generation (D_Q_Qdot_P_partial_Q_partial_oe_func.c) */
