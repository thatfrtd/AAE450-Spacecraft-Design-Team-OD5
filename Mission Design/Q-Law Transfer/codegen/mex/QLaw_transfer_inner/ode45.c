/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * ode45.c
 *
 * Code generation for function 'ode45'
 *
 */

/* Include files */
#include "ode45.h"
#include "QLaw_transfer_inner.h"
#include "QLaw_transfer_inner_data.h"
#include "QLaw_transfer_inner_emxutil.h"
#include "QLaw_transfer_inner_mexutil.h"
#include "QLaw_transfer_inner_types.h"
#include "odearguments.h"
#include "rt_nonfinite.h"
#include "warning.h"
#include "mwmathutil.h"
#include <emmintrin.h>
#include <math.h>
#include <string.h>

/* Variable Definitions */
static emlrtRSInfo uc_emlrtRSI = {
    17,      /* lineNo */
    "ode45", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2025b\\toolbox\\eml\\lib\\matlab\\funfun\\ode45.m" /* pathName
                                                                         */
};

static emlrtRSInfo vc_emlrtRSI = {
    563,                  /* lineNo */
    "explicitRungeKutta", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2025b\\toolbox\\eml\\lib\\matlab\\funfun\\private\\explici"
    "tRungeKutta.m" /* pathName */
};

static emlrtRSInfo wc_emlrtRSI = {
    358,                  /* lineNo */
    "explicitRungeKutta", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2025b\\toolbox\\eml\\lib\\matlab\\funfun\\private\\explici"
    "tRungeKutta.m" /* pathName */
};

static emlrtRSInfo xc_emlrtRSI = {
    345,                  /* lineNo */
    "explicitRungeKutta", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2025b\\toolbox\\eml\\lib\\matlab\\funfun\\private\\explici"
    "tRungeKutta.m" /* pathName */
};

static emlrtRSInfo yc_emlrtRSI = {
    348,                  /* lineNo */
    "explicitRungeKutta", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2025b\\toolbox\\eml\\lib\\matlab\\funfun\\private\\explici"
    "tRungeKutta.m" /* pathName */
};

static emlrtRSInfo ad_emlrtRSI = {
    347,                  /* lineNo */
    "explicitRungeKutta", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2025b\\toolbox\\eml\\lib\\matlab\\funfun\\private\\explici"
    "tRungeKutta.m" /* pathName */
};

static emlrtRSInfo bd_emlrtRSI = {
    305,                  /* lineNo */
    "explicitRungeKutta", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2025b\\toolbox\\eml\\lib\\matlab\\funfun\\private\\explici"
    "tRungeKutta.m" /* pathName */
};

static emlrtRSInfo cd_emlrtRSI = {
    298,                  /* lineNo */
    "explicitRungeKutta", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2025b\\toolbox\\eml\\lib\\matlab\\funfun\\private\\explici"
    "tRungeKutta.m" /* pathName */
};

static emlrtRSInfo dd_emlrtRSI = {
    73,                   /* lineNo */
    "explicitRungeKutta", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2025b\\toolbox\\eml\\lib\\matlab\\funfun\\private\\explici"
    "tRungeKutta.m" /* pathName */
};

static emlrtRSInfo
    md_emlrtRSI =
        {
            15,        /* lineNo */
            "num2str", /* fcnName */
            "C:\\Program "
            "Files\\MATLAB\\R2025b\\toolbox\\eml\\eml\\+coder\\+"
            "internal\\num2str.m" /* pathName */
};

static emlrtMCInfo c_emlrtMCI = {
    53,        /* lineNo */
    19,        /* colNo */
    "flt2str", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2025b\\toolbox\\shared\\coder\\coder\\lib\\+coder\\+"
    "internal\\flt2str.m" /* pName */
};

static emlrtRTEInfo r_emlrtRTEI = {
    185,                  /* lineNo */
    5,                    /* colNo */
    "explicitRungeKutta", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2025b\\toolbox\\eml\\lib\\matlab\\funfun\\private\\explici"
    "tRungeKutta.m" /* pName */
};

static emlrtRTEInfo s_emlrtRTEI = {
    186,                  /* lineNo */
    5,                    /* colNo */
    "explicitRungeKutta", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2025b\\toolbox\\eml\\lib\\matlab\\funfun\\private\\explici"
    "tRungeKutta.m" /* pName */
};

static emlrtRTEInfo t_emlrtRTEI = {
    588,                  /* lineNo */
    1,                    /* colNo */
    "explicitRungeKutta", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2025b\\toolbox\\eml\\lib\\matlab\\funfun\\private\\explici"
    "tRungeKutta.m" /* pName */
};

static emlrtRTEInfo u_emlrtRTEI = {
    589,                  /* lineNo */
    1,                    /* colNo */
    "explicitRungeKutta", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2025b\\toolbox\\eml\\lib\\matlab\\funfun\\private\\explici"
    "tRungeKutta.m" /* pName */
};

static emlrtRTEInfo v_emlrtRTEI = {
    14,                  /* lineNo */
    9,                   /* colNo */
    "appendZeroColumns", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2025b\\toolbox\\eml\\lib\\matlab\\funfun\\private\\appendZ"
    "eroColumns.m" /* pName */
};

static emlrtRSInfo od_emlrtRSI = {
    53,        /* lineNo */
    "flt2str", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2025b\\toolbox\\shared\\coder\\coder\\lib\\+coder\\+"
    "internal\\flt2str.m" /* pathName */
};

/* Function Declarations */
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               char_T y[23]);

static const mxArray *b_sprintf(const emlrtStack *sp, const mxArray *m,
                                const mxArray *m1, emlrtMCInfo *location);

static int32_T div_nde_s32_floor(int32_T numerator);

static void emlrt_marshallIn(const emlrtStack *sp,
                             const mxArray *a__output_of_sprintf_,
                             const char_T *identifier, char_T y[23]);

static void k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                               const emlrtMsgIdentifier *msgId, char_T ret[23]);

/* Function Definitions */
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId, char_T y[23])
{
  k_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static const mxArray *b_sprintf(const emlrtStack *sp, const mxArray *m,
                                const mxArray *m1, emlrtMCInfo *location)
{
  const mxArray *pArrays[2];
  const mxArray *m2;
  pArrays[0] = m;
  pArrays[1] = m1;
  return emlrtCallMATLABR2012b((emlrtConstCTX)sp, 1, &m2, 2, &pArrays[0],
                               "sprintf", true, location);
}

static int32_T div_nde_s32_floor(int32_T numerator)
{
  int32_T quotient;
  if ((numerator < 0) && (numerator % 7 != 0)) {
    quotient = -1;
  } else {
    quotient = 0;
  }
  quotient += numerator / 7;
  return quotient;
}

static void emlrt_marshallIn(const emlrtStack *sp,
                             const mxArray *a__output_of_sprintf_,
                             const char_T *identifier, char_T y[23])
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  b_emlrt_marshallIn(sp, emlrtAlias(a__output_of_sprintf_), &thisId, y);
  emlrtDestroyArray(&a__output_of_sprintf_);
}

static void k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                               const emlrtMsgIdentifier *msgId, char_T ret[23])
{
  static const int32_T dims[2] = {1, 23};
  emlrtCheckBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "char", false, 2U,
                          (const void *)&dims[0]);
  emlrtImportCharArrayR2015b((emlrtConstCTX)sp, src, &ret[0], 23);
  emlrtDestroyArray(&src);
}

void ode45(const emlrtStack *sp,
           const real_T c_ode_workspace_a_control_works[3],
           real_T ode_workspace_mdot, const real_T tspan[2],
           const real_T b_y0[7], real_T options_AbsTol,
           const real_T options_InitialStep_data[],
           const real_T options_MaxStep_data[], real_T options_RelTol,
           emxArray_real_T *varargout_1, emxArray_real_T *varargout_2)
{
  static const real_T x[21] = {0.2,
                               0.075,
                               0.225,
                               0.97777777777777775,
                               -3.7333333333333334,
                               3.5555555555555554,
                               2.9525986892242035,
                               -11.595793324188385,
                               9.8228928516994358,
                               -0.29080932784636487,
                               2.8462752525252526,
                               -10.757575757575758,
                               8.9064227177434727,
                               0.27840909090909088,
                               -0.2735313036020583,
                               0.091145833333333329,
                               0.0,
                               0.44923629829290207,
                               0.65104166666666663,
                               -0.322376179245283,
                               0.13095238095238096};
  static const real_T b[7] = {0.0012326388888888888,
                              0.0,
                              -0.0042527702905061394,
                              0.036979166666666667,
                              -0.05086379716981132,
                              0.0419047619047619,
                              -0.025};
  static const real_T b_b[7] = {-2.859375,
                                0.0,
                                4.0431266846361185,
                                -3.90625,
                                2.7939268867924527,
                                -1.5714285714285714,
                                1.5};
  static const real_T c_b[7] = {3.0833333333333335,
                                0.0,
                                -6.2893081761006293,
                                10.416666666666666,
                                -6.8773584905660377,
                                3.6666666666666665,
                                -4.0};
  static const real_T d_b[7] = {-1.1328125,
                                0.0,
                                2.6954177897574123,
                                -5.859375,
                                3.7610554245283021,
                                -1.9642857142857142,
                                2.5};
  static const int32_T iv[2] = {1, 7};
  static const int32_T iv1[2] = {1, 7};
  static const char_T rfmt[7] = {'%', '2', '3', '.', '1', '5', 'e'};
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack f_st;
  emlrtStack st;
  emxArray_real_T *tout;
  emxArray_real_T *yout;
  const mxArray *b_y;
  const mxArray *c_y;
  const mxArray *m;
  real_T f[49];
  real_T youtnew[28];
  real_T fhBI2[7];
  real_T fhBI3[7];
  real_T fhBI4[7];
  real_T y[7];
  real_T ystage[7];
  real_T toutnew[4];
  real_T tref[3];
  real_T dv[2];
  real_T absh;
  real_T absx;
  real_T d1;
  real_T hmax;
  real_T rtol;
  real_T t;
  real_T tdir;
  real_T tfinal;
  real_T threshold;
  real_T *tout_data;
  real_T *varargout_1_data;
  real_T *yout_data;
  int32_T exponent;
  int32_T i;
  int32_T ia;
  int32_T j;
  int32_T nout;
  int32_T tout_tmp;
  int32_T ynew_tmp;
  boolean_T Done;
  boolean_T MinStepExit;
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
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  st.site = &uc_emlrtRSI;
  for (i = 0; i < 7; i++) {
    y[i] = b_y0[i];
  }
  real_T b_tspan[2];
  int32_T htry_size[2];
  b_st.site = &dd_emlrtRSI;
  t = odearguments(&b_st, c_ode_workspace_a_control_works, ode_workspace_mdot,
                   tspan, y, options_AbsTol, options_InitialStep_data,
                   options_MaxStep_data, options_RelTol, b_tspan, &tfinal,
                   &tdir, ystage, &threshold, &rtol, &hmax, (real_T *)&d1,
                   htry_size, &absx);
  emxInit_real_T(&st, &tout, 2, &r_emlrtRTEI);
  tout_tmp = tout->size[0] * tout->size[1];
  tout->size[0] = 1;
  tout->size[1] = 200;
  emxEnsureCapacity_real_T(&st, tout, tout_tmp, &r_emlrtRTEI);
  tout_data = tout->data;
  for (i = 0; i < 200; i++) {
    tout_data[i] = 0.0;
  }
  emxInit_real_T(&st, &yout, 2, &s_emlrtRTEI);
  tout_tmp = yout->size[0] * yout->size[1];
  yout->size[0] = 7;
  yout->size[1] = 200;
  emxEnsureCapacity_real_T(&st, yout, tout_tmp, &s_emlrtRTEI);
  yout_data = yout->data;
  for (i = 0; i < 1400; i++) {
    yout_data[i] = 0.0;
  }
  nout = 0;
  tout_data[0] = t;
  for (i = 0; i < 7; i++) {
    yout_data[i] = y[i];
  }
  absx = muDoubleScalarAbs(t);
  if (muDoubleScalarIsInf(absx) || muDoubleScalarIsNaN(absx)) {
    absx = rtNaN;
  } else if (absx < 4.4501477170144028E-308) {
    absx = 4.94065645841247E-324;
  } else {
    frexp(absx, &ynew_tmp);
    absx = ldexp(1.0, ynew_tmp - 53);
  }
  absh = muDoubleScalarMin(
      hmax, muDoubleScalarMax(muDoubleScalarMax(16.0 * absx, 0.0), d1));
  memset(&f[0], 0, 49U * sizeof(real_T));
  for (i = 0; i < 7; i++) {
    f[i] = ystage[i];
  }
  MinStepExit = false;
  Done = false;
  int32_T exitg1;
  do {
    __m128d r;
    __m128d r1;
    __m128d r2;
    real_T ynew[7];
    real_T d2;
    real_T err;
    real_T h;
    real_T hmin;
    real_T mxerr;
    real_T tnew;
    int32_T Bcolidx;
    int32_T ncols;
    boolean_T NoFailedAttempts;
    exitg1 = 0;
    absx = muDoubleScalarAbs(t);
    if (muDoubleScalarIsInf(absx) || muDoubleScalarIsNaN(absx)) {
      absx = rtNaN;
    } else if (absx < 4.4501477170144028E-308) {
      absx = 4.94065645841247E-324;
    } else {
      frexp(absx, &exponent);
      absx = ldexp(1.0, exponent - 53);
    }
    hmin = muDoubleScalarMax(16.0 * absx, 0.0);
    absh = muDoubleScalarMin(hmax, muDoubleScalarMax(hmin, absh));
    h = tdir * absh;
    absx = tfinal - t;
    d1 = muDoubleScalarAbs(absx);
    if (1.1 * absh >= d1) {
      h = absx;
      absh = d1;
      Done = true;
    }
    NoFailedAttempts = true;
    int32_T exitg2;
    do {
      exitg2 = 0;
      Bcolidx = 0;
      for (j = 0; j < 5; j++) {
        Bcolidx += j;
        for (i = 0; i < 7; i++) {
          ystage[i] = y[i];
        }
        if (!(h == 0.0)) {
          ynew_tmp = 7 * j + 1;
          for (i = 1; i <= ynew_tmp; i += 7) {
            absx = h * x[Bcolidx + div_nde_s32_floor(i - 1)];
            tout_tmp = i + 6;
            for (ia = i; ia <= tout_tmp; ia++) {
              ncols = ia - i;
              ystage[ncols] += f[ia - 1] * absx;
            }
          }
        }
        b_st.site = &cd_emlrtRSI;
        c_st.site = &wb_emlrtRSI;
        d_st.site = &gd_emlrtRSI;
        e_st.site = &hd_emlrtRSI;
        f_st.site = &wb_emlrtRSI;
        QLaw_transfer_inner_anonFcn2(&f_st, c_ode_workspace_a_control_works,
                                     ode_workspace_mdot, ystage,
                                     &f[7 * (j + 1)]);
      }
      tnew = t + h;
      for (i = 0; i < 7; i++) {
        ynew[i] = y[i];
      }
      if (!(h == 0.0)) {
        for (i = 0; i <= 35; i += 7) {
          absx = h * x[(Bcolidx + div_nde_s32_floor(i)) + 5];
          tout_tmp = i + 7;
          for (ia = i + 1; ia <= tout_tmp; ia++) {
            ynew_tmp = (ia - i) - 1;
            ynew[ynew_tmp] += f[ia - 1] * absx;
          }
        }
      }
      b_st.site = &bd_emlrtRSI;
      c_st.site = &wb_emlrtRSI;
      d_st.site = &gd_emlrtRSI;
      e_st.site = &hd_emlrtRSI;
      f_st.site = &wb_emlrtRSI;
      QLaw_transfer_inner_anonFcn2(&f_st, c_ode_workspace_a_control_works,
                                   ode_workspace_mdot, ynew, &f[42]);
      memset(&ystage[0], 0, 7U * sizeof(real_T));
      for (i = 0; i < 7; i++) {
        r = _mm_loadu_pd(&f[7 * i]);
        r1 = _mm_loadu_pd(&ystage[0]);
        absx = b[i];
        r2 = _mm_set1_pd(absx);
        _mm_storeu_pd(&ystage[0], _mm_add_pd(r1, _mm_mul_pd(r, r2)));
        r = _mm_loadu_pd(&f[7 * i + 2]);
        r1 = _mm_loadu_pd(&ystage[2]);
        _mm_storeu_pd(&ystage[2], _mm_add_pd(r1, _mm_mul_pd(r, r2)));
        r = _mm_loadu_pd(&f[7 * i + 4]);
        r1 = _mm_loadu_pd(&ystage[4]);
        _mm_storeu_pd(&ystage[4], _mm_add_pd(r1, _mm_mul_pd(r, r2)));
        ystage[6] += f[7 * i + 6] * absx;
      }
      if (Done) {
        tnew = tfinal;
      }
      h = tnew - t;
      mxerr = 0.0;
      for (i = 0; i < 7; i++) {
        absx = muDoubleScalarAbs(ystage[i]);
        d1 = muDoubleScalarAbs(y[i]);
        d2 = muDoubleScalarAbs(ynew[i]);
        if ((d1 > d2) || muDoubleScalarIsNaN(d2)) {
          if (d1 > threshold) {
            absx /= d1;
          } else {
            absx /= threshold;
          }
        } else if (d2 > threshold) {
          absx /= d2;
        } else {
          absx /= threshold;
        }
        if ((absx > mxerr) || muDoubleScalarIsNaN(absx)) {
          mxerr = absx;
        }
      }
      err = absh * mxerr;
      if (!(err <= rtol)) {
        if (absh <= hmin) {
          char_T b_str[23];
          char_T str[23];
          b_st.site = &ad_emlrtRSI;
          c_st.site = &md_emlrtRSI;
          b_y = NULL;
          m = emlrtCreateCharArray(2, &iv[0]);
          emlrtInitCharArrayR2013a(&c_st, 7, m, &rfmt[0]);
          emlrtAssign(&b_y, m);
          d_st.site = &od_emlrtRSI;
          emlrt_marshallIn(
              &d_st, b_sprintf(&d_st, b_y, emlrt_marshallOut(t), &c_emlrtMCI),
              "<output of sprintf>", str);
          b_st.site = &yc_emlrtRSI;
          c_st.site = &md_emlrtRSI;
          c_y = NULL;
          m = emlrtCreateCharArray(2, &iv1[0]);
          emlrtInitCharArrayR2013a(&c_st, 7, m, &rfmt[0]);
          emlrtAssign(&c_y, m);
          d_st.site = &od_emlrtRSI;
          emlrt_marshallIn(
              &d_st,
              b_sprintf(&d_st, c_y, emlrt_marshallOut(hmin), &c_emlrtMCI),
              "<output of sprintf>", b_str);
          b_st.site = &xc_emlrtRSI;
          b_warning(&b_st, str, b_str);
          MinStepExit = true;
          exitg2 = 1;
        } else {
          if (NoFailedAttempts) {
            NoFailedAttempts = false;
            b_st.site = &wc_emlrtRSI;
            absx = rtol / err;
            c_st.site = &h_emlrtRSI;
            d_st.site = &i_emlrtRSI;
            if (absx < 0.0) {
              emlrtErrorWithMessageIdR2018a(
                  &d_st, &e_emlrtRTEI, "Coder:toolbox:power_domainError",
                  "Coder:toolbox:power_domainError", 0);
            }
            absh = muDoubleScalarMax(
                hmin, absh * muDoubleScalarMax(
                                 0.1, 0.8 * muDoubleScalarPower(absx, 0.2)));
          } else {
            absh = muDoubleScalarMax(hmin, 0.5 * absh);
          }
          h = tdir * absh;
          Done = false;
        }
      } else {
        exitg2 = 1;
      }
    } while (exitg2 == 0);
    if (MinStepExit) {
      exitg1 = 1;
    } else {
      __m128d r3;
      __m128d r4;
      Bcolidx = nout + 1;
      dv[0] = 0.0;
      dv[1] = 1.0;
      r = _mm_loadu_pd(&dv[0]);
      r1 = _mm_set1_pd(0.25);
      r = _mm_add_pd(
          _mm_set1_pd(t),
          _mm_mul_pd(_mm_set1_pd(h), _mm_add_pd(r1, _mm_mul_pd(r1, r))));
      _mm_storeu_pd(&tref[0], r);
      _mm_storeu_pd(&toutnew[0], r);
      absx = t + h * 0.75;
      tref[2] = absx;
      toutnew[2] = absx;
      toutnew[3] = tnew;
      memset(&fhBI2[0], 0, 7U * sizeof(real_T));
      memset(&fhBI3[0], 0, 7U * sizeof(real_T));
      memset(&fhBI4[0], 0, 7U * sizeof(real_T));
      for (i = 0; i < 7; i++) {
        ystage[i] = f[i] * h;
        absx = h * b_b[i];
        d1 = h * c_b[i];
        d2 = h * d_b[i];
        r1 = _mm_loadu_pd(&f[7 * i]);
        r = _mm_loadu_pd(&fhBI2[0]);
        r2 = _mm_set1_pd(absx);
        _mm_storeu_pd(&fhBI2[0], _mm_add_pd(r, _mm_mul_pd(r1, r2)));
        r = _mm_loadu_pd(&fhBI3[0]);
        r3 = _mm_set1_pd(d1);
        _mm_storeu_pd(&fhBI3[0], _mm_add_pd(r, _mm_mul_pd(r1, r3)));
        r = _mm_loadu_pd(&fhBI4[0]);
        r4 = _mm_set1_pd(d2);
        _mm_storeu_pd(&fhBI4[0], _mm_add_pd(r, _mm_mul_pd(r1, r4)));
        r1 = _mm_loadu_pd(&f[7 * i + 2]);
        r = _mm_loadu_pd(&fhBI2[2]);
        _mm_storeu_pd(&fhBI2[2], _mm_add_pd(r, _mm_mul_pd(r1, r2)));
        r = _mm_loadu_pd(&fhBI3[2]);
        _mm_storeu_pd(&fhBI3[2], _mm_add_pd(r, _mm_mul_pd(r1, r3)));
        r = _mm_loadu_pd(&fhBI4[2]);
        _mm_storeu_pd(&fhBI4[2], _mm_add_pd(r, _mm_mul_pd(r1, r4)));
        r1 = _mm_loadu_pd(&f[7 * i + 4]);
        r = _mm_loadu_pd(&fhBI2[4]);
        _mm_storeu_pd(&fhBI2[4], _mm_add_pd(r, _mm_mul_pd(r1, r2)));
        r = _mm_loadu_pd(&fhBI3[4]);
        _mm_storeu_pd(&fhBI3[4], _mm_add_pd(r, _mm_mul_pd(r1, r3)));
        r = _mm_loadu_pd(&fhBI4[4]);
        _mm_storeu_pd(&fhBI4[4], _mm_add_pd(r, _mm_mul_pd(r1, r4)));
        mxerr = f[7 * i + 6];
        fhBI2[6] += mxerr * absx;
        fhBI3[6] += mxerr * d1;
        fhBI4[6] += mxerr * d2;
      }
      for (i = 0; i < 3; i++) {
        __m128d r5;
        absx = (tref[i] - t) / h;
        r = _mm_loadu_pd(&fhBI4[0]);
        r1 = _mm_loadu_pd(&fhBI3[0]);
        r2 = _mm_loadu_pd(&fhBI2[0]);
        r3 = _mm_loadu_pd(&ystage[0]);
        r4 = _mm_loadu_pd(&y[0]);
        r5 = _mm_set1_pd(absx);
        _mm_storeu_pd(
            &youtnew[7 * i],
            _mm_add_pd(
                _mm_mul_pd(
                    _mm_add_pd(
                        _mm_mul_pd(
                            _mm_add_pd(
                                _mm_mul_pd(_mm_add_pd(_mm_mul_pd(r, r5), r1),
                                           r5),
                                r2),
                            r5),
                        r3),
                    r5),
                r4));
        r = _mm_loadu_pd(&fhBI4[2]);
        r1 = _mm_loadu_pd(&fhBI3[2]);
        r2 = _mm_loadu_pd(&fhBI2[2]);
        r3 = _mm_loadu_pd(&ystage[2]);
        r4 = _mm_loadu_pd(&y[2]);
        _mm_storeu_pd(
            &youtnew[7 * i + 2],
            _mm_add_pd(
                _mm_mul_pd(
                    _mm_add_pd(
                        _mm_mul_pd(
                            _mm_add_pd(
                                _mm_mul_pd(_mm_add_pd(_mm_mul_pd(r, r5), r1),
                                           r5),
                                r2),
                            r5),
                        r3),
                    r5),
                r4));
        r = _mm_loadu_pd(&fhBI4[4]);
        r1 = _mm_loadu_pd(&fhBI3[4]);
        r2 = _mm_loadu_pd(&fhBI2[4]);
        r3 = _mm_loadu_pd(&ystage[4]);
        r4 = _mm_loadu_pd(&y[4]);
        _mm_storeu_pd(
            &youtnew[7 * i + 4],
            _mm_add_pd(
                _mm_mul_pd(
                    _mm_add_pd(
                        _mm_mul_pd(
                            _mm_add_pd(
                                _mm_mul_pd(_mm_add_pd(_mm_mul_pd(r, r5), r1),
                                           r5),
                                r2),
                            r5),
                        r3),
                    r5),
                r4));
        youtnew[7 * i + 6] =
            (((fhBI4[6] * absx + fhBI3[6]) * absx + fhBI2[6]) * absx +
             ystage[6]) *
                absx +
            y[6];
      }
      for (i = 0; i < 7; i++) {
        youtnew[i + 21] = ynew[i];
      }
      nout += 4;
      if (nout + 1 > tout->size[1]) {
        ncols = tout->size[1];
        tout_tmp = tout->size[0] * tout->size[1];
        tout->size[0] = 1;
        ynew_tmp = tout->size[1] + 200;
        tout->size[1] += 200;
        emxEnsureCapacity_real_T(&st, tout, tout_tmp, &v_emlrtRTEI);
        tout_data = tout->data;
        tout_tmp = yout->size[0] * yout->size[1];
        yout->size[0] = 7;
        yout->size[1] = ynew_tmp;
        emxEnsureCapacity_real_T(&st, yout, tout_tmp, &v_emlrtRTEI);
        yout_data = yout->data;
        for (i = 0; i < 200; i++) {
          tout_tmp = ncols + i;
          tout_data[tout_tmp] = 0.0;
          for (ia = 0; ia < 7; ia++) {
            yout_data[ia + 7 * tout_tmp] = 0.0;
          }
        }
      }
      for (i = 0; i < 4; i++) {
        tout_tmp = i + Bcolidx;
        tout_data[tout_tmp] = toutnew[i];
        for (ia = 0; ia < 7; ia++) {
          yout_data[ia + 7 * tout_tmp] = youtnew[ia + 7 * i];
        }
      }
      if (Done) {
        exitg1 = 1;
      } else {
        if (NoFailedAttempts) {
          b_st.site = &vc_emlrtRSI;
          absx = err / rtol;
          c_st.site = &h_emlrtRSI;
          d_st.site = &i_emlrtRSI;
          if (absx < 0.0) {
            emlrtErrorWithMessageIdR2018a(&d_st, &e_emlrtRTEI,
                                          "Coder:toolbox:power_domainError",
                                          "Coder:toolbox:power_domainError", 0);
          }
          absx = 1.25 * muDoubleScalarPower(absx, 0.2);
          if (absx > 0.2) {
            absh /= absx;
          } else {
            absh *= 5.0;
          }
        }
        t = tnew;
        for (i = 0; i < 7; i++) {
          y[i] = ynew[i];
          f[i] = f[i + 42];
        }
      }
    }
  } while (exitg1 == 0);
  if (nout + 1 < 1) {
    nout = -1;
  }
  tout_tmp = varargout_1->size[0];
  varargout_1->size[0] = nout + 1;
  emxEnsureCapacity_real_T(&st, varargout_1, tout_tmp, &t_emlrtRTEI);
  varargout_1_data = varargout_1->data;
  for (i = 0; i <= nout; i++) {
    varargout_1_data[i] = tout_data[i];
  }
  emxFree_real_T(&st, &tout);
  tout_tmp = varargout_2->size[0] * varargout_2->size[1];
  varargout_2->size[0] = nout + 1;
  varargout_2->size[1] = 7;
  emxEnsureCapacity_real_T(&st, varargout_2, tout_tmp, &u_emlrtRTEI);
  tout_data = varargout_2->data;
  for (i = 0; i < 7; i++) {
    for (ia = 0; ia <= nout; ia++) {
      tout_data[ia + varargout_2->size[0] * i] = yout_data[i + 7 * ia];
    }
  }
  emxFree_real_T(&st, &yout);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

/* End of code generation (ode45.c) */
