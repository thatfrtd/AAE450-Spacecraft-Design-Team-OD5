/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * linspace.c
 *
 * Code generation for function 'linspace'
 *
 */

/* Include files */
#include "linspace.h"
#include "QLaw_transfer_inner_emxutil.h"
#include "QLaw_transfer_inner_types.h"
#include "rt_nonfinite.h"
#include "mwmathutil.h"
#include <emmintrin.h>

/* Variable Definitions */
static emlrtRTEInfo f_emlrtRTEI = {
    31,         /* lineNo */
    33,         /* colNo */
    "linspace", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2025b\\toolbox\\eml\\lib\\matlab\\elmat\\linspace.m" /* pName
                                                                           */
};

static emlrtRTEInfo q_emlrtRTEI = {
    45,         /* lineNo */
    20,         /* colNo */
    "linspace", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2025b\\toolbox\\eml\\lib\\matlab\\elmat\\linspace.m" /* pName
                                                                           */
};

/* Function Definitions */
void linspace(const emlrtStack *sp, real_T N, emxArray_real_T *y)
{
  real_T dv[2];
  real_T *y_data;
  int32_T k;
  if (muDoubleScalarIsNaN(N)) {
    emlrtErrorWithMessageIdR2018a(sp, &f_emlrtRTEI,
                                  "Coder:toolbox:MustNotBeNaN",
                                  "Coder:toolbox:MustNotBeNaN", 3, 4, 1, "N");
  }
  if (N < 1.0) {
    y->size[0] = 1;
    y->size[1] = 0;
  } else {
    int32_T i;
    int32_T scalarLB;
    scalarLB = y->size[0] * y->size[1];
    y->size[0] = 1;
    i = (int32_T)muDoubleScalarFloor(N);
    y->size[1] = i;
    emxEnsureCapacity_real_T(sp, y, scalarLB, &q_emlrtRTEI);
    y_data = y->data;
    y_data[i - 1] = 6.2831853071795862;
    if (y->size[1] >= 2) {
      y_data[0] = 0.0;
      if (y->size[1] >= 3) {
        real_T delta1;
        int32_T vectorUB;
        delta1 = 6.2831853071795862 / ((real_T)y->size[1] - 1.0);
        scalarLB = ((y->size[1] - 2) / 2) << 1;
        vectorUB = scalarLB - 2;
        for (k = 0; k <= vectorUB; k += 2) {
          __m128d r;
          dv[0] = k + 1;
          dv[1] = k + 2;
          r = _mm_loadu_pd(&dv[0]);
          _mm_storeu_pd(&y_data[k + 1], _mm_mul_pd(r, _mm_set1_pd(delta1)));
        }
        for (k = scalarLB; k <= i - 3; k++) {
          y_data[k + 1] = (real_T)(k + 1) * delta1;
        }
      }
    }
  }
}

/* End of code generation (linspace.c) */
