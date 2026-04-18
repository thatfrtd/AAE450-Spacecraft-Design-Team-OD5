/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * warning.c
 *
 * Code generation for function 'warning'
 *
 */

/* Include files */
#include "warning.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtMCInfo emlrtMCI = {
    14,        /* lineNo */
    25,        /* colNo */
    "warning", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2025b\\toolbox\\shared\\coder\\coder\\lib\\+coder\\+"
    "internal\\warning.m" /* pName */
};

static emlrtMCInfo b_emlrtMCI = {
    14,        /* lineNo */
    9,         /* colNo */
    "warning", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2025b\\toolbox\\shared\\coder\\coder\\lib\\+coder\\+"
    "internal\\warning.m" /* pName */
};

static emlrtRSInfo nd_emlrtRSI = {
    14,        /* lineNo */
    "warning", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2025b\\toolbox\\shared\\coder\\coder\\lib\\+coder\\+"
    "internal\\warning.m" /* pathName */
};

/* Function Declarations */
static void b_feval(const emlrtStack *sp, const mxArray *m, const mxArray *m1,
                    emlrtMCInfo *location);

static const mxArray *c_feval(const emlrtStack *sp, const mxArray *m,
                              const mxArray *m1, const mxArray *m2,
                              const mxArray *m3, emlrtMCInfo *location);

static const mxArray *feval(const emlrtStack *sp, const mxArray *m,
                            const mxArray *m1, const mxArray *m2,
                            emlrtMCInfo *location);

/* Function Definitions */
static void b_feval(const emlrtStack *sp, const mxArray *m, const mxArray *m1,
                    emlrtMCInfo *location)
{
  const mxArray *pArrays[2];
  pArrays[0] = m;
  pArrays[1] = m1;
  emlrtCallMATLABR2012b((emlrtConstCTX)sp, 0, NULL, 2, &pArrays[0], "feval",
                        true, location);
}

static const mxArray *c_feval(const emlrtStack *sp, const mxArray *m,
                              const mxArray *m1, const mxArray *m2,
                              const mxArray *m3, emlrtMCInfo *location)
{
  const mxArray *pArrays[4];
  const mxArray *m4;
  pArrays[0] = m;
  pArrays[1] = m1;
  pArrays[2] = m2;
  pArrays[3] = m3;
  return emlrtCallMATLABR2012b((emlrtConstCTX)sp, 1, &m4, 4, &pArrays[0],
                               "feval", true, location);
}

static const mxArray *feval(const emlrtStack *sp, const mxArray *m,
                            const mxArray *m1, const mxArray *m2,
                            emlrtMCInfo *location)
{
  const mxArray *pArrays[3];
  const mxArray *m3;
  pArrays[0] = m;
  pArrays[1] = m1;
  pArrays[2] = m2;
  return emlrtCallMATLABR2012b((emlrtConstCTX)sp, 1, &m3, 3, &pArrays[0],
                               "feval", true, location);
}

void b_warning(const emlrtStack *sp, const char_T varargin_1[23],
               const char_T varargin_2[23])
{
  static const int32_T iv[2] = {1, 7};
  static const int32_T iv1[2] = {1, 7};
  static const int32_T iv2[2] = {1, 33};
  static const int32_T iv3[2] = {1, 23};
  static const int32_T iv4[2] = {1, 23};
  static const char_T msgID[33] = {'M', 'A', 'T', 'L', 'A', 'B', ':', 'o', 'd',
                                   'e', '4', '5', ':', 'I', 'n', 't', 'e', 'g',
                                   'r', 'a', 't', 'i', 'o', 'n', 'T', 'o', 'l',
                                   'N', 'o', 't', 'M', 'e', 't'};
  static const char_T b_u[7] = {'m', 'e', 's', 's', 'a', 'g', 'e'};
  static const char_T u[7] = {'w', 'a', 'r', 'n', 'i', 'n', 'g'};
  emlrtStack st;
  const mxArray *b_y;
  const mxArray *c_y;
  const mxArray *d_y;
  const mxArray *e_y;
  const mxArray *m;
  const mxArray *y;
  st.prev = sp;
  st.tls = sp->tls;
  y = NULL;
  m = emlrtCreateCharArray(2, &iv[0]);
  emlrtInitCharArrayR2013a((emlrtConstCTX)sp, 7, m, &u[0]);
  emlrtAssign(&y, m);
  b_y = NULL;
  m = emlrtCreateCharArray(2, &iv1[0]);
  emlrtInitCharArrayR2013a((emlrtConstCTX)sp, 7, m, &b_u[0]);
  emlrtAssign(&b_y, m);
  c_y = NULL;
  m = emlrtCreateCharArray(2, &iv2[0]);
  emlrtInitCharArrayR2013a((emlrtConstCTX)sp, 33, m, &msgID[0]);
  emlrtAssign(&c_y, m);
  d_y = NULL;
  m = emlrtCreateCharArray(2, &iv3[0]);
  emlrtInitCharArrayR2013a((emlrtConstCTX)sp, 23, m, &varargin_1[0]);
  emlrtAssign(&d_y, m);
  e_y = NULL;
  m = emlrtCreateCharArray(2, &iv4[0]);
  emlrtInitCharArrayR2013a((emlrtConstCTX)sp, 23, m, &varargin_2[0]);
  emlrtAssign(&e_y, m);
  st.site = &nd_emlrtRSI;
  b_feval(&st, y, c_feval(&st, b_y, c_y, d_y, e_y, &emlrtMCI), &b_emlrtMCI);
}

void warning(const emlrtStack *sp)
{
  static const int32_T iv[2] = {1, 7};
  static const int32_T iv1[2] = {1, 7};
  static const int32_T iv2[2] = {1, 34};
  static const int32_T iv3[2] = {1, 23};
  static const char_T msgID[34] = {'M', 'A', 'T', 'L', 'A', 'B', ':', 'o', 'd',
                                   'e', 'a', 'r', 'g', 'u', 'm', 'e', 'n', 't',
                                   's', ':', 'R', 'e', 'l', 'T', 'o', 'l', 'I',
                                   'n', 'c', 'r', 'e', 'a', 's', 'e'};
  static const char_T varargin_1[23] = {' ', ' ', '2', '.', '2', '2', '0', '4',
                                        '4', '6', '0', '4', '9', '2', '5', '0',
                                        '3', '1', '3', 'e', '-', '1', '4'};
  static const char_T b_u[7] = {'m', 'e', 's', 's', 'a', 'g', 'e'};
  static const char_T u[7] = {'w', 'a', 'r', 'n', 'i', 'n', 'g'};
  emlrtStack st;
  const mxArray *b_y;
  const mxArray *c_y;
  const mxArray *d_y;
  const mxArray *m;
  const mxArray *y;
  st.prev = sp;
  st.tls = sp->tls;
  y = NULL;
  m = emlrtCreateCharArray(2, &iv[0]);
  emlrtInitCharArrayR2013a((emlrtConstCTX)sp, 7, m, &u[0]);
  emlrtAssign(&y, m);
  b_y = NULL;
  m = emlrtCreateCharArray(2, &iv1[0]);
  emlrtInitCharArrayR2013a((emlrtConstCTX)sp, 7, m, &b_u[0]);
  emlrtAssign(&b_y, m);
  c_y = NULL;
  m = emlrtCreateCharArray(2, &iv2[0]);
  emlrtInitCharArrayR2013a((emlrtConstCTX)sp, 34, m, &msgID[0]);
  emlrtAssign(&c_y, m);
  d_y = NULL;
  m = emlrtCreateCharArray(2, &iv3[0]);
  emlrtInitCharArrayR2013a((emlrtConstCTX)sp, 23, m, &varargin_1[0]);
  emlrtAssign(&d_y, m);
  st.site = &nd_emlrtRSI;
  b_feval(&st, y, feval(&st, b_y, c_y, d_y, &emlrtMCI), &b_emlrtMCI);
}

/* End of code generation (warning.c) */
