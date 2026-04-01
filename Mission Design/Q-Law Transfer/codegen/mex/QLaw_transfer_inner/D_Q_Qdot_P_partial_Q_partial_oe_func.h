/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * D_Q_Qdot_P_partial_Q_partial_oe_func.h
 *
 * Code generation for function 'D_Q_Qdot_P_partial_Q_partial_oe_func'
 *
 */

#pragma once

/* Include files */
#include "rtwtypes.h"
#include "emlrt.h"
#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Function Declarations */
real_T c_D_Q_Qdot_P_partial_Q_partial_(
    const emlrtStack *sp, const real_T in1[5], real_T L1, const real_T in3[5],
    const real_T in4[5], real_T m1, real_T n1, real_T r1, real_T F_max1,
    real_T W_p1, real_T r_p_min1, real_T k1, real_T D[3], real_T *Qdot,
    real_T *P, real_T partial_Q_partial_oe[5]);

/* End of code generation (D_Q_Qdot_P_partial_Q_partial_oe_func.h) */
