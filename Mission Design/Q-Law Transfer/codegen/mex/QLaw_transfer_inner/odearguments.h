/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * odearguments.h
 *
 * Code generation for function 'odearguments'
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
real_T odearguments(const emlrtStack *sp,
                    const real_T c_odeFun_workspace_ODEFunction_[3],
                    real_T d_odeFun_workspace_ODEFunction_,
                    const real_T tspan[2], const real_T b_y0[7],
                    real_T options_AbsTol,
                    const real_T options_InitialStep_data[],
                    const real_T options_MaxStep_data[], real_T options_RelTol,
                    real_T tspanOut[2], real_T *tfinal, real_T *tdir,
                    real_T f0[7], real_T *threshold, real_T *rtol, real_T *hmax,
                    real_T htry_data[], int32_T htry_size[2], real_T *htspan);

/* End of code generation (odearguments.h) */
