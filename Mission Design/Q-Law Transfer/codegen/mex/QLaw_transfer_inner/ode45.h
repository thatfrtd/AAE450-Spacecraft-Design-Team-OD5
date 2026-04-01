/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * ode45.h
 *
 * Code generation for function 'ode45'
 *
 */

#pragma once

/* Include files */
#include "QLaw_transfer_inner_types.h"
#include "rtwtypes.h"
#include "emlrt.h"
#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Function Declarations */
void ode45(const emlrtStack *sp,
           const real_T c_ode_workspace_a_control_works[3],
           real_T ode_workspace_mdot, const real_T tspan[2],
           const real_T b_y0[7], real_T options_AbsTol,
           const real_T options_InitialStep_data[],
           const real_T options_MaxStep_data[], real_T options_RelTol,
           emxArray_real_T *varargout_1, emxArray_real_T *varargout_2);

/* End of code generation (ode45.h) */
