/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * QLaw_transfer_inner.h
 *
 * Code generation for function 'QLaw_transfer_inner'
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
    real_T u_nd[3], boolean_T *not_coast);

void QLaw_transfer_inner_anonFcn2(const emlrtStack *sp,
                                  const real_T a_control_workspace_u_nd[3],
                                  real_T mdot, const real_T x[7],
                                  real_T varargout_1[7]);

/* End of code generation (QLaw_transfer_inner.h) */
