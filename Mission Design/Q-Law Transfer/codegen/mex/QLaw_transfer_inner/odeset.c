/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * odeset.c
 *
 * Code generation for function 'odeset'
 *
 */

/* Include files */
#include "odeset.h"
#include "QLaw_transfer_inner_internal_types.h"
#include "rt_nonfinite.h"

/* Function Definitions */
struct_T odeset(real_T varargin_2, real_T varargin_4)
{
  struct_T options;
  options.RelTol = varargin_2;
  options.AbsTol = varargin_4;
  return options;
}

/* End of code generation (odeset.c) */
