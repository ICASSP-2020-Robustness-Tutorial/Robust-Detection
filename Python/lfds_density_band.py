#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 10:55:32 2020

@author: mfauss
"""

import sys
import os

import numpy as np
from scipy import optimize

sys.path.append(os.path.abspath("../../Helper_Functions"))
from robust_detection_helpers import *


def lfds_density_band(p_min, p_max, dx, alpha=0.0, q0_init=np.nan, q1_init=np.nan, p_norm=np.inf, tol=1e-6, itmax=100):

    # sanity checks
    if p_min.shape[0] == p_min.shape[0] == 2 and p_min.shape[1] == p_max.shape[1]:
        K = p_min.shape[1]
        p0_min = p_min[0, :]
        p0_max = p_max[0, :]
        p1_min = p_min[1, :]
        p1_max = p_max[1, :]
    else:
        raise ValueError("'p_min' and 'p_max' must be of size 2xK")

    if not is_nonnegative_scalar(alpha):
        raise ValueError("The parameter 'alpha' must be a nonegative scalar")
        
    # if not is_valid_density_band(p0_min, p0_max, dx):
    #     raise ValueError('Invalid density band under H0.')
        
    # if not is_valid_density_band(p1_min, p1_max, dx):
    #     raise ValueError('Invalid density band under H1.')
    
    dist = np.inf
    nit = 0

    # while dist > tol and nit < itmax
        
    #     # assigne updated lfds
    #     q0 = q0_new;
    #     q1 = q1_new;
        
    #     % update q0
    #     def func0(c0):
    #         np.sum(np.minimum(p0_max, np.maximum(c0*(alpha*q0 + q1), p0_min)))*dx - 1;
        
    #     c0 = fzero(func0,c0);
    #     q0_new = min(p0_max, max(c0*(alpha*q0 + q1), p0_min));
            
    #     % update q1 using q0_new (!)
    #     func1 = @(c1) sum(min(p1_max, max(c1*(q0_new + alpha*q1), p1_min))) - 1/dx;
    #     c1 = fzero(func1, c1);
    #     q1_new = min(p1_max, max(c1*(q0_new + alpha*q1), p1_min));
        
    #     % calculate sup-norm
    #     dist = max( vecnorm(q0_new-q0, p), vecnorm(q1_new-q1, p) );
          
    #     % count iterations
    #     nit = nit+1;
           
    # end