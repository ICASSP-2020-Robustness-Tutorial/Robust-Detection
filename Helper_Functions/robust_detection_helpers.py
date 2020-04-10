#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 16:57:16 2020

@author: mfauss
"""

import numpy as np


def is_nonnegative_scalar(x):
    return np.isscalar(x) and np.isreal(x) and x >= 0


def is_positive_scalar(x):
    return np.isscalar(x) and np.isreal(x) and x > 0


def is_valid_density_band(p_min, p_max, dx):
    return (
        is_positive_scalar(dx)
        and np.all(np.isreal(p_min))
        and p_min.ndim == 1
        and np.all(np.isreal(p_max))
        and p_max.ndim == 1
        and p_min.size == p_max.size
        and np.sum(p_min) * dx <= 1
        and np.sum(p_max) * dx >= 1
        and np.all(p_min >= 0)
        and np.all(p_min <= p_max)
    )
