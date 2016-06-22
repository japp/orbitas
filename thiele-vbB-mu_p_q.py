# -*- coding: utf-8 -*-
"""
Created on Mon Sep 21 19:51:12 2015

@author: japp
"""

import numpy as np
from scipy.optimize import broyden1, brute

L12 = 9.18
L23 = 6.98
L13 = 47.22



def equations(params):
    mu, p, q = params
    
    eqs = [mu*L12 - p + np.sin(p),
           mu*L23 - q + np.sin(q),
           mu*L13 - p - q + np.sin(p + q)
           ]
    return eqs

mu1, p1, q1 = np.array([0.03, 1.2, 1.])

res_mu0, res_p0, res_q0 = np.array([0.01, 0.01, 0.01])

for mu0 in np.arange(0.01, 0.08, 0.01):
    try:
        for p0 in np.arange(1, 2, 0.1):
            for q0 in np.arange(1, 2, 0.1):

                mu1, p1, q1 = broyden1(equations, (mu0, p0, q0))
                res_mu, res_p, res_q = equations((mu1, p1, q1))
                
                if (abs(np.array([res_mu, res_p, res_q])) <
                    abs(np.array([res_mu0, res_p0, res_q0]))).all():
                    res_mu0, res_p0, res_q0 = res_mu, res_p, res_q
                
                    mu, p, q = mu1, p1, q1
    except:
        continue


print(mu, p, q)

