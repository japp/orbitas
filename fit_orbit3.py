#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 20 Dic 2014

@author: japp

Ajusta la órbita de una binaria usando medidas astrometricas (rho, theta)
Necesita valores iniciales de los parametros orbitales
Ajusta un modelo analitico de la orbita por Decamps (2005)

Decamps 2005
Celestial Mechanics and Dynamical Astronomy (2005) 92:381–402
DOI 10.1007/s10569-005-2629-8

"""

from scipy import optimize

import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.lines import Line2D
from pylab import *
from astropy.io import ascii
import binarylib as bl

""" Torres
P = 15.643
T = 1992.297
a = 0.926
e = 0.295
i = radians(103.000)
O = radians(143.480)
w = radians(347.200)
"""

""" SRH
P = 15.184596
T = 1992.297
a = 0.94600318
e =  0.32667310
i =  radians(76.950024)
O = radians(322.77072%180.)
w = radians(191.80015)

# Elem Orb iniciales para el ajuste
P0 = 15.0
T0 = 1991.7
a0 = 0.9
e0 = 0.3
i0 = radians(102.0)
O0 = radians(143.0)
w0 = radians(347.)
"""

ranges = (
(0.1, 1.3),                          # a
(70., 89.),       # i
(300., 359.),      # omega
(80., 89.),      # Omega
(0.1,0.5)                            # e
)

# Observaciones
obs = ascii.read('binaries_data/gl747_obs.dat')

rhos = obs['rho']
thetas = radians(obs['theta'])
times = obs['epoch']
rho_err = obs['rho_err']

oe = bl.load_elements("binaries_data/GJ747.eo.txt")
#p0 = a0, i0, w0, O0, e0 #, P0, T0
p0 = oe['a'], oe['i'], oe['w'], oe['O'], oe['e']


# Modelo orbital de Decamps
def G(theta, i, O, w):
    numer = cos(i)*sin(theta+O)*cos(w) + sin(w)*cos(theta+O)
    denom = sqrt(cos(theta+O)**2 + sin(theta+O)**2*cos(i)**2)
    return numer/denom


def rho_model(p, theta):
    G_val = G(theta, p[1], p[3], p[2])

    numer = p[0]*(1-p[4]**2)*cos(i)
    denom = (1 + p[4]*G_val)*sqrt(cos(theta+p[3])**2 + sin(theta+p[3])**2 * cos(i)**2)
    return abs(numer/denom)


""" """
ef = lambda p, rho, theta: (((rho_model(p, theta) - rho)/rho_err)**2).sum()

# Downhill Simplex
p2 = optimize.fmin(ef, p0, args=(rhos, thetas))

# Simplex con limites (bounds). Por ahora no funciona
#p2 = optimize.fmin_slsqp(ef, p0, args=(rhos,thetas),  bounds=ranges)

# Fuerza bruta
#p2 = optimize.brute(ef, ranges, args=(rhos, thetas))

# Simulated annealing
#p2, k = optimize.anneal(ef, p0, args=(rhos,thetas))

elements = {
    "a": p2[0],
    "i": degrees(p2[1])%180,
    "w": degrees(p2[2])%360,
    "O": degrees(p2[3])%360,
    "e": p2[4],
    "P": oe['P'],
    "T": oe['T']
}

print("\n------------------------------------")
for k, v in elements.iteritems():
    print("%s: %f" % (k, v))
print("------------------------------------\n")


# Datos
xs, ys = bl.rt2xy(rhos, obs['theta'])
plot(xs, ys, 'ro')

# Modelo
"""
ts = arange(0, 2*pi, 0.01)
xd = cos(ts)*rho_model(p2, ts)
yd = sin(ts)*rho_model(p2, ts)
"""
xd, yd = bl.orbit_model_xy(oe)
plot(xd, yd)


show()

"""
from openopt import DFP


lb = [0.7,radians(99.),radians(340.), radians(135.),0.1]
ub = [1.1, radians(110.), radians(359.), radians(145.),0.5]
p = DFP(rho_model, p0, rhos, thetas, lb=lb, ub=ub)


r = p.solve('nlp:ralg', plot=0, iprint = 10)
print('solution: '+str(r.xf)+'\n||residuals||^2 = '+str(r.ff)+'\nresiduals: '+str([f(p.xf, rhos[i])-thetas[i] for i in xrange(len(thetas))]))

"""
