# -*- coding: utf-8 -*-
"""
Created on 20 Dic 2014
Revision: 2015-07-19

@author: japp
"""

import fitellipse
import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA
from scipy.optimize import minimize

def plotEllipse(a, b, phi=0, x0=0, y0=0, n=100):
    """
    Genera puntos de una elipse de lados a y b, angulo phi,
    centrado en x0, y0 con n puntos.
    """
    interv = (2*np.pi)/(n-1)
    phi = np.radians(phi)
    th = np.arange(0, 2*np.pi + interv, interv)
    x = a*np.cos(th)
    y = b*np.sin(th)
    c = np.cos(phi)
    s = np.sin(phi)
    th = x*c-y*s+x0
    y = x*s+y*c+y0
    x = th
    return x, y


def fitEllipse2(x, y):
    x = x[:, np.newaxis]
    y = y[:, np.newaxis]
    D = np.hstack((x*x, x*y, y*y, x, y, np.ones_like(x)))
    S = np.dot(D.T, D)
    C = np.zeros([6, 6])
    C[0, 2] = C[2, 0] = 2; C[1, 1] = -1
    E, V = LA.eig(np.dot(LA.inv(S), C))
    n = np.argmax(np.abs(E))
    a = V[:, n]
    return a

def rt2xy(rho, theta):
    """Pasa de coordenadas rho y theta a (x,y) y rota 90º
       para tener N arriba
    """
    x = rho*np.cos(np.radians(theta+90))
    y = rho*np.sin(np.radians(theta+90))
    return x, y


# observaciones
t, th, rh = np.loadtxt('../orbitas/binaries_data/GJ1245Aab_obs.dat',
                       unpack=True, usecols=(0, 1, 2), skiprows=4)

# Polares a cartesianas
ox, oy = rt2xy(rh, th)
xs = np.zeros((2, len(ox)))
xs[0] = ox
xs[1] = oy

z, a, b, alpha = fitellipse.fitellipse(xs)

print ("Resultado del ajuste")
print (" a = %.3f\n b = %.3f\n alfa = %.3f" % (a, b, alpha))
print (" Centro: x={:.3f} y={:.3f}".format(*z))

print("fitEllipse2 ----------------")
ellipse_params0 = fitEllipse2(ox, oy)
print(ellipse_params0)


def ellipse_fn(p):
    # Ax² + 2Hxy + By² + 2Gx + 2Fy + 1 = 0

    return np.sum((p[0]*x**2 + 2*p[1]*x*y + p[2]*y**2 + 2*p[3]*x + 2*p[4]*y + p[5])**2.0, axis=0)

#ellipse_params_NM = minimize(ellipse_fn, ellipse_params0, method='Nelder-Mead')
#print(ellipse_params_NM)

# Centro de la elipse (primaria)
plt.axhline(0, c='k', alpha=0.5)
plt.axvline(0, c='k', alpha=0.5)
plt.plot(0, 0, '+k', ms=25, mew=3)

# Elipse ajustada
x, y = plotEllipse(a, b, np.degrees(alpha), z[0], z[1])
plt.plot(x, y, 'k-')

# Medidas
plt.plot(xs[0], xs[1], 'ro')

plt.show()

