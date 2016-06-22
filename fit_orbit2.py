# -*- coding: utf-8 -*-
"""
Created on Sun Feb  5 17:43:33 2012

@author: japp
"""

from scipy import optimize

import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.lines import Line2D
from pylab import *

""" Torres"""
P = 15.643
T = 1992.297
a = 0.926
e = 0.295
i = radians(103.000)
O = radians(143.480)
w = radians(347.200)


""" SRH
P = 15.184596
T = 1992.297
a = 0.94600318
e =  0.32667310
i =  radians(76.950024)
O = radians(322.77072%180.)
w = radians(191.80015)
"""
# Elem Orb iniciales para el ajuste
P0 = 15.0
T0 = 1991.7
a0 = 0.9
e0 = 0.15
i0 = radians(102.000)
O0 = radians(141.480)
w0 = radians(349.200)
# p0 = a0, i0, w0, O0, e0
ranges = (
(0.7,1.1),                           # a 
(radians(99.), radians(110.)),       # i
(radians(340.), radians(359.)),      # omega
(radians(135.), radians(145.)),      # Omega
(0.1,0.5)                            # e
)

# Observaciones
obs = array([
[1.098, 322.7, 1986.4499],
[0.220, 278.0, 1990.04],
[0.179, 248., 1990.181],
[0.177, 246., 1990.183],
[0.177, 226., 1990.359],
[0.186, 210.5, 1990.4399],
[0.43, 160., 1991.237],
[0.436, 161.6, 1991.422],
[0.602, 135.8, 1993.3514],
[0.2359, 26.04, 1995.5733],
[0.416, 353.2, 1996.2923],
[0.4639, 350.53, 1996.4606],
[0.7344,338.44, 1997.4241],
[0.9245, 333.46, 1998.2510],
[0.424, 160.6, 2007.150],
[0.631, 145.9, 2008.027],
[0.638, 143.83, 2008.305],
[0.648, 142.57, 2008.406],
[0.604, 136.1, 2009.039],
[0.579, 132.6, 2009.212],
[0.548, 131.5, 2009.351],
[0.231, 78.7,  2010.479]
])


rhos = obs[:,0]
thetas = radians(obs[:,1]) 
times = obs[:,2]


p0 = a0, i0, w0, O0, e0 #, P0, T0

def t_anom(O, w, i, theta):
    true_anomaly = arctan(tan(theta-O)/cos(i)) - w
    return true_anomaly
    
    

def rho_model(p, theta):
    v = t_anom(p[3], p[2], p[1], theta)
    num = p[0]*(1-p[4]**2)*cos(v+p[2])
    denom = (1 + p[4]*cos(v))*cos(theta - p[3])
    return abs(num/denom)


""" """
ef = lambda p,rho,theta: ( (rho_model(p,theta) - rho)**2 ).sum()

# Downhill Simplex
#p2 = optimize.fmin(ef, p0, args=(rhos,thetas), maxiter=100000, maxfun=10000)

# Fuerza bruta
p2 = optimize.brute(ef, ranges, args=(rhos,thetas), Ns=5)

elements = {
"Semimajor axis (a)": p2[0],
"Inclination (i)": degrees(p2[1]),
"Periastron Longitude  (omega)": degrees(p2[2]),
"Ascending node position (Omega)": degrees(p2[3]),
"Eccentricity (e)": p2[4],
#"Period (P)": P0,
#"Epoch of periastron (T)" : T0
}

print("\n------------------------------------")
for k, v in elements.iteritems():
    print("%s: %f" % (k, v))
print("------------------------------------\n")


# Datos
xs = cos(thetas)*rhos
ys = sin(thetas)*rhos
plot(xs,ys, 'ro')

# Modelo
ts = arange(0,2*pi,0.01)
xd = cos(ts)*rho_model(p2, ts)
yd = sin(ts)*rho_model(p2, ts)
plot(xd,yd)

show()

"""
from openopt import DFP


lb = [0.7,radians(99.),radians(340.), radians(135.),0.1]
ub = [1.1, radians(110.), radians(359.), radians(145.),0.5]
p = DFP(rho_model, p0, rhos, thetas, lb=lb, ub=ub)


r = p.solve('nlp:ralg', plot=0, iprint = 10)
print('solution: '+str(r.xf)+'\n||residuals||^2 = '+str(r.ff)+'\nresiduals: '+str([f(p.xf, rhos[i])-thetas[i] for i in xrange(len(thetas))]))

"""
