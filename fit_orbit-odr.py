# -*- coding: utf-8 -*-
"""
Created on Sun Feb  5 17:43:33 2012

@author: japp
"""

from scipy.odr import *
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


def implicit_fcn(p,rho, theta):
     return (rho_model(p, theta) - rho)     
     

implicit_mod = Model(implicit_fcn, implicit=1, meta=dict(name='Sample Implicit Model', ref='ODRPACK UG, pg. 49'),)

implicit_dat = Data([
            rhos,
            thetas],
            1,
        )
    
implicit_odr = ODR(implicit_dat, implicit_mod, beta0=[a0, i0, w0, O0, e0])

out = implicit_odr.run()
out.pprint()

