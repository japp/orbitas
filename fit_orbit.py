# -*- coding: utf-8 -*-
"""

@author: japp
@date: Mon Feb 22 17:43:07 2010
@version: 0.0
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
a0 = 0.8
e0 = 0.15
i0 = radians(102.000)
O0 = radians(141.480)
w0 = radians(349.200)

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

   
T_O = thetas - O  
signs = where((T_O < 3*pi/2) * (T_O >= pi/2),+1,-1)

# Anomalia verdadera

quad = (ceil((thetas-O)/(2*pi)) + 8)%4
v0 = ceil(quad/2)*pi

v = signs*(thetas-O)*arctan(tan(thetas - O)/cos(i)) - w #+ v0
#v = v + 4*pi % (2*pi)

# Anomalia excentrica

quad_v = (ceil(v/(2*pi)) + 8 )%4
E0 = ceil(quad_v/2) * pi


#E = 2.0*arctan(sqrt((1-e)/(1+e))*tan(v/2.0)) 
E = arccos((e + cos(v))/(e*cos(v) + 1)) #+ E0
#E = E + 4*pi % (2*pi)

orbital = 2*pi*((times - min(times))/15.9)

# Anomalia media
M = E - e*sin(E) - orbital

ps = polyfit(times, M, 1)



p0 = a0, i0, w0, O0, e0 #, P0, T0

 QUAD = (ceil((theta-A[3])/!pi*2) + 8) mod 4
 V_0 = ceil(quad/2) * !pi
 V = atan(tan(theta - A[3])/cos(A[2])) - A[4] + v_0
 V = v + 4 * !pi mod (2 * !pi)

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
p2 = optimize.fmin(ef, p0, args=(rhos,thetas), maxiter=20000, maxfun=10000)

elements = {
"Semimajor axis (a)": p2[0],
"Inclination (i)": degrees(p2[1]),
"Periastron Longitude  (omega)": degrees(p2[2]%360),
"Ascending node position (Omega)": degrees(p2[3])%360,
"Eccentricity (e)": p2[4],
#"Period (P)": P0,
#"Epoch of periastron (T)" : T0
}

print("\n------------------------------------")
for k, v in elements.iteritems():
    print("%s: %f" % (k, v))
print("------------------------------------\n")





