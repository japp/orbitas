# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 11:50:05 2010

@author: japp
@description: calcula los elementos orbitales de una binaria a partir de medidas
(rho, theta) usando el metodo de Asada et al. (1994, 1997)
@version: 1.0
"""

from numpy import *

def u(a,b, x, y):
    return arctan2((a*y),(b*x))
    
def rt2xy(rho, theta):
    """Pasa de coordenadas rho y theta a (x,y) y rota 90º para tener N arriba"""
    x = rho*cos((theta)*(pi/180.))
    y = rho*sin((theta)*(pi/180.))
    return x, y

x0, y0, b, a, alpha = (
0.119274, 0.19294178, 
0.18852648378457651, 
0.91132208933263736, 
2.5296087217674175
)
"""
# Asada test
z, a,b, alpha = (
(-0.01173414,  0.02960827),
0.97682619131519521,
0.77347746049757993,
-2.6143520153986577
)
"""
# Puntos base
cpa = array([
[0.602, 135.8, 1993.3514],
[0.4639, 350.53, 1996.4606],
[0.7344,338.44, 1997.4241],
[0.9245, 333.46, 1998.2510]
])

x1, y1 = rt2xy(*cpa[0,:2])
x2, y2 = rt2xy(*cpa[1,:2])
x3, y3 = rt2xy(*cpa[2,:2])
x4, y4 = rt2xy(*cpa[3,:2])

t1, t2, t3, t4 = cpa[:,2]

"""
# Asada test
t1, x1, y1 = 1,0.372003,0.838658  
t2, x2, y2 =2,-0.0648698,0.831542  
t3, x3, y3 =3,-0.404177,0.696231  
t4, x4, y4 =4,-0.646827,0.509083  
t5, x5, y5 =5,-0.809209,0.304280  
"""
t21 = t2 - t1
t31 = t3 - t1
t32 = t3 - t2
t42 = t4 - t2
t43 = t4 - t3

u1, u2, u3, u4 =  u(a,b, x1, y1), u(a,b, x2, y2), u(a,b, x3, y3), u(a,b, x4, y4)

#u1, u2, u3, u4 = arccos(x1/a), arccos(x2/a), arccos(x3/a), arccos(x4/a)

A1 = t21*sin(u3) + t32*sin(u1) - t31*sin(u2)
A2 = t21*cos(u3) + t32*cos(u1) - t31*cos(u2)
A3 = t21*u3 + t32*u1 - t31*u2

B1 = t32*sin(u4) + t43*sin(u2) - t42*sin(u3)
B2 = t32*cos(u4) + t43*cos(u2) - t42*cos(u3)
B3 = t32*u4 + t43*u2 - t42*u3


xe = -a *(A2*B3 - A3*B2)/(A1*B2 - A2*B1)
ye =  b *(A3*B1 - A1*B3)/(A1*B2 - A2*B1)

e = sqrt((xe/a)**2 + (ye/b)**2)

print xe, ye, e


