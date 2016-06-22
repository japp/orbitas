# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 18:04:59 2010

@author: japp
"""

from numpy import *

def par2canonic(a, b, alpha, z):
    h, k = z # ellipse center

    A = (b*cos(alpha))**2 + (a*sin(alpha))**2
    B = -2*cos(alpha)*sin(alpha)*(a**2-b**2)
    C = (b*sin(alpha))**2 + (a*cos(alpha))**2
    D = -2*A*h - k*B
    E = -2*C*k - h*B
    F = -(a*b)**2 + A*h**2 + B*h*k + C*k**2

    return A, B, C, D, E, F

#canonical_params = array([14., 23., 18, 2, -3, -31, -100])

a, b, c, d, f = array([ 0.14, 0.18, -0.23/2,  -0.03/2, -0.31/2 ])

G = mat([[a,c], [c,b]])

row = mat([d,f]).transpose()

r, s = array(G.I*row)
r, = r
s, = s

Omega = 0.5*arctan2((2*c),(a-b))

I = a*r**2 + b*s**2 + 2*c*r*s
A = a*cos(Omega)**2 + b*sin(Omega)**2 + c*sin(2*Omega)
B = a*sin(Omega)**2 + b*cos(Omega)**2 - c*sin(2*Omega)

ap = sqrt((1+I)/A)
bp = sqrt((1+I)/B)

print(
"""a: %f
b: %f
phi: %f
x0: %f
y0: %f
""" % (ap, bp, degrees(Omega), r, s))

print par2canonic(ap, bp, Omega, [r, s])

"""
A, B, C, D, E = [-2.71494563,  1.3022704,   0.01109172,  0.00505402, -0.09340989]
A, B, C, D, E = [-1.0, -3.0, 0.09, 0.02, 0.08]

a = C
b = D
c = E
d = -(A*C+B*D)
f = -(A*D+B*E)
g = A**2*C + 2*A*B*D + B**2*E -1.0


# Abreviaturas
q = (b**2-a*c)
r = sqrt((a-c)**2+4.0*b**2)
p = 2.0*(a*f**2 +c*d**2 + g*b**2 -2.0*b*d*f - a*c*g)

# Parametros canonicos
x0 = (c*d-b*f)/q
y0 = (a*f-b*d)/q
ap = sqrt(p/( q*(r - (a-c))))
bp = sqrt(p/( q*(-r - (a-c))))


if (b==0 and a<c):
    phi = 0
elif (b==0 and a>c):
    phi = pi/2.0
elif (b!=0 and a<c):
    phi = 0.5* arctan(2*b/(c-a))
elif (b!=0 and a>c):
    phi = pi/2.0 + 0.5* arctan(2*b/(c-a))


#print("num: %f, q: %f, root: %f" % (num, q, root))
print(
a: %f
b: %f
phi: %f
x0: %f
y0: %f
 % (ap, bp, phi*(180./pi), x0, y0))
"""