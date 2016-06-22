#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "J A Perez Prieto <japp@denebola.org>"
__version__ = "$ Revision: Sept 2008 $"

"""
Calcula efemerides de estrellas binarias
a partir de los elementos orbitales
"""

from scipy import pi, sin, cos, tan, arctan, sqrt
from numpy import arange, zeros
from pylab import *
"""  GJ 473 """

P = 15.643                  # Periodo (años)
T = 1992.297                # Paso por el periastro (años)
e = 0.295                  # Ecentricidad
a = 0.926                  # Semieje mayor (seg de arco)
i = 103.0                  # Inclinacion de la orbita (grados)
O_node = 143.48              # Posicion del nodo ascendente (Omega, grados)
o_peri = 347.2            # Argumento del periastro (omega, grados)

"""
P = 41.623                  # Periodo (años)
T = 1934.008                # Paso por el periastro (años)
e = 0.2763                  # Ecentricidad
a = 0.907                  # Semieje mayor (seg de arco)
i = 59.025                  # Inclinacion de la orbita (grados)
O_node = 23.717              # Posicion del nodo ascendente (Omega, grados)
o_peri = 219.907            # Argumento del periastro (omega, grados)
"""
#01 / 01 /2007 		    0.4352 	 159.96
#31 / 12 /2009 		    0.3631 	 116.52
t = 2009.99


print_out = 1   # imprime o no la salida en tabla
def binary_ephem(P, T, e, a, i, O_node, o_peri, t):
   # Grados a radianes
   d2rad = pi/180.
   rad2d = 180./pi
   i = i*d2rad
   O_node = (O_node*d2rad)%(2*pi)
   o_peri = (o_peri*d2rad)%(2*pi)
 
   # Anomalia media
   M = ((2.0*pi)/P)*(t - T)  # radianes
    
   if M >2*pi: M = M - 2*pi 
   M=M%(2*pi)

   # Anomalia excentrica (1ra aproximacion)
   E0 = M  + e*sin(M) + (e**2/M) * sin(2.0*M)

   for itera in range(15):
      M0 = E0 - e*sin(E0)
      E0 = E0 + (M-M0)/(1-e*cos(E0))

   true_anom = 2.0*arctan(sqrt((1+e)/(1-e))*tan(E0/2.0))
 
   #radius = (a*(1-e**2))/(1+e*cos(true_anom))
   radius = a*(1-e*cos(E0))

   theta = arctan( tan(true_anom + o_peri)*cos(i) ) + O_node
   rho = radius * (cos(true_anom + o_peri)/cos(theta - O_node))
   
   # revuelve rho ("), theta (grad), Anomalia excentrica (grad), Anomalia verdadera (grad)
   return rho, (theta*rad2d)%360. #, E0*rad2d, M*rad2d, true_anom*rad2d

dates = arange(2008.0, 2014.0, 0.2)
date_count = 0
coords = zeros([2, len(dates)])
"""
for t in dates:
   rho, theta, E, M, true_anom = binary_ephem(P, T, e, a, i, O_node, o_peri, t)
   #print rho, theta, t 
   coords[0][date_count] = rho * sin(theta*(pi/180.))
   coords[1][date_count] = rho * cos(theta*(pi/180.))
   #print coords[0][date_count], coords[1][date_count], date
   date_count=date_count+1
"""

plot_out = False
if (plot_out): 
   plot(coords[0,:], coords[1,:])
   show()

rho, theta, E, M, true_anom = binary_ephem(P, T, e, a, i, O_node, o_peri, t)
if (print_out): 
   print "\n------------------------------------------------"
   print "Fecha: ", t
   print "Rho (\"):", rho
   print "Theta (deg):", theta
   print "Anomalia ecentrica, E (deg):", E
   print "Anomalia media, M (deg):", M
   print "Anomalia verdadera, nu (deg):", true_anom
   print "------------------------------------------------\n"

