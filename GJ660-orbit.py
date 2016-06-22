# -*- coding: utf-8 -*-
"""
Created on Sat Sep 12 18:51:21 2015

@author: japp

Masas dinámicas de Malkov 2012 para GJ 660
mag1   mag2 SpT plx   Ext       Per        a      e      M1    M2
11.00 11.20  M 98.19  0.016 	   33.5 	0.69 	0.2100 	0.31 	0.74

Doble velocidad areolar:
 (2*pi*0.371*0.625)/33.5 = 0.043489958

"""

import matplotlib.pyplot as plt
import binarylib as bl

plt.style.use('jjplot')

fig = plt.figure(figsize=(8, 5))
ax = plt.subplot()

oe = bl.load_elements("GJ660")

# Observaciones
t = bl.read_binary_obs('GJ660')
x, y = bl.rt2xy(t['rho'], t['theta'])

# Modelo
xd, yd = bl.orbit_model_xy(**oe)
plt.plot(xd, yd, 'k')

# Medidas
OC_lines = bl.OC_lines(ax, t['rho'], t['theta'], t['date'], oe, color='k')
plt.plot(x, y, 'o')



# Origen primaria
bl.draw_origin()

# Linea de nodos
#bl.draw_nodes_line(oe, ax)


# Leyenda NE con direccion
#bl.draw_compass(ax, direct=True)

plt.grid(False)

bl.set_plot_axis(xfmt=':.1f', yfmt=':.2f', half_yticks=True)

plt.gca().set_aspect('equal', adjustable='box')
plt.tight_layout()


#bl.plot_OC(obs, oe_new)

plt.show()