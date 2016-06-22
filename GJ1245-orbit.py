# -*- coding: utf-8 -*-
"""
Created on Sat Sep 12 18:51:21 2015

@author: japp

Orbita de la binaria GJ 1245 A (MCY   3Aa,Ab)
"""

import matplotlib.pyplot as plt
import binarylib as bl

plt.style.use('jjplot')

#oe = bl.load_WDS_elements('19539+4425')  # Nombre de GJ1245 en WDS
#oe['O'] = 80 +180
#oe = bl.load_elements("GJ1245A-law")
oe = bl.load_elements("GJ1245A-japp")

# Modelo segun elementos orbitales
xd, yd = bl.orbit_model_xy(oe)

# Observaciones
t = bl.read_binary_obs('GJ1245Aab')
x, y = bl.rt2xy(t['rho'], t['theta'])

fig = plt.figure(figsize=(8, 8))
ax = plt.subplot()

# Lineas O-C
OC_lines = bl.OC_lines(ax, t['rho'], t['theta'], t['date'],
                       oe, color='k', alpha=0.3)
# Observaciones
for i, date in enumerate(t):

    if t['obs'][i] == "Mcy":
        marker = 'ys'
    if t['obs'][i] == "Shd":
        marker = 'g^'
    if t['obs'][i] == "Law":
        marker = 'b+'
    if t['obs'][i] == "Jod":
        marker = 'ro'

    pl = plt.plot(x[i], y[i], marker)


# Modelo publicado
plt.plot(xd, yd, 'k-', color='gray')

# Origen primaria
bl.draw_origin()
bl.draw_scale_circle(ax)

# Linea de nodos
bl.draw_nodes_line(oe, ax)

# Leyenda NE con direccion
bl.draw_compass(ax, direct=False, y=-0.65, l=0.1)

plt.grid(False)

bl.set_plot_axis(xfmt=':.1f', yfmt=':.1f', half_yticks=False)

plt.gca().set_aspect('equal', adjustable='box')
plt.tight_layout()

#bl.plot_OC(obs, oe_new)

plt.show()