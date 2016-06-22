# -*- coding: utf-8 -*-
"""
Created on Sat Sep 12 18:51:21 2015

@author: japp
"""

import matplotlib.pyplot as plt
import binarylib as bl

plt.style.use('jjplot')

oe = bl.load_WDS_elements('19074+3230')  # Nombre de GJ747 en WDS
# Modelo segun elementos orbitales
xd, yd = bl.orbit_model_xy(oe)

oe_new = bl.load_elements("GJ747-new")
# Modelo segun NUEVOS elementos orbitales
xd_new, yd_new = bl.orbit_model_xy(oe_new)

# Observaciones
t = bl.read_binary_obs('GJ747')
x, y = bl.rt2xy(t['rho'], t['theta'])

fig = plt.figure(figsize=(18, 5))
ax = plt.subplot()

# Medidas
OC_lines = bl.OC_lines(ax, t['rho'], t['theta'], t['date'], oe_new, color='k')
plt.plot(x, y, 'o')

# Modelo Ségransan
plt.plot(xd, yd, 'k:', color='gray')


# Modelo nuevo
plt.plot(xd_new, yd_new, 'k')

# Origen primaria
bl.draw_origin()

# Linea de nodos
bl.draw_nodes_line(oe_new, ax)


# Leyenda NE con direccion
bl.draw_compass(ax, direct=True)

plt.grid(False)

bl.set_plot_axis(xfmt=':.1f', yfmt=':.2f', half_yticks=True)

plt.gca().set_aspect('equal', adjustable='box')
plt.tight_layout()


#bl.plot_OC(obs, oe_new)

plt.show()