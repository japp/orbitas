# -*- coding: utf-8 -*-
"""
Created on Sat Sep 12 18:51:21 2015

@author: japp


"""

import matplotlib.pyplot as plt
import binarylib as bl

plt.style.use('jjplot')

fig = plt.figure(figsize=(8, 5))
ax = plt.subplot()

oe = bl.load_elements("GJ804-new")

# Observaciones
t = bl.read_binary_obs('GJ804')
x, y = bl.rt2xy(t['rho'], t['theta'])


# Medidas
#OC_lines = bl.OC_lines(ax, t['rho'], t['theta'], t['date'], oe, color='k')
plt.plot(x, y, 'o')

# Modelo
xd, yd = bl.orbit_model_xy(oe)
plt.plot(xd, yd, 'k', zorder=-1)

# Origen primaria
bl.draw_origin()
bl.draw_scale_circle(ax)

# Linea de nodos
bl.draw_nodes_line(oe, ax)


# Leyenda NE con direccion
#bl.draw_compass(ax, direct=True)

plt.grid(False)

bl.set_plot_axis(xfmt=':.1f', yfmt=':.2f', half_yticks=True)

plt.gca().set_aspect('equal', adjustable='box')
plt.tight_layout()


#bl.plot_OC(obs, oe_new)

plt.show()