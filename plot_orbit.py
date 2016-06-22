# -*- coding: utf-8 -*-
"""
Created on Sat Aug  1 11:45:55 2015

@author: japp

Dibuja la orbita de una binaria a partir de sus elementos orbitales

"""

#import numpy as np
import matplotlib.pyplot as plt
#from matplotlib.lines import Line2D
import binarylib as bl

oe = bl.load_elements("binaries_data/GJ473.eo.txt")

ax = plt.subplot()



# Orbita
xd, yd = bl.orbit_model_xy(oe)
plt.plot(xd, yd, 'k')

# Origen primaria
bl.draw_origin()

# Linea de nodos
bl.draw_nodes_line(oe, ax)

# Leyenda NE con direccion
bl.draw_compass(ax,  0.75, -0.5, 0.1, direct=False)

plt.show()
