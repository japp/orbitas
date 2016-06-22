# -*- coding: utf-8 -*-

"""
@author: japp
@date:  Feb 11  2011
@version: 0.0

Calculo de periodo y paso por el periastro de una binaria conocidos
los elementos geometricos de la orbita.

"""

import numpy as np
import binarylib as bl
import matplotlib.pyplot as plt

# Observaciones
obs = array([
#[1.098, 322.7, 1986.4499],
#[0.220, 278.0, 1990.04],
#[0.179, 248., 1990.181],
#[0.177, 246., 1990.183],
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


obs = bl.read_binary_obs('GJ660')

rhos = obs['rho']
thetas = np.radians(obs['theta'])
times = obs['date']

oe = bl.load_elements('GJ660')


def true_anomaly_eq(theta, Omega, inc, omega):
    """Anomalia verdadera"""

    return np.arctan2(np.tan(theta - Omega), np.cos(inc)) - omega


def ecc_anomaly_eq(v, e):
    """Anomalia excentrica"""
    en = np.sqrt((1.0 + e)/(1.0 - e))
    return 2.0*np.arctan(np.tan(v/2)/en)

v = true_anomaly_eq(thetas, oe['O'], oe['i'], oe['w'])

E = ecc_anomaly_eq(v, oe['e'])

# Anomalia media
M = E - oe['e']*np.sin(E)

M2 = np.array([])
for M0 in M:
    if M0 > 0:
        M0 = M0 - 2*np.pi
    M2 = np.append(M2, M0)
  
p = np.polyfit(times, M2, 1)
ts = np.arange(times.min() - 2, times.max() + 2)

print(p)

# Periodo y paso por el periastro
P = 2*np.pi/p[0]
T0 = -p[1]/p[0]

print("Periodo: {:.1f}".format(P))
print("T0: {:.1f}".format(T0))

plt.plot(times, M2, 'o')
plt.plot(ts, np.polyval(p, ts), '-k')

plt.show()