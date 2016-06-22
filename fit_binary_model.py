# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 15:59:56 2015

@author: japp


Ajusta la orbita de una binaria astrométrica minimizando
el O-C usando la formula de efemérides.


"""

import binarylib as bl
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize, fmin, anneal, basinhopping


def binary_position(t, oe):
    """
    Da la posición (rho, theta) en un tiempo t, de una binaria con
    elementos orbitales oe

    t: Año besseliano YYYY.YYYY
    oe: Lista de elementos orbitales, con angulos en grados



    """

    #P, T, e, a, i, O, w = oe['P'], oe['T'], oe['e'], oe['a'], oe['i'], oe['O'], oe['w']
    P, T, e, a, i, O, w = oe

    # Grados a radianes
    """
    i = np.radians(i)
    O = np.radians(O) % (2*np.pi)
    w = np.radians(w) % (2*np.pi)
    """

    # Anomalia media
    M = ((2.0*np.pi)/P)*(t - T)  # radianes

    if M > 2*np.pi:
        M = M - 2*np.pi
    M = M % (2*np.pi)

    # Anomalia excentrica (1ra aproximacion)
    # Para evitar dividir entre cero, si M es cero (o < 1e-8)
    # le doy valor 1e-8 como primera proximación
    if M < 1e-8:
        M = 1e-8
    E0 = M + e*np.sin(M) + (e**2/M) * np.sin(2.0*M)

    # Itero para calcular la anomalia excentrica
    for itera in range(15):
        M0 = E0 - e*np.sin(E0)
        E0 = E0 + (M-M0)/(1-e*np.cos(E0))

    # Anomalia verdadera
    true_anom = 2.0*np.arctan(np.sqrt((1 + e)/(1 - e))*np.tan(E0 / 2.0))

    # radius = (a*(1-e**2))/(1+e*cos(true_anom))
    radius = a*(1.0 - e*np.cos(E0))

    # No uso esta formula porque se pierde la informacion de los
    # cuadrantes, hay que usar atan2
    # theta = arctan(sin(true_anom + w)*cos(i)/cos(true_anom + w)) + O

    cos_vo = np.cos(true_anom + w)
    sin_vo = np.sin(true_anom + w)

    num = sin_vo * np.cos(i)
    # Uso arctan2 en lugar de cos/sin para no perder la
    # la informacion de los cuadrantes.
    theta = np.arctan2(num, cos_vo) + O

    rho = radius * (np.cos(true_anom + w)/np.cos(theta - O))

    if theta < 0:
        theta += 2*np.pi
    return rho, theta


def err_func(oe, rhos, thetas, times ):
    """
    Funcion de error O-C a minimizar

    theta en radianes
    """

    rho0 = np.array([])
    theta0 = np.array([])

    for time in times:
        rho, theta = binary_position(time, oe)
        rho0 = np.append(rho0, rho)
        theta0 = np.append(theta0, theta)

    return sum((rhos - rho0)**2 + (thetas - theta0)**2)


# Observaciones y elementos orbitales iniciales
obs = bl.read_binary_obs('GJ804')
obs['theta'] = np.radians(obs['theta'])

oe = bl.load_elements('GJ804')
oe0 = oe['P'], oe['T'], oe['e'], oe['a'], np.radians(oe['i']),\
    np.radians(oe['O']), np.radians(oe['w'])

bounds = [
    (13.0, 16.0),       # P
    (1990., 2010),      # T
    (0.3, 1.0),         # e
    (0.25, 0.6),        # a
    (np.pi, 2*np.pi),   # i  Mov retrogrado,  90 < i < 180
    (np.radians(10.), np.radians(350.)),      # OMEGA, O
    (np.radians(0.), np.radians(320.)),      # omega, w
]

res = basinhopping(err_func,
                   oe0,
                   niter=100,
                   stepsize=0.1,
                   minimizer_kwargs={'method': 'L-BFGS-B',
                                     'args': (obs['rho'],
                                              obs['theta'],
                                              obs['date']),
                                     'bounds': bounds
                                     }
                )

oe = {
    'WDS': ' GJ 804',
    'P': "{:.2f}".format(res['x'][0]),
    'T': "{:.2f}".format(res['x'][1]),
    'e': "{:.3f}".format(res['x'][2]),
    'a': "{:.3f}".format(res['x'][3]),
    'i': "{:.1f}".format(np.degrees(res['x'][4])),
    'O': "{:.1f}".format(np.degrees(res['x'][5])),
    'w': "{:.1f}".format(np.degrees(res['x'][6]))
}

print(res['message'])
print(bl.pprint(oe))

