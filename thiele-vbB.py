# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 19:55:33 2012

@author: japp
"""

import numpy as np
from scipy.optimize import broyden1


# --- Datos de entrada ------
# Tres medidas separadas que cumplan bien
# la ley de las areas
# ---------------------------


# Estimacion del doble de la velocidad areolar
# Signo menos si es retrogrado
c =+0.04348
        
times = np.array([1986.4475, 2006.4488, 2015.494])    # Argyle

t1, t2, t3 = times

rhos = np.array([0.650, 0.891, 0.706])   # Argyle

rho1, rho2, rho3 = rhos

thetas = np.radians(np.array([156.2, 188.9, 271.0]))    # Argyle

theta1, theta2, theta3 = thetas


def areas_constants(rhos, thetas, times):

    dt = np.array([])
    for i in range(len(times)-1):
        dt = np.append(dt, times[i+1] - times[i])

    dS = np.array([])
    for i in range(len(times)-1):
        dS = np.append(dS,0.5*doubleArea(rhos[i], rhos[i+1], thetas[i], thetas[i+1]))

    return dt, dS


def neg_to_pos_angle(angle, degrees=False):
    """
    Cambia un angulo negativo (es decir, contando hacia W),
    si lo es, a positivo

    angle : float
        Angulo a camibiar, en radianes
    degrees : bool
        Si True, acepta grados en lugar de radianes

    return:
        angle : float
        Angulo cambiado a positivo en radianes (grados, si degrees=True)
    """

    if degrees:
        circunference = 2*360.0
    else:
        circunference = 2*np.pi

    if angle < 0:
        angle = circunference + angle

    return angle


def mu_p_q(L12, L23, L13):
    """
    Calcula los valores de mu, p y q por fuerza bruta
    """

    def equations(params):
        mu, p, q = params

        eqs = [mu*L12 - p + np.sin(p),
               mu*L23 - q + np.sin(q),
               mu*L13 - p - q + np.sin(p + q)
               ]
        return eqs

    res_mu0, res_p0, res_q0 = np.array([0.01, 0.01, 0.01])

    for mu0 in np.arange(0.01, 5, 0.01):
        try:
            for p0 in np.arange(0.01, 10, 0.01):
                for q0 in np.arange(0.01, 10, 0.01):
    
                    mu1, p1, q1 = broyden1(equations, (mu0, p0, q0))
                    res_mu, res_p, res_q = equations((mu1, p1, q1))
                    
                    if (abs(np.array([res_mu, res_p, res_q])) <
                        abs(np.array([res_mu0, res_p0, res_q0]))).all()\
                        and mu1 > 0.0001:
                        res_mu0, res_p0, res_q0 = res_mu, res_p, res_q
                    
                        mu, p, q = mu1, p1, q1
        except:
            continue

    try:
        return np.array([mu, p, q])
    except:
        print(mu1, p1, q1)
        raise ValueError("No se encontraron valores para mu, p, q")


x1, x2, x3 = rhos*np.cos(thetas)
y1, y2, y3 = rhos*np.sin(thetas)


def doubleArea(rho1, rho2, theta1, theta2):
    """
    Doble area del triangulo formado por dos posiciones (rho, theta)
    """
    return rho1*rho2*np.sin(theta2 - theta1)


delta12 = doubleArea(rho1, rho2, theta1, theta2)
delta23 = doubleArea(rho2, rho3, theta2, theta3)
delta13 = doubleArea(rho1, rho3, theta1, theta3)


print("delta12={delta12:.3f}  delta23={delta23:.3f}  delta13={delta13:.3f}".\
      format(delta12=delta12, delta23=delta23, delta13=delta13))


L12 = t2 - t1 - delta12/c
L23 = t3 - t2 - delta23/c
L13 = t3 - t1 - delta13/c


print("L12={L12:.3f}  L23={L23:.3f}  L13={L13:.3f}".\
    format(L12=L12, L23=L23, L13=L13))


mu, p, q =  mu_p_q(L12, L23, L13)

print("mu={mu:.3f}  p={p:.3f}  q={q:.3f}".format(mu=mu, p=p, q=q))


esinE2 = (delta23*np.sin(p) - delta12*np.sin(q))/(delta23 + delta12 - delta13)
ecosE2 = (delta23*np.cos(p) + delta12*np.cos(q) - delta13)/(delta23 + delta12 - delta13)

E2 = np.arctan2(esinE2, ecosE2)

E2 = neg_to_pos_angle(E2)

e1 = esinE2/np.sin(E2)
e2 = ecosE2/np.cos(E2)

# Excentricidad
e = (e1 + e2)/2

print("Ecentricidad: {:.3f}".format(e))

E1 = E2 - p
E3 = E2 + q

M1 = E1 - e*np.sin(E1)
M2 = E2 - e*np.sin(E2)
M3 = E3 - e*np.sin(E3)

T1 = t1 - M1/mu
T2 = t2 - M2/mu
T3 = t3 - M3/mu

# Tiempo de paso por el periastro
T = np.mean([T1, T2, T3])

print("Tiempo de paso por el periastro (T): {:.3f}".format(T))


Es = np.array([E1, E2, E3])

# Coordenadas reducidas de la compañera
X1, X2, X3 = np.cos(Es) - e
Y1, Y2, Y3 = np.sqrt(1 - e**2)*np.sin(Es)

a1 = np.array([[X1, Y1], [X2, Y2]])
b1 = np.array([x1, x2])
A, F = np.linalg.solve(a1, b1)


a2 = np.array([[X1, Y1], [X2, Y2]])
b2 = np.array([y1, y2])
B, G = np.linalg.solve(a2, b2)
# Hay que repetirlo para (X1,Y1)+(X3,Y3) t sacar la media de A, B, F, G

u = (A**2 + B**2 + F**2 + G**2)/2
v = A*G - B*F

a = np.sqrt(u + np.sqrt((u+v)*(u-v)))

print("Semieje mayor (a): {:.3f}".format(a))


i = np.arccos(v/a**2)

i = np.degrees(neg_to_pos_angle(i))

print("Inclinacion (i): {:.1f}".format(i))

w_plus_O = np.arctan2((B - F), (A + G))
w_minus_O = np.arctan2((B + F), (A - G))

a3 = np.array([[+1.0, +1.0], [1.0, -1.0]])
b3 = np.array([w_plus_O, w_minus_O])
O, w = np.linalg.solve(a3, b3)

O = neg_to_pos_angle(O)
w = neg_to_pos_angle(w)

print("Omega: (O) {:.1f}".format(np.degrees(O)))
print("omega (w): {:.1f}".format(np.degrees(w)))
