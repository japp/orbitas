# -*- coding: utf-8 -*-

"""
Codigo para convertir los coeficientes de una elipse en representación canónica a los elementos orbitales geométricos por el metodo
de Kowalski.

Primero se debe ajustar los datos a una elipse y representarla de
for canónica:

Ax² + 2Hxy + By² + 2Gx + 2Fy + 1 = 0

Astrophysics Through Computation
Koberlein & Meisel, pp. 5


"""
import numpy as np
from scipy.optimize import fsolve

# Parametros de la elipse: A H B G F Z (Z debe ser +1)


# Elipse de Koberlein & Meisel
"""
------------------------
Semieje mayor:   2.489
Inclinacion:     57.2
Excentricidad:   -0.384
Omega:           102.1
omega:           50.5
------------------------
"""
#params = np.array([-0.65684, -0.251766/2.0, -0.242333, -0.552722/2.0, 0.1166782/2.0, +1])

#Gl473
#params = np.array([-0.57288251,  0.7584859/2.0,  -0.30767696, -0.01162644/2.0,  0.03013993/2.0,  0.0282147 ])/0.0282147

#Gl747
#params = np.array([0.08700319,  0.29248078/2.0,  0.95226477, -0.00478779/2.0,  0.00541168/2.0, -0.00500368])/-0.00500368


# GJ 804
#params = np.array([-0.42060675, -0.55834131/2.0, -0.70447197,  0.04170526/2.0, 0.11426085/2.0,  0.01641067])/0.01641067
#params = np.array([-0.54913466, -0.48391244/2.0, -0.67859635,  0.00497844/2.0, 0.05274851/2.0,  0.03129993])/0.03129993


# GJ 660
#params = np.array([-0.49547874,  0.44958611/2.0, -0.44151997,  0.40289302/2.0, -0.44023686/2.0, -0.03608297])/-0.03608297
#params = np.array([-0.75257952, -0.20733189/2.0, -0.52650233,  0.05736625/2.0, -0.01585286/2.0,  0.33149759])/0.33149759

# GJ 330
params = np.array([-0.39858254, -0.74397134, -0.52816025,  0.09302244,  0.00191116, -0.00534103])/-0.00534103


retrograde_orbit = True


def Omega_eq1(params):
    """Longitud del nodo ascendente"""
    A, H, B, G, F, Z = params

    result = np.arctan2(F**2 - G**2 + A - B, -2.0*(F*G - H))/2.0
    return result


def Omega_eq2(x, params):
    """Longitud del nodo ascendente"""
    A, H, B, G, F, Z = params

    result, = (F**2 - G**2 + A - B)*np.sin(2*x) + 2*(F*G - H)*np.cos(2*x)

    return result

# Omega, radianes
#Omega, = fsolve(Omega_eq2, 2.0, params)
Omega = Omega_eq1(params)

# Para que coincida N con Y+ y E con -X
#Omega =  Omega + np.pi/2.0
Omega =  Omega%np.pi

# Inclinacion y semilatus rectum
# Par de ecuaciones que relaciona i y p
# Esta funcion devuelve las ecuaciones en una
# lista para resolver con solve

def incl_semilatus(x, params):
    """
    NO LO ESTOY USANDO - REVISAR
    """

    A, H, B, G, F, Z = params

    eq1 = F*G - H + (np.sin(2.0*Omega) * np.tan(x[0])**2)/(2.0 * x[1] ** 2)

    eq2 = (F**2 - G**2 + A - B) - np.tan(x[0])**2/x[1]**2 - 2.0/x[1]**2

    return [eq1, eq2]


def tanip_eq(params):
    """
    Calcula la relación entre inclinacion (i) y el semilatus-rectum,
    para luego obtener i usando Omega

    tan(i)²/p²

    """
    A, H, B, G, F, Z = params

    tanip2 = (F**2 - G**2 + A - B)/(np.cos(2*Omega))
    #tanip3 = -2*(F*G-H)/np.sin(2*Omega)

    return abs(tanip2) #abs

# Valor de tan(i)²/p²
tanip = tanip_eq(params)


def semilatus2(params, Omega):
    A, H, B, G, F, Z = params

    semilat = (F**2 +G**2) - (A + B)/(F*G - H) * (1./np.sin(2*Omega))
    return 1./np.sqrt(semilat)


def semilatus(params, tanip):
    A, H, B, G, F, Z = params

    return np.sqrt(2.0/(F**2 + G**2 - A - B - tanip))

p  = semilatus(params, tanip)
#p = semilatus2(params, Omega)


def inclination(tanip, p):
    """Inclinacion de la orbita
    """
    return np.arctan(np.sqrt(tanip)*p)


def inclination2(params, Omega, p):
    """Inclinacion de la orbita"""
    A, H, B, G, F, Z = params

    return np.arctan(np.sqrt((2*p**2 *(H - F*G))/(np.sin(2*Omega))))


# Inclinacion, radianes
inc = inclination(tanip, p)

# Si el movimiento es retrógrado, la inclinación debe ser >90º
# si no lo es, debe ser <90
if retrograde_orbit and inc < np.pi:
    inc = np.pi - inc


def omega_eq(params, Omega, i):
    """Argumento del periastro

    Aitken, equation IV.31
    """

    A, H, B, G, F, Z = params

    num = (F*np.cos(Omega) - G*np.sin(Omega)) * np.cos(i)
    denom = F*np.sin(Omega) + G*np.cos(Omega)
    omega = np.arctan2(num, denom)

    return omega


omega = omega_eq(params, Omega, inc)


def ecc_eq(params, i, omega, p):
    """Excentricidad de la orbita

    Si i>90º la eccentricidad será negativa
    """

    A, H, B, G, F, Z = params

    return -((G*np.sin(Omega) + F*np.cos(Omega))*p)/np.cos(omega)


def ecc_eq2(params, i, omega, p):
    """Inclinacion de la orbita

    Si i>90º la eccentricidad será negativa
    """

    A, H, B, G, F, Z = params

    return ((G*np.sin(Omega) - F*np.cos(Omega))*p*np.cos(i))/np.sin(omega)


# Excentricidad
e = ecc_eq2(params, inc, omega, p)

def a_eq(e, p):
    """Semieje mayor de la orbita en segundos de arco"""

    return p/(1 - e**2)

# Semieje mayor
a = a_eq(e, p)


def true_anomaly_eq(theta, Omega, inc, omega):
    """Anomalia verdadera"""

    return np.arctan2(np.tan(theta - Omega), np.cos(inc)) - omega


def ecc_anomaly_eq(v, e):
    """Anomalia excentrica"""
    en = np.sqrt((1.0 + e)/(1.0 - e))
    return 2.0*np.arctan2(np.tan(v/2.0), en)

result = """
------------------------
Semieje mayor:   {a:.3f}
Inclinacion:     {Inc:.1f}
Excentricidad:   {eccen:.3f}
Omega:           {Omega:.1f}
omega:           {omega:.1f}
------------------------
"""

result = result.format(Omega=np.degrees(Omega),
                       Inc=np.degrees(inc),
                       omega=np.degrees(omega),
                       eccen=e,
                       a=a)

print(result)
