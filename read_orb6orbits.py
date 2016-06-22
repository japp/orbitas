# -*- coding: utf-8 -*-
"""
Created on Sun Sep  6 20:29:25 2015

@author: japp

Funciones para leer los elementos orbitales del
6th Catalog of Orbits.

Este script es de pruebas, la version estable 
está en binarylib


"""

from astropy.table import Table
from astropy.time import Time


def get_oe(name):
    """
    Lee los elementos orbitales del Sixth Catalog of Orbits of Visual Binary Stars
    que estan guardados en disco. La versión en disco es la misma que en la web
    pero con la cabecera cambiada, que es mas clara:
    
    Nueva cabecera
    RA|DEC|WDS|DD|ADS|HD|HIP|V1|V_flag|V2|V2_flag|P|P_flag|P_e|a|a_flag|a_e|i|i_e|
    O|O_flag|O_e|T|T_flag|T_e|e|e_e|w|w_e|EQNX|LAST|G|N|REF|PNGFILE
    
    Algunos parametros como T, P o a pueden tener distintas unidades y por eso
    antes se conververten a años (P, T) o arcsec para el semieje mayor (a)
    El flag de cada elemento suele indicar las unidades en que estan y
    se usa para identificar.
    
    A veces no se error en un parámetro y en ese caso el campo esta enmascarado (mask)
    si es asi no se hace conversión de unides de P, T, a para que no de error.
    
    Parameters
    ----------
    name : str or tuple
           Nombre de la estrellas en WDS '19512-7248' en string o
           de otro catalogo como tupla ('HIP', '97690')
        
    
    Returns
    -------
    oe : dict
        Diccionario con los elementos orbitales de la estrella
        
        P, T en años
        i, O y w en grados
        a en arcsec

    Examples
    --------
    oe = get_oe('19512-7248')  # Busca por nombre WDS
    oe = get_oe( ('HIP', '97690') )   # Busca por nombre HIP  (DD, ADS, HIP, HD)  
    
    
    Reference
    ---------
    
    Sixth Catalog webpage
    http://ad.usno.navy.mil/wds/orb6.html

    Catalogo orginal limitado por "|"
    http://ad.usno.navy.mil/wds/orb6/orb6orbits.sql    
    
    Catalogo original en txt con columnas fijas
    http://ad.usno.navy.mil/wds/orb6/orb6format.txt
    
    """

    # Tabla de elementos orbitales
    table = read_orb6_data()    
    
    # Estrella a buscar
    # Por defecto busca el nombre WDS (tipo 00026-0829)
    # Si se da una tupla, busca por otro catalogo (DD, ADS, HIP, HD) 
    if type(name) == tuple:   
        oe = table[table[name[0]] == name[1]]
    else:
        oe = table[table['WDS']] == name
    
    # De Table a dict
    oe_dict = dict(zip(oe[0].colnames, oe[0].data))
   
    # ----- Cambio de unidades de P a años  -----
   
    # Posibles unidades de P
    P_units2year = {
        'm' : 365.242199*24*3600,  # minutes
        'd' : 365.242199, # days
        'y' : 1,          # year  
        'c' : 1/100.0    # centuries
    }
    
    # Unidades en las que esta P
    # Unidades en las que esta P
    P_u = oe_dict['P_flag']
    
    oe_dict['P'] = oe_dict['P']/P_units2year[P_u]
    oe_dict['P_e'] = oe_dict['P_e']/P_units2year[P_u]
    
    
    # -- Cambio de unidades de T a años    
    
    oe_dict['T_e'] = float(oe_dict['T_e'])
    
    # Julian date  (-2 400 000 days)
    if oe_dict['T_flag'] == 'd':
        t = Time(oe_dict['T'] + 2400000, format='jd')
        if oe_dict['T_e']:
            t_err = Time(oe_dict['T_e'] + 2400000, format='jd')
    # MJD
    elif oe_dict['T_flag'] == 'm':
        t = Time(oe_dict['T'], format='mjd')
        if oe_dict['T_e']:
            t_err = Time(oe_dict['T_e'], format='mjd')
    # fractional Besselian year
    elif oe_dict['T_flag'] == 'y':
        t = Time(oe_dict['T'], format='byear')
        if oe_dict['T_e']:
            t_err = Time(oe_dict['T_e'], format='byear')  
   
    oe_dict['T'] = t.byear
    
    if oe_dict['T_e']:
        oe_dict['T_e'] = t_err.byear
    
    # ----- Cambio de unidades de 'a' a arcseconds   ----- 
    
    # Posibles unidades de a
    a_units2as = {
        'a' : 1,  # arcseconds
        'm' : 1000.0, # milliarcseconds  (mas)
        'u' : 1.e+6,          #  microarcseconds (uas - not yet used)  
    }
    
    # Unidades en las que esta a (normalmente arcsec)
    a_u = oe_dict['a_flag']
    
    oe_dict['a'] = oe_dict['a']/a_units2as[a_u]
    
    if not oe_dict['a_e'].mask:
        oe_dict['a_e'] = oe_dict['a_e']/a_units2as[a_u]

    return oe_dict


def read_orb6_data():
    """
    Lee el fichero de elementos orbitales del 6th orbit catalog
    
    Parameters
    ----------
    None
        
    
    Returns
    -------
    table : astropy.Table type object


    """
    binary_data_path = "/home/japp/codigo/py/orbitas/binaries_data/"
    filename = binary_data_path + "orb6orbits.sql"
    table = Table.read(filename, format="ascii", delimiter='|', header_start=1, data_start=2)
    
    return table


def export_oe(wds_name, alias=None):
    
    binary_data_path = "/home/japp/codigo/py/orbitas/binaries_data/"
    
    oe = get_oe(wds_name)
    if alias:
        oe['WDS'].data[0] = alias
    
    content = """
        name        {WDS}
        P           {P}                 
        T           {T}              
        e           {e}                 
        a           {a}                 
        i           {i}                 
        O           {O}             
        w           {w}
        equinox     {EQNX}
        """.format(**oe)
    fout = open(binary_data_path + "{}.oe.txt".format(oe['WDS'].data[0]), 'w')
    fout.write(content)
    fout.close()

