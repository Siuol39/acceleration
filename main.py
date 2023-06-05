# -*- coding: utf-8 -*-
"""
Created on Thu May 25 10:05:32 2023

@author: HYENNE
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from math import sin, pi

KILO = 1000
TONNE = 1000
WATT = 1

# MASSE = 180 * TONNE # Regiolis 4 caisses
MASSE = 275 * TONNE # Regiolis 6 caisses

# PUISSANCE = 0.692 * 1800 * KILO * WATT # Regiolis 4 caisses
PUISSANCE = 0.655 * 2700 * KILO * WATT # Regiolis 6 caisses

# PUISSANCE_FREINAGE = -3500 * KILO * WATT # Regiolis 4 caisses
PUISSANCE_FREINAGE = -5000 * KILO * WATT # Regiolis 6 caisses

RHO = 1.2
S = 13
CX = 0.6
ALPHA = 0.5 * RHO * S * CX # 1/2 rho S Cx Regiolis

CRR = 0.003

VITESSE_CIBLE = 160/3.6
ECART = 1/3.6

DATA = pd.read_csv("data.csv", sep=";", header=0, decimal=",")

def get_index(pk):
    """
    Donne l'index i tel que le point kilométrique pk est entre la ligne i - 1 et i.
    Retourne None si pk trop grand.

    Parameters
    ----------
    pk : float
        Point kilométrique en m.

    Returns
    -------
    i : int

    """
    pk_km = pk / 1000
    i = 0
    while i < DATA.count()["PK"] and DATA.iloc[i]["PK"] <= pk_km:
        i += 1
    if i >= DATA.count()["PK"]:
        return None
    return i

def pente(pk, croissant = True):
    """
    Calcule la pente au point spécifié.

    Parameters
    ----------
    pk : flot
        Point kilométrique en m
    croissant : bool
        vrai si la ligne est parcourue dans le sens des PK croissants

    Returns
    -------
    float
        Pente en m/m

    """
    i = get_index(pk)
    if i == None:
        return 0
    longueur = DATA.iloc[i]["PK"] - DATA.iloc[i - 1]["PK"]
    deniv = (DATA.iloc[i]["alt"] - DATA.iloc[i - 1]["alt"]) / 1000
    if croissant:
        return deniv / longueur
    else:
        return -deniv / longueur

def altitude(pk):
    """
    Calcule l'altitude au point donné.

    Parameters
    ----------
    pk : float
        PK en m.

    Returns
    -------
    float
        Altitude en m.

    """
    i = get_index(pk)
    if i == None:
        return 0
    a = DATA.iloc[i - 1]["PK"]
    b = DATA.iloc[i]["PK"]
    aa = DATA.iloc[i - 1]["alt"]
    ba = DATA.iloc[i]["alt"]
    return (pk / 1000 - a) * (ba - aa) / (b - a) + aa

def resistance(vitesse, pente):
    """
    Calcule la résistance subie à une vitesse et pour une pente donnée.
    Le résultat peut être négatif si la pente est négative.

    Parameters
    ----------
    vitesse : float
        Vitesse en m/s.
    pente : float
        pente en mm/m.

    Returns
    -------
    float
        Résistance en N

    """
    v = vitesse
    R_aero = ALPHA * v ** 2
    R_elevation = pente * MASSE * 9.81 # m * g * v verticale / v
    R_roulement = CRR * MASSE * 9.81
    
    return R_aero + R_elevation + R_roulement

resistance_pk_croissant = lambda vitesse, pk : resistance(vitesse, pente(pk))
resistance_pk_decroissant = \
    lambda vitesse, pk : resistance(vitesse, pente(pk, False))
    
def force_traction(vitesse, vitesse_cible = VITESSE_CIBLE):
    """
    Force de traction en fonction de la vitesse.
    
    En dessous de VITESSE_CIBLE - ECART : pleine puissance
    Entre VITESSE_CIBLE - ECART et VITESSE_CIBLE : puissance décroissante vers
        0 en suivant une sinusoïde
    Entre VITESSE_CIBLE et VITESSE_CIBLE + ECART : puissance croissante depuis
        0 en suivant une sinusoïde
    Au-delà de VITESSE_CIBLE + ECART : pleine puissance de freinage
    
    Vitesse traitée en valeur absolue

    Parameters
    ----------
    vitesse : float
        m/s

    Returns
    -------
    float
        N

    """
    vitesse = abs(vitesse)
    if vitesse <= vitesse_cible - ECART:
        p = PUISSANCE
    elif vitesse <= vitesse_cible:
        #p = PUISSANCE * sin(pi / 2 * (vitesse - VITESSE_CIBLE) / ECART - pi)
        p = (vitesse_cible - vitesse) * PUISSANCE / ECART
    elif vitesse < vitesse_cible + ECART:
        #p = PUISSANCE_FREINAGE * sin(pi / 2 * (vitesse - VITESSE_CIBLE) / ECART)
         p = (vitesse - vitesse_cible) * PUISSANCE_FREINAGE / ECART
    else:
        p = PUISSANCE_FREINAGE
    return p / vitesse

def f_integree_plat(t, y, v_cible):
    return [force_traction(y[0], v_cible) / MASSE - resistance(y[0], 0) / MASSE, y[0]]

def f_integree_croissant(t, y):
    return [force_traction(y[0]) / MASSE - resistance_pk_croissant(y[0], y[1]) / MASSE, y[0]]

def f_integree_decroissant(t, y):
    return [-force_traction(y[0]) / MASSE + resistance_pk_decroissant(y[0], y[1]) / MASSE, y[0]]

t0 = 0
tmax = 300
range_dist = np.linspace(0, 19000, 100)


# for v in np.linspace(50, 160, 111):
#     v_ini = v / 3.6
#     v_cible = 1 / 3.6
#     sol_frein = solve_ivp(f_integree_plat, [t0, tmax], [v_ini, 0],
#                     t_eval = np.linspace(t0, tmax, 200),
#                     args = [v_cible],
#                     events = lambda t, y, v_c : y[0] - 1.02 * v_c)
    
#     print(v_ini * 3.6, sol_frein.y_events[0][0][1], sol_frein.t_events[0][0])
    
for v in np.linspace(50, 160, 111):
    v_ini = 1e-3
    v_cible = v / 3.6
    sol_acc = solve_ivp(f_integree_plat, [t0, tmax], [v_ini, 0],
                    t_eval = np.linspace(t0, tmax, 200),
                    args = [v_cible],
                    events = lambda t, y, v_c : y[0] - 0.99 * v_c)
    
    print(v_cible * 3.6, sol_acc.y_events[0][0][1], sol_acc.t_events[0][0])
    

# sol = solve_ivp(f_integree_plat, [t0, tmax], [1e-3, 0],
#                 t_eval = np.linspace(t0, tmax, 200),
#                 args = [VITESSE_CIBLE])

# fig = plt.figure(0, dpi=300)
# plt.title("Vitesse du train (accélération à plat)")
# ax = plt.gca()
# ax.plot(sol.y[1] / 1000, 3.6 * sol.y[0])
# ax.set_xlabel("distance (km)")
# ax.set_ylabel("vitesse (km/h)")
# # plt.savefig(".png")
# plt.show()

# sol = solve_ivp(f_integree_croissant, [t0, tmax], [0.01, 0],
#                 t_eval = np.linspace(t0, tmax, 100))
    
# fig = plt.figure(1, dpi=300)
# plt.title("Vitesse du train (Annecy -> La Roche)")
# ax1 = plt.gca()
# ax2 = ax1.twinx()
# ax1.plot(sol.y[1] / 1000, 3.6 * sol.y[0])
# ax2.plot(range_dist / 1000, [altitude(pk) for pk in range_dist], 'r')
# ax1.set_xlabel("distance (km)")
# ax1.set_ylabel("vitesse (km/h)")
# ax2.set_ylabel("altitude (m)")
# plt.savefig("vitesse A-LR.png")
# plt.show()

# sol = solve_ivp(f_integree_decroissant, [t0, tmax], [-0.01, 19400],
#                 t_eval = np.linspace(t0, tmax, 100))
    
# fig = plt.figure(2, dpi=300)
# plt.title("Vitesse du train (La Roche -> Annecy)")
# ax1 = plt.gca()
# ax2 = ax1.twinx()
# ax1.plot(sol.y[1] / 1000, -3.6 * sol.y[0])
# ax2.plot(range_dist / 1000, [altitude(pk) for pk in range_dist], 'r')
# ax1.set_xlabel("distance (km)")
# ax1.set_ylabel("vitesse (km/h)")
# ax2.set_ylabel("altitude (m)")
# plt.savefig("vitesse LR-A.png")
# plt.show()








































