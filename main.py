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

vehicules = {
        "citadis": {
            "nom": "Alstom Citadis",
            "masse": 75 * TONNE,
            "puissance": 0.7 * 960 * KILO * WATT,
            "puissance f": -1500 * KILO * WATT,
            "surface": 2.5 * 3,
            "cx": 0.6,
            "crr": 0.003
            },
        "citadis+": {
            "nom": "Alstom Citadis plus pouissant",
            "masse": 75 * TONNE,
            "puissance": 0.7 * 1500 * KILO * WATT,
            "puissance f": -2000 * KILO * WATT,
            "surface": 2.5 * 3,
            "cx": 0.6,
            "crr": 0.003
            },
        "regiolis4r": {
            "nom": "Regiolis 4 caisses « régional »",
            "masse": 145.6 * TONNE,
            "puissance": 0.692 * 1800 * KILO * WATT,
            "puissance f": -3500 * KILO * WATT,
            "surface": 2.85 * 4.29,
            "cx": 0.6,
            "crr": 0.003
            },
        "regiolis6r": {
            "nom": "Regiolis 6 caisses « régional »",
            "masse": 244.6 * TONNE,
            "puissance": 0.655 * 2700 * KILO * WATT,
            "puissance f": -3500 * KILO * WATT,
            "surface": 2.85 * 4.29,
            "cx": 0.6,
            "crr": 0.003
            },
        "regio2nTETnormandie": {
            "nom": "Regio2N « TET Normandie »",
            "masse": 363.4 * TONNE,
            "puissance": 0.5 * 3200 * KILO * WATT,
            "puissance f": -4000 * KILO * WATT,
            "surface": 3.05 * 4.32,
            "cx": 0.6,
            "crr": 0.003
            },
        "regio2n_8_3200": {
            "nom": "Regio2N version L (8 caisses) 3200 kW",
            "masse": 328 * TONNE,
            "puissance": 0.7 * 3200 * KILO * WATT,
            "puissance f": -4000 * KILO * WATT,
            "surface": 3.05 * 4.32,
            "cx": 0.6,
            "crr": 0.003
            },
        "flirt_521_523_CFF": {
            "nom": "Flirt 521-523 CFF",
            "masse": 157 * TONNE,
            "puissance": 0.7 * 2000 * KILO * WATT,
            "puissance f": -3500 * KILO * WATT,
            "surface": 2.88 * 4.15,
            "cx": 0.6,
            "crr": 0.003
            },
        "fretBB26000": {
            "nom": "Fret 3000 t BB 26000",
            "masse": 3000 * TONNE,
            "puissance": 0.95 * 5600 * KILO * WATT,
            "puissance f": -5000 * KILO * WATT,
            "surface": 2.93 * 4.27,
            "cx": 0.6,
            "crr": 0.003
            },
        "fretBB27000": {
            "nom": "Fret 2000 t BB 27000",
            "masse": 2050 * TONNE,
            "puissance": 0.7 * 4200 * KILO * WATT,
            "puissance f": -5000 * KILO * WATT,
            "surface": 2.93 * 4.27,
            "cx": 0.6,
            "crr": 0.003
            },
        }

vehicule = vehicules["flirt_521_523_CFF"]


RHO = 1.2

ALPHA = 0.5 * RHO * vehicule["surface"] * vehicule["cx"]

VITESSE_INI = 160/3.6
VITESSE_CIBLE = 160/3.6
ECART = 1/3.6
PENTE = 5.35/1000

# DATA = pd.read_csv("profil autoroute La Roche Annecy.csv", sep=";", header=0, decimal=",")
DATA = pd.read_csv("A48 Bourgouin-Voreppe.csv", sep=";", header=0, decimal=",")

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
        pente en m/m.

    Returns
    -------
    float
        Résistance en N

    """
    v = vitesse
    R_aero = ALPHA * v ** 2
    R_elevation = pente * vehicule["masse"] * 9.81 # m * g * v verticale / v
    R_roulement = vehicule["crr"] * vehicule["masse"] * 9.81
    
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
        p = vehicule["puissance"]
    elif vitesse <= vitesse_cible:
        #p = PUISSANCE * sin(pi / 2 * (vitesse - VITESSE_CIBLE) / ECART - pi)
        p = (vitesse_cible - vitesse) * vehicule["puissance"] / ECART
    elif vitesse < vitesse_cible + ECART:
        #p = PUISSANCE_FREINAGE * sin(pi / 2 * (vitesse - VITESSE_CIBLE) / ECART)
         p = (vitesse - vitesse_cible) * vehicule["puissance f"] / ECART
    else:
        p = vehicule["puissance f"]
    return p / vitesse

def f_integree_pente(t, y, v_cible, pente):
    return [force_traction(y[0], v_cible) / vehicule["masse"] - resistance(y[0], pente) / vehicule["masse"], y[0]]

def f_integree_croissant(t, y):
    return [force_traction(y[0]) / vehicule["masse"] - resistance_pk_croissant(y[0], y[1]) / vehicule["masse"], y[0]]

def f_integree_decroissant(t, y):
    return [-force_traction(y[0]) / vehicule["masse"] + resistance_pk_decroissant(y[0], y[1]) / vehicule["masse"], y[0]]

t0 = 0
tmax = 1200

# for v in np.linspace(20, 100, 81):
#     v_ini = v / 3.6
#     v_cible = 1 / 3.6
#     sol_frein = solve_ivp(f_integree_pente, [t0, tmax], [v_ini, 0],
#                     t_eval = np.linspace(t0, tmax, 200),
#                     args = [v_cible, 0/1000],
#                     events = lambda t, y, v_c, p : y[0] - 1.02 * v_c)
    
#     print("{0:0.0f}".format(v_ini * 3.6),
#           "{0:0.5f}".format(sol_frein.y_events[0][0][1]).replace('.', ","),
#           "{0:0.5f}".format(sol_frein.t_events[0][0]).replace('.', ","))
    
# for v in np.linspace(20, 70, 51):
#     v_ini = 1e-3
#     v_cible = v / 3.6
#     sol_acc = solve_ivp(f_integree_pente, [t0, tmax], [v_ini, 0],
#                     t_eval = np.linspace(t0, tmax, 200),
#                     args = [v_cible, 0/1000],
#                     events = lambda t, y, v_c, p : y[0] - 0.98 * v_c)
    
#     print("{0:0.0f}".format(v_cible * 3.6),
#           "{0:0.5f}".format(sol_acc.y_events[0][0][1]).replace('.', ","),
#           "{0:0.5f}".format(sol_acc.t_events[0][0]).replace('.', ","))

# sol = solve_ivp(f_integree_pente, [t0, tmax], [VITESSE_INI, 0],
#                 t_eval = np.linspace(t0, tmax, 200),
#                 args = [VITESSE_CIBLE, PENTE],
#                 events = [lambda t, y, v_c, p : y[0] - 50/3.6])

# fig = plt.figure(0, dpi=80)
# plt.title("Vitesse du {0} (pente {1:0.0f} mm/m)".format(vehicule["nom"], 1000 * PENTE))
# ax = plt.gca()
# ax.plot(sol.y[1] / 1000, 3.6 * sol.y[0])
# ax.set_xlabel("distance (km)")
# ax.set_ylabel("vitesse (km/h)")
# # plt.savefig(".png")
# plt.show()

# print(sol.y_events[0][0][0], sol.y_events[0][0][1], sol.t_events[0][0])
# print(sol.y_events[0][0][0] / sol.t_events[0][0])

range_dist = np.linspace(0, 45000, 100)    


# sol = solve_ivp(f_integree_croissant, [t0, tmax], [VITESSE_INI, 0],
#                 t_eval = np.linspace(t0, tmax, 100))
    
# fig = plt.figure(1, dpi=300)
# plt.title("Vitesse du train (aller)")
# ax1 = plt.gca()
# ax2 = ax1.twinx()
# ax1.plot(sol.y[1] / 1000, 3.6 * sol.y[0])
# ax2.plot(range_dist / 1000, [altitude(pk) for pk in range_dist], 'r')
# ax1.set_xlabel("distance (km)")
# ax1.set_ylabel("vitesse (km/h)")
# ax2.set_ylabel("altitude (m)")
# plt.savefig("vitesse A-LR.png")
# plt.show()

sol = solve_ivp(f_integree_decroissant, [t0, tmax], [-VITESSE_INI, 42030],
                t_eval = np.linspace(t0, tmax, 100))
    
fig = plt.figure(2, dpi=300)
plt.title("Vitesse du train (retour)")
ax1 = plt.gca()
ax2 = ax1.twinx()
ax1.plot(sol.y[1] / 1000, -3.6 * sol.y[0])
ax2.plot(range_dist / 1000, [altitude(pk) for pk in range_dist], 'r')
ax1.set_xlabel("distance (km)")
ax1.set_ylabel("vitesse (km/h)")
ax2.set_ylabel("altitude (m)")
plt.savefig("vitesse LR-A.png")
plt.show()








































