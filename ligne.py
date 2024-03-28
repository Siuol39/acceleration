# -*- coding: utf-8 -*-

import json
import numpy as np
from scipy.integrate import solve_ivp

class Ligne:
    
    _g = 9.81
    
    def __init__(self, filename):
        with open(filename, 'r') as f:
            data = json.load(f)
        self.points_altitude = data["points_altitude"]
        self.points_vitesse = data["points_vitesse"]
        self.points_arret = data["points_arret"]
    
    def v_max(self):
        return max([point["vitesse"] for point in self.points_vitesse])
    
    def longueur(self):
        return self.points_altitude[-1]["pk"]
    
    def _get_index(self, pk):
        """
        Donne l'index i tel que le point kilométrique pk est entre la ligne i - 1 et i.
        Retourne None si pk trop grand.
    
        Parameters
        ----------
        pk : float
            Point kilométrique en km.
    
        Returns
        -------
        i : int
    
        """
        i = 0
        while i < len(self.points_altitude) and self.points_altitude[i]["pk"] <= pk:
            i += 1
        if i >= len(self.points_altitude):
            return None
        return i
    
    def pente(self, pk, croissant = True):
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
        i = self._get_index(pk)
        if i == None:
            return 0
        longueur = self.points_altitude[i]["pk"] \
            - self.points_altitude[i-1]["pk"]
        deniv = (self.points_altitude[i]["altitude"] \
            - self.points_altitude[i-1]["altitude"]) / 1000
        if croissant:
            return deniv / longueur
        else:
            return -deniv / longueur
    
    def altitude(self, pk):
        """
        Calcule l'altitude au point donné.
    
        Parameters
        ----------
        pk : float
            PK en km.
    
        Returns
        -------
        float
            Altitude en m.
    
        """
        i = self.get_index(pk)
        if i == None:
            return 0
        a = self.points_altitude[i-1]["pk"]
        b = self.points_altitude[i]["pk"]
        aa = self.points_altitude[i-1]["altitude"]
        ba = self.points_altitude[i]["altitude"]
        return (pk - a) * (ba - aa) / (b - a) + aa


    def acceleration(self, vehicule, pk, vitesse, croissant):
        return (vehicule.force_traction(vitesse) \
            - vehicule.resistance(vitesse, self.pente(pk, croissant))) \
            / vehicule.masse()
    
    def deceleration(self, vehicule, pk, croissant):
        return vehicule.deceleration + self._g * self.pente(pk, croissant)
    
    def profil_vitesse_acceleration(self, vehicule, pk_ini, v_ini, croissant):
        if croissant:
            sgn = 1
        else:
            sgn = -1
        
        def f_integree(t, y):
            v, x = y
            return [sgn * self.acceleration(vehicule, pk_ini + sgn * x, v, croissant), v]
        
        t0 = 0
        tmax = 10000
        sol = solve_ivp(f_integree, [t0, tmax], [v_ini, pk_ini],
                    t_eval = np.linspace(t0, tmax, 200))
        
        return sol


    def profil_vitesse_deceleration(self, vehicule, pk_fin, v_fin, croissant):
        pass































