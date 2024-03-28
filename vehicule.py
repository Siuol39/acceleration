# -*- coding: utf-8 -*-

import json

class Vehicule:
    
    _rho = 1.2
    
    def __init__(self, filename):
        with open(filename, 'r') as f:
            data = json.load(f)
        self.nom = data["nom"]
        self.masse_vide = data["masse_vide"] * 1000
        self.charge_utile = data["charge_utile"] * 1000
        self.puissance = data["puissance"] * 1000
        self.ratio_puissance = 0.7
        self.traction_max = data["traction_max"] * 1000
        self.hauteur = data["hauteur"]
        self.largeur = data["largeur"]
        self.longueur = data["longueur"]
        self.cx = data["cx"]
        self.crr = data["crr"]
        self.deceleration = data["deceleration"]
        self.v_traction_max = data["v_traction_max"] /3.6
    
    def surface(self):
        return self.hauteur * self.largeur
    
    def masse(self):
        return self.masse_vide + self.charge_utile
    
    def resistance(self, vitesse, pente):
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
        ALPHA = 0.5 * self._rho * self.surface() * self.cx
        R_aero = ALPHA * vitesse ** 2
        R_elevation = pente * self.masse() * 9.81 # m * g * v verticale / v
        R_roulement = self.crr * self.masse() * 9.81
        
        return R_aero + R_elevation + R_roulement
    
    def _force_traction_max(self, vitesse):
        t_max = self.traction_max
        t_v0 = t_max / 10
        v_t_max = self.v_traction_max
        
        if vitesse > v_t_max:
            return t_max
        else:
            return vitesse * (t_max - t_v0) / v_t_max + t_v0
    
    def force_traction(self, vitesse):
        if vitesse == 0:
            return self._force_traction_max(vitesse)
        
        return min(self.ratio_puissance * self.puissance / vitesse,
                   self._force_traction_max(vitesse))
        
        





























