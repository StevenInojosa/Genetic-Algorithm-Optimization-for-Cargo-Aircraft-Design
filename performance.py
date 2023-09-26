#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This module contains functions to compute the Performance of an Aircraft.

(c) Copyright 2023, Steven Inojosa

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import numpy as np
import math
from sympy import symbols, solve, diff

from sklearn.linear_model import LinearRegression

class Engine:
    def __init__(self, name, a, b, c, d, density, altitude):        
        """
        Initialize a new instance of the class to calculate parameters a, b, c, d, and density for a given altitude.

        Parameters:
        -----------
        name : str
            The name of the instance.
        a : list
            The available power coefficient values for the corresponding altitudes.
        b : list
            The available power coefficient values for the corresponding altitudes.
        c : list
            The available power coefficient values for the corresponding altitudes.
        d : list
            The available power coefficient values for the corresponding altitudes.
        density : list
            The density [kg/m^3] for the corresponding altitudes in International Standard Atmosphere (ISA).
        altitude : float
            The altitude [m] for which the parameters are to be calculated in International Standard Atmosphere (ISA).
        """
        self.__name = name
        
        a = np.array(a).reshape(-1, 1)
        b = np.array(b).reshape(-1, 1)
        c = np.array(c).reshape(-1, 1)
        d = np.array(d).reshape(-1, 1)
        density = np.array(density).reshape(-1, 1)
        altitude = np.array(altitude).reshape(-1, 1)
        
        self.__a = LinearRegression()
        self.__a.fit(altitude, a) 
        
        self.__b = LinearRegression()
        self.__b.fit(altitude, b) 
        
        self.__c = LinearRegression()
        self.__c.fit(altitude, c) 
        
        self.__d = LinearRegression()
        self.__d.fit(altitude, d) 
        
        self.__density = LinearRegression()
        self.__density.fit(altitude, density) 
        
    def get_a(self, altitude):
        return self.__a.predict(np.array(altitude).reshape(-1, 1))[0][0]
    
    def get_b(self, altitude):
        return self.__b.predict(np.array(altitude).reshape(-1, 1))[0][0]
    
    def get_c(self, altitude):
        return self.__c.predict(np.array(altitude).reshape(-1, 1))[0][0]
    
    def get_d(self, altitude):
        return self.__d.predict(np.array(altitude).reshape(-1, 1))[0][0]
    
    def get_density(self, altitude):
        return self.__density.predict(np.array(altitude).reshape(-1, 1))[0][0]
        
OS061_12x6 = Engine('OS061_12x6',
                    a = [0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02],
                    b = [-1.42, -1.40, -1.37, -1.35, -1.32, -1.30, -1.27, -1.25, -1.22, -1.20, -1.17, -1.15, -1.12, -1.0734],
                    c = [54.13, 53.33, 52.54, 51.75, 50.97, 50.19, 49.41, 48.63, 47.86, 47.09, 46.33, 45.56, 44.81, 43.2335],
                    d = [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],
                    density = [1.225225683, 1.213506234, 1.201872738, 1.190324756, 1.178861852, 1.16748359, 1.156189537, \
                               1.144979259, 1.133852325, 1.122808306, 1.111846771, 1.100967292, 1.09016944, 1.06881694],
                    altitude = [0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1400] )

class Performance:
    
    def __init__(self):
        pass

    def ground_effect(self, H, CL):
        
        """
        Calculate the parameters for the ground effect  correction.
    
        Parameters:
            H (float): Height of the wing above the ground [m].
            CL (float): Lift coefficient [-].
            
        Returns:
            G (float): Ground effect parameters.
            G2 (float): Ground effect parameters.
        """
        
        wingspan = self.get_Bref()
        S = self.get_Sref()
        N = self.get_Nwings()

        Ra = wingspan**2/(S/N)
        dl = 1-2.25*(0.3**0.00273-0.997)*(Ra**0.717+13.6)

        G = (1-2/np.pi+(16/np.pi*H/wingspan)**2)/(1+(16/np.pi*H/wingspan)**2)
        G2 = (1+dl*288*(H/wingspan)**0.787*math.exp(-9.14*(H/wingspan)**0.327)/Ra**0.882) \
            /(1+0.269*CL**1.45/Ra**3.18*(H/wingspan)*1.12)
            
        return G, G2
    
    def set_engine(self, engine = OS061_12x6):
        """
        Set the engine used by the team.
    
        Parameters:
            engine (Engine): The engine used by the team.
 
        """
        
        self.engine = engine
        
    def takeoff(self, altitude, Sg_max, Sa_max = 0, Ur = 0.13, H = 0.2, mass_inf =  1,
                mass_sup = 30, dt = 0.0001, Error = 0.001):
        
        """
        Calculate the maximum take-off mass for a given altitude.
    
        Parameters:
            altitude (float): Altitude [m], International Standard Atmosphere (ISA).
            Sg_max (float): Maximum ground roll [m].
            Sa_max (float): Maximum horizontal distance While Airborne to Clear an Obstacle [m].
            Ur (float): Friction coefficient.
            H (float): Height of the wing above the ground [m].
            mass_inf (float): Minimum mass for bisection method [m].
            mass_sup (float): Maximum mass for bisection method [m].
            dt (float): Temporal step [s].
            Error (float): Absolute admisible error for mass [kg].
            
        Returns:
            mass (float): Maximum take-off mass in kg.   
        """
        
        density = self.engine.get_density(altitude)
        a = self.engine.get_a(altitude)
        b = self.engine.get_b(altitude)
        c = self.engine.get_c(altitude)
        d = self.engine.get_d(altitude) 
        
        wingspan = self.get_Bref()
        S = self.get_Sref()
        N = self.get_Nwings()
        
        
        CLmin = self.get_CLmin()
        CLmax = self.get_CLmax()
        CL0 = self.get_CL0()
        CD0 = self.get_CD0()
        K = self.get_K()
        
        g     = 9.8066          
        mass  = (mass_inf + mass_sup)/2
            
        CL = CL0
        G, G2 = self.ground_effect(H, CL)
        
        if all([CLmax<0,
                CL0<0,
                CD0<0,
                K<0]):
            return 0
        
        # Algorithm for the calculation of the maximum take-off mass
        while (mass_sup-mass_inf)/2 > Error:
            t  = 0               # Inicializacion tiempo
            V0 = 0.02            # Inicializacion velocidad
            Sg = 0               # Inicializacion distancia recorrida
            
            W = mass*g
            Ka = -density/( 2*(W/S) )*(CD0 + G*K*((CL-CLmin)*G2)**2 - Ur*CL*G2)
            V_stall = ( W / (0.5*density*CLmax*S) )**0.5
            V_LOF    = 1.2*V_stall
            
            
            
            # Ground Roll
            while V0 < V_LOF:
                # Aceleración
                dv_dt = g*( (a/W+Ka)*V0**2 + b/W*V0 + c/W-Ur + (d/W)/V0 )
                V1 = V0 + dv_dt*dt
                V = 0.5*(V0+V1)
                Sg = Sg + V*dt
                V0 = V1
                t += dt;
                
                if t > 300:
                    mass = 0.001
                    return mass
                
            # We reset the limits
            if Sg > Sg_max:
                mass_sup = mass
            else:
                mass_inf = mass
        

            mass = (mass_sup + mass_inf) / 2    

        return mass           

    
    def landing(self, altitude, mass, Ur = 0.13, H = 0.2, dt = 0.01, Error = 0.001):
        
        """
        Calculate the ground roll in landing for a given altitude.
    
        Parameters:
            altitude (float): Altitude [m], International Standard Atmosphere (ISA).
            Ur (float): Friction coefficient.
            H (float): Height of the wing above the ground [m].
            dt (float): Temporal step [s].
            Error (float): Absolute admisible error for mass [kg].
            
        Returns:
            distance (float): Ground roll for landing.   
        """
                                  
        g         = 9.8066                                                
        distance  = [float("+inf"), 0]     
        
        density = self.engine.get_density(altitude)
        CLmax = self.get_CLmax()
        CL0 = self.get_CL0()
        CD0 = self.get_CD0()
        K = self.get_K()
        S = self.get_Sref()
        
        CL    = CL0
        G, G2 = self.ground_effect(H, CL)
        
        if all([CLmax<0,
                CL0<0,
                CD0<0,
                K<0]):
            return 0
        
        # Algorithm for calculating landing distance
        while abs( distance[0] - distance[1] ) > Error:
        
            distance[1] = distance[0]               
            distance[0] = 0                          
            vstall = (mass*g/(density*S*CLmax/2))**(1/2) 
            v=1.1*vstall                             
            t = 0                                    
            
            while v > 0:                          
                CD = CD0 + G*K*(CL*G2)**2                  
                D  = 0.5 * density * v**2 * S * CD        
                L  = 0.5 * density * v**2 * S * CL * G2   
                dv_dt = -D/mass - Ur*(g-L/mass)                
                dx = v*dt + dv_dt*dt**2/2                  
                distance[0] = distance[0]+dx           
                t = t+dt                                  
                v = v+dt*dv_dt                           

            dt = dt/2                                

        return distance[0]          
    
    def speeds(self, altitude, W):
        
        """
        Calculate performance speeds for a given altitude.
    
        Parameters:
            altitude (float): Altitude [m], International Standard Atmosphere (ISA).
            W (float): Weight [kg].
            
        Returns:
            V_min (float): Minimum speed [m/s].   
            V_max (float): Maximum speed [m/s].  
            V_carson (float): Carson speed [m/s].  
            V_cruise (float): Cruise speed [m/s].  
            V_loiter (float): Loiter speed [m/s].  
        """
        
        density = self.engine.get_density(altitude)
        a = self.engine.get_a(altitude)
        b = self.engine.get_b(altitude)
        c = self.engine.get_c(altitude)
        d = self.engine.get_d(altitude) 
        CD0 = self.get_CD0()
        S = self.get_Sref()
        K = self.get_K()
        
        
        # Define symbols used in the equation
        V = symbols('V')
        
        # # Define the Rate of Climb as Equal to 0
        # dH_dt = (a*V**3 + b*V**2 + c*V + d)/W - \
        #     density*V**3*S/(2*W) * (CD0+K*(2*W/(density*V**2*S))**2)
        
        # # Find the zeros of the equation
        # min_max = solve(dH_dt, V)
        # min_max = [v for v in min_max if v >= 0]
        
        # # Obtain the minimum and maximum values of velocity
        # V_min = min(min_max)
        # V_max = max(min_max)
        
        # Obtain the Carson, Cruise and Loiter speeds
        V_carson = (( 2/density) * (3*K/CD0)**0.5 * W/S )**0.5
        V_cruise = (( 2/density) * (  K/CD0)**0.5 * W/S )**0.5 
        V_loiter = (( 2/density) * (K/(3*CD0))**0.5 * W/S )**0.5
        
        Pav = []
        Preq = []
        for vel in np.linspace(start = 1, stop = 30, num = 16):
            Pav.append( [vel, a*vel**3 + b*vel**2 + c*vel + d] )
            Preq.append( [vel, density*vel**3*S/2 * (CD0+K*(2*W/(density*vel**2*S))**2)] )
            
        # V_min, V_max ,
        
        return V_carson ,V_cruise, V_loiter, Pav, Preq  
