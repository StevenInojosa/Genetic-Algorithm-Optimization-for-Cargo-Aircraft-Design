#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: stevenstefanoinojosasiso
"""

from aerodynamics import *
from performance import *
from stability import *
from structures import *

class Aircraft(Aerodynamics, Performance, Stability, Structures):
    """
    A class that represents an aircraft.

    Parameters
    ----------
    name : str
        The name of the aircraft.
    Sref : float, optional
        The reference area of the aircraft in square meters.
    Cref : float, optional
        The reference chord of the aircraft in meters.
    Bref : float, optional
        The reference wingspan of the aircraft in meters.
    Xref : float, optional
        The longitudinal reference position of the aircraft in meters.
    Yref : float, optional
        The lateral reference position of the aircraft in meters.
    Zref : float, optional
        The vertical reference position of the aircraft in meters.
    Nwings : int, optional
        The number of wings of the aircraft.

    Attributes
    ----------
    __name : str
        The name of the aircraft.
    __Sref : float
        The reference area of the aircraft in square meters.
    __Cref : float
        The reference chord of the aircraft in meters.
    __Bref : float
        The reference wingspan of the aircraft in meters.
    __Xref : float
        The longitudinal reference position of the aircraft in meters.
    __Yref : float
        The lateral reference position of the aircraft in meters.
    __Zref : float
        The vertical reference position of the aircraft in meters.
    __Nwings : int
        The number of wings of the aircraft.
    __surfaces : dict
        A dictionary containing the surfaces of the aircraft.
    """
    def __init__(self, name, Sref = 0, Cref = 0, Bref = 0, Xref = 0, Yref = 0,
                 Zref = 0, Nwings = 1): 
        """
        Initializes an instance of the Aircraft class.

        Parameters
        ----------
        name : str
            The name of the aircraft.
        Sref : float, optional
            The reference area of the aircraft in square meters.
        Cref : float, optional
            The reference chord of the aircraft in meters.
        Bref : float, optional
            The reference wingspan of the aircraft in meters.
        Xref : float, optional
            The longitudinal reference position of the aircraft in meters.
        Yref : float, optional
            The lateral reference position of the aircraft in meters.
        Zref : float, optional
            The vertical reference position of the aircraft in meters.
        Nwings : int, optional
            The number of wings of the aircraft.
        """
        self.__name = name
        self.__Sref = Sref
        self.__Cref = Cref
        self.__Bref = Bref
        self.__Xref = Xref
        self.__Yref = Yref
        self.__Zref = Zref
        
        self.__Nwings = Nwings
        
        self.__surfaces = {}
        
    def update_dim_ref():
        pass
        
    def add_surface(self, name, surface):
        self.__surfaces[name] = surface
        
    def get_surface(self, name = ""): 
        """
        Returns a dictionary containing the surfaces of the aircraft, or
        the surface with the given name.

        :param name: A string indicating the name of the surface to return.
        :return: A dictionary containing the surfaces of the aircraft, or
            the surface with the given name.
        """
        if name == "":
            return self.__surfaces
        else:
            return self.__surfaces[name]
        
    def get_name(self):
        """
        Returns the name of the aircraft.

        :return: A string indicating the name of the aircraft.
        """
        return self.__name
    
    def get_Sref(self):
        """
        Returns the reference area of the aircraft.

        :return: A float indicating the reference area of the aircraft.
        """
        return self.__Sref
    
    def get_Cref(self):
        """
        Returns the reference chord of the aircraft.

        :return: A float indicating the reference chord of the aircraft.
        """
        return self.__Cref
    
    def get_Bref(self):
        """
        Returns the reference span of the aircraft.

        :return: A float indicating the reference span of the aircraft.
        """
        return self.__Bref
    
    def get_Xref(self):
        """
        Returns the reference X coordinate of the aircraft.

        :return: A float indicating the reference X coordinate of the aircraft.
        """
        return self.__Xref
    
    def get_Yref(self):
        """
        Returns the reference Y coordinate of the aircraft.

        :return: A float indicating the reference Y coordinate of the aircraft.
        """
        return self.__Yref
    
    def get_Zref(self):
        """
        Returns the reference Z coordinate of the aircraft.

        :return: A float indicating the reference Z coordinate of the aircraft.
        """
        return self.__Zref
    
    def get_Nwings(self):
        """
        Returns the number of wings of the aircraft.

        :return: An integer indicating the number of wings of the aircraft.
        """
        return self.__Nwings
        
class Surface:
    def __init__(self, translate, dAinc):
        """
        Initializes an instance of the Surface class with the given parameters.

        :param translate: A list of three floats indicating the translation
            vector of the surface.
        :param dAinc: A float indicating the incremental surface area of the
            surface.
        """
        self.__translate = translate
        self.__dAinc = dAinc        

        self.__sections = []
        
    def add_section(self, section):
        self.__sections.append(section)
        
        try:
            self.chords()
        except:
            pass
        
    def get_sections(self):
        """
        Returns the list of sections of the surface.

        :return: A list of Section objects.
        """
        return self.__sections
        
    def get_translate(self, index = -1):
        """
        Returns the translation vector of the surface, or the component
        at the given index.

        :param index: An integer indicating the index of the component to return.
            If not provided, the whole translation vector
        """
        if index == -1:
            return self.__translate
        else:
            return self.__translate[index]
    
    def get_dAinc(self):
        """
        Returns the incremental change in the wing's lift due to a change in angle of incidence.
        """
        return self.__dAinc
            
    def chords(self):
        """
        Computes the wing's mean aerodynamic chord (MAC), its reference area (S), and its reference area divided by its 
        mean chord (c).
        """
        try:
            y = sp.symbols(('y'))
            I_mac = 0
            S = 0
            for ii in range(len(self.get_sections())-1):
                c0 = self.get_sections()[ii].get_Chord()
                c1 = self.get_sections()[ii+1].get_Chord()
                y0 = self.get_sections()[ii].get_Yle()
                y1 = self.get_sections()[ii+1].get_Yle()
                
                if y0 == y1:
                    break
                
                c_y = (c1-c0)/(y1-y0)*y + c0
                
                I_mac += sp.integrate(c_y**2, (y,0,y1-y0))
                S += 2*(y1-y0)*(c1+c0)/2
        
            self.__mac = (2/S)*I_mac
            self.__smc = S/(2*self.get_sections()[-1].get_Yle())
            self.__S = S
        except: 
            self.__mac = self.__smc = ( self.get_sections()[0].get_Chord() + 
                                        self.get_sections()[1].get_Chord() ) / 2
            bv = self.get_Yle()
            self.__S = SMCv * self.__smc
        
    def get_mac(self):
        """
        Returns the wing's mean aerodynamic chord (MAC).
        """
        return self.__mac
    
    def get_smc(self):
        """
        Returns the wing's reference area (S) divided by its mean chord (c).
        """
        return self.__smc
       
    def get_S(self):
        """
        Returns the wing's reference area.
        """
        return self.__S

class Section:
    """A class to represent a section of an aircraft wing.

    Attributes:
        Xle (float): X-coordinate of the leading edge of the section
        Yle (float): Y-coordinate of the leading edge of the section
        Zle (float): Z-coordinate of the leading edge of the section
        Chord (float): Chord length of the section
        Ainc (float): Angle of incidence of the section
        airfoil (tuple): Tuple containing the name of the airfoil and its thickness
        Nspan (int): Number of spanwise vortices
        Sspace (float): Spanwise distance between vortices
    
    Methods:
        get_Xle(): Returns the X-coordinate of the leading edge
        get_Yle(): Returns the Y-coordinate of the leading edge
        get_Zle(): Returns the Z-coordinate of the leading edge
        get_Chord(): Returns the chord length
        get_Ainc(): Returns the angle of incidence
        get_Nspan(): Returns the number of spanwise vortices
        get_Sspace(): Returns the spanwise distance between vortices
        get_airfoil(get='name'): Returns the name or thickness of the airfoil
        
    """
    def __init__(self, Xle, Yle, Zle, Chord, Ainc, airfoil,
                 Nspan=12, Sspace=-1):
        """
        Constructor for the Section class.
        
        Parameters:
        -----------
        Xle : float
            X-coordinate of the leading edge.
        Yle : float
            Y-coordinate of the leading edge.
        Zle : float
            Z-coordinate of the leading edge.
        Chord : float
            Chord length.
        Ainc : float
            Angle of incidence (in degrees).
        airfoil : tuple
            Tuple containing the name of the airfoil file and thickness value.
        Nspan : int, optional
            Number of spanwise vortices. Default is 12.
        Sspace : float, optional
            Spanwise distance between vortices. Default is -1, which means automatic spacing based on Nspan.
        """
        self.__Xle = Xle
        self.__Yle = Yle
        self.__Zle = Zle
        self.__Chord = Chord
        self.__Ainc = Ainc 
        self.__Nspan = Nspan
        self.__Sspace = Sspace
        self.__airfoil = airfoil
        
    def get_Xle(self):
        """Returns the X-coordinate of the leading edge"""
        return self.__Xle 
   
    def get_Yle(self):
        """Returns the Y-coordinate of the leading edge"""
        return self.__Yle 
   
    def get_Zle(self):
        """Returns the Z-coordinate of the leading edge"""
        return self.__Zle 
   
    def get_Chord(self):
        """Returns the chord length"""
        return self.__Chord 
   
    def get_Ainc(self):
        """Returns the angle of incidence"""
        return self.__Ainc
   
    def get_Nspan(self):
        """Returns the number of spanwise vortices"""
        return self.__Nspan 
   
    def get_Sspace(self):
        """Returns the spanwise distance between vortices"""
        return self.__Sspace
   
    def get_airfoil(self, get='name'):
        """Returns the name or thickness of the airfoil

        Args:
            get (str): Specifies whether to return the name or thickness. Defaults to 'name'.

        Returns:
            str: The name or thickness of the airfoil
        """
        if get == 'name':
            return self.__airfoil[0]
        elif get == 't':
            return self.__airfoil[1]
