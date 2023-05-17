#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

@author: stevenstefanoinojosasiso
"""

# It is necessary to delete all files at the start, say if is MacOS or Windows

import os
import subprocess as sp
import pandas as pd
import numpy as np
import time

class AVL_python:
    """
    A class for generating AVL geometries from an Aircraft object.
    """
    
    def __init__(self, Alpha, Roll, Pitch, Yaw, Beta, Aircraft):
        """
        Initializes a new AVL_python object.

        Args:
            Alpha (float): Angle of attack in degrees.
            Roll (float): Roll angle in radians.
            Pitch (float): Pitch angle in radians.
            Yaw (float): Yaw angle in radians.
            Beta (float): Sideslip angle in degrees.
            Aircraft (Aircraft): An Aircraft object containing the aircraft configuration.

        Returns:
            None.
        """
        self.__Aircraft = Aircraft
        
        self.__Alpha = Alpha
        self.__Roll  = Roll
        self.__Pitch = Pitch
        self.__Yaw   = Yaw
        self.__Beta  = Beta
        
        self.__cwd = os.getcwd() 
        
        if '\\' in self.__cwd:       # Is Windows
            self.__cwd = ''
            self.__folder1 = '\\Geometries\\' 
            self.__folder2 = '\\Airfoils\\' 
            self.__airfoil
        else:                        # Is MacOS or Linux
            self.__folder1 = '/Geometries/' 
            self.__folder2 = '/Airfoils/' 
            
    def get_Aircraft(self):
        """
        Returns the current aircraft object.
        
        Returns:
        Aircraft: An Aircraft object representing the current aircraft.
        """
        return self.__Aircraft
    
    def get_Alpha(self):
        """
        Returns a list of angles of attack for simulation.
        
        Returns:
        list: A list of angles of attack for simulation.
        """
        return self.__Alpha
    
    def get_Roll(self):
        """
        Returns a list of roll angles for simulation.
        
        Returns:
        list: A list of roll angles for simulation.
        """
        return self.__Roll
    
    def get_Pitch(self):
        """
        Returns a list of pitch angles for simulation.
        
        Returns:
        list: A list of pitch angles for simulation.
        """
        return self.__Pitch
    
    def get_Yaw(self):
        """
        Returns a list of yaw angles for simulation.
        
        Returns:
        list: A list of yaw angles for simulation.
        """
        return self.__Yaw
    
    def get_Beta(self):
        """
        Returns a list of sideslip angles for simulation.
        
        Returns:
        list: A list of sideslip angles for simulation.
        """
        return self.__Beta
    
    def get_cwd(self):
        """
        Returns the current working directory.
        
        Returns:
        str: The current working directory.
        """
        return self.__cwd
    
    def get_folder1(self):
        """
        Returns the name of the first folder.
        
        Returns:
        str: The name of the first folder.
        """
        return self.__folder1
    
    def get_folder2(self):
        """
        Returns the name of the second folder.
        
        Returns:
        str: The name of the second folder.
        """
        return self.__folder2


    def geometry_creator(self, text = ''):
        """
        Creates the AVL geometry file content for the aircraft.
    
        Args:
        - self: instance of the AVL_python class
        - text: (optional) string to add to the geometry file (default: '')
    
        Returns:
        - text: string containing the geometry file content
        """
        Aircraft = self.get_Aircraft()
        
        text += '# Este archivo fue creado usando un programa de Python \n'
        text += '# Fecha: %s \n\n' %(time.ctime())
        text += '%s \n\n' % Aircraft.get_name()
        text += '0.0                             | Mach \n'
        text += '0       0        0.0            | iYsym  iZsym  Zsym \n'
        text += '%.3f    %.3f    %.3f         | Sref   Cref   Bref \n' %(
            Aircraft.get_Sref(), Aircraft.get_Cref(), Aircraft.get_Bref() )
        text += '%.3f    %.3f    %.3f         | Xref   Yref   Zref \n' %(
            Aircraft.get_Xref(), Aircraft.get_Yref(), Aircraft.get_Zref() )
        
        try:
            text += '%.3f                           | CDp    \n\n\n' % Aircraft.get_CDp('Aircraft')
        except:
            text += '%.3f                           | CDp    \n\n\n' %(0.0)
        

        for name, surface in Aircraft.get_surface().items():
            text += '#================================================= \n'
            text += 'SURFACE                         | (keyword) \n'
            text += '%s \n' % name
            text += '#Nchord   Cspace     [ Nspan Sspace ] \n'
            text += '10        1.0 \n\n'

            text += 'INDEX                           | (keyword) \n'
            text += '31621                           | Lsurf \n\n'

            text += 'YDUPLICATE \n'
            text += '0.0 \n\n'

            text += 'SCALE \n'
            text += '1.0  1.0  1.0 \n\n'

            text += 'TRANSLATE \n'
            text += '%.3f    %.3f    %.3f \n\n' %(
                surface.get_translate(0), surface.get_translate(1), surface.get_translate(2))

            text += 'ANGLE \n'
            text += '  %.3f                         | dAinc \n\n' % surface.get_dAinc()
            
            for section in surface.get_sections(): #surface.get_sections():
                text += '#______________ \n'
                text += 'SECTION                         |  (keyword) \n\n'
                text += '# Xle    Yle      Zle      Chord    Ainc    [Nspan Sspace] \n'
                text += '%.3f    %.3f    %.3f    %.3f    %.3f    %.0f    %.0f \n\n' %(
                    section.get_Xle(), section.get_Yle(), section.get_Zle(),
                    section.get_Chord(), section.get_Ainc(),
                    section.get_Nspan(), section.get_Sspace() )  
                text += 'AFIL 0.0 1.0 \n'
                text += self.get_cwd() + self.get_folder2() + section.get_airfoil() + '.dat' + '\n\n'
                
        f = open(self.get_cwd() + self.get_folder1() + Aircraft.get_name() + '.avl', 'w') 
        f.write(text)
        f.close()            
        return text
        
    def cmd(self, cmd = ""):
        """
        Generates the command string to send to the AVL executable based on the aircraft properties
        and simulation parameters stored in the AVLModel object.
    
        Args:
            cmd (str): A string containing a command to append to the generated command string. Default is an empty string.
    
        Returns:
            The encoded command string to send to the AVL executable.
        """
        Aircraft = self.get_Aircraft()
               
        cmd += 'LOAD ' + self.get_cwd() + self.get_folder1() + Aircraft.get_name() + '\n' #
        cmd += 'OPER' + '\n'
        
        # Alpha simulation
        for ii in self.get_Alpha():
            cmd += 'A A ' + str(ii) + '\n X \n'
        cmd += 'A A 0 \n X \n' # Initialize
        
        # Roll simulation
        for ii in self.get_Roll():
            cmd += 'R R ' + str(ii) + '\n X \n'
        cmd += 'R R 0 \n X \n' # Initialize
        
        # Pitch simulation
        for ii in self.get_Pitch():
            cmd += 'P P ' + str(ii) + '\n X \n'
        cmd += 'P P 0 \n X \n' # Initialize

        # Yaw simulation
        for ii in self.get_Yaw():
            cmd += 'Y Y ' + str(ii) + '\n X \n'
        cmd += 'Y Y 0 \n X \n' # Initialize
        
        # Beta simulation
        for ii in self.get_Beta():
            cmd += 'B B ' + str(ii) + '\n X \n'
        cmd += 'B B 0 \n X \n' # Initialize
            
        return cmd.encode()
    
    def run(self):
        """
        Runs the AVL executable and parses its output to extract aerodynamic coefficients.
    
        Returns:
            A tuple containing the following elements:
            - out: A string with the raw output generated by the AVL executable.
            - cmd: The command that was sent to the executable as input.
            - alpha_table: A dictionary containing the following keys and values:
              - 'Deg': A list of float values with the angles of attack (in degrees) for which data was generated.
              - 'Cm': A list of float values with the corresponding pitching moments.
              - 'CL': A list of float values with the corresponding lift coefficients.
              - 'CLff': A list of float values with the corresponding lift coefficients for the flat plate airfoil.
              - 'CDtot': A list of float values with the corresponding drag coefficients (total).
              - 'CDind': A list of float values with the corresponding drag coefficients (induced).
              - 'CDff': A list of float values with the corresponding drag coefficients for the flat plate airfoil.
              - 'CL^2': A list of float values with the square of the lift coefficients.
              - 'CLff^2': A list of float values with the square of the lift coefficients for the flat plate airfoil.
            - roll_table: A dictionary containing the following keys and values:
              - 'Deg': A list of float values with the roll angles (in degrees) for which data was generated.
              - 'rad': A list of float values with the corresponding roll angles in radians.
              - 'Cy': A list of float values with the corresponding side force coefficients.
              - 'Cn': A list of float values with the corresponding yawing moment coefficients.
              - 'Cl': A list of float values with the corresponding rolling moment coefficients.
            - pitch_table: A dictionary containing the following keys and values:
              - 'Deg': A list of float values with the pitch angles (in degrees) for which data was generated.
              - 'rad': A list of float values with the corresponding pitch angles in radians.
              - 'Cz': A list of float values with the corresponding axial force coefficients.
              - 'Cm': A list of float values with the corresponding pitching moment coefficients.
            - yaw_table: A dictionary containing the following keys and values:
              - 'Deg': A list of float values with the yaw angles (in degrees) for which data was generated.
              - 'rad': A list of float values with the corresponding yaw angles in radians.
              - 'Cy': A list of float values with the corresponding side force coefficients.
              - 'Cn': A list of float values with the corresponding yawing moment coefficients.
              - 'Cl': A list of float values with the corresponding rolling moment coefficients.
            - beta_table: A dictionary containing the following keys and values:
              - 'Deg': A list of float values with the sideslip angles (in degrees) for which data was generated.
              - 'rad': A list of float values with the corresponding sideslip angles in radians.
              - 'Cy': A list of float values with the corresponding side force coefficients.
              - 'Cn': A list of float values with the corresponding yawing moment coefficients.
              - 'Cl': A list of float values with the corresponding rolling moment coefficients.
        """
        
        p = sp.Popen([self.get_cwd() + '/AVL'], stdin=sp.PIPE, stdout=sp.PIPE, shell= True) # stdout=sp.PIPE,
        
        cmd = self.cmd()
        out, err = p.communicate(input =  cmd)
        out = out.decode('utf-8').split()
        
        # We obtain the index for every Alpha in AVL
        cases = [index for index, txt in enumerate(out) if txt == 'Alpha']
        cases_validator = len(self.get_Alpha())*['alpha'] + [''] + \
                          len(self.get_Roll())*['roll'] + [''] + \
                          len(self.get_Pitch())*['pitch'] + [''] + \
                          len(self.get_Yaw())*['yaw'] + [''] + \
                          len(self.get_Beta())*['beta'] + [''] 
        

        counter = 0    
        aa = 0
        rr = 0
        pp = 0
        yy = 0
        bb = 0
        
        alpha_table = {'Deg':[], 'Cm':[], 'CL':[], 'CLff':[], 'CDtot':[],
                       'CDind':[], 'CDff':[], 'CL^2':[], 'CLff^2':[]}
        roll_table = {'Deg':[], 'rad':[], 'Cy':[], 'Cn':[], 'Cl':[]}
        pitch_table = {'Deg':[], 'rad':[], 'Cz':[], 'Cm':[]}
        yaw_table = {'Deg':[], 'rad':[], 'Cy':[], 'Cn':[], 'Cl':[]}
        beta_table = {'Deg':[], 'rad':[], 'Cy':[], 'Cn':[], 'Cl':[]}
        
        for index in cases:
            alpha = float(out[ index + 2  ])
            roll  = float(out[ index + 5  ])
            pitch = float(out[ index + 14 ])
            yaw   = float(out[ index + 20 ])
            beta  = float(out[ index + 11 ])
            
            Cl    = float(out[ index + 29 ])
            Cm    = float(out[ index + 38 ])
            Cn    = float(out[ index + 44 ])
            
            Cy    = float(out[ index + 35 ])
            Cz    = float(out[ index + 41 ])
            
            CL    = float(out[ index + 50 ])
            CLff  = float(out[ index + 62 ])
            CDtot = float(out[ index + 53 ])
            CDind = float(out[ index + 59 ])    
            CDff  = float(out[ index + 65 ])
            
            if cases_validator[counter] == 'alpha':
                # Deg Cm	CL	CLff	Cdind	cdff	CL^2	CLff^2
                alpha_table['Deg'].append(alpha)
                alpha_table['Cm'].append(Cm)
                alpha_table['CL'].append(CL)
                alpha_table['CLff'].append(CLff)
                alpha_table['CDtot'].append(CDtot)
                alpha_table['CDind'].append(CDind)
                alpha_table['CDff'].append(CDff)
                alpha_table['CL^2'].append(CL**2)
                alpha_table['CLff^2'].append(CLff**2)
                
            elif cases_validator[counter] == 'roll':
                # Deg	rad	Cy	Cn	Cl
                roll_table['Deg'].append(roll*180/np.pi)
                roll_table['rad'].append(roll)
                roll_table['Cy'].append(Cy)
                roll_table['Cn'].append(Cn)
                roll_table['Cl'].append(Cl)
            
            elif cases_validator[counter] == 'pitch':
                # Deg	rad	Cz	Cm
                pitch_table['Deg'].append(pitch*180/np.pi)
                pitch_table['rad'].append(pitch)
                pitch_table['Cz'].append(Cz)
                pitch_table['Cm'].append(Cm)

            elif cases_validator[counter] == 'yaw':
                # Deg	rad	Cy	Cn	Cl
                yaw_table['Deg'].append(yaw*180/np.pi)
                yaw_table['rad'].append(yaw)
                yaw_table['Cy'].append(Cy)
                yaw_table['Cn'].append(Cn)
                yaw_table['Cl'].append(Cl)
            
            elif cases_validator[counter] == 'beta':
                # B	Cy	Cl	Cn
                beta_table['Deg'].append(beta)
                beta_table['rad'].append(beta*np.pi/180)
                beta_table['Cy'].append(Cy)
                beta_table['Cn'].append(Cn)
                beta_table['Cl'].append(Cl)

            counter += 1
        
            
        return out, cmd, alpha_table, roll_table, pitch_table, yaw_table, beta_table