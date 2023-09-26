#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: stevenstefanoinojosasiso
"""

import numpy as np
import math

class Stability:
    
    def static_stability(self):
        """
        Determine the static stability of the aircraft in pitch, roll, and yaw directions.
        
        Sets instance variables for longitudinal, lateral, and directional static stability, where
        True indicates stable and False indicates unstable.
        
        Returns:
            None
        """
        self.__lon_static_stability = self.get_Cma() < 0 
        self.__lat_static_stability = self.get_CnB() > 0
        self.__dir_static_stability = self.get_ClB() < 0
        
    def is_point_below_line(self, p1, p2, p):
        """
        Checks if the y-coordinate of a given point is less than the y-coordinate on a line
        passing through two given points, corresponding to the point's x-coordinate.
    
        Args:
            p1 (tuple): The coordinates (x1, y1) of the first point on the line.
            p2 (tuple): The coordinates (x2, y2) of the second point on the line.
            p (tuple): The coordinates (x, y) of the point to compare the y-coordinate with.
    
        Returns:
            bool: True if the y-coordinate of the point is less than the corresponding y-coordinate on the line,
                  False otherwise.
    
        Formula:
            The function uses the slope-intercept form of a linear equation (y = mx + c),
            where m is the slope of the line and c is the y-intercept.
    
            The slope (m) is calculated as (y2 - y1) / (x2 - x1).
    
            The y-coordinate on the line is determined using the formula:
            y = m(x - x1) + y1
    
            The function compares the y-coordinate of the point with the calculated y-coordinate on the line
            and returns True if the point's y-coordinate is less than the line's y-coordinate, and False otherwise.
    
        """
        x1, y1 = p1
        x2, y2 = p2
        x, y = p
    
        return y < (y2 - y1) / (x2 - x1) * (x - x1) + y1


    def phase_B1(self, short_period):
        """
        Determines the phase B classification based on the short period characteristics.
    
        Args:
            short_period (dict): A dictionary containing the short period characteristics.
                                 It should have keys 'n/alfa' and 'w' representing the values of
                                 the parameter 'n/alfa' and the undamped natural frequency 'w' respectively.
    
        Returns:
            int: The phase B classification value:
                 - 4 if the point is below the line defined by (1, 0.2) and (100, 2).
                 - 2 if the point is below the line defined by (1, 0.3) and (100, 3).
                 - 1 if the point is below the line defined by (1, 2) and (100, 11).
                 - 2 if the point is below the line defined by (1, 3) and (100, 12).
                 - 3 if none of the above conditions are met.
    
        """
        x = short_period['n/alfa']
        y = short_period['w']
        
        if self.is_point_below_line((1, 0.2), (100, 2), (x, y)):
            return 4
        elif self.is_point_below_line((1, 0.3), (100, 3), (x, y)):
            return 2
        elif self.is_point_below_line((1, 2), (100, 11), (x, y)):
            return 1
        elif self.is_point_below_line((1, 3), (100, 12), (x, y)):
            return 2
        else:
            return 3

    def phase_C1(self, short_period):
        """
        Determines the phase C classification for short period characteristics.
    
        Args:
            short_period (dict): A dictionary containing the short period characteristics.
                                 It should have keys 'n/alfa' and 'w' representing the values of
                                 the parameter 'n/alfa' and the undamped natural frequency 'w' respectively.
    
        Returns:
            int: The phase C classification value   
        """
        
        x = short_period['n/alfa']
        y = short_period['w']
        
        if x < 1.81:
            if self.is_point_below_line( (1, 0.6), (1.81, 0.6), (x, y) ):
                return 4
            else:
                return 3
        
        elif 1.81 <= x < 2.82:
            if self.is_point_below_line( (1.81, 0.60), (2.82, 0.60), (x, y) ):
                return 4
            elif self.is_point_below_line( (1.81, 4.07), (100, 30), (x, y) ):
                return 2
            else:
                return 3
                
        elif 2.82 <= x < 3.87:
            if self.is_point_below_line( (2.82, 0.60), (3.87, 0.60), (x, y) ):
                return 4
            elif self.is_point_below_line( (2.82, 0.84), (3.87, 0.84), (x, y) ):
                return 2
            elif self.is_point_below_line( (2.82, 3.27), (100, 20), (x, y) ):
                return 1
            elif self.is_point_below_line( (1.81, 4.07), (100, 30), (x, y) ):
                return 2
            else:
                return 3
            
        
        elif 3.87 <= x < 4.93:
            if self.is_point_below_line( (3.87, 0.60), (100, 3), (x, y) ):
                return 4
            elif self.is_point_below_line( (3.87, 0.84), (4.93, 0.84), (x, y) ):
                return 2
            elif self.is_point_below_line( (2.82, 3.27), (100, 20), (x, y) ):
                return 1
            elif self.is_point_below_line( (1.81, 4.07), (100, 30), (x, y) ):
                return 2
            else:
                return 3
    
        else:
            if self.is_point_below_line( (3.87, 0.60), (100, 3), (x, y) ):
                return 4
            elif self.is_point_below_line( (4.93, 0.84), (100, 4), (x, y) ):
                return 2
            elif self.is_point_below_line( (2.82, 3.27), (100, 20), (x, y) ):
                return 1
            elif self.is_point_below_line( (1.81, 4.07), (100, 30), (x, y) ):
                return 2
            else:
                return 3
            
    def phase_B2(self, short_period):
        """
        Determines the phase B classification for short period characteristics.
    
        Args:
            short_period (dict): A dictionary containing the short period characteristics.
                                 It should have keys 'zeta' and 'CAP' representing the values of
                                 the damping ratio 'zeta' and the control anticipation parameter 'CAP' respectively.
    
        Returns:
            int: The phase B classification value   
        """
        
        x = short_period['zeta']
        y = short_period['CAP']
        
        if x < 0.16:
            return 4
        
        elif x < 0.20:
            if y < 0.27:
                return 4
            else:
                return 3
            
        elif x < 0.30:
            if y < 0.27:
                return 4
            elif y < 10:
                return 2
            else:
                return 3
            
        elif x < 2:
            if y < 0.27:
                return 4
            elif y < 0.72:
                return 2
            elif y < 3.61:
                return 1
            elif y < 10:
                return 2
            else:
                return 3
            
        else:
            if y < 0.27:
                return 4
            else:
                return 3
            
    def phase_C2(self, short_period):
        """
        Determines the phase C classification for short period characteristics.
    
        Args:
            short_period (dict): A dictionary containing the short period characteristics.
                                 It should have keys 'zeta' and 'CAP' representing the values of
                                 the damping ratio 'zeta' and the control anticipation parameter 'CAP' respectively.
    
        Returns:
            int: The phase C classification value   
        """
        
        x = short_period['zeta']
        y = short_period['CAP']
        
        if x < 0.15:
            return 4
        
        elif x < 0.25:
            if y < 0.1:
                return 4
            else:
                return 3
        
        elif x < 0.35:
            if y < 0.1:
                return 4
            elif y < 10:
                return 2
            else:
                return 3
            
        
        elif x < 1.15:
            if y < 0.1:
                return 4
            elif y < 0.15:
                return 2
            elif y < 1.36:
                return 1
            elif y < 10:
                return 2
            else:
                return 3
        
        elif x < 2:
            if y < 0.1:
                return 4
            elif y < 10:
                return 2
            else:
                return 3
        
        else:
            if y < 0.1:
                return 4
            else:
                return 3
            
    def pughoid_level(self, pughoid):
        """
        Determines the level of the pughoid based on the given pughoid dictionary.
    
        Parameters:
            pughoid (dict): Dictionary containing the pughoid parameters.
                - 'zeta' (float): The value of the pughoid damping ratio.
    
        Returns:
            level (int): The level of the pughoid.
    
        Description:
            This function determines the level of the pughoid based on the damping ratio value provided in the pughoid dictionary.
            
            The levels are determined as follows:
            - If the damping ratio (zeta) is greater than 0.04, the pughoid level is set to 1.
            - If the damping ratio (zeta) is greater than 0, but less than or equal to 0.04, the pughoid level is set to 2.
            - If the damping ratio (zeta) is less than or equal to 0, the pughoid level is set to -1.
    
            The calculated pughoid level is returned as the output.
        """
        if pughoid['zeta'] > 0.04:
            return 1
        elif pughoid['zeta'] > 0:
            return 2
        else:
            return -1
        
    def dutch_roll_level_A(self, dutch_roll):
        """
        Determines the level of Dutch Roll phase A based on the given Dutch Roll parameters.
    
        Parameters:
            dutch_roll (dict): Dictionary containing the Dutch Roll parameters.
                - 'zeta' (float): The damping ratio of Dutch Roll.
                - 'w' (float): The undamped natural frequency.
    
        Returns:
            level (int): The level of Dutch Roll phase A.
    
        Description:
            This function determines the level of Dutch Roll phase A based on the provided damping ratio (zeta) and natural frequency (w) values.
            
            The levels are determined as follows:
            - If the damping ratio (zeta) is greater than 0.19, zeta*w is greater than 0.35, and w is greater than 0.4, the Dutch Roll A level is set to 1.
            - If the damping ratio (zeta) is greater than 0.02, zeta*w is greater than 0.05, and w is greater than 0.4, the Dutch Roll A level is set to 2.
            - If the damping ratio (zeta) is greater than 0 and w is greater than 0.4, the Dutch Roll phase A level is set to 3.
            - Otherwise, if none of the above conditions are met, the Dutch Roll A level is set to -1.
    
            The calculated Dutch Roll phase A level is returned as the output.
        """
        if dutch_roll == False:
            return -1
        
        zeta = dutch_roll['zeta']
        w = dutch_roll['w']
        
        if zeta > 0.19 and zeta*w > 0.35 and w > 0.4:
            return 1
        elif zeta > 0.02 and zeta*w > 0.05 and w > 0.4:
            return 2
        elif zeta > 0.00 and w > 0.4:
            return 3
        else:
            return -1
        
    def dutch_roll_level_B(self, dutch_roll):
        """
        Determines the level of Dutch Roll phase B based on the given Dutch Roll parameters.
    
        Parameters:
            dutch_roll (dict): Dictionary containing the Dutch Roll parameters.
                - 'zeta' (float): The damping ratio of Dutch Roll.
                - 'w' (float): The undamped natural frequency.
    
        Returns:
            level (int): The level of Dutch Roll phase B.
    
        Description:
            This function determines the level of Dutch Roll phase B based on the provided damping ratio (zeta) and natural frequency (w) values.
            
            The levels are determined as follows:
            - If the damping ratio (zeta) is greater than 0.08, zeta*w is greater than 0.15, and w is greater than 0.4, the Dutch Roll B level is set to 1.
            - If the damping ratio (zeta) is greater than 0.02, zeta*w is greater than 0.05, and w is greater than 0.4, the Dutch Roll B level is set to 2.
            - If the damping ratio (zeta) is greater than 0 and w is greater than 0.4, the Dutch Roll phase B level is set to 3.
            - Otherwise, if none of the above conditions are met, the Dutch Roll B level is set to -1.
    
            The calculated Dutch Roll phase B level is returned as the output.
        """
        if dutch_roll == False:
            return -1
        
        zeta = dutch_roll['zeta']
        w = dutch_roll['w']
        
        if zeta > 0.08 and zeta*w > 0.15 and w > 0.4:
            return 1
        elif zeta > 0.02 and zeta*w > 0.05 and w > 0.4:
            return 2
        elif zeta > 0.00 and w > 0.4:
            return 3
        else:
            return -1
        
    def dutch_roll_level_C(self, dutch_roll):
        """
        Determines the level of Dutch Roll phase C based on the given Dutch Roll parameters.
        
        Parameters:
            dutch_roll (dict): Dictionary containing the Dutch Roll parameters.
                - 'zeta' (float): The damping ratio of Dutch Roll.
                - 'w' (float): The undamped natural frequency.
        
        Returns:
            level (int): The level of Dutch Roll phase C.
        
        Description:
            This function determines the level of Dutch Roll phase C based on the provided damping ratio (zeta) and natural frequency (w) values.
            
            The levels are determined as follows:
            - If the damping ratio (zeta) is greater than 0.08, zeta*w is greater than 0.15, and w is greater than 1.0, the Dutch Roll C level is set to 1.
            - If the damping ratio (zeta) is greater than 0.02, zeta*w is greater than 0.05, and w is greater than 0.4, the Dutch Roll C level is set to 2.
            - If the damping ratio (zeta) is greater than 0 and w is greater than 0.4, the Dutch Roll phase A level is set to 3.
            - Otherwise, if none of the above conditions are met, the Dutch Roll C level is set to -1.
        
            The calculated Dutch Roll phase C level is returned as the output.
        """
        
        if dutch_roll == False:
            return -1
        
        zeta = dutch_roll['zeta']
        w = dutch_roll['w']
        
        if zeta > 0.08 and zeta*w > 0.15 and w > 1.0:
            return 1
        elif zeta > 0.02 and zeta*w > 0.05 and w > 0.4:
            return 2
        elif zeta > 0.00 and w > 0.4:
            return 3
        else:
            return -1
            
    def roll_level(self, roll):
        """
        Determines the level of Roll based on the provided Roll parameter.
    
        Parameters:
            roll (dict): Dictionary containing the Roll parameter.
                - 'sigma' (float): The Roll parameter value.
    
        Returns:
            level (int): The level of Roll.
    
        Description:
            This function determines the level of Roll based on the provided Roll parameter value (sigma).
            
            The levels are determined as follows:
            - If the inverse of the Roll parameter (1/sigma) is less than 1.4, the Roll level is set to 1.
            - If the inverse of the Roll parameter (1/sigma) is less than 3, the Roll level is set to 2.
            - If the inverse of the Roll parameter (1/sigma) is less than 10, the Roll level is set to 3.
            - Otherwise, if none of the above conditions are met, the Roll level is set to -1.
    
            The calculated Roll level is returned as the output.
        """
        if roll == False:
            return -1
        
        sigma = roll['sigma']
        
        if 1/sigma < 1.4:
            return 1
        elif 1/sigma < 3:
            return 2
        elif 1/sigma < 10:
            return 3
        else:
            return -1
    
    def spiral_level(self, spiral):
        """
        Determines the level of Spiral based on the provided Spiral parameter.
    
        Parameters:
            spiral (dict): Dictionary containing the Spiral parameter.
                - 'sigma' (float): The Spiral parameter value.
    
        Returns:
            level (int): The level of Spiral.
    
        Description:
            This function determines the level of Spiral based on the provided Spiral parameter value.
            
            The levels are determined as follows:
            - If the calculated value of t2 (log base 10 of 0.5 divided by the Spiral parameter sigma) is greater than 20, the Spiral level is set to 1.
            - If the calculated value of t2 is greater than 12, the Spiral level is set to 2.
            - If the calculated value of t2 is greater than 4, the Spiral level is set to 3.
            - Otherwise, if none of the above conditions are met, the Spiral level is set to -1.
    
            The calculated Spiral level is returned as the output.
        """
        if spiral == False:
            return -1
        
        t2 = math.log10(0.5) / spiral['sigma']
        
        if t2 > 20:
            return 1
        elif t2 > 12:
            return 2
        elif t2 > 4:
            return 3
        else:
            return -1

    def dynamic_stability(self, altitude, mass, Ix, Iy, Iz, u0 = 20,
                            CDu = 0, CLu = 0, Cmu = 0):
        """
        Calculates the dynamic stability characteristics of an aircraft.
    
        Args:
            altitude (float): The altitude in meters.
            mass (float): The mass of the aircraft in kilograms.
            Ix (float): The moment of inertia around the x-axis.
            Iy (float): The moment of inertia around the y-axis.
            Iz (float): The moment of inertia around the z-axis.
            u0 (float, optional): The reference velocity in meters per second. Defaults to 20.
            CDu (float, optional): The drag coefficient due to the control deflection u. Defaults to 0.
            CLu (float, optional): The lift coefficient due to the control deflection u. Defaults to 0.
            Cmu (float, optional): The pitching moment coefficient due to the control deflection u. Defaults to 0.
    
        Returns:
            tuple: A tuple containing the following stability characteristics:
                - short_period (dict): A dictionary containing the characteristics of the short period mode.
                    - 't12': Time to reach half amplitude.
                    - 'T': Period of the mode.
                    - 'N12': Number of cycles to half amplitude.
                    - 'zeta': Damping ratio of the mode.
                    - 'w': Undamped natural frequency.
                    - 'n/alfa': Acceleration sensitivity.
                    - 'CAP': Control acceleration parameter.
                - pughoid (dict): A dictionary containing the characteristics of the phugoid mode.
                    - 't12': Time to reach half amplitude.
                    - 'T': Period of the mode.
                    - 'N12': Number of cycles to half amplitude.
                    - 'zeta': Damping ratio of the mode.
                    - 'w': Undamped natural frequency.
                - dutch_roll (dict): A dictionary containing the characteristics of the Dutch roll mode.
                    - 't12': Time to reach half amplitude.
                    - 'zeta': Damping ratio of the mode.
                    - 'sigma': Damping ratio.
                    - 'w': Undamped natural frequency.
                    - 'T': Period of the mode.
                - roll (dict): A dictionary containing the characteristics of the roll mode.
                    - 't12': Time to reach half amplitude.
                    - 'sigma': Damping ratio.
                - spiral (dict): A dictionary containing the characteristics of the spiral mode.
                    - 't12': Time to reach half amplitude.
                    - 'sigma': Damping ratio.
        """
        
        density = self.engine.get_density(altitude)
        chord = self.get_surface("wing").get_smc()
        Bref = self.get_Bref()
        

        CD0 = self.get_CD0()
        CL0 = self.get_CL0()
        CDa = self.get_CDa()*180/np.pi  # 1/Deg --> 1/rad
        CLa = self.get_CLa()*180/np.pi  # 1/Deg --> 1/rad
        Cma = self.get_Cma()*180/np.pi  # 1/Deg --> 1/rad
        Cmq = self.get_Cmq()
        
        CYB = self.get_CyB() 
        CYp = self.get_Cyp()
        CYr = self.get_Cyr()
        
        ClB = self.get_ClB()
        Clp = self.get_Clp()
        Clr = self.get_Clr()
        
        CnB = self.get_CnB()
        Cnp = self.get_Cnp()
        Cnr = self.get_Cnr()
        
        CLda = 0#1.918 # 0 # 1.918
        Cmda = 0#-3.36 # 0 # -3.36
        
        Q = (1/2)*density*u0**2 
        S = self.get_Sref()
        
        # Longitudinal derivatives 
        Xu = -(CDu+2*CD0)*Q*S/(mass*u0)
        Zu = -(CLu+2*CL0)*Q*S/(mass*u0)
        Mu = Cmu*(Q*S*chord)/(u0*Iy)
        Xw = -(CDa-CL0)*Q*S/(mass*u0)
        Zw = -(CLa+CD0)*Q*S/(mass*u0)
        Mw = Cma*(Q*S*chord)/(u0*Iy)
        # Zdw = CLda*chord/(2*u0)*(Q*S)/(u0*mass)
        Mdw = Cmda*chord/(2*u0)*(Q*S*chord)/(u0*Iy)
        # Zq = CLq*(chord/(2*u0))*(Q*S/mass)
        Mq = Cmq*(chord/(2*u0))*((Q*S*chord)/Iy)
        # Ma = u0*Mw
        # Mda = u0*Mdw
        Za = u0*Zw
        
        # Latero-Directional derivatives
        YB = Q*S*CYB/mass
        Yp = Q*S*Bref*CYp/(2*mass*u0)
        Yr = Q*S*Bref*CYr/(2*mass*u0)
        LB = Q*S*Bref*ClB/Ix
        Lp = Q*S*Bref**2*Clp/(2*Ix*u0)
        Lr = Q*S*Bref**2*Clr/(2*Ix*u0)
        NB = Q*S*Bref*CnB/Iz
        Np = Q*S*Bref**2*Cnp/(2*Iz*u0)
        Nr = Q*S*Bref**2*Cnr/(2*Iz*u0)
        
        g = 9.8066  
         
        # Longitudinal 
        A_long = np.array([[Xu, Xw, 0., -g],
                           [Zu, Zw, u0, 0.],
                           [Mu+Mdw*Zu, Mw+Mdw*Zw, Mq+Mdw*u0, 0.],
                           [0., 0., 1, 0.]])

        eig_long = np.linalg.eig(A_long)[0]
        
        # print('\n\n', eig_long, '\n\n')
            
        ind_sp = abs(eig_long.real).argmax() 
        ind_pg = abs(eig_long.real).argmin() 
        
        eig_sp = [abs(eig_long[ind_sp].real), abs(eig_long[ind_sp].imag)]
        eig_pg = [abs(eig_long[ind_pg].real), abs(eig_long[ind_pg].imag)]
        
        short_period = {}
        short_period['t12'] = 0.69/abs(eig_sp[0])
        short_period['T'] = 2*np.pi/eig_sp[1]
        short_period['N12'] = short_period['t12']/short_period['T']
        short_period['zeta'] = (1/((eig_sp[1]/eig_sp[0])**2 + 1))**(1/2)
        short_period['w'] = eig_sp[0]/(short_period['zeta'])
        short_period['n/alfa'] = -Za/g
        short_period['CAP'] = short_period['w']**2/(-Za/g)
        
        pughoid = {}
        pughoid['t12'] = 0.69/abs(eig_pg[0])
        pughoid['T'] = 2*np.pi/eig_pg[1]
        pughoid['N12'] = pughoid['t12']/pughoid['T']
        pughoid['zeta'] = (1/((eig_pg[1]/eig_pg[0])**2 + 1))**(1/2)
        pughoid['w'] = eig_pg[0]/(-pughoid['zeta'])
        
        # Latero - Directional
        A_latd = np.array([[YB/u0, Yp/u0, -(1-Yr/u0), g/u0],
                           [LB, Lp, Lr, 0],
                           [NB, Np, Nr, 0],
                           [0, 1, 0, 0]])
        
        eig_latd = np.linalg.eig(A_latd)[0]
        # print('\n\n\n', eig_latd, '\n\n\n')
        Roll_Spiral = []
        
        for eig in eig_latd:
            if eig.imag == 0:
                Roll_Spiral.append(eig.real)
            else:
                eig_DR = [eig.real, abs(eig.imag)]
                
        try:
            eig_roll = [min(Roll_Spiral),0]
            eig_spiral = [max(Roll_Spiral),0]
                
            dutch_roll = {}
            dutch_roll['t12'] = 0.69/abs(eig_DR[0])
            dutch_roll['zeta'] = (1/((eig_DR[1]/eig_DR[0])**2 + 1))**(1/2)
            dutch_roll['sigma'] = -eig_DR[0]*(2*u0/chord)
            dutch_roll['w'] = eig_DR[0]/(-dutch_roll['zeta'])
            dutch_roll['T'] = 2*np.pi/dutch_roll['w']
            
            roll = {}
            roll['t12'] = 0.69/abs(eig_roll[0])
            roll['sigma'] = -eig_roll[0]*(2*u0/chord)
            
            
            spiral = {}
            spiral['t12'] = 0.69/abs(eig_spiral[0])
            spiral['sigma'] = -eig_spiral[0]*(2*u0/chord)
            
            # print('\n\\n',eig_latd,'\n\n')
        except:
            dutch_roll = False 
            roll = False
            spiral = False
        
        return short_period, pughoid, dutch_roll, roll, spiral
        
    def get_static(self, name='all'):
        """
        Returns the static stability values of the aircraft.

        Parameters
        ----------
        name : str
            The name of the stability value to retrieve. Default is 'all'.

        Returns
        -------
        tuple or boolean
            The static stability values of the aircraft. If 'name' is given, returns the specified value as a boolean.
            If 'all' is given, returns a tuple containing all static stability values.
        """

        if name == 'all':
            return self.__lon_static_stability, self.__lat_static_stability, self.__dir_static_stability
        elif name == 'longitudinal':
            return self.__lon_static_stability
        elif name == 'lateral':
            return self.__lat_static_stability
        elif name == 'directional':
            return self.__dir_static_stability

    def get_dynamic(self, name='all'):
        """
        Returns the dynamic stability values of the aircraft.

        Parameters
        ----------
        name : str
            The name of the stability value to retrieve. Default is 'all'.

        Returns
        -------
        tuple or boolean
            The dynamic stability values of the aircraft. If 'name' is given, returns the specified value as a boolean.
            If 'all' is given, returns a tuple containing all dynamic stability values.
        """

        if name == 'all':
            return self.__lon_dynamic_stability, self.__lat_dynamic_stability, self.__dir_dynamic_stability
        elif name == 'longitudinal':
            return self.__lon_dynamic_stability
        elif name == 'lateral':
            return self.__lat_dynamic_stability
        elif name == 'directional':
            return self.__dir_dynamic_stability
        
