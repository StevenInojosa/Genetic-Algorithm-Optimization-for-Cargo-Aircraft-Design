#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

@author: stevenstefanoinojosasiso
"""

from AVL_python import *
import numpy as np
import sympy as sp
import warnings 
import copy
import math

class Aerodynamics:
    
    def __init__(self, Aircraft):
        pass
         
    def coefficients(self, 
                     Alpha  = [aa for aa in np.arange(-8,12,2)],
                     Roll   = [-6*np.pi/180],
                     Pitch  = [aa for aa in np.arange(-6,6,2)*np.pi/180],
                     Yaw    = [-6*np.pi/180],
                     Beta   = [-6] ):
        """
        Computes the aerodynamic coefficients of the aircraft at the specified flight conditions.
    
        Args:
            Alpha (list): List of angles of attack [deg] to consider. Default is [-8, -6, -4, -2, 0, 2, 4, 6, 8, 10].
            Roll (list): List of roll angles [rad] to consider. Default is [-6*pi/180].
            Pitch (list): List of pitch angles [rad] to consider. Default is [-6*pi/180, -4*pi/180, -2*pi/180, 0, 2*pi/180, 4*pi/180].
            Yaw (list): List of yaw angles [rad] to consider. Default is [-6*pi/180].
            Beta (list): List of sideslip angles [deg] to consider. Default is [-6].
    
        Returns:
            None
    
        The computed aerodynamic coefficients are stored in the instance variables of the class:
            self.__CLa (float): Lift slope [1/deg].
            self.__CL0 (float): Lift coefficient at zero angle of attack.
            self.__Cma (float): Moment coefficient about the aerodynamic center [1/deg].
            self.__Cm0 (float): Moment coefficient about the aerodynamic center at zero angle of attack.
            self.__K (float): Drag coefficient k factor.
            self.__CD0 (float): Drag coefficient at zero lift.
            self.__CDa (float): Lift-dependent drag coefficient slope [1/deg].
            self.__Cyp (float): Roll damping coefficient [1/rad].
            self.__Cnp (float): Yaw damping coefficient [1/rad].
            self.__Clp (float): Roll moment due to roll rate [1/rad].
            self.__Czq (float): Pitch damping coefficient [1/rad].
            self.__Cmq (float): Pitch moment due to pitch rate [1/rad].
            self.__Cyr (float): Side force due to yaw rate [1/rad].
            self.__Cnr (float): Yaw moment due to yaw rate [1/rad].
            self.__Clr (float): Roll moment due to yaw rate [1/rad].
            self.__CyB (float): Side force due to sideslip angle [1/rad].
            self.__CnB (float): Yaw moment due to sideslip angle [1/rad].
            self.__ClB (float): Roll moment due to sideslip angle [1/rad].
            self.__CLmin (float): Minimum lift coefficient.
    
        In addition, the following pandas DataFrames are created:
            self.alpha_pd (DataFrame): Table of aerodynamic coefficients vs. angle of attack.
            self.roll_pd (DataFrame): Table of aerodynamic coefficients vs. roll angle.
            self.pitch_pd (DataFrame): Table of aerodynamic coefficients vs. pitch angle.
            self.yaw_pd (DataFrame): Table of aerodynamic coefficients vs. yaw angle.
            self.beta_pd (DataFrame): Table of aerodynamic coefficients vs. sideslip angle.
        """
    
        
        avl_turpial = AVL_python(Alpha, Roll, Pitch, Yaw, Beta, self)
        
        avl_turpial.geometry_creator()
        avl_turpial.cmd()
        
        out, cmd, alpha_table, roll_table, pitch_table, yaw_table, beta_table = avl_turpial.run()
        
        # Alpha Table
        self.__CLa, self.__CL0 = np.polyfit(alpha_table['Deg'],
                                            alpha_table['CL'],
                                            1)
        self.__Cma, self.__Cm0 = np.polyfit(alpha_table['Deg'],
                                            alpha_table['Cm'],
                                            1)
        self.__K, self.__CD0 = np.polyfit(alpha_table['CL^2'],
                                          alpha_table['CDtot'],
                                          1)
        self.__CDa, CD_a0 = np.polyfit(alpha_table['Deg'],
                                       alpha_table['CDtot'],
                                       1)
        
        # Roll Table
        self.__Cyp, Cy_p0 = np.polyfit(roll_table['rad'],
                                       roll_table['Cy'],
                                       1)
        self.__Cnp, Cn_p0 = np.polyfit(roll_table['rad'],
                                       roll_table['Cn'],
                                       1)
        self.__Clp, Cl_p0 = np.polyfit(roll_table['rad'],
                                       roll_table['Cl'],
                                       1)
        
        # Pitch Table
        self.__Czq, Cz_q0 = np.polyfit(pitch_table['rad'],
                                       pitch_table['Cz'],
                                       1)
        self.__Cmq, Cm_q0 = np.polyfit(pitch_table['rad'],
                                       pitch_table['Cm'],
                                       1)
        
        # Yaw Table
        self.__Cyr, Cy_r0 = np.polyfit(yaw_table['rad'],
                                       yaw_table['Cy'],
                                       1)
        self.__Cnr, Cn_r0 = np.polyfit(yaw_table['rad'],
                                       yaw_table['Cn'],
                                       1)
        self.__Clr, Cl_r0 = np.polyfit(yaw_table['rad'],
                                       yaw_table['Cl'],
                                       1)
        
        # Beta Table
        self.__CyB, Cy_B0 = np.polyfit(beta_table['rad'],
                                       beta_table['Cy'],
                                       1)
        self.__CnB, Cn_B0 = np.polyfit(beta_table['rad'],
                                       beta_table['Cn'],
                                       1)
        self.__ClB, Cl_B0 = np.polyfit(beta_table['rad'],
                                       beta_table['Cl'],
                                       1)

        self.alpha_table = alpha_table
        self.roll_table = roll_table
        self.pitch_table = pitch_table
        self.yaw_table = yaw_table
        self.beta_table = beta_table 
           
        self.alpha_pd = pd.DataFrame(
            {
                "Deg": alpha_table['Deg'],
                "Cm": alpha_table['Cm'],
                "CL": alpha_table['CL'],
                "CD": alpha_table['CDtot'],
                "CLff": alpha_table['CLff'],
                "CL^2": alpha_table['CL^2'],
                "CLff^2": alpha_table['CLff^2']
            }
        )
            
        self.roll_pd = pd.DataFrame(
            {
                "Deg": roll_table['Deg'],
                "rad": roll_table['rad'],
                "Cy": roll_table['Cy'],
                "Cn": roll_table['Cn'],
                "Cl": roll_table['Cl']
            }
        )
            
        self.pitch_pd = pd.DataFrame(
            {
                "Deg": pitch_table['Deg'],
                "rad": pitch_table['rad'],
                "Cz": pitch_table['Cz'],
                "Cm": pitch_table['Cm']
            }
        )
            
        self.yaw_pd = pd.DataFrame(
            {
                "Deg": yaw_table['Deg'],
                "rad": yaw_table['rad'],
                "Cy": yaw_table['Cy'],
                "Cn": yaw_table['Cn'],
                "Cl": yaw_table['Cl']
            }
        )
            
        self.beta_pd = pd.DataFrame(
            {
                "Deg": beta_table['Deg'],
                "Cy": beta_table['Cy'],
                "Cn": beta_table['Cn'],
                "Cl": beta_table['Cl']
            }
        )

        self.__CLmin = 0
        
    def CDp_streamlines(self, name, V = 10, ro = 1.1226, visc = 1.91*10**(-5) ):
        """
        Calculates the parasite drag coefficient for a given surface section at a specific velocity.
        
        Parameters:
            name (str): Name of the surface section.
            V (float): Velocity of the airflow in meters per second. Defaults to 10.
            ro (float): Density of the air in kilograms per cubic meter. Defaults to 1.1226.
            visc (float): Viscosity of the air in pascal seconds. Defaults to 1.91*10^(-5).
            
        Returns:
            CDp (float): Profile drag coefficient for the surface section at the given velocity.
        """
        # Get surface important data
        surface = self.get_surface(name)
        t = surface.get_sections()[0].get_airfoil('t') 
        c = surface.get_mac()
        
        Re = float(ro*V*c/visc)
        Cf = 0.455/(np.log10(Re)**2.58)
        
        return 2*Cf*(1+2.2*(t/c)+60*(t/c)**4)
    
    def CDp_interaction(self, A1, A2, t, c):
        """
        Calculates the drag coefficient due to interactions.
    
        Parameters:
            A1 (float): Planform area 1 [m^2].
            A2 (float): Planform area 2 [m^2].
            t (float): Maximum thickness of airfoil [m].
            c (float): Mean aerodynamic chord of wing [m].
    
        Returns:
            CDp (float): Drag coefficient due to wing-to-wing interaction.
        """
        return (0.75*(t/c)-0.0003/(t/c)**2)*(t*c)**2/(A1+A2)
        
    def CDp_Aircraft(self, streamlines = ['wing', 'horizontal', 'vertical']): 
        '''
        Computes the CDp for the entire aircraft, including streamlines and interactions, 
        based on the given wing and horizontal surfaces.
        
        Args:
        - streamlines: a list of the names of the surfaces to be considered for streamlines. 
                       The default is ['wing', 'horizontal', 'vertical'].
                       
        Returns:
        - A dictionary with the CDp values for the aircraft, the wing, the horizontal surface, the fuselage, 
          the wing-fuselage interaction, the horizontal-vertical tail interaction, and various components 
          such as wheels, motors, and propellers.
        '''
        Sref = self.get_Sref()
        self.__CDp = {'Aircraft':0}
        
        # Compute the CDp for streamlines shapes
        for name in streamlines: 
            try:
                self.__CDp[name] = self.CDp_streamlines(name)*self.get_surface(name).get_S()/Sref
            except:
                vertical = self.get_surface(name)
                SMCv = ( vertical.get_sections()[0].get_Chord() + 
                         vertical.get_sections()[1].get_Chord() ) / 2
                bv = vertical.get_sections()[1].get_Yle()
                Sv = SMCv * bv
                self.__CDp[name] = self.CDp_streamlines(name)*Sv/Sref
    
        Sh = self.get_surface('horizontal').get_S()
        th = self.get_surface('horizontal').get_sections()[0].get_airfoil('t')
        Ch = self.get_surface('horizontal').get_mac()
        
        tw = self.get_surface('wing').get_sections()[0].get_airfoil('t')
        Cw = self.get_surface('wing').get_mac()
        
        
        Af = 0.08*0.08 # Change if you wanted to consider fuselaje
        
        if Af != 0:
            self.__CDp['Fuselage'] = 0.8*(Af/Sref)
        else:
            self.__CDp['Fuselage'] = 0
        
        # Compute Landing Gear
        self.__CDp['Tren Aterrizaje'] = 0.00195/Sref
        
        
        
        # Compute the CDp for interactions 
        self.__CDp['Wing/Fuselage'] = self.CDp_interaction(Sref,Af,tw,Cw)/Sref;
        
        try:
            Sv = self.get_surface('vertical').get_S()
            Nv = 2
            self.__CDp['Horizontal/Vertical'] = self.CDp_interaction(Sh,Sv,th,Ch)/Sref*Nv
        except:
            print("There is no vertical stabilizer")
              
        self.__CDp['Rueda/T. Principal']= 0.0079/Sref                 
        self.__CDp['Rueda/T. Nariz'] = 0.0021/Sref                     
        self.__CDp['Motor/Fuselage'] = 0.0011/Sref                            
        self.__CDp['Motor/Helice'] = 0.0001/Sref      
        
        for name, CDp in copy.deepcopy(self.__CDp).items():
            self.__CDp['Aircraft'] += CDp       
        
    def get_CDp(self, name='all'):
        """
        Get the coefficient of parasite drag for a specific name, or all if no name is given.
    
        Parameters:
        -----------
        name : str, optional
            The name of the coefficient of drag to get. If 'all', returns all coefficients.
            Defaults to 'all'.
    
        Returns:
        --------
        CDp : float or dict
            The coefficient of drag value for the specified name, or a dictionary containing
            all coefficient of drag values if 'all' is given for the name.
        """
        if name == 'all':
            return self.__CDp
        else:
            return self.__CDp[name]
        
    def CLmax(self, AR, RT, b, reltrap, admissible_error = 1e-5):
        """
        Calculate the CLmax using the Lifting-line theory
            
        Returns:
            CLmax (float): Maximum Lift coefficient [-]
      
        """
        CLmax = [float('inf'), -float('inf')]
        if reltrap < 50:
            RT = 1
    
        Cla = 0.0957 * 180 / math.pi    # Pendiente del coeficiente de sustentacion en funcion del angulo de ataque (perfil 1/rad)
        Clmax = 2.2375                  # Coeficiente de sustentacion maxima (perfil)

    
        n = 20
        omega = 0.0 * math.pi / 180  # alabeo en radianes
        
        while abs(CLmax[-1] - CLmax[-2]) > admissible_error:
            o = (90 * math.pi / 180) / n     # initial angle for the Fourier series
            o1 = o                           # initial angle for the Fourier series (to be used in calculations)
            
            a = np.zeros((n, n))
            
            # Filling the matrix An and Bn
            for n2 in range(n):
                v = 1
                
                for n1 in range(n):
                    if v < 2*n:
                        # Eqs 10 and 46.
                        a[n2, n1] = (((2 * AR * (1 + RT)) / (Cla * (1 - (1 - RT) * (abs(math.cos(o1)))))) + v / (math.sin(o1))) * math.sin(v * o1)
                    v += 2
                o1 += o
                
            # vectores que multiplicarán la matrices
            v1 = np.ones((n, 1))
            v2p = np.arange(o, o * (n + 1), o)
            v21 = np.cos(v2p)
            v2 = v21.T
            
            # vectores resultantes de operar las matrices (incluye la inversa de la matriz)
            af = np.linalg.inv(a) @ v1
            bf = np.linalg.inv(a) @ v2
                
            CLa = math.pi * AR * af[0]                                 # Eq. 16
            theta_max = math.acos(1 - RT)                              # Eq. 48
            CL_Clmax = math.pi * math.sin(theta_max) / (2 * (1 + RT))  # Eq. 47
            
            KLs = 1 + (0.0042 * AR - 0.068) * (1 + 2.3 * CLa * omega / Clmax) # Eq. 52
            
            CL_Clmax_omega0  = 0
            for ii in np.arange(2, n, 2):
                CL_Clmax_omega0 += (-1)**(ii+1) * (af[ii] / af[0])
            CL_Clmax_omega0 = np.pi/4 * (1 + CL_Clmax_omega0)
            
            if omega == 0:
                CLmax.append( CL_Clmax_omega0 * KLs * (Clmax) )
            else:
                CLmax.append( CL_Clmax_omega0 * KLs * KL_Lambda * (Clmax - KLo*CLa*omega) )
                
            n *= 2

        self.__CLmax = CLmax[-1][0]
        
        return self.__CLmax
        
    def get_CLmin(self):
        """
        Returns:
            CLmin (float): Minimum Lift coefficient [-]
      
        """
        return self.__CLmin
    
    def get_CLmax(self):
        """
        Returns:
            CLmax (float): Maximum Lift coefficient [-]
      
        """
        return self.__CLmax
    
    def get_CLa(self):
        """
        Returns the lift coefficient slope with respect to angle of attack for the aircraft.
    
        Returns:
        -------
        float:
            Lift coefficient slope with respect to angle of attack [1/deg].
        """
        return self.__CLa
    
    def get_CL0(self):
        """
        Returns the lift coefficient when the angle of attack is zero for the aircraft.
    
        Returns:
        -------
        float:
            Lift coefficient at zero angle of attack.
        """
        return self.__CL0
    
    def get_Cma(self):
        """
        Returns the pitching moment coefficient slope with respect to angle of attack for the aircraft.
    
        Returns:
        -------
        float:
            Pitching moment coefficient slope with respect to angle of attack [1/deg].
        """
        return self.__Cma
    
    def get_Cm0(self):
        """
        Returns the pitching moment coefficient when the angle of attack is zero for the aircraft.
    
        Returns:
        -------
        float:
            Pitching moment coefficient at zero angle of attack.
        """
        return self.__Cm0
    
    def get_K(self):
        """
        Returns the Drag coefficient k factor for the aircraft.
    
        Returns:
        -------
        float:
            Drag coefficient k factor.
        """
        return self.__K
    
    def get_CD0(self):
        """
        Returns the zero lift drag coefficient for the aircraft.
    
        Returns:
        -------
        float:
            Zero lift drag coefficient.
        """
        return self.__CD0
    
    def get_CDa(self):
        """
        Returns the drag coefficient slope with respect to angle of attack for the aircraft.
    
        Returns:
        -------
        float:
            Drag coefficient slope with respect to angle of attack [1/deg].
        """
        return self.__CDa
    
    def get_Cyp(self):
        """
        Returns the roll damping derivative with respect to yaw rate for the aircraft.
    
        Returns:
        -------
        float:
            Roll damping derivative with respect to yaw rate [1/rad].
        """
        return self.__Cyp
    
    def get_Cnp(self):
        """
        Returns the yaw damping derivative with respect to yaw rate for the aircraft.
    
        Returns:
        -------
        float:
            Yaw damping derivative with respect to yaw rate [1/rad].
        """
        return self.__Cnp
    
    def get_Clp(self):
        """
        Returns the roll damping derivative with respect to roll rate for the aircraft.
    
        Returns:
        -------
        float:
            Roll damping derivative with respect to roll rate [1/rad].
        """
        return self.__Clp
    
    def get_Czq(self):
        """
        Returns the pitch damping derivative with respect to pitch rate for the aircraft.
    
        Returns:
        -------
        float:
            Pitch damping derivative with respect to pitch rate [1/rad].
        """
        return self.__Czq
    
    def get_Cmq(self):
        """
        Returns the pitch damping derivative with respect to pitch rate for the aircraft.
    
        Returns:
        -------
        float:
            Pitch damping derivative with respect to pitch rate [1/rad].
        """
        return self.__Cmq
    
    def get_Cyr(self):
        """
        Returns the roll damping derivative with respect to yaw rate for the aircraft.
    
        Returns:
        -------
        float:
            Roll damping derivative with respect to yaw rate [1/rad].
        """
        return self.__Cyr
    
    def get_Cnr(self):
        """
        Returns the yaw damping derivative with respect to yaw rate for the aircraft.
    
        Returns:
        -------
        float:
            Yaw damping derivative with respect to yaw rate [1/rad].
        """
        return self.__Cnr
    
    def get_Clr(self):
        """
        Returns the roll damping derivative with respect to roll rate for the aircraft.
    
        Returns:
        -------
        float:
            Roll damping derivative with respect to roll rate [1/rad].
        """
        return self.__Clr
    
    def get_CyB(self):
        """
        Returns the side force coefficient slope with respect to sideslip angle for the aircraft.
    
        Returns:
        -------
        float:
            Side force coefficient slope with respect to sideslip angle [1/rad].
        """
        return self.__CyB
    
    def get_CnB(self):
        """
        Returns the yawing moment coefficient slope with respect to sideslip angle for the aircraft.
    
        Returns:
        -------
        float:
            Yawing moment coefficient slope with respect to sideslip angle [1/rad].
        """
        return self.__CnB
    
    def get_ClB(self):
        """
        Returns the rolling moment coefficient slope with respect to sideslip angle for the aircraft.
    
        Returns:
        -------
        float:
            Rolling moment coefficient slope with respect to sideslip angle [1/rad].
        """
        return self.__ClB
    
