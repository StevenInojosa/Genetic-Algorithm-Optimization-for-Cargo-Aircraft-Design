#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 12:44:43 2023

@author: stevenstefanoinojosasiso
"""
import numpy as np

class Structures:
    
    def rib_inertia(self, chord, Ix = 541e-9, Iy = 23457e-9 , Iz = 22920e-9,
                    xcm = 142e-3, ycm = 1.5e-3, zcm = 27e-3, 
                    og_chord = 400e-3, mass = 2.83e-3):
        """
        Inertia [g/mm^3] in reference to the center of mass. eading edge.
        """
        ra = chord/og_chord
        return [ra**4*Ix, ra**4*Iy, ra**4*Iz], [ra*xcm, ra*ycm, ra*zcm], ra**2*mass
        
    def parallel_axis(self, Icm, mass, d):
        Ix = Icm[0] + mass* (d[1]**2 + d[2]**2)
        Iy = Icm[1] + mass* (d[0]**2 + d[2]**2)
        Iz = Icm[2] + mass* (d[0]**2 + d[1]**2)
        
        return [Ix, Iy, Iz]
    
    def inertia_surface(self, name, separation=0.08, R=0.02, e=0.001, density_beam=1730):
        Xref = self.get_Xref()
        Yref = self.get_Yref()
        Zref = self.get_Zref()
    
        surface = self.get_surface(name)
        translate = surface.get_translate()
    
        Xle_list = []
        Yle_list = []
        Zle_list = []
        Chord_list = []
    
        for section in surface.get_sections():
            Xle_list.append(section.get_Xle() + translate[0])
            Yle_list.append(section.get_Yle() + translate[1])
            Zle_list.append(section.get_Zle() + translate[2])
            Chord_list.append(section.get_Chord())
    
        Ix_ribs = 0
        Iy_ribs = 0
        Iz_ribs = 0
        Ix_beam = 0
        Iy_beam = 0
        Iz_beam = 0
        mass_ribs_total = 0
    
        for ii in range(len(Yle_list) - 1):
            x1 = Xle_list[ii]
            x2 = Xle_list[ii + 1]
            y1 = Yle_list[ii]
            y2 = Yle_list[ii + 1]
            z1 = Zle_list[ii]
            z2 = Zle_list[ii + 1]
    
            Chord1 = Chord_list[ii]
            Chord2 = Chord_list[ii + 1]
    
            n = int(np.ceil((y2 - y1) / separation))
    
            # Compute the Ribs inertia and mass
            for yle in np.linspace(y1, y2, n):
                xle = (x2 - x1) / (y2 - y1) * (yle - y1) + x1
                zle = (z2 - z1) / (y2 - y1) * (yle - y1) + z1
                chord = (Chord2 - Chord1) / (y2 - y1) * (yle - y1) + Chord1
    
                # Compute the inertia and mass on the center of mass of the rib.
                Icm, rcm, mass_rib = self.rib_inertia(chord)
    
                # Compute the inertia in reference to the Aerodynamic center.
                d_rib = [Xref - xle - rcm[0],
                         Yref - yle - rcm[1],
                         Zref - zle - rcm[2]]
    
                Inertia_rib = self.parallel_axis(Icm, mass_rib, d_rib)
    
                Ix_ribs += Inertia_rib[0]
                Iy_ribs += Inertia_rib[1]
                Iz_ribs += Inertia_rib[2]
    
                mass_ribs_total += mass_rib
    
    
        # Consider the symmetry
        Ix_ribs *= 2
        Iy_ribs *= 2
        Iz_ribs *= 2
        mass_ribs_total *= 2
    
        # Add the inertia of the beam
        r = R - e
        H = 2 * Yle_list[-1]
        mass_beam = density_beam * np.pi * (R**2 - r**2) * H
        Ix_beam = Iz_beam = mass_beam * (1 / 4) * (R**2 + r**2 + H**2 / 3)
        Iy_beam = mass_beam * (1 / 2) * (R**2 + r**2)
    
        d_beam = [Xref - Xle_list[0] - Chord_list[0] / 4,
                  Yref - Yle_list[0],
                  Zref - Zle_list[0]]
        
        Icm_beam = [Ix_beam, Iy_beam, Iz_beam]
        
        # print(Ix_ribs, Iy_ribs, Iz_ribs)
        # print(mass_ribs_total)
        # print('\n')
        # print(Ix_beam, Iy_beam, Iz_beam)
        # print(mass_beam)
        # print('\n')

        
        Ix_beam, Iy_beam, Iz_beam = self.parallel_axis(Icm_beam, mass_beam, d_beam)
    
        Ix = Ix_beam + Ix_ribs
        Iy = Iy_beam + Iy_ribs
        Iz = Iz_beam + Iz_ribs
    
        mass = mass_ribs_total + mass_beam
        
        # print(Ix_ribs, Iy_ribs, Iz_ribs)
        # print(mass_ribs_total)
        # print('\n')
        # print(Ix_beam, Iy_beam, Iz_beam)
        # print(mass_beam)
        # print('\n')
        # print(mass, d_beam, Xref, )
        # print('\n')

        return [Ix, Iy, Iz], mass

    
    def inertia_engine(self, xle, yle, zle, mass = 0.750, R = 0.12):
        Xref = self.get_Xref()
        Yref = self.get_Yref()
        Zref = self.get_Zref()
        
        Icm = [2/5 * mass * R**2, 2/5 * mass * R**2, 2/5 * mass * R**2]
        d = [Xref - (xle + R/2),
             Yref - yle,
             Zref - zle]
        
        Inertia = self.parallel_axis(Icm, mass, d)
        
        return Inertia, mass

    def inertia_main_gear(self, mass = 0.175, separation = 0.2, R = 0.2, name = 'wing'):
        
        Xref = self.get_Xref()
        Yref = self.get_Yref()
        Zref = self.get_Zref()
        
        surface = self.get_surface(name)
        surface.get_sections()[0]
        translate = surface.get_translate()
            
        section_root = surface.get_sections()[0]
        
        Xler = section_root.get_Xle() + translate[0]
        Yer = section_root.get_Yle() + translate[1]
        Zler = section_root.get_Zle() + translate[2]
        
        chord = chord =  section_root.get_Chord()
        
        Icm = [(1/4)*mass*R**2 - mass * (4*R/(3*np.pi))**2 ,
               (1/2)*mass*R**2 - mass * (4*R/(3*np.pi))**2,
               (1/4)*mass*R**2]
        
        d = [Xref - Xler - chord/2,
             Yref - separation,
             Zref - Zler]
        
        Ix, Iy, Iz = self.parallel_axis(Icm, mass, d)
        
        # Consider the two gears
        Inertia = [2*Ix, 2*Iy, 2*Iz]
        
        mass = 2*mass
        
        return Inertia, mass

    
    def inertia_nose_gear(self, xle, yle, zle, mass = 0.150, R = 1/100, H = 10/100):
        Icm = [(1/4)*mass*(R**2 + H**2/3),
               (1/4)*mass*(R**2 + H**2/3),
               (1/2)*mass*R**2]
        
        Xref = self.get_Xref()
        Yref = self.get_Yref()
        Zref = self.get_Zref()
        
        d = [Xref - (xle + R/2),
             Yref - yle,
             Zref - (zle - H/2)]
        
        Inertia = self.parallel_axis(Icm, mass, d)
        
        return Inertia, mass
    
    def inertia_fuselage(self):
        pass
    
    def inertia_payload(self, payload, payload_density = 7800, name = 'wing',
                        ly = 8/100, lz = 8/100):
        
        Xref = self.get_Xref()
        Yref = self.get_Yref()
        Zref = self.get_Zref()
        
        surface = self.get_surface(name)
        surface.get_sections()[0]
        translate = surface.get_translate()
            
        section_root = surface.get_sections()[0]
        
        Xler = section_root.get_Xle() + translate[0]
        Yler = section_root.get_Yle() + translate[1]
        Zler = section_root.get_Zle() + translate[2]
        
        lx = payload/(payload_density*ly*lz)
        
        mass =  payload
        
        Icm  = [(mass/12) * (ly**2 + lz**2),
                (mass/12) * (lx**2 + lz**2),
                (mass/12) * (lx**2 + ly**2)]
        
        d = [Xref - Xler - lx/2,
             Yref,
             Zref - Zler]
        
        Inertia = self.parallel_axis(Icm, mass, d)

        return Inertia, mass
        
        
        
                

        
        
        
        
        
    