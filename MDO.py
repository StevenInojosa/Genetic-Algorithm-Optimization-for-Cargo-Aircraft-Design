#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: stevenstefanoinojosasiso
"""
from Aircraft import *
from sympy import *
import numpy as np

import pygad
import time

# Record the start time
start_time = time.time()


def dim_ref(xle, yle, chord):
    """
    Compute the reference dimensions.

    Parameters:
        xle (list): x distance of the section [m].
        yle (list): x distance of the section [m].
        chord (list): chord of every section [m].
                
    Returns:
        Sref (float): Reference area.
        Cref (float): Mean aerodynamic chord.
        Bref (float): Wing span.
        Xref (float): X-coordinate of the reference point.
        Yref (float): Y-coordinate of the reference point.
        Zref (float): Z-coordinate of the reference point.
    """
    
    I_mac = 0
    I_xc_4 = 0
    Sref = 0

    for ii in range(len(chord) - 1):
        y = symbols('y')
        c0 = float(chord[ii])
        c1 = float(chord[ii + 1])
        x0 = float(xle[ii])
        x1 = float(xle[ii + 1])
        y0 = float(yle[ii])
        y1 = float(yle[ii + 1])
        
        Sref += (y1 - y0) * (c1 + c0)
        c_y = (c1 - c0) / (y1 - y0) * y + c0
        x_c4 = x0 + (x1 - x0) / (y1 - y0) * y + c_y / 4
        
        I_mac += integrate(c_y**2, (y, 0, y1 - y0))
        I_xc_4 += integrate(x_c4 * c_y, (y, 0, y1 - y0))
        
    Cref = (2 / Sref) * I_mac
    Bref = yle[-1] * 2
    Xref = (2 / Sref) * I_xc_4
    Yref = 0
    Zref = 0
    
    return Sref, Cref, Bref, Xref, Yref, Zref


def ADN_to_Aircraft(ADN): 
    """Create an Aircraft object from the given genetic code (ADN).

    Parameters:
        ADN (list): Genetic code representing the aircraft design.
            ADN = [Cr, RT, Y1, Y2, tx, tz, Ch, Yh, Zv] 

    Returns:
        Aircraft: An aircraft object with the specified geometry.

    Legend:
        Cr: Wing root chord.
        RT: Wing taper ratio.
        Y1: Y-coordinate of the wing second section.
        Y2: Y-coordinate of the wing third section.
        tx: Horizontal distance bewteen wing and horizontal stabilizer.
        tz: Vertical distance bewteen wing and horizontal stabilizer.
        Ch: Horizontal stabilizer chord.
        Yh: Y-coordinate of the horizontal stabilizer .
        Zv: Z-coordinate of the vertical stabilizer.
    """

    C0 = ADN[0]
    RT = ADN[1] 
    Y1 = ADN[2] 
    Y2 = ADN[2] + ADN[3]
    tx = ADN[4]
    tz = ADN[5]
    Ch = ADN[6] 
    Yh = ADN[7]
    Zv = ADN[8]
    
    # Create the wing (Surface)
    wing = Surface(translate=[0, 0, 0], dAinc=3)
    
    wing.add_section(Section(Xle=0, Yle=0, Zle=0, Chord=C0, Ainc=0, airfoil=["s1223", 0.121]))
    wing.add_section(Section(Xle=0, Yle=Y1, Zle=0, Chord=C0, Ainc=0, airfoil=["s1223", 0.121]))
    wing.add_section(Section(Xle=0, Yle=Y2, Zle=0, Chord=C0 * RT, Ainc=0, airfoil=["s1223", 0.121]))
    
    # Create the horizontal stabilizer (Surface)
    horizontal = Surface(translate=[tx + C0, 0, tz], dAinc=0)

    horizontal.add_section(Section(Xle=0, Yle=0, Zle=0, Chord=Ch, Ainc=0, airfoil=["naca0012", 0.12]))
    horizontal.add_section(Section(Xle=0, Yle=Yh, Zle=0, Chord=Ch, Ainc=0, airfoil=["naca0012", 0.12]))
        
    # Create the vertical stabilizer (Surface)
    vertical = Surface(translate=[tx + C0, 0, tz], dAinc=0)

    vertical.add_section(Section(Xle=0, Yle=Yh, Zle=-Zv, Chord=Ch, Ainc=0, airfoil=["naca0012", 0.12]))
    vertical.add_section(Section(Xle=0, Yle=Yh, Zle=Zv, Chord=Ch, Ainc=0, airfoil=["naca0012", 0.12]))
    
    # Compute the Ref. Dimentions  
    Sref, Cref, Bref, Xref, Yref, Zref = dim_ref(xle = [0, 0, 0], yle = [0, Y1, Y2], chord = [C0, C0, C0*RT])
    
    # Create the aircraft (Aircraft)
    aircraft = Aircraft("Turpial", Sref, Cref, Bref, Xref, Yref, Zref)

    aircraft.add_surface(name = "wing", surface = wing)
    aircraft.add_surface(name = "horizontal", surface = horizontal)
    aircraft.add_surface(name="vertical", surface=vertical)
    
    return aircraft

def check_geometry_2018(aircraft, l_engine):
    Bref = aircraft.get_Bref()
    Lref = aircraft.get_surface('horizontal').get_translate(0) + \
           aircraft.get_surface('horizontal').get_sections()[0].get_Chord() + \
           aircraft.get_surface('wing').get_sections()[0].get_Xle() + \
           max(-l_engine,0)
    
    # Check geometry conditions    
    Sd = Bref + Lref + abs(l_engine) 
    H = aircraft.get_surface('horizontal').get_sections()[0].get_Zle()
    
    return (Sd<=3.5) 
    
def score(ADN, ADN_idx):
    """
    Calculate the fitness score for an individual in a genetic algorithm using its ADN (Aircraft Design Network).

    The function evaluates the aircraft design's stability, aerodynamics, performance, and structure to compute the
    fitness score based on certain constraints and design criteria. It also considers the stability of short-period
    motion and pughoid motion to ensure a stable aircraft design. The geometry and design constraints are checked to
    ensure the design is within acceptable limits. The function returns the computed fitness score, which is used in
    the genetic algorithm to select and evolve individuals for better aircraft designs.
    
    Parameters:
        ADN (list): The Aircraft Design Network containing design parameters for an individual.
        ADN_idx (int): Index of the ADN (used in the genetic algorithm).

    Returns:
        score (float): The computed fitness score for the individual's aircraft design.
    """
    # Required Data
    NR = 105
    NRmax = 185
    Sgmax_landing = 120
    l_engine = -0.12
    altitude = 900
    
    # Create aircraft.
    aircraft = ADN_to_Aircraft(ADN)
    
    Bref = aircraft.get_Bref()
    Sref = aircraft.get_Sref()
    
    AR = Bref**2 / Sref
    RT = ADN[1]
    Y1 = ADN[2]
    Y2 = ADN[3]
    reltrap = 100 * (Y2-Y1)/Y2
    
    # Check geometry conditions        
    if not check_geometry_2018(aircraft, l_engine):
        return 0
    
    # Compute Aerodynamics
    aircraft.CDp_Aircraft()
    try:
        aircraft.coefficients()
    except:
        return 0
    
    aircraft.CLmax(AR, RT, Bref, reltrap)
    
    # Compute performance
    aircraft.set_engine()
    mmax = aircraft.takeoff(altitude = 900, Sg_max = 50)
    dmin = aircraft.landing(altitude = 900, mass = mmax)
    
    # Compute Structures
    I_wing, mass_wing = aircraft.inertia_surface('wing')
    I_horizontal, mass_horizontal= aircraft.inertia_surface('horizontal')
    I_engine, mass_engine = aircraft.inertia_engine(xle = l_engine , yle = 0, zle = 0.05)
    I_main, mass_main = aircraft.inertia_main_gear()
    I_nose, mass_nose = aircraft.inertia_nose_gear(xle = l_engine , yle = 0, zle = 0)
    I_payload, mass_payload = aircraft.inertia_payload(mmax)
    
    Ix = float(I_wing[0] + I_horizontal[0] + I_engine[0] + I_main[0] + I_nose[0] + I_payload[0])
    Iy = float(I_wing[1] + I_horizontal[1] + I_engine[1] + I_main[1] + I_nose[1] + I_payload[1])
    Iz = float(I_wing[2] + I_horizontal[2] + I_engine[2] + I_main[2] + I_nose[2] + I_payload[2])
    
    
    # Compute Stability
    aircraft.static_stability()
    
    short_period, pughoid, dutch_roll, roll, spiral = aircraft.dynamic_stability(altitude, mmax, Ix, Iy, Iz )
    
    phase_B1 = aircraft.phase_B1(short_period)
    phase_C1 = aircraft.phase_C1(short_period)
    phase_B2 = aircraft.phase_B2(short_period)
    phase_C2 = aircraft.phase_C2(short_period)
    
    if (1<= phase_B1 <=3) and (1<= phase_C1 <=3) and (1<= phase_B2 <=3) and (1<= phase_C2 <=3):
        short_period_stability = True
    else:
        short_period_stability = False
    
    pughoid_level = aircraft.pughoid_level(pughoid)
    
    if (1<= pughoid_level <=2):
        pughoid_stability = True
    else:
        pughoid_stability = False
        
    dutch_roll_level_A = aircraft.dutch_roll_level_A(dutch_roll)
    dutch_roll_level_B = aircraft.dutch_roll_level_B(dutch_roll)
    dutch_roll_level_C = aircraft.dutch_roll_level_C(dutch_roll)
    
    if (1<= dutch_roll_level_A <=2) and (1<= dutch_roll_level_B <=2) and (1<= dutch_roll_level_C <=2):
        dutch_roll_stability = True
    else:
        dutch_roll_stability = False
        
    roll_level = aircraft.roll_level(roll)
    
    if (1<= roll_level <=2):
        roll_level_stability = True
    else:
        roll_level_stability = False
    
    spiral_level = aircraft.spiral_level(spiral)
    spiral_stability = (1<=spiral_level <= 3)
    
    stability = short_period_stability and pughoid_stability and dutch_roll_stability and roll_level_stability
    
    # Compute Score
    mass_tank = 0.064 + 0.143
    mass_electric_system = 0.25
    mass_epoxy = 0.1
    mass_monokote = 0.2
    mass_others = mass_tank + mass_electric_system + mass_epoxy + mass_monokote
    PVprev = mass_wing + mass_horizontal + mass_engine + \
        mass_main + mass_nose + mass_others
    PVreal = PVprev
    CPreal = mmax - PVreal
    FPV = max( 1.10 - 15*((PVprev - PVreal)/PVprev)**2, 0.8 )
    FPR = min(1.0, 0.5 + 0.9*NR/NRmax)
    Pcp = 12.5*CPreal
    Pvoo = FPV*FPR*Pcp 
    Bpo = mmax * (dmin <= Sgmax_landing)
    
    # Check design criteria
    Cma = aircraft.get_Cma() 
    design = (4<=AR)
    
    return (Pvoo + Bpo) * design * stability

def ADN_to_ref(ADN):
    l_engine = -0.12
    aircraft = ADN_to_Aircraft(ADN)
    Sref = aircraft.get_Sref()
    Bref = aircraft.get_Bref()
    Lref = aircraft.get_surface('horizontal').get_sections()[0].get_Xle() + \
           aircraft.get_surface('horizontal').get_sections()[0].get_Chord() + \
           aircraft.get_surface('wing').get_sections()[0].get_Xle() + \
           aircraft.get_surface('wing').get_sections()[0].get_Chord() + \
            max(-l_engine,0)
    Sd = Bref + Lref
    AR = Bref**2 / Sref
    print('Sref: ', Sref)
    print('Bref: ', Bref)
    print('AR:   ', AR)
    print('Sd:   ', Sd)
    

Cr_min = 0.2
Cr_max = 0.8
RT_min = 0.3
RT_max = 1.0
Y1_min = 0.2
Y1_max = 1.0
Y2_min = 0.2
Y2_max = 1.0
tx_min = 0.1
tx_max = 0.6
tz_min = 0.1
tz_max = 0.6
Ch_min = 0.15
Ch_max = 0.4
Yh_min = 0.2
Yh_max = 1.6
Zv_min = 0.1
Zv_max = 0.5

ga_instance = pygad.GA(num_generations=10,
                       num_parents_mating=20,
                       fitness_func=score,
                       sol_per_pop=20,
                       num_genes=9,
                       gene_space = [list(np.linspace(Cr_min, Cr_max, round((Cr_max-Cr_min)/0.02))),
                                     list(np.linspace(RT_min, RT_max, round((RT_max-RT_min)/0.1))),
                                     list(np.linspace(Y1_min, Y1_max, round((Y1_max-Y1_min)/0.02))),
                                     list(np.linspace(Y2_min, Y2_max, round((Y2_max-Y2_min)/0.02))),
                                     list(np.linspace(tx_min, tx_max, round((tx_max-tx_min)/0.02))),
                                     list(np.linspace(tz_min, tz_max, round((tz_max-tz_min)/0.02))),
                                     list(np.linspace(Ch_min, Ch_max, round((Ch_max-Ch_min)/0.02))),
                                     list(np.linspace(Yh_min, Yh_max, round((Yh_max-Yh_min)/0.02))),
                                     list(np.linspace(Zv_min, Zv_max, round((Zv_max-Zv_min)/0.02)))
                                     ]              
                       )

def run_10_instances(name):
    # Running Genetic Algorithm and Saving Data
    ga_instance.run()
    ga_instance.save(filename = os.getcwd() + '/Results/' + name + '_10')
    
    ga_instance.run()
    ga_instance.save(filename = os.getcwd() + '/Results/' + name + '_20')
    
    ga_instance.run()
    ga_instance.save(filename = os.getcwd() + '/Results/' + name + '_30')
    
    ga_instance.run()
    ga_instance.save(filename = os.getcwd() + '/Results/' + name + '_40')
    
    ga_instance.run()
    ga_instance.save(filename = os.getcwd() + '/Results/' + name + '_50')
    
    # ga_instance.run()
    ga_instance.save(filename = os.getcwd() + '/Results/' + name + '_60')
    
    ga_instance.run()
    ga_instance.save(filename = os.getcwd() + '/Results/' + name + '_70')
    
    ga_instance.run()
    ga_instance.save(filename = os.getcwd() + '/Results/' + name + '_80')
    
    ga_instance.run()
    ga_instance.save(filename = os.getcwd() + '/Results/' + name + '_90')
    
    ga_instance.run()
    ga_instance.save(filename = os.getcwd() + '/Results/' + name + '_100')
    
    solution, solution_fitness, solution_idx = ga_instance.best_solution()
    ga_instance.plot_fitness()

def saving_instances(filename, name):
    # Saving Data to CSV
    result_0 = ga_instance.initial_population
    result_1 = pygad.load(filename=os.getcwd() + '/Results/' + name + '_10').population
    result_2 = pygad.load(filename=os.getcwd() + '/Results/' + name + '_20').population
    result_3 = pygad.load(filename=os.getcwd() + '/Results/' + name + '_30').population
    result_4 = pygad.load(filename=os.getcwd() + '/Results/' + name + '_40').population
    result_6 = pygad.load(filename=os.getcwd() + '/Results/' + name + '_60').population
    result_8 = pygad.load(filename=os.getcwd() + '/Results/' + name + '_80').population
    result_F = pygad.load(filename=os.getcwd() + '/Results/' + name + '_100').population
    
    fitness_0 = []
    for genes in result_0:
        fitness_0.append( score(list(genes), 1))
    
    fitness_1 = []
    for genes in result_1:
        fitness_1.append( score(list(genes), 1))
        
    fitness_2 = []
    for genes in result_2:
        fitness_2.append( score(genes, 1))
        
    fitness_3 = []
    for genes in result_3:
        fitness_3.append( score(genes, 1))
        
    fitness_4 = []
    for genes in result_4:
        fitness_4.append( score(genes, 1))
        
    fitness_6 = []
    for genes in result_3:
        fitness_6.append( score(genes, 1))
        
    fitness_8 = []
    for genes in result_4:
        fitness_8.append( score(genes, 1))
        
    fitness_F = []
    for genes in result_F:
        fitness_F.append( score(genes, 1))
      
    df0 = pd.DataFrame(result_0)
    df1 = pd.DataFrame(result_1)
    df2 = pd.DataFrame(result_2)
    df3 = pd.DataFrame(result_3)
    df4 = pd.DataFrame(result_4)
    df6 = pd.DataFrame(result_6)
    df8 = pd.DataFrame(result_8)
    dfF = pd.DataFrame(result_F)
    
    r0 = pd.DataFrame(fitness_0)
    r1 = pd.DataFrame(fitness_1)
    r2 = pd.DataFrame(fitness_2)
    r3 = pd.DataFrame(fitness_3)
    r4 = pd.DataFrame(fitness_4)
    r6 = pd.DataFrame(fitness_6)
    r8 = pd.DataFrame(fitness_8)
    rF = pd.DataFrame(fitness_F)
    
    with pd.ExcelWriter(filename, engine='openpyxl', mode='w') as writer:
        df0.to_excel(writer, sheet_name='R0', header=False, index=False, startcol=0, startrow=0)
        df1.to_excel(writer, sheet_name='R1', header=False, index=False, startcol=0, startrow=0)
        df2.to_excel(writer, sheet_name='R2', header=False, index=False, startcol=0, startrow=0)
        df3.to_excel(writer, sheet_name='R3', header=False, index=False, startcol=0, startrow=0)
        df4.to_excel(writer, sheet_name='R4', header=False, index=False, startcol=0, startrow=0)
        df6.to_excel(writer, sheet_name='R6', header=False, index=False, startcol=0, startrow=0)
        df8.to_excel(writer, sheet_name='R8', header=False, index=False, startcol=0, startrow=0)
        dfF.to_excel(writer, sheet_name='RF', header=False, index=False, startcol=0, startrow=0)
        
        r0.to_excel(writer, sheet_name='R0', header=False, index=False, startcol=12, startrow=0)
        r1.to_excel(writer, sheet_name='R1', header=False, index=False, startcol=12, startrow=0)
        r2.to_excel(writer, sheet_name='R2', header=False, index=False, startcol=12, startrow=0)
        r3.to_excel(writer, sheet_name='R3', header=False, index=False, startcol=12, startrow=0)
        r4.to_excel(writer, sheet_name='R4', header=False, index=False, startcol=12, startrow=0)
        r6.to_excel(writer, sheet_name='R6', header=False, index=False, startcol=12, startrow=0)
        r8.to_excel(writer, sheet_name='R8', header=False, index=False, startcol=12, startrow=0)
        rF.to_excel(writer, sheet_name='RF', header=False, index=False, startcol=12, startrow=0)


# filename = 'Resultados.xlsx'
# name = "GA"

# run_10_instances(name)

# # Record the end time
# end_time = time.time()

# # Calculate the elapsed time
# elapsed_time = end_time - start_time

# saving_instances(filename, name)

# ga_instance.plot_fitness()

# print(f"Elapsed time: {elapsed_time:.6f} seconds")