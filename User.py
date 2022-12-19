# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This is the file where the user can change the different parameters for the simulation.
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

import math
import numpy as np
from pathlib import Path

#Own functions and classes
import Grain

#-------------------------------------------------------------------------------
#User
#-------------------------------------------------------------------------------

def All_parameters():
    '''
    This function is called in main.py to have all the parameters needed in the simulation.

        Input :
            Nothing
        Output :
            an algorithm dictionnary (a dictionnary)
            a material dictionnary (a dictionnary)
            a sample dictionnary (a dictionnary)
            a sollicitation dictionnary (a dictionnary)
    '''
    #---------------------------------------------------------------------------
    #Geometry parameters

    N_grain = 300 #total number of grains

    R0 = 350 #µm radius to compute the grain distribution
    L_R = [1.2*R_mean,1.1*R_mean,0.9*R_mean,0.8*R_mean] #from larger to smaller
    L_percentage_R = [1/6,1/3,1/3,1/6] #distribution of the different radius
    #Recompute the mean radius
    R_mean = 0
    for i in range(len(L_R)):
        R_mean = R_mean + L_R[i]*L_percentage_R[i]

    #write dict
    dict_geometry = {
    'N_grain' : N_grain,
    'R_mean' : R_mean,
    'L_R' : L_R,
    'L_percentage_R' : L_percentage_R
    }

    #---------------------------------------------------------------------------
    #Sample parameters

    #Box définition
    x_box_min = -230 #µm
    x_box_max = 230 #µm
    y_box_min = -130 #µm

    #spatial discretisation
    x_min = -230
    x_max = 230
    nx = 180
    x_L = np.linspace(x_min,x_max,nx)
    y_min = -130
    y_max = 130
    ny = 100
    y_L = np.linspace(y_min,y_max,ny)

    #approximatively the number of vertices for one grain during DEM simulation
    grain_discretisation = 20

    dict_sample = {
    'x_L' : x_L,
    'y_L' : y_L,
    'grain_discretisation' : grain_discretisation,
    'Emec_M' : np.array(np.zeros((len(y_L),len(x_L)))),
    'Eche_M' : np.array(np.zeros((len(y_L),len(x_L)))),
    'Ed_M' : np.array(np.zeros((len(y_L),len(x_L)))),
    'x_box_min' : x_box_min,
    'x_box_max' : x_box_max,
    'y_box_min' : y_box_min
    }

    #---------------------------------------------------------------------------
    #Algorithm parameters

    np_proc = 4 #number of processor used
    n_t_PFDEM = 200 #number of cycle PF-DEM

    #Time step for phase field
    n_t_PF = 8
    dt_PF_init = 0.3
    dt_PF_level1 = dt_PF_init/2
    dt_PF_level2 = dt_PF_level1/2
    dt_PF_level3 = dt_PF_level2/2
    #criteria to switch level
    Ed_level1 = 10
    Ed_level2 = 20
    Ed_level3 = 25

    #DEM parameters
    dt_DEM_crit = math.pi*min(min(L_Dimension)/2,min(L_R))/(0.16*nu+0.88)*math.sqrt(rho*(2+2*nu)/Y) #s critical time step from O'Sullivan 2011
    dt_DEM = dt_DEM_crit/8 #s time step during DEM simulation
    factor_neighborhood = 1.5 #margin to detect a grain into a neighborhood
    i_update_neighborhoods = 100 #the frequency of the update of the neighborhood of the grains and the walls
    #Stop criteria of the DEM
    i_DEM_stop = 3000 #maximum iteration for one DEM simulation
    Ecin_ratio = 0.0002
    n_window_stop = 50
    dy_box_max_stop = 0.5

    #Margin for sphericity study
    sphericity_margin = 0.05
    #Discretisation to find the inscribing (number of nodes in one direction)
    n_spatial_inscribing = 100

    #List of plot to do
    Debug = True #plot configuration before and after DEM simulation
    Debug_DEM = False #plot configuration inside DEM
    i_print_plot = 200 #frenquency of the print and plot (if Debug_DEM) in DEM step
    # Config, C_at_P, Diff_Solute, dt, Ed, Eta_c, Init_Current_Shape, Kc, Movie (need Config to work), Sint_MinEtai, Sphericity, sum_Ed
    L_flag_plot = ['Config', 'C_at_P', 'Diff_Solute', 'dt', 'Sphericity', 'Kc', 'Movie', 'sum_Ed']
    #Visual parameters (for plot Config)
    c_min = 0
    c_max = 0.08

    #Save the simulation
    SaveData = True #Save data or not
    clean_memory = True #delete Data, Input, Output at the end of the simulation
    foldername = 'Data_2G_ACS' #name of the folder where data are saved
    template = 'PS_Long_Run' #template of the name of the simulation
    if SaveData :
        i_run = 1
        folderpath = Path('../'+foldername+'/'+template+'_'+str(i_run))
        while folderpath.exists():
            i_run = i_run + 1
            folderpath = Path('../'+foldername+'/'+template+'_'+str(i_run))
        namefile = template+'_'+str(i_run)
    else :
        namefile = template

    dict_algorithm = {
    'c_min' : c_min,
    'c_max' : c_max,
    'sphericity_margin' : sphericity_margin,
    'n_spatial_inscribing' : n_spatial_inscribing,
    'np_proc' : np_proc,
    'Debug' : Debug,
    'Debug_DEM' : Debug_DEM,
    'i_print_plot' : i_print_plot,
    'factor_neighborhood' : factor_neighborhood,
    'factor_distribution_etai' : factor_distribution_etai,
    'clean_memory' : clean_memory,
    'SaveData' : SaveData,
    'namefile' : namefile,
    'dt_PF' : dt_PF_init,
    'dt_PF_init' : dt_PF_init,
    'dt_PF_level1' : dt_PF_level1,
    'dt_PF_level2' : dt_PF_level2,
    'dt_PF_level3' : dt_PF_level3,
    'Ed_level1' : Ed_level1,
    'Ed_level2' : Ed_level2,
    'Ed_level3' : Ed_level3,
    'n_t_PF' : n_t_PF,
    'dt_DEM' : dt_DEM,
    'i_update_neighborhoods': i_update_neighborhoods,
    'i_DEM_stop' : i_DEM_stop,
    'Ecin_ratio' : Ecin_ratio,
    'n_window_stop' : n_window_stop,
    'dy_box_max_stop' : dy_box_max_stop,
    'foldername' : foldername,
    'n_t_PFDEM' : n_t_PFDEM,
    'L_flag_plot' : L_flag_plot
    }

    #---------------------------------------------------------------------------
    #External sollicitation parameters

    overlap_target = 10 #overlap verified before each phase field iteration
    chi = 0.5 #coefficient applied to the chemical energy

    dict_sollicitation = {
    'overlap_target' : overlap_target,
    'chi' : chi
    }

    #---------------------------------------------------------------------------
    #Material parameters

    #phase field
    width_int = math.sqrt((x_L[4]-x_L[0])**2+(y_L[4]-y_L[0])**2)
    Mobility = 1 #L1, L2 in .i
    #Diffusion of etai
    kappa_eta = 0.01
    #Energy barrier of etai
    Energy_barrier = 20*kappa_eta/(width_int)**2
    #Diffusion of c, the solute generated by the dissolution
    kappa_c = 50
    method_to_compute_kc = 'wfd' #dilation, wfd or interpolation
    #Definition of the evolution of the penalty term inside a grain (for interpolation)
    tau_kappa_c = 5 #kc = kc0 e(-(R_i-d)/R_i/tau_kappa_c)
    #DEM
    Y = 70*(10**9)*(10**6)*(10**(-12)) #Young Modulus µN/µm2
    nu = 0.3 #Poisson's ratio
    rho = 2500*10**(-6*3) #density kg/µm3
    rho_surf = 4/3*rho*R_mean #kg/µm2
    mu_friction_gg = 0.5 #grain-grain
    mu_friction_gw = 0 #grain-wall
    coeff_restitution = 0.2 #1 is perfect elastic

    dict_material = {
    'w' : width_int,
    'M' : Mobility,
    'kappa_eta' : kappa_eta,
    'kappa_c' : kappa_c,
    'tau_kappa_c' : tau_kappa_c,
    'method_to_compute_kc' : method_to_compute_kc,
    'Energy_barrier' : Energy_barrier,
    'Y' : Y,
    'nu' : nu,
    'rho' : rho,
    'rho_surf' : rho_surf,
    'mu_friction_gg' : mu_friction_gg,
    'mu_friction_gw' : mu_friction_gw,
    'coeff_restitution' : coeff_restitution
    }

    #---------------------------------------------------------------------------
    #Initial condition parameters

    n_generation = 2 #number of grains generation
    factor_ymax_box = 2.5 #margin to generate grains
    N_test_max = 5000 # maximum number of tries to generate a grain without overlap
    i_DEM_stop_IC = 2000 #stop criteria for DEM during IC
    Debug_DEM_IC = False #plot configuration inside DEM during IC
    i_print_plot_IC = 200 #frenquency of the print and plot (if Debug_DEM_IC) for IC
    dt_DEM_IC = dt_DEM_crit/5 #s time step during IC
    Ecin_ratio_IC = 0.0005
    factor_neighborhood_IC = 1.5 #margin to detect a grain into a neighborhood
    i_update_neighborhoods_gen = 5 #the frequency of the update of the neighborhood of the grains and the walls during IC generations
    i_update_neighborhoods_com = 100 #the frequency of the update of the neighborhood of the grains and the walls during IC combination

    #write dict
    dict_ic = {
    'n_generation' : n_generation,
    'i_update_neighborhoods_gen': i_update_neighborhoods_gen,
    'i_update_neighborhoods_com': i_update_neighborhoods_com,
    'factor_ymax_box' : factor_ymax_box,
    'i_DEM_stop_IC' : i_DEM_stop_IC,
    'Debug_DEM' : Debug_DEM_IC,
    'dt_DEM_IC' : dt_DEM_IC,
    'Ecin_ratio_IC' : Ecin_ratio_IC,
    'i_print_plot_IC' : i_print_plot_IC,
    'factor_neighborhood_IC' : factor_neighborhood_IC,
    'N_test_max' : N_test_max
    }


    #---------------------------------------------------------------------------

    return dict_algorithm, dict_material, dict_sample, dict_sollicitation

#-------------------------------------------------------------------------------

def Add_2grains(dict_material,dict_sample):
    '''
    Generate two grains in the sample.

        Input :
            a material dictionnary (a dictionnary)
            a sample dictionnary (a dictionnary)
        Output :
            Nothing but the sample dictionnary gets a new grains configuration (a list)
    '''
    #grain 1
    radius = 100
    center = np.array([np.mean(dict_sample['x_L'])-radius,np.mean(dict_sample['y_L'])])
    grain_1 = Grain.Grain(1,radius,center,dict_material,dict_sample)

    #grain 2
    radius = 100
    center = np.array([np.mean(dict_sample['x_L'])+radius,np.mean(dict_sample['y_L'])])
    grain_2 = Grain.Grain(2,radius,center,dict_material,dict_sample)

    #add element in dict
    dict_sample['L_g'] = [grain_1, grain_2]

#-------------------------------------------------------------------------------

def Add_S0(dict_sample, dict_sollicitation):
    '''
    Compute the initial intersection surface between two disk particles.
    See https://calculis.net/q/aire-intersection-disques-35

        Input :
            a sample dictionnary (a dictionnary)
        Output :
            Nothing but the sample dictionnary gets a new entry (a float)
    '''
    #notation from the forum
    R = dict_sample['L_g'][0].r_mean
    R_prime = dict_sample['L_g'][1].r_mean
    D = R + R_prime - dict_sollicitation['overlap_target']
    d  = (R**2 + D**2 - R_prime**2)/2/D
    d_prime = D - d
    #Part from particle 1
    S1 = R**2*math.acos(d/R)-d*math.sqrt(R**2-d**2)
    #Part from particle 2
    S2 = R_prime**2*math.acos(d_prime/R_prime)-d_prime*math.sqrt(R_prime**2-d_prime**2)

    #add element in dictionnary
    dict_sample['S_int_0'] = S1 + S2

#-------------------------------------------------------------------------------

def Add_solute(dict_sample):
    '''
    Generate solute in the sample.

        Input :
            a sample dictionnary (a dictionnary)
        Output :
            Nothing but the sample dictionnary gets a new solute configuration (a n_y x n_x numpy array)
    '''
    #solute
    solute_M = np.zeros((len(dict_sample['y_L']),len(dict_sample['x_L'])))

    #add element in dict
    dict_sample['solute_M'] = solute_M

#-------------------------------------------------------------------------------

def Add_alpha_emec(dict_sample, dict_sollicitation):
    '''
    Compute the coefficient alpha applied to the mechanical energy.

        Input :
            a sample dictionnary (a dictionnary)
            a sollicitation dictionnary (a dictionnary)
        Output :
            Nothing but the sollicitation dictionnary gets a new attribut (a float)
    '''
    alpha =  0.2*dict_sample['S_int_0']

    #update dict
    dict_sollicitation['alpha'] = alpha

#-------------------------------------------------------------------------------

def Criteria_StopSimulation(dict_algorithm):
    '''
    Define a stop criteria for the PFDEM simulation.

    The simulation stops if the number of iterations reaches a maximum value.

        Input :
            an algorithm dictionnary (a dictionnary)
        Output :
            The result depends on the fact if the stop criteria is reached or not (a bool)
    '''
    Criteria_Verified = False
    if dict_algorithm['i_PFDEM'] >= dict_algorithm['n_t_PFDEM']:
        Criteria_Verified = True
    return Criteria_Verified
