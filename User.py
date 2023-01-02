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

    N_grain = 100 #total number of grains

    R0 = 350 #µm radius to compute the grain distribution
    L_R = [1.2*R0,1.1*R0,0.9*R0,0.8*R0] #from larger to smaller
    L_percentage_R = [1/6,1/3,1/3,1/6] #distribution of the different radius
    #Recompute the mean radius
    R_mean = 0
    for i in range(len(L_R)):
        R_mean = R_mean + L_R[i]*L_percentage_R[i]

    #write dict
    dict_geometry = {
    'R_mean' : R_mean,
    'L_R' : L_R,
    'L_percentage_R' : L_percentage_R,
    'N_grain' : N_grain
    }

    #---------------------------------------------------------------------------
    #Sample parameters

    #Box définition
    x_box_min = 0 #µm
    x_box_max = 2*R_mean*math.sqrt(N_grain/0.6) #µm
    y_box_min = 0 #µm

    #spatial discretisation
    nx = 180
    ny = 100

    #approximatively the number of vertices for one grain during DEM simulation
    grain_discretisation = 20

    dict_sample = {
    'nx' : nx,
    'ny' : ny,
    'grain_discretisation' : grain_discretisation,
    'x_box_min' : x_box_min,
    'x_box_max' : x_box_max,
    'y_box_min' : y_box_min
    }

    #---------------------------------------------------------------------------
    #Material parameters

    #phase field
    Mobility = 1 #L1, L2 in .i
    #Diffusion of etai
    kappa_eta = 0.01
    #Diffusion of c, the solute generated by the dissolution
    kappa_c = 50
    method_to_compute_kc = 'dilation' #dilation, wfd or interpolation
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
    'M' : Mobility,
    'kappa_eta' : kappa_eta,
    'kappa_c' : kappa_c,
    'tau_kappa_c' : tau_kappa_c,
    'method_to_compute_kc' : method_to_compute_kc,
    'Y' : Y,
    'nu' : nu,
    'rho' : rho,
    'rho_surf' : rho_surf,
    'mu_friction_gg' : mu_friction_gg,
    'mu_friction_gw' : mu_friction_gw,
    'coeff_restitution' : coeff_restitution
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
    dt_DEM_crit = math.pi*min(L_R)/(0.16*nu+0.88)*math.sqrt(rho*(2+2*nu)/Y) #s critical time step from O'Sullivan 2011
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
    L_flag_plot = ['Config', 'Sphericity', 'Movie', 'sum_Ed']
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
    gravity = 0
    Vertical_Confinement_Linear_Force = Y*2*R_mean/1000 #µN/µm used to compute the Vertical_Confinement_Force
    Vertical_Confinement_Force = Vertical_Confinement_Linear_Force*(x_box_max-x_box_min) #µN

    dict_sollicitation = {
    'overlap_target' : overlap_target,
    'chi' : chi,
    'gravity' : gravity,
    'Vertical_Confinement_Force' : Vertical_Confinement_Force
    }

    #---------------------------------------------------------------------------
    #Initial condition parameters

    n_generation = 2 #number of grains generation
    factor_ymax_box = 2.5 #margin to generate grains
    N_test_max = 5000 # maximum number of tries to generate a grain without overlap
    i_DEM_stop_IC = 2000 #stop criteria for DEM during IC
    Debug_DEM_IC = True #plot configuration inside DEM during IC
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

    return dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitation

#-------------------------------------------------------------------------------

def Add_mesh(dict_geometry, dict_sample):
    '''
    Generate the mesh in the sample.

        Input :
            a geometry dictionnary (a dict)
            a sample dictionnary (a dictionnary)
        Output :
            Nothing but the sample dictionnary gets a new grains configuration (a list)
    '''
    #spatial discretisation
    x_min = dict_sample['x_box_min'] - 0.1*dict_geometry['R_mean']
    x_max = dict_sample['x_box_max'] + 0.1*dict_geometry['R_mean']
    x_L = np.linspace(x_min,x_max,dict_sample['nx'])
    y_min = dict_sample['x_box_min'] - 0.1*dict_geometry['R_mean']
    y_max = dict_sample['x_box_max'] + 0.1*dict_geometry['R_mean']
    y_L = np.linspace(y_min,y_max,dict_sample['nx'])

    #add element in dict
    dict_sample['x_L'] = x_L
    dict_sample['y_L'] = y_L

#-------------------------------------------------------------------------------

def Add_variables_needed(dict_material, dict_sample):
    '''
    Generate the mesh in the sample.

        Input :
            a material dictionnary (a dict)
            a sample dictionnary (a dictionnary)
        Output :
            Nothing but the sample dictionnary gets new variables
    '''
    #Arrays with only 0
    dict_sample['Emec_M'] = np.array(np.zeros((len(dict_sample['y_L']),len(dict_sample['x_L']))))
    dict_sample['Eche_M'] = np.array(np.zeros((len(dict_sample['y_L']),len(dict_sample['x_L'])))),
    dict_sample['Ed_M'] = np.array(np.zeros((len(dict_sample['y_L']),len(dict_sample['x_L'])))),

    #compute the phase field width interface
    dict_material['w'] = math.sqrt((dict_sample['x_L'][4]-dict_sample['x_L'][0])**2+(dict_sample['y_L'][4]-dict_sample['y_L'][0])**2)

    #Energy barrier of etai
    dict_material['Energy_barrier'] = 20*dict_material['kappa_eta']/(dict_material['w'])**2

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
