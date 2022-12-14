# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This file contains the different functions used to compute parameters or variables in the simulation.
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

import numpy as np
import math
import random
from scipy.ndimage import binary_dilation

#-------------------------------------------------------------------------------

def Compute_S_int(dict_sample):
    '''
    Searching Surface,
    Monte Carlo Method
    A box is defined, we take a random point and we look if it is inside or outside the grain
    Properties are the statistic times the box properties


        Input :
            a sample dictionnary (a dict)
        Output :
            Nothing but the dictionnary gets an updated value for the intersection surface (a float)
    '''
    #Defining the limit
    box_min_x = min(dict_sample['L_g'][1].l_border_x)
    box_max_x = max(dict_sample['L_g'][0].l_border_x)
    box_min_y = min(dict_sample['L_g'][0].l_border_y)
    box_max_y = max(dict_sample['L_g'][0].l_border_y)

    #Compute the intersection surface
    N_MonteCarlo = 5000 #The larger it is, the more accurate it is
    sigma = 1
    M_Mass = 0

    for i in range(N_MonteCarlo):
        P = np.array([random.uniform(box_min_x,box_max_x),random.uniform(box_min_y,box_max_y)])
        if dict_sample['L_g'][0].P_is_inside(P) and dict_sample['L_g'][1].P_is_inside(P):
            M_Mass = M_Mass + sigma

    Mass = (box_max_x-box_min_x)*(box_max_y-box_min_y)/N_MonteCarlo*M_Mass
    Surface = Mass/sigma

    #Update element in dictionnary
    dict_sample['S_int'] = Surface

#-------------------------------------------------------------------------------

def Compute_sum_min_etai(dict_sample, dict_sollicitation):
    '''
    Compute the sum over the sample of the minimum of etai.

    This variable is used for Compute_Emec().

        Input :
            a sample dictionnary (a dict)
        Output :
            Nothing but the dictionnary gets an updated value (a float)
            If it is the first call of the function, dictionnaries gets new value (2 floats)
    '''
    #Initialisation
    sum_min_etai = 0
    #compute the sum over the sample of the minimum of etai
    for l in range(len(dict_sample['y_L'])):
        for c in range(len(dict_sample['x_L'])):
            sum_min_etai = sum_min_etai + min(dict_sample['L_g'][0].etai_M[-1-l][c],dict_sample['L_g'][1].etai_M[-1-l][c])

    #Update element in dictionnary
    dict_sample['sum_min_etai'] = sum_min_etai
    #create element in dictionnary if not already created
    if 'sum_min_etai0' not in dict_sample.keys():
        dict_sample['sum_min_etai0'] = sum_min_etai
        dict_sollicitation['alpha'] = 0.2*sum_min_etai

#-------------------------------------------------------------------------------

def Compute_Emec(dict_sample, dict_sollicitation):
    '''
    Compute the mechanical energy field in the sample.

        Input :
            a sample dictionnary (a dict)
        Output :
            Nothing but the dictionnary gets an updated value for the mechanical energy map (a nx x ny numpy array)
    '''
    #Initialisation
    Emec_M = np.array(np.zeros((len(dict_sample['y_L']),len(dict_sample['x_L']))))
    #compute the variable e_mec
    e_mec = dict_sollicitation['alpha']/dict_sample['sum_min_etai']
    #compute the distribution of the mechanical energy
    for l in range(len(dict_sample['y_L'])):
        for c in range(len(dict_sample['x_L'])):
            Emec_M[-1-l][c] = e_mec*min(dict_sample['L_g'][0].etai_M[-1-l][c],dict_sample['L_g'][1].etai_M[-1-l][c])

    #Update element in dictionnary
    dict_sample['Emec_M'] = Emec_M

#-------------------------------------------------------------------------------

def Compute_kc_dil(dict_material, dict_sample):
    '''
    Compute the solute diffusion coefficient field in the sample.

    Here, a dilation method is applied. For all node, a Boolean variable is defined.
    This variable is True if eta_i and eta_j are greater than 0.5 (in the contact zone).
                  is True if eta_i and eta_j are lower than 0.5 (in the pore zone).
                  is False else.

    A dilation method is applied, the size of the structural element is the main case.

    The diffusion map is built on the Boolean map. If the variable is True, the diffusion is kc, else 0.

        Input :
            a material dictionnary (a dict)
            a sample dictionnary (a dict)
        Output :
            Nothing but the dictionnary gets an updated value for the solute diffusion coefficient map (a nx x ny numpy array)
    '''
    #Initialisation
    on_off_M = np.array(np.zeros((len(dict_sample['y_L']),len(dict_sample['x_L']))), dtype = bool)

    #compute the on off map
    for l in range(len(dict_sample['y_L'])):
        for c in range(len(dict_sample['x_L'])):
            #at the contact
            if dict_sample['L_g'][0].etai_M[-1-l][c] > 0.5 and dict_sample['L_g'][1].etai_M[-1-l][c] > 0.5:
                on_off_M[-l-1][c] = True
            #in the pore space
            elif dict_sample['L_g'][0].etai_M[-1-l][c] < 0.5 and dict_sample['L_g'][1].etai_M[-1-l][c] < 0.5:
                on_off_M[-l-1][c] = True

    #dilatation
    struct_element = np.array(np.ones((8,8)), dtype = bool)
    dilated_M = binary_dilation(on_off_M, struct_element)

    #compute the map of the solute diffusion coefficient
    kc_M = np.array(np.zeros((len(dict_sample['y_L']),len(dict_sample['x_L']))))
    for l in range(len(dict_sample['y_L'])):
        for c in range(len(dict_sample['x_L'])):
            if dilated_M[-1-l][c] :
                kc_M[-1-l][c] = dict_material['kappa_c']

    #Update element in dictionnary
    dict_sample['kc_M'] = kc_M

#-------------------------------------------------------------------------------

def Compute_kc_wfd(dict_material, dict_sample):
    '''
    Compute the solute diffusion coefficient field in the sample.

    Here, a water film diffusion is assumed. For all node, a Boolean variable is defined.
    This variable is True if eta_i and eta_j are greater than 0.5 (in the contact zone).
                  is True if eta_i and eta_j are lower than 0.5 (in the pore zone).
                  is False else.

    Then, the nodes inside the water film diffusion are searched. This film is assumed centered on x = 0 (approximately the contact center).
    All nodes in this film get a variable True.

    The diffusion map is built on the Boolean map. If the variable is True, the diffusion is kc, else 0.

        Input :
            a material dictionnary (a dict)
            a sample dictionnary (a dict)
        Output :
            Nothing but the dictionnary gets an updated value for the solute diffusion coefficient map (a nx x ny numpy array)
    '''
    #Initialisation
    on_off_M = np.array(np.zeros((len(dict_sample['y_L']),len(dict_sample['x_L']))), dtype = bool)

    #compute the on off map
    for l in range(len(dict_sample['y_L'])):
        for c in range(len(dict_sample['x_L'])):
            #at the contact
            if dict_sample['L_g'][0].etai_M[-1-l][c] > 0.5 and dict_sample['L_g'][1].etai_M[-1-l][c] > 0.5:
                on_off_M[-l-1][c] = True
            #in the pore space
            elif dict_sample['L_g'][0].etai_M[-1-l][c] < 0.5 and dict_sample['L_g'][1].etai_M[-1-l][c] < 0.5:
                on_off_M[-l-1][c] = True

    #look for nodes delimiting the wfd
    c = 0
    search_limits_wfd = True
    search_first = True
    while search_limits_wfd :
        if search_first and dict_sample['x_L'][c] > -dict_material['width_wfd']/2:
            c_start = c - 1
            search_first = False
        elif dict_sample['x_L'][c] > dict_material['width_wfd']/2:
            c_end = c
            search_limits_wfd = False
        c = c + 1

    #force diffusion in the wter film diffusion
    for c in range(c_start,c_end+1):
        for l in range(len(dict_sample['y_L'])):
            on_off_M[-l-1][c] = True

    #compute the map of the solute diffusion coefficient
    kc_M = np.array(np.zeros((len(dict_sample['y_L']),len(dict_sample['x_L']))))
    for l in range(len(dict_sample['y_L'])):
        for c in range(len(dict_sample['x_L'])):
            if dilated_M[-1-l][c] :
                kc_M[-1-l][c] = dict_material['kappa_c']

    #Update element in dictionnary
    dict_sample['kc_M'] = kc_M

#-------------------------------------------------------------------------------

def Compute_kc_int(dict_material, dict_sample):
    '''
    Compute the solute diffusion coefficient field in the sample.

    For all nodes in the mesh a diffusion coefficient is computed.
    If eta_i and eta_j are greater than 0.5 (in the contact zone), the diffusion is kc0.
    If eta_i and eta_j are lower than 0.5 (in the pore zone), the diffusion is kc0.
    If eta_i (resp. eta_j) is greater than 0.5 and eta_j (resp. eta_i) is lower than 0.5 (in one grain but not the other), an interpolated diffusion is used.

    This interpolated diffusion is kc = kc0*exp(tau*(d_to_center_i-radius_i)/radius_i).

        Input :
            a material dictionnary (a dict)
            a sample dictionnary (a dict)
        Output :
            Nothing but the dictionnary gets an updated value for the solute diffusion coefficient map (a nx x ny numpy array)
    '''
    #Initialisation
    kc_M = np.array(np.zeros((len(dict_sample['y_L']),len(dict_sample['x_L']))))

    #compute the distribution of the solute diffusion coefficient
    for l in range(len(dict_sample['y_L'])):
        for c in range(len(dict_sample['x_L'])):
            #at the contact
            if dict_sample['L_g'][0].etai_M[-1-l][c] > 0.5 and dict_sample['L_g'][1].etai_M[-1-l][c] > 0.5:
                kc_M[-l-1][c] = dict_material['kappa_c']
            #inside g1 and not g2
            elif dict_sample['L_g'][0].etai_M[-1-l][c] > 0.5 and dict_sample['L_g'][1].etai_M[-1-l][c] < 0.5:
                P = np.array([dict_sample['x_L'][c], dict_sample['y_L'][-1-l]])
                Distance = np.linalg.norm(P - dict_sample['L_g'][0].center)
                #exponential decrease
                kappa_c_trans = dict_material['kappa_c']*math.exp(-(dict_sample['L_g'][0].r_mean-Distance)/(dict_sample['L_g'][0].r_mean/dict_material['tau_kappa_c']))
                kc_M[-l-1][c] = kappa_c_trans
            #inside g2 and not g1
            elif dict_sample['L_g'][0].etai_M[-1-l][c] < 0.5 and dict_sample['L_g'][1].etai_M[-1-l][c] > 0.5:
                #compute the distance to g1
                P = np.array([dict_sample['x_L'][c], dict_sample['y_L'][-1-l]])
                Distance = np.linalg.norm(P - dict_sample['L_g'][1].center)
                #exponential decrease
                kappa_c_trans = dict_material['kappa_c']*math.exp(-(dict_sample['L_g'][1].r_mean-Distance)/(dict_sample['L_g'][1].r_mean/dict_material['tau_kappa_c']))
                kc_M[-l-1][c] = kappa_c_trans
            #outside
            else :
                kc_M[-l-1][c] = dict_material['kappa_c']

    #Update element in dictionnary
    dict_sample['kc_M'] = kc_M

#-------------------------------------------------------------------------------

def Compute_sum_c(dict_sample):
    '''
    Compute the quantity of the solute.

        Input :
            a sample dictionnary (a dict)
        Output :
            Nothing but the dictionnary gets an updated value for the intersection surface (a float)
    '''
    sum_c = 0
    for l in range(len(dict_sample['y_L'])):
        for c in range(len(dict_sample['x_L'])):
            sum_c = sum_c + dict_sample['solute_M'][l][c]

    #update element in dict
    dict_sample['sum_c'] = sum_c

#-------------------------------------------------------------------------------

def Compute_sum_eta(dict_sample):
    '''
    Compute the quantity of the grain.

        Input :
            a sample dictionnary (a dict)
        Output :
            Nothing but the dictionnary gets an updated value for the intersection surface (a float)
    '''
    sum_eta = 0
    for grain in dict_sample['L_g']:
        for l in range(len(dict_sample['y_L'])):
            for c in range(len(dict_sample['x_L'])):
                sum_eta = sum_eta + grain.etai_M[l][c]

    #update element in dict
    dict_sample['sum_eta'] = sum_eta

#-------------------------------------------------------------------------------

def Compute_sum_Ed_plus_minus(dict_sample, dict_sollicitation):
    '''
    Compute the total energy in the sample. This energy is divided in a plus term and a minus term.

        Input :
            a sample dictionnary (a dict)
        Output :
            Nothing but the dictionnary gets an updated value for energy inside the sample (three floats)
    '''
    sum_ed_plus = 0
    sum_ed_minus = 0
    sum_ed = 0
    sum_Ed_mec = 0
    sum_Ed_che = 0
    for l in range(len(dict_sample['y_L'])):
        for c in range(len(dict_sample['x_L'])):

            #Emec
            Ed_mec = dict_sample['Emec_M'][-1-l][c]

            #Eche
            Ed_che = dict_sollicitation['chi']*dict_sample['solute_M'][-1-l][c]*(3*dict_sample['L_g'][0].etai_M[-1-l][c]**2-2*dict_sample['L_g'][0].etai_M[-1-l][c]**3+\
                                                                                 3*dict_sample['L_g'][1].etai_M[-1-l][c]**2-2*dict_sample['L_g'][1].etai_M[-1-l][c]**3)

            #Ed
            Ed = Ed_mec - Ed_che

            #sum actualisation
            sum_ed = sum_ed + Ed
            sum_Ed_mec = sum_Ed_mec + Ed_mec
            sum_Ed_che = sum_Ed_che + Ed_che
            if Ed > 0 :
                sum_ed_plus = sum_ed_plus + Ed
            else :
                sum_ed_minus = sum_ed_minus - Ed

    #update elements in dict
    dict_sample['sum_ed'] = sum_ed
    dict_sample['sum_Ed_mec'] = sum_Ed_mec
    dict_sample['sum_Ed_che'] = sum_Ed_che
    dict_sample['sum_ed_plus'] = sum_ed_plus
    dict_sample['sum_ed_minus'] = sum_ed_minus
