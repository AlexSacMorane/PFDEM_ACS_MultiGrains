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
import matplotlib.pyplot as plt #to delete

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

def Compute_Emec(dict_material, dict_sample, dict_sollicitation):
    '''
    Compute the mechanical energy over the sample.

    There is an iteration over all the contacts detected (grain-grain and grain-wall). First, the sum of the minimal grain phase variables is computed.
    Then, the mechanical energy is computed.

        Input :
            a material dictionnary (a dict)
            a sample dictionnary (a dict)
            a sollicitation dictionnary (a dict)
        Output :
            Nothing, but the sample dictionnary gets an updated value for the mechanical term (a nx x ny numpy array)
    '''
    #Initialisation
    Emec_M = np.array(np.zeros((len(dict_sample['y_L']),len(dict_sample['x_L']))))
    #contact grain-grain part
    for contact in dict_sample['L_contact']:
        #extract a spatial zone
        x_min = min(min(contact.g1.l_border_x),min(contact.g2.l_border_x))-dict_material['w']
        x_max = max(max(contact.g1.l_border_x),max(contact.g2.l_border_x))+dict_material['w']
        y_min = min(min(contact.g1.l_border_y),min(contact.g2.l_border_y))-dict_material['w']
        y_max = max(max(contact.g1.l_border_y),max(contact.g2.l_border_y))+dict_material['w']
        #look for this part inside the global mesh
        #create search list
        x_L_search_min = abs(np.array(dict_sample['x_L'])-x_min)
        x_L_search_max = abs(np.array(dict_sample['x_L'])-x_max)
        y_L_search_min = abs(np.array(dict_sample['y_L'])-y_min)
        y_L_search_max = abs(np.array(dict_sample['y_L'])-y_max)
        #get index
        i_x_min = list(x_L_search_min).index(min(x_L_search_min))
        i_x_max = list(x_L_search_max).index(min(x_L_search_max))
        i_y_min = list(y_L_search_min).index(min(y_L_search_min))
        i_y_max = list(y_L_search_max).index(min(y_L_search_max))
        #Initialisation
        sum_min_etai = 0
        #compute the sum over the sample of the minimum of etai
        for l in range(i_y_min, i_y_max):
            for c in range(i_x_min, i_x_max):
                sum_min_etai = sum_min_etai + min(contact.g1.etai_M[-1-l][c],contact.g2.etai_M[-1-l][c])
        if sum_min_etai != 0 :
            #compute the variable e_mec
            e_mec = dict_sollicitation['alpha']/sum_min_etai
        else :
            e_mec = 0
        #compute the distribution of the mechanical energy
        for l in range(i_y_min, i_y_max):
            for c in range(i_x_min, i_x_max):
                Emec_M[-1-l][c] = Emec_M[-1-l][c] + e_mec*min(contact.g1.etai_M[-1-l][c],contact.g2.etai_M[-1-l][c])
    #contact grain-wall part
    if dict_sollicitation['contact_gw_for_Emec']:
        for contact in dict_sample['L_contact_gw']:
            #extract a spatial zone
            x_min = min(contact.g.l_border_x)-dict_material['w']
            x_max = max(contact.g.l_border_x)+dict_material['w']
            y_min = min(contact.g.l_border_y)-dict_material['w']
            y_max = max(contact.g.l_border_y)+dict_material['w']
            #look for this part inside the global mesh
            #create search list
            x_L_search_min = abs(np.array(dict_sample['x_L'])-x_min)
            x_L_search_max = abs(np.array(dict_sample['x_L'])-x_max)
            y_L_search_min = abs(np.array(dict_sample['y_L'])-y_min)
            y_L_search_max = abs(np.array(dict_sample['y_L'])-y_max)
            #get index
            i_x_min = list(x_L_search_min).index(min(x_L_search_min))
            i_x_max = list(x_L_search_max).index(min(x_L_search_max))
            i_y_min = list(y_L_search_min).index(min(y_L_search_min))
            i_y_max = list(y_L_search_max).index(min(y_L_search_max))
            #Initialisation
            sum_min_etai = 0
            #compute the sum over the sample of the etai outside the box
            for l in range(i_y_min, i_y_max):
                for c in range(i_x_min, i_x_max):
                    if contact.nature == 'gwy_min' and dict_sample['y_L'][l] < contact.limit:
                        sum_min_etai = sum_min_etai + contact.g.etai_M[-1-l][c]
                    elif contact.nature == 'gwy_max' and dict_sample['y_L'][l] > contact.limit:
                        sum_min_etai = sum_min_etai + contact.g.etai_M[-1-l][c]
                    elif contact.nature == 'gwx_min' and dict_sample['x_L'][c] < contact.limit:
                        sum_min_etai = sum_min_etai + contact.g.etai_M[-1-l][c]
                    elif contact.nature == 'gwx_max' and dict_sample['x_L'][c] > contact.limit:
                        sum_min_etai = sum_min_etai + contact.g.etai_M[-1-l][c]
            #compute the variable e_mec
            e_mec = dict_sollicitation['alpha']/sum_min_etai*5*2*math.sqrt(2) #the term 5 is related to the factor applied to the spring grain - wall
                                                                              #the term 2*math.sqrt(2) is related to the equivalent radius and young modulus
            #compute the distribution of the mechanical energy
            for l in range(i_y_min, i_y_max):
                for c in range(i_x_min, i_x_max):
                    if contact.nature == 'gwy_min' and dict_sample['y_L'][l] < contact.limit:
                        Emec_M[-1-l][c] = Emec_M[-1-l][c] + e_mec*contact.g.etai_M[-1-l][c]
                    elif contact.nature == 'gwy_max' and dict_sample['y_L'][l] > contact.limit:
                        Emec_M[-1-l][c] = Emec_M[-1-l][c] + e_mec*contact.g.etai_M[-1-l][c]
                    elif contact.nature == 'gwx_min' and dict_sample['x_L'][c] < contact.limit:
                        Emec_M[-1-l][c] = Emec_M[-1-l][c] + e_mec*contact.g.etai_M[-1-l][c]
                    elif contact.nature == 'gwx_max' and dict_sample['x_L'][c] > contact.limit:
                        Emec_M[-1-l][c] = Emec_M[-1-l][c] + e_mec*contact.g.etai_M[-1-l][c]

    #Update element in dictionnary
    dict_sample['Emec_M'] = Emec_M

#-------------------------------------------------------------------------------

def Compute_kc_dil(dict_algorithm, dict_material, dict_sample):
    '''
    Compute the solute diffusion coefficient field in the sample.

    Here, a dilation method is applied. For all node, a Boolean variable is defined.
    This variable is True at least 2 eta_i are greater than 0.5 (in the contact zone).
                  is True all etai are lower than 0.5 (in the pore zone).
                  is False else.

    A dilation method is applied, the size of the structural element is the main case.

    The diffusion map is built on the Boolean map. If the variable is True, the diffusion is kc, else 0.

        Input :
            an algorithm dictionnary (a dict)
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
            n_etai_over_0_5 = 0 #counter
            for etai in dict_sample['L_etai'] :
                if etai.etai_M[-1-l][c] > 0.5 :
                    n_etai_over_0_5 = n_etai_over_0_5 + 1
            #at a contact
            if n_etai_over_0_5 >= 2:
                on_off_M[-l-1][c] = True
            #in the pore space
            elif n_etai_over_0_5 == 0:
                on_off_M[-l-1][c] = True

    #dilatation
    dilated_M = binary_dilation(on_off_M, dict_algorithm['struct_element'])

    #compute the map of the solute diffusion coefficient
    kc_M = np.array(np.zeros((len(dict_sample['y_L']),len(dict_sample['x_L']))))
    for l in range(len(dict_sample['y_L'])):
        for c in range(len(dict_sample['x_L'])):
            if dilated_M[-1-l][c] :
                kc_M[-1-l][c] = dict_material['kappa_c']

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
    for etai in dict_sample['L_etai']:
        for l in range(len(dict_sample['y_L'])):
            for c in range(len(dict_sample['x_L'])):
                sum_eta = sum_eta + etai.etai_M[l][c]

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
            Ed_mec = 0
            for etai in dict_sample['L_etai']:
                Ed_mec = Ed_mec + Ed_mec*(3*etai.etai_M[-1-l][c]**2-2*etai.etai_M[-1-l][c]**3)

            #Eche
            Ed_che = 0
            for etai in dict_sample['L_etai']:
                Ed_che = Ed_che + dict_sollicitation['chi']*dict_sample['solute_M'][-1-l][c]*(3*etai.etai_M[-1-l][c]**2-2*etai.etai_M[-1-l][c]**3)

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

#-------------------------------------------------------------------------------

def Compute_mean_sphericity(dict_algorithm, dict_sample):
    '''
    Compute the mean sphericities of the grains in the sample.

        Input :
            a sample dictionnary (a dict)
        Output :
            an area sphericity (a float)
            a diameter sphericity (a float)
            a circle ratio sphericity (a float)
            a perimeter sphericity (a float)
            a width to length ration sphericity (a float)
    '''
    L_area_sphericity = []
    L_diameter_sphericity = []
    L_circle_ratio_sphericity = []
    L_perimeter_sphericity = []
    L_width_to_length_ratio_sphericity = []
    for grain in dict_sample['L_g']:
        grain.geometric_study(dict_sample)
        grain.Compute_sphericity(dict_algorithm)
        #sphericities
        L_area_sphericity.append(grain.area_sphericity)
        L_diameter_sphericity.append(grain.diameter_sphericity)
        L_circle_ratio_sphericity.append(grain.circle_ratio_sphericity)
        L_perimeter_sphericity.append(grain.perimeter_sphericity)
        L_width_to_length_ratio_sphericity.append(grain.width_to_length_ratio_sphericity)

    return np.mean(L_area_sphericity), np.mean(L_diameter_sphericity), np.mean(L_circle_ratio_sphericity), np.mean(L_perimeter_sphericity), np.mean(L_width_to_length_ratio_sphericity)
