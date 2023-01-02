# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

The goal of this file is to define a new class.
The new class is about the unconserved variables used in phase-field simulation.
"""

#-------------------------------------------------------------------------------
#Libs
#-------------------------------------------------------------------------------

import numpy as np
import random

#Own
import Grain

#-------------------------------------------------------------------------------
#Class
#-------------------------------------------------------------------------------

class Etai:

    #-------------------------------------------------------------------------------
    def __init__(self, ID, L_ig, L_g):
        '''
        Defining the phase variables etai.

        One eta represents several grains.

            Input :
                an id (an integer)
                a list of grain id represented by the variable (a list)
                the list of the grains (a list)
            Output :
                an eta (an eta)
        '''
        self.id = ID
        self.l_ig = L_ig
        self.update_etai_M(L_g)

    #-------------------------------------------------------------------------------

    def update_etai_M(self,L_g):
        '''
        Sum variables field of all the grain assigned to this eta.

            Input :
                itself (an eta)
                the list of grains (a list)
            Output :
                Nothing, but the phase field attribute is updated (a nx x ny numpy array)
        '''
        self.etai_M = np.array(L_g[self.l_ig[0]].etai_M.copy())
        for i in range(1,len(self.l_ig)):
            self.etai_M = self.etai_M + np.array(L_g[self.l_ig[i]].etai_M.copy())

#-------------------------------------------------------------------------------
#Function
#-------------------------------------------------------------------------------

def etai_distribution(dict_algorithm, dict_sample):
    '''
    Assign grains to etai.

    A minimal distance between grains with same eta is computed from the factor variable

        Input :
            an algorithm dictionnary (a dict)
            a sample dictionnary (a dict)
        Output :
            Nothing, but the sample dictionnary is updated with the list of the etai (a list)
    '''
    #initialisation
    for grain in dict_sample['L_g']:
        grain.ig_near_L = []
        grain.eta_near_L = []
    #look for grains near
    for i_grain1 in range(1,len(dict_sample['L_g'])):
        grain1 = dict_sample['L_g'][i_grain1]
        for i_grain2 in range(0,i_grain1):
            grain2 = dict_sample['L_g'][i_grain2]
            if np.linalg.norm(grain1.center - grain2.center) < dict_algorithm['factor_etai'] * (grain1.r_max+grain2.r_max):
                grain1.ig_near_L.append(i_grain2)
                grain2.ig_near_L.append(i_grain1)
    #first try
    n_etai_L = [1]
    L_etai = [0]
    L_ig_etai = [[0]]
    dict_sample['L_g'][0].id_eta = 0
    for grain in dict_sample['L_g']:
        if 0 in grain.ig_near_L:
            grain.eta_near_L.append(0)
    for i_grain in range(1,len(dict_sample['L_g'])):
        grain = dict_sample['L_g'][i_grain]
        etai_defined = False
        etai = 0
        while etai < len(L_etai) and not etai_defined:
            if etai not in grain.eta_near_L:
                grain.id_eta = etai
                etai_defined = True
                n_etai_L[etai] = n_etai_L[etai] + 1
                L_ig_etai[etai].append(i_grain)
                for grain2 in dict_sample['L_g']:
                    if i_grain in grain2.ig_near_L:
                        grain2.eta_near_L.append(etai)
            etai = etai + 1
        if etai == len(L_etai) and not etai_defined:
            grain.id_eta = etai
            L_etai.append(etai)
            n_etai_L.append(1)
            L_ig_etai.append([i_grain])
            for grain2 in dict_sample['L_g']:
                if i_grain in grain2.ig_near_L:
                    grain2.eta_near_L.append(etai)
    #adaptation (try to have the same number of grain assigned to all eta)
    n_etai_mean = np.mean(n_etai_L)
    adaptation_done = False
    adaptation_i = 0
    while not adaptation_done :
        adaptation_i = adaptation_i + 1
        L_ig_over =  L_ig_etai[n_etai_L.index(max(n_etai_L))]
        id_g_to_work = L_ig_over[random.randint(0,len(L_ig_over)-1)]
        etai_over = n_etai_L.index(max(n_etai_L))
        etai_under = n_etai_L.index(min(n_etai_L))
        grain = dict_sample['L_g'][id_g_to_work]
        if etai_under not in grain.eta_near_L:
            grain.id_eta = etai_under
            n_etai_L[etai_over] = n_etai_L[etai_over] - 1
            n_etai_L[etai_under] = n_etai_L[etai_under] + 1
            L_ig_etai[etai_over].remove(id_g_to_work)
            L_ig_etai[etai_under].append(id_g_to_work)
            for grain2 in dict_sample['L_g']:
                if id_g_to_work in grain2.ig_near_L:
                    grain2.eta_near_L.remove(etai_over)
                    grain2.eta_near_L.append(etai_under)
        #check the quality
        adaptation_done = True
        for n_etai in n_etai_L :
            if n_etai_mean-2 < n_etai and n_etai < n_etai_mean+2:
                adaptation_done = False
        if adaptation_i > len(dict_sample['L_g']):
            adaptation_done = True

    #Create etai
    L_etai = []
    for i in range(len(L_ig_etai)) :
        etai = Etai(i, L_ig_etai, dict_sample['L_g'])
        L_etai.append(etai)

    #update the dict
    dict_sample['L_etai'] = L_etai

#-------------------------------------------------------------------------------
