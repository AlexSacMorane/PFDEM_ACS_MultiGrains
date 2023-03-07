# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This file contains the different functions used to save in the simulation.
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

import os
import pickle

#-------------------------------------------------------------------------------

def save_dicts_tempo(dict_algorithm, dict_material, dict_sample, dict_sollicitation, dict_tracker, simulation_report):
    '''
    Save dictionnaries at the end of each PFDEM interations.

        Input :
            an algorithm dictionnary (a dictionnary)
            a material dictionnary (a dictionnary)
            a sample dictionnary (a dictionnary)
            a sollicitation dictionnary (a dictionnary)
            a tracker dictionnary (a dictionnary)
        Output :
            Nothing but a save file is generated (a file)
    '''
    outfile = open(dict_algorithm['namefile']+'_save_tempo','wb')
    dict_save = {}
    dict_save['algorithm'] = dict_algorithm
    dict_save['material'] = dict_material
    dict_save['sample'] = dict_sample
    dict_save['sollicitation'] = dict_sollicitation
    dict_save['tracker'] = dict_tracker
    dict_save['report'] = simulation_report
    pickle.dump(dict_save,outfile)
    outfile.close()

#-------------------------------------------------------------------------------

def save_dicts_ite(dict_algorithm, dict_material, dict_sample, dict_sollicitation, dict_tracker, simulation_report):
    '''
    Save dictionnaries before the dem simulation.

        Input :
            an algorithm dictionnary (a dictionnary)
            a material dictionnary (a dictionnary)
            a sample dictionnary (a dictionnary)
            a sollicitation dictionnary (a dictionnary)
            a tracker dictionnary (a dictionnary)
        Output :
            Nothing but a save file is generated (a file)
    '''
    outfile = open('Dicts/'+dict_algorithm['namefile']+'_'+str(dict_algorithm['i_PFDEM'])+'_save_tempo','wb')
    dict_save = {}
    dict_save['algorithm'] = dict_algorithm
    dict_save['material'] = dict_material
    dict_save['sample'] = dict_sample
    dict_save['sollicitation'] = dict_sollicitation
    dict_save['tracker'] = dict_tracker
    dict_save['report'] = simulation_report
    pickle.dump(dict_save,outfile)
    outfile.close()

#-------------------------------------------------------------------------------

def save_dicts_before_pf(dict_algorithm, dict_material, dict_sample, dict_sollicitation, dict_tracker, simulation_report):
    '''
    Save dictionnaries before the phase-field simulation.

        Input :
            an algorithm dictionnary (a dictionnary)
            a material dictionnary (a dictionnary)
            a sample dictionnary (a dictionnary)
            a sollicitation dictionnary (a dictionnary)
            a tracker dictionnary (a dictionnary)
        Output :
            Nothing but a save file is generated (a file)
    '''
    outfile = open(dict_algorithm['namefile']+'_save_tempo_before_pf','wb')
    dict_save = {}
    dict_save['algorithm'] = dict_algorithm
    dict_save['material'] = dict_material
    dict_save['sample'] = dict_sample
    dict_save['sollicitation'] = dict_sollicitation
    dict_save['tracker'] = dict_tracker
    dict_save['report'] = simulation_report
    pickle.dump(dict_save,outfile)
    outfile.close()

#-------------------------------------------------------------------------------

def save_dicts_final(dict_algorithm, dict_material, dict_sample, dict_sollicitation, dict_tracker, simulation_report):
    '''
    Save dictionnaries at the end of the simulation.

        Input :
            an algorithm dictionnary (a dictionnary)
            a material dictionnary (a dictionnary)
            a sample dictionnary (a dictionnary)
            a sollicitation dictionnary (a dictionnary)
            a tracker dictionnary (a dictionnary)
        Output :
            Nothing but a save file is generated (a file)
    '''
    outfile = open('../'+dict_algorithm['foldername']+'/'+dict_algorithm['namefile']+'_save','wb')
    dict_save = {}
    dict_save['algorithm'] = dict_algorithm
    dict_save['material'] = dict_material
    dict_save['sample'] = dict_sample
    dict_save['sollicitation'] = dict_sollicitation
    dict_save['tracker'] = dict_tracker
    dict_save['report'] = simulation_report
    pickle.dump(dict_save,outfile)
    outfile.close()
    os.remove(dict_algorithm['namefile']+'_save_tempo')
    os.remove(dict_algorithm['namefile']+'_save_tempo_before_pf')
