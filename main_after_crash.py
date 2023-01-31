# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This is the file to restart a simulation after a crash.
There is a save at the end of each PFDEM iteration.
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

from pathlib import Path
from datetime import datetime
import numpy as np
import os
import shutil
import math
import pickle

#Own functions and classes
import Grain
import User
import Report
import main

#-------------------------------------------------------------------------------
#User
#-------------------------------------------------------------------------------

name_to_load = '../Data_MG_ACS/PS_MG_1_before_pf'

#-------------------------------------------------------------------------------
#load data
#-------------------------------------------------------------------------------

toload = open(name_to_load,'rb')
dict_save = pickle.load(toload,encoding = 'bytes')
toload.close()
dict_algorithm = dict_save['algorithm']
dict_material = dict_save['material']
dict_sample = dict_save['sample']
dict_sollicitation = dict_save['sollicitation']
dict_tracker = dict_save['tracker']
simulation_report = dict_save['report']

#-------------------------------------------------------------------------------
#Plan the simulation
#-------------------------------------------------------------------------------

simulation_report.write('\nA crash occurs...\n\n')

#delete last folder
if Path('Output/Ite_'+str(dict_algorithm['i_PFDEM']+1)).exists():
    shutil.rmtree('Output/Ite_'+str(dict_algorithm['i_PFDEM']+1))
if 'Diff_Solute' in dict_algorithm['L_flag_plot']:
    if Path('Debug/Diff_Solute/Ite_'+str(dict_algorithm['i_PFDEM']+1)).exists():
        shutil.rmtree('Debug/Diff_Solute/Ite_'+str(dict_algorithm['i_PFDEM']+1))

#-------------------------------------------------------------------------------
#main
#-------------------------------------------------------------------------------

#crash during the phase-field simulation
if name_to_load[-10:] == '_before_pf':
    main.iteration_main_from_pf(dict_algorithm, dict_material, dict_sample, dict_sollicitation, dict_tracker, simulation_report)

while not User.Criteria_StopSimulation(dict_algorithm):
    main.iteration_main_until_pf(dict_algorithm, dict_material, dict_sample, dict_sollicitation, dict_tracker, simulation_report)
    main.iteration_main_from_pf(dict_algorithm, dict_material, dict_sample, dict_sollicitation, dict_tracker, simulation_report)

#-------------------------------------------------------------------------------
#close simulation
#-------------------------------------------------------------------------------

main.close_main(dict_algorithm, dict_material, dict_sample, dict_sollicitation, dict_tracker, simulation_report)
