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

name_to_load = '../Data_2G_ACS/PS_Adaptative_dt_PF_3_save'
name_report = 'Debug/Report_after_crash'

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

#-------------------------------------------------------------------------------
#Plan the simulation
#-------------------------------------------------------------------------------

#create a simulation report
simulation_report = Report.Report(name_report,datetime.now())

#delete last folder
if Path('Output/Ite_'+str(dict_algorithm['i_PFDEM']+1)).exists():
    shutil.rmtree('Output/Ite_'+str(dict_algorithm['i_PFDEM']+1))
if 'Diff_Solute' in dict_algorithm['L_flag_plot']:
    if Path('Debug/Diff_Solute/Ite_'+str(dict_algorithm['i_PFDEM']+1)).exists():
        shutil.rmtree('Debug/Diff_Solute/Ite_'+str(dict_algorithm['i_PFDEM']+1))

#-------------------------------------------------------------------------------
#main
#-------------------------------------------------------------------------------

while not User.Criteria_StopSimulation(dict_algorithm):

    main.iteration_main(dict_algorithm, dict_material, dict_sample, dict_sollicitation, dict_tracker, simulation_report)

#-------------------------------------------------------------------------------
#close simulation
#-------------------------------------------------------------------------------

main.close_main(dict_algorithm, dict_material, dict_sample, dict_sollicitation, dict_tracker, simulation_report)
