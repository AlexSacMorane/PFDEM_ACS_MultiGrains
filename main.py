# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This is the main file.
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

#Own functions and classes
import Grain
import Owntools
import Owntools.Compute
import Owntools.PFtoDEM_Multi
import Owntools.Plot
import Owntools.Save
import Owntools.Write
import User
import Report

#-------------------------------------------------------------------------------

def iteration_main(dict_algorithm, dict_material, dict_sample, dict_sollicitation, dict_tracker, simulation_report):
    '''
    Description of one PDEM iteration.

    The iteration is composed by a DEM step (to obtain a steady state configuration) and a PF step (to obtain dissolution and precipitation).

        Input :
            an algorithm dictionnary (a dict)
            a material dictionnary (a dict)
            a sample dictionnary (a dict)
            a sollicitation dictionnary (a dict)
            a tracker dictionnary (a dict)
            a simulation report (a Report)
        Output :
            Nothing but the dictionnaies and the report are updated
    '''
    #---------------------------------------------------------------------------
    #prepare iteration
    #---------------------------------------------------------------------------
    simulation_report.tic_tempo(datetime.now())
    dict_algorithm['i_PFDEM'] = dict_algorithm['i_PFDEM'] + 1
    simulation_report.write_and_print(f"\nITERATION {dict_algorithm['i_PFDEM']} / {dict_algorithm['n_t_PFDEM']}\n\n",f"\nITERATION {dict_algorithm['i_PFDEM']} / {dict_algorithm['n_t_PFDEM']}\n")
    os.mkdir('Output/Ite_'+str(dict_algorithm['i_PFDEM']))

    #---------------------------------------------------------------------------
    #move grains
    #---------------------------------------------------------------------------

    Grain.Compute_overlap_2_grains(dict_sample)
    Grain.Apply_overlap_target(dict_material,dict_sample,dict_sollicitation,dict_tracker)

    #---------------------------------------------------------------------------
    #prepare phase field simulation
    #---------------------------------------------------------------------------

    #Compute parameters needed
    Owntools.Compute.Compute_S_int(dict_sample) #the intersection surface
    Owntools.Compute.Compute_sum_min_etai(dict_sample, dict_sollicitation) #the sum of the minimum of etai
    Owntools.Compute.Compute_Emec(dict_sample, dict_sollicitation) #the mechanical energy
    if dict_material['method_to_compute_kc'] == 'dilation':
        Owntools.Compute.Compute_kc_dil(dict_material, dict_sample) #the solute diffusion
    elif dict_material['method_to_compute_kc'] == 'wfd':
        Owntools.Compute.Compute_kc_wfd(dict_material, dict_sample) #the solute diffusion
    if dict_material['method_to_compute_kc'] == 'interpolation':
        Owntools.Compute.Compute_kc_int(dict_material, dict_sample) #the solute diffusion
    dict_tracker['S_int_L'].append(dict_sample['S_int'])
    dict_tracker['sum_min_etai_L'].append(dict_sample['sum_min_etai'])

    #compute for total energy in the sample and track the value
    Owntools.Compute.Compute_sum_Ed_plus_minus(dict_sample, dict_sollicitation)
    dict_tracker['sum_ed_L'].append(dict_sample['sum_ed'])
    dict_tracker['sum_Ed_che_L'].append(dict_sample['sum_Ed_che'])
    dict_tracker['sum_Ed_mec_L'].append(dict_sample['sum_Ed_mec'])
    dict_tracker['sum_ed_plus_L'].append(dict_sample['sum_ed_plus'])
    dict_tracker['sum_ed_minus_L'].append(dict_sample['sum_ed_minus'])

    #Adaptative time step
    if abs(dict_sample['sum_ed']) < dict_algorithm['Ed_level1']:
        dict_algorithm['dt_PF'] = dict_algorithm['dt_PF_init']
    elif dict_algorithm['Ed_level1'] <= abs(dict_sample['sum_ed']) and abs(dict_sample['sum_ed']) < dict_algorithm['Ed_level2']:
        dict_algorithm['dt_PF'] = dict_algorithm['dt_PF_level1']
    elif dict_algorithm['Ed_level2'] <= abs(dict_sample['sum_ed']) and abs(dict_sample['sum_ed']) < dict_algorithm['Ed_level3']:
        dict_algorithm['dt_PF'] = dict_algorithm['dt_PF_level2']
    elif dict_algorithm['Ed_level3'] <= abs(dict_sample['sum_ed']) :
        dict_algorithm['dt_PF'] = dict_algorithm['dt_PF_level3']

    #plot
    if 'Config' in dict_algorithm['L_flag_plot']:
        Owntools.Plot.Plot_config(dict_algorithm, dict_sample)
    if 'Kc' in dict_algorithm['L_flag_plot']:
        Owntools.Plot.Plot_kc(dict_sample)

    #write data
    Owntools.Write.Write_eta_txt(dict_algorithm, dict_sample)
    Owntools.Write.Write_solute_txt(dict_algorithm, dict_sample)
    Owntools.Write.Write_kc_txt(dict_algorithm, dict_sample)
    Owntools.Write.Write_Emec_txt(dict_algorithm, dict_sample)

    #plot the difference of solute conentration in the case of a pure diffusion problem
    if 'Diff_Solute' in dict_algorithm['L_flag_plot']:
        os.mkdir('Debug/Diff_Solute/Ite_'+str(dict_algorithm['i_PFDEM']))
        Owntools.Plot.Plot_Diffusion_Solute(dict_algorithm, dict_material, dict_sample)

    #create i
    Owntools.Write.Write_i(dict_algorithm, dict_material, dict_sample, dict_sollicitation)

    simulation_report.tac_tempo(datetime.now(),f"Iteration {dict_algorithm['i_PFDEM']}: preparation of the pf simulation")
    simulation_report.tic_tempo(datetime.now())

    #---------------------------------------------------------------------------
    #PF simulation
    #---------------------------------------------------------------------------

    #run
    os.system('mpiexec -n '+str(dict_algorithm['np_proc'])+' ~/projects/moose/modules/combined/combined-opt -i '+dict_algorithm['namefile']+'_'+str(dict_algorithm['i_PFDEM'])+'.i')

    simulation_report.tac_tempo(datetime.now(),f"Iteration {dict_algorithm['i_PFDEM']}: pf simulation")
    simulation_report.tic_tempo(datetime.now())

    #sorting files
    j_str = Owntools.Sort_Files(dict_algorithm)

    #---------------------------------------------------------------------------
    #PF to DEM
    #---------------------------------------------------------------------------

    #look for the new grains shape
    for grain in dict_sample['L_g']:
        grain.PFtoDEM_Multi('Output/Ite_'+str(dict_algorithm['i_PFDEM'])+'/'+dict_algorithm['namefile']+'_'+str(dict_algorithm['i_PFDEM'])+'_other_'+j_str,dict_algorithm,dict_sample)
        grain.geometric_study(dict_sample)
    #look for the new solute shape
    Owntools.PFtoDEM_Multi.solute_PFtoDEM_Multi('Output/Ite_'+str(dict_algorithm['i_PFDEM'])+'/'+dict_algorithm['namefile']+'_'+str(dict_algorithm['i_PFDEM'])+'_other_'+j_str,dict_algorithm,dict_sample)
    #look for the initial external energy sources
    Owntools.PFtoDEM_Multi.Ed_PFtoDEM_Multi('Output/Ite_'+str(dict_algorithm['i_PFDEM'])+'/'+dict_algorithm['namefile']+'_'+str(dict_algorithm['i_PFDEM'])+'_other_000',dict_algorithm,dict_sample)

    #plot
    if 'Config' in dict_algorithm['L_flag_plot']:
        Owntools.Plot.Plot_config(dict_algorithm, dict_sample)
    if 'Init_Current_Shape' in dict_algorithm['L_flag_plot']:
        Owntools.Plot.Plot_init_current_shape(dict_sample)
    if 'Ed' in dict_algorithm['L_flag_plot']:
        Owntools.Plot.Plot_Ed(dict_sample)

    #---------------------------------------------------------------------------
    #postprocess
    #---------------------------------------------------------------------------

    #Compute the sphericity for the first grain
    dict_sample['L_g'][0].Compute_sphericity(dict_algorithm)

    #compute the mass of grain
    Owntools.Compute.Compute_sum_eta(dict_sample)

    #compute the mass of the solute
    Owntools.Compute.Compute_sum_c(dict_sample)

    #---------------------------------------------------------------------------
    #tracker
    #---------------------------------------------------------------------------

    dict_tracker['L_t'].append(dict_tracker['L_t'][-1]+dict_algorithm['dt_PF']*dict_algorithm['n_t_PF'])
    dict_tracker['L_dt'].append(dict_algorithm['dt_PF'])
    dict_tracker['L_sum_solute'].append(dict_sample['sum_c'])
    dict_tracker['L_sum_eta'].append(dict_sample['sum_eta'])
    dict_tracker['L_sum_total'].append(dict_sample['sum_c']+dict_sample['sum_eta'])
    dict_tracker['L_area_sphericity_g0'].append(dict_sample['L_g'][0].area_sphericity)
    dict_tracker['L_diameter_sphericity_g0'].append(dict_sample['L_g'][0].diameter_sphericity)
    dict_tracker['L_circle_ratio_sphericity_g0'].append(dict_sample['L_g'][0].circle_ratio_sphericity)
    dict_tracker['L_perimeter_sphericity_g0'].append(dict_sample['L_g'][0].perimeter_sphericity)
    dict_tracker['L_width_to_length_ratio_sphericity_g0'].append(dict_sample['L_g'][0].width_to_length_ratio_sphericity)
    dict_tracker['c_at_the_center'].append(Owntools.Extract_solute_at_p(dict_sample,(int(len(dict_sample['y_L'])/2),int(len(dict_sample['x_L'])/2))))

    #Plot trackers
    if 'Eta_c' in dict_algorithm['L_flag_plot'] :
        Owntools.Plot.Plot_sum_eta_c(dict_tracker)
    if 'Sphericity' in dict_algorithm['L_flag_plot'] :
        Owntools.Plot.Plot_sphericity(dict_tracker)
    if 'C_at_P' in dict_algorithm['L_flag_plot'] :
        Owntools.Plot.Plot_c_at_p(dict_tracker)
    if 'sum_Ed' in dict_algorithm['L_flag_plot'] :
        Owntools.Plot.Plot_sum_Ed(dict_tracker)
    if 'Sint_MinEtai' in dict_algorithm['L_flag_plot'] :
        Owntools.Plot.Plot_Sint_SumMinEtai(dict_tracker)
    if 'dt' in dict_algorithm['L_flag_plot'] :
        Owntools.Plot.Plot_dt_used(dict_tracker)

    #---------------------------------------------------------------------------
    #tempo save
    #---------------------------------------------------------------------------

    if dict_algorithm['SaveData']:
        Owntools.Save.save_dicts_tempo(dict_algorithm, dict_material, dict_sample, dict_sollicitation, dict_tracker)
        shutil.copy('Debug/Report.txt','../'+dict_algorithm['foldername']+'/Report_'+dict_algorithm['namefile']+'_tempo.txt')

    simulation_report.tac_tempo(datetime.now(),f"Iteration {dict_algorithm['i_PFDEM']}: from pf to dem")

#-------------------------------------------------------------------------------

def close_main(dict_algorithm, dict_material, dict_sample, dict_sollicitation, dict_tracker, simulation_report):
    '''
    Close the PFDEM.

        Input :
            an algorithm dictionnary (a dict)
            a material dictionnary (a dict)
            a sample dictionnary (a dict)
            a sollicitation dictionnary (a dict)
            a tracker dictionnary (a dict)
            a simulation report (a Report)
        Output :
            Nothing but the dictionnaries and the report are updated
    '''
    #make movie of the different configuration
    if 'Movie' in dict_algorithm['L_flag_plot'] and 'Config' in dict_algorithm['L_flag_plot']:
        Owntools.Plot.Plot_mp4('Debug/Configuration/Configuration_','Debug/Configuration.mp4')

    simulation_report.end(datetime.now())

    #final save
    if dict_algorithm['SaveData']:
        print()
        print('Copying data, it can take long times...')

        Owntools.Save.save_dicts_final(dict_algorithm, dict_material, dict_sample, dict_sollicitation, dict_tracker)
        name_actual_folder = os.path.dirname(os.path.realpath(__file__))
        shutil.copytree(name_actual_folder, '../'+dict_algorithm['foldername']+'/'+dict_algorithm['namefile'])
        os.remove('../'+dict_algorithm['foldername']+'/User_'+dict_algorithm['namefile']+'_tempo.txt')
        os.remove('../'+dict_algorithm['foldername']+'/Report_'+dict_algorithm['namefile']+'_tempo.txt')

#-------------------------------------------------------------------------------

if '__main__' == __name__:
    #-------------------------------------------------------------------------------
    #Plan simulation
    #-------------------------------------------------------------------------------

    if Path('Input').exists():
        shutil.rmtree('Input')
    os.mkdir('Input')
    if Path('Output').exists():
        shutil.rmtree('Output')
    os.mkdir('Output')
    if Path('Data').exists():
        shutil.rmtree('Data')
    os.mkdir('Data')
    if Path('Debug').exists():
        shutil.rmtree('Debug')
    os.mkdir('Debug')

    #-------------------------------------------------------------------------------
    #Create a simulation
    #-------------------------------------------------------------------------------

    #create a simulation report
    simulation_report = Report.Report('Debug/Report',datetime.now())
    simulation_report.tic_tempo(datetime.now())

    #general parameters
    dict_algorithm, dict_material, dict_sample, dict_sollicitation = User.All_parameters()
    if dict_algorithm['SaveData']:
        if not Path('../'+dict_algorithm['foldername']).exists():
            os.mkdir('../'+dict_algorithm['foldername'])
        #tempo save of the user file
        shutil.copy('User.py','../'+dict_algorithm['foldername']+'/User_'+dict_algorithm['namefile']+'_tempo.txt')

    #prepare plot
    if 'Config' in dict_algorithm['L_flag_plot']:
        os.mkdir('Debug/Configuration')
    if 'Init_Current_Shape' in dict_algorithm['L_flag_plot']:
        os.mkdir('Debug/Comparison_Init_Current')
    if 'Ed' in dict_algorithm['L_flag_plot']:
        os.mkdir('Debug/Ed')
    if 'Kc' in dict_algorithm['L_flag_plot']:
        os.mkdir('Debug/Kc')
    if 'Diff_Solute' in dict_algorithm['L_flag_plot']:
        os.mkdir('Debug/Diff_Solute')

    #create the two grains
    User.Add_2grains(dict_material,dict_sample)
    #Compute initial sum_eta
    Owntools.Compute.Compute_sum_eta(dict_sample)
    #Compute the sphericity initially for the first grain
    dict_sample['L_g'][0].geometric_study(dict_sample)
    dict_sample['L_g'][0].Compute_sphericity(dict_algorithm)
    #create the solute
    User.Add_solute(dict_sample)

    #plot
    if 'Config' in dict_algorithm['L_flag_plot']:
        Owntools.Plot.Plot_config(dict_algorithm, dict_sample)

    simulation_report.tac_tempo(datetime.now(),'Initialisation')

    #-------------------------------------------------------------------------------
    #trackers
    #-------------------------------------------------------------------------------

    dict_tracker = {
    'L_t' : [0],
    'L_dt' : [],
    'L_displacement' : [0],
    'L_int_displacement' : [0],
    'L_sum_solute' : [0],
    'L_sum_eta' : [dict_sample['sum_eta']],
    'L_sum_total' : [dict_sample['sum_eta']],
    'L_area_sphericity_g0' : [dict_sample['L_g'][0].area_sphericity],
    'L_diameter_sphericity_g0' : [dict_sample['L_g'][0].diameter_sphericity],
    'L_circle_ratio_sphericity_g0' : [dict_sample['L_g'][0].circle_ratio_sphericity],
    'L_perimeter_sphericity_g0' : [dict_sample['L_g'][0].perimeter_sphericity],
    'L_width_to_length_ratio_sphericity_g0' : [dict_sample['L_g'][0].width_to_length_ratio_sphericity],
    'c_at_the_center' : [Owntools.Extract_solute_at_p(dict_sample,(int(len(dict_sample['y_L'])/2),int(len(dict_sample['x_L'])/2)))],
    'sum_ed_L': [],
    'sum_Ed_che_L': [],
    'sum_Ed_mec_L': [],
    'sum_ed_plus_L' : [],
    'sum_ed_minus_L' : [],
    'S_int_L' : [],
    'sum_min_etai_L' : []
    }

    #-------------------------------------------------------------------------------
    #main
    #-------------------------------------------------------------------------------

    dict_algorithm['i_PFDEM'] = 0
    while not User.Criteria_StopSimulation(dict_algorithm):

        iteration_main(dict_algorithm, dict_material, dict_sample, dict_sollicitation, dict_tracker, simulation_report)

    #-------------------------------------------------------------------------------
    #close simulation
    #-------------------------------------------------------------------------------

    close_main(dict_algorithm, dict_material, dict_sample, dict_sollicitation, dict_tracker, simulation_report)
