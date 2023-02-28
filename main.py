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
import Contact_gg
import Contact_gw
import Owntools
import Owntools.Compute
import Owntools.PFtoDEM_Multi
import Owntools.Plot
import Owntools.Save
import Owntools.Write
import Create_IC
import Create_IC.Grain_ic
import Create_IC.Contact_gg_ic
import Create_IC.Contact_gw_ic
import User
import Report
import Etai

#-------------------------------------------------------------------------------

def iteration_main_until_pf(dict_algorithm, dict_material, dict_sample, dict_sollicitation, dict_tracker, simulation_report):
    '''
    Description of one PDEM iteration until the phase-field simulation.

    The iteration is composed by a DEM step (to obtain a steady state configuration).

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
    if dict_algorithm['Debug_DEM']:
        os.mkdir('Debug/Configuration/PFDEM_'+str(dict_algorithm['i_PFDEM']))
        os.mkdir('Debug/txt/PFDEM_'+str(dict_algorithm['i_PFDEM']))

    # Saving center and angle to compute a rigid body motion
    L_center_g = []
    for grain in dict_sample['L_g']:
        L_center_g.append(grain.center.copy())
    L_theta_g = []
    for grain in dict_sample['L_g']:
        L_theta_g.append(grain.theta)

    #initial value
    DEM_loop_statut = True
    dict_algorithm['i_DEM'] = 0

    # Compute kinetic energy criteria and update element in dict
    Ecin_stop = 0
    for grain in dict_sample['L_g']:
        Ecin_stop = Ecin_stop + 0.5*grain.mass*(dict_algorithm['Ecin_ratio']*grain.r_mean/dict_algorithm['dt_DEM'])**2/len(dict_sample['L_g'])
    dict_algorithm['Ecin_stop'] = Ecin_stop

    # Initialize the DEM step with 0 speed and 0 angular speed
    for grain in dict_sample['L_g']:
        grain.v = np.array([0, 0])
        grain.w = 0

    #Trackers and add element in dict
    dict_tracker['Ecin'] = []
    dict_tracker['y_box_max_DEM'] = [dict_sample['y_box_max']]
    dict_tracker['Force_on_upper_wall'] = []

    #---------------------------------------------------------------------------
    #DEM to move grains
    #---------------------------------------------------------------------------

    while DEM_loop_statut :
        # update element in dict
        dict_algorithm['i_DEM'] = dict_algorithm['i_DEM'] + 1

        #Initialize the grain kinematic
        for grain in dict_sample['L_g']:
            grain.init_f_control(dict_sollicitation)

        # Detection of contacts between grains
        if (dict_algorithm['i_DEM']-1) % dict_algorithm['i_update_neighborhoods']  == 0:
            Contact_gg.Update_Neighborhoods(dict_algorithm,dict_sample)
        Contact_gg.Grains_Polyhedral_contact_Neighborhoods(dict_material,dict_sample)

        # Detection of contacts between grain and walls
        if (dict_algorithm['i_DEM']-1)% dict_algorithm['i_update_neighborhoods']  == 0:
            Contact_gw.Update_wall_Neighborhoods(dict_algorithm, dict_sample)
        Contact_gw.Grains_Polyhedral_Wall_contact_Neighborhood(dict_material,dict_sample)

        #Compute contact interactions (g-g and g-w)
        for contact in dict_sample['L_contact']:
            contact.DEM_gg_Polyhedral_normal()
            contact.DEM_gg_Polyhedral_tangential(dict_algorithm['dt_DEM'])
        for contact in dict_sample['L_contact_gw'] :
            contact.DEM_gw_Polyhedral_normal()
            contact.DEM_gw_Polyhedral_tangential(dict_algorithm['dt_DEM'])

        #Move particles and trackers
        #Semi implicit euler scheme
        Ecin = 0
        Force_applied = 0
        for grain in dict_sample['L_g']:
            a_i = grain.f/grain.mass
            v_i = grain.v + a_i*dict_algorithm['dt_DEM']
            dw_i = grain.mz/grain.inertia
            w_i = grain.w + dw_i*dict_algorithm['dt_DEM']
            grain.update_geometry_kinetic(v_i,a_i,w_i,dict_algorithm['dt_DEM']) #Move grains
            Ecin = Ecin + 0.5*grain.mass*np.linalg.norm(grain.v)**2/len(dict_sample['L_g'])

        #Control the y_max to verify vertical confinement
        Owntools.Control_y_max_NR(dict_sample,dict_sollicitation)
        #trackers
        dict_tracker['Ecin'].append(Ecin)
        dict_tracker['y_box_max_DEM'].append(dict_sample['y_box_max'])
        dict_tracker['Force_on_upper_wall'].append(dict_sollicitation['Force_on_upper_wall'])

        if (dict_algorithm['i_DEM']-1) % dict_algorithm['i_print_plot'] == 0:
            print('\nPF '+str(dict_algorithm['i_PFDEM'])+' -> i_DEM '+str(dict_algorithm['i_DEM'])+' / '+str(dict_algorithm['i_DEM_stop'])+' (max)')
            print('Ecin',int(Ecin),'/',int(dict_algorithm['Ecin_stop']),'('+str(int(100*Ecin/dict_algorithm['Ecin_stop'])),' %)')
            print('F_confinement',int(dict_sollicitation['Force_on_upper_wall']),'/',int(dict_sollicitation['Vertical_Confinement_Force']),'('+str(int(100*dict_sollicitation['Force_on_upper_wall']/dict_sollicitation['Vertical_Confinement_Force'])),' %)')

            if dict_algorithm['Debug_DEM'] :
                Owntools.Plot.Plot_config_DEM(dict_algorithm, dict_sample)
                Owntools.Write.Write_DEM_txt_DEM(dict_algorithm,dict_sample)

        #-----------------------------------------------------------------------------
        # Stop conditions
        #-----------------------------------------------------------------------------

        if dict_algorithm['i_DEM'] >= dict_algorithm['i_DEM_stop'] :
            DEM_loop_statut = False
            print("DEM loop stopped by too many iterations.")
            simulation_report.write('/!\ End of DEM steps with '+str(dict_algorithm['i_DEM']+1)+' iterations / '+str(dict_algorithm['i_DEM_stop']+1)+'/!\ \n')
        if dict_algorithm['i_DEM'] > max(0.1*dict_algorithm['i_DEM_stop'],dict_algorithm['n_window_stop']) and Ecin < dict_algorithm['Ecin_stop'] and (dict_sollicitation['Vertical_Confinement_Force']*0.95<dict_sollicitation['Force_on_upper_wall'] and dict_sollicitation['Force_on_upper_wall']<dict_sollicitation['Vertical_Confinement_Force']*1.05):
            y_box_max_window = dict_tracker['y_box_max_DEM'][dict_algorithm['i_DEM']+1-dict_algorithm['n_window_stop']:dict_algorithm['i_DEM']+1]
            if max(y_box_max_window) - min(y_box_max_window) < dict_algorithm['dy_box_max_stop']:
                DEM_loop_statut = False
                print("DEM loop stopped by steady state reached.")
                simulation_report.write("DEM loop stopped by steady state reached with "+str(dict_algorithm['i_DEM']+1)+' iterations / '+str(dict_algorithm['i_DEM_stop']+1)+"\n")

    #---------------------------------------------------------------------------
    #Close DEM step
    #---------------------------------------------------------------------------

    #Write DEM data in a .txt
    if 'DEM_txt' in dict_algorithm['L_flag_plot']:
        Owntools.Write.Write_DEM_txt(dict_algorithm,dict_sample)

    simulation_report.tac_tempo(datetime.now(),f"Iteration {dict_algorithm['i_PFDEM']}: dem")
    simulation_report.tic_tempo(datetime.now())

    #---------------------------------------------------------------------------
    #Compute and apply rigid boby motion
    #---------------------------------------------------------------------------

    mean_delta_sum_eta = 0
    for i_grain in range(len(dict_sample['L_g'])):
        U = dict_sample['L_g'][i_grain].center - L_center_g[i_grain] #translation
        dtheta = dict_sample['L_g'][i_grain].theta - L_theta_g[i_grain] #rotation
        dict_sample['L_g'][i_grain].move_grain_interpolation(U, dtheta, dict_material, dict_sample)
        mean_delta_sum_eta = mean_delta_sum_eta + abs(dict_sample['L_g'][i_grain].delta_sum_eta)
    mean_delta_sum_eta = mean_delta_sum_eta/len(dict_sample['L_g'])
    simulation_report.write_and_print('Mean Delta sum eta = '+str(round(mean_delta_sum_eta,2))+' %\n','Mean Delta sum eta = '+str(round(mean_delta_sum_eta,2))+' %')

    #---------------------------------------------------------------------------
    #Recreate the etai
    #---------------------------------------------------------------------------

    for etai in dict_sample['L_etai'] :
        etai.update_etai_M(dict_sample['L_g'])

    #---------------------------------------------------------------------------
    #prepare phase field simulation
    #---------------------------------------------------------------------------

    #plot
    if 'Config' in dict_algorithm['L_flag_plot']:
        Owntools.Plot.Plot_config(dict_algorithm, dict_sample)

    #Move out solute in grains
    Owntools.Interpolate_solute_out_grains(dict_algorithm, dict_sample)

    #Compute the mechanical energy term
    Owntools.Compute.Compute_Emec(dict_material, dict_sample, dict_sollicitation)

    #Compute the solute diffusion
    Owntools.Compute.Compute_kc_dil(dict_algorithm, dict_material, dict_sample) #the solute diffusion

    #compute for total energy in the sample and track the value
    Owntools.Compute.Compute_sum_Ed_plus_minus(dict_sample, dict_sollicitation)
    dict_tracker['sum_ed_L'].append(dict_sample['sum_ed'])
    dict_tracker['sum_Ed_che_L'].append(dict_sample['sum_Ed_che'])
    dict_tracker['sum_Ed_mec_L'].append(dict_sample['sum_Ed_mec'])
    dict_tracker['sum_ed_plus_L'].append(dict_sample['sum_ed_plus'])
    dict_tracker['sum_ed_minus_L'].append(dict_sample['sum_ed_minus'])

    #compute absolute total energy per contact node
    Owntools.Compute.Compute_Ed_abs_node_contact(dict_sample, dict_sollicitation)
    dict_tracker['sum_ed_abs_node_L'].append(dict_sample['sum_ed_abs_node'])

    #Adaptative time step
    if abs(dict_sample['sum_ed_abs_node']) < dict_algorithm['Ed_level1']:
        dict_algorithm['dt_PF'] = dict_algorithm['dt_PF_init']
    elif dict_algorithm['Ed_level1'] <= abs(dict_sample['sum_ed_abs_node']) and abs(dict_sample['sum_ed_abs_node']) < dict_algorithm['Ed_level2']:
        dict_algorithm['dt_PF'] = dict_algorithm['dt_PF_level1']
    elif dict_algorithm['Ed_level2'] <= abs(dict_sample['sum_ed_abs_node']) and abs(dict_sample['sum_ed_abs_node']) < dict_algorithm['Ed_level3']:
        dict_algorithm['dt_PF'] = dict_algorithm['dt_PF_level2']
    elif dict_algorithm['Ed_level3'] <= abs(dict_sample['sum_ed_abs_node']) :
        dict_algorithm['dt_PF'] = dict_algorithm['dt_PF_level3']

    #plot
    if 'Config' in dict_algorithm['L_flag_plot']:
        Owntools.Plot.Plot_config(dict_algorithm, dict_sample)
    if 'Kc' in dict_algorithm['L_flag_plot']:
        Owntools.Plot.Plot_kc(dict_sample)
    if 'DEM_tracker' in dict_algorithm['L_flag_plot']:
        Owntools.Plot.Plot_DEM_tracker(dict_tracker)
    if 'YBoxMax' in dict_algorithm['L_flag_plot']:
         Owntools.Plot.Plot_yboxmax(dict_tracker)

    #write data
    Owntools.Write.Write_eta_txt(dict_algorithm, dict_sample)
    Owntools.Write.Write_solute_txt(dict_algorithm, dict_sample)
    Owntools.Write.Write_kc_txt(dict_algorithm, dict_sample)
    Owntools.Write.Write_Emec_txt(dict_algorithm, dict_sample)

    #plot the difference of solute concentration in the case of a pure diffusion problem
    if 'Diff_Solute' in dict_algorithm['L_flag_plot']:
        raise ValueError('This plot is not adapted for multiple grains... WIP !')
        os.mkdir('Debug/Diff_Solute/Ite_'+str(dict_algorithm['i_PFDEM']))
        Owntools.Plot.Plot_Diffusion_Solute(dict_algorithm, dict_material, dict_sample)
    #look for the initial external energy sources
    if 'Ed' in dict_algorithm['L_flag_plot']:
        Owntools.Plot.Plot_Ed(dict_sample, dict_sollicitation)

    #create i
    Owntools.Write.Write_i(dict_algorithm, dict_material, dict_sample, dict_sollicitation)

    simulation_report.tac_tempo(datetime.now(),f"Iteration {dict_algorithm['i_PFDEM']}: preparation of the pf simulation")
    simulation_report.tic_tempo(datetime.now())

    #save
    Owntools.Save.save_dicts_before_pf(dict_algorithm, dict_material, dict_sample, dict_sollicitation, dict_tracker, simulation_report)

#-------------------------------------------------------------------------------

def iteration_main_from_pf(dict_algorithm, dict_material, dict_sample, dict_sollicitation, dict_tracker, simulation_report):
    '''
    Description of one PDEM iteration from the phase field simulation.

    The iteration is composed by a PF step (to obtain dissolution and precipitation).

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
    Etai.PFtoDEM_Multi('Output/Ite_'+str(dict_algorithm['i_PFDEM'])+'/'+dict_algorithm['namefile']+'_'+str(dict_algorithm['i_PFDEM'])+'_other_'+j_str, dict_algorithm, dict_material, dict_sample)

    #look for the new solute shape
    Owntools.PFtoDEM_Multi.solute_PFtoDEM_Multi('Output/Ite_'+str(dict_algorithm['i_PFDEM'])+'/'+dict_algorithm['namefile']+'_'+str(dict_algorithm['i_PFDEM'])+'_other_'+j_str,dict_algorithm,dict_sample)

    #plot
    if 'Config' in dict_algorithm['L_flag_plot']:
        Owntools.Plot.Plot_config(dict_algorithm, dict_sample)
    if 'Init_Current_Shape' in dict_algorithm['L_flag_plot']:
        Owntools.Plot.Plot_init_current_shape(dict_sample)

    #---------------------------------------------------------------------------
    #postprocess
    #---------------------------------------------------------------------------

    #porosity
    Owntools.Compute.Compute_porosity(dict_sample)

    #Compute the mean sphericities
    area_sphericity_mean, diameter_sphericity_mean, circle_ratio_sphericity_mean, perimeter_sphericity_mean, width_to_length_ratio_sphericity_mean = Owntools.Compute.Compute_mean_sphericity(dict_algorithm, dict_sample)

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
    dict_tracker['L_area_sphericity_mean'].append(area_sphericity_mean)
    dict_tracker['L_diameter_sphericity_mean'].append(diameter_sphericity_mean)
    dict_tracker['L_circle_ratio_sphericity_mean'].append(circle_ratio_sphericity_mean)
    dict_tracker['L_perimeter_sphericity_mean'].append(perimeter_sphericity_mean)
    dict_tracker['L_width_to_length_ratio_sphericity_mean'].append(width_to_length_ratio_sphericity_mean)
    dict_tracker['L_area_sphericity_g0'].append(dict_sample['L_g'][0].area_sphericity)
    dict_tracker['L_diameter_sphericity_g0'].append(dict_sample['L_g'][0].diameter_sphericity)
    dict_tracker['L_circle_ratio_sphericity_g0'].append(dict_sample['L_g'][0].circle_ratio_sphericity)
    dict_tracker['L_perimeter_sphericity_g0'].append(dict_sample['L_g'][0].perimeter_sphericity)
    dict_tracker['L_width_to_length_ratio_sphericity_g0'].append(dict_sample['L_g'][0].width_to_length_ratio_sphericity)
    dict_tracker['L_y_box_max'].append(dict_sample['y_box_max'])
    dict_tracker['L_porosity'].append(dict_sample['porosity'])

    #Plot trackers
    if 'Eta_c' in dict_algorithm['L_flag_plot'] :
        Owntools.Plot.Plot_sum_eta_c(dict_tracker)
    if 'Sphericity' in dict_algorithm['L_flag_plot'] :
        Owntools.Plot.Plot_sphericity(dict_tracker)
    if 'sum_Ed' in dict_algorithm['L_flag_plot'] :
        Owntools.Plot.Plot_sum_Ed(dict_tracker)
    if 'dt' in dict_algorithm['L_flag_plot'] :
        Owntools.Plot.Plot_dt_used(dict_sample, dict_tracker)
    if 'Porosity' in dict_algorithm['L_flag_plot'] :
        Owntools.Plot.Plot_porosity(dict_tracker)

    #---------------------------------------------------------------------------
    #tempo save
    #---------------------------------------------------------------------------

    Owntools.Save.save_dicts_tempo(dict_algorithm, dict_material, dict_sample, dict_sollicitation, dict_tracker, simulation_report)
    if dict_algorithm['SaveData']:
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

        #clean memory
        if dict_algorithm['clean_memory']:
            shutil.rmtree('Data')
            shutil.rmtree('Input')
            shutil.rmtree('Output')

        print()
        print('Copying data, it can take long times...')

        Owntools.Save.save_dicts_final(dict_algorithm, dict_material, dict_sample, dict_sollicitation, dict_tracker, simulation_report)
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
    dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitation = User.All_parameters()
    if dict_algorithm['SaveData']:
        if not Path('../'+dict_algorithm['foldername']).exists():
            os.mkdir('../'+dict_algorithm['foldername'])
        #tempo save of the user file
        shutil.copy('User.py','../'+dict_algorithm['foldername']+'/User_'+dict_algorithm['namefile']+'_tempo.txt')

    #prepare plot
    if 'Config' in dict_algorithm['L_flag_plot'] or dict_algorithm['Debug_DEM'] or dict_ic['Debug_DEM']:
        os.mkdir('Debug/Configuration')
        if dict_ic['Debug_DEM'] :
            os.mkdir('Debug/Configuration/Init')
    if 'DEM_txt' in dict_algorithm['L_flag_plot'] or dict_algorithm['Debug_DEM']:
        os.mkdir('Debug/txt')
    if 'Init_Current_Shape' in dict_algorithm['L_flag_plot']:
        os.mkdir('Debug/Comparison_Init_Current')
    if 'Ed' in dict_algorithm['L_flag_plot']:
        os.mkdir('Debug/Ed')
    if 'Kc' in dict_algorithm['L_flag_plot']:
        os.mkdir('Debug/Kc')
    if 'Diff_Solute' in dict_algorithm['L_flag_plot']:
        os.mkdir('Debug/Diff_Solute')
    if 'DEM_tracker' in dict_algorithm['L_flag_plot']:
        os.mkdir('Debug/DEM_tracker')

    #create the initial configuration
    Create_IC.LG_tempo(dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitation, simulation_report)

    #Define the mesh
    User.Add_mesh(dict_geometry, dict_sample)
    simulation_report.write('\nThe mesh used is '+str(len(dict_sample['x_L']))+' x '+str(len(dict_sample['y_L']))+' ( '+str(len(dict_sample['x_L'])*len(dict_sample['y_L']))+' nodes)\n\n')

    #Add needed variables
    User.Add_variables_needed(dict_geometry, dict_material, dict_sample, dict_sollicitation)

    #Convert the tempo grain to real grain
    Create_IC.From_LG_tempo_to_usable(dict_ic, dict_material, dict_sample)

    #plot mesh
    if 'Mesh' in dict_algorithm['L_flag_plot'] :
        Owntools.Plot.Plot_mesh(dict_sample)

    #Distribution of the etai and plot
    Etai.etai_distribution(dict_algorithm, dict_sample, simulation_report)
    if 'Etai_distribution' in dict_algorithm['L_flag_plot'] :
        Owntools.Plot.Plot_etai_distribution(dict_sample)

    #Compute initial sum_eta
    Owntools.Compute.Compute_sum_eta(dict_sample)
    #Compute the mean sphericities
    area_sphericity_mean, diameter_sphericity_mean, circle_ratio_sphericity_mean, perimeter_sphericity_mean, width_to_length_ratio_sphericity_mean = Owntools.Compute.Compute_mean_sphericity(dict_algorithm, dict_sample)
    #create the solute
    User.Add_solute(dict_sample)
    #compute the porosity
    Owntools.Compute.Compute_porosity(dict_sample)

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
    'L_y_box_max' : [dict_sample['y_box_max']],
    'L_porosity' : [dict_sample['porosity']],
    'L_sum_solute' : [0],
    'L_sum_eta' : [dict_sample['sum_eta']],
    'L_sum_total' : [dict_sample['sum_eta']],
    'L_area_sphericity_g0' : [dict_sample['L_g'][0].area_sphericity],
    'L_diameter_sphericity_g0' : [dict_sample['L_g'][0].diameter_sphericity],
    'L_circle_ratio_sphericity_g0' : [dict_sample['L_g'][0].circle_ratio_sphericity],
    'L_perimeter_sphericity_g0' : [dict_sample['L_g'][0].perimeter_sphericity],
    'L_width_to_length_ratio_sphericity_g0' : [dict_sample['L_g'][0].width_to_length_ratio_sphericity],
    'L_area_sphericity_mean' : [area_sphericity_mean],
    'L_diameter_sphericity_mean' : [diameter_sphericity_mean],
    'L_circle_ratio_sphericity_mean' : [circle_ratio_sphericity_mean],
    'L_perimeter_sphericity_mean' : [perimeter_sphericity_mean],
    'L_width_to_length_ratio_sphericity_mean' : [width_to_length_ratio_sphericity_mean],
    'sum_Ed_che_L': [],
    'sum_Ed_mec_L': [],
    'sum_ed_L': [],
    'sum_ed_plus_L' : [],
    'sum_ed_minus_L' : [],
    'sum_ed_abs_node_L' : []
    }

    #-------------------------------------------------------------------------------
    #main
    #-------------------------------------------------------------------------------

    # Preparation and add elements in dicts
    dict_algorithm['i_PFDEM'] = 0
    dict_sample['L_contact_gw'] = []
    dict_sample['L_ij_contact_gw'] = []
    dict_sample['id_contact_gw'] = 0
    dict_sample['L_contact'] = []
    dict_sample['L_ij_contact'] = []
    dict_sample['id_contact'] = 0

    while not User.Criteria_StopSimulation(dict_algorithm):

        iteration_main_until_pf(dict_algorithm, dict_material, dict_sample, dict_sollicitation, dict_tracker, simulation_report)
        iteration_main_from_pf(dict_algorithm, dict_material, dict_sample, dict_sollicitation, dict_tracker, simulation_report)

    #-------------------------------------------------------------------------------
    #close simulation
    #-------------------------------------------------------------------------------

    close_main(dict_algorithm, dict_material, dict_sample, dict_sollicitation, dict_tracker, simulation_report)
