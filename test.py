# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This is the test file to verify if the code is running well.
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

from pathlib import Path
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import unittest
import os
import shutil
import math

#own functions and classes
import User
import Owntools
import Owntools.Compute
import Owntools.Plot
import Owntools.Write
import Grain
import Report
from main import iteration_main

#-------------------------------------------------------------------------------
#Test
#-------------------------------------------------------------------------------

class TestExample(unittest.TestCase):
    '''
    This is a template to write new tests.
    '''
    def test_example(self):
        '''
        This is a template to write new tests.

            Output :
                The result is always True (a bool)
        '''
        self.assertTrue(True)

#-------------------------------------------------------------------------------

class TestGlobal(unittest.TestCase):
    '''
    Test global.
    '''
    def test_all_files_here(self):
        '''
        This test verifies that all files are in the directory.

            Output :
                The result is True if all files needed for a simulation are here (a bool)
        '''
        L_files = ['Grain.py',
                   'main_after_crash.py',
                   'main.py',
                   'Owntools/__init__.py',
                   'Owntools/Debug_Diff_Solute_base.i',
                   'Owntools/PFtoDEM_Multi.py',
                   'Owntools/Plot.py',
                   'Owntools/Save.py',
                   'Owntools/Write.py',
                   'PF_ACS_base.i',
                   'Report.py',
                   'User.py']
        AllHere = True
        MissingFiles = ''
        for file in L_files :
            plotpath = Path(file)
            if not plotpath.exists():
                MissingFiles = MissingFiles + ' '+ file
                AllHere = False
        self.assertTrue(AllHere,'Some files are missing :'+MissingFiles)

    #---------------------------------------------------------------------------

    def test_main_iteration_main(self):
        '''
        This test verifies that one PFDEM iteration can be done.

            Output :
                No result, except if an error appears
        '''
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
        simulation_report.tac_tempo(datetime.now(),'Initialisation')

        #trackers
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

        #Try to run one iteration
        dict_algorithm['i_PFDEM'] = 0
        iteration_main(dict_algorithm, dict_material, dict_sample, dict_sollicitation, dict_tracker, simulation_report)

#-------------------------------------------------------------------------------

class TestReport(unittest.TestCase):
    '''Test functions from Report.py.'''
    def test_Report(self):
        '''
        Try to generate a report.txt file.

            Output :
                The result depends on the fact if the file is well generated or not (a bool)
        '''
        #try to create a report.txt file
        simulation_report = Report.Report('Report',datetime.now())
        #check if the .txt has been created
        self.assertTrue(Path('Report.txt').is_file(),"The file Report.txt has not been created by the function Report.Report()!")
        os.remove('Report.txt')

#-------------------------------------------------------------------------------

class TestUser(unittest.TestCase):
    '''Test functions from User.py.'''
    def test_All_parameters(self):
        '''
        Try to acquire data from the function User.All_parameters().

            Output :
                The result depends on the fact if data are well acquired or not (a bool)
        '''
        #try to acquire data
        dict_algorithm, dict_material, dict_sample, dict_sollicitation = User.All_parameters()
        #check if all data are dictionnaries
        self.assertTrue((type(dict_algorithm) is dict) and (type(dict_material) is dict) and (type(dict_sample) is dict) and (type(dict_sollicitation) is dict),'Outputs from User.All_parameters() are not dictionaries!')

    #---------------------------------------------------------------------------

    def test_Add_2grains(self):
        '''
        Try to generate two grains with the function User.Add_2grains().

            Output :
                The result depends on the fact if grains are well generated or not (a bool)
        '''
        #Acquire data
        dict_algorithm, dict_material, dict_sample, dict_sollicitation = User.All_parameters()
        #try to create 2 grains
        User.Add_2grains(dict_material,dict_sample)
        #check if there are 2 grains
        self.assertTrue(len(dict_sample['L_g'])==2,'The function User.Add_2grains does not create 2 grains!')

    #---------------------------------------------------------------------------

    def test_Add_solute(self):
        '''
        Try to generate a solute with the function User.Add_solute().

            Output :
                The result depends on the fact if solute is well generated or not (a bool)
        '''
        #Acquire data
        dict_algorithm, dict_material, dict_sample, dict_sollicitation = User.All_parameters()
        #try to create 2 grains
        User.Add_solute(dict_sample)
        #check if there are 2 grains
        self.assertTrue('solute_M' in dict_sample.keys(),'The function User.Add_solute does not create a solute!')

#-------------------------------------------------------------------------------

class TestOwntools(unittest.TestCase):
    '''Test functions from Owntools.py.'''
    def test_is_PF_ACS_base_here(self):
        '''
        Verify if the template PF_ACS_base.i is in the directory.

        This file is used to generated MOOSE simulation input file.

            Output :
                The result depends on the fact if the file is here or not (a bool)
        '''
        #Check if the file PF_ACS_base.i is in the directory
        self.assertTrue(Path('PF_ACS_base.i').is_file(),"The file PF_CH_AC_base.i should exists!")

    #---------------------------------------------------------------------------

    def test_Write_i(self):
        '''
        Try to create a MOOSE simulation input file with Owntools.Write.Write_i().

        This file is created from the template PF_ACS_base.i.

            Output :
                The result depends on the fact if the file is well created or not (a bool)
        '''
        #Acquire data
        dict_algorithm, dict_material, dict_sample, dict_sollicitation = User.All_parameters()
        dict_algorithm['i_PFDEM'] = 0
        #try to create a .i file
        Owntools.Write.Write_i(dict_algorithm, dict_material, dict_sample, dict_sollicitation)
        #Check if the file PF_CH_AC_base.i is in the directory
        self.assertTrue(Path(dict_algorithm['namefile']+'_'+str(dict_algorithm['i_PFDEM'])+'.i').is_file(),"The file namefile.i has not been created!")
        os.remove(dict_algorithm['namefile']+'_'+str(dict_algorithm['i_PFDEM'])+'.i')

    #---------------------------------------------------------------------------

    def test_index_to_str(self):
        '''
        Try to convert an integer into a string with Owntools.index_to_str().

        The string is composed by three elements.

            Output :
                The result depends on the fact if the conversions are well done or not (a bool)
        '''
        #check if the function works well in different configurations
        self.assertTrue(Owntools.index_to_str(7)=='007','The conversion index_to_str() seems to do not for 00x...')
        self.assertTrue(Owntools.index_to_str(26)=='026','The conversion index_to_str() seems to do not for 0xx...')
        self.assertTrue(Owntools.index_to_str(666)=='666','The conversion index_to_str() seems to do not for xxx...')

    #---------------------------------------------------------------------------

    def test_Plot_config(self):
        '''
        Try to plot a sample configuration with Owntools.Plot_config().

            Output :
                The result depends on the fact if the .png file is well created or not (a bool)
        '''
        #Acquire data
        dict_algorithm, dict_material, dict_sample, dict_sollicitation = User.All_parameters()
        #create the two grains
        User.Add_2grains(dict_material,dict_sample)
        #create the solute
        User.Add_solute(dict_sample)
        #create a folder
        if Path('Debug').exists():
            shutil.rmtree('Debug')
        os.mkdir('Debug')
        os.mkdir('Debug/Configuration')
        #try to plot the configuration
        Owntools.Plot.Plot_config(dict_algorithm, dict_sample)
        #check if the .png has been created
        self.assertTrue(Path('Debug/Configuration/Configuration_0.png').is_file(),"The image Debug/Configuration/Configuration_0.png has not been created by the function Owntools.Plot_config()!")
        shutil.rmtree('Debug')

    #---------------------------------------------------------------------------

    def test_Cosine_Profile(self):
        '''
        Try to compute a phase variable with Owntools.Cosine_Profile().

        Three cases are considered : inside the grain, at the interface, outside the grain.

            Output :
                The result depends on the fact if the phase variable is well computed or not (a bool)
        '''
        #check if the function works well in different configurations
        self.assertTrue(Owntools.Cosine_Profile(1,0,0.5)==1,'The Owntools.Cosine_Profile() seems to do not for a point inside the grain...')
        self.assertTrue(Owntools.Cosine_Profile(1,1,0.5)==0.5*(1 + math.cos(math.pi/2)),'The Owntools.Cosine_Profile() seems to do not for a point inside the interface...')
        self.assertTrue(Owntools.Cosine_Profile(1,2,0.5)==0,'The Owntools.Cosine_Profile() seems to do not for a point outside the grain...')

    #---------------------------------------------------------------------------

    def test_Write_eta_txt(self):
        '''
        Try to create file needed for MOOSE simulation with Owntools.Write_eta_txt().

        This file is about etai variables.
        A sample with two grains is assumed.

            Output :
                The result depends on the fact if the file is well created or not (a bool)
        '''
        #Acquire data
        dict_algorithm, dict_material, dict_sample, dict_sollicitation = User.All_parameters()
        dict_algorithm['i_PFDEM'] = 0
        #create the two grains
        User.Add_2grains(dict_material,dict_sample)
        #create a folder
        if Path('Data').exists():
            shutil.rmtree('Data')
        os.mkdir('Data')
        #try to create .txt files
        Owntools.Write.Write_eta_txt(dict_algorithm, dict_sample)
        #Check if the files are in the directory
        self.assertTrue(Path('Data/eta1_0.txt').is_file(),"The file Data/eta1_0.txt has not been created!")
        self.assertTrue(Path('Data/eta2_0.txt').is_file(),"The file Data/eta2_0.txt has not been created!")
        shutil.rmtree('Data')

    #---------------------------------------------------------------------------

    def test_Write_solute_txt(self):
        '''
        Try to create file needed for MOOSE simulation with Owntools.Write_solute_txt().

        This file is about c variable.

            Output :
                The result depends on the fact if the file is well created or not (a bool)
        '''
        #Acquire data
        dict_algorithm, dict_material, dict_sample, dict_sollicitation = User.All_parameters()
        dict_algorithm['i_PFDEM'] = 0
        #create the two grains
        User.Add_solute(dict_sample)
        #create a folder
        if Path('Data').exists():
            shutil.rmtree('Data')
        os.mkdir('Data')
        #try to create .txt files
        Owntools.Write.Write_solute_txt(dict_algorithm, dict_sample)
        #Check if the files are in the directory
        self.assertTrue(Path('Data/c_0.txt').is_file(),"The file Data/c_0.txt has not been created!")
        shutil.rmtree('Data')

    #---------------------------------------------------------------------------

    def test_Write_Emec_txt(self):
        '''
        Try to create file needed for MOOSE simulation with Owntools.Write_Emec_txt().

        This file is about the dissolution field due to mechanical loading.

            Output :
                The result depends on the fact if the file is well created or not (a bool)
        '''
        #Acquire data
        dict_algorithm, dict_material, dict_sample, dict_sollicitation = User.All_parameters()
        dict_algorithm['i_PFDEM'] = 0
        #create the two grains
        User.Add_2grains(dict_material,dict_sample)
        #create a folder
        if Path('Data').exists():
            shutil.rmtree('Data')
        os.mkdir('Data')
        #try to create .txt files
        Owntools.Write.Write_Emec_txt(dict_algorithm, dict_sample)
        #Check if the files are in the directory
        self.assertTrue(Path('Data/ep_0.txt').is_file(),"The file Data/ep_0.txt has not been created!")
        shutil.rmtree('Data')

    #---------------------------------------------------------------------------

    def test_Write_kc_txt(self):
        '''
        Try to create file needed for MOOSE simulation with Owntools.Write_kc_txt().

        This file is about the diffusion of the solute.

            Output :
                The result depends on the fact if the file is well created or not (a bool)
        '''
        #Acquire data
        dict_algorithm, dict_material, dict_sample, dict_sollicitation = User.All_parameters()
        dict_algorithm['i_PFDEM'] = 0
        #create the two grains
        User.Add_2grains(dict_material,dict_sample)
        #create a folder
        if Path('Data').exists():
            shutil.rmtree('Data')
        os.mkdir('Data')
        #Compute the diffusion coefficient
        Owntools.Compute.Compute_kc_dil(dict_material,dict_sample)
        #try to create .txt files
        Owntools.Write.Write_kc_txt(dict_algorithm, dict_sample)
        #Check if the files are in the directory
        self.assertTrue(Path('Data/kc_0.txt').is_file(),"The file Data/kc_0.txt has not been created!")
        shutil.rmtree('Data')

#-------------------------------------------------------------------------------

class TestGrain(unittest.TestCase):
    '''Test functions from Grain.py.'''
    def test_geometric_study(self):
        '''
        Try to study the grain geometry with Grain.geometric_study().

        As Monte Carlo is used, some noise can occur. It is advised to test multiple times.

            Output :
                The result depends on the fact if the center is well located or not (a bool)
                The result depends on the fact if the surface is well assumed or not (a bool)
        '''
        #Acquire data
        dict_algorithm, dict_material, dict_sample, dict_sollicitation = User.All_parameters()
        #Create one grain
        grain = Grain.Grain(0,10,np.array([np.mean(dict_sample['x_L']),np.mean(dict_sample['y_L'])]),dict_material,dict_sample)
        #try to do the geometric study of the grain
        grain.geometric_study(dict_sample)
        #define margine of the surface and the center (the study is done with a Monte Carlo method, some noise can be introduced)
        margin_center = 0.03*10 #10 is the radius see line upper
        margin_surface = 0.03*math.pi*10**2 #10 is the radius see line upper
        #check if the center computed is near the analytical one
        self.assertTrue(np.linalg.norm(np.array([np.mean(dict_sample['x_L']),np.mean(dict_sample['y_L'])])-grain.center)<margin_center,'The estimation of the center position seems false (because of the Monte Carlo Method try to rerun or increase the margin)...')
        #check if the surface computed is near the analytical one
        self.assertTrue(abs(math.pi*10**2-grain.surface)<margin_surface,'The estimation of the surface seems false (because of the Monte Carlo Method try to rerun or increase the margin)...')

    #---------------------------------------------------------------------------

    def test_P_is_inside(self):
        '''
        Try to determine if a point is inside a grain geometry with Grain.P_is_inside().

        Two cases are considered : the point is inside, the point is outside.

            Output :
                The result depends on the fact if points are well located or not (a bool)
        '''
        #Acquire data
        dict_algorithm, dict_material, dict_sample, dict_sollicitation = User.All_parameters()
        #Create one grain
        grain = Grain.Grain(0,10,np.array([np.mean(dict_sample['x_L']),np.mean(dict_sample['y_L'])]),dict_material,dict_sample)
        #check if the function works well in different configurations
        self.assertFalse(grain.P_is_inside(np.array([np.mean(dict_sample['x_L'])+11,np.mean(dict_sample['y_L'])])),'An outside point is detected as inside by Grain.P_is_inside()...')
        self.assertTrue(grain.P_is_inside(np.array([np.mean(dict_sample['x_L']),np.mean(dict_sample['y_L'])+1])),'An inside point is detected as outside by Grain.P_is_inside()...')

    #---------------------------------------------------------------------------

    def test_move_grain_rebuild(self):
        '''
        Try to move a grain by deconstruction and rebuild with Grain.move_grain_rebuild().

        As Monte Carlo is used, some noise can occur. It is advised to test multiple times.

            Output :
                The result depends on the fact if the grain is well moved or not (a bool)
        '''
        #Acquire data
        dict_algorithm, dict_material, dict_sample, dict_sollicitation = User.All_parameters()
        #Create one grain
        grain = Grain.Grain(0,10,np.array([np.mean(dict_sample['x_L']),np.mean(dict_sample['y_L'])]),dict_material,dict_sample)
        #try to move the grain
        grain.move_grain_rebuild(np.array([5,0]),dict_material,dict_sample)
        #check if the center computed is near the analytical one
        self.assertTrue(np.linalg.norm(np.array([np.mean(dict_sample['x_L'])+5,np.mean(dict_sample['y_L'])])-grain.center)==0,'The grain has not been well moved!')
        #Study the geometric of the grain
        grain.geometric_study(dict_sample)
        #define margine of the new center (the study is done with a Monte Carlo method, some noise can be introduced)
        margin_center = 0.03*10 #10 is the radius see line upper
        #check if the center computed is near the analytical one
        self.assertTrue(np.linalg.norm(np.array([np.mean(dict_sample['x_L'])+5,np.mean(dict_sample['y_L'])])-grain.center)<margin_center,'The displacement of the grain seems false after the etai_M rebuild (because of the Monte Carlo Method try to rerun or increase the margin)...')

    #---------------------------------------------------------------------------

    def test_move_grain_interpolation(self):
        '''
        Try to move a grain by interpolation with Grain.move_grain_interpolation().

        As Monte Carlo is used, some noise can occur. It is advised to test multiple times.

            Output :
                The result depends on the fact if the grain is well moved or not (a bool)
        '''
        #Acquire data
        dict_algorithm, dict_material, dict_sample, dict_sollicitation = User.All_parameters()
        #Create one grain
        grain = Grain.Grain(0,10,np.array([np.mean(dict_sample['x_L']),np.mean(dict_sample['y_L'])]),dict_material,dict_sample)
        #try to move the grain
        grain.move_grain_interpolation(np.array([5,0]),dict_sample)
        #check if the center computed is near the analytical one
        self.assertTrue(np.linalg.norm(np.array([np.mean(dict_sample['x_L'])+5,np.mean(dict_sample['y_L'])])-grain.center)==0,'The grain has not been well moved!')
        #Study the geometric of the grain
        grain.geometric_study(dict_sample)
        #define margine of the new center (the study is done with a Monte Carlo method, some noise can be introduced)
        margin_center = 0.03*10 #10 is the radius see line upper
        #check if the center computed is near the analytical one
        self.assertTrue(np.linalg.norm(np.array([np.mean(dict_sample['x_L'])+5,np.mean(dict_sample['y_L'])])-grain.center)<margin_center,'The displacement of the grain seems false after the etai_M rebuild (because of the Monte Carlo Method try to rerun or increase the margin)...')

    #---------------------------------------------------------------------------

    def test_Compute_overlap_2_grains(self):
        '''
        Try to compute the overlap between two grains with Grain.Compute_overlap_2_grains().

        Two cases are considered : 2 grains in contact, 2 grains not in contact.

            Output :
                The result depends on the fact if the overlap are well computed or not (a bool)
        '''
        #Acquire data
        dict_algorithm, dict_material, dict_sample, dict_sollicitation = User.All_parameters()
        dict_sample['grain_discretisation'] = 20
        #Create two grain
        g1 = Grain.Grain(0,10,np.array([np.mean(dict_sample['x_L'])-5,np.mean(dict_sample['y_L'])]),dict_material,dict_sample)
        g2 = Grain.Grain(0,10,np.array([np.mean(dict_sample['x_L'])+5,np.mean(dict_sample['y_L'])]),dict_material,dict_sample)
        dict_sample['L_g'] = [g1,g2]
        #Check if there is an overlap between those grains
        Grain.Compute_overlap_2_grains(dict_sample)
        self.assertTrue(dict_sample['overlap']==10,'The overlap between two grains in contact is not well computed!')
        #Create two grain
        g1 = Grain.Grain(0,10,np.array([np.mean(dict_sample['x_L'])-11,np.mean(dict_sample['y_L'])]),dict_material,dict_sample)
        g2 = Grain.Grain(0,10,np.array([np.mean(dict_sample['x_L'])+11,np.mean(dict_sample['y_L'])]),dict_material,dict_sample)
        dict_sample['L_g'] = [g1,g2]
        #Check if there is not an overlap between those grains
        Grain.Compute_overlap_2_grains(dict_sample)
        self.assertTrue(dict_sample['overlap']==-2,'The overlap between two grains not in contact is not well computed!')

    #---------------------------------------------------------------------------

    def test_Apply_overlap_target(self):
        '''
        Try to apply the overlap between two grains with Grain.Apply_overlap_target().

            Output :
                The result depends on the fact if the overlap are well applied or not (a bool)
        '''
        #Acquire data
        dict_algorithm, dict_material, dict_sample, dict_sollicitation = User.All_parameters()
        #Create two grain
        User.Add_2grains(dict_material,dict_sample)
        #Compute the initial overlap
        Grain.Compute_overlap_2_grains(dict_sample)
        #Create the dict_tracker
        dict_tracker = {'L_displacement': [0], 'L_int_displacement' : [0]}
        #try to apply a target overlap
        Grain.Apply_overlap_target(dict_material,dict_sample,dict_sollicitation,dict_tracker)
        #Compute the current overlap
        Grain.Compute_overlap_2_grains(dict_sample)
        #Check if the overlap target is well applied
        self.assertTrue(dict_sample['overlap']==dict_sollicitation['overlap_target'],'The target overlap between two grains is not well applied!')

#-------------------------------------------------------------------------------
#main
#-------------------------------------------------------------------------------

unittest.main(verbosity = 2)
