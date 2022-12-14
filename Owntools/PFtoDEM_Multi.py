# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This file contains the different functions used to translate phase field data in DEM data in the simulation.
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

import numpy as np

#-------------------------------------------------------------------------------

def solute_PFtoDEM_Multi(FileToRead,dict_algorithm,dict_sample):
    '''
    Read file from MOOSE simulation to reconstruct the phase field of the solute.

        Input :
            the name of the file to read (a string)
            an algorithm dictionnary (a dictionnary)
            a sample dictionnary (a dictionnary)
        Output :
            Nothing but the sample dictionnary gets an updated attribute (a n_y x n_x numpy array)
    '''
    #---------------------------------------------------------------------------
    #Global parameters
    #---------------------------------------------------------------------------

    dict_sample['solute_M'] = np.array(np.zeros((len(dict_sample['y_L']),len(dict_sample['x_L']))))

    id_L = None
    c_selector_len = len('        <DataArray type="Float64" Name="c')
    end_len = len('        </DataArray>')
    XYZ_selector_len = len('        <DataArray type="Float64" Name="Points"')
    data_jump_len = len('          ')

    for i_proc in range(dict_algorithm['np_proc']):

        L_Work = [[], #X
                  [], #Y
                  []] #c

    #---------------------------------------------------------------------------
    #Reading file
    #---------------------------------------------------------------------------

        f = open(f'{FileToRead}_{i_proc}.vtu','r')
        data = f.read()
        f.close()
        lines = data.splitlines()

        #iterations on line
        for line in lines:

            if line[0:c_selector_len] == '        <DataArray type="Float64" Name="c':
                id_L = 2

            elif line[0:XYZ_selector_len] == '        <DataArray type="Float64" Name="Points"':
                id_L = 0

            elif (line[0:end_len] == '        </DataArray>' or  line[0:len('          <InformationKey')] == '          <InformationKey') and id_L != None:
                id_L = None

            elif line[0:data_jump_len] == '          ' and id_L == 2: #Read c
                line = line[data_jump_len:]
                c_start = 0
                for c_i in range(0,len(line)):
                    if line[c_i]==' ':
                        c_end = c_i
                        L_Work[id_L].append(float(line[c_start:c_end]))
                        c_start = c_i+1
                L_Work[id_L].append(float(line[c_start:]))

            elif line[0:data_jump_len] == '          ' and id_L == 0: #Read [X, Y, Z]
                line = line[data_jump_len:]
                XYZ_temp = []
                c_start = 0
                for c_i in range(0,len(line)):
                    if line[c_i]==' ':
                        c_end = c_i
                        XYZ_temp.append(float(line[c_start:c_end]))
                        if len(XYZ_temp)==3:
                            L_Work[0].append(XYZ_temp[0])
                            L_Work[1].append(XYZ_temp[1])
                            XYZ_temp = []
                        c_start = c_i+1
                XYZ_temp.append(float(line[c_start:]))
                L_Work[0].append(XYZ_temp[0])
                L_Work[1].append(XYZ_temp[1])

        #Adaptating data and update of solute_M
        for i in range(len(L_Work[0])):
            #Interpolation method
            L_dy = []
            for y_i in dict_sample['y_L'] :
                L_dy.append(abs(y_i - L_Work[1][i]))
            L_dx = []
            for x_i in dict_sample['x_L'] :
                L_dx.append(abs(x_i - L_Work[0][i]))
            dict_sample['solute_M'][-1-list(L_dy).index(min(L_dy))][list(L_dx).index(min(L_dx))] = L_Work[2][i]

#---------------------------------------------------------------------------

def Ed_PFtoDEM_Multi(FileToRead,dict_algorithm,dict_sample):
    '''
    Read file from MOOSE simulation to follow the external energy used in the phase field formulation.

        Input :
            the name of the file to read (a string)
            an algorithm dictionnary (a dictionnary)
            a sample dictionnary (a dictionnary)
        Output :
            Nothing but the sample dictionnary gets updated attributes (three n_y x n_x numpy array)
    '''
    #---------------------------------------------------------------------------
    #Global parameters
    #---------------------------------------------------------------------------

    id_L = None
    Emec_selector_len = len('        <DataArray type="Float64" Name="Ed_mec')
    Eche_selector_len = len('        <DataArray type="Float64" Name="Ed_pre')
    end_len = len('        </DataArray>')
    XYZ_selector_len = len('        <DataArray type="Float64" Name="Points"')
    data_jump_len = len('          ')

    for i_proc in range(dict_algorithm['np_proc']):

        L_Work = [[], #X
                  [], #Y
                  [], #Emec
                  []] #Eche

    #---------------------------------------------------------------------------
    #Reading file
    #---------------------------------------------------------------------------

        f = open(f'{FileToRead}_{i_proc}.vtu','r')
        data = f.read()
        f.close()
        lines = data.splitlines()

        #iterations on line
        for line in lines:

            if line[0:Emec_selector_len] == '        <DataArray type="Float64" Name="Ed_mec':
                id_L = 2

            elif line[0:Eche_selector_len] == '        <DataArray type="Float64" Name="Ed_pre':
                id_L = 3

            elif line[0:XYZ_selector_len] == '        <DataArray type="Float64" Name="Points"':
                id_L = 0

            elif (line[0:end_len] == '        </DataArray>' or  line[0:len('          <InformationKey')] == '          <InformationKey') and id_L != None:
                id_L = None

            elif line[0:data_jump_len] == '          ' and (id_L == 2 or id_L == 3): #Read Ed_mec or Ed_pre
                line = line[data_jump_len:]
                c_start = 0
                for c_i in range(0,len(line)):
                    if line[c_i]==' ':
                        c_end = c_i
                        L_Work[id_L].append(float(line[c_start:c_end]))
                        c_start = c_i+1
                L_Work[id_L].append(float(line[c_start:]))

            elif line[0:data_jump_len] == '          ' and id_L == 0: #Read [X, Y, Z]
                line = line[data_jump_len:]
                XYZ_temp = []
                c_start = 0
                for c_i in range(0,len(line)):
                    if line[c_i]==' ':
                        c_end = c_i
                        XYZ_temp.append(float(line[c_start:c_end]))
                        if len(XYZ_temp)==3:
                            L_Work[0].append(XYZ_temp[0])
                            L_Work[1].append(XYZ_temp[1])
                            XYZ_temp = []
                        c_start = c_i+1
                XYZ_temp.append(float(line[c_start:]))
                L_Work[0].append(XYZ_temp[0])
                L_Work[1].append(XYZ_temp[1])

        #Adaptating data and update of Ed_M, Emec_M and Eche_M
        for i in range(len(L_Work[0])):
            #Interpolation method
            L_dy = []
            for y_i in dict_sample['y_L'] :
                L_dy.append(abs(y_i - L_Work[1][i]))
            L_dx = []
            for x_i in dict_sample['x_L'] :
                L_dx.append(abs(x_i - L_Work[0][i]))
            dict_sample['Emec_M'][-1-list(L_dy).index(min(L_dy))][list(L_dx).index(min(L_dx))] = L_Work[2][i]
            dict_sample['Eche_M'][-1-list(L_dy).index(min(L_dy))][list(L_dx).index(min(L_dx))] = L_Work[3][i]
            dict_sample['Ed_M'][-1-list(L_dy).index(min(L_dy))][list(L_dx).index(min(L_dx))] = L_Work[2][i] - L_Work[3][i]
