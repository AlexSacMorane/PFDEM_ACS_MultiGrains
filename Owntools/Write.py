# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This file contains the different functions to write files used in the simulation.
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------

def Write_i(dict_algorithm, dict_material, dict_sample, dict_sollicitation):
  '''
  Create the .i file to run MOOSE simulation.

  The file is generated from a template nammed PF_ACS_base.i

    Input :
        a algorithm dictionnary (a dictionnary)
        a sample dictionnary (a dictionnary)
        a material dictionnary (a dictionnary)
    Output :
        Nothing but a txt .i file is generated (a file)
  '''
  file_to_write = open(dict_algorithm['namefile']+'_'+str(dict_algorithm['i_PFDEM'])+'.i','w')
  file_to_read = open('PF_ACS_base.i','r')
  lines = file_to_read.readlines()
  file_to_read.close()

  j = 0
  for line in lines :
    j = j + 1
    if j == 4:
      line = line[:-1] + ' ' + str(len(dict_sample['x_L'])-1)+'\n'
    elif j == 5:
      line = line[:-1] + ' ' + str(len(dict_sample['y_L'])-1)+'\n'
    elif j == 7:
      line = line[:-1] + ' ' + str(min(dict_sample['x_L']))+'\n'
    elif j == 8:
      line = line[:-1] + ' ' + str(max(dict_sample['x_L']))+'\n'
    elif j == 9:
      line = line[:-1] + ' ' + str(min(dict_sample['y_L']))+'\n'
    elif j == 10:
      line = line[:-1] + ' ' + str(max(dict_sample['y_L']))+'\n'
    elif j == 116:
      line = line[:-1] + "'"+str(dict_material['M'])+' '+str(dict_material['kappa_eta'])+"'\n"
    elif j == 136:
      line = line[:-1] + ' ' + str(dict_material['Energy_barrier'])+"'\n"
    elif j == 157 :
      line = line[:-1] + str(dict_sollicitation['chi'])+"'\n"
    elif j == 177 or j == 181 or j == 185 or j == 189 or j == 193:
      line = line[:-1] + str(dict_algorithm['i_PFDEM']) + '.txt\n'
    elif j == 223:
      line = line[:-1] + ' ' + str(dict_algorithm['dt_PF']*dict_algorithm['n_t_PF']) +'\n'
    elif j == 227:
      line = line[:-1] + ' ' + str(dict_algorithm['dt_PF']) +'\n'
    file_to_write.write(line)

  file_to_write.close()

#-------------------------------------------------------------------------------

def Write_eta_txt(dict_algorithm, dict_sample):
    '''
    Write a .txt file needed for MOOSE simulation.

    The variables etai are transmitted to the MOOSE simulation.

        Input :
            an algorithm dictionnary (a dict)
            an sample dictionnary (a dict)
        Output :
            Nothing but a .txt file is generated (a file)
    '''
    for etai in dict_sample['L_etai'] :
        file_to_write = open('Data/eta'+str(int(etai.id+1))+'_'+str(dict_algorithm['i_PFDEM'])+'.txt','w')
        file_to_write.write('AXIS X\n')
        line = ''
        for x in dict_sample['x_L']:
            line = line + str(x)+ ' '
        line = line + '\n'
        file_to_write.write(line)

        file_to_write.write('AXIS Y\n')
        line = ''
        for y in dict_sample['y_L']:
            line = line + str(y)+ ' '
        line = line + '\n'
        file_to_write.write(line)

        file_to_write.write('DATA\n')
        for l in range(len(dict_sample['y_L'])):
            for c in range(len(dict_sample['x_L'])):
                file_to_write.write(str(etai.etai_M[-1-l][c])+'\n')

        file_to_write.close()

#-------------------------------------------------------------------------------

def Write_solute_txt(dict_algorithm, dict_sample):
    '''
    Write a .txt file needed for MOOSE simulation.

    The variable c is transmitted to the MOOSE simulation.

        Input :
            an algorithm dictionnary (a dict)
            an sample dictionnary (a dict)
        Output :
            Nothing but a .txt file is generated (a file)
    '''
    file_to_write = open('Data/c_'+str(dict_algorithm['i_PFDEM'])+'.txt','w')
    file_to_write.write('AXIS X\n')
    line = ''
    for x in dict_sample['x_L']:
        line = line + str(x)+ ' '
    line = line + '\n'
    file_to_write.write(line)

    file_to_write.write('AXIS Y\n')
    line = ''
    for y in dict_sample['y_L']:
        line = line + str(y)+ ' '
    line = line + '\n'
    file_to_write.write(line)

    file_to_write.write('DATA\n')
    for l in range(len(dict_sample['y_L'])):
        for c in range(len(dict_sample['x_L'])):
            file_to_write.write(str(dict_sample['solute_M'][-1-l][c])+'\n')

    file_to_write.close()

#-------------------------------------------------------------------------------

def Write_Emec_txt(dict_algorithm, dict_sample):
    '''
    Write a .txt file needed for MOOSE simulation.

    The variable ep is transmitted to the MOOSE simulation.
    This variable is the dissolution field due to mechanical energy at the contact level.

        Input :
            an algorithm dictionnary (a dict)
            a sample dictionnary (a dict)
            a sollicitation dictionnary (a dict)
        Output :
            Nothing but a .txt file is generated (a file)
    '''
    #write data
    file_to_write = open('Data/ep_'+str(dict_algorithm['i_PFDEM'])+'.txt','w')
    file_to_write.write('AXIS X\n')
    line = ''
    for x in dict_sample['x_L']:
        line = line + str(x)+ ' '
    line = line + '\n'
    file_to_write.write(line)

    file_to_write.write('AXIS Y\n')
    line = ''
    for y in dict_sample['y_L']:
        line = line + str(y)+ ' '
    line = line + '\n'
    file_to_write.write(line)

    file_to_write.write('DATA\n')
    for l in range(len(dict_sample['y_L'])):
        for c in range(len(dict_sample['x_L'])):
            file_to_write.write(str(dict_sample['Emec_M'][-1-l][c])+'\n')

    file_to_write.close()

#-------------------------------------------------------------------------------

def Write_kc_txt(dict_algorithm, dict_sample):
    '''
    Write a .txt file needed for MOOSE simulation.

    The variable kc is transmitted to the MOOSE simulation.
    This variable is the diffusion coefficient of the solute.
    It takes the value 0 if the point is inside one grain and not in the other.
    Else it takes an user defined value.

        Input :
            an algorithm dictionnary (a dict)
            an material dictionnary (a dict)
            an sample dictionnary (a dict)
        Output :
            Nothing but a .txt file is generated (a file)
    '''

    file_to_write = open('Data/kc_'+str(dict_algorithm['i_PFDEM'])+'.txt','w')
    file_to_write.write('AXIS X\n')
    line = ''
    for x in dict_sample['x_L']:
        line = line + str(x)+ ' '
    line = line + '\n'
    file_to_write.write(line)

    file_to_write.write('AXIS Y\n')
    line = ''
    for y in dict_sample['y_L']:
        line = line + str(y)+ ' '
    line = line + '\n'
    file_to_write.write(line)

    file_to_write.write('DATA\n')
    for l in range(len(dict_sample['y_L'])):
        for c in range(len(dict_sample['x_L'])):
            file_to_write.write(str(dict_sample['kc_M'][-l-1][c])+'\n')

    file_to_write.close()
