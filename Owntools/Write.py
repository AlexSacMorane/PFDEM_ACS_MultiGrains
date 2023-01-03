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
    elif j == 20:
      line = ''
      for etai in dict_sample['L_etai']:
          line = line + '\t[./eta'+str(int(etai.id+1))+']\n'+\
                        '\t\torder = FIRST\n'+\
                        '\t\tfamily = LAGRANGE\n'+\
                        '\t\t[./InitialCondition]\n'+\
                        '\t\t\ttype = FunctionIC\n'+\
                        '\t\t\tfunction = eta'+str(int(etai.id+1))+'_txt\n'+\
                        '\t\t[../]\n'+\
                        '\t[../]\n'
    elif j == 30:
      line = ''
      for etai in dict_sample['L_etai']:
          line = line + '\t#\n'+\
                        '\t# Order parameter eta'+str(int(etai.id+1))+'\n'+\
                        '\t#\n'+\
                        '\t[./deta'+str(int(etai.id+1))+'dt]\n'+\
                        '\t\ttype = TimeDerivative\n'+\
                        '\t\tvariable = eta'+str(int(etai.id+1))+'\n'+\
                        '\t[../]\n'+\
                        '\t[./ACBulk'+str(int(etai.id+1))+']\n'+\
                        '\t\ttype = AllenCahn\n'+\
                        '\t\tvariable = eta'+str(int(etai.id+1))+'\n'
          args_str = ''
          for etaj in dict_sample['L_etai'] :
              if etai.id != etaj.id :
                  args_str = args_str + 'eta'+str(int(etaj.id+1))+' '
          line = line + "\t\targs = '"+args_str[:-1]+"'\n"+\
                        '\t\tmob_name = L\n'+\
                        '\t\tf_name = F_total\n'+\
                        '\t[../]\n'+\
                        '\t[./ACInterface'+str(int(etai.id+1))+']\n'+\
                        '\t\ttype = ACInterface\n'+\
                        '\t\tvariable = eta'+str(int(etai.id+1))+'\n'+\
                        '\t\tmob_name = L\n'+\
                        "\t\tkappa_name = 'kappa_eta'\n"+\
                        '\t[../]\n'
    elif j == 38:
      line = ''
      for etai in dict_sample['L_etai']:
          line = line + '\t[./eta'+str(int(etai.id+1))+'_c]\n'+\
                        '\t\ttype = CoupledTimeDerivative\n'+\
                        "\t\tv = 'eta"+str(int(etai.id+1))+"'\n"+\
                        '\t\tvariable = c\n'+\
                        '\t[../]\n'
    elif j == 50:
      line = "\t\tprop_values ='"+str(dict_material['M'])+' '+str(dict_material['kappa_eta'])+"'\n"
    elif j == 68:
        args_str = ''
        for etai in dict_sample['L_etai'] :
            args_str = args_str + 'eta'+str(int(etai.id+1))+' '
        line = "\t\targs = '"+args_str[:-1]+"'\n"
    elif j == 70:
        line = "\t\tconstant_expressions = '"+str(dict_material['Energy_barrier'])+"'\n"
    elif j == 71:
        line = "\t\tfunction = 'h*"
        etai_str = ''
        for etai in dict_sample['L_etai']:
            etai_str = etai_str + 'eta'+str(int(etai.id+1))+'^2*(1-eta'+str(int(etai.id+1))+')^2+'
        line = line + etai_str[:-1] + ")'\n"
    elif j == 79:
        line = "\t\targs = '"
        etai_str = ''
        for etai in dict_sample['L_etai']:
            etai_str = etai_str + 'eta'+str(int(etai.id+1))+' '
        line = line + etai_str[:-1] + "'\n"
    elif j == 81:
        line =  "\t\tfunction = '"
        etai_str = ''
        for etai in dict_sample['L_etai']:
            etai_str = etai_str + 'ep*3*eta'+str(int(etai.id+1))+'^2-ep*2*eta'+str(int(etai.id+1))+'^3+'
        line = line + etai_str[:-1] + "'\n"
    elif j == 89:
        line =  "\t\targs = 'c "
        etai_str = ''
        for etai in dict_sample['L_etai']:
            etai_str = etai_str + 'eta'+str(int(etai.id+1))+' '
        line = line + etai_str[:-1] + "'\n"
    elif j == 91 :
      line = "\t\tconstant_expressions = '"+str(dict_sollicitation['chi'])+"'\n"
    elif j == 92:
        line =  "\t\tfunction = '"
        etai_str = ''
        for etai in dict_sample['L_etai']:
            etai_str = etai_str + 'chi*c*(3*eta'+str(int(etai.id+1))+'^2-2*eta'+str(int(etai.id+1))+'^3)+'
        line = line + etai_str[:-1] + "'\n"
    elif j == 100:
        line =  "\t\targs = '"
        etai_str = ''
        for etai in dict_sample['L_etai']:
            etai_str = etai_str + 'eta'+str(int(etai.id+1))+' '
        line = line + etai_str[:-1] + "'\n"
    elif j == 101:
        etai_str = ''
        for etai in dict_sample['L_etai']:
            etai_str = etai_str + 'eta'+str(int(etai.id+1))+','
        line =  "\t\tmaterial_property_names = 'F("+etai_str[:-1]+') Ed_mec('+etai_str[:-1]+') Ed_pre('+etai_str[:-1]+")'\n"
    elif j == 109:
        line = ''
        for etai in dict_sample['L_etai']:
            line = line + '\t[eta'+str(int(etai.id+1))+'_txt]\n'+\
                          '\t\ttype = PiecewiseMultilinear\n'+\
                          '\t\tdata_file = Data/eta'+str(int(etai.id+1))+'_'+str(dict_algorithm['i_PFDEM'])+'.txt\n'+\
                          '\t[]\n'
    elif j == 112:
        line = '\t\tdata_file = Data/c_'+str(dict_algorithm['i_PFDEM'])+'.txt\n'
    elif j == 116:
        line = '\t\tdata_file = Data/ep_'+str(dict_algorithm['i_PFDEM'])+'.txt\n'
	elif j == 120:
        line = '\t\tdata_file = Data/kc_'+str(dict_algorithm['i_PFDEM'])+'.txt\n'
    elif j == 150:
        line =  '\t\tend_time = '+str(dict_algorithm['dt_PF']*dict_algorithm['n_t_PF'])+'\n'
    elif j == 154:
        line =  '\t\tdt = '+str(dict_algorithm['dt_PF']) +'\n'

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
