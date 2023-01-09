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
        line = "\t\tfunction = 'h*("
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

#-------------------------------------------------------------------------------

def Write_DEM_txt(dict_algorithm,dict_sample):
    """
    Write a .txt file to give information about grains and contacts.

        Input :
            an algorithm dictionnary (a dict)
            a sample dictionnary (a dict)
        Output :
            Nothing, but a .txt file is generated (a .txt)
    """
    file_to_write = open('Debug/txt/PF_'+str(dict_algorithm['i_PFDEM'])+'_ite_'+str(dict_algorithm['i_DEM'])+'.txt','w')
    file_to_write.write('Iteration PF : '+str(dict_algorithm['i_PFDEM'])+'\n'+\
                        'Iteration DEM : '+str(dict_algorithm['i_DEM'])+'\n')
    file_to_write.write('\n')
    file_to_write.write('GRAINS LIST\n')
    file_to_write.write('\n')
    for grain in dict_sample['L_g']:
        file_to_write.write('<grain_o>\n')
        file_to_write.write('\tid : '+str(grain.id)+'\n')
        file_to_write.write('\tCenter : ['+str(int(grain.center[0]))+', '+str(int(grain.center[1]))+']\n')
        file_to_write.write('\tSpeed : ['+str(round(grain.v[0],2))+', '+str(round(grain.v[1],2))+']\n')
        file_to_write.write('\tForce applied : ['+str(round(grain.f[0],1))+', '+str(round(grain.f[1],1))+']\n')
        file_to_write.write('\tOmega : '+str(grain.w)+'\n')
        file_to_write.write('\tSurface : '+str(int(grain.surface))+'\n')
        file_to_write.write('\tCoordinate X of the border : '+str(grain.l_border_x)+'\n')
        file_to_write.write('\tCoordinate Y of the border : '+str(grain.l_border_y)+'\n')
        file_to_write.write('<grain_c>\n')
    file_to_write.write('\n')
    file_to_write.write('CONTACTS LIST\n')
    file_to_write.write('\n')
    for contact in dict_sample['L_contact']:
        file_to_write.write('<contact_o>\n')
        file_to_write.write('\tid : '+str(contact.id)+'\n')
        file_to_write.write('\tGrains : '+str(contact.g1.id)+'-'+str(contact.g2.id)+'\n')
        file_to_write.write('\tNormal : ['+str(round(contact.pc_normal[0],2))+', '+str(round(contact.pc_normal[1],2))+']\n')
        file_to_write.write('\tNormal overlap : '+str(round(contact.overlap_normal,2))+'\n')
        file_to_write.write('\tNormal reaction : '+str(int(contact.F_2_1_n))+'\n')
        file_to_write.write('\tNormal damping : '+str(int(contact.F_2_1_damp))+'\n')
        file_to_write.write('\tTangential : ['+str(round(contact.pc_tangential[0],2))+', '+str(round(contact.pc_tangential[1],2))+']\n')
        file_to_write.write('\tTangential overlap : '+str(round(contact.overlap_tangential,2))+'\n')
        file_to_write.write('\tTangential reaction : '+str(int(contact.ft))+'\n')
        file_to_write.write('\tTangential damping : '+str(int(contact.ft_damp))+'\n')
        file_to_write.write('<contact_c>\n')
    file_to_write.write('\n')
    file_to_write.write('CONTACTS WITH WALL LIST\n')
    file_to_write.write('\n')
    for contact in dict_sample['L_contact_gw']:
        file_to_write.write('<contact_w_o>\n')
        file_to_write.write('\tid : '+str(contact.id)+'\n')
        file_to_write.write('\tType : '+str(contact.nature)+'\n')
        file_to_write.write('\tGrain : '+str(contact.g.id)+'\n')
        file_to_write.write('\tLimit : '+str(contact.limit)+'\n')
        file_to_write.write('\tNormal : ['+str(round(contact.nwg[0],2))+', '+str(round(contact.nwg[1],2))+']\n')
        file_to_write.write('\tNormal overlap : '+str(round(contact.overlap,2))+'\n')
        file_to_write.write('\tNormal reaction : '+str(int(contact.Fwg_n))+'\n')
        file_to_write.write('\tNormal damping : '+str(int(contact.Fwg_damp_n))+'\n')
        file_to_write.write('\tTangential : ['+str(round(contact.twg[0],2))+', '+str(round(contact.twg[1],2))+']\n')
        file_to_write.write('\tTangential overlap : '+str(round(contact.overlap_tangential,2))+'\n')
        file_to_write.write('\tTangential reaction : '+str(int(contact.ft))+'\n')
        file_to_write.write('<contact_w_c>\n')
    file_to_write.close()
