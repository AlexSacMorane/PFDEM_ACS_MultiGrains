# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This file contains the different functions used in the simulation.
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

from pathlib import Path
import os
import math

#-------------------------------------------------------------------------------

def index_to_str(j):
  '''
  An integer is converted to a float with 3 components

    Input :
        a number (a float)
    Output :
        a string with 3 components (a string)
  '''
  if j < 10:
      j_str = '00'+str(j)
  elif 10 <= j and j < 100:
      j_str = '0'+str(j)
  else :
      j_str = str(j)
  return j_str

#-------------------------------------------------------------------------------

def Sort_Files(dict_algorithm):
     '''
     Sort files generated by MOOSE to different directories

        Input :
            an algorithm dictionnary (a dict)
        Output :
            Nothing but files are sorted
     '''

     os.rename(dict_algorithm['namefile']+'_'+str(dict_algorithm['i_PFDEM'])+'_out.e','Output/'+dict_algorithm['namefile']+'_'+str(dict_algorithm['i_PFDEM'])+'_out.e')
     os.rename(dict_algorithm['namefile']+'_'+str(dict_algorithm['i_PFDEM'])+'.i','Input/'+dict_algorithm['namefile']+'_'+str(dict_algorithm['i_PFDEM'])+'.i')
     j = 0
     j_str = index_to_str(j)
     filepath = Path(dict_algorithm['namefile']+'_'+str(dict_algorithm['i_PFDEM'])+'_other_'+j_str+'.pvtu')
     while filepath.exists():
         for i_proc in range(dict_algorithm['np_proc']):
            os.rename(dict_algorithm['namefile']+'_'+str(dict_algorithm['i_PFDEM'])+'_other_'+j_str+'_'+str(i_proc)+'.vtu','Output/Ite_'+str(dict_algorithm['i_PFDEM'])+'/'+dict_algorithm['namefile']+'_'+str(dict_algorithm['i_PFDEM'])+'_other_'+j_str+'_'+str(i_proc)+'.vtu')
         os.rename(dict_algorithm['namefile']+'_'+str(dict_algorithm['i_PFDEM'])+'_other_'+j_str+'.pvtu','Output/Ite_'+str(dict_algorithm['i_PFDEM'])+'/'+dict_algorithm['namefile']+'_'+str(dict_algorithm['i_PFDEM'])+'_other_'+j_str+'.pvtu')
         j = j + 1
         j_str = index_to_str(j)
         filepath = Path(dict_algorithm['namefile']+'_'+str(dict_algorithm['i_PFDEM'])+'_other_'+j_str+'.pvtu')

     return index_to_str(j-1)

#-------------------------------------------------------------------------------

def Extract_solute_at_p(dict_sample,ij_p):
    '''
    Extract the value of the solute concentration at a given point.

     Input :
         a sample dictionnary (a dict)
         a coordinate of the point (a tuple of int)
     Output :
         the value of the solute concentration (a float)
    '''
    return dict_sample['solute_M'][-1-ij_p[0]][ij_p[1]]

#-------------------------------------------------------------------------------

def Cosine_Profile(R,r,w):
    '''
    Compute the phase field variable at some point.

    A cosine profile is assumed (see https://mooseframework.inl.gov/source/ics/SmoothCircleIC.html).

    Input :
        the radius R of the grain in the direction (a float)
        the distance r between the current point and the center (a float)
        the width w of the interface (a float)
    Output :
        the value of the phase field variable (a float)
    '''
    #inside the grain
    if r<R-w/2:
        return 1
    #outside the grain
    elif r>R+w/2:
        return 0
    #inside the interface
    else :
        return 0.5*(1 + math.cos(math.pi*(r-R+w/2)/w))

#-------------------------------------------------------------------------------

def error_on_ymax_f(dy,overlap_L,k_L,Force_target) :
    """
    Compute the function f to control the upper wall. It is the difference between the force applied and the target value.

        Input :
            an increment of the upper wall position (a float)
            a list of overlap for contact between grain and upper wall (a list)
            a list of spring for contact between grain and upper wall (a list)
            a confinement force (a float)
        Output :
            the difference between the force applied and the confinement (a float)
    """
    f = Force_target
    for i in range(len(overlap_L)):
        f = f - k_L[i]*(max(overlap_L[i]-dy,0))**(3/2)
    return f

#-------------------------------------------------------------------------------

def error_on_ymax_df(dy,overlap_L,k_L) :
    """
    Compute the derivative function df to control the upper wall (error_on_ymax_f()).

        Input :
            an increment of the upper wall position (a float)
            a list of overlap for contact between grain and upper wall (a list)
            a list of spring for contact between grain and upper wall (a list)
        Output :
            the derivative of error_on_ymax_f() (a float)
    """
    df = 0
    for i in range(len(overlap_L)):
        df = df + 3/2*k_L[i]*(max(overlap_L[i]-dy,0))**(1/2)
    return df

#-------------------------------------------------------------------------------

def Control_y_max_NR(dict_sample,dict_sollicitation):
    """
    Control the upper wall to apply force.

    A Newton-Raphson method is applied to verify the confinement.
        Input :
            a sample dictionnary (a dict)
            a sollicitations dictionnary (a dict)
        Output :
            Nothing, but the sample dictionnary is updated, concerning the upper wall position and force applied (two floats)
    """
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    #load data needed
    Force_target = dict_sollicitation['Vertical_Confinement_Force']
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

    F = 0
    overlap_L = []
    k_L = []
    for contact in dict_sample['L_contact_gw']:
        if contact.nature == 'gwy_max':
            F = F + contact.Fwg_n
            overlap_L.append(contact.overlap)
            k_L.append(contact.k)
            #compute force applied, save contact overlap and spring

    if overlap_L != []:
        i_NR = 0
        dy = 0
        ite_criteria = True
        if -0.01*Force_target<error_on_ymax_f(dy,overlap_L,k_L,Force_target) and error_on_ymax_f(dy,overlap_L,k_L,Force_target)<0.01*Force_target:
            ite_criteria = False
        while ite_criteria :
            i_NR = i_NR + 1
            dy = dy - error_on_ymax_f(dy,overlap_L,k_L,Force_target)/error_on_ymax_df(dy,overlap_L,k_L)
            if i_NR > 100:
                ite_criteria = False
            if -0.01*Force_target<error_on_ymax_f(dy,overlap_L,k_L,Force_target) and error_on_ymax_f(dy,overlap_L,k_L,Force_target)<0.01*Force_target:
                ite_criteria = False
        dict_sample['y_box_max'] = dict_sample['y_box_max'] + dy

    else :
        #if there is no contact with the upper wall, the wall is reset
        dict_sample['y_box_max'] = Reset_y_max(dict_sample['L_g'],Force_target)

    for contact in dict_sample['L_contact_gw']:
        if contact.nature == 'gwy_max':
            #reactualisation
            contact.limit = dict_sample['y_box_max']

    #Update dict
    dict_sollicitation['Force_on_upper_wall'] = F

#-------------------------------------------------------------------------------

def Reset_y_max(L_g,Force):
    """
    The upper wall is located as a single contact verify the target value.

        Input :
            the list of temporary grains (a list)
            the confinement force (a float)
        Output :
            the upper wall position (a float)
    """
    print('Reset of the y_max on the upper grain')
    y_max = None
    id_grain_max = None
    for id_grain in range(len(L_g)):
        grain = L_g[id_grain]
        y_max_grain = max(grain.l_border_y)

        if y_max != None and y_max_grain > y_max:
            y_max = y_max_grain
            id_grain_max = id_grain
        elif y_max == None:
            y_max = y_max_grain
            id_grain_max = id_grain

    k = 5*4/3*L_g[id_grain_max].y/(1-L_g[id_grain_max].nu*L_g[id_grain_max].nu)*math.sqrt(L_g[id_grain_max].r_mean)
    y_max = y_max - (Force/k)**(2/3)

    return y_max
