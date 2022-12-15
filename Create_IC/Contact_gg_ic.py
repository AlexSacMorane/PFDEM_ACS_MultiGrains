# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This file contains ??.nt functions used in the simulation.
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

import numpy as np

#Own
import Grain_ic

#-------------------------------------------------------------------------------
#Class
#-------------------------------------------------------------------------------

class Contact_Tempo:
  """
  A temporary contact grain - grain used to generated an initial condition.
  """

#-------------------------------------------------------------------------------

  def __init__(self, ID, G1, G2, dict_material):
    """
    Defining the contact grain-grain.

        Input :
            itself (a contact_tempo)
            an id (a int)
            two grains (two grain_tempos)
            a material dictionnary (a dict)
        Output :
            Nothing, but the contact grain - grain is generated (a contact_tempo)
    """
    self.id = ID
    self.g1 = G1
    self.g2 = G2
    self.ft = 0
    self.mu = 0
    self.coeff_restitution = dict_material['coeff_restitution']
    self.tangential_old_statut = False
    self.overlap_tangential = 0

#-------------------------------------------------------------------------------

  def normal(self):
    """
    Compute the normal reaction of a contact grain-grain.

    Here a pontual spring is considered.

        Input :
            itself (a contact_tempo)
        Output :
            Nothing, but attributes are updated
    """
    # Compute the normal and tangential planes
    PC_normal = (self.g2.center-self.g1.center)/np.linalg.norm(self.g2.center-self.g1.center)
    self.pc_normal = PC_normal #n12
    self.pc_tangential = np.array([-PC_normal[1],PC_normal[0]])

    # Compute the overlap
    overlap = self.g1.radius + self.g2.radius - np.linalg.norm(self.g1.center - self.g2.center)
    self.overlap_normal = overlap

    if overlap > 0:
    #-----------------------------------------------------------------------------
    # Compute the reaction
    #-----------------------------------------------------------------------------

        #Spring term
        Y_eq = 1/((1-self.g1.nu*self.g1.nu)/self.g1.y+(1-self.g2.nu*self.g2.nu)/self.g2.y)
        R_eq = 1/(1/np.linalg.norm(self.g1.center-self.g1.l_border[self.ij_min[0]])+1/np.linalg.norm(self.g2.center-self.g2.l_border[self.ij_min[1]]))
        k = 4/3*Y_eq*math.sqrt(R_eq)
        F_2_1_n = -k * overlap**(3/2)  #unlinear spring
        F_2_1 = F_2_1_n * PC_normal
        self.F_2_1_n = F_2_1_n
        self.Ep_n = 2/5 * k * overlap**(5/2) #-dEp/dx = F_2_1_n
        self.g1.add_F( F_2_1, self.g1.center + self.g1.radius*self.pc_normal)
        self.g2.add_F(-F_2_1, self.g2.center - self.g2.radius*self.pc_normal)

        #Damping term
        gamma = -math.log(self.coeff_restitution)/math.sqrt(math.pi**2+math.log(self.coeff_restitution)**2)
        mass_eq = self.g1.mass*self.g2.mass/(self.g1.mass+self.g2.mass)
        eta = 2 * gamma * math.sqrt(mass_eq*k)
        F_2_1_damp_n = np.dot(self.g2.v - self.g1.v,PC_normal)*eta
        F_2_1_damp = F_2_1_damp_n *PC_normal
        self.F_2_1_damp = F_2_1_damp_n
        self.g1.add_F( F_2_1_damp, self.g1.center + self.g1.radius*self.pc_normal)
        self.g2.add_F(-F_2_1_damp, self.g2.center - self.g2.radius*self.pc_normal)

    #no contact finally
    else :
        self.F_2_1_n = 0
        self.F_2_1_damp = 0
        self.Ep_n = 0

#-------------------------------------------------------------------------------

  def tangential(self,dt_DEM):
    """
    Compute the tangential reaction of a contact grain-grain.

    Here a pontual spring is considered

        Input :
            itself (a contact_tempo)
            a time step (a float)
        Output :
            Nothing, but attributes are updated
    """
    if self.overlap_normal > 0 and self.mu > 0:

        if self.tangential_old_statut:
          #if a reaction has been already computed
          #need to project the tangential reaction on the new tangential plane
          self.ft = self.ft*np.dot(self.tangential_old,self.pc_tangential)
        else:
          self.tangential_old_statut = True

        G_eq = 1/((1-self.g1.nu)/self.g1.g+(1-self.g2.nu)/self.g2.g)
        R_eq = 1/(1/self.g1.radius+1/self.g2.radius)
        kt0 = 8 * G_eq *math.sqrt(R_eq*abs(self.overlap_normal))
        kt = kt0*math.sqrt(max(1-2/3*kt0*abs(self.overlap_tangential)/self.mu/abs(self.F_2_1_n),0)) #not linear spring

        r1 = self.g1.radius - self.overlap_normal/2
        r2 = self.g2.radius - self.overlap_normal/2
        Delta_Us = (np.dot(self.g1.v-self.g2.v,self.pc_tangential) + r1*self.g1.w + r2*self.g2.w)*dt_DEM
        self.overlap_tangential = self.overlap_tangential + Delta_Us
        self.ft = self.ft - kt*Delta_Us
        self.tangential_old = self.pc_tangential
        if abs(self.ft) > abs(self.mu*self.F_2_1_n) or kt == 0: #Coulomb criteria
            self.ft = self.mu * abs(self.F_2_1_n) * np.sign(self.ft)

        self.g1.add_F( self.ft*self.pc_tangential, self.g1.center + self.g1.radius*self.pc_normal)
        self.g2.add_F(-self.ft*self.pc_tangential, self.g2.center - self.g2.radius*self.pc_normal)

        #Damping term
        gamma = -math.log(self.coeff_restitution)/math.sqrt(math.pi**2+math.log(self.coeff_restitution)**2)
        mass_eq = self.g1.mass*self.g2.mass/(self.g1.mass+self.g2.mass)
        eta = 2 * gamma * math.sqrt(mass_eq*kt)
        F_2_1_damp_t = -Delta_Us/dt_DEM*eta/2
        F_2_1_damp = F_2_1_damp_t *self.pc_tangential
        self.ft_damp = F_2_1_damp_t
        self.g1.add_F( F_2_1_damp, self.g1.center + self.g1.radius*self.pc_normal)
        self.g2.add_F(-F_2_1_damp, self.g2.center - self.g2.radius*self.pc_normal)

    #no contact finally
    else :
        tangential_old_statut = False
        self.overlap_tangential = 0
        self.ft = 0
        self.ft_damp = 0

#-------------------------------------------------------------------------------
#Function
#-------------------------------------------------------------------------------

def Grains_contact_f(g1,g2):
  """
  Detect the contact grain-grain.

    Input :
        two temporary grains (two grain_tempos)
    Output :
        a Boolean, True if there is contact between the two grains (a Boolean)
  """
  return np.linalg.norm(g1.center-g2.center) < g1.radius+g2.radius:

#-------------------------------------------------------------------------------

def Update_Neighborhoods(dict_ic):
    """
    Determine a neighborhood for each grain.

    This function is called every x time step. The contact is determined by Grains_contact_Neighborhoods().
    Notice that if there is a potential contact between grain_i and grain_j, grain_i is not in the neighborhood of grain_j.
    Whereas grain_j is in the neighborhood of grain_i. With i_grain < j_grain.

        Input :
            an initial condition dictionnary (a dict)
        Output :
            Nothing, but the neighborhood of the temporary grains is updated
    """
    for i_grain in range(len(dict_ic['L_g_tempo'])-1) :
        neighborhood = []
        for j_grain in range(i_grain+1,len(dict_ic['L_g_tempo'])):
            if np.linalg.norm(dict_ic['L_g_tempo'][i_grain].center-dict_ic['L_g_tempo'][j_grain].center) < dict_ic['factor_neighborhood_IC']*(dict_ic['L_g_tempo'][i_grain].radius+dict_ic['L_g_tempo'][j_grain].radius):
                neighborhood.append(dict_ic['L_g_tempo'][j_grain])
        dict_ic['L_g_tempo'][i_grain].neighbourood = neighborhood


#-------------------------------------------------------------------------------

def Grains_contact_Neighborhoods(dict_ic,dict_material):
    """
    Detect contact between a grain and grains from its neighborhood.

    The neighborhood is updated with Update_Neighborhoods().

        Input :
            an initial condition dictionnary (a dict)
            a material dictionnary (a dict)
        Output :
            Nothing, but the initial condition dictionnary is updated with grain - grain contacts
    """
    for i_grain in range(len(dict_ic['L_g_tempo'])-1) :
        grain_i = dict_ic['L_g_tempo'][i_grain]
        for neighbor in dict_ic['L_g_tempo'][i_grain].neighbourood:
            j_grain = neighbor.id
            grain_j = neighbor
            if Grains_contact_f(grain_i,grain_j):
                if (i_grain,j_grain) not in dict_ic['L_contact_ij']:  #contact not detected previously
                   #creation of contact
                   dict_ic['L_contact_ij'].append((i_grain,j_grain))
                   dict_ic['L_contact'].append(Contact_Tempo(dict_ic['id_contact'], grain_i, grain_j, dict_material))
                   dict_ic['id_contact'] = dict_ic['id_contact'] + 1

            else :
                if (i_grain,j_grain) in dict_ic['L_contact_ij'] : #contact detected previously is not anymore
                       dict_ic['L_contact'].pop(dict_ic['L_contact_ij'].index((i_grain,j_grain)))
                       dict_ic['L_contact_ij'].remove((i_grain,j_grain))
