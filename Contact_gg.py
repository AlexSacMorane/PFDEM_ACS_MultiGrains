# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

The goal of this file is to define a new class.
The new class is about the contact between two grains.
"""

#-------------------------------------------------------------------------------
#Librairies
#-------------------------------------------------------------------------------

import numpy as np
import math
import Grain

#-------------------------------------------------------------------------------
#Class
#-------------------------------------------------------------------------------

class Contact_gg:

#-------------------------------------------------------------------------------

  def __init__(self, ID, G1, G2, dict_material):
    """
    Defining the contact grain-grain.

        Input :
            itself (a contact grain - grain)
            an id (a int)
            two grains (two grains)
            a material dictionnary (a dict)
        Output :
            Nothing, but the contact grain - grain is generated
    """
    self.id = ID
    self.g1 = G1
    self.g2 = G2
    self.ft = 0
    self.mu = dict_material['mu_friction_gg']
    self.coeff_restitution = dict_material['coeff_restitution']
    self.tangential_old_statut = False
    self.overlap_tangential = 0

#-------------------------------------------------------------------------------

  def init_contact(self,L_g):
    """
    Initialize the contact.

    The grains are updated, the tangetial reaction is set to 0 and a boolean attribute is set to False (new contact grain - grain).

        Input :
            itself (a contact grain - grain)
            a list of grains (a list)
        Output :
            Nothing, but attributes are updated (two grains, two floats and a boolean)
    """
    self.g1 = L_g[self.g1.id]
    self.g2 = L_g[self.g2.id]
    self.ft = 0
    self.tangential_old_statut = False
    self.overlap_tangential = 0

#-------------------------------------------------------------------------------

  def DEM_gg_Polyhedral_normal(self):
    """
    Compute the normal reaction of a contact grain-grain

    Here a pontual spring is considered

        Input :
            itself (a contact grain - grain)
        Output :
            Nothing, but attributes are updated
    """
    #looking for the nearest nodes
    d_virtual = max(self.g1.r_max,self.g2.r_max)
    ij_min = [0,0]
    d_ij_min = 100*d_virtual #Large
    for i in range(len(self.g1.l_border[:-1])):
        for j in range(len(self.g2.l_border[:-1])):
            d_ij = np.linalg.norm(self.g2.l_border[:-1][j]-self.g1.l_border[:-1][i]+d_virtual*(self.g2.center-self.g1.center)/np.linalg.norm(self.g2.center-self.g1.center))
            if d_ij < d_ij_min :
                d_ij_min = d_ij
                ij_min = [i,j]
    self.ij_min = ij_min

    #-----------------------------------------------------------------------------
    #Computing CP
    #-----------------------------------------------------------------------------

    M = (self.g1.l_border[:-1][ij_min[0]]+self.g2.l_border[:-1][ij_min[1]])/2
    # 5 candidates for CP
    N = np.array([self.g1.l_border[:-1][ij_min[0]][0] - self.g2.l_border[:-1][ij_min[1]][0],
                  self.g1.l_border[:-1][ij_min[0]][1] - self.g2.l_border[:-1][ij_min[1]][1]])
    N = N/np.linalg.norm(N)
    PB = np.array([-N[1] ,N[0]])
    PB = PB/np.linalg.norm(PB)

    #candidats from grain 1
    if ij_min[0] <len(self.g1.l_border[:-1]) - 1:
        M1 = self.g1.l_border[:-1][ij_min[0]+1]-self.g1.l_border[:-1][ij_min[0]]
    else :
        M1 = self.g1.l_border[:-1][0]-self.g1.l_border[:-1][ij_min[0]]
    M1 = M1/np.linalg.norm(M1)
    M3 = self.g1.l_border[:-1][ij_min[0]-1]-self.g1.l_border[:-1][ij_min[0]]
    M3 = M3/np.linalg.norm(M3)
    #reorganize the candidats
    if np.dot(M1,PB) < 0:
        Mtempo = M1.copy()
        M1 = M3.copy()
        M3 = Mtempo.copy()

    #candidats from grain 2
    if ij_min[1] <len(self.g2.l_border[:-1]) - 1:
        M2 = self.g2.l_border[:-1][ij_min[1]+1]-self.g2.l_border[:-1][ij_min[1]]
    else :
        M2 = self.g2.l_border[:-1][0]-self.g2.l_border[:-1][ij_min[1]]
    M2 = M2/np.linalg.norm(M2)
    M4 = self.g2.l_border[:-1][ij_min[1]-1]-self.g2.l_border[:-1][ij_min[1]]
    M4 = M4/np.linalg.norm(M4)
    #reorganize the candidats
    if np.dot(M2,PB) < 0:
      Mtempo = M2.copy()
      M2 = M4.copy()
      M4 = Mtempo.copy()

    #compute the different angles
    theta_PB = math.pi/2
    theta_M1 =  math.acos(np.dot(M1,N))
    theta_M2 =  math.acos(np.dot(M2,N))
    theta_M3 = -math.acos(np.dot(M3,N))
    theta_M4 = -math.acos(np.dot(M4,N))

    #find the PC
    if theta_M2 < theta_PB and theta_PB < theta_M1\
       and theta_M3 < -theta_PB and -theta_PB < theta_M4:
       PC = PB
    else:
      L_Mi = [M1,M2,M3,M4]
      L_theta_Mi_PB=[theta_M1-theta_PB, theta_PB-theta_M2, -theta_M3-theta_PB, theta_PB+theta_M4]
      PC = L_Mi[L_theta_Mi_PB.index(min(L_theta_Mi_PB))]

    #-----------------------------------------------------------------------------
    # Compute the normal and tangential planes
    #-----------------------------------------------------------------------------

    PC_normal = np.array([PC[1],-PC[0]])
    if np.dot(PC_normal,(self.g2.center-self.g1.center)/np.linalg.norm(self.g2.center-self.g1.center))<0 :
        PC_normal = np.array([-PC[1],PC[0]])
    self.pc_normal = PC_normal #n12
    self.pc_tangential = np.array([-PC_normal[1],PC_normal[0]])

    #-----------------------------------------------------------------------------
    # Compute the overlap
    #-----------------------------------------------------------------------------

    d_b = np.dot(M-self.g2.l_border[:-1][ij_min[1]],PC_normal)
    d_a = np.dot(M-self.g1.l_border[:-1][ij_min[0]],PC_normal)
    overlap = d_b - d_a
    self.overlap_normal = overlap

    if overlap > 0:
    #-----------------------------------------------------------------------------
    # Compute the reaction
    #-----------------------------------------------------------------------------

        #Spring term
        Y_eq = 1/((1-self.g1.nu*self.g1.nu)/self.g1.y+(1-self.g2.nu*self.g2.nu)/self.g2.y)
        R_eq = 1/(1/self.g1.r_mean+1/self.g2.r_mean)
        k = 4/3*Y_eq*math.sqrt(R_eq)
        F_2_1_n = -k * overlap**(3/2)  #unlinear spring
        F_2_1 = F_2_1_n * PC_normal
        self.F_2_1_n = F_2_1_n
        self.Ep_n = 2/5 * k * overlap**(5/2) #-dEp/dx = F_2_1_n
        self.g1.update_f( F_2_1[0],  F_2_1[1], self.g1.l_border[:-1][ij_min[0]])
        self.g2.update_f(-F_2_1[0], -F_2_1[1], self.g2.l_border[:-1][ij_min[1]])

        #Damping term
        gamma = -math.log(self.coeff_restitution)/math.sqrt(math.pi**2+math.log(self.coeff_restitution)**2)
        mass_eq = self.g1.mass*self.g2.mmass/(self.g1.mmass+self.g2.mmass)
        eta = 2 * gamma * math.sqrt(mass_eq*k)
        F_2_1_damp_n = np.dot(self.g2.v - self.g1.v,PC_normal)*eta
        F_2_1_damp = F_2_1_damp_n *PC_normal
        self.F_2_1_damp = F_2_1_damp_n
        self.g1.update_f( F_2_1_damp[0],  F_2_1_damp[1], self.g1.l_border[:-1][ij_min[0]])
        self.g2.update_f(-F_2_1_damp[0], -F_2_1_damp[1], self.g2.l_border[:-1][ij_min[1]])

    #no contact finally
    else :
        self.F_2_1_n = 0
        self.F_2_1_damp = 0
        self.Ep_n = 0

#-------------------------------------------------------------------------------

  def DEM_gg_Polyhedral_tangential(self,dt_DEM):
    """
    Compute the tangential reaction of a contact grain-grain.

    Here a pontual spring is considered.

        Input :
            itself (a contact grain - grain)
            a time step (a float)
        Output :
            Nothing, but attributes are updated
    """
    if self.overlap_normal > 0:

        if self.tangential_old_statut:
          #if a reaction has been already computed
          #need to project the tangential reaction on the new tangential plane
          self.ft = self.ft*np.dot(self.tangential_old,self.pc_tangential)
        else:
          self.tangential_old_statut = True

        G_eq = 1/((1-self.g1.nu)/self.g1.g+(1-self.g2.nu)/self.g2.g)
        R_eq = 1/(1/self.g1.r_mean+1/self.g2.r_mean)
        kt0 = 8 * G_eq *math.sqrt(R_eq*abs(self.overlap_normal))
        kt = kt0*math.sqrt(max(1-2/3*kt0*abs(self.overlap_tangential)/self.mu/abs(self.F_2_1_n),0)) #not linear spring

        r1 = np.linalg.norm(self.g1.l_border[:-1][self.ij_min[0]] - self.g1.center) - self.overlap_normal/2
        r2 = np.linalg.norm(self.g2.l_border[:-1][self.ij_min[1]] - self.g2.center) - self.overlap_normal/2
        Delta_Us = (np.dot(self.g1.v-self.g2.v,self.pc_tangential) + r1*self.g1.w + r2*self.g2.w)*dt_DEM
        self.overlap_tangential = self.overlap_tangential + Delta_Us
        self.ft = self.ft - kt*Delta_Us
        self.tangential_old = self.pc_tangential
        if abs(self.ft) > abs(self.mu*self.F_2_1_n) or kt == 0: #Coulomb criteria
            self.ft = self.mu * abs(self.F_2_1_n) * np.sign(self.ft)

        self.g1.update_f( self.ft*self.pc_tangential[0],  self.ft*self.pc_tangential[1], self.g1.l_border[:-1][self.ij_min[0]])
        self.g2.update_f(-self.ft*self.pc_tangential[0], -self.ft*self.pc_tangential[1], self.g2.l_border[:-1][self.ij_min[1]])

        #Damping term
        gamma = -math.log(self.coeff_restitution)/math.sqrt(math.pi**2+math.log(self.coeff_restitution)**2)
        mass_eq = self.g1.mass*self.g2.mass/(self.g1.mass+self.g2.mass)
        eta = 2 * gamma * math.sqrt(mass_eq*kt)
        F_2_1_damp_t = -Delta_Us/dt_DEM*eta/2
        F_2_1_damp = F_2_1_damp_t *self.pc_tangential
        self.ft_damp = F_2_1_damp_t
        self.g1.update_f( F_2_1_damp[0],  F_2_1_damp[1], self.g1.l_border[:-1][self.ij_min[0]])
        self.g2.update_f(-F_2_1_damp[0], -F_2_1_damp[1], self.g2.l_border[:-1][self.ij_min[1]])

    #no contact finally
    else :
        tangential_old_statut = False
        self.overlap_tangential = 0
        self.ft = 0
        self.ft_damp = 0

#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------

def Update_Neighborhoods(dict_algorithm,dict_sample):
    """
    Determine a neighborhood for each grain.

    This function is called every x time step. The contact is determined by Grains_Polyhedral_contact_Neighborhoods().

    Notice that if there is a potential contact between grain_i and grain_j, grain_i is not in the neighborhood of grain_j.
    Whereas grain_j is in the neighbourood of grain_i, with i_grain < j_grain

        Input :
            an algorithm dictionnary (a dict)
            a sample dictionnary (a dict)
        Output :
            Nothing, but attribut is updated for each grains (a list)
    """
    for i_grain in range(len(dict_sample['L_g'])-1) :
        neighborhood = []
        for j_grain in range(i_grain+1,len(dict_sample['L_g'])):
            if np.linalg.norm(dict_sample['L_g'][i_grain].center-dict_sample['L_g'][j_grain].center) < dict_algorithm['factor_neighborhood']*(dict_sample['L_g'][i_grain].r_max+dict_sample['L_g'][j_grain].r_max):
                neighborhood.append(dict_sample['L_g'][j_grain])
        dict_sample['L_g'][i_grain].neighbourood = neighborhood

#-------------------------------------------------------------------------------

def Grains_Polyhedral_contact_Neighborhoods_bool(g1,g2):
  """
  Detect contact grain-grain.

    Input :
        two grains (two grains)
    Output :
        a Boolean (True if there is a contact)
  """
  #looking for the nearest nodes
  d_virtual = max(g1.r_max,g2.r_max)
  ij_min = [0,0]
  d_ij_min = 100*d_virtual #Large
  for i in range(len(g1.l_border[:-1])):
    for j in range(len(g2.l_border[:-1])):
        d_ij = np.linalg.norm(g2.l_border[:-1][j]-g1.l_border[:-1][i]+d_virtual*(g2.center-g1.center)/np.linalg.norm(g2.center-g1.center))
        if d_ij < d_ij_min :
            d_ij_min = d_ij
            ij_min = [i,j]

  d_ij_min = np.dot(g2.l_border[:-1][ij_min[1]]-g1.l_border[:-1][ij_min[0]],-(g2.center-g1.center)/np.linalg.norm(g2.center-g1.center))
  return d_ij_min > 0

#-------------------------------------------------------------------------------

def Grains_Polyhedral_contact_Neighborhoods(dict_material,dict_sample):
    """
    Detect contact between a grain and grains from its neighborhood.

    The neighbourood is updated with Update_Neighborhoods_f()

        Input :
            a material dictionnary (a dict)
            a sample dictionnary (a dict)
        Output :
            Nothing, but the sample dictionnary is updated with a list of contacts grain - grain
    """
    for i_grain in range(len(dict_sample['L_g'])-1) :
        for neighbour in dict_sample['L_g'][i_grain].neighbourood:
            j_grain = neighbour.id
            if Grains_Polyhedral_contact_Neighborhoods_bool(dict_sample['L_g'][i_grain],dict_sample['L_g'][j_grain]):
                if (i_grain,j_grain) not in dict_sample['L_ij_contact']:  #contact not detected previously
                   #creation of contact
                   dict_sample['L_ij_contact'].append((i_grain,j_grain))
                   dict_sample['L_contact'].append(Contact_gg(dict_sample['id_contact'], dict_sample['L_g'][i_grain], dict_sample['L_g'][j_grain], dict_material))
                   dict_sample['id_contact'] = dict_sample['id_contact'] + 1

            else :
                if (i_grain,j_grain) in dict_sample['L_ij_contact'] : #contact detected previously is not anymore
                       dict_sample['L_contact'].pop(dict_sample['L_ij_contact'].index((i_grain,j_grain)))
                       dict_sample['L_ij_contact'].remove((i_grain,j_grain))
