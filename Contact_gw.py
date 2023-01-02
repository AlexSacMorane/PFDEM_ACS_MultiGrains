# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

The goal of this file is to define a new class.
The new class is about the contact between a wall and a grain
"""

#-------------------------------------------------------------------------------
#Librairies
#-------------------------------------------------------------------------------

import numpy as np
import math

#-------------------------------------------------------------------------------
#Class
#-------------------------------------------------------------------------------

class Contact_gw:

#-------------------------------------------------------------------------------

  def __init__(self, ID, G, dict_material, Nature, Limit, Overlap):
    """
    Defining the contact grain - wall.

        Input :
            itself (a contact grain - wall)
            an id (a int)
            a grain (a grain)
            a material dictionnary (a dict)
            a nature to identify the wall (a string)
            the coordinate of the wall (a float)
            an overlap (a float)
        Output :
            a contact grain - wall is generated
    """
    self.id = ID
    self.g = G
    factor = 5 #factor just to increase the stiffness
    self.k = factor*4/3*self.g.y/(1-self.g.nu*self.g.nu)*math.sqrt(self.g.r_mean) #Hertz law
    self.kt = 0
    self.ft = 0
    self.limit = Limit
    self.nature = Nature
    self.mu = dict_material['mu_friction_gw']
    self.coeff_restitution = dict_material['coeff_restitution']
    self.overlap = Overlap
    self.tangential_old_statut = False
    self.overlap_tangential = 0

#-------------------------------------------------------------------------------

  def update_overlap(self,new_overlap):
    """
    Update the overlap of a contact already created.

        Input :
            itself (a contact grain - wall)
            a new overlap (a float)
        Output :
            Nothing, but the attribute is updated (a float)
    """
    self.overlap = new_overlap

#-------------------------------------------------------------------------------

  def init_contact_gw(self,L_g):
    """
    Initialize the contact.

    The grain is updated, the tangential reaction is set at 0 and a boolean is set at False (new contact grain-wall)

        Input :
            itself (a contact grain - wall)
            the list of the grain (a list)
        Output :
            Nothing, but the attributes are updated (a grain, a float and a boolean)
    """
    self.g = L_g[self.g.id]
    self.ft = 0
    self.tangential_old_statut = False

#-------------------------------------------------------------------------------

  def  DEM_gw_Polyhedral_normal(self):
    """
    Compute the normal reaction of a contact grain-wall.

    Here a pontual spring is considered.

        Input :
            itself (a contact grain - wall)
        Output :
            Nothing, but the interaction between grain and wall is computed. Differents attributes are updated.
    """
    #conditions "if" are defined and same for each wall nature
    if self.nature == 'gwy_min':
        #unlinear stiffness
        nwg = np.array([0,1])
        self.nwg = nwg
        Fwg_n = self.k*self.overlap**(3/2)
        Fwg = Fwg_n*nwg
        self.Fwg_n = Fwg_n
        self.g.update_f(Fwg[0],Fwg[1],self.g.l_border[self.g.l_border_y.index(min(self.g.l_border_y))])
        #damping
        gamma = -math.log(self.coeff_restitution)/math.sqrt(math.pi**2+math.log(self.coeff_restitution)**2)
        mass_eq = self.g.mass
        eta = 2 * gamma * math.sqrt(mass_eq*self.k)
        Fwg_damp_n = -np.dot(self.g.v,nwg)*eta
        Fwg_damp = Fwg_damp_n*nwg
        self.Fwg_damp_n = Fwg_damp_n
        self.g.update_f(Fwg_damp[0],Fwg_damp[1],self.g.l_border[self.g.l_border_y.index(min(self.g.l_border_y))])

    elif self.nature == 'gwy_max':
        #unlinear stiffness
        nwg = np.array([0,-1])
        self.nwg = nwg
        Fwg_n = self.k*self.overlap**(3/2)
        Fwg = Fwg_n*nwg
        self.Fwg_n = Fwg_n
        self.g.update_f(Fwg[0],Fwg[1],self.g.l_border[self.g.l_border_y.index(max(self.g.l_border_y))])
        #damping
        Fwg_damp_n = 0
        self.Fwg_damp_n = Fwg_damp_n

    elif self.nature == 'gwx_min':
        #unlinear stiffness
        nwg = np.array([1,0])
        self.nwg = nwg
        Fwg_n = self.k*self.overlap**(3/2)
        Fwg = Fwg_n*nwg
        self.Fwg_n = Fwg_n
        self.g.update_f(Fwg[0],Fwg[1],self.g.l_border[self.g.l_border_x.index(min(self.g.l_border_x))])
        #damping
        gamma = -math.log(self.coeff_restitution)/math.sqrt(math.pi**2+math.log(self.coeff_restitution)**2)
        mass_eq = self.g.mass
        eta = 2 * gamma * math.sqrt(mass_eq*self.k)
        Fwg_damp_n = -np.dot(self.g.v,nwg)*eta
        Fwg_damp = Fwg_damp_n*nwg
        self.Fwg_damp_n = Fwg_damp_n
        self.g.update_f(Fwg_damp[0],Fwg_damp[1],self.g.l_border[self.g.l_border_x.index(min(self.g.l_border_x))])

    elif self.nature == 'gwx_max':
        #unlinear stiffness
        nwg = np.array([-1,0])
        self.nwg = nwg
        Fwg_n = self.k*self.overlap**(3/2)
        Fwg = Fwg_n*nwg
        self.Fwg_n = Fwg_n
        self.g.update_f(Fwg[0],Fwg[1],self.g.l_border[self.g.l_border_x.index(max(self.g.l_border_x))])
        #damping
        gamma = -math.log(self.coeff_restitution)/math.sqrt(math.pi**2+math.log(self.coeff_restitution)**2)
        mass_eq = self.g.mass
        eta = 2 * gamma * math.sqrt(mass_eq*self.k)
        Fwg_damp_n = -np.dot(self.g.v,nwg)*eta
        Fwg_damp = Fwg_damp_n*nwg
        self.Fwg_damp_n = Fwg_damp_n
        self.g.update_f(Fwg_damp[0],Fwg_damp[1],self.g.l_border[self.g.l_border_x.index(max(self.g.l_border_x))])

#-------------------------------------------------------------------------------

  def DEM_gw_Polyhedral_tangential(self, dt_DEM):
   """
   Compute the tangential reaction of a contact grain-wall.

   Here a pontual spring is considered.

        Input :
            itself (a contact grain - wall)
        Output :
            Nothing, but the interaction between grain and wall is computed. Differents attributes are updated.
   """
   #conditions "if" are defined and same for each wall nature
   if self.nature == 'gwy_min':
       #unlinear stiffness
       twg = np.array([-1, 0])
       self.twg = twg
       r = np.linalg.norm(self.g.l_border[:-1][self.g.l_border_y.index(min(self.g.l_border_y))] - self.g.center) - self.overlap
       Delta_Us = (np.dot(self.g.v,self.twg) - r*self.g.w) * dt_DEM
       self.overlap_tangential = self.overlap_tangential + Delta_Us
       self.ft = self.ft - self.kt*Delta_Us
       if abs(self.ft) > abs(self.mu*self.Fwg_n) :
           self.ft = self.mu * abs(self.Fwg_n) * np.sign(self.ft)
       Fwg = self.ft*twg
       self.g.update_f(Fwg[0],Fwg[1],self.g.l_border[self.g.l_border_y.index(min(self.g.l_border_y))])

   elif self.nature == 'gwy_max':
       #unlinear stiffness
       twg = np.array([1, 0])
       self.twg = twg
       r = np.linalg.norm(self.g.l_border[:-1][self.g.l_border_y.index(max(self.g.l_border_y))] - self.g.center) - self.overlap
       Delta_Us = (np.dot(self.g.v,self.twg) - r*self.g.w) * dt_DEM
       self.overlap_tangential = self.overlap_tangential + Delta_Us
       self.ft = self.ft - self.kt*Delta_Us
       if abs(self.ft) > abs(self.mu*self.Fwg_n) :
           self.ft = self.mu * abs(self.Fwg_n) * np.sign(self.ft)
       Fwg = self.ft*twg
       self.g.update_f(Fwg[0],Fwg[1],self.g.l_border[self.g.l_border_y.index(max(self.g.l_border_y))])

   elif self.nature == 'gwx_min':
       #unlinear stiffness
       twg = np.array([0, 1])
       self.twg = twg
       r = np.linalg.norm(self.g.l_border[:-1][self.g.l_border_x.index(min(self.g.l_border_x))] - self.g.center) - self.overlap
       Delta_Us = (np.dot(self.g.v,self.twg) - r*self.g.w) * dt_DEM
       self.overlap_tangential = self.overlap_tangential + Delta_Us
       self.ft = self.ft - self.kt*Delta_Us
       if abs(self.ft) > abs(self.mu*self.Fwg_n) :
           self.ft = self.mu * abs(self.Fwg_n) * np.sign(self.ft)
       Fwg = self.ft*twg
       self.g.update_f(Fwg[0],Fwg[1],self.g.l_border[self.g.l_border_x.index(min(self.g.l_border_x))])

   elif self.nature == 'gwx_max':
       #linear stiffness
       twg = np.array([0, -1])
       self.twg = twg
       r = np.linalg.norm(self.g.l_border[:-1][self.g.l_border_x.index(max(self.g.l_border_x))] - self.g.center) - self.overlap
       Delta_Us = (np.dot(self.g.v,self.twg) - r*self.g.w) * dt_DEM
       self.overlap_tangential = self.overlap_tangential + Delta_Us
       self.ft = self.ft - self.kt*Delta_Us
       if abs(self.ft) > abs(self.mu*self.Fwg_n) :
           self.ft = self.mu * abs(self.Fwg_n) * np.sign(self.ft)
       Fwg = self.ft*twg
       self.g.update_f(Fwg[0],Fwg[1],self.g.l_border[self.g.l_border_x.index(max(self.g.l_border_x))])

#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------

def Grains_Polyhedral_Wall_contact_Neighborhood(dict_material,dict_sample):
  """
  Detect contact grain in the neighborhood of the wall and  the wall.

  The neighborhood is updated with Update_wall_Neighborhoods().
  An iteration over the grain list and a comparison with the coordinate of the different walls are done.

    Input :
        a material dictionnary (a dict)
        a sample dict (a dict)
    Output :
        Nothing, but the sample dict is updated with the different grain - wall in contact
  """
  for grain in dict_sample['wall_neighborhood']:

      p_x_min = min(grain.l_border_x)
      p_x_max = max(grain.l_border_x)
      p_y_min = min(grain.l_border_y)
      p_y_max = max(grain.l_border_y)

      #grain-wall x_min
      if p_x_min < dict_sample['x_box_min'] and (grain.id,-1) not in dict_sample['L_ij_contact_gw']:
          overlap = dict_sample['x_box_min'] - p_x_min
          dict_sample['L_contact_gw'].append(Contact_gw(dict_sample['id_contact_gw'], grain, dict_material, 'gwx_min', dict_sample['x_box_min'], overlap))
          dict_sample['L_ij_contact_gw'].append((grain.id,-1))
          dict_sample['id_contact_gw'] = dict_sample['id_contact_gw'] + 1
      elif p_x_min < dict_sample['x_box_min'] and (grain.id,-1) in dict_sample['L_ij_contact_gw']:
          overlap = dict_sample['x_box_min'] - p_x_min
          dict_sample['L_contact_gw'][dict_sample['L_ij_contact_gw'].index((grain.id,-1))].update_overlap(overlap)
      elif p_x_min > dict_sample['x_box_min'] and (grain.id,-1) in dict_sample['L_ij_contact_gw']:
          i_contact = dict_sample['L_ij_contact_gw'].index((grain.id,-1))
          dict_sample['L_contact_gw'].pop(i_contact)
          dict_sample['L_ij_contact_gw'].pop(i_contact)
      #grain-wall x_max
      if p_x_max > dict_sample['x_box_max'] and (grain.id,-2) not in dict_sample['L_ij_contact_gw']:
          overlap = p_x_max - dict_sample['x_box_max']
          dict_sample['L_contact_gw'].append(Contact_gw(dict_sample['id_contact_gw'], grain, dict_material, 'gwx_max', dict_sample['x_box_max'], overlap))
          dict_sample['L_ij_contact_gw'].append((grain.id,-2))
          dict_sample['id_contact_gw'] = dict_sample['id_contact_gw'] + 1
      elif p_x_max > dict_sample['x_box_max'] and (grain.id,-2) in dict_sample['L_ij_contact_gw']:
          overlap = p_x_max - dict_sample['x_box_max']
          dict_sample['L_contact_gw'][dict_sample['L_ij_contact_gw'].index((grain.id,-2))].update_overlap(overlap)
      elif p_x_max < dict_sample['x_box_max'] and (grain.id,-2) in dict_sample['L_ij_contact_gw']:
          i_contact = dict_sample['L_ij_contact_gw'].index((grain.id,-2))
          dict_sample['L_contact_gw'].pop(i_contact)
          dict_sample['L_ij_contact_gw'].pop(i_contact)
      #grain-wall y_min
      if p_y_min < dict_sample['y_box_min'] and (grain.id,-3) not in dict_sample['L_ij_contact_gw']:
          overlap = dict_sample['y_box_min'] - p_y_min
          dict_sample['L_contact_gw'].append(Contact_gw(dict_sample['id_contact_gw'], grain, dict_material, 'gwy_min', dict_sample['y_box_min'], overlap))
          dict_sample['L_ij_contact_gw'].append((grain.id,-3))
          dict_sample['id_contact_gw'] = dict_sample['id_contact_gw'] + 1
      elif p_y_min < dict_sample['y_box_min'] and (grain.id,-3) in dict_sample['L_ij_contact_gw']:
          overlap = dict_sample['y_box_min'] - p_y_min
          dict_sample['L_contact_gw'][dict_sample['L_ij_contact_gw'].index((grain.id,-3))].update_overlap(overlap)
      elif p_y_min > dict_sample['y_box_min'] and (grain.id,-3) in dict_sample['L_ij_contact_gw']:
          i_contact = dict_sample['L_ij_contact_gw'].index((grain.id,-3))
          dict_sample['L_contact_gw'].pop(i_contact)
          dict_sample['L_ij_contact_gw'].pop(i_contact)
      #grain-wall y_max
      if p_y_max > dict_sample['y_box_max'] and (grain.id,-4) not in dict_sample['L_ij_contact_gw']:
          overlap = p_y_max - dict_sample['y_box_max']
          dict_sample['L_contact_gw'].append(Contact_gw(dict_sample['id_contact_gw'], grain, dict_material, 'gwy_max', dict_sample['y_box_max'], overlap))
          dict_sample['L_ij_contact_gw'].append((grain.id,-4))
          dict_sample['id_contact_gw'] = dict_sample['id_contact_gw'] + 1
      elif p_y_max > dict_sample['y_box_max'] and (grain.id,-4) in dict_sample['L_ij_contact_gw']:
          overlap = p_y_max - dict_sample['y_box_max']
          dict_sample['L_contact_gw'][dict_sample['L_ij_contact_gw'].index((grain.id,-4))].update_overlap(overlap)
      elif p_y_max < dict_sample['y_box_max'] and (grain.id,-4) in dict_sample['L_ij_contact_gw']:
          i_contact = dict_sample['L_ij_contact_gw'].index((grain.id,-4))
          dict_sample['L_contact_gw'].pop(i_contact)
          dict_sample['L_ij_contact_gw'].pop(i_contact)

#-------------------------------------------------------------------------------

def Update_wall_Neighborhoods(dict_algorithm, dict_sample):
    """
    Determine a neighborhoods for wall.

    This function is called every x time step. A contact grain - wall is determined by Grains_Polyhedral_Wall_contact_Neighborhood().
    A factor determines the size of the neighborhood window.
    """
    wall_neighborhood = []
    for grain in dict_sample['L_g']:

        p_x_min = min(grain.l_border_x)
        p_x_max = max(grain.l_border_x)
        p_y_min = min(grain.l_border_y)
        p_y_max = max(grain.l_border_y)

        #grain-wall x_min
        if abs(p_x_min-dict_sample['x_box_min']) < dict_algorithm['factor_neighborhood']*grain.r_max :
            wall_neighborhood.append(grain)
        #grain-wall x_max
        if abs(p_x_max-dict_sample['x_box_max']) < dict_algorithm['factor_neighborhood']*grain.r_max :
            wall_neighborhood.append(grain)
        #grain-wall y_min
        if abs(p_y_min-dict_sample['y_box_min']) < dict_algorithm['factor_neighborhood']*grain.r_max :
            wall_neighborhood.append(grain)
        #grain-wall y_max
        if abs(p_y_max-dict_sample['y_box_max']) < dict_algorithm['factor_neighborhood']*grain.r_max :
            wall_neighborhood.append(grain)

    dict_sample['wall_neighborhood'] = wall_neighborhood
