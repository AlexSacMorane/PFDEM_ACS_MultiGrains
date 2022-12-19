# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This file contains ??.nt functions used in the simulation.
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

import numpy  as np
import math

#Own
import Create_IC.Grain_ic

#-------------------------------------------------------------------------------
#Class
#-------------------------------------------------------------------------------

class Contact_gw_Tempo:
  """
  A temporary contact grain - wall used to generated an initial condition.
  """

#-------------------------------------------------------------------------------

  def __init__(self, ID, G, dict_material, Nature, Limit, Overlap):
    """
    Defining the contact grain-wall.

        Input :
             itself (a contact_gw_tempo)
             an id (a int)
             a grain (a grain_tempo)
             a material dictionnary (a dict)
             the nature of the wall (a string)
             the coordinate of the wall (a float)
             an overlap (a float)
    """
    self.id = ID
    self.g = G
    factor = 5 #factor just to increase the stiffness
    self.k = factor*4/3*self.g.y/(1-self.g.nu*self.g.nu)*math.sqrt(self.g.radius) #Hertz law
    self.kt = 0
    self.ft = 0
    self.limit = Limit
    self.nature = Nature
    self.mu = 0
    self.coeff_restitution = dict_material['coeff_restitution']
    self.overlap = Overlap
    self.tangential_old_statut = False
    self.overlap_tangential = 0

#-------------------------------------------------------------------------------

  def update_overlap(self,new_overlap):
    '''
    Update the overlap of a contact already created.

        Input :
            itself (a contact_gw_tempo)
            an overlap (a float)
        Output :
            Nothing, but the attribut concerning the overlap is updated (a float)
    '''
    self.overlap = new_overlap

#-------------------------------------------------------------------------------

  def  normal(self):
    """
    Compute the normal reaction of a contact grain-wall.

    Here a pontual spring is considered

        Input :
            itself (a contact_gw_tempo)
        Output :
            Nothing, but attributes are updated
    """
    #conditions "if" are defined and same for each wall nature
    if self.nature == 'gwy_min':
        #unlinear stiffness
        nwg = np.array([0,1])
        self.nwg = nwg
        Fwg_n = self.k*self.overlap**(3/2)
        Fwg = Fwg_n*nwg
        self.Fwg_n = Fwg_n
        self.g.add_F(Fwg, self.g.center - self.g.radius*self.nwg)
        #damping
        gamma = -math.log(self.coeff_restitution)/math.sqrt(math.pi**2+math.log(self.coeff_restitution)**2)
        mass_eq = self.g.mass
        eta = 2 * gamma * math.sqrt(mass_eq*self.k)
        Fwg_damp_n = -np.dot(self.g.v,nwg)*eta
        Fwg_damp = Fwg_damp_n*nwg
        self.Fwg_damp_n = Fwg_damp_n
        self.g.add_F(Fwg_damp, self.g.center - self.g.radius*self.nwg)

    elif self.nature == 'gwy_max':
        #unlinear stiffness
        nwg = np.array([0,-1])
        self.nwg = nwg
        Fwg_n = self.k*self.overlap**(3/2)
        Fwg = Fwg_n*nwg
        self.Fwg_n = Fwg_n
        self.g.add_F(Fwg, self.g.center - self.g.radius*self.nwg)
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
        self.g.add_F(Fwg, self.g.center - self.g.radius*self.nwg)
        #damping
        gamma = -math.log(self.coeff_restitution)/math.sqrt(math.pi**2+math.log(self.coeff_restitution)**2)
        mass_eq = self.g.mass
        eta = 2 * gamma * math.sqrt(mass_eq*self.k)
        Fwg_damp_n = -np.dot(self.g.v,nwg)*eta
        Fwg_damp = Fwg_damp_n*nwg
        self.Fwg_damp_n = Fwg_damp_n
        self.g.add_F(Fwg_damp, self.g.center - self.g.radius*self.nwg)

    elif self.nature == 'gwx_max':
        #unlinear stiffness
        nwg = np.array([-1,0])
        self.nwg = nwg
        Fwg_n = self.k*self.overlap**(3/2)
        Fwg = Fwg_n*nwg
        self.Fwg_n = Fwg_n
        self.g.add_F(Fwg, self.g.center - self.g.radius*self.nwg)
        #damping
        gamma = -math.log(self.coeff_restitution)/math.sqrt(math.pi**2+math.log(self.coeff_restitution)**2)
        mass_eq = self.g.mass
        eta = 2 * gamma * math.sqrt(mass_eq*self.k)
        Fwg_damp_n = -np.dot(self.g.v,nwg)*eta
        Fwg_damp = Fwg_damp_n*nwg
        self.Fwg_damp_n = Fwg_damp_n
        self.g.add_F(Fwg_damp, self.g.center - self.g.radius*self.nwg)

#-------------------------------------------------------------------------------

  def tangential(self, dt_DEM):
   """
   Compute the tangential reaction of a contact grain-wall.

   Here a pontual spring is considered.

        Input :
            itself (a contact_gw_tempo)
            a time step (a float)
        Output :
            Nothing, but attributes are updated
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
       self.g.add_F(Fwg, self.g.center - self.g.radius*self.nwg)

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
       self.g.add_F(Fwg, self.g.center - self.g.radius*self.nwg)

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
       self.g.add_F(Fwg, self.g.center - self.g.radius*self.nwg)

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
       self.g.add_F(Fwg, self.g.center - self.g.radius*self.nwg)

#-------------------------------------------------------------------------------
#Function
#-------------------------------------------------------------------------------

def Update_wall_Neighborhoods(L_g_tempo,factor_neighborhood_IC,x_min,x_max,y_min,y_max):
    """
    Determine a neighborhood for wall.

    This function is called every x time step. The grain - wall contact is determined by Grains_Polyhedral_Wall_contact_Neighborhood().
    A factor determines the size of the neighborhood window.

        Input :
            a list of temporary grains (a list)
            a factor to determine the neighborhood window (a float)
            the coordinates of the left, right, lower, upper walls (four floats)
        Output :
            a list of temporary grains in the neighborhood of the walls (a list)
    """
    wall_neighborhood = []
    for grain in L_g_tempo:

        p_x_min = min(grain.l_border_x)
        p_x_max = max(grain.l_border_x)
        p_y_min = min(grain.l_border_y)
        p_y_max = max(grain.l_border_y)

        #grain-wall x_min
        if abs(p_x_min-x_min) < factor_neighborhood_IC*grain.radius :
            wall_neighborhood.append(grain)
        #grain-wall x_max
        if abs(p_x_max-x_max) < factor_neighborhood_IC*grain.radius :
            wall_neighborhood.append(grain)
        #grain-wall y_min
        if abs(p_y_min-y_min) < factor_neighborhood_IC*grain.radius :
            wall_neighborhood.append(grain)
        #grain-wall y_max
        if abs(p_y_max-y_max) < factor_neighborhood_IC*grain.radius :
            wall_neighborhood.append(grain)

    return wall_neighborhood

#-------------------------------------------------------------------------------

def Grains_Polyhedral_Wall_contact_Neighborhood(wall_neighborhood,x_box_min,x_box_max,y_box_min,y_box_max, dict_ic, dict_material):
  """
  Detect contact grain in the neighborhood of the wall and the wall.

  The neighbourood is updated with Update_wall_Neighborhoods(). An iteration over the grains in the wall neighborhood is done. A comparison is done with the coordinates of the wall to determine if there is a contact.

        Input :
            a walls neighborhood (a list)
            the coordinates of the walls (four floats)
            an initial condition dictionnary (a dict)
            a material dictionnary (a dict)
        Output :
            Nothing, but the initial condition dictionnary is updated with the contact grain - walls.
  """
  #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
  #load data needed
  L_ij_contact_gw = dict_ic['L_contact_gw_ij']
  L_contact_gw = dict_ic['L_contact_gw']
  id_contact = dict_ic['id_contact']
  #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-

  for grain in wall_neighborhood:

      p_x_min = grain.center[0] - grain.radius
      p_x_max = grain.center[0] + grain.radius
      p_y_min = grain.center[1] - grain.radius
      p_y_max = grain.center[1] + grain.radius

      #grain-wall x_min
      if p_x_min < x_box_min and (grain.id,-1) not in L_ij_contact_gw:
          overlap = x_box_min - p_x_min
          L_contact_gw.append(Contact_gw_Tempo(id_contact, grain, dict_material, 'gwx_min', x_box_min, overlap))
          L_ij_contact_gw.append((grain.id,-1))
          id_contact = id_contact + 1
      elif p_x_min < x_box_min and (grain.id,-1) in L_ij_contact_gw:
          overlap = x_box_min - p_x_min
          L_contact_gw[L_ij_contact_gw.index((grain.id,-1))].update_overlap(overlap)
      elif p_x_min > x_box_min and (grain.id,-1) in L_ij_contact_gw:
          i_contact = L_ij_contact_gw.index((grain.id,-1))
          L_contact_gw.pop(i_contact)
          L_ij_contact_gw.pop(i_contact)
      #grain-wall x_max
      if p_x_max > x_box_max and (grain.id,-2) not in L_ij_contact_gw:
          overlap = p_x_max - x_box_max
          L_contact_gw.append(Contact_gw_Tempo(id_contact, grain, dict_material, 'gwx_max', x_box_max, overlap))
          L_ij_contact_gw.append((grain.id,-2))
          id_contact = id_contact + 1
      elif p_x_max > x_box_max and (grain.id,-2) in L_ij_contact_gw:
          overlap = p_x_max - x_box_max
          L_contact_gw[L_ij_contact_gw.index((grain.id,-2))].update_overlap(overlap)
      elif p_x_max < x_box_max and (grain.id,-2) in L_ij_contact_gw:
          i_contact = L_ij_contact_gw.index((grain.id,-2))
          L_contact_gw.pop(i_contact)
          L_ij_contact_gw.pop(i_contact)
      #grain-wall y_min
      if p_y_min < y_box_min and (grain.id,-3) not in L_ij_contact_gw:
          overlap = y_box_min - p_y_min
          L_contact_gw.append(Contact_gw_Tempo(id_contact, grain, dict_material, 'gwy_min', y_box_min, overlap))
          L_ij_contact_gw.append((grain.id,-3))
          id_contact = id_contact + 1
      elif p_y_min < y_box_min and (grain.id,-3) in L_ij_contact_gw:
          overlap = y_box_min - p_y_min
          L_contact_gw[L_ij_contact_gw.index((grain.id,-3))].update_overlap(overlap)
      elif p_y_min > y_box_min and (grain.id,-3) in L_ij_contact_gw:
          i_contact = L_ij_contact_gw.index((grain.id,-3))
          L_contact_gw.pop(i_contact)
          L_ij_contact_gw.pop(i_contact)
      #grain-wall y_max
      if p_y_max > y_box_max and (grain.id,-4) not in L_ij_contact_gw:
          overlap = p_y_max - y_box_max
          L_contact_gw.append(Contact_gw_Tempo(id_contact, grain, dict_material, 'gwy_max', y_box_max, overlap))
          L_ij_contact_gw.append((grain.id,-4))
          id_contact = id_contact + 1
      elif p_y_max > y_box_max and (grain.id,-4) in L_ij_contact_gw:
          overlap = p_y_max - y_box_max
          L_contact_gw[L_ij_contact_gw.index((grain.id,-4))].update_overlap(overlap)
      elif p_y_max < y_box_max and (grain.id,-4) in L_ij_contact_gw:
          i_contact = L_ij_contact_gw.index((grain.id,-4))
          L_contact_gw.pop(i_contact)
          L_ij_contact_gw.pop(i_contact)

      #Update dict
      dict_ic['L_contact_gw_ij'] = L_ij_contact_gw
      dict_ic['L_contact_gw'] = L_contact_gw
      dict_ic['id_contact'] = id_contact
