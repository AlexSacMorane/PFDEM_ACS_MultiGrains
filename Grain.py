# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

The goal of this file is to define a new class.
The new class is about the grains
"""

#-------------------------------------------------------------------------------
#Libs
#-------------------------------------------------------------------------------

import numpy as np
import math
import random
from scipy.ndimage import binary_dilation

#Own  functions and classes
import Owntools

#-------------------------------------------------------------------------------
#Class
#-------------------------------------------------------------------------------

class Grain:

    #---------------------------------------------------------------------------

    def __init__(self, dict_ic_to_real, dict_material, dict_sample):
        '''
        Defining a disk grain

            Input :
                an id (an integer)
                a radius (a float)
                a center (a 2 x 1 numpy array)
                a material dictionnary (a dictionnary)
                a sample dictionnary (a dictionnary)
            Output :
                a grain (a grain)
        '''
        self.id = dict_ic_to_real['Id']
        self.center = dict_ic_to_real['Center']
        #save description
        self.r_mean = np.mean(dict_ic_to_real['L_r'])
        self.r_min = min(dict_ic_to_real['L_r'])
        self.r_max = np.max(dict_ic_to_real['L_r'])
        self.l_r = dict_ic_to_real['L_r']
        self.l_theta_r = dict_ic_to_real['L_theta_r']
        self.l_border = dict_ic_to_real['L_border']
        self.l_border_x = dict_ic_to_real['L_border_x']
        self.l_border_y = dict_ic_to_real['L_border_y']
        self.surface = dict_ic_to_real['Surface']
        self.inertia = dict_ic_to_real['Inertia']
        #save initial
        self.center_init = self.center.copy()
        self.l_border_x_init = self.l_border_x.copy()
        self.l_border_y_init = self.l_border_y.copy()
        #material
        self.y = dict_material['Y']
        self.nu = dict_material['nu']
        self.g = self.y/2/(1+self.nu) #shear modulus
        self.rho_surf = dict_ic_to_real['Rho_surf']
        self.mass = dict_ic_to_real['Mass']
        #kinetic
        self.theta = 0
        self.v = np.array([0, 0])
        self.w = 0

        self.build_etai_M(dict_material,dict_sample)

    #---------------------------------------------------------------------------

    def __eq__(self, other):
        """
        Define the equality between two grains.

        To be equal two grains must have the same id.

          Input :
              two grains (grain)
          Return
              A Boolean (bool)
        """
        return self.id == other.id

    #---------------------------------------------------------------------------

    def build_etai_M(self,dict_material,dict_sample):
        '''
        Build the phase field for one grain.

        A cosine profile is assumed (see https://mooseframework.inl.gov/source/ics/SmoothCircleIC.html).

            Input :
                itself (a grain)
                a material dictionnary (a dictionnary)
                a sample dictionnary (a dictionnary)
            Output :
                Nothing but the grain gets a new attribute (a n_y x n_x numpy array)
        '''
        #initialization
        self.etai_M = np.array(np.zeros((len(dict_sample['y_L']),len(dict_sample['x_L']))))

        #extract a spatial zone
        x_min = min(self.l_border_x)-dict_material['w']
        x_max = max(self.l_border_x)+dict_material['w']
        y_min = min(self.l_border_y)-dict_material['w']
        y_max = max(self.l_border_y)+dict_material['w']

        #look for this part inside the global mesh
        #create search list
        x_L_search_min = abs(np.array(dict_sample['x_L'])-x_min)
        x_L_search_max = abs(np.array(dict_sample['x_L'])-x_max)
        y_L_search_min = abs(np.array(dict_sample['y_L'])-y_min)
        y_L_search_max = abs(np.array(dict_sample['y_L'])-y_max)

        #get index
        i_x_min = list(x_L_search_min).index(min(x_L_search_min))
        i_x_max = list(x_L_search_max).index(min(x_L_search_max))
        i_y_min = list(y_L_search_min).index(min(y_L_search_min))
        i_y_max = list(y_L_search_max).index(min(y_L_search_max))

        for l in range(i_y_min,i_y_max+1):
            for c in range(i_x_min,i_x_max+1):
                y = dict_sample['y_L'][l]
                x = dict_sample['x_L'][c]
                p = np.array([x,y])
                r = np.linalg.norm(self.center - p)
                #look for the radius on this direction
                if p[1]>self.center[1]:
                    theta = math.acos((p[0]-self.center[0])/np.linalg.norm(self.center-p))
                else :
                    theta= 2*math.pi - math.acos((p[0]-self.center[0])/np.linalg.norm(self.center-p))
                L_theta_R_i = list(abs(np.array(self.l_theta_r)-theta))
                R = self.l_r[L_theta_R_i.index(min(L_theta_R_i))]
                #build etai_M
                self.etai_M[-1-l][c] = Owntools.Cosine_Profile(R,r,dict_material['w'])

    #---------------------------------------------------------------------------

    def geometric_study(self,dict_sample):
      '''
      Searching limits of the grain

      Not best method but first approach
      We iterate on y constant, we look for a value under and over 0.5
      If both conditions are verified, there is a limit at this y
      Same with iteration on x constant

      Once the border of the grain is defined, a Monte Carlo method is used to computed some geometrical properties.

        Input :
            itself (a grain)
            a sample dictionnary (a dictionnary)
        Output :
            Nothing but the grain gets new attributes
                r_min : the minimum radius of the grain (a float)
                r_max : the maximum radius of the grain (a float)
                r_mean : the mean radius of the grain (a float)
                l_r : a list of radius of the grain, work with l_theta_r (a list)
                l_theta_r : a list of angle to see the distribution of the radius of the grain, work with l_r (a list)
                surface : the surface of the grain (a float)
                center : the coordinate of the grain center (a 2 x 1 numpy array)
                l_border_x : the list of the coordinate x of the grain vertices (a list)
                l_border_y : the list of the coordinate y of the grain vertices (a list)
                l_border : the list of the coordinate [x,y] of the grain vertices (a list)
      '''
      #-------------------------------------------------------------------------
      #load data needed
      n = dict_sample['grain_discretisation']
      x_L = dict_sample['x_L']
      y_L = dict_sample['y_L']
      #-------------------------------------------------------------------------

      L_border_old = []
      for y_i in range(len(y_L)):
          L_extract_x = self.etai_M[y_i][:]
          if max(L_extract_x)>0.5 and min(L_extract_x)<0.5:
              y_intersect = y_L[len(y_L)-1-y_i]
              for x_i in range(len(x_L)-1):
                  if (L_extract_x[x_i]-0.5)*(L_extract_x[x_i+1]-0.5)<0:
                      x_intersect = (0.5-L_extract_x[x_i])/(L_extract_x[x_i+1]-L_extract_x[x_i])*\
                                  (x_L[x_i+1]-x_L[x_i]) + x_L[x_i]
                      L_border_old.append(np.array([x_intersect,y_intersect]))

      for x_i in range(len(x_L)):
          L_extract_y = []
          for y_i in range(len(y_L)):
              L_extract_y.append(self.etai_M[y_i][x_i])
          if max(L_extract_y)>0.5 and min(L_extract_y)<0.5:
              x_intersect = x_L[x_i]
              for y_i in range(len(y_L)-1):
                  if (L_extract_y[y_i]-0.5)*(L_extract_y[y_i+1]-0.5)<0:
                      y_intersect = (0.5-L_extract_y[y_i])/(L_extract_y[y_i+1]-L_extract_y[y_i])*\
                                  (y_L[len(y_L)-1-y_i-1]-y_L[len(y_L)-1-y_i]) + y_L[len(y_L)-1-y_i]
                      L_border_old.append(np.array([x_intersect,y_intersect]))

      #Adaptating
      L_id_used = [0]
      L_border = [L_border_old[0]]
      HighValue = 100000000 #Large

      current_node = L_border_old[0]
      for j in range(1,len(L_border_old)):
          L_d = list(np.zeros(len(L_border_old)))
          for i in range(0,len(L_border_old)):
              node = L_border_old[i]
              if  i not in L_id_used:
                  d = np.linalg.norm(node - current_node)
                  L_d[i] = d
              else :
                  L_d[i] = HighValue #Value need to be larger than potential distance between node

          index_nearest_node = L_d.index(min(L_d))
          nearest_node = L_border_old[index_nearest_node]
          current_node = nearest_node
          L_border.append(nearest_node)
          L_id_used.append(index_nearest_node)

      #Correcting
      L_d_final = []
      for i in range(len(L_border)-1):
          L_d_final.append(np.linalg.norm(L_border[i+1] - L_border[i]))

      #look for really far points, we assume the first point is accurate
      d_final_mean = np.mean(L_d_final)
      while np.max(L_d_final) > 5 * d_final_mean : #5 here is an user choixe value
          i_error = L_d_final.index(np.max(L_d_final))+1
          #simulation_report.write('Point '+str(L_border[i_error])+' is deleted because it is detected as an error\n')
          L_border.pop(i_error)
          L_id_used.pop(i_error)
          L_d_final = []
          for i in range(len(L_border)-1):
              L_d_final.append(np.linalg.norm(L_border[i+1] - L_border[i]))

      #-------------------------------------------------------------------------------
      #Reduce the number of nodes for a grain
      #-------------------------------------------------------------------------------

      Perimeter = 0
      for i_p in range(len(L_border)-1):
          Perimeter = Perimeter + np.linalg.norm(L_border[i_p+1]-L_border[i_p])
      Perimeter = Perimeter + np.linalg.norm(L_border[-1]-L_border[0])
      distance_min = Perimeter/n
      L_border_adapted = [L_border[0]]
      for p in L_border[1:]:
          distance = np.linalg.norm(p-L_border_adapted[-1])
          if distance >= distance_min:
              L_border_adapted.append(p)
      L_border = L_border_adapted
      L_border.append(L_border[0])
      self.l_border = L_border

      #-------------------------------------------------------------------------------
      #Searching Surface, Center of mass and Inertia.
      #Monte Carlo Method
      #A box is defined, we take a random point and we look if it is inside or outside the grain
      #Properties are the statistic times the box properties
      #-------------------------------------------------------------------------------

      min_max_defined = False
      for p in L_border[:-1] :
          if not min_max_defined:
              box_min_x = p[0]
              box_max_x = p[0]
              box_min_y = p[1]
              box_max_y = p[1]
              min_max_defined = True
          else:
              if p[0] < box_min_x:
                  box_min_x = p[0]
              elif p[0] > box_max_x:
                  box_max_x = p[0]
              if p[1] < box_min_y:
                  box_min_y = p[1]
              elif p[1] > box_max_y:
                  box_max_y = p[1]

      N_MonteCarlo = 3000 #The larger it is, the more accurate it is
      sigma = self.rho_surf
      M_Mass = 0
      M_Center_Mass = np.array([0,0])
      M_Inertia = 0

      for i in range(N_MonteCarlo):
          P = np.array([random.uniform(box_min_x,box_max_x),random.uniform(box_min_y,box_max_y)])
          if self.P_is_inside(P):
              M_Mass = M_Mass + sigma
              M_Center_Mass = M_Center_Mass + sigma*P
              M_Inertia = M_Inertia + sigma*np.dot(P,P)

      Mass = (box_max_x-box_min_x)*(box_max_y-box_min_y)/N_MonteCarlo*M_Mass
      Center_Mass = (box_max_x-box_min_x)*(box_max_y-box_min_y)/N_MonteCarlo*M_Center_Mass/Mass
      Inertia = (box_max_x-box_min_x)*(box_max_y-box_min_y)/N_MonteCarlo*M_Inertia-Mass*np.dot(Center_Mass,Center_Mass)

      #-------------------------------------------------------------------------------
      #Updating the grain geometry and properties
      #-------------------------------------------------------------------------------

      L_R = []
      L_theta_R = []
      L_border_x = []
      L_border_y = []
      for p in L_border[:-1]:
          L_R.append(np.linalg.norm(p-Center_Mass))
          L_border_x.append(p[0])
          L_border_y.append(p[1])
          if (p-Center_Mass)[1] > 0:
              theta = math.acos((p-Center_Mass)[0]/np.linalg.norm(p-Center_Mass))
          else :
              theta = 2*math.pi - math.acos((p-Center_Mass)[0]/np.linalg.norm(p-Center_Mass))
          L_theta_R.append(theta)
      L_border_x.append(L_border_x[0])
      L_border_y.append(L_border_y[0])
      #reorganize lists
      i_theta = 0
      while L_theta_R[i_theta] < 2*math.pi/10 or 2*math.pi*9/10 < L_theta_R[i_theta]:
          i_theta = i_theta + 1
      if L_theta_R[i_theta + 1] < L_theta_R[i_theta] :
          L_R.reverse()
          L_theta_R.reverse()
          L_border_x.reverse()
          L_border_y.reverse()
          L_border.reverse()

      self.r_min = np.min(L_R)
      self.r_max = np.max(L_R)
      self.r_mean = np.mean(L_R)
      self.l_r = L_R
      self.l_theta_r = L_theta_R
      self.mass = Mass
      self.surface = Mass/sigma
      self.inertia = Inertia
      self.center = Center_Mass
      self.l_border_x = L_border_x
      self.l_border_y = L_border_y
      self.l_border = L_border

    #-------------------------------------------------------------------------------

    def P_is_inside(self,P):
      '''Determine if a point P is inside of a grain

      Make a slide on constant y. Every time a border is crossed, the point switches between in and out.
      see Franklin 1994, see Alonso-Marroquin 2009

          Input :
              itself (a grain)
              a point (a 2 x 1 numpy array)
          Output :
              True or False, depending on the fact that the point is inside the grain or not (a bool)
      '''
      counter = 0
      for i_p_border in range(len(self.l_border)-1):
          #consider only points if the coordinates frame the y-coordinate of the point
          if (self.l_border[i_p_border][1]-P[1])*(self.l_border[i_p_border+1][1]-P[1]) < 0 :
            x_border = self.l_border[i_p_border][0] + (self.l_border[i_p_border+1][0]-self.l_border[i_p_border][0])*(P[1]-self.l_border[i_p_border][1])/(self.l_border[i_p_border+1][1]-self.l_border[i_p_border][1])
            if x_border > P[0] :
                counter = counter + 1
      if counter % 2 == 0:
        return False
      else :
        return True

    #-------------------------------------------------------------------------------

    def Compute_sphericity(self, dict_algorithm):
      '''Compute sphericity of the particle with five parameters.

      The parameters used are the area, the diameter, the circle ratio, the perimeter and the width to length ratio sphericity.
      See Zheng, J., Hryciw, R.D. (2015) Traditional soil particle sphericity, roundness and surface roughness by computational geometry, Geotechnique, Vol 65

          Input :
              itself (a grain)
              an algorithm dictionnary (a dict)
          Output :
              Nothing, but the grain gets updated attributes (five floats)
      '''
      #Find the minimum circumscribing circle
      #look for the two farthest and nearest points
      MaxDistance = 0
      for i_p in range(0,len(self.l_border)-2):
          for j_p in range(i_p+1,len(self.l_border)-1):
              Distance = np.linalg.norm(self.l_border[i_p]-self.l_border[j_p])
              if Distance > MaxDistance :
                  ij_farthest = (i_p,j_p)
                  MaxDistance = Distance

      #Trial circle
      center_circumscribing = (self.l_border[ij_farthest[0]]+self.l_border[ij_farthest[1]])/2
      radius_circumscribing = MaxDistance/2
      Circumscribing_Found = True
      Max_outside_distance = radius_circumscribing
      for i_p in range(len(self.l_border)-1):
          #there is a margin here because of the numerical approximation
          if np.linalg.norm(self.l_border[i_p]-center_circumscribing) > (1+dict_algorithm['sphericity_margin'])*radius_circumscribing and i_p not in ij_farthest: #vertex outside the trial circle
            Circumscribing_Found = False
            if np.linalg.norm(self.l_border[i_p]-center_circumscribing) > Max_outside_distance:
                k_outside_farthest = i_p
                Max_outside_distance = np.linalg.norm(self.l_border[i_p]-center_circumscribing)
      #The trial guess does not work
      if not Circumscribing_Found:
          L_ijk_circumscribing = [ij_farthest[0],ij_farthest[1],k_outside_farthest]
          center_circumscribing, radius_circumscribing = FindCircleFromThreePoints(self.l_border[L_ijk_circumscribing[0]],self.l_border[L_ijk_circumscribing[1]],self.l_border[L_ijk_circumscribing[2]])
          Circumscribing_Found = True
          for i_p in range(len(self.l_border)-1):
              #there is a margin here because of the numerical approximation
              if np.linalg.norm(self.l_border[i_p]-center_circumscribing) > (1+dict_algorithm['sphericity_margin'])*radius_circumscribing and i_p not in L_ijk_circumscribing: #vertex outside the circle computed
                Circumscribing_Found = False
          #see article for other case
          if not Circumscribing_Found:
              print('This algorithm is not developped for this case...')

      #look for length and width
      length = MaxDistance
      u_maxDistance = (self.l_border[ij_farthest[0]]-self.l_border[ij_farthest[1]])/np.linalg.norm(self.l_border[ij_farthest[0]]-self.l_border[ij_farthest[1]])
      v_maxDistance = np.array([u_maxDistance[1], -u_maxDistance[0]])
      MaxWidth = 0
      for i_p in range(0,len(self.l_border)-2):
        for j_p in range(i_p+1,len(self.l_border)-1):
            Distance = abs(np.dot(self.l_border[i_p]-self.l_border[j_p],v_maxDistance))
            if Distance > MaxWidth :
                ij_width = (i_p,j_p)
                MaxWidth = Distance
      width = MaxWidth

      #look for maximum inscribed circle
      #discretisation of the grain
      l_x_inscribing = np.linspace(min(self.l_border_x),max(self.l_border_x),dict_algorithm['n_spatial_inscribing'])
      l_y_inscribing = np.linspace(min(self.l_border_y),max(self.l_border_y),dict_algorithm['n_spatial_inscribing'])
      #creation of an Euclidean distance map to the nearest boundary vertex
      map_inscribing = np.zeros((dict_algorithm['n_spatial_inscribing'],dict_algorithm['n_spatial_inscribing']))
      #compute the map
      for i_x in range(dict_algorithm['n_spatial_inscribing']):
          for i_y in range(dict_algorithm['n_spatial_inscribing']):
              p = np.array([l_x_inscribing[i_x], l_y_inscribing[-1-i_y]])
              #work only if the point is inside the grain
              if self.P_is_inside(p):
                  #look for the nearest vertex
                  MinDistance = None
                  for q in self.l_border[:-1]:
                      Distance = np.linalg.norm(p-q)
                      if MinDistance == None or Distance < MinDistance:
                          MinDistance = Distance
                  map_inscribing[-1-i_y][i_x] = MinDistance
              else :
                  map_inscribing[-1-i_y][i_x] = 0
      #look for the peak of the map
      index_max = np.argmax(map_inscribing)
      l = index_max//dict_algorithm['n_spatial_inscribing']
      c = index_max%dict_algorithm['n_spatial_inscribing']
      radius_inscribing = map_inscribing[l][c]

      #Area Sphericity
      SurfaceParticle = self.surface
      SurfaceCircumscribing = math.pi*radius_circumscribing**2
      AreaSphericity = SurfaceParticle / SurfaceCircumscribing
      if Circumscribing_Found: #else, same value
          self.area_sphericity = AreaSphericity

      #Diameter Sphericity
      DiameterSameAreaParticle = 2*math.sqrt(self.surface/math.pi)
      DiameterCircumscribing = radius_circumscribing*2
      DiameterSphericity = DiameterSameAreaParticle / DiameterCircumscribing
      if Circumscribing_Found: #else, same value
          self.diameter_sphericity = DiameterSphericity

      #Circle Ratio Sphericity
      DiameterInscribing = radius_inscribing*2
      CircleRatioSphericity = DiameterInscribing / DiameterCircumscribing
      self.circle_ratio_sphericity = CircleRatioSphericity

      #Perimeter Sphericity
      PerimeterSameAreaParticle = 2*math.sqrt(self.surface*math.pi)
      PerimeterParticle = 0
      for i in range(len(self.l_border)-1):
          PerimeterParticle = PerimeterParticle + np.linalg.norm(self.l_border[i+1]-self.l_border[i])
      PerimeterSphericity = PerimeterSameAreaParticle / PerimeterParticle
      self.perimeter_sphericity = PerimeterSphericity

      #Width to length ratio Spericity
      WidthToLengthRatioSpericity = width / length
      self.width_to_length_ratio_sphericity = WidthToLengthRatioSpericity

    #---------------------------------------------------------------------------

    def PFtoDEM_Multi(self,FileToRead,dict_algorithm,dict_sample):
        '''
        Read file from MOOSE simulation to reconstruct the phase field of the grain.

            Input :
                itself (a grain)
                the name of the file to read (a string)
                an algorithm dictionnary (a dictionnary)
                a sample dictionnary (a dictionnary)
            Output :
                Nothing but the grain gets an updated attribute (a n_y x n_x numpy array)
        '''
        #--------------------------------------------------------------------------
        #Global parameters
        #---------------------------------------------------------------------------

        self.etai_M = np.array(np.zeros((len(dict_sample['y_L']),len(dict_sample['x_L']))))

        id_L = None
        eta_selector_len = len('        <DataArray type="Float64" Name="etai')
        end_len = len('        </DataArray>')
        XYZ_selector_len = len('        <DataArray type="Float64" Name="Points"')
        data_jump_len = len('          ')

        for i_proc in range(dict_algorithm['np_proc']):

            L_Work = [[], #X
                      [], #Y
                      []] #etai

        #---------------------------------------------------------------------------
        #Reading file
        #---------------------------------------------------------------------------

            f = open(f'{FileToRead}_{i_proc}.vtu','r')
            data = f.read()
            f.close
            lines = data.splitlines()

            #iterations on line
            for line in lines:

                if line[0:eta_selector_len] == '        <DataArray type="Float64" Name="eta'+str(self.id):
                    id_L = 2

                elif line[0:XYZ_selector_len] == '        <DataArray type="Float64" Name="Points"':
                    id_L = 0

                elif (line[0:end_len] == '        </DataArray>' or  line[0:len('          <InformationKey')] == '          <InformationKey') and id_L != None:
                    id_L = None

                elif line[0:data_jump_len] == '          ' and id_L == 2: #Read etai
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

            #Adaptating data and update of etai_M
            for i in range(len(L_Work[0])):
                #Interpolation method
                L_dy = []
                for y_i in dict_sample['y_L'] :
                    L_dy.append(abs(y_i - L_Work[1][i]))
                L_dx = []
                for x_i in dict_sample['x_L'] :
                    L_dx.append(abs(x_i - L_Work[0][i]))
                self.etai_M[-1-list(L_dy).index(min(L_dy))][list(L_dx).index(min(L_dx))] = L_Work[2][i]

    #---------------------------------------------------------------------------

    def ExtractPF_from_Eta(self, L_etai_M, dict_algorithm, dict_material, dict_sample):
        '''
        Extract from the total phase field the variable associated with the grain.

            Input :
                itself (a grain)
                a list of phase field (a neta x nx x ny numpy array)
                an algorithm dictionnary (a dict)
                a material dictionnary (a dict)
                a sample dictionnary (a dict)
            Output :
                Nothing but the grain gets an updated attribute (a n_y x n_x numpy array)
        '''
        etai_M = L_etai_M[self.etai].copy()
        #extract a spatial zone
        x_min = self.center[0] - 1.2*self.r_max
        x_max = max(self.l_border_x)+dict_material['w']
        x_max = self.center[0] + 1.2*self.r_max
        y_min = min(self.l_border_y)-dict_material['w']
        y_min = self.center[1] - 1.2*self.r_max
        y_max = max(self.l_border_y)+dict_material['w']
        y_max = self.center[1] + 1.2*self.r_max

        #look for this part inside the global mesh
        #create search list
        x_L_search_min = abs(np.array(dict_sample['x_L'])-x_min)
        x_L_search_max = abs(np.array(dict_sample['x_L'])-x_max)
        y_L_search_min = abs(np.array(dict_sample['y_L'])-y_min)
        y_L_search_max = abs(np.array(dict_sample['y_L'])-y_max)

        #get index
        i_x_min = list(x_L_search_min).index(min(x_L_search_min))
        i_x_max = list(x_L_search_max).index(min(x_L_search_max))
        i_y_min = list(y_L_search_min).index(min(y_L_search_min))
        i_y_max = list(y_L_search_max).index(min(y_L_search_max))

        #Filter
        filter = np.array(np.zeros((len(dict_sample['y_L']),len(dict_sample['x_L']))), dtype = bool)
        for l in range(i_y_min,i_y_max+1):
            for c in range(i_x_min,i_x_max+1):
                if self.etai_M[-1-l][c] > 0.1:
                    filter[-1-l][c] = True
        struc_element = np.array(np.ones(5,5),dtype = bool)
        filter = binary_dilation(filter, struc_element)

        self.etai_M = np.array(np.zeros((len(dict_sample['y_L']),len(dict_sample['x_L']))))
        for l in range(i_y_min,i_y_max+1):
            for c in range(i_x_min,i_x_max+1):
                if filter[-1-l][c] :
                    self.etai_M[-1-l][c] = etai_M[-1-l][c]

    #---------------------------------------------------------------------------

    def move_grain_interpolation(self, U, dtheta, dict_material, dict_sample):
        '''
        Move the grain by updating the phase field of the grain.

        A bilinear interpolation on the phase field is done. See https://en.wikipedia.org/wiki/Bilinear_interpolation

        The mass conservation is better than with move_grain_rebuild().

            Input :
                itself (a grain)
                the displacement asked (a 2 x 1 numpy array)
                the rotation asked (a float)
                a material dictionnary (a dict)
                a sample dictionnary (a dictionnary)
            Output :
                Nothing but the grain gets an updated attribute (a n_y x n_x numpy array)
        '''
        #initialization
        etai_M_init = self.etai_M.copy()
        sum_eta_before = 0

        #extract a spatial zone (#consider the rotation ?)
        if U[0] >= 0 :
            x_min = min(self.l_border_x)-dict_material['w'] - U[0]
            x_max = max(self.l_border_x)+dict_material['w']
        else :
            x_min = min(self.l_border_x)-dict_material['w']
            x_max = max(self.l_border_x)+dict_material['w'] - U[0]
        if U[1] >= 0:
            y_min = min(self.l_border_y)-dict_material['w'] - U[1]
            y_max = max(self.l_border_y)+dict_material['w']
        else:
            y_min = min(self.l_border_y)-dict_material['w']
            y_max = max(self.l_border_y)+dict_material['w'] - U[1]
        #look for this part inside the global mesh
        #create search list
        x_L_search_min = abs(np.array(dict_sample['x_L'])-x_min)
        x_L_search_max = abs(np.array(dict_sample['x_L'])-x_max)
        y_L_search_min = abs(np.array(dict_sample['y_L'])-y_min)
        y_L_search_max = abs(np.array(dict_sample['y_L'])-y_max)
        #get index
        i_x_min = list(x_L_search_min).index(min(x_L_search_min))
        i_x_max = list(x_L_search_max).index(min(x_L_search_max))
        i_y_min = list(y_L_search_min).index(min(y_L_search_min))
        i_y_max = list(y_L_search_max).index(min(y_L_search_max))

        #iteration on the mesh
        for l in range(i_y_min, i_y_max+1):
            for c in range(i_x_min, i_x_max+1):
                sum_eta_before = sum_eta_before + etai_M_init[-1-l][c]

                #Interpolation of the rotation
                p = np.array([dict_sample['x_L'][c], dict_sample['y_L'][l]])
                pp = p - (self.center - U) #take center before the rigib body motion
                M_rot = np.array([[math.cos(dtheta), -math.sin(dtheta)],
                                  [math.sin(dtheta),  math.cos(dtheta)]])
                pp = np.dot(M_rot, pp)
                p = pp + self.center - U #take center before the rigib body motion
                #check the new position in the mesh
                if p[0] <= min(dict_sample['x_L']) or max(dict_sample['x_L']) <= p[0] or p[1] <= min(dict_sample['y_L']) or max(dict_sample['y_L']) <= p[1] :
                    self.etai_M[-1-l][c] = 0 #no information because out of the mesh, it is 0
                else : #look for the nearest nodes
                    #look on x axis
                    i_x = 1
                    while not (dict_sample['x_L'][i_x-1] <= p[0] and p[0] < dict_sample['x_L'][i_x])  :
                        i_x = i_x + 1
                    #look on y axis
                    i_y = 1
                    while not (dict_sample['y_L'][i_y-1] <= p[1] and p[1] < dict_sample['y_L'][i_y])  :
                        i_y = i_y + 1
                    #definition of the nearest nodes
                    p1  = np.array([dict_sample['x_L'][i_x-1], dict_sample['y_L'][i_y-1]])
                    p2  = np.array([dict_sample['x_L'][i_x]  , dict_sample['y_L'][i_y-1]])
                    p12 = np.array([p[0]                     , dict_sample['y_L'][i_y-1]])
                    p3  = np.array([dict_sample['x_L'][i_x-1], dict_sample['y_L'][i_y]  ])
                    p4  = np.array([dict_sample['x_L'][i_x]  , dict_sample['y_L'][i_y]  ])
                    p34 = np.array([p[0]                     , dict_sample['y_L'][i_y]  ])
                    #definition of value at those nearest nodes
                    q1  = etai_M_init[-(i_y-1)-1][i_x-1]
                    q2  = etai_M_init[-(i_y-1)-1][i_x]
                    q3  = etai_M_init[-i_y-1][i_x-1]
                    q4  = etai_M_init[-i_y-1][i_x]
                    #first interpolations, compute intermediate nodes
                    q12 = (q1*(p2[0]-p[0]) + q2*(p[0]-p1[0]))/(p2[0]-p1[0])
                    q34 = (q3*(p4[0]-p[0]) + q4*(p[0]-p3[0]))/(p4[0]-p1[0])
                    #interpolation of p
                    q = (q12*(p[1]-p34[1]) + q34*(p12[1]-p[1]))/(p12[1]-p34[1])
                    #update etai_M
                    self.etai_M[-1-l][c] = q

        #reinitialization
        etai_M_init = self.etai_M.copy()
        sum_eta_after = 0

        #extract a spatial zone (#consider the rotation ?)
        if U[0] >= 0 :
            x_min = min(self.l_border_x)-dict_material['w'] - U[0]
            x_max = max(self.l_border_x)+dict_material['w']
        else :
            x_min = min(self.l_border_x)-dict_material['w']
            x_max = max(self.l_border_x)+dict_material['w'] - U[0]
        if U[1] >= 0:
            y_min = min(self.l_border_y)-dict_material['w'] - U[1]
            y_max = max(self.l_border_y)+dict_material['w']
        else:
            y_min = min(self.l_border_y)-dict_material['w']
            y_max = max(self.l_border_y)+dict_material['w'] - U[1]
        #look for this part inside the global mesh
        #create search list
        x_L_search_min = abs(np.array(dict_sample['x_L'])-x_min)
        x_L_search_max = abs(np.array(dict_sample['x_L'])-x_max)
        y_L_search_min = abs(np.array(dict_sample['y_L'])-y_min)
        y_L_search_max = abs(np.array(dict_sample['y_L'])-y_max)
        #get index
        i_x_min = list(x_L_search_min).index(min(x_L_search_min))
        i_x_max = list(x_L_search_max).index(min(x_L_search_max))
        i_y_min = list(y_L_search_min).index(min(y_L_search_min))
        i_y_max = list(y_L_search_max).index(min(y_L_search_max))

        #iteration on the mesh
        for l in range(i_y_min, i_y_max+1):
            for c in range(i_x_min, i_x_max+1):
                #Interpolation of the translation
                p = np.array([dict_sample['x_L'][c], dict_sample['y_L'][l]])
                p = p - U
                #check the new position in the mesh
                if p[0] <= min(dict_sample['x_L']) or max(dict_sample['x_L']) <= p[0] or p[1] <= min(dict_sample['y_L']) or max(dict_sample['y_L']) <= p[1] :
                    self.etai_M[-1-l][c] = 0 #no information because out of the mesh, it is 0
                else : #look for the nearest nodes
                    #look on x axis
                    i_x = 1
                    while not (dict_sample['x_L'][i_x-1] <= p[0] and p[0] < dict_sample['x_L'][i_x])  :
                        i_x = i_x + 1
                    #look on y axis
                    i_y = 1
                    while not (dict_sample['y_L'][i_y-1] <= p[1] and p[1] < dict_sample['y_L'][i_y])  :
                        i_y = i_y + 1
                    #definition of the nearest nodes
                    p1  = np.array([dict_sample['x_L'][i_x-1], dict_sample['y_L'][i_y-1]])
                    p2  = np.array([dict_sample['x_L'][i_x]  , dict_sample['y_L'][i_y-1]])
                    p12 = np.array([p[0]                     , dict_sample['y_L'][i_y-1]])
                    p3  = np.array([dict_sample['x_L'][i_x-1], dict_sample['y_L'][i_y]  ])
                    p4  = np.array([dict_sample['x_L'][i_x]  , dict_sample['y_L'][i_y]  ])
                    p34 = np.array([p[0]                     , dict_sample['y_L'][i_y]  ])
                    #definition of value at those nearest nodes
                    q1  = etai_M_init[-(i_y-1)-1][i_x-1]
                    q2  = etai_M_init[-(i_y-1)-1][i_x]
                    q3  = etai_M_init[-i_y-1][i_x-1]
                    q4  = etai_M_init[-i_y-1][i_x]
                    #first interpolations, compute intermediate nodes
                    q12 = (q1*(p2[0]-p[0]) + q2*(p[0]-p1[0]))/(p2[0]-p1[0])
                    q34 = (q3*(p4[0]-p[0]) + q4*(p[0]-p3[0]))/(p4[0]-p1[0])
                    #interpolation of p
                    q = (q12*(p[1]-p34[1]) + q34*(p12[1]-p[1]))/(p12[1]-p34[1])
                    #update etai_M
                    self.etai_M[-1-l][c] = q
                    sum_eta_after = sum_eta_after + q

        self.delta_sum_eta = (sum_eta_before - sum_eta_after)/sum_eta_before*100

    #---------------------------------------------------------------------------

    def move_grain_rebuild(self,dict_material,dict_sample):
        '''
        Move the grain by updating the phase field of the grain.

        A rebuild on the phase field is done. The mass conservation is verified by producing solute.

            Input :
                itself (a grain)
                a material dictionnary (a dict)
                a sample dictionnary (a dictionnary)
            Output :
                Nothing but the grain gets an updated attribute (a n_y x n_x numpy array)
        '''
        #save previous phase field
        save_etai_M = self.etai_M.copy()
        #compute the new phase field
        self.build_etai_M(dict_material,dict_sample)
        #compute the delta etai
        sum_eta_before = 0
        sum_eta_after = 0
        for l in range(len(dict_sample['y_L'])):
            for c in range(len(dict_sample['x_L'])):
                sum_eta_before = sum_eta_before + save_etai_M[l][c]
                sum_eta_after = sum_eta_after + self.etai_M[l][c]
        self.delta_sum_eta = (sum_eta_before - sum_eta_after)/sum_eta_before*100

    #-------------------------------------------------------------------------------

    def update_geometry_kinetic(self, V, A, W, DT):
        """
        Update the acceleration and the velocity of a grain. Update geometrical parameters as border and center nodes.

            Input :
                itself (a grain)
                a speed (a 1 x 2 numpy array)
                an acceleration (a 1 x 2 numpy array)
                an angular speed (a float)
                a time step (a float)
            Ouput :
                Nothing, but the position of the grain is updated
        """
        #translation
        #max_speed
        if np.linalg.norm(V) > self.r_mean / (DT*1000) :
            V = V * self.r_mean / (DT*1000) / np.linalg.norm(V)
        self.v = V
        self.a = A
        for i in range(len(self.l_border)):
            self.l_border[i] = self.l_border[i] + self.v*DT
            self.l_border_x[i] = self.l_border_x[i] + self.v[0]*DT
            self.l_border_y[i] = self.l_border_y[i] + self.v[1]*DT
        self.center = self.center + self.v*DT

        #rotation
        self.w = W
        self.theta = self.theta + self.w*DT

        for i_theta_r in range(len(self.l_theta_r)) :
            theta_r = self.l_theta_r[i_theta_r]
            theta_r = theta_r + self.w*DT
            while theta_r >= 2*math.pi:
                theta_r = theta_r - 2*math.pi
            while theta_r < 0 :
                theta_r = theta_r + 2*math.pi
            self.l_theta_r[i_theta_r] = theta_r

        for i in range(len(self.l_border)):
            p = self.l_border[i] - self.center
            Rot_Matrix = np.array([[math.cos(self.w*DT), -math.sin(self.w*DT)],
                                   [math.sin(self.w*DT),  math.cos(self.w*DT)]])
            p = np.dot(Rot_Matrix,p)
            self.l_border[i] = p + self.center
            self.l_border_x[i] = p[0] + self.center[0]
            self.l_border_y[i] = p[1] + self.center[1]

    #-------------------------------------------------------------------------------

    def init_f_control(self,dict_sollicitation):
        """
        Initialize the force applied to the grain.

        A gravity of g is applied.

            Input :
                itself (a grain)
                a sollicitations dictionnary (a dict)
            Ouput :
                Nothing, but the force applied on the grain is initialized
        """
        self.fx = 0
        self.fy = -dict_sollicitation['gravity']*self.mass
        self.f = np.array([self.fx,self.fy])
        self.mz = 0

    #-------------------------------------------------------------------------------

    def update_f(self, Fx, Fy, p_application):
        """
        Add a force to the grain.

            Input :
                itself (a grain)
                the value x and y of the force (two float)
                an applicaiton point (a 1 x 2 numpy array)
            Output :
                Nothing, but a force is applied to the grain
        """
        self.fx = self.fx + Fx
        self.fy = self.fy + Fy
        self.f = np.array([self.fx,self.fy])

        v1 = np.array([p_application[0]-self.center[0], p_application[1]-self.center[1], 0])
        v2 = np.array([Fx, Fy, 0])
        self.mz = self.mz + np.cross(v1,v2)[2]

#-------------------------------------------------------------------------------
#Functions
#-------------------------------------------------------------------------------

def Compute_overlap_2_grains(dict_sample):
    '''
    Compute the current overlap between two grains.

    It is assumed the sample is composed  by only 2 grains.

        Input :
            a sample dictionnary (a dictionnary)
        Output :
            Nothing but the sample dictionnary gets an updated value (a float)
    '''
    #2 grains
    g1 = dict_sample['L_g'][0]
    g2 = dict_sample['L_g'][1]

    #compute overlap
    #assume the normal n12 is +x axis
    overlap = max(g1.l_border_x) - min(g2.l_border_x)

    #Add element in dict
    dict_sample['overlap'] = overlap

#-------------------------------------------------------------------------------

def FindCircleFromThreePoints(P1,P2,P3):
    '''
    Compute the circumscribing circle of a triangle defined by three points.

    https://www.geeksforgeeks.org/program-find-circumcenter-triangle-2/

        Input :
            three points (a 2 x 1 numpy array)
        Output :
            a center (a 2 x 1 numpy array)
            a radius (a float)
    '''
    # Line P1P2 is represented as ax + by = c and line P2P3 is represented as ex + fy = g
    a, b, c = lineFromPoints(P1, P2)
    e, f, g = lineFromPoints(P2, P3)

    # Converting lines P1P2 and P2P3 to perpendicular bisectors.
    #After this, L : ax + by = c and M : ex + fy = g
    a, b, c = perpendicularBisectorFromLine(P1, P2, a, b, c)
    e, f, g = perpendicularBisectorFromLine(P2, P3, e, f, g)

    # The point of intersection of L and M gives the circumcenter
    circumcenter = lineLineIntersection(a, b, c, e, f, g)

    if np.linalg.norm(circumcenter - np.array([10**9,10**9])) == 0:
        raise ValueError('The given points do not form a triangle and are collinear...')
    else :
        #compute the radius
        radius = max([np.linalg.norm(P1-circumcenter), np.linalg.norm(P2-circumcenter), np.linalg.norm(P3-circumcenter)])

    return circumcenter, radius

#-------------------------------------------------------------------------------

def lineFromPoints(P, Q):
    '''
    Function to find the line given two points

    Used in FindCircleFromThreePoints().
    The equation is c = ax + by.
    https://www.geeksforgeeks.org/program-find-circumcenter-triangle-2/

        Input :
            two points (a 2 x 1 numpy array)
        Output :
            three characteristic of the line (three floats)
    '''
    a = Q[1] - P[1]
    b = P[0] - Q[0]
    c = a * (P[0]) + b * (P[1])
    return a, b, c

#-------------------------------------------------------------------------------

def perpendicularBisectorFromLine(P, Q, a, b, c):
    '''
    Function which converts the input line to its perpendicular bisector.

    Used in FindCircleFromThreePoints().
    The equation is c = ax + by.
    https://www.geeksforgeeks.org/program-find-circumcenter-triangle-2/

        Input :
            two points (a 2 x 1 numpy array)
            three characteristic of the line (three floats)
        Output :
            three characteristic of the perpendicular bisector (three floats)
    '''
    mid_point = [(P[0] + Q[0])//2, (P[1] + Q[1])//2]
    # c = -bx + ay
    c = -b * (mid_point[0]) + a * (mid_point[1])
    temp = a
    a = -b
    b = temp
    return a, b, c

#-------------------------------------------------------------------------------

def lineLineIntersection(a1, b1, c1, a2, b2, c2):
    '''
    Returns the intersection point of two lines.

    Used in FindCircleFromThreePoints().
    https://www.geeksforgeeks.org/program-find-circumcenter-triangle-2/

        Input :
            six characteristics of the line 1 and 2 (six floats)
        Output :
            the intersection point (a 2 x 1 numpy array)
    '''
    determinant = a1 * b2 - a2 * b1
    if (determinant == 0):
        # The lines are parallel.
        return np.array([10**9,10**9])
    else:
        x = (b2 * c1 - b1 * c2)//determinant
        y = (a1 * c2 - a2 * c1)//determinant
        return np.array([x, y])
