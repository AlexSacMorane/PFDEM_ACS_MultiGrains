# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This file contains ??.nt functions used in the simulation.
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#Function
#-------------------------------------------------------------------------------

def LG_tempo(dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitations, simulation_report):
    """
    Create an initial condition

        Input :
            an algorithm dictionnary (a dict)
            a geometry dictionnary (a dict)
            an initial condition dictionnary (a dict)
            a material dictionnary (a dict)
            a sample dictionnary (a dict)
            a sollicitations dictionnary (a dict)
            a simulation report (a report)
        Output :
            Nothing, but dictionnaries are updated
    """
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    #load data needed
    n_generation = dict_ic['n_generation']
    if n_generation != 2:
        simulation_report.write('n_generation must be equal to 2 !')
        raise ValueError('n_generation must be equal to 2 !')
    factor = dict_ic['factor_ymax_box']
    N_grain_disk = dict_geometry['N_grain_disk']
    N_grain_square = dict_geometry['N_grain_square']
    N_grain = dict_geometry['N_grain_disk'] + dict_geometry['N_grain_square']
    L_radius = dict_geometry['L_R']
    L_percentage_radius = dict_geometry['L_percentage_R']
    L_dimension = dict_geometry['L_Dimension']
    L_percentage_dimension = dict_geometry['L_percentage_Dimension']
    x_min = dict_sample['x_box_min']
    x_max = dict_sample['x_box_max']
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

    #define the y_max for the grains generation
    radius_mean = 0
    for i in range(len(L_radius)):
        radius_mean = radius_mean + L_radius[i]*L_percentage_radius[i]
    dimension_mean = 0
    for i in range(len(L_dimension)):
        dimension_mean = dimension_mean + L_dimension[i]*L_percentage_dimension[i]
    dy_creation = N_grain_disk/n_generation*factor*(2*radius_mean)**2/(x_max-x_min) + N_grain_square/n_generation*factor*(dimension_mean)**2/(x_max-x_min)

    #plan the grains generation
    L_n_grain_radius_try_one = []
    L_n_grain_radius = []
    L_n_grain_radius_done = []
    for percentage in L_percentage_radius:
        L_n_grain_radius_try_one.append(int(N_grain_disk*percentage/n_generation))
        L_n_grain_radius.append(int(N_grain_disk*percentage))
        L_n_grain_radius_done.append(0)
    L_n_grain_dimension_try_one = []
    L_n_grain_dimension = []
    L_n_grain_dimension_done = []
    for percentage in L_percentage_dimension:
        L_n_grain_dimension_try_one.append(int(N_grain_square*percentage/n_generation))
        L_n_grain_dimension.append(int(N_grain_square*percentage))
        L_n_grain_dimension_done.append(0)

    #Creation of grains
    #grains generation is decomposed in several steps (creation of grain then settlement)
    i_DEM = 0
    L_L_g_tempo = []

    #---------------------------------------------------------------------------

    print('First generation of grains')
    L_g_tempo = []

    #add elements in dicts
    dict_ic['L_g_tempo'] = L_g_tempo
    dict_ic['L_L_g_tempo'] = L_L_g_tempo
    dict_ic['i_DEM_IC'] = i_DEM
    dict_ic['L_n_grain_radius_try_one'] = L_n_grain_radius_try_one
    dict_ic['L_n_grain_radius'] = L_n_grain_radius
    dict_ic['L_n_grain_radius_done'] = L_n_grain_radius_done
    dict_ic['L_n_grain_dimension_try_one'] = L_n_grain_dimension_try_one
    dict_ic['L_n_grain_dimension'] = L_n_grain_dimension
    dict_ic['L_n_grain_dimension_done'] = L_n_grain_dimension_done
    dict_ic['last_id'] = 0
    dict_sample['y_box_min_ic'] = dict_sample['y_box_min']
    dict_sample['dy_creation'] = dy_creation

    Create_grains(dict_ic, dict_geometry, dict_sample, dict_material, 1, simulation_report)

    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    #load data needed
    L_g_tempo = dict_ic['L_g_tempo']
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

    #DEM to find the steady-state configuration after loading
    #find the maximum y (center+radius)
    y_max = dict_sample['y_box_min_ic']
    for grain in L_g_tempo:
        if max(grain.l_border_y) > y_max:
            y_max = max(grain.l_border_y)

    #add element in dict
    dict_sample['y_box_max'] = y_max

    DEM_loading(dict_ic, dict_material, dict_sample, dict_sollicitations, False, simulation_report)

    #---------------------------------------------------------------------------

    print('Second generation of grains')

    #Update dict
    L_g_tempo = []
    dict_ic['L_g_tempo'] = L_g_tempo

    Create_grains(dict_ic, dict_geometry, dict_sample, dict_material, 2, simulation_report)

    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    #load data needed
    L_g_tempo = dict_ic['L_g_tempo']
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

    #DEM to find the steady-state configuration after loading
    #find the maximum y (center+radius)
    y_max = dict_sample['y_box_min_ic']
    for grain in L_g_tempo:
        if max(grain.l_border_y) > y_max:
            y_max = max(grain.l_border_y)

    DEM_loading(dict_ic, dict_material, dict_sample, dict_sollicitations, False, simulation_report)

    #---------------------------------------------------------------------------

    print('Combine generations of grains')

    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    #load data needed
    L_L_g_tempo = dict_ic['L_L_g_tempo']
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

    L_g = []
    for L_g_tempo in L_L_g_tempo:
        for g_tempo in L_g_tempo:
            L_g.append(g_tempo)

    #update element in dict
    dict_ic['L_g_tempo'] = L_g

    DEM_loading(dict_ic, dict_material, dict_sample, dict_sollicitations, True, simulation_report)

    #update element in dict
    dict_sample['y_box_max'] = dict_ic['y_box_max']

    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    #load data needed
    L_g_tempo = dict_ic['L_g_tempo']
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

    simulation_report.write_and_print(str(len(L_g_tempo))+' / '+str(N_grain)+' grains have been created\n','\n'+str(len(L_g_tempo))+' / '+str(N_grain)+' grains have been created\n')

#-------------------------------------------------------------------------------

def DEM_loading(dict_ic, dict_material, dict_sample, dict_sollicitations, multi_generation, simulation_report):
    """
    Loading the granular system.

        Input :
            an initial condition dictionnary (a dict)
            a material dictionnary (a dict)
            a smaple dictionnary (a dict)
            a sollicitations dictionnary (a dict)
            a Boolean to determine if multiple generations are considered (a Boolean)
            a simultion report (a report)
        Output :
            Nothing, but initial condition dictionnary is updated
    """
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    #load data needed
    L_g_tempo = dict_ic['L_g_tempo']
    dt_DEM = dict_ic['dt_DEM_IC']
    i_DEM_stop = dict_ic['i_DEM_stop_IC']
    i_DEM = dict_ic['i_DEM_IC']
    Ecin_ratio_IC = dict_ic['Ecin_ratio_IC']
    i_print_plot_IC = dict_ic['i_print_plot_IC']
    factor_neighborhood_IC = dict_ic['factor_neighborhood_IC']
    if multi_generation :
        i_update_neighbouroods = dict_ic['i_update_neighborhoods_com']
        y_min = dict_sample['y_box_min']
    else :
        i_update_neighbouroods = dict_ic['i_update_neighborhoods_gen']
        y_min = dict_sample['y_box_min_ic']
    x_min = dict_sample['x_box_min']
    x_max = dict_sample['x_box_max']
    y_max = dict_sample['y_box_max']
    Forcev_target = dict_sollicitations['Vertical_Confinement_Force']
    gravity = dict_sollicitations['gravity']
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

    i_DEM_0 = i_DEM
    DEM_loop_statut = True

    #Initialisation
    dict_ic['L_contact'] = []
    dict_ic['L_contact_ij'] = []
    dict_ic['L_contact_gw'] = []
    dict_ic['L_contact_gw_ij'] = []
    dict_ic['id_contact'] = 0

    #trackers and stop conditions
    Force_tracker = []
    Force_stop = 0
    Ecin_tracker = []
    Ecin_stop = 0
    Ymax_tracker = []
    Ymax_stop = 0
    for grain in L_g_tempo:
        Force_stop = Force_stop + 0.5*grain.mass*gravity
        Ecin_stop = Ecin_stop + 0.5*grain.mass*(Ecin_ratio_IC*grain.r_max/dt_DEM)**2

    while DEM_loop_statut :

        i_DEM = i_DEM + 1

        #Contact detection
        if (i_DEM-i_DEM_0-1) % i_update_neighbouroods  == 0:
            Update_Neighbouroods(dict_ic)
        Grains_Polyhedral_contact_Neighbouroods(dict_ic,dict_material)
        # Detection of contacts between grain and walls
        if (i_DEM-i_DEM_0-1) % i_update_neighbouroods  == 0:
            wall_neighborhood = Update_wall_Neighborhoods(L_g_tempo,factor_neighborhood_IC,x_min,x_max,y_min,y_max)
        Grains_Polyhedral_Wall_contact_Neighborhood(wall_neighborhood,x_min,x_max,y_min,y_max, dict_ic, dict_material)

        #Sollicitation computation
        for grain in dict_ic['L_g_tempo']:
             grain.init_F_control(gravity)
        for contact in  dict_ic['L_contact']+dict_ic['L_contact_gw']:
            contact.normal()
            contact.tangential(dt_DEM)

        #Move grains
        for grain in dict_ic['L_g_tempo']:
            grain.euler_semi_implicite(dt_DEM,10*Ecin_ratio_IC)

        #check if some grains are outside of the study box
        L_ig_to_delete = []
        for id_grain in range(len(dict_ic['L_g_tempo'])):
            if dict_ic['L_g_tempo'][id_grain].center[0] < x_min :
                L_ig_to_delete.append(id_grain)
            elif dict_ic['L_g_tempo'][id_grain].center[0] > x_max :
                L_ig_to_delete.append(id_grain)
            elif dict_ic['L_g_tempo'][id_grain].center[1] < y_min :
                L_ig_to_delete.append(id_grain)
            elif dict_ic['L_g_tempo'][id_grain].center[1] > y_max :
                L_ig_to_delete.append(id_grain)
        L_ig_to_delete.reverse()
        for id_grain in L_ig_to_delete:
            simulation_report.write_and_print('Grain '+str(dict_ic['L_g_tempo'][id_grain].id)+' has been deleted because it is out of the box\n','Grain '+str(dict_ic['L_g_tempo'][id_grain].id)+' has been deleted because it is out of the box')
            dict_ic['L_g_tempo'].pop(id_grain)

        #Control the y_max to have the pressure target
        y_max, Fv = Control_y_max_NR(y_max,Forcev_target,dict_ic['L_contact_gw'],dict_ic['L_g_tempo'])

        #Tracker
        F = F_total(dict_ic['L_g_tempo'])
        Ecin = E_cin_total(dict_ic['L_g_tempo'])
        Force_tracker.append(F)
        Ecin_tracker.append(Ecin)
        Ymax_tracker.append(y_max)

        if i_DEM % i_print_plot_IC ==0:
            if gravity > 0 :
                print('i_DEM',i_DEM,'and Ecin',int(100*Ecin/Ecin_stop),'% and Force',int(100*F/Force_stop),'% and Confinement',int(100*Fv/Forcev_target),'%')
            else :
                print('i_DEM',i_DEM,'and Ecin',int(100*Ecin/Ecin_stop),'% and Confinement',int(100*Fv/Forcev_target),'%')
            if dict_ic['Debug_DEM']:
                Plot_Config_Loaded(dict_ic['L_g_tempo'],x_min,x_max,y_min,y_max,i_DEM)

        #Check stop conditions for DEM
        if i_DEM >= i_DEM_stop + i_DEM_0:
             DEM_loop_statut = False
        if gravity > 0:
            if Ecin < Ecin_stop and F < Force_stop and (0.95*Forcev_target<Fv and Fv<1.05*Forcev_target):
                  DEM_loop_statut = False
        else:
            if Ecin < Ecin_stop and i_DEM >= i_DEM_stop*0.1 + i_DEM_0 and (0.95*Forcev_target<Fv and Fv<1.05*Forcev_target):
                DEM_loop_statut = False
        if dict_ic['L_g_tempo'] == []:
            DEM_loop_statut = False

    #Update dict
    dict_ic['L_g_tempo'] = L_g_tempo
    L_L_g_tempo = dict_ic['L_L_g_tempo']
    L_L_g_tempo.append(L_g_tempo.copy())
    dict_ic['L_L_g_tempo'] = L_L_g_tempo
    dict_ic['y_box_max'] = y_max
    dict_ic['i_DEM_IC'] = i_DEM

#-------------------------------------------------------------------------------

def Create_grains(dict_ic, dict_geometry, dict_sample, dict_material, id_generation, simulation_report):
    """
    Generate the grains.

    A position is tried. This position must not create overlap with already creared temporary grain. If there is no overlap, a temporary grai nis created.

        Input :
            an initial condition dictionnary (a dict)
            a geometry dictionnary (a dict)
            a sample dictionnary (a dict)
            a material dictionnary (a dict)
            a generation id (a int)
            a simulation report (a report)
        Output :
            Nothing, but initial configuration dictionnary is updated
    """
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    #load data needed
    L_g_tempo = dict_ic['L_g_tempo']
    N_test_max = dict_ic['N_test_max']
    if id_generation == 1:
        L_n_grain_radius = dict_ic['L_n_grain_radius_try_one']
        L_n_grain_dimension = dict_ic['L_n_grain_dimension_try_one']
    else :
        L_n_grain_radius = dict_ic['L_n_grain_radius']
        L_n_grain_dimension = dict_ic['L_n_grain_dimension']
    L_n_grain_radius_done = dict_ic['L_n_grain_radius_done']
    L_radius = dict_geometry['L_R']
    L_n_grain_dimension_done = dict_ic['L_n_grain_dimension_done']
    L_dimension = dict_geometry['L_Dimension']
    x_min = dict_sample['x_box_min']
    x_max = dict_sample['x_box_max']
    y_min = dict_sample['y_box_min_ic']
    dy_creation = dict_sample['dy_creation']
    y_max = y_min + dy_creation
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

    #Parameters for the method
    n_not_created = 0

    for i in range(len(L_radius)):
        radius = L_radius[i]
        n_grain = L_n_grain_radius[i]
        n_grain_done = L_n_grain_radius_done[i]
        last_id_grain_created = dict_ic['last_id']
        for id_grain in range(last_id_grain_created, last_id_grain_created + n_grain - n_grain_done):
            i_test = 0
            grain_created = False
            while (not grain_created) and i_test < N_test_max:
                i_test = i_test + 1
                center = np.array([random.uniform(x_min+1.1*radius,x_max-1.1*radius),random.uniform(y_min+1.1*radius,y_max)])
                g_tempo = Grain_Tempo(id_grain-n_not_created,center,radius,dict_material,'Disk')
                grain_created = True
                for grain in L_g_tempo:
                    if Grains_Polyhedral_contact_f(g_tempo,grain):
                        grain_created = False
            if i_test == N_test_max and not grain_created:
                n_not_created = n_not_created + 1
                simulation_report.write_and_print('Grain '+str(id_grain)+' has not been created after '+str(i_test)+' tries\n','Grain '+str(id_grain)+' has not been created after '+str(i_test)+' tries')
            else :
                L_g_tempo.append(g_tempo)
                L_n_grain_radius_done[i] = L_n_grain_radius_done[i] + 1
                dict_ic['last_id'] = dict_ic['last_id'] + 1
    for i in range(len(L_dimension)):
        dimension = L_dimension[i]
        n_grain = L_n_grain_dimension[i]
        n_grain_done = L_n_grain_dimension_done[i]
        last_id_grain_created = dict_ic['last_id']
        for id_grain in range(last_id_grain_created, last_id_grain_created + n_grain - n_grain_done):
            i_test = 0
            grain_created = False
            while (not grain_created) and i_test < N_test_max:
                i_test = i_test + 1
                center = np.array([random.uniform(x_min+1.1*dimension/2,x_max-1.1*dimension/2),random.uniform(y_min+1.1*dimension/2,y_max)])
                g_tempo = Grain_Tempo(id_grain-n_not_created,center,dimension,dict_material,'Square')
                grain_created = True
                for grain in L_g_tempo:
                    if Grains_Polyhedral_contact_f(g_tempo,grain):
                        grain_created = False
            if i_test == N_test_max and not grain_created:
                n_not_created = n_not_created + 1
                simulation_report.write_and_print('Grain '+str(id_grain)+' has not been created after '+str(i_test)+' tries\n','Grain '+str(id_grain)+' has not been created after '+str(i_test)+' tries')
            else :
                L_g_tempo.append(g_tempo)
                L_n_grain_dimension_done[i] = L_n_grain_dimension_done[i] + 1
                dict_ic['last_id'] = dict_ic['last_id'] + 1

    #Update dict
    dict_ic['L_g_tempo'] = L_g_tempo
    dict_ic['L_n_grain_dimension_done'] = L_n_grain_dimension_done
    dict_ic['L_n_grain_radius_done'] = L_n_grain_radius_done

#-------------------------------------------------------------------------------

def Grains_Polyhedral_contact_f(g1,g2):
  """
  Detect the contact grain-grain.

    Input :
        two temporary grains (two grain_tempos)
    Output :
        a Boolean, True if there is contaxct between the twwo grains (a Boolean)
  """
  if np.linalg.norm(g1.center-g2.center) < 1.5*(g1.r_max+g2.r_max):

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

  else:
    return False

#-------------------------------------------------------------------------------

def E_cin_total(L_g):
    """
    Compute total kinetic energy.

        Input :
            a list of temporary grains (a list)
        Output :
            the total kinetic energy (a float)
    """
    Ecin = 0
    for grain in L_g:
        Ecin = Ecin + 1/2*grain.mass*np.dot(grain.v,grain.v)
    return Ecin

#-------------------------------------------------------------------------------

def F_total(L_g):
    """
    Compute total force applied on grains in the sample.

        Input :
            a list of temporary grains (a list)
        Output :
            the total force applied (a float)
    """
    F = 0
    for grain in L_g:
        F = F + np.linalg.norm([grain.fx, grain.fy])
    return F

#-------------------------------------------------------------------------------

def Control_y_max_NR(y_max,Force_target,L_contact_gw,L_g):
    """
    Control the upper wall to apply force.

    A Newton-Raphson method is applied to verify the confinement.
        Input :
            a coordinate of the upper wall (a float)
            a confinement value (a float)
            a list of contact grain - wall (a list)
            a list of temporary grain (a list)
        Output :
            the coordinate of the upper wall (a float)
            a force applied on the upper wall before control (a float)
    """
    F = 0
    overlap_L = []
    k_L = []
    for contact in L_contact_gw:
        if contact.nature == 'gwy_max':
            F = F + contact.Fwg_n
            overlap_L.append(contact.overlap)
            k_L.append(contact.k)
            #compute force applied, save contact overlap and spring

    if overlap_L != []:
        i_NR = 0
        dy = 0
        ite_criteria = True
        #control the upper wall
        if -0.01*Force_target<error_on_ymax_f(dy,overlap_L,k_L,Force_target) and error_on_ymax_f(dy,overlap_L,k_L,Force_target)<0.01*Force_target:
            ite_criteria = False
        while ite_criteria :
            i_NR = i_NR + 1
            dy = dy - error_on_ymax_f(dy,overlap_L,k_L,Force_target)/error_on_ymax_df(dy,overlap_L,k_L)
            if i_NR > 100: #Maximum try
                ite_criteria = False
            if -0.01*Force_target<error_on_ymax_f(dy,overlap_L,k_L,Force_target) and error_on_ymax_f(dy,overlap_L,k_L,Force_target)<0.01*Force_target:
                ite_criteria = False
        y_max = y_max + dy

    else :
        #if there is no contact with the upper wall, the wall is reset
        y_max = Reset_y_max(L_g,Force_target)

    for contact in L_contact_gw:
        if contact.nature == 'gwy_max':
            #reactualisation
            contact.limit = y_max

    return y_max, F

#-------------------------------------------------------------------------------

def error_on_ymax_f(dy,overlap_L,k_L,Force_target) :
    """
    Compute the function f to control the upper wall. It is the difference between the force applied and the target value.

        Input :
            an increment of the upper wall position (a float)
            a list of overlap for contact between temporary grain and upper wall (a list)
            a list of spring for contact between temporary grain and upper wall (a list)
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
            a list of overlap for contact between temporary grain and upper wall (a list)
            a list of spring for contact between temporary grain and upper wall (a list)
        Output :
            the derivative of error_on_ymax_f() (a float)
    """
    df = 0
    for i in range(len(overlap_L)):
        df = df + 3/2*k_L[i]*(max(overlap_L[i]-dy,0))**(1/2)
    return df

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

    factor = 5
    k = factor*4/3*L_g[id_grain_max].y/(1-L_g[id_grain_max].nu*L_g[id_grain_max].nu)*math.sqrt(L_g[id_grain_max].r_max)
    y_max = y_max - (Force/k)**(2/3)

    return y_max

#-------------------------------------------------------------------------------

def Grains_Polyhedral_contact_Neighborhoods_f(g1,g2):
  """
  Detect contact grain - grain.

    Input :
        two temporary grains (two grain_tempo)
    Output :
        a Boolean, True if there is a contact (a Boolean)
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

def Update_Neighbouroods(dict_ic):
    """
    Determine a neighborhood for each grain.

    This function is called every x time step. The contact is determined by Grains_Polyhedral_contact_Neighbouroods().
    Notice that if there is a potential contact between grain_i and grain_j, grain_i is not in the neighborhood of grain_j.
    Whereas grain_j is in the neighbourood of grain_i. With i_grain < j_grain.

        Input :
            an initial condition dictionnary (a dict)
        Output :
            Nothing, but the neighborhood of the temporary grains is updated
    """
    for i_grain in range(len(dict_ic['L_g_tempo'])-1) :
        neighbourood = []
        for j_grain in range(i_grain+1,len(dict_ic['L_g_tempo'])):
            if np.linalg.norm(dict_ic['L_g_tempo'][i_grain].center-dict_ic['L_g_tempo'][j_grain].center) < dict_ic['factor_neighborhood_IC']*(dict_ic['L_g_tempo'][i_grain].r_max+dict_ic['L_g_tempo'][j_grain].r_max):
                neighbourood.append(dict_ic['L_g_tempo'][j_grain])
        dict_ic['L_g_tempo'][i_grain].neighbourood = neighbourood

#-------------------------------------------------------------------------------

def Grains_Polyhedral_contact_Neighbouroods(dict_ic,dict_material):
    """
    Detect contact between a grain and grains from its neighbourood.

    The neighbourood is updated with Update_Neighbouroods().

        Input :
            an initial condition dictionnary (a dict)
            a material dictionnary (a dict)
        Output :
            Nothing, but the initial condition dictionnary is updated with grain - grain contacts
    """
    for i_grain in range(len(dict_ic['L_g_tempo'])-1) :
        grain_i = dict_ic['L_g_tempo'][i_grain]
        for neighbour in dict_ic['L_g_tempo'][i_grain].neighbourood:
            j_grain = neighbour.id
            grain_j = neighbour
            if Grains_Polyhedral_contact_Neighborhoods_f(grain_i,grain_j):
                if (i_grain,j_grain) not in dict_ic['L_contact_ij']:  #contact not detected previously
                   #creation of contact
                   dict_ic['L_contact_ij'].append((i_grain,j_grain))
                   dict_ic['L_contact'].append(Contact_Tempo(dict_ic['id_contact'], grain_i, grain_j, dict_material))
                   dict_ic['id_contact'] = dict_ic['id_contact'] + 1

            else :
                if (i_grain,j_grain) in dict_ic['L_contact_ij'] : #contact detected previously is not anymore
                       dict_ic['L_contact'].pop(dict_ic['L_contact_ij'].index((i_grain,j_grain)))
                       dict_ic['L_contact_ij'].remove((i_grain,j_grain))

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
        if abs(p_x_min-x_min) < factor_neighborhood_IC*grain.r_max :
            wall_neighborhood.append(grain)
        #grain-wall x_max
        if abs(p_x_max-x_max) < factor_neighborhood_IC*grain.r_max :
            wall_neighborhood.append(grain)
        #grain-wall y_min
        if abs(p_y_min-y_min) < factor_neighborhood_IC*grain.r_max :
            wall_neighborhood.append(grain)
        #grain-wall y_max
        if abs(p_y_max-y_max) < factor_neighborhood_IC*grain.r_max :
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

      p_x_min = min(grain.l_border_x)
      p_x_max = max(grain.l_border_x)
      p_y_min = min(grain.l_border_y)
      p_y_max = max(grain.l_border_y)

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

#-------------------------------------------------------------------------------

def Plot_Config_Loaded(L_g,x_min,x_max,y_min,y_max,i):
    """
    Plot loaded configuration.

        Input :
            a list of temporary grain (a list)
            the coordinates of the walls (four floats)
            an iteration (a int)
        Output :
            Nothing, but a .png file is generated (a file)
    """
    plt.figure(1,figsize=(16,9))
    L_x = []
    L_y = []
    L_u = []
    L_v = []
    for grain in L_g:
        plt.plot(grain.l_border_x,grain.l_border_y,'k')
        plt.plot(grain.center[0],grain.center[1],'xk')
        L_x.append(grain.center[0])
        L_y.append(grain.center[1])
        L_u.append(grain.fx)
        L_v.append(grain.fy)
    plt.plot([x_min,x_min,x_max,x_max,x_min],[y_max,y_min,y_min,y_max,y_max],'k')
    plt.axis('equal')
    plt.savefig('Debug/DEM_ite/Init/Config_Loaded_'+str(i)+'.png')
    plt.close(1)

#-------------------------------------------------------------------------------

def Plot_Config_Loaded_End(L_g,x_min,x_max,y_min,y_max):
    """
    Plot loaded configuration at the end of the initial configuration.

        Input :
            a list of temporary grain (a list)
            the coordinates of the walls (four floats)
            an iteration (a int)
        Output :
            Nothing, but a .png file is generated (a file)
    """
    plt.figure(1,figsize=(16,9))
    L_x = []
    L_y = []
    L_u = []
    L_v = []
    for grain in L_g:
        plt.plot(grain.l_border_x,grain.l_border_y,'k')
        plt.plot(grain.center[0],grain.center[1],'xk')
        L_x.append(grain.center[0])
        L_y.append(grain.center[1])
        L_u.append(grain.fx)
        L_v.append(grain.fy)
    plt.plot([x_min,x_min,x_max,x_max,x_min],[y_max,y_min,y_min,y_max,y_max],'k')
    plt.axis('equal')
    plt.savefig('Debug/ConfigLoaded.png')
    plt.close(1)

#-------------------------------------------------------------------------------

def From_LG_tempo_to_usable(dict_ic, dict_geometry, dict_material, dict_sample):
    """
    Create a rea lgrain from a temporary grain.

    The phase variable is built. The distance between the point of the mesh and the particle center determines the value of the variable.
    A cosine profile is applied inside the interface.

        Input :
            an initial condition dictionnary (a dict)
            a geometry dictionnary (a dict)
            a material dictionnary (a dict)
            a sample dictionnary (a dict)
        Output :
            Nothing, but the sample dictionnary is updated with the list of real grains
    """
    L_g = []
    for grain_tempo in dict_ic['L_g_tempo']:

        dict_ic_to_real = {
        'Id' : grain_tempo.id,
        'Type' : grain_tempo.type,
        'Y' : grain_tempo.y,
        'Nu' : grain_tempo.nu,
        'Rho_surf' : grain_tempo.rho_surf,
        'Center' : grain_tempo.center,
        'L_border' : grain_tempo.l_border,
        'L_border_x' : grain_tempo.l_border_x,
        'L_border_y' : grain_tempo.l_border_y,
        'L_r' : grain_tempo.l_r,
        'L_theta_r' : grain_tempo.l_theta_r,
        'R_min' : min(grain_tempo.l_r),
        'R_max' : max(grain_tempo.l_r),
        'R_mean' : np.mean(grain_tempo.l_r),
        'Surface' : grain_tempo.surface,
        'Mass' : grain_tempo.mass,
        'Inertia' : grain_tempo.inertia
        }
        #create real grain
        L_g.append(Grain.Grain(dict_ic_to_real))

    #Add element in dict
    dict_sample['L_g'] = L_g

#-------------------------------------------------------------------------------
