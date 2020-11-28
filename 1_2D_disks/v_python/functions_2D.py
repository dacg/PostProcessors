# Local functions for 2D postprocessing

import os, sys, copy
import numpy as np

# Global variable
# In case of periodic interactions with polyg or clusters, it's better to explicitely 
# define the periodicity
global len_periodic
len_periodic = 1592.9464


############################################################
# Defining classes 
############################################################

# Creating an object disk 
class body_disk():
  def __init__(self, id_body, center_ini, radius):
    self.shape      = 'disk'
    self.id_body    = id_body
    self.center_ini = center_ini
    self.radius     = radius
    self.center_ii  = np.array([0.,0.])   # Default initialization
    self.veloc      = np.array([0.,0.,0]) # Default initialization
    self.rotation   = 0.        # Default initialization

    self.area       = np.pi*radius*radius
  #
#

# Creating an object polygon 
class body_poly():
  def __init__(self, center, radius):
    self.shape = 'poly'
    self.center = center
    self.radius = radius
  #
#

# Creating an object wall 
class body_wall():
  def __init__(self, id_body, center, radius, ax1, ax2):
    self.shape = 'wall'
    self.id_body    = id_body
    self.center_ini = center
    self.center_ii  = [0.,0.,0]
    self.radius     = radius
    self.ax1        = ax1
    self.ax2        = ax2
  #
#

# Creating an object cluster 
class body_cluster():
  def __init__(self, id_body, center_ini, radius,vertices,radius_c1,radius_c2, center_c1, center_c2):
    self.shape         = 'cluster'
    self.id_body       = id_body
    self.center_ini    = center_ini
    self.radius        = radius
    self.radius_c1     = radius_c1
    self.radius_c2     = radius_c2
    self.center_ii     = np.array([0.,0.]) # Default initialization
    self.center_c1     = center_c1
    self.center_c2     = center_c2
    self.center_c1_ii  = np.zeros_like(center_c1)
    self.center_c2_ii  = np.zeros_like(center_c1)
    self.veloc         = np.array([0.,0.,0]) # Default initialization
    self.rotation      = 0.        # Default initialization
    self.vertices_ini  = vertices
    self.vertices_ii   = np.zeros_like(vertices) # Default initialization

    # Computing area
    self.area = 0.

    self.area += radius_c1*radius_c1*np.pi
    area_rect = abs(vertices[0,0]*vertices[1,1] - vertices[0,1]*vertices[1,0] + \
                    vertices[1,0]*vertices[2,1] - vertices[1,1]*vertices[2,0] + \
                    vertices[2,0]*vertices[3,1] - vertices[2,1]*vertices[3,0] + \
                    vertices[3,0]*vertices[0,1] - vertices[3,1]*vertices[0,0])
    area_rect = area_rect/2.
    self.area += area_rect
  #
#

# Creating an object single contact  
class contact():
  def __init__(self, inter_type, id_ctc, id_cd, id_an, status, 
               i_law, rn, rt, vln, vlt, gap, n_frame, coor_ctc):
    self.inter_type = inter_type
    self.id_ctc     = id_ctc
    self.id_cd      = id_cd
    self.id_an      = id_an
    self.status     = status
    self.i_law      = i_law
    self.rn         = rn
    self.rt         = rt
    self.vln        = vln
    self.vlt        = vlt
    self.gap        = gap
    self.n_frame    = n_frame
    self.t_frame    = np.zeros_like(n_frame)
    self.coor_ctc   = coor_ctc
    self.n_interact = 1

    # Setting the tangential frame
    self.t_frame[0] =  n_frame[1]
    self.t_frame[1] = -n_frame[0]
  #
#

############################################################
# Function that initialized the bodies
############################################################
def init_bodies():

  # Initalizing a container for several bodies
  bodies = []

  # Opening the file of BODIES
  file_bodies=open("OUTBOX/BODIES.OUT", "r")

  # Initializing number of every type of body. THESE ARE GLOBAL VARIABLES
  global n_disks
  global n_polyg
  global n_walls
  global n_pt2dx
  global n_clust
  global n_bodies

  n_disks  = 0
  n_polyg  = 0
  n_walls  = 0
  n_pt2dx  = 0
  n_clust  = 0
  n_bodies = 0

  # Looking for keywords of body shapes
  while (True):
    # Reading line by lines
    curr_line = file_bodies.readline()
    
    # Checking if this is the end of the file. Note: End of file is an empty string
    if (len(curr_line) == 0):
      break
    #

    # Identifying the type of body
    if (curr_line[1:6] == 'DISKx'):
      n_disks+=1
    #
    if (curr_line[1:6] == 'POLYG'):
      # Checking if this is a cluster of a simple polyg
      elements_line = curr_line.split()
      n_sides = int(elements_line[-1])

      # Skipping the vertices info
      for kk in range (0,n_sides,1): 
        curr_line = file_bodies.readline()
      #

      # Getting the next line to check if this is a cluster
      curr_line = file_bodies.readline()

      if (curr_line[1:6] == 'DISKb'):
        n_clust += 1
      else :
        n_polyg+=1
      #
    #
    if (curr_line[1:6] == 'JONCx'):
      n_walls+=1
    #
    if (curr_line[1:6] == 'PTPT2'):
      n_pt2dx+=1
    #
  #

  print ('Number of disks: ',    n_disks)
  print ('Number of polygons: ', n_polyg)
  print ('Number of walls: ',    n_walls)
  print ('Number of points: ',   n_pt2dx)
  print ('Number of clusters: ', n_clust)

  # Total number of bodies
  n_bodies = n_disks+n_polyg+n_walls+n_pt2dx+n_clust
  print ('Total number of objects: ', n_bodies)

  # Rewinding the file
  file_bodies.seek(0)

  # Iterating the lines looking for common data among the objects
  while (True):
    # Reading line by line
    curr_line = file_bodies.readline()
    # Checking if this is the end of the file. Note: End of file is an empty string
    if (len(curr_line) == 0):
      break
    #

    # Identifying the type of body
    if (curr_line[0:6] == '$bdyty'):
      curr_line = file_bodies.readline()
      elements_line = curr_line.split()

      # Looking for the id of the body
      id_body = int(elements_line[1])
      
      curr_line = file_bodies.readline() # Skipable line
      
      # Getting the radius of the body
      curr_line = file_bodies.readline()
      elements_line = curr_line.split()
      radius = float(elements_line[5].replace('D','E'))

      curr_line = file_bodies.readline() # Skipable like

      # Getting the initial center of the body
      curr_line = file_bodies.readline()
      # Careful. Coordinates can be negative. So replacing to make an easy split
      curr_line = curr_line.replace('=',' ')
      elements_line = curr_line.split()
      center = np.array([float(elements_line[3].replace('D','E')), float(elements_line[5].replace('D','E'))])

      curr_line = file_bodies.readline() # Skipable like

      # Identifying the contactor name and id
      curr_line = file_bodies.readline()
      elements_line = curr_line.split()

      if (elements_line[0] == 'DISKx'):
        # Creating body and appending to the list
        bodies.append(body_disk(id_body,center,radius))
      if (elements_line[0] == 'POLYG'):
        # Checking if this is a cluster of a simple polyg
        n_sides = int(elements_line[-1])

        # Creating array with the vertices information
        curr_vertices = np.zeros((n_sides,2))

        # Reading the vertices coordinations
        for kk in range (0,n_sides,1): 
          curr_line = file_bodies.readline()
          # Careful. Coordinates can be negative. So replacing to make an easy split
          curr_line = curr_line.replace('=',' ')
          elements_line = curr_line.split()
          curr_vertices[kk,0] = float(elements_line[1].replace('D','E'))
          curr_vertices[kk,1] = float(elements_line[3].replace('D','E'))
        #

        # Getting the next line to check if this is a cluster
        curr_line = file_bodies.readline()
        elements_line = curr_line.split()
        if (elements_line[0] == 'DISKb'):
          # Getting the radius of C1
          radius_c1 = float(elements_line[-1].replace('D','E'))
          # Getting the center of C1
          center_c1 = np.array([0.,0.])
          curr_line = file_bodies.readline()
          curr_line = curr_line.replace('=',' ')
          elements_line = curr_line.split()
          center_c1[0] = float(elements_line[1].replace('D','E'))
          center_c1[1] = float(elements_line[3].replace('D','E'))

          # Getting the radius for C2
          curr_line = file_bodies.readline()
          radius_c2 = float(elements_line[-1].replace('D','E'))
          # Getting the center of C2
          center_c2 = np.array([0.,0.])
          curr_line = file_bodies.readline()
          curr_line = curr_line.replace('=',' ')
          elements_line = curr_line.split()
          center_c2[0] = float(elements_line[1].replace('D','E'))
          center_c2[1] = float(elements_line[3].replace('D','E'))

          # Building the object cluster
          bodies.append(body_cluster(id_body,center,radius,curr_vertices,radius_c1,radius_c2, center_c1, center_c2))
        else :
          print ('error particle case')
          sys.exit()
        #
      #
      if (elements_line[0] == 'JONCx'):
        ax1 = float(elements_line[6].replace('D','E'))
        ax2 = float(elements_line[9].replace('D','E'))
        bodies.append(body_wall(id_body,center,radius,ax1,ax2))
      #
    #
  #

  # Closing the file of BODIES
  file_bodies.close()

  return bodies
#

############################################################
# Function to read the historic files
############################################################
def read_historic(frame_ii,bodies):

  # Initializing container of contacts
  # 
  contacts_raw = []
  contacts     = []
  
  # Updating the bodies to current configuration
  # Opening the DOF
  file_dof_ii=open("OUTBOX/DOF.OUT." + str(frame_ii), "r")

  # Skipable lines
  curr_line = file_dof_ii.readline()
  curr_line = file_dof_ii.readline()
  curr_line = file_dof_ii.readline()

  # Extracting the time. THIS IS A GLOBAL VARIABLE
  curr_line = file_dof_ii.readline()
  elements_line = curr_line.split()

  global curr_time 

  curr_time = float(elements_line[3].replace('D','E'))

  # Iterating the lines looking for common data among the objects
  while (True):
    # Reading line by line
    curr_line = file_dof_ii.readline()

    # Checking if this is the end of the file. Note: End of file is an empty string
    if (len(curr_line) == 0):
      break
    #

    # Identifying the type of body
    if (curr_line[0:6] == '$bdyty'):
      curr_line = file_dof_ii.readline()
      elements_line = curr_line.split()
      # Looking for the id of the body
      id_body = int(elements_line[1])

      curr_line = file_dof_ii.readline() # Skipable like

      # Careful. Displacements can be negative. So replacing to make an easy split
      curr_line = file_dof_ii.readline()
      curr_line = curr_line.replace('=',' ')
      elements_line = curr_line.split()
      # Storing the particle displacement
      disp_ii = np.array([0.,0.])
      disp_ii[0] = float(elements_line[3].replace('D','E'))
      disp_ii[1] = float(elements_line[5].replace('D','E'))

      bodies[id_body-1].center_ii[0]=bodies[id_body-1].center_ini[0]+disp_ii[0]
      bodies[id_body-1].center_ii[1]=bodies[id_body-1].center_ini[1]+disp_ii[1]

      # Storing the rotation
      bodies[id_body-1].rotation=float(elements_line[7].replace('D','E'))

      # Reading the velocity
      curr_line = file_dof_ii.readline()
      curr_line = curr_line.replace('=',' ')

      # Storing the linear and angular velocit
      bodies[id_body-1].veloc = np.array([float(elements_line[3].replace('D','E')), 
                               float(elements_line[5].replace('D','E')),
                               float(elements_line[7].replace('D','E'))])

      # Building the rotation matrix
      mat_rot = np.zeros((2,2))
      mat_rot[0,0] =  np.cos(bodies[id_body-1].rotation)
      mat_rot[1,0] =  np.sin(bodies[id_body-1].rotation)
      mat_rot[1,0] = -np.sin(bodies[id_body-1].rotation)
      mat_rot[1,1] =  np.cos(bodies[id_body-1].rotation)

      # Updating vertices if necessary
      if (bodies[id_body-1].shape == 'cluster'):
        for kk in range (0, len(bodies[id_body-1].vertices_ii),1):
          bodies[id_body-1].vertices_ii[kk,:]  = mat_rot.dot(bodies[id_body-1].vertices_ini[kk,:]) + bodies[id_body-1].center_ii[:]
        #
        bodies[id_body-1].center_c1_ii = mat_rot.dot(bodies[id_body-1].center_c1) + bodies[id_body-1].center_ii
        bodies[id_body-1].center_c2_ii = mat_rot.dot(bodies[id_body-1].center_c2) + bodies[id_body-1].center_ii
      #
    #
  #

  # Closing the DOF file
  file_dof_ii.close()

  # Contacts information
  
  # Initializing the number of contact types
  n_dkdk = 0
  n_plpl = 0
  n_dkpl = 0
  n_dkjc = 0
  n_pljc = 0

  # Opening the Vloc
  file_name = "OUTBOX/Vloc_Rloc.OUT." + str(frame_ii)
  file_vr_ii=open(file_name, "r")

  # Iterating the lines looking for common data among the objects
  while (True):
    # Reading line by line
    curr_line = file_vr_ii.readline()

    # Checking if this is the end of the file. Note: End of file is an empty string
    if (len(curr_line) == 0):
      break
    #

    # Spliting line
    elements_line = curr_line.split()

    # Vloc_Rloc files can have empty lines. In order to skip so: 
    if (len(elements_line)<2): continue

    # Reading contact structure
    if (elements_line[0] == '$icdan'):
      # Reading the type of contact
      inter_type = elements_line[1]

      if (inter_type == 'DKPLx'): 
        n_dkpl += 1
      #
      if (inter_type == 'DKDKx'): 
        n_dkdk += 1
      #
      if (inter_type == 'PLPLx'): 
        n_plpl += 1
      #
      if (inter_type == 'PLJCx'): 
        n_pljc += 1
      #
      if (inter_type == 'DKJCx'): 
        n_dkjc += 1
      #
      id_ctc = int(elements_line[2])
      curr_line = file_vr_ii.readline() # Skipable line

      # Getting info from the bodies in interaction
      curr_line = file_vr_ii.readline()
      elements_line = curr_line.split()

      # The bodies in contact. Note: PLJC and PLPL have different structure
      if (inter_type == 'PLPLx'): 
        id_cd = int(elements_line[1])
        id_an = int(elements_line[8])
        i_law = elements_line[6]
        status = elements_line[13]
      elif(inter_type == 'PLJCx'): 
        id_cd = int(elements_line[1])
        id_an = int(elements_line[8])
        i_law = elements_line[6]
        status = elements_line[11]
      else :
        id_cd = int(elements_line[1])
        id_an = int(elements_line[6])
        i_law = elements_line[4]
        status = elements_line[9]
      #

      # Reading forces
      curr_line = file_vr_ii.readline()
      curr_line = curr_line.replace('H',' ')
      elements_line = curr_line.split()

      rn = float(elements_line[3].replace('D','E'))
      rt = float(elements_line[1].replace('D','E'))

      # Reading relative velocities
      curr_line = file_vr_ii.readline()
      curr_line = curr_line.replace('=',' ')
      elements_line = curr_line.split()

      vln = float(elements_line[3].replace('D','E'))
      vlt = float(elements_line[1].replace('D','E'))

      # Reading the gap
      curr_line = file_vr_ii.readline()
      curr_line = curr_line.replace('TT','  ')
      elements_line = curr_line.split()

      gap = float(elements_line[1].replace('D','E'))

      # Reading the local frame
      n_frame = np.array([0.,0.])

      curr_line = file_vr_ii.readline()
      curr_line = curr_line.replace('=',' ')
      elements_line = curr_line.split()

      n_frame[0] = float(elements_line[1].replace('D','E'))
      n_frame[1] = float(elements_line[3].replace('D','E'))

      # Reading the local frame
      coor_ctc = [0.,0.]

      curr_line = file_vr_ii.readline()
      curr_line = curr_line.replace('=',' ')
      elements_line = curr_line.split()

      coor_ctc[0] = float(elements_line[1].replace('D','E'))
      coor_ctc[1] = float(elements_line[3].replace('D','E'))

      # Building the object
      contacts_raw.append(contact(inter_type, id_ctc, id_cd, id_an, status, i_law, rn, rt, vln, vlt, gap, n_frame, coor_ctc))
    #
  #

  # Closing the Vloc_Rloc file
  file_vr_ii.close()

  # To count the number of double interactions
  n_double = 0

  # Checking if the contact list has more than simple interactions and summarizing if necessary
  for kk in range (0, len(contacts_raw), 1):
    # We can verify up to the contact before the last 
    if (kk < len(contacts_raw)-1):
      # Those that can be douple are PLPL or PLJC
      if (contacts_raw[kk].inter_type == 'PLPLx' or contacts_raw[kk].inter_type == 'PLJCx'):
        # Checking is this is a double 
        ## They should be the same type
        if (contacts_raw[kk].inter_type == contacts_raw[kk].inter_type):
          # They should not have been already counted
          if (contacts_raw[kk].n_interact == 1): 
            # They should have the same bodies in interaction
            if (contacts_raw[kk].id_cd == contacts_raw[kk+1].id_cd and contacts_raw[kk].id_an == contacts_raw[kk+1].id_an):
              # Then we can add the double interaction
              # The sum of forces
              rn_total = contacts_raw[kk].rn + contacts_raw[kk+1].rn
              rt_total = contacts_raw[kk].rt + contacts_raw[kk+1].rt
              # The average velocity
              vln_total = (contacts_raw[kk].vln + contacts_raw[kk+1].vln)/2.
              vlt_total = (contacts_raw[kk].vlt + contacts_raw[kk+1].vlt)/2.
              # The average gap
              gap_total = (contacts_raw[kk].gap + contacts_raw[kk+1].gap)/2
              # The average coordinations of the contact point
              coor_ctc_total = np.zeros_like(contacts_raw[kk].coor_ctc)
              coor_ctc_total[0] = (contacts_raw[kk].coor_ctc[0] + contacts_raw[kk+1].coor_ctc[0])/2.
              coor_ctc_total[1] = (contacts_raw[kk].coor_ctc[1] + contacts_raw[kk+1].coor_ctc[1])/2.
              #
              # We change the number of interactions in the raw list
              contacts_raw[kk].n_interact = 2
              contacts_raw[kk+1].n_interact = 2
              #
              # We count the number of double interactions
              n_double += 1
              #
              # The common data for both raw interactions 
              # The contact frame is the same
              n_frame_total = contacts_raw[kk].n_frame
              # We keep the id_ctc of the kk interaction
              id_ctc = contacts_raw[kk].id_ctc
              #
              id_cd      = contacts_raw[kk].id_cd
              id_an      = contacts_raw[kk].id_an
              inter_type = contacts_raw[kk].inter_type
              status     = contacts_raw[kk].status
              i_law      = contacts_raw[kk].i_law
              #
              # We build the new summarized object and we add it to the clean list
              contacts.append(contact(inter_type, id_ctc, id_cd, id_an, status,i_law,
                              rn_total, rt_total, vln_total, vlt_total, gap_total, n_frame_total, coor_ctc_total ))
            #
          #
        #
      else:
        contacts.append(contacts_raw[kk])
      #
    else:
      # In the case is not double, we add it to the clean list of contacts
      contacts.append(contacts_raw[kk])
    #
  #

  print ('Number of dkdk: ',    n_dkdk)
  print ('Number of dkpl: ',    n_dkpl)
  print ('Number of plpl: ',    n_plpl)
  print ('Number of pljc: ',    n_pljc)
  print ('Number of dkjc: ',    n_dkjc)
  print ('Total number of raw contacts: ', n_dkdk+n_plpl+n_dkpl+n_pljc+n_dkjc)
  print ('Total number of double contacts: ', n_double)
  return bodies, contacts
#

############################################################
# Function that computes the position of the walls
############################################################
def walls_position(ii,init_frame,last_frame,bodies,contacts):

  # Checking if we need to open the file
  if (ii == init_frame): 
    global file_w_pos
    file_w_pos = open("POSTPRO/WALLS_POSITION.DAT","w")
    # Then we print the heading of the file
    file_w_pos.write('#   Time    ')
    for kk in range (0,n_walls,1):
      file_w_pos.write('   Wall X-' + str(kk+1) + '  ' + '   Wall Y-' + str(kk+1) + '  ')
    #
    # Changing line
    file_w_pos.write('\n')
  #

  # Printing current time
  file_w_pos.write(format_e(curr_time))

  # Identifying the walls
  for kk in range(0,len(bodies),1):
    if (bodies[kk].shape == 'wall'):
      file_w_pos.write(format_e(bodies[kk].center_ii[0]))
      file_w_pos.write(format_e(bodies[kk].center_ii[1]))
    #
  #

  # Ending this line
  file_w_pos.write('\n')

  # Checking if we have to close the file
  if (ii==last_frame):
    file_w_pos.close()
  #
#

############################################################
# Function that computes the forces on the walls
############################################################
def walls_forces(ii,init_frame,last_frame,bodies,contacts):

  # Checking if we need to open the file
  if (ii == init_frame): 
    global file_w_frc
    file_w_frc = open("POSTPRO/WALLS_FORCES.DAT","w")
    # Then we print the heading of the file
    file_w_frc.write('#   Time    ')
    for kk in range (0,n_walls,1):
      file_w_frc.write('   Wall X-' + str(kk+1) + '  ' + '   Wall Y-' + str(kk+1) + '  ')
    #
    # Changing line
    file_w_frc.write('\n')
  #

  # Printing current time
  file_w_frc.write(format_e(curr_time))

  # Preparing an array to sum the forces. Fx, Fy
  forces_sum = np.zeros((n_walls,2))

  # Looping the contacts
  for kk in range(0,len(contacts),1):
    # Identifying contacts with the walls
    if (contacts[kk].inter_type == 'DKJCx' or contacts[kk].inter_type == 'PLJCx'):
      # Identifying the ID of the wall. It is always the antagonist!
      id_wall = contacts[kk].id_an

      # For simplicity
      n_no_wall = n_bodies - n_walls

      # The force vector
      f_vec = contacts[kk].rn * contacts[kk].n_frame + contacts[kk].rt*contacts[kk].t_frame
      
      # Adding the force
      forces_sum[id_wall-n_no_wall-1,0] += f_vec[0] # Component X
      forces_sum[id_wall-n_no_wall-1,1] += f_vec[1] # Component Y
    #
  #

  # Writing the resume
  for kk in range (0,len(forces_sum),1):
    file_w_frc.write(format_e(forces_sum[kk,0]))
    file_w_frc.write(format_e(forces_sum[kk,1]))
  #

  # Ending this line
  file_w_frc.write('\n')

  # Checking if we have to close the file
  if (ii==last_frame):
    file_w_frc.close()
  #
#

############################################################
# Function that computes the position of the walls
############################################################
def compacity(ii,init_frame,last_frame,bodies,contacts):

  # Checking if we need to open the file
  if (ii == init_frame): 
    global file_compacity
    file_compacity = open("POSTPRO/COMPACITY.DAT","w")
    # Then we print the heading of the file
    file_compacity.write('#   Time    ' + '    Height   ' + '     Width   ' + '   Compacity ')
    # Changing line
    file_compacity.write('\n')
  #

  # Printing current time
  file_compacity.write(format_e(curr_time))

  # We need to know the size of the box
  box_size(bodies)

  v_solid = 0

  # The solid volume of particles
  for kk in range(0,len(bodies),1):
    if (bodies[kk].shape == 'wall'): continue
    v_solid += bodies[kk].area
  #

  # Writing in the file
  file_compacity.write(format_e(box_height))
  file_compacity.write(format_e(box_width))
  file_compacity.write(format_e(v_solid/(box_height*box_width)))

  # Ending this line
  file_compacity.write('\n')

  # Checking if we have to close the file
  if (ii==last_frame):
    file_compacity.close()
  #
#

############################################################
# Function that computes the stress tensor
############################################################
def qoverp(ii,init_frame,last_frame,bodies,contacts):

  # 
  box_size(bodies)

  # Checking if we need to open the file
  if (ii == init_frame): 
    global file_qoverp
    file_qoverp = open("POSTPRO/QOVERP.DAT","w")
    # Then we print the heading of the file
    file_qoverp.write('#   Time    ' + '     S1     ' + '      S2     ' + '     Q/P    ')
    # Changing line
    file_qoverp.write('\n')
  #

  # Printing current time
  file_qoverp.write(format_e(curr_time))

  # We build the stress tensor

  stress_tensor = np.zeros((2,2))

  # Over all the contacts 
  for kk in range(0,len(contacts),1):
    # We skip contacts with the walls
    if (contacts[kk].inter_type == 'DKJCx' or contacts[kk].inter_type == 'PLJCx'): continue
    # Active contacts
    if (contacts[kk].rn < 1e-8): continue
    
    # Force vector
    f_vec = contacts[kk].rn*contacts[kk].n_frame + contacts[kk].rt*contacts[kk].t_frame

    # The branch vector
    cd = contacts[kk].id_cd
    an = contacts[kk].id_an

    b_vec = bodies[cd-1].center_ii - bodies[an-1].center_ii

    if ((b_vec[0]**2 + b_vec[1]**2)**0.5 > box_width/2.):
      b_vec[0] = b_vec[0] - len_periodic*(b_vec[0]/abs(b_vec[0]))
    #

    #stress_tensor += np.outer(f_vec,b_vec)
    stress_tensor[0,0] += f_vec[0]*b_vec[0]
    stress_tensor[0,1] += f_vec[0]*b_vec[1]
    stress_tensor[1,0] += f_vec[1]*b_vec[0]
    stress_tensor[1,1] += f_vec[1]*b_vec[1]
  #

  # From moment to stress
  stress_tensor = stress_tensor/(box_width*box_height)

  eig_val, eig_vec = np.linalg.eig(stress_tensor)

  S1 = max(eig_val)
  S2 = min(eig_val)

  p = (S1 + S2)/2.
  q = (S1 - S2)/2.

  # Writing in the file
  file_qoverp.write(format_e(S1))
  file_qoverp.write(format_e(S2))
  file_qoverp.write(format_e(q/p))

  # Ending this line
  file_qoverp.write('\n')

  # Checking if we have to close the file
  if (ii==last_frame):
    file_qoverp.close()
  #
#

############################################################
# Function that computes the contact orientation anisotropy
############################################################
def c_anisotropy(ii,init_frame,last_frame,bodies,contacts):

  # 
  box_size(bodies)

  # Checking if we need to open the file
  if (ii == init_frame): 
    global file_c_aniso
    file_c_aniso = open("POSTPRO/CTC_ANISOTROPY.DAT","w")
    # Then we print the heading of the file
    file_c_aniso.write('#   Time    ' + '     ac     ')
    # Changing line
    file_c_aniso.write('\n')
  #

  # Printing current time
  file_c_aniso.write(format_e(curr_time))

  # We build the fabric tensor
  fabric_tensor = np.zeros((2,2))
  valid_ctc = 0

  # Over all the contacts 
  for kk in range(0,len(contacts),1):
    # We skip contacts with the walls
    if (contacts[kk].inter_type == 'DKJCx' or contacts[kk].inter_type == 'PLJCx'): continue
    # Active contacts
    if (contacts[kk].rn < 1e-8): continue
    
    # The contact frame
    n_frame = contacts[kk].n_frame

    #stress_tensor += np.outer(f_vec,b_vec)
    fabric_tensor[0,0] += n_frame[0]*n_frame[0]
    fabric_tensor[0,1] += n_frame[0]*n_frame[1]
    fabric_tensor[1,0] += n_frame[1]*n_frame[0]
    fabric_tensor[1,1] += n_frame[1]*n_frame[1]

    valid_ctc += 1.
  #

  # The average
  fabric_tensor = fabric_tensor/valid_ctc

  eig_val, eig_vec = np.linalg.eig(fabric_tensor)

  S1 = max(eig_val)
  S2 = min(eig_val)

  ac = 2*(S1 - S2)

  # Writing in the file
  file_c_aniso.write(format_e(ac))

  # Ending this line
  file_c_aniso.write('\n')

  # Checking if we have to close the file
  if (ii==last_frame):
    file_c_aniso.close()
  #
#

############################################################
# Function that computes the force magnitude anisotropy
############################################################
def f_anisotropy(ii,init_frame,last_frame,bodies,contacts):

  # 
  box_size(bodies)

  # Checking if we need to open the file
  if (ii == init_frame): 
    global file_f_aniso
    file_f_aniso = open("POSTPRO/FRC_ANISOTROPY.DAT","w")
    # Then we print the heading of the file
    file_f_aniso.write('#   Time    ' + '   afn+ac    '+ '  aft+afn+2ac  ')
    # Changing line
    file_f_aniso.write('\n')
  #

  # Printing current time
  file_f_aniso.write(format_e(curr_time))

  # We build the fabric tensor
  fn_tensor = np.zeros((2,2))
  ft_tensor = np.zeros((2,2))
  f_tensor = np.zeros((2,2))
  
  average_force = 0.
  valid_ctc = 0

  # Over all the contacts 
  for kk in range(0,len(contacts),1):
    # We skip contacts with the walls
    if (contacts[kk].inter_type == 'DKJCx' or contacts[kk].inter_type == 'PLJCx'): continue
    # Active contacts
    if (contacts[kk].rn < 1e-8): continue
    
    # The force vector 
    f_vec = contacts[kk].rn*contacts[kk].n_frame + contacts[kk].rt*contacts[kk].t_frame

    # Adding the force for the average
    average_force += contacts[kk].rn
    valid_ctc += 1

    # The contact frame
    n_frame = contacts[kk].n_frame
    t_frame = contacts[kk].t_frame

    #stress_tensor += np.outer(f_vec,b_vec)
    fn_tensor[0,0] += contacts[kk].rn*n_frame[0]*n_frame[0]
    fn_tensor[0,1] += contacts[kk].rn*n_frame[0]*n_frame[1]
    fn_tensor[1,0] += contacts[kk].rn*n_frame[1]*n_frame[0]
    fn_tensor[1,1] += contacts[kk].rn*n_frame[1]*n_frame[1]

    ft_tensor[0,0] += contacts[kk].rt*n_frame[0]*t_frame[0]
    ft_tensor[0,1] += contacts[kk].rt*n_frame[0]*t_frame[1]
    ft_tensor[1,0] += contacts[kk].rt*n_frame[1]*t_frame[0]
    ft_tensor[1,1] += contacts[kk].rt*n_frame[1]*t_frame[1]
  #

  fn_tensor = fn_tensor/(average_force/valid_ctc)
  ft_tensor = ft_tensor/(average_force/valid_ctc)

  f_tensor = fn_tensor + ft_tensor

  # For the normal force anisotropy
  eig_val_n, eig_vec_n = np.linalg.eig(fn_tensor)

  S1_n = max(eig_val_n)
  S2_n = min(eig_val_n)

  afn = 2*(S1_n - S2_n)/(S1_n + S2_n)

  # Writing in the file
  file_f_aniso.write(format_e(afn))

  # For the tangential force anisotropy
  eig_val_t, eig_vec_t = np.linalg.eig(f_tensor)

  S1_t = max(eig_val_t)
  S2_t = min(eig_val_t)

  aft = 2*(S1_t - S2_t)/(S1_t + S2_t)

  # Writing in the file
  file_f_aniso.write(format_e(aft))

  # Ending this line
  file_f_aniso.write('\n')

  # Checking if we have to close the file
  if (ii==last_frame):
    file_f_aniso.close()
  #
#

############################################################
# Function that computes the force magnitude anisotropy
############################################################
def b_anisotropy(ii,init_frame,last_frame,bodies,contacts):

  # 
  box_size(bodies)

  # Checking if we need to open the file
  if (ii == init_frame): 
    global file_b_aniso
    file_b_aniso = open("POSTPRO/FRC_ANISOTROPY.DAT","w")
    # Then we print the heading of the file
    file_b_aniso.write('#   Time    ' + '   aln+ac    '+ '  alt+aln+2ac  ')
    # Changing line
    file_b_aniso.write('\n')
  #

  # Printing current time
  file_b_aniso.write(format_e(curr_time))

  # We build the fabric tensor
  bn_tensor = np.zeros((2,2))
  bt_tensor = np.zeros((2,2))
  b_tensor = np.zeros((2,2))
  
  average_branch = 0.
  valid_ctc = 0

  # Over all the contacts 
  for kk in range(0,len(contacts),1):
    # We skip contacts with the walls
    if (contacts[kk].inter_type == 'DKJCx' or contacts[kk].inter_type == 'PLJCx'): continue
    # Active contacts
    if (contacts[kk].rn < 1e-8): continue
    
    # The branch vector
    cd = contacts[kk].id_cd
    an = contacts[kk].id_an

    b_vec = bodies[cd-1].center_ii - bodies[an-1].center_ii

    if ((b_vec[0]**2 + b_vec[1]**2)**0.5 > box_width/2.):
      b_vec[0] = b_vec[0] - len_periodic*(b_vec[0]/abs(b_vec[0]))
    #

    # The contact frame
    n_frame = contacts[kk].n_frame
    t_frame = contacts[kk].t_frame

    # Adding the branch for the average
    proy_normal = b_vec[0]*n_frame[0] + b_vec[1]*n_frame[1]
    proy_tangen = b_vec[0]*t_frame[0] + b_vec[1]*t_frame[1]

    average_branch += proy_normal
    valid_ctc += 1

    #stress_tensor += np.outer(f_vec,b_vec)
    bn_tensor[0,0] += proy_normal*n_frame[0]*n_frame[0]
    bn_tensor[0,1] += proy_normal*n_frame[0]*n_frame[1]
    bn_tensor[1,0] += proy_normal*n_frame[1]*n_frame[0]
    bn_tensor[1,1] += proy_normal*n_frame[1]*n_frame[1]

    bt_tensor[0,0] += proy_tangen*n_frame[0]*t_frame[0]
    bt_tensor[0,1] += proy_tangen*n_frame[0]*t_frame[1]
    bt_tensor[1,0] += proy_tangen*n_frame[1]*t_frame[0]
    bt_tensor[1,1] += proy_tangen*n_frame[1]*t_frame[1]
  #

  bn_tensor = bn_tensor/(average_branch/valid_ctc)
  bt_tensor = bt_tensor/(average_branch/valid_ctc)

  b_tensor = bn_tensor + bt_tensor

  # For the normal force anisotropy
  eig_val_n, eig_vec_n = np.linalg.eig(bn_tensor)

  S1_n = max(eig_val_n)
  S2_n = min(eig_val_n)

  aln = 2*(S1_n - S2_n)/(S1_n + S2_n)

  # Writing in the file
  file_b_aniso.write(format_e(aln))

  # For the tangential force anisotropy
  eig_val_t, eig_vec_t = np.linalg.eig(b_tensor)

  S1_t = max(eig_val_t)
  S2_t = min(eig_val_t)

  alt = 2*(S1_t - S2_t)/(S1_t + S2_t)

  # Writing in the file
  file_b_aniso.write(format_e(alt))

  # Ending this line
  file_b_aniso.write('\n')

  # Checking if we have to close the file
  if (ii==last_frame):
    file_b_aniso.close()
  #
#

############################################################
# Function that sets the scientific formating for writing
############################################################
def format_e(text):
  return (str("{:12.5E}".format(text))+' ')
#

############################################################
# Function that defines the size of the box using bodies
############################################################
def box_size(bodies):
  # Initializing variables
  x_min =  9999.
  x_max = -9999.
  y_min =  9999.
  y_max = -9999.

  global box_width
  global box_height

  for kk in range(0,len(bodies),1):
    # Omiting walls and points
    if (bodies[kk].shape == 'wall'): continue

    # In case of disks, we check with the radius
    if (bodies[kk].shape == 'disk'):
      x_min = min(x_min,bodies[kk].center_ii[0] - bodies[kk].radius)
      x_max = max(x_max,bodies[kk].center_ii[0] + bodies[kk].radius)
      y_min = min(y_min,bodies[kk].center_ii[1] - bodies[kk].radius)
      y_max = max(y_max,bodies[kk].center_ii[1] + bodies[kk].radius)
    #
    if (bodies[kk].shape == 'cluster'):
      n_sides = len(bodies[kk].vertices_ii)
      for ll in range (0,n_sides,1):
        x_min = min(x_min,bodies[kk].vertices_ii[ll,0])
        x_max = max(x_max,bodies[kk].vertices_ii[ll,0])
        y_min = min(y_min,bodies[kk].vertices_ii[ll,1])
        y_max = max(y_max,bodies[kk].vertices_ii[ll,1])
      #
      # for diskp 
      x_min = min(x_min,bodies[kk].center_c1_ii[0] - bodies[kk].radius_c1)
      x_max = max(x_max,bodies[kk].center_c1_ii[0] + bodies[kk].radius_c1)
      y_min = min(y_min,bodies[kk].center_c1_ii[1] - bodies[kk].radius_c1)
      y_max = max(y_max,bodies[kk].center_c1_ii[1] + bodies[kk].radius_c1)

      x_min = min(x_min,bodies[kk].center_c2_ii[0] - bodies[kk].radius_c2)
      x_max = max(x_max,bodies[kk].center_c2_ii[0] + bodies[kk].radius_c2)
      y_min = min(y_min,bodies[kk].center_c2_ii[1] - bodies[kk].radius_c2)
      y_max = max(y_max,bodies[kk].center_c2_ii[1] + bodies[kk].radius_c2)
    #
  #
  box_width  = x_max - x_min
  box_height = y_max - y_min
#