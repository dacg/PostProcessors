# Local functions for dkdk postprocessing

import os, copy

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
    self.center_ii  = None
    self.veloc  = None
    self.rotation   = None
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
  def __init__(self, center, radius):
    self.shape = 'wall'
    self.center = center
    self.radius = radius
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
  # Getting the number of line
  n_lines_file = sum(1 for line in file_bodies)
  # Rewind the cursor in the file
  file_bodies.seek(0)

  # Initializing variables
  n_disks = 0
  n_polyg = 0
  n_walls = 0
  n_pt2dx = 0

  # Looking for keywords of body shapes
  for jj in range(0,n_lines_file,1):
    # Reading line by lines
    curr_line = file_bodies.readline()

    # Identifying the type of body
    if (curr_line[1:6] == 'DISKx'):
      n_disks+=1
    #
    if (curr_line[1:6] == 'POLYG'):
      n_polyg+=1
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
  # Rewinding the file
  file_bodies.seek(0)

  # Iterating the lines looking for common data among the objects
  for jj in range(0,n_lines_file,1):
    # Reading line by line
    curr_line = file_bodies.readline()

    # Identifying the type of body
    if (curr_line[0:6] == '$bdyty'):
      curr_line = file_bodies.readline()
      elements_line = curr_line.split()
      # Looking for the id of the body
      id_body = int(elements_line[1])
      
      curr_line = file_bodies.readline() # Skipable like
      
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
      center = [float(elements_line[3].replace('D','E')), float(elements_line[5].replace('D','E'))]

      curr_line = file_bodies.readline() # Skipable like

      # Identifying the contactor name and id
      curr_line = file_bodies.readline()
      elements_line = curr_line.split()

      # The ID
      id_body = int(elements_line[1])

      if (elements_line[0] == 'DISKx'):
        # Creating body and appending to the list
        bodies.append(body_disk(id_body,center,radius))
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

  # Initializing containers
  contacts=[]
  
  # Updating the bodies to current configuration
  # Opening the DOF
  file_dof_ii=open("OUTBOX/DOF.OUT." + str(frame_ii), "r")

  # Getting the number of lines
  n_lines_file = sum(1 for line in file_dof_ii)
  # Rewind the cursor in the file
  file_dof_ii.seek(0)

    # Iterating the lines looking for common data among the objects
  for jj in range(0,n_lines_file,1):
    # Reading line by line
    curr_line = file_dof_ii.readline()

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
      print (elements_line)
      bodies[id_body].center_ii[0]=bodies[id_body].center_ini[0]+float(elements_line[3].replace('D','E'))
      bodies[id_body].center_ii[1]=bodies[id_body].center_ini[1]+float(elements_line[5].replace('D','E'))

      # Storing the rotation
      bodies[id_body].rotation=float(elements_line[7].replace('D','E'))

      # Reading the velocity
      curr_line = file_dof_ii.readline()
      curr_line = curr_line.replace('=',' ')

      # Storing the velocity
      bodies[id_body].veloc = [float(elements_line[3].replace('D','E')),float(elements_line[5].replace('D','E')),0.]

      # Storing the angular velocity
      bodies[id_body].veloc[2] = float(elements_line[7].replace('D','E'))

      # Deep copying the body

      # Reading the 


  # Opening the Vloc
  file_vr_ii=open("OUTBOX/Vloc_Rloc.OUT." + str(frame_ii), "r")

  return bodies,contacts
#

