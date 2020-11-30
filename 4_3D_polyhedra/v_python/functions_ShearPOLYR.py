# Local functions for ShearPOLYR postprocessing

import os, copy

############################################################
# Defining classes 
############################################################

# Creating an object POLYR 
class body_POLYR():
  def __init__(self, id_body, center_ini, radius, vertices, faces, material, vertex, connectivity):
    self.shape      = 'POLYR'
    self.id_body    = id_body
    self.center_ini = center_ini
    self.radius     = radius
    self.vertices   = vertices
    self.faces      = faces
    self.material   = material
    self.center_ii  = [0., 0., 0.]
    self.veloc      = None
    self.rotation   = None
    self.vertex 	= vertex
    self.connectivity= connectivity
    #self.chi		= chi # Particle's deviation
    #self.psi 		= psi # Sphericity
    #self.alpha 	= alpha # Mean angularity
  #
#

# Creating an object WALL 
class body_PLAN():
  def __init__(self, id_body, center_ini, radius, material):
    self.shape      = 'PLAN'
    self.id_body    = id_body
    self.center_ini = center_ini
    self.radius     = radius
    self.material   = material
    self.center_ii  = [0., 0., 0.]
    self.veloc      = None
    self.rotation   = None
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
  n_POLYR = 0
  n_PLANs = 0

  # Looking for keywords of body shapes
  for jj in range(0,n_lines_file,1):
    # Reading line by lines
    curr_line = file_bodies.readline()
    # Identifying the type of body
    if (curr_line[1:6] == 'POLYR'):
      n_POLYR+=1
    #
    if (curr_line[1:6] == 'PLANx'):
      n_PLANs+=1
    #

  print ('Number of polyhedra: ', n_POLYR)
  print ('Number of plans: ', n_PLANs)
  # Rewinding the file

  # Reouverture du fichier de lecture
  file_bodies.seek(0)
  # Iterating the lines looking for common data among the objects
  for jj in range(0,n_lines_file,1):
    # Reading line by line
    curr_line = file_bodies.readline() # Seems to jump implicitly to the next line.
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
        material = str(elements_line[4])
        radius = float(elements_line[5].replace('D','E'))

        curr_line = file_bodies.readline() # Skipable like
        curr_line = file_bodies.readline() # Skipable like

      # Getting the initial center of the body
        curr_line = file_bodies.readline()
      # Careful. Coordinates can be negative. So replacing to make an easy split
        curr_line = curr_line.replace('=',' ')
        elements_line = curr_line.split()
        center = [float(elements_line[3].replace('D','E')), float(elements_line[5].replace('D','E')), float(elements_line[7].replace('D','E'))]

        curr_line = file_bodies.readline() # Skipable like
        curr_line = file_bodies.readline() # Skipable like

      # Identifying the contactor name and id
        curr_line = file_bodies.readline()
        elements_line = curr_line.split()
      
      # Pour l'instant pas fan de cette partie que je ne comprends pas: pourquoi effacer le precedent id_body ?
      ## The ID
      ##id_body = int(elements_line[1])
        if (elements_line[0] == 'POLYR'):
        # Creating body and appending to the list
          # Adding specific parameters...
          vertices = int(elements_line[5])
          faces = int(elements_line[7])
          
          vertex = []
          for i in range(vertices):
          	curr_line = file_bodies.readline() 
          	curr_line = curr_line.replace('=',' ')
          	elements_line = curr_line.split()
          	vertex += [[float(elements_line[1].replace('D','E')), float(elements_line[3].replace('D','E')), float(elements_line[5].replace('D','E'))]]
          connectivity = []
          for i in range(faces):
            curr_line = file_bodies.readline()
            elements_line = curr_line.split()
            connectivity += [[int(elements_line[1]), int(elements_line[3]), int(elements_line[5])]]

          bodies.append(body_POLYR(id_body,center,radius,vertices,faces,material,vertex,connectivity))
        else :
      	  bodies.append(body_PLAN(id_body,center,radius,material))
      #
    #
  #
  # Closing the file of BODIES
  print(len(bodies))
  return bodies
  file_bodies.close()
#

############################################################
# Function to read the historic files
############################################################
def read_historic(frame_ii,bodies):

  ## Initializing containers
  ## contacts=[]
  
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
      id_body = int(elements_line[1])-1
      #print(id_body)

      curr_line = file_dof_ii.readline() # Skipable like

      # Careful. Displacements can be negative. So replacing to make an easy split
      curr_line = file_dof_ii.readline()
      curr_line = curr_line.replace('=',' ')
      elements_line = curr_line.split()
      # Storing the particle displacement
      #####print(elements_line)
      bodies[id_body].center_ii[0]=bodies[id_body].center_ini[0]+float(elements_line[3].replace('D','E'))
      bodies[id_body].center_ii[1]=bodies[id_body].center_ini[1]+float(elements_line[5].replace('D','E'))
      bodies[id_body].center_ii[2]=bodies[id_body].center_ini[2]+float(elements_line[7].replace('D','E'))

      # Pas adapte aux fichiers generes pour les modeles 3D.
      ## Storing the rotation
      ## bodies[id_body].rotation=float(elements_line[7].replace('D','E'))

      curr_line = file_dof_ii.readline() # Skipable like

      # Reading the velocity
      curr_line = file_dof_ii.readline()
      curr_line = curr_line.replace('=',' ')

      # Storing the velocity
      bodies[id_body].veloc = [float(elements_line[3].replace('D','E')),float(elements_line[5].replace('D','E')),float(elements_line[7].replace('D','E'))]

      # Pas adapte aux fichiers generes pour les modeles 3D.
      ## Storing the angular velocity
      ## bodies[id_body].veloc[2] = float(elements_line[7].replace('D','E'))

      # Deep copying the body

      # Reading the 

  # Pour l'intant pas besoin de compute des infos liees aux contacts.
  ## Opening the Vloc
  ## file_vr_ii=open("OUTBOX/Vloc_Rloc.OUT." + str(frame_ii), "r")

  return bodies #,contacts
  file_dof_ii.close()
#

