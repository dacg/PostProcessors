############################################################
########### POSTPROCESSOR FOR ASSEMBLIES OF DKDK ###########
############################################################

# Importing functions
import sys, glob

# Local functions
from functions_dkdk import *

############################################################
# Setting up the postprocessor
############################################################

# Defining the frames to analyze
h_frame = 'all'                 # Options: 'all': for all frames in the OUTBOX. 
                                #         [initial_frame, last_frame]

# Defining the parameters to compute
c_wall_position = True
c_wall_forces   = True
c_compacity     = True
c_shape_ratio   = True

############################################################
# The postprocessor
############################################################

# Defining the initial frame and last frame to analyse
# If all the frames, we look for the number using glob
if (h_frame == 'all'):
  # We suppose the list begins at 1
  init_frame = 1
  # We look for the number of files in the OUTBOX
  last_frame = len(glob.glob1('OUTBOX',"DOF.OUT.*"))
else:
  # Otherwise we read the entries
  init_frame = h_frame[0]
  last_frame = h_frame[1]
#

print ('Analysing frames: ', init_frame, ' to ', last_frame)

# One time functions
# Initializing bodies 
bodies = init_bodies()

# Looping the desired frames
for ii in range(init_frame,last_frame,1):

  # Reading the historic file. This allocates contacts and updates variables in bodies for current configuration
  bodies, contacts = read_historic(ii,bodies)

  print (bodies[0].veloc)

  # Computing the parameters that are needed
  #if (c_wall_position): wall_position(ii,last_frame,bodies_current,contacts)
  #if (c_wall_forces): wall_forces(ii,last_frame,bodies_current,contacts)
  #if (c_compacity): compacity(ii,last_frame,bodies_current,contacts)
  #

#

