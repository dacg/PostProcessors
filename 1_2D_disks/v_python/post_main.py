############################################################
###########  POSTPROCESSOR FOR ASSEMBLIES OF 2D  ###########
############################################################

# Importing functions
import sys, glob

# Local functions
from functions_2D import *

############################################################
# Setting up the postprocessor
############################################################

# Defining the frames to analyze
h_frame = 'all'                 # Options: 'all': for all frames in the OUTBOX. 
                                #         [initial_frame, last_frame]

# Defining the parameters to compute
c_walls_position = False #True
c_walls_forces   = False #True
c_compacity      = False #True
c_qoverp         = False #True
c_c_anisotropy   = False #True
c_f_anisotropy   = True #True
c_b_anisotropy   = True #True


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

print ('Frames to analyze: ', init_frame, ' to ', last_frame)

# One time functions
# Initializing bodies 
bodies = init_bodies()

# Looping the desired frames
for ii in range(init_frame,last_frame,1):

  # Printing info on this step
  print ('################################')
  print ('Analyzing frame: ', str(ii))

  # Reading the historic file. This allocates contacts and updates variables in bodies for current configuration
  bodies, contacts = read_historic(ii,bodies)

  # Computing the parameters that are needed
  if (c_walls_position): walls_position(ii,init_frame,last_frame,bodies,contacts)
  if (c_walls_forces)  : walls_forces(ii, init_frame,last_frame,bodies,contacts)
  if (c_compacity)     : compacity(ii,init_frame,last_frame,bodies,contacts)
  if (c_qoverp)        : qoverp(ii,init_frame,last_frame,bodies,contacts)
  if (c_c_anisotropy)  : c_anisotropy(ii,init_frame,last_frame,bodies,contacts)
  if (c_f_anisotropy)  : f_anisotropy(ii,init_frame,last_frame,bodies,contacts)
  if (c_b_anisotropy)  : b_anisotropy(ii,init_frame,last_frame,bodies,contacts)
  #

#

