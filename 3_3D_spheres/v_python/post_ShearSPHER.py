############################################################
# POSTPROCESSOR FOR SIMPLE SHEARING OF SPHERES ASSEMBLIES ##
############################################################

# Importing functions
import math, os, numpy, string, random, shutil,sys,time, glob

# Local functions
from functions_ShearSPHER import *
dim = 3

############################################################
# Setting up the postprocessor
############################################################

# Defining the frames to analyze
h_frame = 'all'                 # Options: 'all': for all frames in the OUTBOX. 
                                #         [initial_frame, last_frame]

# Defining the parameters to compute
c_packing_frac  = True
c_Lacey	        = True ; mLS = 1.; size_ratio = 3.; radius_min = 0.0006
#c_wall_forces   = True
#c_compacity     = True
#c_shape_ratio   = True

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
result_file = open("POSTPRO/RESULT.DAT","w")
# Looping the desired frames
for ii in range(init_frame,last_frame+1,1):
  result_file.write('%s \t %d \n'%('Frame:',ii))
  result_file.write('%s \n'%('------'))
  # Reading the historic file. This allocates contacts and updates variables in bodies for current configuration
  bodies = read_historic(ii,bodies)
  vol_asmbly = 0.
  # Computing the parameters that are needed
  if (c_packing_frac):
    result_file.write('%s \t %s \t\t %s \t\t %s \n'%('Body:','CoorX','CoorY','CoorZ'))
    for i in range(0,len(bodies),1):
      if bodies[i].shape == 'PLAN':
        result_file.write('%d \t %10.4E \t %10.4E \t %10.4E \n'%(i+1,bodies[i].center_ii[0],bodies[i].center_ii[1],bodies[i].center_ii[2]))
      else:
        vol_asmbly += 4./3.*math.pi*bodies[i].radius**3
    h0 = bodies[1].center_ii[2]+bodies[0].center_ii[2]
    h1 = bodies[3].center_ii[0]+bodies[2].center_ii[0]
    h2 = bodies[5].center_ii[1]+bodies[4].center_ii[1]
    vol_box = h0*h1*h2
    packing = vol_asmbly/vol_box
    result_file.write('%s \n'%('------'))
    result_file.write('%s \t %10.4E \n'%('Volume BOX=',vol_box))
    result_file.write('%s \n'%('------'))
    result_file.write('%s \t %10.4E \n'%('Volume PART=',vol_asmbly))
    result_file.write('%s \n'%('------'))
    result_file.write('%s \t %10.4E \n'%('Packing fraction=',packing))
    result_file.write('%s \n'%('======'))

  if(c_Lacey):
    p = mLS/(mLS+1.)
    count_box = 0
    small_lowcomp = 0
    small_highcomp = 0
    coarse_lowcomp = 0
    coarse_highcomp = 0
    h0 = bodies[1].center_ii[2]+bodies[0].center_ii[2]
    for i in range(0,len(bodies),1):
      if bodies[i].shape == 'PLAN':
      	count_box += 1
      else :
      	if bodies[i].radius < radius_min:
      		if bodies[i].center_ii[2] < h0/2. :
      		  small_lowcomp += 1
      		else :
      		  small_highcomp += 1
      	else :
      		if bodies[i].center_ii[2] < h0/2. :
      		  coarse_lowcomp += 1
      		else :
      		  coarse_highcomp += 1
    S02 = p*(1-p)
    Sr2 = p*(1-p)*2./(len(bodies)-count_box)
    sigmaS = numpy.array([small_lowcomp/(small_lowcomp+coarse_lowcomp*size_ratio**dim),small_highcomp/(small_highcomp+coarse_highcomp*size_ratio**dim)])
    sigmaC = numpy.array([coarse_lowcomp*size_ratio**dim/(small_lowcomp+coarse_lowcomp*size_ratio**dim),coarse_highcomp*size_ratio**dim/(small_highcomp+coarse_highcomp*size_ratio**dim)])
    smallM = (Sr2 - numpy.var(sigmaS))/(Sr2-S02)
    coarseM = (Sr2 - numpy.var(sigmaC))/(Sr2-S02)
    result_file.write('%s \t %10.4E \n'%('small Lacey=',smallM))
    result_file.write('%s \t %10.4E \n'%('coarse Lacey=',coarseM))

  #if (c_packing_frac):
  #if (c_wall_forces): wall_forces(ii,last_frame,bodies_current,contacts)
  #if (c_compacity): compacity(ii,last_frame,bodies_current,contacts)
  #
  result_file.write('%s \n'%('#################################################'))
  result_file.write('\n')
result_file.close()
#

