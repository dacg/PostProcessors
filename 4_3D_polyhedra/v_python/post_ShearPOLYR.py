############################################################
## POSTPROCESSOR FOR SIMPLE SHEARING OF POLYR ASSEMBLIES ###
############################################################

# Importing functions
import math, os, numpy, numpy.linalg, string, random, shutil,sys,time, glob

# Local functions
from functions_ShearPOLYR import *
dim = 3

############################################################
# Setting up the postprocessor
############################################################

# Defining the frames to analyze
h_frame = 'all'                 # Options: 'all': for all frames in the OUTBOX. 
                                #         [initial_frame, last_frame]

# Defining the parameters to compute
c_packing_frac  = True
c_Lacey	        = True ; mLS = 1.; size_ratio = 3.; radius_min = 0.6
c_eccentricty   = False ; small_vertices = 60
c_sphericiy     = False
c_angularity    = False
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
    result_file.write('%s \n'%('======'))

  if ii == 1:
    if (c_eccentricty):
      avg_small = numpy.array([])
      avg_coarse = numpy.array([])
      for i in range(0,len(bodies),1):
        if bodies[i].shape == 'PLAN':
         continue
        else:
         dist_x = numpy.array([])
         dist_y = numpy.array([])
         dist_z = numpy.array([])
         for v in bodies[i].vertex:
           dist_x = numpy.append(dist_x,numpy.array([v[0]]))
           dist_y = numpy.append(dist_y,numpy.array([v[1]]))
           dist_z = numpy.append(dist_z,numpy.array([v[2]]))
         A = numpy.amax(dist_x)-numpy.amin(dist_x)
         B = numpy.amax(dist_y)-numpy.amin(dist_y)
         C = numpy.amax(dist_z)-numpy.amin(dist_z)
         if A < B :
           if B < C :
             LS = numpy.array([A/C, B/C])
             LS = numpy.array([[numpy.amax(LS), numpy.amin(LS)]])
           else :
             LS = numpy.array([A/B, C/B])
             LS = numpy.array([[numpy.amax(LS), numpy.amin(LS)]])
         else :
           if A < C:
             LS = numpy.array([A/C, B/C])
             LS = numpy.array([[numpy.amax(LS), numpy.amin(LS)]])
           else :
             LS = numpy.array([B/A, C/A])
             LS = numpy.array([[numpy.amax(LS), numpy.amin(LS)]])
         if bodies[i].vertices == small_vertices:
           if len(avg_small) == 0:
             avg_small = LS
           else :
             avg_small = numpy.append(avg_small, LS, axis=0)
         else:
           if len(avg_coarse) == 0:
             avg_coarse = LS
           else :
             avg_coarse = numpy.append(avg_coarse, LS, axis=0)
        small_chi = numpy.mean(avg_small, axis=0)
        coarse_chi = numpy.mean(avg_coarse, axis=0)
        result_file.write('%s \t %10.4E \t %s \t %10.4E \n'%('small particles deviation=',small_chi[0],';',small_chi[1]))
        result_file.write('%s \t %10.4E \t %s \t %10.4E \n'%('coarse particles deviation=',coarse_chi[0],';',coarse_chi[1]))

    if (c_sphericiy):
      avg_small = numpy.array([])
      avg_coarse = numpy.array([])
      for i in range(0,len(bodies),1):
        if bodies[i].shape == 'PLAN':
         continue
        else:
         area = 0.
         for f in bodies[i].connectivity:
           u = [bodies[i].vertex[f[1]-1][0]-bodies[i].vertex[f[0]-1][0], bodies[i].vertex[f[1]-1][1]-bodies[i].vertex[f[0]-1][1], bodies[i].vertex[f[1]-1][2]-bodies[i].vertex[f[0]-1][2]]
           v = [bodies[i].vertex[f[2]-1][0]-bodies[i].vertex[f[0]-1][0], bodies[i].vertex[f[2]-1][1]-bodies[i].vertex[f[0]-1][1], bodies[i].vertex[f[2]-1][2]-bodies[i].vertex[f[0]-1][2]]
           area += 1./2.*numpy.linalg.norm(numpy.cross(u,v))
         psi = (numpy.pi**(1./3.))*((6*4./3.*numpy.pi*bodies[i].radius**3)**(2./3.))/area
         if bodies[i].vertices == small_vertices:
           avg_small = numpy.append(avg_small, psi)
         else:
           avg_coarse = numpy.append(avg_coarse, psi)
      small_chi = numpy.mean(avg_small)
      coarse_chi = numpy.mean(avg_coarse)
      result_file.write('%s \t %10.4E \n'%('small particles sphericity=',small_chi))
      result_file.write('%s \t %10.4E \n'%('coarse particles sphericity=',coarse_chi))

    if (c_angularity):
      avg_small = numpy.array([])
      avg_coarse = numpy.array([])
      for i in range(0,len(bodies),1):
        alpha_mean = numpy.array([])
        if bodies[i].shape == 'PLAN':
         continue
        else:
         for f in range(len(bodies[i].connectivity)):
           for c in bodies[i].connectivity[(f+1):]:
             common = numpy.array([])
             sum = numpy.array([])
             for o in bodies[i].connectivity[f] :
               sum = numpy.append(sum, o)
             for el in c :
               if el in sum :
                 common = numpy.append(common, el)
               else :
                 sum = numpy.append(sum, el)
             print(sum)
             print(common)
             excl = numpy.setdiff1d(sum,common)
             print(excl)
             if len(common)==2 :
               u = numpy.array([bodies[i].vertex[common[1]-1][0]-bodies[i].vertex[common[0]-1][0], bodies[i].vertex[common[1]-1][1]-bodies[i].vertex[common[0]-1][1], bodies[i].vertex[common[1]-1][2]-bodies[i].vertex[common[0]-1][2]])
               v1 = numpy.array([bodies[i].vertex[excl[0]-1][0]-bodies[i].vertex[common[0]-1][0], bodies[i].vertex[excl[0]-1][1]-bodies[i].vertex[common[0]-1][1], bodies[i].vertex[excl[0]-1][2]-bodies[i].vertex[common[0]-1][2]])
               v2 = numpy.array([bodies[i].vertex[excl[1]-1][0]-bodies[i].vertex[common[0]-1][0], bodies[i].vertex[excl[1]-1][1]-bodies[i].vertex[common[0]-1][1], bodies[i].vertex[excl[1]-1][2]-bodies[i].vertex[common[0]-1][2]])
               n1 = numpy.cross(u,v1)
               n2 = numpy.cross(u,v2)
               alpha = numpy.arcsin(numpy.linalg.norm(numpy.cross(n1,n2))/(numpy.linalg.norm(n1)*numpy.linalg.norm(n2)))
               alpha_mean = numpy.append(alpha_mean, alpha)
             else :
               continue
        mean_ang = numpy.mean(alpha_mean)
        if bodies[i].vertices == small_vertices:
          avg_small = numpy.append(avg_small, mean_ang)
        else:
          avg_coarse = numpy.append(avg_coarse, mean_ang)
      small_alpha = numpy.mean(avg_small)
      coarse_alpha = numpy.mean(avg_coarse)
      result_file.write('%s \t %10.4E \n'%('small particles angularity=',small_alpha))
      result_file.write('%s \t %10.4E \n'%('coarse particles angularity=',coarse_alpha))


  #if (c_packing_frac):
  #if (c_wall_forces): wall_forces(ii,last_frame,bodies_current,contacts)
  #if (c_compacity): compacity(ii,last_frame,bodies_current,contacts)
  #
  result_file.write('%s \n'%('#################################################'))
  result_file.write('\n')
result_file.close()
#

