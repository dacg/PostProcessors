#################################################
###############    Command 3D   #################
############### D. Cantor 21/01/15 ##############
#################################################

# importing modules
from pylmgc90.chipy import *
import numpy as np

Initialize()                                # Initializing
checkDirectories()                          # checking/creating mandatory subfolders
utilities_DisableLogMes()                   # logMes

################ Variables #################

p_step = 10                                 # Print time step every p_step

# time evolution parameters
dt = 1e-6                                   # Time step
nb_steps = 2000                             # Steps

theta = 0.5                                 # Theta integrator parameter

deformable = 0                              # Deformable  yes=1, no=0

freq_detect = 1                             # Detection frequency
Rloc_tol = 5.e-2                            # Ask!!!!!!!!!!!!!!!!!!!!!!!!

PRPRx_UseCpCundallDetection(40)            # Detection method. Arg: # of iteretions
#PRPRx_UseCpF2fExplicitDetection(1e-2)

PRPRx_LowSizeArrayPolyr(1000)              # Memory parameter

# The 'correct' working of the Mohr Coulomb law need the following functions
PRPRx_ShrinkPolyrFaces(0.01)
PRPRx_SetVoronoiInitialFaces()
POLYR_SetVoronoiInitialFaces()
nlgs_3D_DiagonalResolution()
nlgs_SetWithQuickScramble()
#nlgs_3D_ScrambleContactOrder()


# nlgs parameters
type='Stored_Delassus_Loops         '
norm = 'QM/16'
tol = 1e-4
relax = 1.0
gs_it1 = 51
gs_it2 = 501

freq_write   = 5                       # Writing frequency

freq_display = 5                       # Display frequency
ref_radius = 1.0e-3

# Triaxial parameters
#sigma_wall = 100                          # Stress on the walls

# 1=mur down, 2=mur right, 3=mur up, 4=mur left, 5=mur front, 6=mur rear
#tri_loads = np.array([3,sigma_wall,4,sigma_wall,5,sigma_wall])

################ Simulation ################

# Set space dimension
SetDimension(3)

# read and load
#
utilities_logMes('INIT TIME STEPPING')
TimeEvolution_SetTimeStep(dt)
Integrator_InitTheta(theta)
#
utilities_logMes('READ BEHAVIOURS')
ReadBehaviours()
if deformable: ReadModels()
#
utilities_logMes('READ BODIES')
ReadBodies()
#
utilities_logMes('LOAD BEHAVIOURS')
LoadBehaviours()
if deformable: LoadModels()
#
utilities_logMes('READ INI DOF')
ReadIniDof()
#
if deformable:
  utilities_logMes('READ INI GPV')
  ReadIniGPV()
#
utilities_logMes('READ DRIVEN DOF')
ReadDrivenDof()
#
utilities_logMes('LOAD TACTORS')
LoadTactors()
#
utilities_logMes('READ INI Vloc Rloc')
ReadIniVlocRloc()

# paranoid writes

utilities_logMes('WRITE BODIES')
WriteBodies()
utilities_logMes('WRITE BEHAVIOURS')
WriteBehaviours()
utilities_logMes('WRITE DRIVEN DOF')
WriteDrivenDof()

#
# open display & postpro
#

utilities_logMes('DISPLAY & WRITE')
OpenDisplayFiles()
OpenPostproFiles()

# since constant, compute elementary mass once
ComputeMass()

nb_rbdy2 = RBDY3_GetNbRBDY3()
wall_down  = nb_rbdy2 - 5
wall_left  = nb_rbdy2 - 4
wall_right = nb_rbdy2 - 3
wall_front = nb_rbdy2 - 2
wall_rear  = nb_rbdy2 - 1
wall_up   = nb_rbdy2

iteration = 0                               # Time steping counter

# time loop
for k in xrange(1,nb_steps+1):
  #
  iteration = iteration + 1
    
  if (iteration == p_step):
    print "Step: " + str(k)
    iteration = 0
  #

  IncrementStep()
  # bulk_behav_SetGravity([x,y,z])
  ComputeFext()
  
  #RBDY3_TriaxialLoading(wall_down, wall_right, wall_up, wall_left, wall_front, wall_rear, 3, tri_loads)
  
  ComputeBulk()
  ComputeFreeVelocity()

  SelectProxTactors(freq_detect)

  RecupRloc(Rloc_tol)

  ExSolver(type, norm, tol, relax, gs_it1, gs_it2)
  UpdateTactBehav()

  StockRloc()

  ComputeDof()

  UpdateStep()

  WriteOutDof(freq_write)
  WriteOutVlocRloc(freq_write)

  WriteDisplayFiles(freq_display,ref_radius)
  WritePostproFiles()

#
# close display & postpro
#
CloseDisplayFiles()
ClosePostproFiles()

# this is the end
Finalize()

os.system("say ya")
