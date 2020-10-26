#################################################
#################### Builder ####################
############### D. Cantor 07/11/14 ##############
#################################################
import os

# Deleting the folders DISPLAY, OUTBOX, and POSTPRO
if (os.path.exists("DISPLAY")):
   os.system("rm -rf ./DISPLAY ./OUTBOX ./POSTPRO")
#

#from chipy import *
from pylmgc90.chipy import *
from numpy import *

checkDirectories()

#################### Inputs #####################

overall_DIME(3,0)                         # A 3D simualtion
utilities_DisableLogMes()                 # Desactivation des messages de log

### computation's parameters definition ### 
utilities_logMes('INIT TIME STEPPING')

dt = 1.e-4                                # time step length
theta = 0.5                               # value of the parameter of the theta-method
nb_steps = 100                              # number of time steps

echo = 0                                  # bavardage de certaines fonctions

p_step = 1                                # Print time step every p_step

### parameters setting ###
freq_detect = 1                           # Detection frequency

#POLYR_SkipAutomaticReorientation()       # Deja de molestar con particulas no convexas
#PRPRx_UseNcF2fExplicitDetection(1.,1e-3) 

#PRPRx_WithNodalContact()
#PRPRx_UseCpF2fExplicitDetection(0.1)
PRPRx_UseCpCundallDetection(30)          # Detection method. OJO preguntar el arg
PRPRx_LowSizeArrayPolyr(100)              # When the polyhedron has a lot of sides
                                          # it helps to control some parameter
                                          # of arrays.
                                          
PRPRx_SetVoronoiInitialFaces()
POLYR_SetVoronoiInitialFaces()

freq_display = 1                         # Visualization frequency
freq_write = 1                           # History files frequency

#   * nlgs solver parameters
type = 'Stored_Delassus_Loops         '
norm = 'QM/16'
tol = 0.1666e-3
relax = 1.0
gs_it1 = 51
gs_it2 = 501

################### Running #####################

#
TimeEvolution_SetTimeStep(dt)
Integrator_InitTheta(theta)
#

### model reading ###
utilities_logMes('READ BODIES')
RBDY3_ReadBodies()
#
PLANx_LoadTactors()

utilities_logMes('READ INI DOF')
TimeEvolution_ReadIniDof()
RBDY3_ReadIniDof()

TimeEvolution_WriteLastDof()
RBDY3_WriteLastDof()

POLYR_LoadTactors()
#
overall_WriteBodies()
RBDY3_WriteBodies()
#
utilities_logMes('READ BEHAVIOURS')
bulk_behav_ReadBehaviours()
tact_behav_ReadBehaviours()

#LOADS
RBDY3_LoadBehaviours()

bulk_behav_WriteBehaviours()
tact_behav_WriteBehaviours()

utilities_logMes('READ INI Vloc Rloc')
TimeEvolution_ReadIniVlocRloc()
PRPRx_ReadIniVlocRloc()
PRPLx_ReadIniVlocRloc()

utilities_logMes('READ DRIVEN DOF')
RBDY3_ReadDrivenDof()

overall_WriteDrivenDof()
RBDY3_WriteDrivenDof()

### compute masses ###
RBDY3_ComputeMass()

### post3D ##
#                         1234567890123456
# definition of fields to be computed by the post3D module
post3D_SetDisplayedField('POSITION')
post3D_SetDisplayedField('AVERAGE VELOCITY')
post3D_SetDisplayedField('STRESS')
# initialization of the post3D module
post3D_Init()

# definition of another fields to be displayed by the display_3D module
display_3D_SetDisplayedField('TACTOR')
display_3D_SetDisplayedField('INTERACTION')
# choose of the file format for the visualiztion
display_3D_SetDisplayFileFormat('VTK')
# initilization of the display_3D module
display_3D_Init(0)

# compute of a first visualization
display_3D_WriteOutDisplayFile(0)

# Counter
iteration = 0
# time loop
for k in range(1, nb_steps + 1):
    
    iteration = iteration + 1
    
    if (iteration == p_step):
       print "Step: " + str(k)
       iteration = 0
    #
    utilities_logMes('itere : '+str(k))
    #
    utilities_logMes('INCREMENT STEP')
    TimeEvolution_IncrementStep()
    RBDY3_IncrementStep()
    #
    utilities_logMes('COMPUTE Fext')
    RBDY3_ComputeFext()
    #
    utilities_logMes('COMPUTE Fint')
    RBDY3_ComputeBulk()
    # 
    utilities_logMes('COMPUTE Free Vlocy')
    RBDY3_ComputeFreeVelocity()
    #
    utilities_logMes('SELECT PROX TACTORS')
    overall_SelectProxTactors(freq_detect)
    
    PRPRx_SelectProxTactors()
    PRPLx_SelectProxTactors()
    #
    PRPRx_RecupRloc()
    PRPLx_RecupRloc()
    
    utilities_logMes('RESOLUTION' )
    
    #utilities_logMes(str(type)+', '+str(norm)+', '+str(tol)+', '+str(relax)+', '+str(gs_it1)+', '+str(gs_it2))
    nlgs_3D_ExSolver(type, norm, tol, relax, gs_it1, gs_it2)
    
    PRPRx_StockRloc()
    PRPLx_StockRloc()
    #
    utilities_logMes('COMPUTE DOF')
    RBDY3_ComputeDof()
    #
    utilities_logMes('UPDATE DOF')
    TimeEvolution_UpdateStep()
    RBDY3_UpdateDof()
    
    ### post3D ###
    post3D_Update()
    
    overall_WriteOutDisplayFile(freq_display)
    display_3D_WriteOutDisplayFile(0)
    
    TimeEvolution_WriteOutDof(freq_write)
    RBDY3_WriteOutDof(-1,9999999)
    
    TimeEvolution_WriteOutVlocRloc(freq_write)
    PRPRx_WriteOutVlocRloc()
    PRPLx_WriteOutVlocRloc()
    
    ### writeout handling ###
    overall_CleanWriteOutFlags()
    
#
TimeEvolution_WriteLastDof()
RBDY3_WriteLastDof()

os.system("say ya")