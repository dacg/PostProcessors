
from pylmgc90.chipy import *

checkDirectories()

# desactivation des messages de log
utilities_DisableLogMes()

####
# info gestion du temps
dt = 1.e-4
theta = 0.5
nb_steps = 120000

# bavardage de certaines fonctions
echo = 0

p_step = 1000                               # Print time step every p_step

# info generation fichier visu
freq_display = 20000
freq_write = 20000

# info contact
freq_detect = 5

#       123456789012345678901234567890
type = 'Stored_Delassus_Loops         '
norm = 'Quad '
tol = 0.1666e-3
relax = 1.0
gs_it1 = 51
gs_it2 = 1001

### definition des parametres du calcul ### 
utilities_logMes('INIT TIME STEPPING')
TimeEvolution_SetTimeStep(dt)
Integrator_InitTheta(theta)

### lecture du modele ###

### model reading ###
utilities_logMes('READ BODIES')
RBDY2_ReadBodies()

utilities_logMes('READ INI DOF')
TimeEvolution_ReadIniDof()
RBDY2_ReadIniDof()

utilities_logMes('READ BEHAVIOURS')
bulk_behav_ReadBehaviours()
tact_behav_ReadBehaviours()

#LOADS
DISKx_LoadTactors()
JONCx_LoadTactors()
RBDY2_LoadBehaviours()

utilities_logMes('READ INI Vloc Rloc')
TimeEvolution_ReadIniVlocRloc()
DKDKx_ReadIniVlocRloc()
DKJCx_ReadIniVlocRloc()

utilities_logMes('READ DRIVEN DOF')
RBDY2_ReadDrivenDof()

### ecriture paranoiaque du modele ###
utilities_logMes('WRITE BODIES')
overall_WriteBodies()
RBDY2_WriteBodies()

utilities_logMes('WRITE BEHAVIOURS')
bulk_behav_WriteBehaviours()
tact_behav_WriteBehaviours()

utilities_logMes('WRITE DRIVEN DOF')
overall_WriteDrivenDof()
RBDY2_WriteDrivenDof()

### post2D ##
post2D_SetDisplayedField('CONTACT POINT   ')
post2D_SetDisplayedField('TACTOR          ')
post2D_Init()

utilities_logMes('COMPUTE MASS')
RBDY2_ComputeMass()

# Counter
iteration = 0

for k in xrange(1, nb_steps + 1, 1):
   
   iteration = iteration + 1
   
   if (iteration == p_step):
      print "Step: " + str(k)
      iteration = 0
   #
   
   #
   utilities_logMes('itere : '+str(k))
   #
   utilities_logMes('INCREMENT STEP')
   TimeEvolution_IncrementStep()
   RBDY2_IncrementStep()

   utilities_logMes('DISPLAY TIMES')
   TimeEvolution_DisplayStep()

   utilities_logMes('COMPUTE Fext')
   RBDY2_ComputeFext()

   utilities_logMes('COMPUTE Fint')
   RBDY2_ComputeBulk()
   
   utilities_logMes('COMPUTE Free Vlocy')
   RBDY2_ComputeFreeVelocity()
   #
   utilities_logMes('SELECT PROX TACTORS')
   overall_SelectProxTactors(freq_detect)
   DKJCx_SelectProxTactors()
   DKDKx_SelectProxTactors()
   #
   DKJCx_RecupRloc()
   DKDKx_RecupRloc()
   nlgs_ExSolver(type, norm, tol, relax, gs_it1, gs_it2)
   DKJCx_StockRloc()
   DKDKx_StockRloc()
   #
   utilities_logMes('COMPUTE DOF')
   RBDY2_ComputeDof()
   #
   utilities_logMes('UPDATE DOF')
   TimeEvolution_UpdateStep()
   RBDY2_UpdateDof()
   #
   utilities_logMes('WRITE LAST DOF')
   TimeEvolution_WriteLastDof()
   RBDY2_WriteLastDof()
   #
   utilities_logMes('WRITE LAST Vloc Rloc')
   TimeEvolution_WriteLastVlocRloc()
   DKDKx_WriteLastVlocRloc()
   DKJCx_WriteLastVlocRloc()
   #
   ### post2D ###
   overall_WriteOutDisplayFile(freq_display)
   post2D_WriteOutDisplayFile(0)
   
   TimeEvolution_WriteOutDof(freq_write)
   RBDY2_WriteOutDof() 
   TimeEvolution_WriteOutVlocRloc(freq_write)
   DKDKx_WriteOutVlocRloc()
   DKJCx_WriteOutVlocRloc()
   
   ### wrtieout handling ###
   overall_CleanWriteOutFlags()
