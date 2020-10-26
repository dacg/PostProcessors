#################################################
############### Post pro Builder ################
############### D. Cantor 24/02/15 ##############
#################################################

import os,sys

import numpy as np
import math

from pylmgc90.pre_lmgc import *

#################### Input ######################

dp = 2600.                                # Density of particles
dw = 1000.                                # Density of walls

# Reading the input file
file_name = "./POSTPRO/VORONOI/0main"
f = open(file_name, "r")
n_particles = int(f.readline())
f.close()

# Reading the support file
file_name = "./POSTPRO/VORONOI/0support"
f = open(file_name, "r")

# Reading the center of the particles
coor_p = np.array([[0.,0.,0.]]*n_particles)
coor_w = np.array([[0.,0.,0.]]*6)
axes_w = np.array([[0.,0.,0.]]*6)

for i in xrange(0,n_particles,1):
   cur_line = f.readline()
   cur_line = cur_line.split("   ")
   
   coor_p[i,0] = float(cur_line[0])
   coor_p[i,1] = float(cur_line[1])
   coor_p[i,2] = float(cur_line[2])
#
# Reading the center of the walls. Note they are supposed 6 walls
for i in xrange(0,6,1):
   cur_line = f.readline()
   cur_line = cur_line.split("   ")
   
   coor_w[i,0] = float(cur_line[0])
   coor_w[i,1] = float(cur_line[1])
   coor_w[i,2] = float(cur_line[2])
#

# Reading the axes of the walls.
for i in xrange(0,6,1):
   cur_line = f.readline()
   cur_line = cur_line.split("   ")
   
   axes_w[i,0] = float(cur_line[0])
   axes_w[i,1] = float(cur_line[1])
   axes_w[i,2] = float(cur_line[2])
#

f.close()

################### Building ####################

# on se place en 3D
dim = 3

# creation des conteneurs
#   * pour les corps
bodies = avatars()
#   * pour les materiaux
mat = materials()

# creation de deux materiaux :
#   * Particles
plex = material(name='PLEXx', type='RIGID', density=dp)
mat.addMaterial(plex)
#   * Walls
tdur = material(name='TDURx', type='RIGID', density=dw)
mat.addMaterial(tdur)

# creation d'un modele rigide 3D
mod = model(name='rigid', type='MECAx', element='Rxx3D', dimension=dim)

# Compilation of the fortran subroutines and tessellation code
# Note that the fortran subroutine is used to create the particles. This is an Emilien's
# code. On the other hand, the tessellation subroutine is written in c++ and uses 
# the package voro++

# Now, it is necessary to read and manipulate the output of the tessellation in order...
# to create the sub-particles

# It will be manipulated each particle in the box
for i in xrange(0,n_particles,1):
      
   # Opening the corresponding output file for this particle
   file_name = "./POSTPRO/VORONOI/vertices" + str(i+1)
   f = open(file_name, "r")
   
   # Finding the final number of subparticles after the tessellation
   nfi_sub_particles = 0
   
   file_temp = "./POSTPRO/VORONOI/post_pro_nei" + str(i+1)
   ftemp = open(file_temp, "r")
   while (ftemp.readline()):
      nfi_sub_particles += 1;
   #
   ftemp.close()
   
   # Reading the information of each sub-particle
   for j in xrange(0,nfi_sub_particles,1):
      
      ac_line = f.readline()
      cur_data = ac_line.split(", ")
      
      # Getting the number of vertices of the sub-particle
      string_n_vertices = cur_data[2]
      n_vertices = int(string_n_vertices[11:])
      
      # Array to save the vertices of the new sub-particle
      sub_vertices = np.array([[0.,0.,0.]])
      
      # Getting the vertices of the sub-particle
      string_sub_vertices = cur_data[3]
      string_sub_vertices = string_sub_vertices[10:(len(string_sub_vertices)-2)]
      
      string_sub_vertices = string_sub_vertices.split(") (")
      
      # Saving each vertex in the array and translating it to the global frame
      # Also, it is performed a filter to ensure that they won't be repeated. 
      # The reason to do so is that sometimes the output of voro++ repeats vertices. 
      # This is done in the local frame
      for k in xrange(0,n_vertices,1):
         s_vert = string_sub_vertices[k].split(",")
         # Writing the first vertex
         if (k==0):
            sub_vertices[0,0] = float(s_vert[0])
            sub_vertices[0,1] = float(s_vert[1])
            sub_vertices[0,2] = float(s_vert[2])
         else:
            # Appending the rest
            sub_vertices = np.append(sub_vertices, [[float(s_vert[0]),
                                                     float(s_vert[1]),
                                                     float(s_vert[2])]], axis=0)
            #
         #
      #
      
      # Translating the sub_vertices to the global frame
      for k in xrange(0,sub_vertices.shape[0],1):
         sub_vertices[k,0]= sub_vertices[k,0] + coor_p[i,0]
         sub_vertices[k,1]= sub_vertices[k,1] + coor_p[i,1]
         sub_vertices[k,2]= sub_vertices[k,2] + coor_p[i,2]
      #
      
      # Building the sub-particle
      body = rigidPolyhedron(model=mod, material=plex, color='BLEUx', 
                             generation_type='vertices',vertices=sub_vertices)
      # Adding the sub-particle to the container 
      bodies += body
      
      # Cleaning varialbles
      del sub_vertices
      #
   #
   f.close()
#

# creation de corps pour les parois
down = rigidPlan(axe1=axes_w[0,0], axe2=axes_w[0,1], axe3=axes_w[0,2],
                 center=numpy.array([coor_w[0,0], coor_w[0,1], coor_w[0,2]]) ,
                 model=mod,material=tdur,color='VERTx')

left = rigidPlan(axe1=axes_w[1,0], axe2=axes_w[1,1], axe3=axes_w[1,2],
                 center=numpy.array([coor_w[1,0], coor_w[1,1], coor_w[1,2]]),
                 model=mod,material=tdur,color='VERTx')

right= rigidPlan(axe1=axes_w[2,0], axe2=axes_w[2,1], axe3=axes_w[2,2],
                 center=numpy.array([coor_w[2,0], coor_w[2,1], coor_w[2,2]]),
                 model=mod,material=tdur,color='VERTx')

front= rigidPlan(axe1=axes_w[3,0], axe2=axes_w[3,1], axe3=axes_w[3,2],
                 center=numpy.array([coor_w[3,0], coor_w[3,1], coor_w[3,2]]),
                 model=mod,material=tdur,color='VERTx')

rear = rigidPlan(axe1=axes_w[4,0], axe2=axes_w[4,1], axe3=axes_w[4,2],
                 center=numpy.array([coor_w[4,0], coor_w[4,1], coor_w[4,2]]),
                 model=mod,material=tdur,color='VERTx')

top  = rigidPlan(axe1=axes_w[5,0], axe2=axes_w[5,1], axe3=axes_w[5,2],
                 center=numpy.array([coor_w[5,0], coor_w[5,1], coor_w[5,2]]),
                 model=mod,material=tdur,color='VERTx')

# on tourne les parois formant les cotes
# rotation autour de l'axe x, par rapport au centre d'inertie, avec un angle
# -pi/2 (parametres: angles d'Euler)
left.rotate(theta=-0.5*math.pi, center=left.nodes[1].coor)
# rotation autour de l'axe x, par rapport au centre d'inertie, avec un angle
# pi/2 (parametres: angles d'Euler)
right.rotate(theta=0.5*math.pi, center=right.nodes[1].coor)
# rotation autour de l'axe y, par rapport au centre d'inertie, avec un angle
# -pi/2 (parametres: axe + angle)  
front.rotate(type='axis', alpha=-0.5*math.pi, axis=[0., 1., 0.], center=front.nodes[1].coor)
# rotation autour de l'axe y, par rapport au centre d'inertie, avec un angle
# pi/2 (parametres: axe + angle)
rear.rotate(type='axis', alpha=0.5*math.pi, axis=[0., 1., 0.], center=rear.nodes[1].coor)

# blocage des parois
down.imposeDrivenDof(component=[1, 2, 3, 4, 5, 6], dofty='vlocy')
left.imposeDrivenDof(component=[1, 2, 3, 4, 5, 6], dofty='vlocy')
right.imposeDrivenDof(component=[1, 2, 3, 4, 5, 6], dofty='vlocy')
front.imposeDrivenDof(component=[1, 2, 3, 4, 5, 6], dofty='vlocy')
rear.imposeDrivenDof(component=[1, 2, 3, 4, 5, 6], dofty='vlocy')
top.imposeDrivenDof(component=[1, 2, 4, 5, 6], dofty='vlocy')

# ajouts de parois au conteneur de corps
bodies.addAvatar(down)
bodies.addAvatar(left)
bodies.addAvatar(right)
bodies.addAvatar(front)
bodies.addAvatar(rear)
bodies.addAvatar(top)

os.system("mkdir ./POSTPRO/VORONOI/DATBOX")
# ecriture des fichiers
writeBodies(bodies, chemin='./POSTPRO/VORONOI/DATBOX/')
writeDrvDof(bodies, chemin='./POSTPRO/VORONOI/DATBOX/')
writeDofIni(bodies, chemin='./POSTPRO/VORONOI/DATBOX/')

# Pre visualization of bodies
visuAvatars(bodies)                   #
