!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!  POST-PROCESSOR !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!! SPHERES !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


program post3D_spheres

implicit none

type  ::  T_CORPS
  real(kind=8)                                  ::  Rmax,volume,Ax1,Ax2,Ax3
  real(kind=8),dimension(3)                     ::  center,center_ref
  real(kind=8),dimension(:,:),pointer           ::  vertex,vertex_ref
  logical                                       ::  inbox
  real(kind=8)                                  ::  Vx,Vy,Vz,Vrx,Vry,Vrz
  real(kind=8)                                  ::  Wx,Wy,Wz
  real(kind=8)                                  ::  Norme_V,Norme_Vr,I1,I2,I3
  real(kind=8),dimension(3,3)                   ::  momentI,Rot,moment_part,moment_part1,moment_part2,&
                                                    moment_part3,moment_part4,moment_part5,moment_part6,moment_part7,&
                                                    moment_part8,moment_part9
  real(kind=8)                                  ::  nctc,nctcslide,nctcstick,nb_ctc,Norm_F,ctc_float, av_fn, av_mob
  character(len=5)                              ::  behav,color
  integer                                       ::  id
  real(kind=8),dimension(3,3)                   ::  Fabric,Fabric_N,Fabric_T,Fabric_T2

end type T_CORPS

type  ::  T_CONTACT
  integer                                        ::  icdent,ianent,type
  real(kind=8),dimension(3)                      ::  n,t,s
  real(kind=8),dimension(3)                      ::  coor_ctc
  real(kind=8)                                   ::  rn,rt,rs,vln,vls,vlt,gap
  integer                                        ::  sect
  character(len=5)                               ::  nature,statut
  logical                                        ::  inbox, contact_perdu
end type T_CONTACT

type(T_CORPS),dimension(:),allocatable          ::  TAB_SPHERE,TAB_PLAN
type(T_CONTACT),dimension(:),allocatable        ::  TAB_CONTACT

!DEFINITIONS Vloc-Dof-Bodies
!================================================================
integer                                         ::  n_walls=0, n_particles=0
integer                                         ::  compteur_clout=0
integer                                         ::  nb_face,nb_sommet
logical                                         ::  the_first_time=.true., fin_post=.false. , first_bodies=.true.
integer                                         ::  step,nb_ligneCONTACT
real(kind=8)                                    ::  time,Mean_total_Normal_force

!DEFINITIONS BOITE
!================================================================
real(kind=8),dimension(3,8)                     ::  tab_sommet_box
real(kind=8)                                    ::  Def,H_ini,Long_ini,Larg_ini,Def_H,step_pdf_start,step_pdf_stop,&
                                                    step_ksi_reseau,stop_ksi_reseau,Ep_avant=0,Eq_avant=0,Ep,Eq,Depsilon_p, &
                                                    Depsilon_q, pressure_in_sample=0,pressure_c_in_sample = 0.D0,&
                                                    sin_de_phi=0,partial_pressure_in_sample=0
logical                                         ::  deformation_ini = .true.,deformation_clout=.false.

!MOYENNE POUR KSI-RESEAU
!================================================================
real(kind=8),dimension(:,:,:),allocatable       ::  Mean_Ksi_reseau
integer                                         ::  pas_ksi = 0
integer                                         ::  quel_pas_de_temps =0

!CALCUL DES PROFILS
!================================================================
integer                                         ::  debut_profils=0,fin=0
logical                                         ::  the_first_time_profils

!CALCUL DU PROFIL DE VITESSE MOYEN
!================================================================
type  ::  T_step_velocity
  real(kind=8),dimension(:,:),allocatable       ::  tab_vmoyen
  real(kind=8),dimension(:,:),allocatable       ::  tab_vrmoyen
end type T_step_velocity

type(T_step_velocity),dimension(:),allocatable  ::  TAB_vitesse

!CALCUL DU PROFIL DE CONTRAINTE MOYEN
!================================================================
type  ::  T_sigma
  real(kind=8),dimension(3,3)                   :: sigma
end type
type  ::  T_step_contrainte
  type(T_sigma),dimension(:),allocatable        ::  tab_sigma
end type

type(T_step_contrainte),dimension(:),allocatable  ::  TAB_contrainte

!CALCUL DU PROFIL DE COMPACITY MOYEN
!================================================================
type  ::  T_step_compacity
  real(kind=8),dimension(:,:),allocatable          ::  tab_compacity_z
end type T_step_compacity
type(T_step_compacity),dimension(:),allocatable    ::  TAB_compacity


!----------------Pour Persiste Contact---------------------------
logical                                         :: first_list = .true.
type(T_CONTACT),dimension(:),pointer            :: TAB_CONTACT_INITIAL
integer                                         :: nb_ligneCONTACT_first=0,debut_persistance=0


!COMMANDES D APPELLE
!================================================================
character(len=30)                               ::  command
real(kind=8)                                    ::  PI = 3.14159265359, &
                                                    width,large,height
integer                                         ::  type_box=0,calcul_vitesse_moyenne=0,calcul_coordination=0, &
                                                    calcul_qoverp=0,calcul_compacity=0, &
                                                    calcul_anisotropy_contact=0, calcul_anisotropy_force=0, &
                                                    calcul_anisotropy_branch=0, &
                                                    calcul_construct_bodies=0, calcul_granulo=0, &
                                                    calcul_walls_positions=0, calcul_walls_forces=0, &
                                                    calcul_contact_dir=0, c_n_ctc_probability=0, &
                                                    c_float_particles=0, c_coordination_class=0, &
                                                    calcul_persiste_contacts=0, c_draw=0, &
                                                    c_list_n_forces=0, c_list_branches=0, c_list_mobilization=0, &
                                                    c_fn_class=0, c_mobilization_class=0, c_ctc_dir=0, &
                                                    c_brc_dir=0, c_frc_dir=0, c_aniso_class=0, c_aniso_part=0, &
                                                    c_list_t_forces=0, c_mobilization = 0, c_radial_dist = 0, &
                                                    c_segregation_index = 0, c_lacey_index=0, c_list_active_par=0, &
                                                    c_float_class=0, c_prop_nctc=0, c_qoverp_class=0, c_stresses=0, &
                                                    c_vloc_profile=0, c_pdf_z=0


! Variables David C.
logical                                          ::  first_over_all = .true.
integer                                          ::  qoverp_option, n_layers

!================================================================
!Lecture du fichier de commande
!================================================================
print*,'-------------------------------------------'
print*,'!                                         !'
print*,'!      Module de Post-traitement 2D       !'
print*,'!         from Emilien Azema Code         !'
print*,'!                 David C.                !'
print*,'-------------------------------------------'
print*,'Reading INPUT'

open(unit=1,file='POST_INPUT.DAT',status='old')
! Reading the commands
do

  read(1,'(A30)') command
  print*,command

  if (command=='VLOC_DOF_BODIES              :') then
    read(1,*) n_walls
    read(1,*) n_particles
    cycle
  end if

  if (command=='COORDINATION                 :') then
    ! Uses port 100
    calcul_coordination=1
    open (unit=100,file='./POSTPRO/COORDINATION.DAT',status='replace')
    cycle
  end if

  if (command=='GRAIN SIZE DISTRIB           :') then
    calcul_granulo=1
    cycle
  end if

  if (command=='VITESSE MOYENNE              :') then
    ! Uses port 101
    calcul_vitesse_moyenne=1
    open (unit=101,file='VITESSE_MOYENNE.DAT',status='replace')
    cycle
  end if

  if (command=='COMPACITY                    :') then
    ! Uses port 102
    calcul_compacity=1
    open (unit=102,file='./POSTPRO/COMPACITY.DAT',status='replace')
    cycle
  end if

  if (command=='QoverP                       :') then
    calcul_qoverp=1
    open (unit=103,file='./POSTPRO/QoverP.DAT',status='replace')
    cycle
  end if

  if (command=='ANISOTROPY CONTACT           :') then
    calcul_anisotropy_contact=1
    open (unit=104,file='./POSTPRO/ANISOTROPY_CONTACT.DAT',status='replace')
    cycle
  end if

  if (command=='ANISOTROPY FORCE             :') then
    calcul_anisotropy_force=1
    open (unit=105,file='./POSTPRO/ANISOTROPY_FORCE.DAT',status='replace')
    cycle
  end if

  if (command=='ANISOTROPY BRANCH            :') then
    calcul_anisotropy_branch=1
    open (unit=106,file='./POSTPRO/ANISOTROPY_BRANCH.DAT',status='replace')
    cycle
  end if

  if (command=='WALLS POSITION               :') then
    calcul_walls_positions=1
    open (unit=108,file='./POSTPRO/WALLSPOS.DAT',status='replace')
    cycle
  end if

  if (command=='WALLS FORCES                 :') then
    calcul_walls_forces=1
    open (unit=109,file='./POSTPRO/WALLSFORCE.DAT',status='replace')
    cycle
  end if

  !if (command=='CONTACTS ORIENTATION         :') then
  !  calcul_contact_dir=1
  !  cycle
    ! Uses port 110
  !end if

  if (command=='CONTACT NUMBER PROBABILITY   :') then
    c_n_ctc_probability = 1
    open (unit=111,file='./POSTPRO/CTC_PROBABILITY.DAT',status='replace')
    cycle
  end if

  if (command=='FLOATING PARTICLES           :') then
    c_float_particles = 1
    open (unit=112,file='./POSTPRO/FLOAT_PARTICLES.DAT',status='replace')
    cycle
  end if

  if (command=='COORDINATION CLASS           :') then
    c_coordination_class = 1
    ! Uses port 113
    cycle
  end if

  if (command=='DRAW                         :') then
    c_draw = 1
    ! Uses port 114
    cycle
  end if

  if (command=='LIST NORMAL FORCES           :') then
    c_list_n_forces = 1
    ! Uses port 115
    cycle
  end if

  if (command=='LIST BRANCHES                :') then
    c_list_branches = 1
    ! Uses port 116
    cycle
  end if

  if (command=='LIST FRICTION MOBILIZATION   :') then
    c_list_mobilization = 1
    ! Uses port 117
    cycle
  end if

  if (command=='NORMAL FORCE CLASS           :') then
    c_fn_class = 1
    ! Uses port 118
    cycle
  end if

  if (command=='MOBILIZATION CLASS           :') then
    c_mobilization_class = 1
    ! Uses port 119
    cycle
  end if

  if (command=='CONTACTS ORIENTATION         :') then
    c_ctc_dir = 1
    ! Uses port 120
    cycle
  end if

  if (command=='BRANCH  DISTRIBUTION         :') then
    c_brc_dir = 1
    ! Uses port 121
    cycle
  end if

  if (command=='FORCES DISTRIBUTION          :') then
    c_frc_dir = 1
    ! Uses port 122
    cycle
  end if

  if (command=='ANISOTROPY CLASS             :') then
    c_aniso_class = 1
    ! Uses port 123
    cycle
  end if

  if (command=='ANISOTROPY CLASS PARTICLES   :') then
    c_aniso_part = 1
    ! Uses port 124
    cycle
  end if

  if (command=='LIST TANGENTIAL FORCES       :') then
    c_list_t_forces = 1
    ! Uses port 125
    cycle
  end if

  if (command=='MOBILIZATION                 :') then
    c_mobilization=1
    open (unit=126,file='./POSTPRO/MOBILIZATION.DAT',status='replace')
    cycle
  end if

  if (command=='RADIAL DISTRIBUTION          :') then
    c_radial_dist = 1
    ! Uses port 127
    cycle
  end if

  if (command=='SEGREGATION                  :') then
    ! Computes the covariance of particle sizes at contact
    c_segregation_index = 1
    open (unit=128,file='./POSTPRO/SEGREGATION.DAT',status='replace')
    cycle
  end if

  if (command=='LACEY                        :') then
    ! Computes the variance of particle positions pondered by volume
    ! It allows to compute the lacey index parameter for mixing
    c_lacey_index = 1
    open (unit=129,file='./POSTPRO/FOR_LACEY.DAT',status='replace')
    cycle
  end if

  if (command=='LIST ACTIVE PARTICLES        :') then
    c_list_active_par = 1
    ! Uses port 130
    cycle
  end if

  if (command=='FLOAT CLASS                  :') then
    c_float_class = 1
    ! Uses port 131
    cycle
  end if

  if (command=='NUM CONTACT PROPORTION       :') then
    ! Uses port 132
    c_prop_nctc = 1
    open (unit=132,file='./POSTPRO/PROPO_NCTC.DAT',status='replace')
    cycle
  end if

  if (command=='QOVERP CLASS                 :') then
    c_qoverp_class = 1
    ! Uses port 133
    cycle
  end if

  if (command=='STRESSES                     :') then
    c_stresses = 1
    ! Uses port 134
    open (unit=134,file='./POSTPRO/STRESSES.DAT',status='replace')
    cycle
  end if

  if (command=='VELOCITY PROFILE             :') then
    c_vloc_profile = 1
    read(1,*) n_layers
    ! Uses port 135
    cycle
  end if

  if (command=='PDF BY COORDINATION NUMBER   :') then
    c_pdf_z = 1
    ! Uses port 136
    cycle
  end if

  !if (command=='CONTACT PERSISTANCE          :') then
  !  calcul_persiste_contacts = 1
  !  read(1,*) debut_persistance
  !  open (unit=134,file='CONTACT_PERSISTANCE.DAT',status='replace')
  !  cycle
  !end if

  if (command=='END                          :') exit
end do
close(1)

!==============================================================================
!Appel des differente ressources lues dans le fichier de commande
!==============================================================================

  do
    call read_Vloc_dof_bodies

    call calcul_deformation

    ! Calling the Subroutines
    if (calcul_coordination         == 1)  call nb_coordination
    if (calcul_walls_positions      == 1)  call walls_position
    if (calcul_walls_forces         == 1)  call walls_force
    !if (calcul_granulo              == 1)  call granulometry
    !if (calcul_contact_dir          == 1)  call contact_direction
    if (calcul_vitesse_moyenne      == 1)  call vitesse_moyenne
    if (calcul_compacity            == 1)  call compacity
    if (calcul_qoverp               == 1)  call qoverp
    if (calcul_anisotropy_contact   == 1)  call anisotropy_contact
    if (calcul_anisotropy_force     == 1)  call anisotropy_force
    if (calcul_anisotropy_branch    == 1)  call anisotropy_branch
    if (c_float_particles           == 1)  call float_particles
    if (c_coordination_class        == 1)  call coordination_class
    if (c_draw                      == 1)  call draw
    if (c_list_n_forces             == 1)  call list_n_forces
    if (c_list_branches             == 1)  call list_branches
    if (c_list_mobilization         == 1)  call list_mobilization
    if (c_fn_class                  == 1)  call fn_class
    if (c_mobilization_class        == 1)  call mobilization_class
    if (c_ctc_dir                   == 1)  call ctc_dir
    if (c_brc_dir                   == 1)  call branch_dir
    if (c_frc_dir                   == 1)  call forces_dir
    if (c_aniso_class               == 1)  call aniso_class
    if (c_aniso_part                == 1)  call aniso_part
    if (c_list_t_forces             == 1)  call list_t_forces
    if (c_mobilization              == 1)  call mobilization
    if (c_radial_dist               == 1)  call radial_dist
    if (c_segregation_index         == 1)  call segregation_index
    if (c_lacey_index               == 1)  call for_lacey
    if (c_list_active_par           == 1)  call list_active_par
    if (c_float_class               == 1)  call float_class
    if (c_prop_nctc                 == 1)  call prop_nctc
    if (c_qoverp_class              == 1)  call qoverp_class
    if (c_stresses                  == 1)  call stresses
    if (c_vloc_profile              == 1)  call vloc_profile
    if (c_pdf_z                     == 1)  call pdf_z
    !if (c_n_ctc_probability         == 1)  call ctc_probability

    if (fin_post) then
      call close_all
      print*,'--> End of post-processing <--'
      exit
    endif
    first_over_all = .false.
  end do


!==============================================================================
!==============================================================================
!==============================================================================
!==============================================================================

 contains

!======================================================================
!Lecture des fichiers Vloc, Dof et BODIES et on commence  remplir notre type
!TAB_SPHER.... c'est un peu bancale je l'accorde volontier :-)
!======================================================================
subroutine read_Vloc_dof_bodies

  implicit none

  integer                              ::  i=0,j=0,n_spsp,n_sppl,num_part,cpt_type,err
  integer                              ::  nb_SPSPx,nb_SPPLx,icdent,ianent
  real(kind=8)                         ::  rn,rs,rt,t1,t2,t3,n1,n2,n3,s1,s2,s3,rayon,vls,vln,vlt,gap
  character(len=6)                     ::  text
  character(len=5)                     ::  statut,behav,color
  character(len=13)                    ::  text2
  real(kind=8)                         ::  center1,center2,center3,ax1r,ax2r,ax3r,coor1,coor2,coor3,&
                                           centercluster1,centercluster2,centercluster3
  real(kind=8)                         ::  coor_ctc1,coor_ctc2,coor_ctc3
  real(kind=8)                         ::  som1,som2,som3,Volume
  real(kind=8),dimension(3)            ::  X,Y,Z,center_gravity
  real(kind=8),dimension(3,3)          ::  rot
  real(kind=8)                         ::  Vrx,Vry,Vrz,vx,vy,vz
  real(kind=8)                         ::  I1,I2,I3,Number_mean_Total_force
  character(len=22)                    ::  clout_DOF
  character(len=28)                    ::  clout_Vloc
  character(len=20)                    ::  clout_Bodies
  
  clout_Bodies = './OUTBOX/BODIES.OUT'
  clout_DOF  =  './OUTBOX/DOF.OUT.    '
  clout_Vloc =  './OUTBOX/Vloc_Rloc.OUT.    '
  
  ! Setting the index of the output files
  compteur_clout = compteur_clout+1
  
  if (compteur_clout<10) then
    WRITE(clout_DOF(18:19),'(I1)')   compteur_clout
    WRITE(clout_Vloc(24:25),'(I1)')  compteur_clout
  else if ( (compteur_clout>=10) .and. (compteur_clout<100) ) then
    WRITE(clout_DOF(18:20),'(I2)')   compteur_clout
    WRITE(clout_Vloc(24:26),'(I2)')  compteur_clout
  else if ( (compteur_clout>=100).and. (compteur_clout<1000) ) then
    WRITE(clout_DOF(18:21),'(I3)')   compteur_clout
    WRITE(clout_Vloc(24:27),'(I3)')  compteur_clout
  end if
  
  ! Reading the BODIES.OUT file
  if (first_over_all) then
    open(unit=2,file=clout_Bodies,status='old')
    print*,'--->',clout_Bodies
    
    allocate(TAB_SPHERE(n_particles))
    allocate(TAB_PLAN(n_walls))
    i=0
    do
      read(2,'(A6)') text
      if (text == '$bdyty') then
        i = i +1
        ! Storing the particles information
        if (i<n_particles+1) then
          TAB_SPHERE(i)%id=i
          read(2,*)
          read(2,*)
          read(2,'(22X,A5)') behav
          TAB_SPHERE(i)%behav = behav
          read(2,'(29X, 3(5x,D14.7,2X))') I1,I2,I3
          read(2,*)
          read(2,'(29X, 3(5x,D14.7,2X))') center1,center2,center3
          TAB_SPHERE(i)%I1=I1
          TAB_SPHERE(i)%I2=I2
          TAB_SPHERE(i)%I3=I3
          TAB_SPHERE(i)%center_ref(1)=center1
          TAB_SPHERE(i)%center_ref(2)=center2
          TAB_SPHERE(i)%center_ref(3)=center3
          read(2,*)
          read(2,*)
          read(2,'(22X,A5,7X,D14.7)') color,rayon
          TAB_SPHERE(i)%color = color
          TAB_SPHERE(i)%Rmax = rayon
          TAB_SPHERE(i)%volume = 4.188790204*rayon**3
          TAB_SPHERE(i)%center =0.D0
          TAB_SPHERE(i)%inbox=.false.
          TAB_SPHERE(i)%momentI=0
          TAB_SPHERE(i)%nctc = 0
          TAB_SPHERE(i)%nctcslide = 0
          TAB_SPHERE(i)%nctcstick = 0
        end if
        ! Storing the walls information
        if (i>n_particles) then
          TAB_PLAN(i-n_particles)%id=i
          read(2,*)
          read(2,*)
          read(2,*)
          read(2,*)
          read(2,*)
          read(2,'(29X, 3(5x,D14.7,2X))') center1,center2,center3
          TAB_PLAN(i-n_particles)%center_ref(1)=center1
          TAB_PLAN(i-n_particles)%center_ref(2)=center2
          TAB_PLAN(i-n_particles)%center_ref(3)=center3
          read(2,*)
          read(2,*)
          !read(2,'(29X, 3(5x,D14.7,2X))') ax1r,ax2r,ax3r
          if (first_over_all) then 
            TAB_PLAN(i-n_particles)%ax1=0.01
            TAB_PLAN(i-n_particles)%ax2=0.01
            TAB_PLAN(i-n_particles)%ax3=0.01
          end if 
          !read(2,*)
          TAB_PLAN(i-n_particles)%center(:)=0.D0
          TAB_PLAN(i-n_particles)%behav = 'WALL_'
          
        end if
        if (i==n_particles+n_walls) exit
      end if
      first_bodies = .false.
    end do
  end if
  
  open(unit=5,file=clout_DOF,iostat=err,status='old')
  
  if (err/=0) then
    fin_post=.true.
  else
    fin_post=.false.
    i = 0
    read(5,*)
    read(5,*)
    read(5,*)
    read(5,'(7X,i9,19X,D14.7)') step,time
    
    print*,'--->',clout_DOF
    read(5,*)
    read(5,*)
    
    do
      read(5,'(A6)') text
      if (text == '!-----') exit
      if (text == '$bdyty') then
        i = i+1
        read(5,'(6X,i9)') num_part
        if (num_part<n_particles+1) then
          read(5,*)
          read(5,'(29X, 3(5x,D14.7,2X))') coor1,coor2,coor3
          TAB_SPHERE(i)%center(1)=coor1+TAB_SPHERE(i)%center_ref(1)
          TAB_SPHERE(i)%center(2)=coor2+TAB_SPHERE(i)%center_ref(2)
          TAB_SPHERE(i)%center(3)=coor3+TAB_SPHERE(i)%center_ref(3)
          read(5,*)
          read(5,'(29X, 3(5x,D14.7,2X))') Vx,Vy,Vz
          read(5,'(29X, 3(5x,D14.7,2X))') Vrx,Vry,Vrz
          TAB_SPHERE(i)%Vx=Vx
          TAB_SPHERE(i)%Vy=Vy
          TAB_SPHERE(i)%Vz=Vz
          TAB_SPHERE(i)%Vrx=Vrx
          TAB_SPHERE(i)%Vry=Vry
          TAB_SPHERE(i)%Vrz=Vrz
          TAB_SPHERE(i)%Norme_V  = sqrt(Vx**2+Vy**2+Vz**2)
          TAB_SPHERE(i)%Norme_Vr = sqrt(Vrx**2+Vry**2+Vrz**2)
          read(5,'(29X, 3(5x,D14.7,2X))') rot(1,1),rot(2,1),rot(3,1)
          read(5,'(29X, 3(5x,D14.7,2X))') rot(1,2),rot(2,2),rot(3,2)
          read(5,'(29X, 3(5x,D14.7,2X))') rot(1,3),rot(2,3),rot(3,3)
          TAB_SPHERE(i)%Rot = rot
          TAB_SPHERE(i)%Wx = Vrx*rot(1,1) + Vry*rot(1,2) + Vrz*rot(1,3)
          TAB_SPHERE(i)%Wy = Vrx*rot(2,1) + Vry*rot(2,2) + Vrz*rot(2,3)
          TAB_SPHERE(i)%Wz = Vrx*rot(3,1) + Vry*rot(3,2) + Vrz*rot(3,3)
          TAB_SPHERE(i)%nctc = 0
        else 
          read(5,*)
          read(5,'(29X, 3(5x,D14.7,2X))') coor1,coor2,coor3
          TAB_PLAN(i-n_particles)%center(1)=coor1+TAB_PLAN(i-n_particles)%center_ref(1)
          TAB_PLAN(i-n_particles)%center(2)=coor2+TAB_PLAN(i-n_particles)%center_ref(2)
          TAB_PLAN(i-n_particles)%center(3)=coor3+TAB_PLAN(i-n_particles)%center_ref(3)
          read(5,*)
          read(5,'(29X, 3(5x,D14.7,2X))') Vx,Vy,Vz
          read(5,'(29X, 3(5x,D14.7,2X))') Vrx,Vry,Vrz
          TAB_PLAN(i-n_particles)%Vx=Vx
          TAB_PLAN(i-n_particles)%Vy=Vy
          TAB_PLAN(i-n_particles)%Vz=Vz
          TAB_PLAN(i-n_particles)%Vrx=Vrx
          TAB_PLAN(i-n_particles)%Vry=Vry
          TAB_PLAN(i-n_particles)%Vrz=Vrz
          TAB_PLAN(i-n_particles)%Norme_V  = sqrt(Vx**2+Vy**2+Vz**2)
          TAB_PLAN(i-n_particles)%Norme_Vr = sqrt(Vrx**2+Vry**2+Vrz**2)
          read(5,'(29X, 3(5x,D14.7,2X))') rot(1,1),rot(2,1),rot(3,1)
          read(5,'(29X, 3(5x,D14.7,2X))') rot(1,2),rot(2,2),rot(3,2)
          read(5,'(29X, 3(5x,D14.7,2X))') rot(1,3),rot(2,3),rot(3,3)
          TAB_PLAN(i-n_particles)%Rot = rot
          TAB_PLAN(i-n_particles)%Wx = Vrx*rot(1,1) + Vry*rot(1,2) + Vrz*rot(1,3)
          TAB_PLAN(i-n_particles)%Wy = Vrx*rot(2,1) + Vry*rot(2,2) + Vrz*rot(2,3)
          TAB_PLAN(i-n_particles)%Wz = Vrx*rot(3,1) + Vry*rot(3,2) + Vrz*rot(3,3)
          TAB_PLAN(i-n_particles)%nctc = 0


                    ! Similar to poly, i'm going to create vertices (8 vertices in each wall)
          allocate(TAB_PLAN(i-n_particles)%vertex_ref(3,8))
          allocate(TAB_PLAN(i-n_particles)%vertex(3,8))
          
          ! First sheet
          TAB_PLAN(i-n_particles)%vertex_ref(1,1) = TAB_PLAN(i-n_particles)%ax1
          TAB_PLAN(i-n_particles)%vertex_ref(2,1) = TAB_PLAN(i-n_particles)%ax2
          TAB_PLAN(i-n_particles)%vertex_ref(3,1) = TAB_PLAN(i-n_particles)%ax3
          
          TAB_PLAN(i-n_particles)%vertex_ref(1,2) = (-TAB_PLAN(i-n_particles)%ax1)
          TAB_PLAN(i-n_particles)%vertex_ref(2,2) = TAB_PLAN(i-n_particles)%ax2
          TAB_PLAN(i-n_particles)%vertex_ref(3,2) = TAB_PLAN(i-n_particles)%ax3
          
          TAB_PLAN(i-n_particles)%vertex_ref(1,3) = (-TAB_PLAN(i-n_particles)%ax1)
          TAB_PLAN(i-n_particles)%vertex_ref(2,3) = (-TAB_PLAN(i-n_particles)%ax2)
          TAB_PLAN(i-n_particles)%vertex_ref(3,3) = TAB_PLAN(i-n_particles)%ax3
          
          TAB_PLAN(i-n_particles)%vertex_ref(1,4) = TAB_PLAN(i-n_particles)%ax1
          TAB_PLAN(i-n_particles)%vertex_ref(2,4) = (-TAB_PLAN(i-n_particles)%ax2)
          TAB_PLAN(i-n_particles)%vertex_ref(3,4) = TAB_PLAN(i-n_particles)%ax3
          
          ! Second sheet
          TAB_PLAN(i-n_particles)%vertex_ref(1,5) =  TAB_PLAN(i-n_particles)%ax1
          TAB_PLAN(i-n_particles)%vertex_ref(2,5) =  TAB_PLAN(i-n_particles)%ax2
          TAB_PLAN(i-n_particles)%vertex_ref(3,5) = (- TAB_PLAN(i-n_particles)%ax3)
          
          TAB_PLAN(i-n_particles)%vertex_ref(1,6) = (- TAB_PLAN(i-n_particles)%ax1)
          TAB_PLAN(i-n_particles)%vertex_ref(2,6) =  TAB_PLAN(i-n_particles)%ax2
          TAB_PLAN(i-n_particles)%vertex_ref(3,6) = (- TAB_PLAN(i-n_particles)%ax3)
          
          TAB_PLAN(i-n_particles)%vertex_ref(1,7) = (- TAB_PLAN(i-n_particles)%ax1)
          TAB_PLAN(i-n_particles)%vertex_ref(2,7) = (- TAB_PLAN(i-n_particles)%ax2)
          TAB_PLAN(i-n_particles)%vertex_ref(3,7) = (- TAB_PLAN(i-n_particles)%ax3)
          
          TAB_PLAN(i-n_particles)%vertex_ref(1,8) =  TAB_PLAN(i-n_particles)%ax1
          TAB_PLAN(i-n_particles)%vertex_ref(2,8) = (- TAB_PLAN(i-n_particles)%ax2)
          TAB_PLAN(i-n_particles)%vertex_ref(3,8) = (- TAB_PLAN(i-n_particles)%ax3)

          do j=1, 8
            TAB_PLAN(i-n_particles)%vertex(1,j) = TAB_PLAN(i-n_particles)%vertex_ref(1,j) *rot(1,1)+&
                                                  TAB_PLAN(i-n_particles)%vertex_ref(2,j) *rot(1,2)+&
                                                  TAB_PLAN(i-n_particles)%vertex_ref(3,j) *rot(1,3)!+&
                                                  !TAB_PLAN(i-n_particles)%center(1)
            
            TAB_PLAN(i-n_particles)%vertex(2,j) = TAB_PLAN(i-n_particles)%vertex_ref(1,j) *rot(2,1)+&
                                                  TAB_PLAN(i-n_particles)%vertex_ref(2,j) *rot(2,2)+&
                                                  TAB_PLAN(i-n_particles)%vertex_ref(3,j) *rot(2,3)!+&
                                                  !TAB_PLAN(i-n_particles)%center(2)
            
            TAB_PLAN(i-n_particles)%vertex(3,j) = TAB_PLAN(i-n_particles)%vertex_ref(1,j) *rot(3,1)+&
                                                  TAB_PLAN(i-n_particles)%vertex_ref(2,j) *rot(3,2)+&
                                                  TAB_PLAN(i-n_particles)%vertex_ref(3,j) *rot(3,3)!+&
                                                  !TAB_PLAN(i-n_particles)%center(3)
          end do

          ! Rotating wall 45 degrees (around z) because weird initial configuration was written by lmgc90
          do j=1, 8
            coor1 = TAB_PLAN(i-n_particles)%vertex(1,j)
            coor2 = TAB_PLAN(i-n_particles)%vertex(2,j)
            coor3 = TAB_PLAN(i-n_particles)%vertex(3,j)

            TAB_PLAN(i-n_particles)%vertex(1,j) = coor1*cos(PI/4) - coor2*sin(PI/4) + TAB_PLAN(i-n_particles)%center(1)
            TAB_PLAN(i-n_particles)%vertex(2,j) = coor1*sin(PI/4) + coor2*cos(PI/4) + TAB_PLAN(i-n_particles)%center(2)
            TAB_PLAN(i-n_particles)%vertex(3,j) = coor3 + TAB_PLAN(i-n_particles)%center(3)
          end do 
        end if
        if (i==n_particles+n_walls) exit
      end if
    end do
  end if
  
  ! Number of contacts between spheres or sphere-wall
  nb_SPSPx=0
  nb_SPPLx=0
  
  open(unit=4,file=clout_Vloc,iostat=err,status='old')
  if (err/=0) then
    ! Void
  else
    print*,'--->',clout_Vloc
    read(4,*)
    read(4,*)
    read(4,*)
    read(4,*)
    read(4,*)
    read(4,*)
    
    do
      read(4,'(A13)',iostat=err) text2
      if (err/=0) then
        close(2)
        close(4)
        close(5)
        exit
      end if
      if (text2 == '!-------------') exit
      if (text2 == '$icdan  SPPLx') then
        nb_SPPLx = nb_SPPLx +1
      end if
      if (text2 == '$icdan  SPSPx') then
        nb_SPSPx = nb_SPSPx +1
      end if
    end do
    
    close(4)
    
    if (allocated(TAB_CONTACT)) deallocate(TAB_CONTACT)
    allocate(TAB_CONTACT(nb_SPSPx+nb_SPPLx))
    
    do i=1,nb_SPSPx+nb_SPPLx
      TAB_CONTACT(i)%icdent=0
      TAB_CONTACT(i)%ianent=0
      TAB_CONTACT(i)%n(1)=0.D0
      TAB_CONTACT(i)%n(2)=0.D0
      TAB_CONTACT(i)%n(3)=0.D0
      TAB_CONTACT(i)%t(1)=0.D0
      TAB_CONTACT(i)%t(2)=0.D0
      TAB_CONTACT(i)%t(3)=0.D0
      TAB_CONTACT(i)%s(1)=0.D0
      TAB_CONTACT(i)%s(2)=0.D0
      TAB_CONTACT(i)%s(3)=0.D0
      TAB_CONTACT(i)%coor_ctc(1)=0.D0
      TAB_CONTACT(i)%coor_ctc(2)=0.D0
      TAB_CONTACT(i)%coor_ctc(3)=0.D0
      TAB_CONTACT(i)%rn=0.D0
      TAB_CONTACT(i)%rt=0.D0
      TAB_CONTACT(i)%rs=0.D0
      TAB_CONTACT(i)%vls=0.D0
      TAB_CONTACT(i)%vln=0.D0
      TAB_CONTACT(i)%vlt=0.D0
      TAB_CONTACT(i)%type=0
      TAB_CONTACT(i)%nature = '    '
    end do

    n_spsp = 0
    n_sppl = 0
    nb_ligneCONTACT=0
    open(unit=4,file=clout_Vloc,status='old')
    read(4,*) 
    read(4,*) 
    read(4,*) 
    read(4,*) 
    read(4,*) 
    read(4,*) 

    do
      read(4,'(A13)',iostat=err) text2
      if (err/=0) then
        close(4)
        exit
      end if
      if (text2 == '!-------------') exit
      if (text2 == '$icdan  SPPLx') then
        read(4,*)
        read(4,'(7X,i6,36X,i6,16X,A5)') icdent,ianent,statut
        read(4,'(29X, 3(5x,D14.7,2X))') rs,rt,rn
        read(4,'(29X, 3(5x,D14.7,2X))') vls,vlt,vln
        read(4,*)
        read(4,'(29X, 3(5x,D14.7,2X))') s1,s2,s3
        read(4,'(29X, 3(5x,D14.7,2X))') t1,t2,t3
        read(4,'(29X, 3(5x,D14.7,2X))') n1,n2,n3
        read(4,'(29X, 3(5x,D14.7,2X))') coor_ctc1,coor_ctc2,coor_ctc3
        
        nb_ligneCONTACT=nb_ligneCONTACT+1
        n_sppl = n_sppl + 1
        TAB_CONTACT(n_sppl+nb_SPSPx)%icdent=icdent
        TAB_CONTACT(n_sppl+nb_SPSPx)%ianent=ianent
        TAB_CONTACT(n_sppl+nb_SPSPx)%statut=statut
        TAB_CONTACT(n_sppl+nb_SPSPx)%n(1)=n1
        TAB_CONTACT(n_sppl+nb_SPSPx)%n(2)=n2
        TAB_CONTACT(n_sppl+nb_SPSPx)%n(3)=n3
        TAB_CONTACT(n_sppl+nb_SPSPx)%t(1)=t1
        TAB_CONTACT(n_sppl+nb_SPSPx)%t(2)=t2
        TAB_CONTACT(n_sppl+nb_SPSPx)%t(3)=t3
        TAB_CONTACT(n_sppl+nb_SPSPx)%s(1)=s1
        TAB_CONTACT(n_sppl+nb_SPSPx)%s(2)=s2
        TAB_CONTACT(n_sppl+nb_SPSPx)%s(3)=s3
        TAB_CONTACT(n_sppl+nb_SPSPx)%coor_ctc(1)=coor_ctc1
        TAB_CONTACT(n_sppl+nb_SPSPx)%coor_ctc(2)=coor_ctc2
        TAB_CONTACT(n_sppl+nb_SPSPx)%coor_ctc(3)=coor_ctc3
        TAB_CONTACT(n_sppl+nb_SPSPx)%rn=rn
        TAB_CONTACT(n_sppl+nb_SPSPx)%rt=rt
        TAB_CONTACT(n_sppl+nb_SPSPx)%rs=rs
        TAB_CONTACT(n_sppl+nb_SPSPx)%type=1
        TAB_CONTACT(n_sppl+nb_SPSPx)%nature = 'SPPLx'
      end if
      
      if (text2 == '$icdan  SPSPx') then
        read(4,*)
        read(4,'(7X,i6,36X,i6,16X,A5)') icdent,ianent,statut
        read(4,'(29X, 3(5x,D14.7,2X))') rs,rt,rn
        read(4,'(29X, 3(5x,D14.7,2X))') vls,vlt,vln
        read(4,'(29X,5x,D14.7)')        gap
        read(4,'(29X, 3(5x,D14.7,2X))') s1,s2,s3
        read(4,'(29X, 3(5x,D14.7,2X))') t1,t2,t3
        read(4,'(29X, 3(5x,D14.7,2X))') n1,n2,n3
        read(4,'(29X, 3(5x,D14.7,2X))') coor_ctc1,coor_ctc2,coor_ctc3
        
        nb_ligneCONTACT=nb_ligneCONTACT+1
        n_spsp = n_spsp + 1 
        TAB_CONTACT(n_spsp)%icdent=icdent
        TAB_CONTACT(n_spsp)%ianent=ianent
        TAB_CONTACT(n_spsp)%statut=statut
        TAB_CONTACT(n_spsp)%n(1)=n1
        TAB_CONTACT(n_spsp)%n(2)=n2
        TAB_CONTACT(n_spsp)%n(3)=n3
        TAB_CONTACT(n_spsp)%t(1)=t1
        TAB_CONTACT(n_spsp)%t(2)=t2
        TAB_CONTACT(n_spsp)%t(3)=t3
        TAB_CONTACT(n_spsp)%s(1)=s1
        TAB_CONTACT(n_spsp)%s(2)=s2
        TAB_CONTACT(n_spsp)%s(3)=s3
        TAB_CONTACT(n_spsp)%coor_ctc(1)=coor_ctc1
        TAB_CONTACT(n_spsp)%coor_ctc(2)=coor_ctc2
        TAB_CONTACT(n_spsp)%coor_ctc(3)=coor_ctc3
        TAB_CONTACT(n_spsp)%rn=rn
        TAB_CONTACT(n_spsp)%rt=rt
        TAB_CONTACT(n_spsp)%rs=rs
        TAB_CONTACT(n_spsp)%vln=vln
        TAB_CONTACT(n_spsp)%vlt=vlt
        TAB_CONTACT(n_spsp)%vls=vls
        TAB_CONTACT(n_spsp)%gap=gap
        TAB_CONTACT(n_spsp)%type=1
        if (ianent .gt. n_particles)then 
          TAB_CONTACT(n_spsp)%nature = 'SPPLx'
        else 
          TAB_CONTACT(n_spsp)%nature = 'SPSPx'
        end if 
      end if
    end do
  end if
  
  Mean_total_Normal_force = 0.D0
  Number_mean_Total_force = 0.D0
  
  do i=1,nb_ligneCONTACT
    if (TAB_CONTACT(i)%rn == 0.D0 ) cycle
    Mean_total_Normal_force = Mean_total_Normal_force + TAB_CONTACT(i)%rn
    Number_mean_Total_force = Number_mean_Total_force + 1.D0
  end do
  
  Mean_total_Normal_force = Mean_total_Normal_force / Number_mean_Total_force

  ! Computing the number of active contacts per particle
  ! Adding the number of contacts per particle
  TAB_SPHERE(:)%nctc = 0.
  TAB_SPHERE(:)%av_mob = 0.
  TAB_SPHERE(:)%av_fn = 0.

  do i=1,nb_ligneCONTACT
    icdent = TAB_CONTACT(i)%icdent
    ianent = TAB_CONTACT(i)%ianent
    ! Contacts between particles
    if (TAB_CONTACT(i)%nature == 'SPPLx') cycle
    if ((icdent .gt. n_particles) .or. (ianent .gt. n_particles)) cycle
    ! Only active contacts
    if (TAB_CONTACT(i)%rn .le. 0.) cycle
    
    TAB_SPHERE(icdent)%nctc = TAB_SPHERE(icdent)%nctc + 1
    TAB_SPHERE(ianent)%nctc = TAB_SPHERE(ianent)%nctc + 1

    if (c_fn_class == 1 .or. c_pdf_z==1) then 
      TAB_SPHERE(icdent)%av_fn = TAB_SPHERE(icdent)%av_fn + TAB_CONTACT(i)%rn
      TAB_SPHERE(ianent)%av_fn = TAB_SPHERE(ianent)%av_fn + TAB_CONTACT(i)%rn
    end if 

    if (c_mobilization_class == 1) then 
      TAB_SPHERE(icdent)%av_mob = TAB_SPHERE(icdent)%av_mob + &
                                  sqrt(TAB_CONTACT(i)%rt**2 + TAB_CONTACT(i)%rs**2)/(0.4*TAB_CONTACT(i)%rn)
      TAB_SPHERE(ianent)%av_mob = TAB_SPHERE(ianent)%av_mob + &
                                  sqrt(TAB_CONTACT(i)%rt**2 + TAB_CONTACT(i)%rs**2)/(0.4*TAB_CONTACT(i)%rn)
    end if 
  end do

  ! Finding averages
  if (c_fn_class == 1 .or. c_pdf_z==1) then 
    do i=1, n_particles
      if (TAB_SPHERE(i)%nctc .gt. 0) then 
        TAB_SPHERE(i)%av_fn = TAB_SPHERE(i)%av_fn/TAB_SPHERE(i)%nctc
      end if 
    end do 
  end if 

  if (c_mobilization_class == 1) then 
    do i=1, n_particles
      if (TAB_SPHERE(i)%nctc .gt. 0) then 
        TAB_SPHERE(i)%av_mob = TAB_SPHERE(i)%av_mob/TAB_SPHERE(i)%nctc
      end if 
    end do 
  end if
  
  print*, 'Nb particules :', n_particles
  print*, 'Nb Contacts :', nb_ligneCONTACT
  
  print*, 'Step:', step
  print*, 'Time:', time
  
end subroutine read_Vloc_dof_bodies


!==============================================================================
!Calcul des deformation
!==============================================================================
subroutine calcul_deformation
  
  implicit none
  
  integer               :: i
  real(kind=8)          :: Xplan_max,Xplan_min,Yplan_max,Yplan_min,Zplan_max,Zplan_min
  real(kind=8)          :: E1,E2,E3,DE1,DE2,DE3
  
  if (deformation_ini) then
    Xplan_max =-10000000
    Xplan_min = 10000000
    Yplan_max =-10000000
    Yplan_min = 10000000
    Zplan_max =-10000000
    Zplan_min = 10000000
    
    ! Finding the size of the box
    do i=1,n_particles
      Xplan_max = max(Xplan_max, TAB_SPHERE(i)%center(1) + TAB_SPHERE(i)%Rmax)
      Xplan_min = min(Xplan_min, TAB_SPHERE(i)%center(1) - TAB_SPHERE(i)%Rmax)
      
      Yplan_max = max(Yplan_max, TAB_SPHERE(i)%center(2) + TAB_SPHERE(i)%Rmax)
      Yplan_min = min(Yplan_min, TAB_SPHERE(i)%center(2) - TAB_SPHERE(i)%Rmax)
      
      Zplan_max = max(Zplan_max, TAB_SPHERE(i)%center(3) + TAB_SPHERE(i)%Rmax)
      Zplan_min = min(Zplan_min, TAB_SPHERE(i)%center(3) - TAB_SPHERE(i)%Rmax)
    end do
    
    Long_ini = Xplan_max-Xplan_min
    Larg_ini = Yplan_max-Yplan_min
    H_ini    = Zplan_max-Zplan_min
    
    Depsilon_p = 0.D0
    Depsilon_q = 0.D0
    
    deformation_ini=.false.
    
  end if
  
  Xplan_max =-10000000
  Xplan_min = 10000000
  Yplan_max =-10000000
  Yplan_min = 10000000
  Zplan_max =-10000000
  Zplan_min = 10000000
  
  do i=1,n_particles
    
    Xplan_max = max(Xplan_max, TAB_SPHERE(i)%center(1) + TAB_SPHERE(i)%Rmax)
    Xplan_min = min(Xplan_min, TAB_SPHERE(i)%center(1) - TAB_SPHERE(i)%Rmax)
    
    Yplan_max = max(Yplan_max, TAB_SPHERE(i)%center(2) + TAB_SPHERE(i)%Rmax)
    Yplan_min = min(Yplan_min, TAB_SPHERE(i)%center(2) - TAB_SPHERE(i)%Rmax)
    
    Zplan_max = max(Zplan_max, TAB_SPHERE(i)%center(3) + TAB_SPHERE(i)%Rmax)
    Zplan_min = min(Zplan_min, TAB_SPHERE(i)%center(3) - TAB_SPHERE(i)%Rmax)
    
  end do
  
  large  = Xplan_max - Xplan_min
  width  = Yplan_max - Yplan_min
  height = Zplan_max - Zplan_min

  ! Adding the size of the walls 
  if (first_over_all) then 
    TAB_PLAN(1)%ax1=large/2.0
    TAB_PLAN(1)%ax2=width/2.0
    TAB_PLAN(1)%ax3=0.01

    TAB_PLAN(2)%ax1=large/2.0
    TAB_PLAN(2)%ax2=width/2.0
    TAB_PLAN(2)%ax3=0.01
  end if 
  
  !E1 = log( Long / Long_ini )
  !E2 = log( Larg / Larg_ini )
  !E3 = log( Haut / H_ini    )
  
  !Depsilon_p = ( DE1 + DE2 + DE3 )
  !Depsilon_q = ( DE3 - DE1 ) 
  
  !Ep = E1+E2+E3
  !Eq = E3-E1
  
  !Def   = 0.047472569*time
  !Def_H = abs(Haut - H_ini )
  
  !print*, 'Deformation = ',Def
  !print*, 'Long =' , Long, 'Larg =' , Larg , 'haut =' , Haut
  
end subroutine calcul_deformation


!==============================================================================
! Walls positions
!==============================================================================
subroutine walls_position
  
  implicit none
  
  if (first_over_all) then
    write(108,*) '#   time     ', '   wall 1-X  ', '   wall 1-Y  ', '   wall 1-Z  ', &
                                  '   wall 2-X  ', '   wall 2-Y  ', '   wall 2-Z  '
  endif
  
  write(108,'(7(1X,E12.5))') time, TAB_PLAN(1)%center(1), TAB_PLAN(1)%center(2), TAB_PLAN(1)%center(3), &
                                   TAB_PLAN(2)%center(1), TAB_PLAN(2)%center(2), TAB_PLAN(2)%center(3)
  
  print*, 'Write Walls position            ---> Ok!'
  
end subroutine walls_position


!==============================================================================
! Walls forces
!==============================================================================
subroutine walls_force
  
  implicit none
  
  integer                              :: i, idw1, idw2
  real*8                               :: fx1, fx2, fx3
  real*8                               :: fy1, fy2, fy3
  real*8                               :: fz1, fz2, fz3
  real*8                               :: rn, rt, rs
  real*8, dimension(3)                 :: nik, tik, sik
  
  ! Initializing variables (n->z, t->x, s->y)
  fx1 = 0.0
  fx2 = 0.0
  
  fy1 = 0.0
  fy2 = 0.0

  fz1 = 0.0
  fz2 = 0.0
  
  ! Ids for the walls
  idw1 = n_particles + 1 
  idw2 = n_particles + 2 
  
  ! Computing 
  do i=1, nb_ligneCONTACT
    if (TAB_CONTACT(i)%nature== 'SPPLx') then

      rn = TAB_CONTACT(i)%rn
      rt = TAB_CONTACT(i)%rt
      rs = TAB_CONTACT(i)%rs
      nik = TAB_CONTACT(i)%n
      tik = TAB_CONTACT(i)%t
      sik = TAB_CONTACT(i)%s

      if(TAB_CONTACT(i)%ianent == idw1) then
        fx1 = fx1 + (rn*nik(1) + rt*tik(1) + rs*sik(1))
        fy1 = fy1 + (rn*nik(2) + rt*tik(2) + rs*sik(2))
        fz1 = fz1 + (rn*nik(3) + rt*tik(3) + rs*sik(3))
      elseif(TAB_CONTACT(i)%ianent == idw2) then
        fx2 = fx2 + (rn*nik(1) + rt*tik(1) + rs*sik(1))
        fy2 = fy2 + (rn*nik(2) + rt*tik(2) + rs*sik(2))
        fz2 = fz2 + (rn*nik(3) + rt*tik(3) + rs*sik(3))
      end if
    end if
  end do

  if (first_over_all) then
    write(109,*) '#   time     ', '     fx1     ', '     fy1     ', '     fz1     ', &
                                  '     fx2     ', '     fy2     ', '     fz2     '
  endif
  
  write(109,'(7(1X,E12.5))') time, fx1, fy1, fz1, fx2, fy2, fz2 
  
  print*, 'Write Walls forces  ---> Ok!'
  
end subroutine walls_force


!==============================================================================
! Computing the compacity
!==============================================================================
subroutine compacity
  
  implicit none
  
  integer                                  :: i, j, k, l, mesh_ele, n_ele_in 
  integer, dimension(:), allocatable       :: status_part
  real*8                                   :: vol_grains, vol_box, l_height, l_width, delta_mesh
  real*8                                   :: d_towall, ave_rad, lx_max, ly_max, lx_min, ly_min
  real*8                                   :: norm_ele
  real*8, dimension(3)                     :: center_ele_ini, center_ele 
  
  real*8                                   :: pos_z_w1, pos_z_w2


  ! Initializing
  vol_grains = 0.D0
  mesh_ele = 75

  pos_z_w1 = TAB_PLAN(1)%center(3)
  pos_z_w2 = TAB_PLAN(2)%center(3)

  if (allocated(status_part)) deallocate(status_part)
  allocate(status_part(n_particles))

  status_part(:) = 0

  lx_min =  9999.0
  lx_max = -9999.0
  ly_min =  9999.0
  ly_max = -9999.0

  ! Computing the average radius of the sample 
  ave_rad = 0.0 
  do i=1, n_particles
    ave_rad = ave_rad + TAB_SPHERE(i)%Rmax
  end do

  ave_rad = ave_rad/n_particles

  ! Separation to the wall 
  d_towall = ave_rad*4

  !print*, "dwall"
  !print*, d_towall
  !print*, "----"
  ! The height of the sub-box 
  l_height = abs(pos_z_w1 - pos_z_w2) - 2*d_towall

  ! Giving status to the particles (1: Inside, 2: Partially inside, 3: Outside).
  do i=1, n_particles
    if (status_part(i) == 3) cycle

    ! Testing the lower border 
    if (TAB_SPHERE(i)%center(3)+TAB_SPHERE(i)%Rmax .gt. (pos_z_w1 + d_towall)) then 
      if (TAB_SPHERE(i)%center(3)-TAB_SPHERE(i)%Rmax .gt. (pos_z_w1 + d_towall)) then 
        status_part(i) = 1
      else 
        status_part(i) = 21
      end if 
    else 
      status_part(i) = 3
    end if

    if (status_part(i) == 3 .or. status_part(i) == 21) cycle

    ! Testing the upper border 
    if (TAB_SPHERE(i)%center(3)-TAB_SPHERE(i)%Rmax .lt. (pos_z_w2 - d_towall)) then 
      if (TAB_SPHERE(i)%center(3)+TAB_SPHERE(i)%Rmax .lt. (pos_z_w2 - d_towall)) then 
        status_part(i) = 1
      else 
        status_part(i) = 22
      end if 
    else 
      status_part(i) = 3
    end if 
  end do

  ! Adding the volume of the particles as function of their status.
  ! Also computing the sample width and large 
  do i=1, n_particles
    !print*, i 
    !print*, status_part(i)
    !print*, vol_grains
    ! The particles outside 
    if (status_part(i) == 3) cycle 

    ! Checking the width and large 
    if (lx_min .gt. TAB_SPHERE(i)%center(1)-TAB_SPHERE(i)%Rmax) lx_min = TAB_SPHERE(i)%center(1)-TAB_SPHERE(i)%Rmax
    if (lx_max .lt. TAB_SPHERE(i)%center(1)+TAB_SPHERE(i)%Rmax) lx_max = TAB_SPHERE(i)%center(1)+TAB_SPHERE(i)%Rmax

    if (ly_min .gt. TAB_SPHERE(i)%center(2)-TAB_SPHERE(i)%Rmax) ly_min = TAB_SPHERE(i)%center(2)-TAB_SPHERE(i)%Rmax
    if (ly_max .lt. TAB_SPHERE(i)%center(2)+TAB_SPHERE(i)%Rmax) ly_max = TAB_SPHERE(i)%center(2)+TAB_SPHERE(i)%Rmax

    ! The particles completely inside 
    if (status_part(i) == 1) then 
      vol_grains = vol_grains + TAB_SPHERE(i)%volume

    ! The particles partially inside towards the lower boundary 
    else if(status_part(i) == 21) then
      ! Size mesh 
      delta_mesh = 2.0*TAB_SPHERE(i)%Rmax/mesh_ele

      n_ele_in = 0
      center_ele_ini(:) = TAB_SPHERE(i)%center(:) - TAB_SPHERE(i)%Rmax
      do j=1, mesh_ele
        do k=1, mesh_ele
          do l=1, mesh_ele
            center_ele(1) = center_ele_ini(1) + delta_mesh*(j-0.5)
            center_ele(2) = center_ele_ini(2) + delta_mesh*(k-0.5)
            center_ele(3) = center_ele_ini(3) + delta_mesh*(l-0.5)

            ! The norm from the sphere center to the element center 
            norm_ele = (TAB_SPHERE(i)%center(1)-center_ele(1))**2 + &
                       (TAB_SPHERE(i)%center(2)-center_ele(2))**2 + &
                       (TAB_SPHERE(i)%center(3)-center_ele(3))**2

            norm_ele = sqrt(norm_ele)

            ! The case it is under the box 
            if (center_ele(3) .lt. (pos_z_w1 + d_towall)) cycle
            ! The case it is not inside the sphere volume 
            if (norm_ele .gt. TAB_SPHERE(i)%Rmax) cycle
            ! Otherwise it is an element of the sphere and we count it   
            n_ele_in = n_ele_in + 1
          end do 
        end do 
      end do

      vol_grains = vol_grains + (n_ele_in*(delta_mesh**3))

    ! The particles partially inside upper the lower boundary 
    else if(status_part(i) == 22) then 
      ! Size mesh 
      delta_mesh = 2.0*TAB_SPHERE(i)%Rmax/mesh_ele

      n_ele_in = 0
      center_ele_ini(:) = TAB_SPHERE(i)%center(:) - TAB_SPHERE(i)%Rmax
      do j=1, mesh_ele
        do k=1, mesh_ele
          do l=1, mesh_ele
            center_ele(1) = center_ele_ini(1) + delta_mesh*(j-0.5)
            center_ele(2) = center_ele_ini(2) + delta_mesh*(k-0.5)
            center_ele(3) = center_ele_ini(3) + delta_mesh*(l-0.5)

            ! The norm from the sphere center to the element center 
            norm_ele = (TAB_SPHERE(i)%center(1)-center_ele(1))**2 + &
                       (TAB_SPHERE(i)%center(2)-center_ele(2))**2 + &
                       (TAB_SPHERE(i)%center(3)-center_ele(3))**2

            norm_ele = sqrt(norm_ele)

            ! The case it is under the box 
            if (center_ele(3) .gt. (pos_z_w2 - d_towall)) cycle
            ! The case it is not inside the sphere volume 
            if (norm_ele .gt. TAB_SPHERE(i)%Rmax) cycle
            ! Otherwise it is an element of the sphere and we count it   
            n_ele_in = n_ele_in + 1
          end do 
        end do 
      end do

      vol_grains = vol_grains + (n_ele_in*(delta_mesh**3))
    end if
  end do 

  vol_box = l_height * (lx_max - lx_min) * (ly_max - ly_min)
  
  if (first_over_all) then
    write(102,*) '#   time     ', '   height    ', '   width    ', '    large    ',' V_grain/V_box '
  end if
  
  write(102,'(5(1X,E12.5))') time, l_height, (lx_max - lx_min), (ly_max - ly_min), vol_grains/vol_box
  
  write (*,*) 'Write Compacity     ---> Ok!'
  
end subroutine compacity

!==============================================================================
! Computing the coordination
!==============================================================================
subroutine nb_coordination

  implicit none

  integer                                  :: i, j, cd, an, float_part, n_l_particles
  real*8                                   :: z,zc,zp
  real*8                                   :: pos_z_w1, pos_z_w2
  integer, dimension(:), allocatable       :: valid_part

  ! Initalizing variables
  z   = 0
  zc  = 0
  zp  = 0
  float_part = 0
  n_l_particles = 0

  pos_z_w1 = TAB_PLAN(1)%center(3)
  pos_z_w2 = TAB_PLAN(2)%center(3)

  if (allocated(valid_part)) deallocate(valid_part)
  allocate(valid_part(n_particles))

  valid_part(:) = 1

  do i=1, nb_ligneCONTACT
    cd = TAB_CONTACT(i)%icdent
    an = TAB_CONTACT(i)%ianent

    if (an .le. n_particles) then
      if(TAB_CONTACT(i)%nature == 'SPPLx' .or. TAB_SPHERE(cd)%nctc .lt. 2) then
        valid_part(cd) = 0
      end if
      if(TAB_CONTACT(i)%nature == 'SPPLx' .or. TAB_SPHERE(an)%nctc .lt. 2) then
        valid_part(an) = 0
      end if
    else
      if(TAB_CONTACT(i)%nature == 'SPPLx' .or. TAB_SPHERE(cd)%nctc .gt. 2) then
        valid_part(cd) = 0.
      end if
    end if
  end do

  ! Computing the number of contacts
  do i=1, nb_ligneCONTACT
    ! Only between particles
    if (TAB_CONTACT(i)%nature == 'SPPLx') cycle
    cd = TAB_CONTACT(i)%icdent
    an = TAB_CONTACT(i)%ianent
    if ((cd .gt. n_particles) .or. (an .gt. n_particles)) cycle
    if (valid_part(an) .lt. 1 .and. valid_part(cd) .lt. 1) cycle
    ! Computing the number of contacts if they are active

    ! Counting a smaller box letting a space to the walls
    !if (TAB_SPHERE(cd)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
    !if (TAB_SPHERE(cd)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle

    !if (TAB_SPHERE(an)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
    !if (TAB_SPHERE(an)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle

    if (TAB_CONTACT(i)%rn .gt. 0) then
      zc=zc+1
    end if
  end do

  ! Computing the number of floating particles
  !do i=1,n_particles
    !if (TAB_SPHERE(i)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
    !if (TAB_SPHERE(i)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle

    !n_l_particles = n_l_particles + 1

    !if (TAB_SPHERE(i)%nctc .lt. 2) float_part = float_part + 1
  !end do

  !zp = n_l_particles - float_part
  do i=1, n_particles
    if (valid_part(i) == 1) then
      zp = zp + 1
    end if
  end do

 !print*, n_l_particles
 !print*, float_part
 !print*, zp
 !print*, zc

  z   = 2*zc/zp

  print*, z

  if (first_over_all) then
    write(100,*) '#   time     ', '      z      '
  end if

  write(100,'(2(1X,E12.5))') time, z

  print*, 'Write Coordination  ---> Ok!'
end subroutine nb_coordination


!==============================================================================
! Calcul de q/p
!==============================================================================
subroutine qoverp

  implicit none

  integer                                  :: i, cd, an
  integer                                  :: ierror, matz, lda
  real*8                                   :: Rtik, Rnik, Rsik
  real*8                                   :: S1,S2,S3,Rayon_max,q,p,dmean
  real*8,dimension(3)                      :: nik,tik,sik,Lik,Fik
  real*8,dimension(3)                      :: wr,wi
  real*8,dimension(3,3)                    :: Moment,M
  real*8,dimension(3,3)                    :: localframe

  real*8                                   :: pos_z_w1, pos_z_w2


  ! Initalizing variables
  Moment(:,:)   = 0.D0
  Lik = 0.0
  Fik = 0.0

  pos_z_w1 = TAB_PLAN(1)%center(3)
  pos_z_w2 = TAB_PLAN(2)%center(3)

  do i=1,nb_ligneCONTACT
    ! Contacts between particles
    if (TAB_CONTACT(i)%nature == 'SPPLx') cycle

    ! Active contacts
    if (TAB_CONTACT(i)%rn .le. 0.0) cycle

    cd   = TAB_CONTACT(i)%icdent
    an   = TAB_CONTACT(i)%ianent

    if ((cd .gt. n_particles) .or. (an .gt. n_particles)) cycle

    ! Counting a smaller box letting a space to the walls
    if (TAB_SPHERE(cd)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
    if (TAB_SPHERE(cd)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle

    if (TAB_SPHERE(an)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
    if (TAB_SPHERE(an)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle

    nik  = TAB_CONTACT(i)%n
    tik  = TAB_CONTACT(i)%t
    sik  = TAB_CONTACT(i)%s
    Rtik = TAB_CONTACT(i)%rt
    Rnik = TAB_CONTACT(i)%rn
    Rsik = TAB_CONTACT(i)%rs

    ! The branch
    Lik(1:3) = TAB_SPHERE(cd)%center(1:3)-TAB_SPHERE(an)%center(1:3)

    if (sqrt(Lik(1)**2 + Lik(2)**2 + Lik(3)**2) > 3*(TAB_SPHERE(cd)%Rmax + TAB_SPHERE(an)%Rmax)) then
      !Rayon_max = sqrt((Lik(1)**2 + Lik(2)**2 + Lik(3)**2))
      !Lik(:) = Lik(:) / Rayon_max
      Lik(:) = abs(TAB_SPHERE(cd)%Rmax+TAB_SPHERE(an)%Rmax)*nik(:)
    end if

    ! The force
    Fik(1:3) = (Rnik*nik(1:3)+Rtik*tik(1:3)+Rsik*sik(1:3))

    Moment(1,1:3) = Fik(1)*Lik(1:3) + Moment(1,1:3)
    Moment(2,1:3) = Fik(2)*Lik(1:3) + Moment(2,1:3)
    Moment(3,1:3) = Fik(3)*Lik(1:3) + Moment(3,1:3)
  end do

  ! The stress tensor
  Moment = Moment / (large * width * height)

  ! Computing the principal stresses
  M = Moment
  lda  = 3
  matz = 1
  call rg ( lda, 3, M, wr, wi, matz, localframe, ierror )

  if ( wr(1)==wr(2) .and. (wr(2)==wr(3)) ) then
    S1=wr(1)
    S2=S1
    S3=S1
  else
    S3 = max( wr(1),max(wr(2),wr(3)))
    S1 = min( wr(1),min(wr(2),wr(3)))
    if (wr(1)==wr(2)) then
      if (S1==wr(1)) S2=S1
      if (S3==wr(1)) S2=S3
    else if (wr(3)==wr(2)) then
      if (S1==wr(2)) S2=S1
      if (S3==wr(3)) S2=S3
    else
      do i=1,3
        if ((wr(i)<S3) .and. (wr(i)>S1)) S2=wr(i)
      end do
    end if
  end if

  q = 0.5*(S3-S1)
  !p = (S1+S2+S3)/3
  p = 0.5*(S1+S3)

  if (first_over_all) then
    write(103,*) '#   time     ', '      S1     ', '      S2     ', '      S3     ', '     q/p     '
  end if

  write(103,'(5(1X,E12.5))') time, S1, S2, S3, q/p

  print*, 'Write qoverp        ---> Ok!'

end subroutine qoverp


!==============================================================================
! Calcul de l'anisotropie des contacts
!==============================================================================
subroutine anisotropy_contact

  implicit none

  integer                                  :: i,cd,an
  integer                                  :: ierror,matz,lda
  real*8                                   :: S1, S2, S3, cpt, Rnik, ac
  real*8,dimension(3,3)                    :: Fabric
  real*8,dimension(3)                      :: nik
  real*8,dimension(3)                      :: wr,wi
  real*8,dimension(3,3)                    :: localframe

  real*8                                   :: pos_z_w1, pos_z_w2

  ! Initializing the number of active contacts
  cpt = 0

  ! Initializing the fabric tensor
  Fabric(:,:) = 0.0

  pos_z_w1 = TAB_PLAN(1)%center(3)
  pos_z_w2 = TAB_PLAN(2)%center(3)

  ! Building the fabric tensor
  do i=1,nb_ligneCONTACT
    ! Checking it is a contact between two bodies
    if(TAB_CONTACT(i)%nature == 'SPPLx') cycle
    cd   = TAB_CONTACT(i)%icdent
    an   = TAB_CONTACT(i)%ianent

    if ((cd .gt. n_particles) .or. (an .gt. n_particles)) cycle

    ! Counting a smaller box letting a space to the walls
    if (TAB_SPHERE(cd)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
    if (TAB_SPHERE(cd)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle

    if (TAB_SPHERE(an)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
    if (TAB_SPHERE(an)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle

    Rnik = TAB_CONTACT(i)%rn

    ! Checking if it is an active contact
    if (Rnik .le. 0.D0) cycle

    nik  = TAB_CONTACT(i)%n

    Fabric(1,1:3) = Fabric(1,1:3) + nik(1)*nik(1:3)
    Fabric(2,1:3) = Fabric(2,1:3) + nik(2)*nik(1:3)
    Fabric(3,1:3) = Fabric(3,1:3) + nik(3)*nik(1:3)

    cpt = cpt + 1
  end do

  ! Normalizing by the number of contacts
  Fabric = Fabric / cpt

  ! Finding the eigen values
  lda  = 3
  matz = 1

  call rg (lda, 3, Fabric, wr, wi, matz, localframe, ierror)

  if ( wr(1)==wr(2) .and. (wr(2)==wr(3)) ) then
    S1=wr(1)
    S2=S1
    S3=S1
  else
    S3 = max( wr(1),max(wr(2),wr(3)))
    S1 = min( wr(1),min(wr(2),wr(3)))
    if (wr(1)==wr(2)) then
      if (S1==wr(1)) S2=S1
      if (S3==wr(1)) S2=S3
    else if (wr(3)==wr(2)) then
      if (S1==wr(2)) S2=S1
      if (S3==wr(3)) S2=S3
    else
      do i=1,3
        if ((wr(i)<S3) .and. (wr(i)>S1)) S2=wr(i)
      end do
    end if
  end if

  ! Computing ac from the principal values of F

  ! We do not need to divide as the trace of F is 1 ... or do we?

  ac = 2*(S3-S1)/(S3+S1)

  if (first_over_all) then
    write(104,*) '#   time     ', '     S1      ', '     S2      ', '     S3      ', '     ac      '
  end if

  write(104,'(5(1X,E12.5))') time, S1, S2, S3, ac

  print*, 'Write Anisotropy Conctacts      ---> Ok!'

end subroutine anisotropy_contact


!==============================================================================
! Calcul de l'anisotropie des contacts
!==============================================================================
subroutine anisotropy_force

  implicit none

  integer                                  :: i,cd,an,j
  integer                                  :: ierror,matz,lda
  real*8                                   :: X1n, X2n, X3n, X1f, X2f, X3f
  real*8                                   :: cpt, Rtik, Rnik, Rsik, av_force, norm_FN, norm_FT
  real*8, dimension(3)                     :: wrn, win, wrf, wif
  real*8, dimension(3)                     :: nik,tik, sik, dirFT
  real*8, dimension(3)                     :: Fik, FikT
  real*8, dimension(3,3)                   :: Fabric_N,Fabric_T, Fabric_F
  real*8, dimension(3,3)                   :: localframeN, localframeF

  real*8                                   :: pos_z_w1, pos_z_w2

  ! Initializing variables
  Fabric_N(:,:) = 0.
  Fabric_T(:,:) = 0.
  Fabric_F(:,:) = 0.
  cpt = 0.
  av_force = 0.
  Fik(:)   = 0.
  FikT(:)  = 0.

  pos_z_w1 = TAB_PLAN(1)%center(3)
  pos_z_w2 = TAB_PLAN(2)%center(3)

  ! Building the fabric forces tensor
  do i=1,nb_ligneCONTACT
    if  (TAB_CONTACT(i)%nature == 'SPPLx') cycle
    cd   = TAB_CONTACT(i)%icdent
    an   = TAB_CONTACT(i)%ianent

    if ((cd .gt. n_particles) .or. (an .gt. n_particles)) cycle

    ! Counting a smaller box letting a space to the walls
    if (TAB_SPHERE(cd)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
    if (TAB_SPHERE(cd)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle

    if (TAB_SPHERE(an)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
    if (TAB_SPHERE(an)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle

    nik  = TAB_CONTACT(i)%n
    tik  = TAB_CONTACT(i)%t
    sik  = TAB_CONTACT(i)%s
    Rnik = TAB_CONTACT(i)%rn
    Rtik = TAB_CONTACT(i)%rt
    Rsik = TAB_CONTACT(i)%rs

    ! Only active contacts
    if (Rnik .le. 0.D0) cycle

    !print*, "First"
    !print*, Rnik

    ! Active contacts
    cpt = cpt + 1

    ! The force vector
    Fik(:) = Rnik*nik(:) + Rtik*tik(:) + Rsik*sik(:)

    ! Average normal force
    !av_force = av_force + Fik(1)*nik(1) + Fik(2)*nik(2) + Fik(3)*nik(3)
    ! Average total force
    av_force = av_force + sqrt(Fik(1)**2 + Fik(2)**2 +Fik(3)**2)

    ! The norm of the normal force
    norm_FN = Fik(1)*nik(1) + Fik(2)*nik(2) + Fik(3)*nik(3)

    !print*, norm_FN

    ! Building the tensor of normal forces
    Fabric_N(1,1:3) = norm_FN*nik(1)*nik(1:3) + Fabric_N(1,1:3)
    Fabric_N(2,1:3) = norm_FN*nik(2)*nik(1:3) + Fabric_N(2,1:3)
    Fabric_N(3,1:3) = norm_FN*nik(3)*nik(1:3) + Fabric_N(3,1:3)

    ! The tangential force vector
    FikT(1) = Fik(1) - Rnik*nik(1)
    FikT(2) = Fik(2) - Rnik*nik(2)
    FikT(3) = Fik(3) - Rnik*nik(3)

    ! The norm of FikT
    norm_FT = sqrt(FikT(1)**2 + FikT(2)**2 + FikT(3)**2)

    ! The direction of FikT
    dirFT = FikT/norm_FT

    ! Building the tensor of tangential forces
    Fabric_T(1,1:3) = norm_FT*nik(1)*dirFT(1:3) + Fabric_T(1,1:3)
    Fabric_T(2,1:3) = norm_FT*nik(2)*dirFT(1:3) + Fabric_T(2,1:3)
    Fabric_T(3,1:3) = norm_FT*nik(3)*dirFT(1:3) + Fabric_T(3,1:3)
  end do

  ! Computing the average normal force
  av_force = av_force / cpt

  ! Computing the fabric forces tensors
  Fabric_N = Fabric_N / cpt !(cpt*av_force)
  Fabric_T = Fabric_T / cpt !(cpt*av_force)

  ! Building the full fabric tensor of forces (Built this always before using rg!!!)
  Fabric_F = Fabric_N + Fabric_T

  !print*, Fabric_N(1,1) + Fabric_N(2,2) + Fabric_N(3,3)

  !print*, Fabric_N
  !print*, 'asd1'
  !print*, Fabric_T
  !print*, 'asd2'
  !print*, Fabric_F
  !print*, 'asd3'

  ! Finding the eigen values of the fabric normal forces tensor
  lda = 3
  matz = 1

  call rg(lda, 3, Fabric_N, wrn, win, matz, localframeN, ierror)

  if ( wrn(1)==wrn(2) .and. (wrn(2)==wrn(3)) ) then
    X1n=wrn(1)
    X2n=X1n
    X3n=X1n
  else
    X3n = max( wrn(1),max(wrn(2),wrn(3)))
    X1n = min( wrn(1),min(wrn(2),wrn(3)))
    if (wrn(1)==wrn(2)) then
      if (X1n==wrn(1)) X2n=X1n
      if (X3n==wrn(1)) X2n=X3n
    else if (wrn(3)==wrn(2)) then
      if (X1n==wrn(2)) X2n=X1n
      if (X3n==wrn(3)) X2n=X3n
    else
      do i=1,3
        if ((wrn(i)<X3n) .and. (wrn(i)>X1n)) X2n=wrn(i)
      end do
    end if
  end if

  ! Finding the eigen values of the full fabric forces tensor
  call rg(lda, 3, Fabric_F, wrf, wif, matz, localframeF, ierror)

  if ( wrf(1)==wrf(2) .and. (wrf(2)==wrf(3)) ) then
    X1f=wrf(1)
    X2f=X1f
    X3f=X1f
  else
    X3f = max(wrf(1),max(wrf(2),wrf(3)))
    X1f = min(wrf(1),min(wrf(2),wrf(3)))
    if (wrf(1)==wrf(2)) then
      if (X1f==wrf(1)) X2f=X1f
      if (X3f==wrf(1)) X2f=X3f
    else if (wrf(3)==wrf(2)) then
      if (X1f==wrf(2)) X2f=X1f
      if (X3f==wrf(3)) X2f=X3f
    else
      do i=1,3
        if ((wrf(i)<X3f) .and. (wrf(i)>X1f)) X2f=wrf(i)
      end do
    end if
  end if

  if (first_over_all) then
     write(105,*) '#   time    ', ' 2(X3-1)/(X1+3)n', ' 2(X3-1)/(X1+3)f '
  end if

  !print*, X1f,  X2f, X3f

  write(105,'(3(1X,E12.5))') time, 2*(X3n-X1n)/(X1n+X3n), 2*(X3f-X1f)/(X1f+X3f)

  print*, 'Write Anisotropy Force          ---> Ok!'
end subroutine anisotropy_force


!==============================================================================
! Calcul de l'anisotropie des branches
!==============================================================================
subroutine anisotropy_branch

  implicit none

  integer                                  :: i,cd,an,j
  integer                                  :: ierror,matz,lda
  real*8                                   :: X1n, X2n, X3n, X1f, X2f, X3f, Rayon_max
  real*8                                   :: cpt, av_length, av_total, Lnik, Ltik, Rnik, Rsik, Rtik
  real*8, dimension(3)                     :: nik, tik, sik, ikt
  real*8, dimension(3)                     :: wrn, win, wrf, wif
  real*8, dimension(3)                     :: Lik
  real*8, dimension(3,3)                   :: Fabric_N,Fabric_T, Fabric_F
  real*8, dimension(3,3)                   :: localframeN, localframeF

  real*8                                   :: pos_z_w1, pos_z_w2

  ! Initializing variables
  Fabric_N(:,:) = 0.
  Fabric_T(:,:) = 0.
  Fabric_F(:,:) = 0.
  cpt = 0.
  av_length = 0.
  av_total = 0.
  Lik(:) = 0.

  pos_z_w1 = TAB_PLAN(1)%center(3)
  pos_z_w2 = TAB_PLAN(2)%center(3)

  ! Building the chi tensor with the length of branches
  do i=1,nb_ligneCONTACT
    if  (TAB_CONTACT(i)%nature == 'SPPLx') cycle
    cd   = TAB_CONTACT(i)%icdent
    an   = TAB_CONTACT(i)%ianent
    nik  = TAB_CONTACT(i)%n
    tik  = TAB_CONTACT(i)%t
    sik  = TAB_CONTACT(i)%s
    Rnik = TAB_CONTACT(i)%rn
    Rsik = TAB_CONTACT(i)%rs
    Rtik = TAB_CONTACT(i)%rt

    ! Only active contacts
    if (Rnik .le. 0.D0) cycle

    if ((cd .gt. n_particles) .or. (an .gt. n_particles)) cycle

    ! Counting a smaller box letting a space to the walls
    if (TAB_SPHERE(cd)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
    if (TAB_SPHERE(cd)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle

    if (TAB_SPHERE(an)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
    if (TAB_SPHERE(an)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle

    ! Active contacts
    cpt = cpt + 1

    ! The branch
    Lik(1:3) = TAB_SPHERE(cd)%center(1:3)-TAB_SPHERE(an)%center(1:3)

    if (sqrt(Lik(1)**2 + Lik(2)**2 + Lik(3)**2) > 3*(TAB_SPHERE(cd)%Rmax + TAB_SPHERE(an)%Rmax)) then
      !Rayon_max = sqrt((Lik(1)**2 + Lik(2)**2 + Lik(3)**2))
      !Lik(:) = Lik(:) / Rayon_max
      Lik(:) = abs(TAB_SPHERE(cd)%Rmax+TAB_SPHERE(an)%Rmax)*nik(:)
    end if

    ! Average branch length in the normal contact direction
    !av_length = av_length + Lik(1)*nik(1) + Lik(2)*nik(2) + Lik(3)*nik(3)

    ! Average total length
    !av_total = av_total + sqrt(Lik(1)**2 + Lik(2)**2 + Lik(3)**2)

    ! The norm of the branch in the normal contact direction
    Lnik = Lik(1)*nik(1) + Lik(2)*nik(2) + Lik(3)*nik(3)

    ! The direction of the resultant for tangent forces and Normalizing
    ikt = tik*Rtik + sik*Rsik

    ! Computing the unitary tangential vector
    ikt = ikt/sqrt(ikt(1)**2 + ikt(2)**2 + ikt(3)**2)

    ! The branch in the tangential resultant
    Ltik = Lik(1)*ikt(1) + Lik(2)*ikt(2) + Lik(3)*ikt(3)

    ! Building the chi_n
    Fabric_N(1,1:3) = Lnik*nik(1)*nik(1:3) + Fabric_N(1,1:3)
    Fabric_N(2,1:3) = Lnik*nik(2)*nik(1:3) + Fabric_N(2,1:3)
    Fabric_N(3,1:3) = Lnik*nik(3)*nik(1:3) + Fabric_N(3,1:3)

    ! Building the tensor of tangential forces
    Fabric_T(1,1:3) = Ltik*nik(1)*ikt(1:3) + Fabric_T(1,1:3)
    Fabric_T(2,1:3) = Ltik*nik(2)*ikt(1:3) + Fabric_T(2,1:3)
    Fabric_T(3,1:3) = Ltik*nik(3)*ikt(1:3) + Fabric_T(3,1:3)
  end do

  ! Computing the average normal force
  !av_length = av_length / cpt

  ! Computing the fabric forces tensors
  !Fabric_N = Fabric_N / (cpt*av_length)
  Fabric_N = Fabric_N / cpt
  !Fabric_T = Fabric_T / (cpt*av_length)
  Fabric_T = Fabric_T / cpt

  ! Building the full fabric forces tensor (Built this always before using rg!!!)
  Fabric_F = Fabric_N + Fabric_T

  ! Finding the eigen values of the fabric normal forces tensor
  lda = 3
  matz = 1

  call rg(lda, 3, Fabric_N, wrn, win, matz, localframeN, ierror)

  if ( wrn(1)==wrn(2) .and. (wrn(2)==wrn(3)) ) then
    X1n=wrn(1)
    X2n=X1n
    X3n=X1n
  else
    X3n = max( wrn(1),max(wrn(2),wrn(3)))
    X1n = min( wrn(1),min(wrn(2),wrn(3)))
    if (wrn(1)==wrn(2)) then
      if (X1n==wrn(1)) X2n=X1n
      if (X3n==wrn(1)) X2n=X3n
    else if (wrn(3)==wrn(2)) then
      if (X1n==wrn(2)) X2n=X1n
      if (X3n==wrn(3)) X2n=X3n
    else
      do i=1,3
        if ((wrn(i)<X3n) .and. (wrn(i)>X1n)) X2n=wrn(i)
      end do
    end if
  end if

  ! Finding the eigen values of the full fabric forces tensor
  call rg(lda, 3, Fabric_F, wrf, wif, matz, localframeF, ierror)

  if ( wrf(1)==wrf(2) .and. (wrf(2)==wrf(3)) ) then
    X1f=wrf(1)
    X2f=X1f
    X3f=X1f
  else
    X3f = max(wrf(1),max(wrf(2),wrf(3)))
    X1f = min(wrf(1),min(wrf(2),wrf(3)))
    if (wrf(1)==wrf(2)) then
      if (X1f==wrf(1)) X2f=X1f
      if (X3f==wrf(1)) X2f=X3f
    else if (wrf(3)==wrf(2)) then
      if (X1f==wrf(2)) X2f=X1f
      if (X3f==wrf(3)) X2f=X3f
    else
      do i=1,3
        if ((wrf(i)<X3f) .and. (wrf(i)>X1f)) X2f=wrf(i)
      end do
    end if
  end if

  if (first_over_all) then
     write(106,*) '#   time     ', '1/2(X3-1)/(X1+3)n', ' 1/2(X3-1)/(X1+3)f '
  end if

  write(106,'(3(1X,E12.5))') time, 2*(X3n-X1n)/(X1n+X3n), 2*(X3f-X1f)/(X1f+X3f)

  print*, 'Write Anisotropy Branch          ---> Ok!'
end subroutine anisotropy_branch


!==============================================================================
! Computing the compacity
!==============================================================================
subroutine float_particles

  implicit none

  integer                                  :: i, j, float_part, real_particles, cd, an
  real*8                                   :: float_ratio
  logical                                  :: has_contact

  real*8                                   :: pos_z_w1, pos_z_w2

  pos_z_w1 = TAB_PLAN(1)%center(3)
  pos_z_w2 = TAB_PLAN(2)%center(3)

  ! Initializing
  float_part = 0
  real_particles = 0

  do i=1, n_particles
    ! Particles far from the walls
    if (TAB_SPHERE(i)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
    if (TAB_SPHERE(i)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle

    real_particles = real_particles + 1

    if (TAB_SPHERE(i)%nctc == 0.0) then
      float_part = float_part + 1
    end if
  end do

  print*, float_part
  print*, real_particles

  ! Normalizing by the number of particles in the box
  float_ratio = real(float_part) / real(real_particles)

  if (first_over_all) then
    write(112,*) '#   time     ', ' float_ratio  '
  end if

  write(112,'(2(1X,E12.5))') time, float_ratio

  write (*,*) 'Write Floating particles ---> Ok!'

end subroutine float_particles

!==============================================================================
! Computing the coordination by class
!==============================================================================
subroutine coordination_class

  implicit none

  integer                                  :: i, j, cd, an, n_bin
  real*8                                   :: r_max, r_min, interval
  real*8, dimension(:,:), allocatable      :: bins
  integer, dimension(:), allocatable       :: valid_part
  character(len=24)                        :: coor_c_file

  real*8                                   :: pos_z_w1, pos_z_w2

  pos_z_w1 = TAB_PLAN(1)%center(3)
  pos_z_w2 = TAB_PLAN(2)%center(3)

  ! The name of the file
  coor_c_file = './POSTPRO/CORCLASS.     '

  ! Preparing the corresponding number of the file
  if (compteur_clout<10) then
    WRITE(coor_c_file(20:21),'(I1)')   compteur_clout
  else if ( (compteur_clout>=10) .and. (compteur_clout<100) ) then
    WRITE(coor_c_file(20:22),'(I2)')   compteur_clout
  else if ( (compteur_clout>=100).and. (compteur_clout<1000) ) then
    WRITE(coor_c_file(20:23),'(I3)')   compteur_clout
  else if ( (compteur_clout>=1000).and. (compteur_clout<10000) ) then
    WRITE(coor_c_file(20:24),'(I4)')   compteur_clout
  else
    print*, "Not enough file names :: Coordination per class"
    stop
  end if

  ! Opening the file
  open(unit=113,file=coor_c_file,status='replace')

  ! Initalizing
  r_max = -9999.
  r_min =  9999.

  ! Finding the max and min radius
  do i=1, n_particles
    ! Particles far from the walls
    if (TAB_SPHERE(i)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
    if (TAB_SPHERE(i)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle

    r_min = min(r_min, TAB_SPHERE(i)%Rmax)
    r_max = max(r_max, TAB_SPHERE(i)%Rmax)
  end do

  if (allocated(valid_part)) deallocate(valid_part)
  allocate(valid_part(n_particles))

  valid_part(:) = 1

  do i=1, nb_ligneCONTACT
    cd = TAB_CONTACT(i)%icdent
    an = TAB_CONTACT(i)%ianent

    if (an .le. n_particles) then
      if(TAB_CONTACT(i)%nature == 'SPPLx' .or. TAB_SPHERE(cd)%nctc .lt. 2) then
        valid_part(cd) = 0
      end if
      if(TAB_CONTACT(i)%nature == 'SPPLx' .or. TAB_SPHERE(an)%nctc .lt. 2) then
        valid_part(an) = 0
      end if
    else
      if(TAB_CONTACT(i)%nature == 'SPPLx' .or. TAB_SPHERE(cd)%nctc .gt. 2) then
        valid_part(cd) = 0.
      end if
    end if
  end do

  ! The number of bins - > Intervals in the size distribution
  n_bin = 36

  ! The size of the interval
  interval = (r_max - r_min)/n_bin

  ! Allocating the bins
  if (allocated(bins)) deallocate(bins)
  allocate(bins(n_bin,3))
  ! initializing
  bins(:,:) = 0

  ! Adding particles to bins as function of their size. First column: average size.
  ! Second: number of contacts. Third: Total number of particles in this interval

  ! Caracteristic size of each bin 
  do i=1, n_bin
    bins(i,1) = r_min + interval*(i - 0.5)
  end do

  do i=1, n_bin
    do j=1, n_particles
      if((TAB_SPHERE(j)%Rmax .ge. (r_min + (i-1)*interval)) .and. (TAB_SPHERE(j)%Rmax*0.999 .lt. (interval*i + r_min))) then 

        if (TAB_SPHERE(j)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle 
        if (TAB_SPHERE(j)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle 

        if (TAB_SPHERE(j)%nctc .gt. 1) then
          bins(i,2) = bins(i,2) + TAB_SPHERE(j)%nctc
          if (valid_part(j) == 0) cycle
          bins(i,3) = bins(i,3) + 1
        end if 
      end if
    end do
  end do

  ! Writing
  write(113,*) '#   Class     ', '      z      '

  print*, bins

  !print*, bins  
  
  do i=1, n_bin
    if (bins(1,3) .gt. 0.) then 
      write(113,'(2(1X,E12.5))') bins(i,1), real(bins(i,2)/bins(i,3))
    end if 
  end do 

  close(113)

  print*, 'Coordination Class  ---> Ok!'

end subroutine coordination_class


!================================================
! Drawing in vtk
!================================================
subroutine draw

  implicit none

  integer                                  ::  i, j, k, n_l_vertices, sphere_n_vertices, sphere_n_faces
  integer                                  ::  n_l_faces, curr_l_faces
  integer                                  ::  n_l_fields, l_counter
  integer, dimension(:,:), allocatable     ::  sphere_connect
  real*8, dimension(3)                     ::  curr_l_vector
  real*8                                   ::  ave_rad
  real*8, dimension(:,:), allocatable      ::  sphere_vertices, local_sphere
  logical                                  ::  dir_vtk
  character(len=8)                         ::  vtk_c_temp, v_n_vertices
  character(:), allocatable                ::  vtk_counter, vtk_part, vtk_n_points

  ! Variables for the forces
  integer                                  ::  l_cdt, l_ant, f_counter, f_counter_periodic
  real*8                                   ::  vtk_ave_force, l_force_scale, l_force
  real*8                                   ::  l_rn, l_rt, l_rs
  real*8, dimension(3)                     ::  vtk_cd_center, vtk_an_center, v_l_normal, v_l_t, v_l_s
  real*8, dimension(3)                     ::  Lik
  character(len=8)                         ::  v_n_v_forces
  character(:), allocatable                ::  vtk_forces, vtk_n_v_forces

  ! Variables for the friction mobilization
  integer                                  ::  real_ctc
  character(len=8)                         ::  v_n_v_friction
  character(:), allocatable                ::  vtk_friction, vtk_n_v_friction

  real*8                                   :: pos_z_w1, pos_z_w2, disp_wall_times
  
  pos_z_w1 = TAB_PLAN(1)%center(3)
  pos_z_w2 = TAB_PLAN(2)%center(3)

  ! Cleaning or creating the folder if necessary
  if (first_over_all) then
    ! Asking if the file already exists
    inquire(file='./POSTPRO/VTK', exist=dir_vtk)
    if(dir_vtk) then
      ! Cleaning
      call system('rm ./POSTPRO/VTK/*')
    else
      ! Creating
      call system('mkdir ./POSTPRO/VTK')
    end if
  end if

  ! Creating the file name for the particles
  write(vtk_c_temp, '(I8)') compteur_clout
  
  if (compteur_clout<10) then
    vtk_counter = vtk_c_temp(8:8)
  else if (compteur_clout >= 10 .and. compteur_clout < 100) then
    vtk_counter = vtk_c_temp(7:8)
  else if (compteur_clout >= 100 .and. compteur_clout<1000) then
    vtk_counter = vtk_c_temp(6:8)
  else if (compteur_clout >= 1000 .and. compteur_clout<10000) then
    vtk_counter = vtk_c_temp(5:8)
  end if

  ! Computing the average radius in the sample 
  ave_rad = 0.
  do i=1, n_particles
    ave_rad = ave_rad + TAB_SPHERE(i)%Rmax
  end do

  ave_rad = ave_rad/n_particles
  
  ! Writing the particles information 
  vtk_part = './POSTPRO/VTK/rigid_' // vtk_counter // '.vtk'
  
  open(unit=114, file=vtk_part, status='replace')
  write(114,'(A)') '# vtk DataFile Version 3.0'
  write(114,'(A,I6)') 'RIGID ', compteur_clout
  write(114,'(A)') 'ASCII'
  write(114,'(A)') 'DATASET POLYDATA'
  
  ! Writing the number of vertices (Each sphere ---> 74 vertices)
  sphere_n_vertices = 74
  sphere_n_faces = 144
  n_l_vertices = n_particles * sphere_n_vertices
  
  ! Adding the vertices of the walls (Periodic --> 2 walls)
  n_l_vertices = n_l_vertices + (2*8)
  
  write(v_n_vertices, '(I8)') n_l_vertices
  
  vtk_n_points = 'POINTS ' // v_n_vertices // ' float'
  write(114,'(A)') vtk_n_points
  
  ! The sphere is written as a polyhedron of several faces
  if (allocated(sphere_vertices)) deallocate(sphere_vertices)
  allocate(sphere_vertices(sphere_n_vertices,3))

  if (allocated(sphere_connect)) deallocate(sphere_connect)
  allocate(sphere_connect(sphere_n_faces,3))

  ! List of vertices of unitary polyhedron 
  i=1
  sphere_vertices(i ,1:3)   = (/  0.00000000e+00,   0.00000000e+00,   1.00000000e+00/); i=i+1
  sphere_vertices(i ,1:3)   = (/  4.33883739e-01,   0.00000000e+00,   9.00968868e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/  3.75754340e-01,   2.16941870e-01,   9.00968868e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/  2.16941870e-01,   3.75754340e-01,   9.00968868e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/  2.65677166e-17,   4.33883739e-01,   9.00968868e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/ -2.16941870e-01,   3.75754340e-01,   9.00968868e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/ -3.75754340e-01,   2.16941870e-01,   9.00968868e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/ -4.33883739e-01,   5.31354332e-17,   9.00968868e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/ -3.75754340e-01,  -2.16941870e-01,   9.00968868e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/ -2.16941870e-01,  -3.75754340e-01,   9.00968868e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/ -7.97031498e-17,  -4.33883739e-01,   9.00968868e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/  2.16941870e-01,  -3.75754340e-01,   9.00968868e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/  3.75754340e-01,  -2.16941870e-01,   9.00968868e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/  7.81831482e-01,   0.00000000e+00,   6.23489802e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/  6.77085925e-01,   3.90915741e-01,   6.23489802e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/  3.90915741e-01,   6.77085925e-01,   6.23489802e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/  4.78733711e-17,   7.81831482e-01,   6.23489802e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/ -3.90915741e-01,   6.77085925e-01,   6.23489802e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/ -6.77085925e-01,   3.90915741e-01,   6.23489802e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/ -7.81831482e-01,   9.57467422e-17,   6.23489802e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/ -6.77085925e-01,  -3.90915741e-01,   6.23489802e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/ -3.90915741e-01,  -6.77085925e-01,   6.23489802e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/ -1.43620113e-16,  -7.81831482e-01,   6.23489802e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/  3.90915741e-01,  -6.77085925e-01,   6.23489802e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/  6.77085925e-01,  -3.90915741e-01,   6.23489802e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/  9.74927912e-01,   0.00000000e+00,   2.22520934e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/  8.44312339e-01,   4.87463956e-01,   2.22520934e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/  4.87463956e-01,   8.44312339e-01,   2.22520934e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/  5.96971174e-17,   9.74927912e-01,   2.22520934e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/ -4.87463956e-01,   8.44312339e-01,   2.22520934e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/ -8.44312339e-01,   4.87463956e-01,   2.22520934e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/ -9.74927912e-01,   1.19394235e-16,   2.22520934e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/ -8.44312339e-01,  -4.87463956e-01,   2.22520934e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/ -4.87463956e-01,  -8.44312339e-01,   2.22520934e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/ -1.79091352e-16,  -9.74927912e-01,   2.22520934e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/  4.87463956e-01,  -8.44312339e-01,   2.22520934e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/  8.44312339e-01,  -4.87463956e-01,   2.22520934e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/  9.74927912e-01,   0.00000000e+00,  -2.22520934e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/  8.44312339e-01,   4.87463956e-01,  -2.22520934e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/  4.87463956e-01,   8.44312339e-01,  -2.22520934e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/  5.96971174e-17,   9.74927912e-01,  -2.22520934e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/ -4.87463956e-01,   8.44312339e-01,  -2.22520934e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/ -8.44312339e-01,   4.87463956e-01,  -2.22520934e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/ -9.74927912e-01,   1.19394235e-16,  -2.22520934e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/ -8.44312339e-01,  -4.87463956e-01,  -2.22520934e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/ -4.87463956e-01,  -8.44312339e-01,  -2.22520934e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/ -1.79091352e-16,  -9.74927912e-01,  -2.22520934e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/  4.87463956e-01,  -8.44312339e-01,  -2.22520934e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/  8.44312339e-01,  -4.87463956e-01,  -2.22520934e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/  7.81831482e-01,   0.00000000e+00,  -6.23489802e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/  6.77085925e-01,   3.90915741e-01,  -6.23489802e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/  3.90915741e-01,   6.77085925e-01,  -6.23489802e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/  4.78733711e-17,   7.81831482e-01,  -6.23489802e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/ -3.90915741e-01,   6.77085925e-01,  -6.23489802e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/ -6.77085925e-01,   3.90915741e-01,  -6.23489802e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/ -7.81831482e-01,   9.57467422e-17,  -6.23489802e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/ -6.77085925e-01,  -3.90915741e-01,  -6.23489802e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/ -3.90915741e-01,  -6.77085925e-01,  -6.23489802e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/ -1.43620113e-16,  -7.81831482e-01,  -6.23489802e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/  3.90915741e-01,  -6.77085925e-01,  -6.23489802e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/  6.77085925e-01,  -3.90915741e-01,  -6.23489802e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/  4.33883739e-01,   0.00000000e+00,  -9.00968868e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/  3.75754340e-01,   2.16941870e-01,  -9.00968868e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/  2.16941870e-01,   3.75754340e-01,  -9.00968868e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/  2.65677166e-17,   4.33883739e-01,  -9.00968868e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/ -2.16941870e-01,   3.75754340e-01,  -9.00968868e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/ -3.75754340e-01,   2.16941870e-01,  -9.00968868e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/ -4.33883739e-01,   5.31354332e-17,  -9.00968868e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/ -3.75754340e-01,  -2.16941870e-01,  -9.00968868e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/ -2.16941870e-01,  -3.75754340e-01,  -9.00968868e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/ -7.97031498e-17,  -4.33883739e-01,  -9.00968868e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/  2.16941870e-01,  -3.75754340e-01,  -9.00968868e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/  3.75754340e-01,  -2.16941870e-01,  -9.00968868e-01/); i=i+1
  sphere_vertices(i ,1:3)   = (/  1.22464680e-16,   0.00000000e+00,  -1.00000000e+00/)

  ! Connectivity of the 144 faced polyhedron 
  i=1
  sphere_connect(i,1:3) = (/1,   2,  3/); i=i+1
  sphere_connect(i,1:3) = (/1,   3,  4/); i=i+1
  sphere_connect(i,1:3) = (/1,   4,  5/); i=i+1
  sphere_connect(i,1:3) = (/1,   5,  6/); i=i+1
  sphere_connect(i,1:3) = (/1,   6,  7/); i=i+1
  sphere_connect(i,1:3) = (/1,   7,  8/); i=i+1
  sphere_connect(i,1:3) = (/1,   8,  9/); i=i+1
  sphere_connect(i,1:3) = (/1,   9, 10/); i=i+1
  sphere_connect(i,1:3) = (/1,  10, 11/); i=i+1
  sphere_connect(i,1:3) = (/1,  11, 12/); i=i+1
  sphere_connect(i,1:3) = (/2,   1, 13/); i=i+1
  sphere_connect(i,1:3) = (/1,  12, 13/); i=i+1
  sphere_connect(i,1:3) = (/3,   2, 14/); i=i+1
  sphere_connect(i,1:3) = (/4,   3, 15/); i=i+1
  sphere_connect(i,1:3) = (/3,  14, 15/); i=i+1
  sphere_connect(i,1:3) = (/5,   4, 16/); i=i+1
  sphere_connect(i,1:3) = (/4,  15, 16/); i=i+1
  sphere_connect(i,1:3) = (/6,   5, 17/); i=i+1
  sphere_connect(i,1:3) = (/5,  16, 17/); i=i+1
  sphere_connect(i,1:3) = (/7,   6, 18/); i=i+1
  sphere_connect(i,1:3) = (/6,  17, 18/); i=i+1
  sphere_connect(i,1:3) = (/8,   7, 19/); i=i+1
  sphere_connect(i,1:3) = (/7,  18, 19/); i=i+1
  sphere_connect(i,1:3) = (/9,   8, 20/); i=i+1
  sphere_connect(i,1:3) = (/8,  19, 20/); i=i+1
  sphere_connect(i,1:3) = (/9,  20, 21/); i=i+1
  sphere_connect(i,1:3) = (/11, 10, 22/); i=i+1
  sphere_connect(i,1:3) = (/10,  9, 22/); i=i+1
  sphere_connect(i,1:3) = (/9,  21, 22/); i=i+1
  sphere_connect(i,1:3) = (/12, 11, 23/); i=i+1
  sphere_connect(i,1:3) = (/11, 22, 23/); i=i+1
  sphere_connect(i,1:3) = (/13, 12, 24/); i=i+1
  sphere_connect(i,1:3) = (/12, 23, 24/); i=i+1
  sphere_connect(i,1:3) = (/2,  13, 25/); i=i+1
  sphere_connect(i,1:3) = (/14,  2, 25/); i=i+1
  sphere_connect(i,1:3) = (/13, 24, 25/); i=i+1
  sphere_connect(i,1:3) = (/14, 25, 26/); i=i+1
  sphere_connect(i,1:3) = (/15, 14, 27/); i=i+1
  sphere_connect(i,1:3) = (/14, 26, 27/); i=i+1
  sphere_connect(i,1:3) = (/16, 15, 27/); i=i+1
  sphere_connect(i,1:3) = (/17, 16, 28/); i=i+1
  sphere_connect(i,1:3) = (/16, 27, 28/); i=i+1
  sphere_connect(i,1:3) = (/17, 28, 29/); i=i+1
  sphere_connect(i,1:3) = (/18, 17, 30/); i=i+1
  sphere_connect(i,1:3) = (/17, 29, 30/); i=i+1
  sphere_connect(i,1:3) = (/19, 18, 31/); i=i+1
  sphere_connect(i,1:3) = (/18, 30, 31/); i=i+1
  sphere_connect(i,1:3) = (/20, 19, 31/); i=i+1
  sphere_connect(i,1:3) = (/20, 31, 32/); i=i+1
  sphere_connect(i,1:3) = (/22, 21, 33/); i=i+1
  sphere_connect(i,1:3) = (/21, 20, 33/); i=i+1
  sphere_connect(i,1:3) = (/20, 32, 33/); i=i+1
  sphere_connect(i,1:3) = (/23, 22, 34/); i=i+1
  sphere_connect(i,1:3) = (/22, 33, 34/); i=i+1
  sphere_connect(i,1:3) = (/23, 34, 35/); i=i+1
  sphere_connect(i,1:3) = (/25, 24, 36/); i=i+1
  sphere_connect(i,1:3) = (/24, 23, 36/); i=i+1
  sphere_connect(i,1:3) = (/23, 35, 36/); i=i+1
  sphere_connect(i,1:3) = (/26, 25, 37/); i=i+1
  sphere_connect(i,1:3) = (/25, 36, 37/); i=i+1
  sphere_connect(i,1:3) = (/27, 26, 38/); i=i+1
  sphere_connect(i,1:3) = (/26, 37, 38/); i=i+1
  sphere_connect(i,1:3) = (/28, 27, 39/); i=i+1
  sphere_connect(i,1:3) = (/27, 38, 39/); i=i+1
  sphere_connect(i,1:3) = (/29, 28, 40/); i=i+1
  sphere_connect(i,1:3) = (/28, 39, 40/); i=i+1
  sphere_connect(i,1:3) = (/30, 29, 41/); i=i+1
  sphere_connect(i,1:3) = (/29, 40, 41/); i=i+1
  sphere_connect(i,1:3) = (/31, 30, 42/); i=i+1
  sphere_connect(i,1:3) = (/30, 41, 42/); i=i+1
  sphere_connect(i,1:3) = (/32, 31, 43/); i=i+1
  sphere_connect(i,1:3) = (/31, 42, 43/); i=i+1
  sphere_connect(i,1:3) = (/33, 32, 44/); i=i+1
  sphere_connect(i,1:3) = (/32, 43, 44/); i=i+1
  sphere_connect(i,1:3) = (/34, 33, 45/); i=i+1
  sphere_connect(i,1:3) = (/33, 44, 45/); i=i+1
  sphere_connect(i,1:3) = (/35, 34, 46/); i=i+1
  sphere_connect(i,1:3) = (/34, 45, 46/); i=i+1
  sphere_connect(i,1:3) = (/36, 35, 47/); i=i+1
  sphere_connect(i,1:3) = (/35, 46, 47/); i=i+1
  sphere_connect(i,1:3) = (/37, 36, 48/); i=i+1
  sphere_connect(i,1:3) = (/36, 47, 48/); i=i+1
  sphere_connect(i,1:3) = (/38, 37, 49/); i=i+1
  sphere_connect(i,1:3) = (/37, 48, 49/); i=i+1
  sphere_connect(i,1:3) = (/38, 49, 50/); i=i+1
  sphere_connect(i,1:3) = (/39, 38, 51/); i=i+1
  sphere_connect(i,1:3) = (/38, 50, 51/); i=i+1
  sphere_connect(i,1:3) = (/40, 39, 52/); i=i+1
  sphere_connect(i,1:3) = (/39, 51, 52/); i=i+1
  sphere_connect(i,1:3) = (/41, 40, 52/); i=i+1
  sphere_connect(i,1:3) = (/41, 52, 53/); i=i+1
  sphere_connect(i,1:3) = (/42, 41, 54/); i=i+1
  sphere_connect(i,1:3) = (/41, 53, 54/); i=i+1
  sphere_connect(i,1:3) = (/43, 42, 54/); i=i+1
  sphere_connect(i,1:3) = (/43, 54, 55/); i=i+1
  sphere_connect(i,1:3) = (/44, 43, 56/); i=i+1
  sphere_connect(i,1:3) = (/43, 55, 56/); i=i+1
  sphere_connect(i,1:3) = (/46, 45, 57/); i=i+1
  sphere_connect(i,1:3) = (/45, 44, 57/); i=i+1
  sphere_connect(i,1:3) = (/44, 56, 57/); i=i+1
  sphere_connect(i,1:3) = (/46, 57, 58/); i=i+1
  sphere_connect(i,1:3) = (/48, 47, 59/); i=i+1
  sphere_connect(i,1:3) = (/47, 46, 59/); i=i+1
  sphere_connect(i,1:3) = (/46, 58, 59/); i=i+1
  sphere_connect(i,1:3) = (/49, 48, 60/); i=i+1
  sphere_connect(i,1:3) = (/48, 59, 60/); i=i+1
  sphere_connect(i,1:3) = (/50, 49, 61/); i=i+1
  sphere_connect(i,1:3) = (/49, 60, 61/); i=i+1
  sphere_connect(i,1:3) = (/51, 50, 62/); i=i+1
  sphere_connect(i,1:3) = (/50, 61, 62/); i=i+1
  sphere_connect(i,1:3) = (/52, 51, 63/); i=i+1
  sphere_connect(i,1:3) = (/51, 62, 63/); i=i+1
  sphere_connect(i,1:3) = (/52, 63, 64/); i=i+1
  sphere_connect(i,1:3) = (/53, 52, 65/); i=i+1
  sphere_connect(i,1:3) = (/52, 64, 65/); i=i+1
  sphere_connect(i,1:3) = (/54, 53, 65/); i=i+1
  sphere_connect(i,1:3) = (/55, 54, 66/); i=i+1
  sphere_connect(i,1:3) = (/54, 65, 66/); i=i+1
  sphere_connect(i,1:3) = (/56, 55, 67/); i=i+1
  sphere_connect(i,1:3) = (/55, 66, 67/); i=i+1
  sphere_connect(i,1:3) = (/57, 56, 68/); i=i+1
  sphere_connect(i,1:3) = (/56, 67, 68/); i=i+1
  sphere_connect(i,1:3) = (/57, 68, 69/); i=i+1
  sphere_connect(i,1:3) = (/59, 58, 70/); i=i+1
  sphere_connect(i,1:3) = (/58, 57, 70/); i=i+1
  sphere_connect(i,1:3) = (/57, 69, 70/); i=i+1
  sphere_connect(i,1:3) = (/60, 59, 71/); i=i+1
  sphere_connect(i,1:3) = (/59, 70, 71/); i=i+1
  sphere_connect(i,1:3) = (/61, 60, 72/); i=i+1
  sphere_connect(i,1:3) = (/60, 71, 72/); i=i+1
  sphere_connect(i,1:3) = (/62, 61, 73/); i=i+1
  sphere_connect(i,1:3) = (/61, 72, 73/); i=i+1
  sphere_connect(i,1:3) = (/63, 62, 74/); i=i+1
  sphere_connect(i,1:3) = (/64, 63, 74/); i=i+1
  sphere_connect(i,1:3) = (/65, 64, 74/); i=i+1
  sphere_connect(i,1:3) = (/66, 65, 74/); i=i+1
  sphere_connect(i,1:3) = (/67, 66, 74/); i=i+1
  sphere_connect(i,1:3) = (/68, 67, 74/); i=i+1
  sphere_connect(i,1:3) = (/69, 68, 74/); i=i+1
  sphere_connect(i,1:3) = (/70, 69, 74/); i=i+1
  sphere_connect(i,1:3) = (/71, 70, 74/); i=i+1
  sphere_connect(i,1:3) = (/72, 71, 74/); i=i+1
  sphere_connect(i,1:3) = (/73, 72, 74/); i=i+1
  sphere_connect(i,1:3) = (/62, 73, 74/)

  ! Writing the coordinates of the vertices of particles and walls
  ! 3 vertices per row
  k=0
  do i=1, n_particles + n_walls
    !
    if (i .le. n_particles) then
      ! The current vertices of this sphere with its corresponding size in global coordinates 
      if (allocated(local_sphere)) deallocate(local_sphere)
      allocate(local_sphere(sphere_n_vertices,3))

      do j=1, sphere_n_vertices
        local_sphere(j,1:3) = (sphere_vertices(j,1:3)*TAB_SPHERE(i)%Rmax) + TAB_SPHERE(i)%center(1:3)
      end do 

      !print*, local_sphere

      do j=1, sphere_n_vertices
        k = k + 1
        if(k .lt. 3) then 
          if (local_sphere(j,1) .lt. 0) then
            write(114,'(F11.8,A)', advance = 'no') local_sphere(j,1), ' '
          else 
            write(114,'(F9.7,A)', advance = 'no') local_sphere(j,1), ' '
          end if
          
          if (local_sphere(j,2) .lt. 0) then
            write(114,'(F11.8,A)', advance = 'no') local_sphere(j,2), ' '
          else 
            write(114,'(F9.7,A)', advance = 'no') local_sphere(j,2), ' '
          end if
          
          if (local_sphere(j,3) .lt. 0) then
            write(114,'(F11.8,A)', advance = 'no') local_sphere(j,3), ' '
          else 
            write(114,'(F9.7,A)', advance = 'no') local_sphere(j,3), ' '
          end if
        else
          if (local_sphere(j,1) .lt. 0) then
            write(114,'(F11.8,A)', advance = 'no') local_sphere(j,1), ' '
          else 
            write(114,'(F9.7,A)', advance = 'no') local_sphere(j,1), ' '
          end if
          
          if (local_sphere(j,2) .lt. 0) then
            write(114,'(F11.8,A)', advance = 'no') local_sphere(j,2), ' '
          else 
            write(114,'(F9.7,A)', advance = 'no') local_sphere(j,2), ' '
          end if
          
          if (local_sphere(j,3) .lt. 0) then
            write(114,'(F11.8,A)') local_sphere(j,3), ' '
          else 
            write(114,'(F9.7,A)') local_sphere(j,3), ' '
          end if
          k = 0
        end if 
      end do
    else
      ! Writing the vertices of the walls
      do j=1, 8
        
        k = k +1
        if(k .lt. 3) then 
          if (TAB_PLAN(i-n_particles)%vertex(1,j) .lt. 0) then
            write(114,'(F11.8,A)', advance = 'no') TAB_PLAN(i-n_particles)%vertex(1,j), ' '
          else 
            write(114,'(F9.7,A)', advance = 'no') TAB_PLAN(i-n_particles)%vertex(1,j), ' '
          end if
          
          if (TAB_PLAN(i-n_particles)%vertex(2,j) .lt. 0) then
            write(114,'(F11.8,A)', advance = 'no') TAB_PLAN(i-n_particles)%vertex(2,j), ' '
          else 
            write(114,'(F9.7,A)', advance = 'no') TAB_PLAN(i-n_particles)%vertex(2,j), ' '
          end if
          
          if (TAB_PLAN(i-n_particles)%vertex(3,j) .lt. 0) then
            write(114,'(F11.8,A)', advance = 'no') TAB_PLAN(i-n_particles)%vertex(3,j), ' '
          else 
            write(114,'(F9.7,A)', advance = 'no') TAB_PLAN(i-n_particles)%vertex(3,j), ' '
          end if
        else
          if (TAB_PLAN(i-n_particles)%vertex(1,j) .lt. 0) then
            write(114,'(F11.8,A)', advance = 'no') TAB_PLAN(i-n_particles)%vertex(1,j), ' '
          else 
            write(114,'(F9.7,A)', advance = 'no') TAB_PLAN(i-n_particles)%vertex(1,j), ' '
          end if
          
          if (TAB_PLAN(i-n_particles)%vertex(2,j) .lt. 0) then
            write(114,'(F11.8,A)', advance = 'no') TAB_PLAN(i-n_particles)%vertex(2,j), ' '
          else 
            write(114,'(F9.7,A)', advance = 'no') TAB_PLAN(i-n_particles)%vertex(2,j), ' '
          end if
          
          if (TAB_PLAN(i-n_particles)%vertex(3,j) .lt. 0) then
            write(114,'(F11.8,A)') TAB_PLAN(i-n_particles)%vertex(3,j), ' '
          else 
            write(114,'(F9.7,A)') TAB_PLAN(i-n_particles)%vertex(3,j), ' '
          end if
          
          k = 0
        end if
      end do
    end if
  end do

  write(114, '(A)') ' '
  
  ! Finding the total number of faces
  n_l_faces = n_particles * sphere_n_faces

  ! Writing the conectivity between vertices
  ! Plus 6 faces each wall (=6*2)
  write(114,'(A)', advance='no') 'POLYGONS '
  write(114, '(2(I8,A))') n_l_faces + 12, ' ' , (n_l_faces*4)+(12*5)
  
  curr_l_faces = 0
  
  do i=1, n_particles
    !write(114,*) 'Particle', i
    do j=1, sphere_n_faces
      ! We have always triangles
      write(114, '(I1,A)', advance= 'no') 3, ' '
      
      ! First face... writing precisely 
      if (sphere_connect(j,1) + curr_l_faces -1 .lt. 10) then
        write(114, '(I1)', advance = 'no') sphere_connect(j,1) + curr_l_faces -1
      else if (sphere_connect(j,1) + curr_l_faces -1 .lt. 100) then
        write(114, '(I2)', advance = 'no') sphere_connect(j,1) + curr_l_faces -1
      else if (sphere_connect(j,1) + curr_l_faces -1 .lt. 1000) then
        write(114, '(I3)', advance = 'no') sphere_connect(j,1) + curr_l_faces -1
      else if (sphere_connect(j,1) + curr_l_faces -1 .lt. 10000) then
        write(114, '(I4)', advance = 'no') sphere_connect(j,1) + curr_l_faces -1
      else if (sphere_connect(j,1) + curr_l_faces -1 .lt. 100000) then
        write(114, '(I5)', advance = 'no') sphere_connect(j,1) + curr_l_faces -1
      else if (sphere_connect(j,1) + curr_l_faces -1 .lt. 1000000) then
        write(114, '(I6)', advance = 'no') sphere_connect(j,1) + curr_l_faces -1
      else if (sphere_connect(j,1) + curr_l_faces -1 .lt. 10000000) then
        write(114, '(I7)', advance = 'no') sphere_connect(j,1) + curr_l_faces -1
      end if
      
      write(114, '(A)', advance='no') ' '
      
      ! Second face
      if (sphere_connect(j,2) + curr_l_faces -1 .lt. 10) then
        write(114, '(I1)', advance = 'no') sphere_connect(j,2) + curr_l_faces -1
      else if (sphere_connect(j,2) + curr_l_faces -1 .lt. 100) then
        write(114, '(I2)', advance = 'no') sphere_connect(j,2) + curr_l_faces -1
      else if (sphere_connect(j,2) + curr_l_faces -1 .lt. 1000) then
        write(114, '(I3)', advance = 'no') sphere_connect(j,2) + curr_l_faces -1
      else if (sphere_connect(j,2) + curr_l_faces -1 .lt. 10000) then
        write(114, '(I4)', advance = 'no') sphere_connect(j,2) + curr_l_faces -1
      else if (sphere_connect(j,2) + curr_l_faces -1 .lt. 100000) then
        write(114, '(I5)', advance = 'no') sphere_connect(j,2) + curr_l_faces -1
      else if (sphere_connect(j,2) + curr_l_faces -1 .lt. 1000000) then
        write(114, '(I6)', advance = 'no') sphere_connect(j,2) + curr_l_faces -1
      else if (sphere_connect(j,2) + curr_l_faces -1 .lt. 10000000) then
        write(114, '(I7)', advance = 'no') sphere_connect(j,2) + curr_l_faces -1
      end if
      
      ! Third face
      write(114, '(A)', advance='no') ' '
      
      if (sphere_connect(j,3) + curr_l_faces -1 .lt. 10) then
        write(114, '(I1)') sphere_connect(j,3) + curr_l_faces -1
      else if (sphere_connect(j,3) + curr_l_faces -1 .lt. 100) then
        write(114, '(I2)') sphere_connect(j,3) + curr_l_faces -1
      else if (sphere_connect(j,3) + curr_l_faces -1 .lt. 1000) then
        write(114, '(I3)') sphere_connect(j,3) + curr_l_faces -1
      else if (sphere_connect(j,3) + curr_l_faces -1 .lt. 10000) then
        write(114, '(I4)') sphere_connect(j,3) + curr_l_faces -1
      else if (sphere_connect(j,3) + curr_l_faces -1 .lt. 100000) then
        write(114, '(I5)') sphere_connect(j,3) + curr_l_faces -1
      else if (sphere_connect(j,3) + curr_l_faces -1 .lt. 1000000) then
        write(114, '(I6)') sphere_connect(j,3) + curr_l_faces -1
      else if (sphere_connect(j,3) + curr_l_faces -1 .lt. 10000000) then
        write(114, '(I7)') sphere_connect(j,3) + curr_l_faces -1
      end if
      
    end do
    curr_l_faces = curr_l_faces + sphere_n_vertices
  end do
  
  ! Writing the conectivity for the walls
  ! For all the walls 
  do i=1, n_walls
    ! Each one with 6 faces
        
    !!!!! 1st (1-2-3-4)
    write(114, '(I1,A)', advance= 'no') 4, ' '
    
    if (curr_l_faces .lt. 10) then
      write(114, '(I1)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 100) then
      write(114, '(I2)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 1000) then
      write(114, '(I3)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 10000) then
      write(114, '(I4)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 100000) then
      write(114, '(I5)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 1000000) then
      write(114, '(I6)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 10000000) then
      write(114, '(I7)', advance = 'no') curr_l_faces
    end if
    
    write(114, '(A)', advance='no') ' '
    
    if (curr_l_faces +1 .lt. 10) then
      write(114, '(I1)', advance = 'no') curr_l_faces +1 
    else if (curr_l_faces +1 .lt. 100) then
      write(114, '(I2)', advance = 'no') curr_l_faces +1 
    else if (curr_l_faces +1 .lt. 1000) then
      write(114, '(I3)', advance = 'no') curr_l_faces +1 
    else if (curr_l_faces +1 .lt. 10000) then
      write(114, '(I4)', advance = 'no') curr_l_faces +1
    else if (curr_l_faces +1 .lt. 100000) then
      write(114, '(I5)', advance = 'no') curr_l_faces +1 
    else if (curr_l_faces +1 .lt. 1000000) then
      write(114, '(I6)', advance = 'no') curr_l_faces +1 
    else if (curr_l_faces +1 .lt. 10000000) then
      write(114, '(I7)', advance = 'no') curr_l_faces +1
    end if
    
    write(114, '(A)', advance='no') ' '
    
    if (curr_l_faces +2 .lt. 10) then
      write(114, '(I1)', advance = 'no') curr_l_faces +2 
    else if (curr_l_faces +2 .lt. 100) then
      write(114, '(I2)', advance = 'no') curr_l_faces +2 
    else if (curr_l_faces +2 .lt. 1000) then
      write(114, '(I3)', advance = 'no') curr_l_faces +2 
    else if (curr_l_faces +2 .lt. 10000) then
      write(114, '(I4)', advance = 'no') curr_l_faces +2
    else if (curr_l_faces +2 .lt. 100000) then
      write(114, '(I5)', advance = 'no') curr_l_faces +2 
    else if (curr_l_faces +2 .lt. 1000000) then
      write(114, '(I6)', advance = 'no') curr_l_faces +2 
    else if (curr_l_faces +2 .lt. 10000000) then
      write(114, '(I7)', advance = 'no') curr_l_faces +2 
    end if
    
    write(114, '(A)', advance='no') ' '
    
    if (curr_l_faces +3 .lt. 10) then
      write(114, '(I1)') curr_l_faces +3
    else if (curr_l_faces +3 .lt. 100) then
      write(114, '(I2)') curr_l_faces +3 
    else if (curr_l_faces +3 .lt. 1000) then
      write(114, '(I3)') curr_l_faces +3 
    else if (curr_l_faces +3 .lt. 10000) then
      write(114, '(I4)') curr_l_faces +3
    else if (curr_l_faces +3 .lt. 100000) then
      write(114, '(I5)') curr_l_faces +3
    else if (curr_l_faces +3 .lt. 1000000) then
      write(114, '(I6)') curr_l_faces +3
    else if (curr_l_faces +3 .lt. 10000000) then
      write(114, '(I7)') curr_l_faces +3
    end if
    
    !!!!! 2nd (5-6-7-8)
    write(114, '(I1,A)', advance= 'no') 4, ' '
    
    if (curr_l_faces + 4 .lt. 10) then
      write(114, '(I1)', advance = 'no') curr_l_faces + 4
    else if (curr_l_faces + 4 .lt. 100) then
      write(114, '(I2)', advance = 'no') curr_l_faces + 4
    else if (curr_l_faces + 4 .lt. 1000) then
      write(114, '(I3)', advance = 'no') curr_l_faces + 4
    else if (curr_l_faces + 4 .lt. 10000) then
      write(114, '(I4)', advance = 'no') curr_l_faces + 4 
    else if (curr_l_faces + 4 .lt. 100000) then
      write(114, '(I5)', advance = 'no') curr_l_faces + 4
    else if (curr_l_faces + 4 .lt. 1000000) then
      write(114, '(I6)', advance = 'no') curr_l_faces + 4 
    else if (curr_l_faces + 4 .lt. 10000000) then
      write(114, '(I7)', advance = 'no') curr_l_faces + 4 
    end if
    
    write(114, '(A)', advance='no') ' '
    
    if (curr_l_faces +5 .lt. 10) then
      write(114, '(I1)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 100) then
      write(114, '(I2)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 1000) then
      write(114, '(I3)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 10000) then
      write(114, '(I4)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 100000) then
      write(114, '(I5)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 1000000) then
      write(114, '(I6)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 10000000) then
      write(114, '(I7)', advance = 'no') curr_l_faces + 5 
    end if
    
    write(114, '(A)', advance='no') ' '
    
    if (curr_l_faces +6 .lt. 10) then
      write(114, '(I1)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 100) then
      write(114, '(I2)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 1000) then
      write(114, '(I3)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 10000) then
      write(114, '(I4)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 100000) then
      write(114, '(I5)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 1000000) then
      write(114, '(I6)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 10000000) then
      write(114, '(I7)', advance = 'no') curr_l_faces + 6
    end if
    
    write(114, '(A)', advance='no') ' '
    
    if (curr_l_faces +7 .lt. 10) then
      write(114, '(I1)') curr_l_faces +7
    else if (curr_l_faces +7 .lt. 100) then
      write(114, '(I2)') curr_l_faces +7 
    else if (curr_l_faces +7 .lt. 1000) then
      write(114, '(I3)') curr_l_faces +7 
    else if (curr_l_faces +7 .lt. 10000) then
      write(114, '(I4)') curr_l_faces +7
    else if (curr_l_faces +7 .lt. 100000) then
      write(114, '(I5)') curr_l_faces +7
    else if (curr_l_faces +7 .lt. 1000000) then
      write(114, '(I6)') curr_l_faces +7
    else if (curr_l_faces +7 .lt. 10000000) then
      write(114, '(I7)') curr_l_faces +7
    end if
    
    !!!!! 3nd (1-5-8-4)
    write(114, '(I1,A)', advance= 'no') 4, ' '
    
    if (curr_l_faces .lt. 10) then
      write(114, '(I1)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 100) then
      write(114, '(I2)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 1000) then
      write(114, '(I3)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 10000) then
      write(114, '(I4)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 100000) then
      write(114, '(I5)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 1000000) then
      write(114, '(I6)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 10000000) then
      write(114, '(I7)', advance = 'no') curr_l_faces
    end if
    
    write(114, '(A)', advance='no') ' '
    
    if (curr_l_faces +4 .lt. 10) then
      write(114, '(I1)', advance = 'no') curr_l_faces + 4
    else if (curr_l_faces +4 .lt. 100) then
      write(114, '(I2)', advance = 'no') curr_l_faces + 4
    else if (curr_l_faces +4 .lt. 1000) then
      write(114, '(I3)', advance = 'no') curr_l_faces + 4
    else if (curr_l_faces +4 .lt. 10000) then
      write(114, '(I4)', advance = 'no') curr_l_faces + 4
    else if (curr_l_faces +4 .lt. 100000) then
      write(114, '(I5)', advance = 'no') curr_l_faces + 4
    else if (curr_l_faces +4 .lt. 1000000) then
      write(114, '(I6)', advance = 'no') curr_l_faces + 4
    else if (curr_l_faces +4 .lt. 10000000) then
      write(114, '(I7)', advance = 'no') curr_l_faces + 4
    end if
    
    write(114, '(A)', advance='no') ' '
    
    if (curr_l_faces +7 .lt. 10) then
      write(114, '(I1)', advance = 'no') curr_l_faces + 7
    else if (curr_l_faces +7 .lt. 100) then
      write(114, '(I2)', advance = 'no') curr_l_faces + 7
    else if (curr_l_faces +7 .lt. 1000) then
      write(114, '(I3)', advance = 'no') curr_l_faces + 7
    else if (curr_l_faces +7 .lt. 10000) then
      write(114, '(I4)', advance = 'no') curr_l_faces + 7
    else if (curr_l_faces +7 .lt. 100000) then
      write(114, '(I5)', advance = 'no') curr_l_faces + 7
    else if (curr_l_faces +7 .lt. 1000000) then
      write(114, '(I6)', advance = 'no') curr_l_faces + 7
    else if (curr_l_faces +7 .lt. 10000000) then
      write(114, '(I7)', advance = 'no') curr_l_faces + 7
    end if
    
    write(114, '(A)', advance='no') ' '
    
    if (curr_l_faces +3 .lt. 10) then
      write(114, '(I1)') curr_l_faces +3
    else if (curr_l_faces +3 .lt. 100) then
      write(114, '(I2)') curr_l_faces +3 
    else if (curr_l_faces +3 .lt. 1000) then
      write(114, '(I3)') curr_l_faces +3 
    else if (curr_l_faces +3 .lt. 10000) then
      write(114, '(I4)') curr_l_faces +3
    else if (curr_l_faces +3 .lt. 100000) then
      write(114, '(I5)') curr_l_faces +3
    else if (curr_l_faces +3 .lt. 1000000) then
      write(114, '(I6)') curr_l_faces +3
    else if (curr_l_faces +3 .lt. 10000000) then
      write(114, '(I7)') curr_l_faces +3
    end if
    
    !!!!! 4th (2-6-7-3)
    write(114, '(I1,A)', advance= 'no') 4, ' '
    
    if (curr_l_faces + 1 .lt. 10) then
      write(114, '(I1)', advance = 'no') curr_l_faces + 1
    else if (curr_l_faces + 1 .lt. 100) then
      write(114, '(I2)', advance = 'no') curr_l_faces + 1
    else if (curr_l_faces + 1 .lt. 1000) then
      write(114, '(I3)', advance = 'no') curr_l_faces + 1
    else if (curr_l_faces + 1 .lt. 10000) then
      write(114, '(I4)', advance = 'no') curr_l_faces + 1
    else if (curr_l_faces + 1 .lt. 100000) then
      write(114, '(I5)', advance = 'no') curr_l_faces + 1
    else if (curr_l_faces + 1 .lt. 1000000) then
      write(114, '(I6)', advance = 'no') curr_l_faces + 1
    else if (curr_l_faces + 1 .lt. 10000000) then
      write(114, '(I7)', advance = 'no') curr_l_faces + 1
    end if
    
    write(114, '(A)', advance='no') ' '
    
    if (curr_l_faces +5 .lt. 10) then
      write(114, '(I1)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 100) then
      write(114, '(I2)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 1000) then
      write(114, '(I3)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 10000) then
      write(114, '(I4)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 100000) then
      write(114, '(I5)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 1000000) then
      write(114, '(I6)', advance = 'no') curr_l_faces + 5 
    else if (curr_l_faces +5 .lt. 10000000) then
      write(114, '(I7)', advance = 'no') curr_l_faces + 5 
    end if
    
    write(114, '(A)', advance='no') ' '
    
    if (curr_l_faces +6 .lt. 10) then
      write(114, '(I1)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 100) then
      write(114, '(I2)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 1000) then
      write(114, '(I3)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 10000) then
      write(114, '(I4)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 100000) then
      write(114, '(I5)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 1000000) then
      write(114, '(I6)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 10000000) then
      write(114, '(I7)', advance = 'no') curr_l_faces + 6
    end if
    
    write(114, '(A)', advance='no') ' '
    
    if (curr_l_faces +2 .lt. 10) then
      write(114, '(I1)') curr_l_faces +2
    else if (curr_l_faces +2 .lt. 100) then
      write(114, '(I2)') curr_l_faces +2
    else if (curr_l_faces +2 .lt. 1000) then
      write(114, '(I3)') curr_l_faces +2
    else if (curr_l_faces +2 .lt. 10000) then
      write(114, '(I4)') curr_l_faces +2
    else if (curr_l_faces +2 .lt. 100000) then
      write(114, '(I5)') curr_l_faces +2
    else if (curr_l_faces +2 .lt. 1000000) then
      write(114, '(I6)') curr_l_faces +2
    else if (curr_l_faces +2 .lt. 10000000) then
      write(114, '(I7)') curr_l_faces +2
    end if
    
    !!!!! 5th (1-2-6-5)
    write(114, '(I1,A)', advance= 'no') 4, ' '
    
    if (curr_l_faces .lt. 10) then
      write(114, '(I1)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 100) then
      write(114, '(I2)', advance = 'no') curr_l_faces
    else if (curr_l_faces  .lt. 1000) then
      write(114, '(I3)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 10000) then
      write(114, '(I4)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 100000) then
      write(114, '(I5)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 1000000) then
      write(114, '(I6)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 10000000) then
      write(114, '(I7)', advance = 'no') curr_l_faces
    end if
    
    write(114, '(A)', advance='no') ' '
    
    if (curr_l_faces +1 .lt. 10) then
      write(114, '(I1)', advance = 'no') curr_l_faces + 1
    else if (curr_l_faces +1 .lt. 100) then
      write(114, '(I2)', advance = 'no') curr_l_faces + 1
    else if (curr_l_faces +1 .lt. 1000) then
      write(114, '(I3)', advance = 'no') curr_l_faces + 1
    else if (curr_l_faces +1 .lt. 10000) then
      write(114, '(I4)', advance = 'no') curr_l_faces + 1
    else if (curr_l_faces +1 .lt. 100000) then
      write(114, '(I5)', advance = 'no') curr_l_faces + 1
    else if (curr_l_faces +1 .lt. 1000000) then
      write(114, '(I6)', advance = 'no') curr_l_faces + 1
    else if (curr_l_faces +1 .lt. 10000000) then
      write(114, '(I7)', advance = 'no') curr_l_faces + 1
    end if
    
    write(114, '(A)', advance='no') ' '
    
    if (curr_l_faces +5 .lt. 10) then
      write(114, '(I1)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 100) then
      write(114, '(I2)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 1000) then
      write(114, '(I3)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 10000) then
      write(114, '(I4)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 100000) then
      write(114, '(I5)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 1000000) then
      write(114, '(I6)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 10000000) then
      write(114, '(I7)', advance = 'no') curr_l_faces + 5
    end if
    
    write(114, '(A)', advance='no') ' '
    
    if (curr_l_faces +4 .lt. 10) then
      write(114, '(I1)') curr_l_faces +4
    else if (curr_l_faces +4 .lt. 100) then
      write(114, '(I2)') curr_l_faces +4
    else if (curr_l_faces +4 .lt. 1000) then
      write(114, '(I3)') curr_l_faces +4
    else if (curr_l_faces +4 .lt. 10000) then
      write(114, '(I4)') curr_l_faces +4
    else if (curr_l_faces +4 .lt. 100000) then
      write(114, '(I5)') curr_l_faces +4
    else if (curr_l_faces +4 .lt. 1000000) then
      write(114, '(I6)') curr_l_faces +4
    else if (curr_l_faces +4 .lt. 10000000) then
      write(114, '(I7)') curr_l_faces +4
    end if
    
    !!!!! 6th (4-3-7-8)
    write(114, '(I1,A)', advance= 'no') 4, ' '
    
    if (curr_l_faces + 3 .lt. 10) then
      write(114, '(I1)', advance = 'no') curr_l_faces + 3
    else if (curr_l_faces + 3 .lt. 100) then
      write(114, '(I2)', advance = 'no') curr_l_faces + 3
    else if (curr_l_faces + 3 .lt. 1000) then
      write(114, '(I3)', advance = 'no') curr_l_faces + 3
    else if (curr_l_faces + 3 .lt. 10000) then
      write(114, '(I4)', advance = 'no') curr_l_faces + 3
    else if (curr_l_faces + 3 .lt. 100000) then
      write(114, '(I5)', advance = 'no') curr_l_faces + 3
    else if (curr_l_faces + 3 .lt. 1000000) then
      write(114, '(I6)', advance = 'no') curr_l_faces + 3
    else if (curr_l_faces + 3 .lt. 10000000) then
      write(114, '(I7)', advance = 'no') curr_l_faces + 3
    end if
    
    write(114, '(A)', advance='no') ' '
    
    if (curr_l_faces +2 .lt. 10) then
      write(114, '(I1)', advance = 'no') curr_l_faces + 2
    else if (curr_l_faces +2 .lt. 100) then
      write(114, '(I2)', advance = 'no') curr_l_faces + 2
    else if (curr_l_faces +2 .lt. 1000) then
      write(114, '(I3)', advance = 'no') curr_l_faces + 2
    else if (curr_l_faces +2 .lt. 10000) then
      write(114, '(I4)', advance = 'no') curr_l_faces + 2
    else if (curr_l_faces +2 .lt. 100000) then
      write(114, '(I5)', advance = 'no') curr_l_faces + 2
    else if (curr_l_faces +2 .lt. 1000000) then
      write(114, '(I6)', advance = 'no') curr_l_faces + 2
    else if (curr_l_faces +2 .lt. 10000000) then
      write(114, '(I7)', advance = 'no') curr_l_faces + 2
    end if
    
    write(114, '(A)', advance='no') ' '
    
    if (curr_l_faces +6 .lt. 10) then
      write(114, '(I1)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 100) then
      write(114, '(I2)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 1000) then
      write(114, '(I3)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 10000) then
      write(114, '(I4)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 100000) then
      write(114, '(I5)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 1000000) then
      write(114, '(I6)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 10000000) then
      write(114, '(I7)', advance = 'no') curr_l_faces + 6
    end if
    
    write(114, '(A)', advance='no') ' '
    
    if (curr_l_faces +7 .lt. 10) then
      write(114, '(I1)') curr_l_faces +7
    else if (curr_l_faces +7 .lt. 100) then
      write(114, '(I2)') curr_l_faces +7 
    else if (curr_l_faces +7 .lt. 1000) then
      write(114, '(I3)') curr_l_faces +7 
    else if (curr_l_faces +7 .lt. 10000) then
      write(114, '(I4)') curr_l_faces +7
    else if (curr_l_faces +7 .lt. 100000) then
      write(114, '(I5)') curr_l_faces +7
    else if (curr_l_faces +7 .lt. 1000000) then
      write(114, '(I6)') curr_l_faces +7
    else if (curr_l_faces +7 .lt. 10000000) then
      write(114, '(I7)') curr_l_faces +7
    end if
    
    curr_l_faces = curr_l_faces + 8
    
  end do

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!! FIELDS !!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! How many fields there will be?
  n_l_fields = 5                          ! They are: Id, Material, Disp, Veloc, Z (coordination)
                                          ! ...
  
  ! A blank space
  write(114,'(A)') ''
  ! The cells begin
  write(114,'(A)', advance='no') 'CELL_DATA'
  
  ! Writing the number of data by field. It corresponds to the same number of faces
  write(114, '(2(I8,A))') n_l_faces + 12
  write(114, '(A,I4)')  'FIELD FieldData', n_l_fields
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Id
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! Naming the field, the dimension of the data, and number of lines (as before FIELD), and data type
  write(114,'(A,I1,I8,A)') 'Id ', 1, (n_l_faces +12), ' float'
  k = 0
  l_counter = 0
  do i=1, n_particles + n_walls
    l_counter = l_counter + 1
    if (i .le. n_particles) then
      do j=1, sphere_n_faces
        k=k+1
        if(k .le. 9) then
          write(114, '(I6)', advance='no') l_counter
        else 
          write(114, '(I6)') l_counter
          k=0
        end if 
      end do
    else
      ! Six faces for each wall 
      do j=1, 6
        k=k+1
        if(k .le. 9) then
          write(114, '(I6)', advance='no') l_counter
        else 
          write(114, '(I6)') l_counter
          k=0
        end if 
      end do
    end if
  end do
  ! And jump
  write(114, '(A)') ' '

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Material
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(114,'(A,I1,I8,A)') 'Material ', 1, (n_l_faces +12), ' float'
  k = 0
  l_counter = 1
  do i=1, n_particles + n_walls
    if (i .le. n_particles) then
      ! Material number 1
      do j=1, sphere_n_faces
        k=k+1
        if(k .le. 9) then
          write(114, '(I6)', advance='no') l_counter
        else 
          write(114, '(I6)') l_counter
          k=0
        end if 
      end do
    else
      ! Material number 2
      if (i==n_particles+1) then
        l_counter = l_counter +1
      end if
      do j=1, 6
        k=k+1
        if(k .le. 9) then
          write(114, '(I6)', advance='no') l_counter
        else 
          write(114, '(I6)') l_counter
          k=0
        end if 
      end do
    end if
  end do
  
  ! And jump
  write(114, '(A)') ' '
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Displacement
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  curr_l_vector(:) = 0.D0
  
  write(114,'(A,I1,I8,A)') 'Disp ', 3, (n_l_faces +12), ' float'
  k = 0
  do i=1, n_particles + n_walls
    if (i .le. n_particles) then
      curr_l_vector(:) = TAB_SPHERE(i)%center(:) - TAB_SPHERE(i)%center_ref(:)
      do j=1, sphere_n_faces
        k=k+3
        if(k .le. 9) then
          write(114, '(3(F12.9,A))', advance='no') curr_l_vector(1), ' ', curr_l_vector(2), ' ', curr_l_vector(3), ' '
        else 
          write(114, '(3(F12.9,A))') curr_l_vector(1), ' ', curr_l_vector(2), ' ', curr_l_vector(3)
          k=0
        end if 
      end do
    else
      curr_l_vector(:) = TAB_PLAN(i-n_particles)%center(:) - TAB_PLAN(i-n_particles)%center_ref
      if (i==n_particles+1) then
        l_counter = l_counter +1
      end if
      do j=1, 6
        k=k+3
        if(k .le. 9) then
          write(114, '(3(F12.9,A))', advance='no') curr_l_vector(1), ' ', curr_l_vector(2), ' ', curr_l_vector(3), ' '
        else 
          write(114, '(3(F12.9,A))') curr_l_vector(1), ' ', curr_l_vector(2), ' ', curr_l_vector(3)
          k=0
        end if 
      end do
    end if
  end do
  
  ! And jump
  write(114, '(A)') ' '
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Velocity
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  curr_l_vector(:) = 0.D0
  
  write(114,'(A,I1,I8,A)') 'Veloc ', 3, (n_l_faces +12), ' float'
  k = 0
  do i=1, n_particles + n_walls
    if (i .le. n_particles) then
      curr_l_vector(1) = TAB_SPHERE(i)%Vx
      curr_l_vector(2) = TAB_SPHERE(i)%Vy
      curr_l_vector(3) = TAB_SPHERE(i)%Vz
      
      do j=1, sphere_n_faces
        k=k+3
        if(k .le. 9) then
          write(114, '(3(F13.9,A))', advance='no') curr_l_vector(1), ' ', curr_l_vector(2), ' ', curr_l_vector(3), ' '
        else 
          write(114, '(3(F13.9,A))') curr_l_vector(1), ' ', curr_l_vector(2), ' ', curr_l_vector(3)
          k=0
        end if 
      end do
    else
      curr_l_vector(1) = TAB_PLAN(i-n_particles)%Vx
      curr_l_vector(2) = TAB_PLAN(i-n_particles)%Vy
      curr_l_vector(3) = TAB_PLAN(i-n_particles)%Vz
      
      if (i==n_particles+1) then
        l_counter = l_counter +1
      end if
      do j=1, 6
        k=k+3
        if(k .le. 9) then
          write(114, '(3(F13.9,A))', advance='no') curr_l_vector(1), ' ', curr_l_vector(2), ' ', curr_l_vector(3), ' '
        else 
          write(114, '(3(F13.9,A))') curr_l_vector(1), ' ', curr_l_vector(2), ' ', curr_l_vector(3)
          k=0
        end if 
      end do
    end if
  end do
  
  ! And jump
  write(114, '(A)') ' '
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Coordination
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  l_counter = 0
  
  write(114,'(A,I1,I8,A)') 'Z ', 1, (n_l_faces +12), ' float'
  k = 0
  do i=1, n_particles + n_walls
    l_counter = 0
    if (i .le. n_particles) then
      ! Counting the number of contacts
      l_counter = TAB_SPHERE(i)%nctc
        
      do j=1, sphere_n_faces
        k=k+1
        if(k .le. 9) then
          write(114, '(I4,A)', advance='no') l_counter, ' '
        else 
          write(114, '(I4,A)') l_counter
          k=0
        end if
      end do
    else
      do j=1, 6
        k=k+1
        if(k .le. 9) then
          write(114, '(I4,A)', advance='no') 0, ' '
        else 
          write(114, '(I4,A)') 0
          k=0
        end if 
      end do
    end if
  end do
  
  ! And jump
  write(114, '(A)') ' '

  ! Closing the files of the particles
  close(114)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!  FORCES   !!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  vtk_forces = './POSTPRO/VTK/forces_' // vtk_counter // '.vtk'
  
  open(unit=114, file=vtk_forces, status='replace')
  write(114,'(A)') '# vtk DataFile Version 3.0'
  write(114,'(A,I6)') 'FORCES ', compteur_clout
  write(114,'(A)') 'ASCII'
  write(114,'(A)') 'DATASET POLYDATA'
  
  ! Writing the number of vertices for the forces. These are parallelepipeds joining the center of 
  ! particles in contact. As the same time we find the average normal force
  ! Eight vertices for each parallelepiped
  f_counter = 0
  f_counter_periodic = 0
  vtk_ave_force = 0.D0

  disp_wall_times = 3.5
  do i=1, nb_ligneCONTACT

    ! Removing contacts with walls 
    if(TAB_CONTACT(i)%nature /= 'SPSPx') cycle
    
    l_cdt = TAB_CONTACT(i)%icdent
    l_ant = TAB_CONTACT(i)%ianent

    if (l_cdt .gt. n_particles .or. l_ant .gt. n_particles) cycle

    ! Far from the walls 
    if (TAB_SPHERE(l_cdt)%center(3) .lt. pos_z_w1 + ave_rad*disp_wall_times) cycle 
    if (TAB_SPHERE(l_cdt)%center(3) .gt. pos_z_w2 - ave_rad*disp_wall_times) cycle 

    if (TAB_SPHERE(l_ant)%center(3) .lt. pos_z_w1 + ave_rad*disp_wall_times) cycle 
    if (TAB_SPHERE(l_ant)%center(3) .gt. pos_z_w2 - ave_rad*disp_wall_times) cycle 

    if (TAB_CONTACT(i)%rn .le. 0) cycle

    Lik(:) = TAB_SPHERE(l_cdt)%center(:)-TAB_SPHERE(l_ant)%center(:)

    ! Forces in the periodic zone should be counted separately 
    if (sqrt(Lik(1)**2 + Lik(2)**2 + Lik(3)**2) .gt. (TAB_SPHERE(l_cdt)%Rmax + TAB_SPHERE(l_ant)%Rmax)) then 
      f_counter_periodic = f_counter_periodic + 1
    else 
      f_counter = f_counter + 1
    end if 
    vtk_ave_force = vtk_ave_force + TAB_CONTACT(i)%rn
  end do
  
  vtk_ave_force = vtk_ave_force/(f_counter+f_counter_periodic)

  write(v_n_v_forces, '(I8)') (f_counter + f_counter_periodic*2) * 8
  
  vtk_n_v_forces = 'POINTS' // v_n_v_forces // ' float'
  write(114,'(A)') vtk_n_v_forces
  
  ! Force scale parameter 
  ! The scale is here!
  ! 4% of average radius
  l_force_scale = (ave_rad)*0.1

  ! Writing
  k=0
  do i=1, nb_ligneCONTACT
    
    if (TAB_CONTACT(i)%nature /= 'SPSPx') cycle
    ! The particles ids
    l_cdt = TAB_CONTACT(i)%icdent
    l_ant = TAB_CONTACT(i)%ianent

    if (l_cdt .gt. n_particles .or. l_ant .gt. n_particles) cycle

    ! Far from the walls 
    if (TAB_SPHERE(l_cdt)%center(3) .lt. pos_z_w1 + ave_rad*disp_wall_times) cycle 
    if (TAB_SPHERE(l_cdt)%center(3) .gt. pos_z_w2 - ave_rad*disp_wall_times) cycle 

    if (TAB_SPHERE(l_ant)%center(3) .lt. pos_z_w1 + ave_rad*disp_wall_times) cycle 
    if (TAB_SPHERE(l_ant)%center(3) .gt. pos_z_w2 - ave_rad*disp_wall_times) cycle 

    if (TAB_CONTACT(i)%rn .le. 0) cycle

    Lik(:) = TAB_SPHERE(l_cdt)%center(:)-TAB_SPHERE(l_ant)%center(:)

    ! Being careful with particles in the periodic borders 
    ! Periodic case 
    if (sqrt(Lik(1)**2 + Lik(2)**2 + Lik(3)**2) .gt. (TAB_SPHERE(l_cdt)%Rmax + TAB_SPHERE(l_ant)%Rmax)) then 
      ! The center of the candidate and computing the center of the Antagonist
      vtk_cd_center(:) = TAB_SPHERE(l_cdt)%center(:)

      ! The normal and tangential vectors
      v_l_normal(:) = TAB_CONTACT(i)%n(:)
      v_l_t(:) = TAB_CONTACT(i)%t(:)
      v_l_s(:) = TAB_CONTACT(i)%s(:)

      ! The center of the Antagonist
      vtk_an_center(:) = TAB_SPHERE(l_cdt)%center(:) - (TAB_SPHERE(l_cdt)%Rmax+TAB_SPHERE(l_ant)%Rmax)*v_l_normal(:)

      ! The forces
      l_rn = TAB_CONTACT(i)%rn
      l_rt = TAB_CONTACT(i)%rt
      l_rs = TAB_CONTACT(i)%rs
    
      l_force = (l_rn/vtk_ave_force)*l_force_scale
    
      ! Write! ... 4 vertices for each line
    
      ! Candidate --- First vertex. 
      write(114, '(3(F15.7,A))', advance='no') vtk_cd_center(1)+l_force*v_l_t(1)+ l_force*v_l_s(1), ' ', &
                                               vtk_cd_center(2)+l_force*v_l_t(2)+ l_force*v_l_s(2), ' ', &
                                               vtk_cd_center(3)+l_force*v_l_t(3)+ l_force*v_l_s(3), ' '
      ! Candidate --- Second vertex. 
      write(114, '(3(F15.7,A))', advance='no') vtk_cd_center(1)+l_force*v_l_t(1)- l_force*v_l_s(1), ' ', &
                                               vtk_cd_center(2)+l_force*v_l_t(2)- l_force*v_l_s(2), ' ', &
                                               vtk_cd_center(3)+l_force*v_l_t(3)- l_force*v_l_s(3), ' '
      ! Candidate --- Third vertex. 
      write(114, '(3(F15.7,A))', advance='no') vtk_cd_center(1)-l_force*v_l_t(1)- l_force*v_l_s(1), ' ', &
                                               vtk_cd_center(2)-l_force*v_l_t(2)- l_force*v_l_s(2), ' ', &
                                               vtk_cd_center(3)-l_force*v_l_t(3)- l_force*v_l_s(3), ' '
      ! Candidate --- 4th vertex. 
      write(114, '(3(F15.7,A))') vtk_cd_center(1)-l_force*v_l_t(1)+ l_force*v_l_s(1), ' ', &
                                 vtk_cd_center(2)-l_force*v_l_t(2)+ l_force*v_l_s(2), ' ', &
                                 vtk_cd_center(3)-l_force*v_l_t(3)+ l_force*v_l_s(3), ' ' 
    
      ! Antagonist --- First vertex. 
      write(114, '(3(F15.7,A))', advance='no') vtk_an_center(1)+l_force*v_l_t(1)+ l_force*v_l_s(1), ' ', &
                                               vtk_an_center(2)+l_force*v_l_t(2)+ l_force*v_l_s(2), ' ', &
                                               vtk_an_center(3)+l_force*v_l_t(3)+ l_force*v_l_s(3), ' '
      ! Antagonist --- Second vertex.          
      write(114, '(3(F15.7,A))', advance='no') vtk_an_center(1)+l_force*v_l_t(1)- l_force*v_l_s(1), ' ', &
                                               vtk_an_center(2)+l_force*v_l_t(2)- l_force*v_l_s(2), ' ', &
                                               vtk_an_center(3)+l_force*v_l_t(3)- l_force*v_l_s(3), ' '
      ! Antagonist --- Third vertex.           
      write(114, '(3(F15.7,A))', advance='no') vtk_an_center(1)-l_force*v_l_t(1)- l_force*v_l_s(1), ' ', &
                                               vtk_an_center(2)-l_force*v_l_t(2)- l_force*v_l_s(2), ' ', &
                                               vtk_an_center(3)-l_force*v_l_t(3)- l_force*v_l_s(3), ' '
      ! Antagonist --- 4th vertex.             
      write(114, '(3(F15.7,A))') vtk_an_center(1)-l_force*v_l_t(1)+ l_force*v_l_s(1), ' ', &
                                 vtk_an_center(2)-l_force*v_l_t(2)+ l_force*v_l_s(2), ' ', &
                                 vtk_an_center(3)-l_force*v_l_t(3)+ l_force*v_l_s(3), ' ' 

      

      ! The center of the antagonist and computing the center of the candidate
      vtk_an_center(:) = TAB_SPHERE(l_ant)%center(:)

      ! The normal and tangential vectors
      v_l_normal(:) = TAB_CONTACT(i)%n(:)
      v_l_t(:) = TAB_CONTACT(i)%t(:)
      v_l_s(:) = TAB_CONTACT(i)%s(:)

      ! The center of the Antagonist
      vtk_cd_center(:) = TAB_SPHERE(l_ant)%center(:) + (TAB_SPHERE(l_cdt)%Rmax+TAB_SPHERE(l_ant)%Rmax)*v_l_normal(:)

      ! The forces
      l_rn = TAB_CONTACT(i)%rn
      l_rt = TAB_CONTACT(i)%rt
      l_rs = TAB_CONTACT(i)%rs
    
      l_force = (l_rn/vtk_ave_force)*l_force_scale
    
      ! Write! ... 4 vertices for each line
    
      ! Candidate --- First vertex. 
      write(114, '(3(F15.7,A))', advance='no') vtk_cd_center(1)+l_force*v_l_t(1)+ l_force*v_l_s(1), ' ', &
                                               vtk_cd_center(2)+l_force*v_l_t(2)+ l_force*v_l_s(2), ' ', &
                                               vtk_cd_center(3)+l_force*v_l_t(3)+ l_force*v_l_s(3), ' '
      ! Candidate --- Second vertex. 
      write(114, '(3(F15.7,A))', advance='no') vtk_cd_center(1)+l_force*v_l_t(1)- l_force*v_l_s(1), ' ', &
                                               vtk_cd_center(2)+l_force*v_l_t(2)- l_force*v_l_s(2), ' ', &
                                               vtk_cd_center(3)+l_force*v_l_t(3)- l_force*v_l_s(3), ' '
      ! Candidate --- Third vertex. 
      write(114, '(3(F15.7,A))', advance='no') vtk_cd_center(1)-l_force*v_l_t(1)- l_force*v_l_s(1), ' ', &
                                               vtk_cd_center(2)-l_force*v_l_t(2)- l_force*v_l_s(2), ' ', &
                                               vtk_cd_center(3)-l_force*v_l_t(3)- l_force*v_l_s(3), ' '
      ! Candidate --- 4th vertex. 
      write(114, '(3(F15.7,A))') vtk_cd_center(1)-l_force*v_l_t(1)+ l_force*v_l_s(1), ' ', &
                                 vtk_cd_center(2)-l_force*v_l_t(2)+ l_force*v_l_s(2), ' ', &
                                 vtk_cd_center(3)-l_force*v_l_t(3)+ l_force*v_l_s(3), ' ' 
    
      ! Antagonist --- First vertex. 
      write(114, '(3(F15.7,A))', advance='no') vtk_an_center(1)+l_force*v_l_t(1)+ l_force*v_l_s(1), ' ', &
                                               vtk_an_center(2)+l_force*v_l_t(2)+ l_force*v_l_s(2), ' ', &
                                               vtk_an_center(3)+l_force*v_l_t(3)+ l_force*v_l_s(3), ' '
      ! Antagonist --- Second vertex.          
      write(114, '(3(F15.7,A))', advance='no') vtk_an_center(1)+l_force*v_l_t(1)- l_force*v_l_s(1), ' ', &
                                               vtk_an_center(2)+l_force*v_l_t(2)- l_force*v_l_s(2), ' ', &
                                               vtk_an_center(3)+l_force*v_l_t(3)- l_force*v_l_s(3), ' '
      ! Antagonist --- Third vertex.           
      write(114, '(3(F15.7,A))', advance='no') vtk_an_center(1)-l_force*v_l_t(1)- l_force*v_l_s(1), ' ', &
                                               vtk_an_center(2)-l_force*v_l_t(2)- l_force*v_l_s(2), ' ', &
                                               vtk_an_center(3)-l_force*v_l_t(3)- l_force*v_l_s(3), ' '
      ! Antagonist --- 4th vertex.             
      write(114, '(3(F15.7,A))') vtk_an_center(1)-l_force*v_l_t(1)+ l_force*v_l_s(1), ' ', &
                                 vtk_an_center(2)-l_force*v_l_t(2)+ l_force*v_l_s(2), ' ', &
                                 vtk_an_center(3)-l_force*v_l_t(3)+ l_force*v_l_s(3), ' ' 

    else 

      ! The center of each particle
      vtk_cd_center(:) = TAB_SPHERE(l_cdt)%center(:)
      vtk_an_center(:) = TAB_SPHERE(l_ant)%center(:)

      ! The normal and tangential vectors
      v_l_normal(:) = TAB_CONTACT(i)%n(:)
      v_l_t(:) = TAB_CONTACT(i)%t(:)
      v_l_s(:) = TAB_CONTACT(i)%s(:)

      ! The forces
      l_rn = TAB_CONTACT(i)%rn
      l_rt = TAB_CONTACT(i)%rt
      l_rs = TAB_CONTACT(i)%rs
    
      l_force = (l_rn/vtk_ave_force)*l_force_scale
    
      ! Write! ... 4 vertices for each line
    
      ! Candidate --- First vertex. 
      write(114, '(3(F15.7,A))', advance='no') vtk_cd_center(1)+l_force*v_l_t(1)+ l_force*v_l_s(1), ' ', &
                                               vtk_cd_center(2)+l_force*v_l_t(2)+ l_force*v_l_s(2), ' ', &
                                               vtk_cd_center(3)+l_force*v_l_t(3)+ l_force*v_l_s(3), ' '
      ! Candidate --- Second vertex. 
      write(114, '(3(F15.7,A))', advance='no') vtk_cd_center(1)+l_force*v_l_t(1)- l_force*v_l_s(1), ' ', &
                                               vtk_cd_center(2)+l_force*v_l_t(2)- l_force*v_l_s(2), ' ', &
                                               vtk_cd_center(3)+l_force*v_l_t(3)- l_force*v_l_s(3), ' '
      ! Candidate --- Third vertex. 
      write(114, '(3(F15.7,A))', advance='no') vtk_cd_center(1)-l_force*v_l_t(1)- l_force*v_l_s(1), ' ', &
                                               vtk_cd_center(2)-l_force*v_l_t(2)- l_force*v_l_s(2), ' ', &
                                               vtk_cd_center(3)-l_force*v_l_t(3)- l_force*v_l_s(3), ' '
      ! Candidate --- 4th vertex. 
      write(114, '(3(F15.7,A))') vtk_cd_center(1)-l_force*v_l_t(1)+ l_force*v_l_s(1), ' ', &
                                 vtk_cd_center(2)-l_force*v_l_t(2)+ l_force*v_l_s(2), ' ', &
                                 vtk_cd_center(3)-l_force*v_l_t(3)+ l_force*v_l_s(3), ' ' 
    
      ! Antagonist --- First vertex. 
      write(114, '(3(F15.7,A))', advance='no') vtk_an_center(1)+l_force*v_l_t(1)+ l_force*v_l_s(1), ' ', &
                                               vtk_an_center(2)+l_force*v_l_t(2)+ l_force*v_l_s(2), ' ', &
                                               vtk_an_center(3)+l_force*v_l_t(3)+ l_force*v_l_s(3), ' '
      ! Antagonist --- Second vertex.          
      write(114, '(3(F15.7,A))', advance='no') vtk_an_center(1)+l_force*v_l_t(1)- l_force*v_l_s(1), ' ', &
                                               vtk_an_center(2)+l_force*v_l_t(2)- l_force*v_l_s(2), ' ', &
                                               vtk_an_center(3)+l_force*v_l_t(3)- l_force*v_l_s(3), ' '
      ! Antagonist --- Third vertex.           
      write(114, '(3(F15.7,A))', advance='no') vtk_an_center(1)-l_force*v_l_t(1)- l_force*v_l_s(1), ' ', &
                                               vtk_an_center(2)-l_force*v_l_t(2)- l_force*v_l_s(2), ' ', &
                                               vtk_an_center(3)-l_force*v_l_t(3)- l_force*v_l_s(3), ' '
      ! Antagonist --- 4th vertex.             
      write(114, '(3(F15.7,A))') vtk_an_center(1)-l_force*v_l_t(1)+ l_force*v_l_s(1), ' ', &
                                 vtk_an_center(2)-l_force*v_l_t(2)+ l_force*v_l_s(2), ' ', &
                                 vtk_an_center(3)-l_force*v_l_t(3)+ l_force*v_l_s(3), ' ' 
    end if 
  end do

! Writing the conectivity
  write(114,'(A)', advance='no') 'POLYGONS '
  ! 6 faces for each parallepiped
  write(114, '(2(I8,A))') (f_counter+f_counter_periodic*2)*6, ' ' , (f_counter+f_counter_periodic*2)*6*5
  
  curr_l_faces = 0
  
  do i=1, (f_counter+f_counter_periodic*2)
    
    !!!!! 1st (1-2-3-4)
    write(114, '(I1,A)', advance= 'no') 4, ' '
    
    if (curr_l_faces .lt. 10) then
      write(114, '(I1)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 100) then
      write(114, '(I2)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 1000) then
      write(114, '(I3)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 10000) then
      write(114, '(I4)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 100000) then
      write(114, '(I5)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 1000000) then
      write(114, '(I6)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 10000000) then
      write(114, '(I7)', advance = 'no') curr_l_faces
    end if
    
    write(114, '(A)', advance='no') ' '
    
    if (curr_l_faces +1 .lt. 10) then
      write(114, '(I1)', advance = 'no') curr_l_faces +1 
    else if (curr_l_faces +1 .lt. 100) then
      write(114, '(I2)', advance = 'no') curr_l_faces +1 
    else if (curr_l_faces +1 .lt. 1000) then
      write(114, '(I3)', advance = 'no') curr_l_faces +1 
    else if (curr_l_faces +1 .lt. 10000) then
      write(114, '(I4)', advance = 'no') curr_l_faces +1
    else if (curr_l_faces +1 .lt. 100000) then
      write(114, '(I5)', advance = 'no') curr_l_faces +1
    else if (curr_l_faces +1 .lt. 1000000) then
      write(114, '(I6)', advance = 'no') curr_l_faces +1 
    else if (curr_l_faces +1 .lt. 10000000) then
      write(114, '(I7)', advance = 'no') curr_l_faces +1 
    end if
    
    write(114, '(A)', advance='no') ' '
    
    if (curr_l_faces +2 .lt. 10) then
      write(114, '(I1)', advance = 'no') curr_l_faces +2 
    else if (curr_l_faces +2 .lt. 100) then
      write(114, '(I2)', advance = 'no') curr_l_faces +2 
    else if (curr_l_faces +2 .lt. 1000) then
      write(114, '(I3)', advance = 'no') curr_l_faces +2 
    else if (curr_l_faces +2 .lt. 10000) then
      write(114, '(I4)', advance = 'no') curr_l_faces +2
    else if (curr_l_faces +2 .lt. 100000) then
      write(114, '(I5)', advance = 'no') curr_l_faces +2
    else if (curr_l_faces +2 .lt. 1000000) then
      write(114, '(I6)', advance = 'no') curr_l_faces +2 
    else if (curr_l_faces +2 .lt. 10000000) then
      write(114, '(I7)', advance = 'no') curr_l_faces +2 
    end if
    
    write(114, '(A)', advance='no') ' '
    
    if (curr_l_faces +3 .lt. 10) then
      write(114, '(I1)') curr_l_faces +3
    else if (curr_l_faces +3 .lt. 100) then
      write(114, '(I2)') curr_l_faces +3 
    else if (curr_l_faces +3 .lt. 1000) then
      write(114, '(I3)') curr_l_faces +3 
    else if (curr_l_faces +3 .lt. 10000) then
      write(114, '(I4)') curr_l_faces +3
    else if (curr_l_faces +3 .lt. 100000) then
      write(114, '(I5)') curr_l_faces +3
    else if (curr_l_faces +3 .lt. 1000000) then
      write(114, '(I6)') curr_l_faces +3
    else if (curr_l_faces +3 .lt. 10000000) then
      write(114, '(I7)') curr_l_faces +3
    end if
    
    !!!!! 2nd (5-6-7-8)
    write(114, '(I1,A)', advance= 'no') 4, ' '
    
    if (curr_l_faces + 4 .lt. 10) then
      write(114, '(I1)', advance = 'no') curr_l_faces + 4
    else if (curr_l_faces + 4 .lt. 100) then
      write(114, '(I2)', advance = 'no') curr_l_faces + 4
    else if (curr_l_faces + 4 .lt. 1000) then
      write(114, '(I3)', advance = 'no') curr_l_faces + 4
    else if (curr_l_faces + 4 .lt. 10000) then
      write(114, '(I4)', advance = 'no') curr_l_faces + 4 
    else if (curr_l_faces + 4 .lt. 100000) then
      write(114, '(I5)', advance = 'no') curr_l_faces + 4
    else if (curr_l_faces + 4 .lt. 1000000) then
      write(114, '(I6)', advance = 'no') curr_l_faces + 4
    else if (curr_l_faces + 4 .lt. 10000000) then
      write(114, '(I7)', advance = 'no') curr_l_faces + 4
    end if
    
    write(114, '(A)', advance='no') ' '
    
    if (curr_l_faces +5 .lt. 10) then
      write(114, '(I1)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 100) then
      write(114, '(I2)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 1000) then
      write(114, '(I3)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 10000) then
      write(114, '(I4)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 100000) then
      write(114, '(I5)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 1000000) then
      write(114, '(I6)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 10000000) then
      write(114, '(I7)', advance = 'no') curr_l_faces + 5
    end if
    
    write(114, '(A)', advance='no') ' '
    
    if (curr_l_faces +6 .lt. 10) then
      write(114, '(I1)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 100) then
      write(114, '(I2)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 1000) then
      write(114, '(I3)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 10000) then
      write(114, '(I4)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 100000) then
      write(114, '(I5)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 1000000) then
      write(114, '(I6)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 10000000) then
      write(114, '(I7)', advance = 'no') curr_l_faces + 6
    end if
    
    write(114, '(A)', advance='no') ' '
    
    if (curr_l_faces +7 .lt. 10) then
      write(114, '(I1)') curr_l_faces +7
    else if (curr_l_faces +7 .lt. 100) then
      write(114, '(I2)') curr_l_faces +7 
    else if (curr_l_faces +7 .lt. 1000) then
      write(114, '(I3)') curr_l_faces +7 
    else if (curr_l_faces +7 .lt. 10000) then
      write(114, '(I4)') curr_l_faces +7
    else if (curr_l_faces +7 .lt. 100000) then
      write(114, '(I5)') curr_l_faces +7
    else if (curr_l_faces +7 .lt. 1000000) then
      write(114, '(I6)') curr_l_faces +7
    else if (curr_l_faces +7 .lt. 10000000) then
      write(114, '(I7)') curr_l_faces +7
    end if
    
    !!!!! 3nd (1-5-8-4)
    write(114, '(I1,A)', advance= 'no') 4, ' '
    
    if (curr_l_faces .lt. 10) then
      write(114, '(I1)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 100) then
      write(114, '(I2)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 1000) then
      write(114, '(I3)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 10000) then
      write(114, '(I4)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 100000) then
      write(114, '(I5)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 1000000) then
      write(114, '(I6)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 10000000) then
      write(114, '(I7)', advance = 'no') curr_l_faces
    end if
    
    write(114, '(A)', advance='no') ' '
    
    if (curr_l_faces +4 .lt. 10) then
      write(114, '(I1)', advance = 'no') curr_l_faces + 4
    else if (curr_l_faces +4 .lt. 100) then
      write(114, '(I2)', advance = 'no') curr_l_faces + 4
    else if (curr_l_faces +4 .lt. 1000) then
      write(114, '(I3)', advance = 'no') curr_l_faces + 4
    else if (curr_l_faces +4 .lt. 10000) then
      write(114, '(I4)', advance = 'no') curr_l_faces + 4
    else if (curr_l_faces +4 .lt. 100000) then
      write(114, '(I5)', advance = 'no') curr_l_faces + 4
    else if (curr_l_faces +4 .lt. 1000000) then
      write(114, '(I6)', advance = 'no') curr_l_faces + 4
    else if (curr_l_faces +4 .lt. 10000000) then
      write(114, '(I7)', advance = 'no') curr_l_faces + 4
    end if
    
    write(114, '(A)', advance='no') ' '
    
    if (curr_l_faces +7 .lt. 10) then
      write(114, '(I1)', advance = 'no') curr_l_faces + 7
    else if (curr_l_faces +7 .lt. 100) then
      write(114, '(I2)', advance = 'no') curr_l_faces + 7
    else if (curr_l_faces +7 .lt. 1000) then
      write(114, '(I3)', advance = 'no') curr_l_faces + 7
    else if (curr_l_faces +7 .lt. 10000) then
      write(114, '(I4)', advance = 'no') curr_l_faces + 7
    else if (curr_l_faces +7 .lt. 100000) then
      write(114, '(I5)', advance = 'no') curr_l_faces + 7
    else if (curr_l_faces +7 .lt. 1000000) then
      write(114, '(I6)', advance = 'no') curr_l_faces + 7
    else if (curr_l_faces +7 .lt. 10000000) then
      write(114, '(I7)', advance = 'no') curr_l_faces + 7
    end if
    
    write(114, '(A)', advance='no') ' '
    
    if (curr_l_faces +3 .lt. 10) then
      write(114, '(I1)') curr_l_faces +3
    else if (curr_l_faces +3 .lt. 100) then
      write(114, '(I2)') curr_l_faces +3 
    else if (curr_l_faces +3 .lt. 1000) then
      write(114, '(I3)') curr_l_faces +3 
    else if (curr_l_faces +3 .lt. 10000) then
      write(114, '(I4)') curr_l_faces +3
    else if (curr_l_faces +3 .lt. 100000) then
      write(114, '(I5)') curr_l_faces +3
    else if (curr_l_faces +3 .lt. 1000000) then
      write(114, '(I6)') curr_l_faces +3
    else if (curr_l_faces +3 .lt. 10000000) then
      write(114, '(I7)') curr_l_faces +3
    end if
    
    
    !!!!! 4th (2-6-7-3)
    write(114, '(I1,A)', advance= 'no') 4, ' '
    
    if (curr_l_faces + 1 .lt. 10) then
      write(114, '(I1)', advance = 'no') curr_l_faces + 1
    else if (curr_l_faces + 1 .lt. 100) then
      write(114, '(I2)', advance = 'no') curr_l_faces + 1
    else if (curr_l_faces + 1 .lt. 1000) then
      write(114, '(I3)', advance = 'no') curr_l_faces + 1
    else if (curr_l_faces + 1 .lt. 10000) then
      write(114, '(I4)', advance = 'no') curr_l_faces + 1
    else if (curr_l_faces + 1 .lt. 100000) then
      write(114, '(I5)', advance = 'no') curr_l_faces + 1
    else if (curr_l_faces + 1 .lt. 1000000) then
      write(114, '(I6)', advance = 'no') curr_l_faces + 1
    else if (curr_l_faces + 1 .lt. 10000000) then
      write(114, '(I7)', advance = 'no') curr_l_faces + 1
    end if
    
    write(114, '(A)', advance='no') ' '
    
    if (curr_l_faces +5 .lt. 10) then
      write(114, '(I1)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 100) then
      write(114, '(I2)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 1000) then
      write(114, '(I3)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 10000) then
      write(114, '(I4)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 100000) then
      write(114, '(I5)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 1000000) then
      write(114, '(I6)', advance = 'no') curr_l_faces + 5 
    else if (curr_l_faces +5 .lt. 10000000) then
      write(114, '(I7)', advance = 'no') curr_l_faces + 5 
    end if
    
    write(114, '(A)', advance='no') ' '
    
    if (curr_l_faces +6 .lt. 10) then
      write(114, '(I1)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 100) then
      write(114, '(I2)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 1000) then
      write(114, '(I3)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 10000) then
      write(114, '(I4)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 100000) then
      write(114, '(I5)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 1000000) then
      write(114, '(I6)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 10000000) then
      write(114, '(I7)', advance = 'no') curr_l_faces + 6
    end if
    
    write(114, '(A)', advance='no') ' '
    
    if (curr_l_faces +2 .lt. 10) then
      write(114, '(I1)') curr_l_faces +2
    else if (curr_l_faces +2 .lt. 100) then
      write(114, '(I2)') curr_l_faces +2
    else if (curr_l_faces +2 .lt. 1000) then
      write(114, '(I3)') curr_l_faces +2
    else if (curr_l_faces +2 .lt. 10000) then
      write(114, '(I4)') curr_l_faces +2
    else if (curr_l_faces +2 .lt. 100000) then
      write(114, '(I5)') curr_l_faces +2
    else if (curr_l_faces +2 .lt. 1000000) then
      write(114, '(I6)') curr_l_faces +2
    else if (curr_l_faces +2 .lt. 10000000) then
      write(114, '(I7)') curr_l_faces +2
    end if
    
    !!!!! 5th (1-2-6-5)
    write(114, '(I1,A)', advance= 'no') 4, ' '
    
    if (curr_l_faces .lt. 10) then
      write(114, '(I1)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 100) then
      write(114, '(I2)', advance = 'no') curr_l_faces
    else if (curr_l_faces  .lt. 1000) then
      write(114, '(I3)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 10000) then
      write(114, '(I4)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 100000) then
      write(114, '(I5)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 1000000) then
      write(114, '(I6)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 10000000) then
      write(114, '(I7)', advance = 'no') curr_l_faces
    end if
    
    write(114, '(A)', advance='no') ' '
    
    if (curr_l_faces +1 .lt. 10) then
      write(114, '(I1)', advance = 'no') curr_l_faces + 1
    else if (curr_l_faces +1 .lt. 100) then
      write(114, '(I2)', advance = 'no') curr_l_faces + 1
    else if (curr_l_faces +1 .lt. 1000) then
      write(114, '(I3)', advance = 'no') curr_l_faces + 1
    else if (curr_l_faces +1 .lt. 10000) then
      write(114, '(I4)', advance = 'no') curr_l_faces + 1
    else if (curr_l_faces +1 .lt. 100000) then
      write(114, '(I5)', advance = 'no') curr_l_faces + 1
    else if (curr_l_faces +1 .lt. 1000000) then
      write(114, '(I6)', advance = 'no') curr_l_faces + 1
    else if (curr_l_faces +1 .lt. 10000000) then
      write(114, '(I7)', advance = 'no') curr_l_faces + 1
    end if
    
    write(114, '(A)', advance='no') ' '
    
    if (curr_l_faces +5 .lt. 10) then
      write(114, '(I1)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 100) then
      write(114, '(I2)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 1000) then
      write(114, '(I3)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 10000) then
      write(114, '(I4)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 100000) then
      write(114, '(I5)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 1000000) then
      write(114, '(I6)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 10000000) then
      write(114, '(I7)', advance = 'no') curr_l_faces + 5
    end if
    
    write(114, '(A)', advance='no') ' '
    
    if (curr_l_faces +4 .lt. 10) then
      write(114, '(I1)') curr_l_faces +4
    else if (curr_l_faces +4 .lt. 100) then
      write(114, '(I2)') curr_l_faces +4
    else if (curr_l_faces +4 .lt. 1000) then
      write(114, '(I3)') curr_l_faces +4
    else if (curr_l_faces +4 .lt. 10000) then
      write(114, '(I4)') curr_l_faces +4
    else if (curr_l_faces +4 .lt. 100000) then
      write(114, '(I5)') curr_l_faces +4
    else if (curr_l_faces +4 .lt. 1000000) then
      write(114, '(I6)') curr_l_faces +4
    else if (curr_l_faces +4 .lt. 10000000) then
      write(114, '(I7)') curr_l_faces +4
    end if
    
    
    !!!!! 6th (4-3-7-8)
    write(114, '(I1,A)', advance= 'no') 4, ' '
    
    if (curr_l_faces + 3 .lt. 10) then
      write(114, '(I1)', advance = 'no') curr_l_faces + 3
    else if (curr_l_faces + 3 .lt. 100) then
      write(114, '(I2)', advance = 'no') curr_l_faces + 3
    else if (curr_l_faces + 3 .lt. 1000) then
      write(114, '(I3)', advance = 'no') curr_l_faces + 3
    else if (curr_l_faces + 3 .lt. 10000) then
      write(114, '(I4)', advance = 'no') curr_l_faces + 3
    else if (curr_l_faces + 3 .lt. 100000) then
      write(114, '(I5)', advance = 'no') curr_l_faces + 3
    else if (curr_l_faces + 3 .lt. 1000000) then
      write(114, '(I6)', advance = 'no') curr_l_faces + 3
    else if (curr_l_faces + 3 .lt. 10000000) then
      write(114, '(I7)', advance = 'no') curr_l_faces + 3
    end if
    
    write(114, '(A)', advance='no') ' '
    
    if (curr_l_faces +2 .lt. 10) then
      write(114, '(I1)', advance = 'no') curr_l_faces + 2
    else if (curr_l_faces +2 .lt. 100) then
      write(114, '(I2)', advance = 'no') curr_l_faces + 2
    else if (curr_l_faces +2 .lt. 1000) then
      write(114, '(I3)', advance = 'no') curr_l_faces + 2
    else if (curr_l_faces +2 .lt. 10000) then
      write(114, '(I4)', advance = 'no') curr_l_faces + 2
    else if (curr_l_faces +2 .lt. 100000) then
      write(114, '(I5)', advance = 'no') curr_l_faces + 2
    else if (curr_l_faces +2 .lt. 1000000) then
      write(114, '(I6)', advance = 'no') curr_l_faces + 2
    else if (curr_l_faces +2 .lt. 10000000) then
      write(114, '(I7)', advance = 'no') curr_l_faces + 2
    end if
    
    write(114, '(A)', advance='no') ' '
    
    if (curr_l_faces +6 .lt. 10) then
      write(114, '(I1)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 100) then
      write(114, '(I2)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 1000) then
      write(114, '(I3)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 10000) then
      write(114, '(I4)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 100000) then
      write(114, '(I5)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 1000000) then
      write(114, '(I6)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 10000000) then
      write(114, '(I7)', advance = 'no') curr_l_faces + 6
    end if
    
    write(114, '(A)', advance='no') ' '
    
    if (curr_l_faces +7 .lt. 10) then
      write(114, '(I1)') curr_l_faces +7
    else if (curr_l_faces +7 .lt. 100) then
      write(114, '(I2)') curr_l_faces +7 
    else if (curr_l_faces +7 .lt. 1000) then
      write(114, '(I3)') curr_l_faces +7 
    else if (curr_l_faces +7 .lt. 10000) then
      write(114, '(I4)') curr_l_faces +7
    else if (curr_l_faces +7 .lt. 100000) then
      write(114, '(I5)') curr_l_faces +7
    else if (curr_l_faces +7 .lt. 1000000) then
      write(114, '(I6)') curr_l_faces +7
    else if (curr_l_faces +7 .lt. 10000000) then
      write(114, '(I7)') curr_l_faces +7
    end if
    
    curr_l_faces = curr_l_faces + 8
  end do

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!! FORCES FIELDS !!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! How many fields there will be?
  n_l_fields = 3                          ! They are: RN, RT, Friction Mobilization Fric_M
                                          ! ...
  
  ! A blank space
  write(114,'(A)') ''
  ! The cells begin
  write(114,'(A)', advance='no') 'CELL_DATA'
  
  ! Writing the number of data by field. It corresponds to the same number of faces
  write(114, '(2(I8,A))') (f_counter+f_counter_periodic*2)*6
  write(114, '(A,I4)')  'FIELD FieldData', n_l_fields
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Normal Force
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! Naming the field, the dimension of the data, and number of lines (as before FIELD), and data type
  write(114,'(A,I1,I8,A)') 'RN ', 1, (f_counter+f_counter_periodic*2)*6, ' float'
  k = 0
  do i=1, nb_ligneCONTACT
    if (TAB_CONTACT(i)%nature /= 'SPSPx') cycle
    l_cdt = TAB_CONTACT(i)%icdent
    l_ant = TAB_CONTACT(i)%ianent

    if (l_cdt .gt. n_particles .or. l_ant .gt. n_particles) cycle

    ! Far from the walls 
    if (TAB_SPHERE(l_cdt)%center(3) .lt. pos_z_w1 + ave_rad*disp_wall_times) cycle 
    if (TAB_SPHERE(l_cdt)%center(3) .gt. pos_z_w2 - ave_rad*disp_wall_times) cycle 

    if (TAB_SPHERE(l_ant)%center(3) .lt. pos_z_w1 + ave_rad*disp_wall_times) cycle 
    if (TAB_SPHERE(l_ant)%center(3) .gt. pos_z_w2 - ave_rad*disp_wall_times) cycle 

    Lik(:) = TAB_SPHERE(l_cdt)%center(:)-TAB_SPHERE(l_ant)%center(:)

    if (TAB_CONTACT(i)%rn .le. 0) cycle
    
    ! Being careful with the periodic case
    if (sqrt(Lik(1)**2 + Lik(2)**2 + Lik(3)**2) .gt. (TAB_SPHERE(l_cdt)%Rmax + TAB_SPHERE(l_ant)%Rmax)) then 
      do j=1, 12
        k=k+1
        if(k .le. 9) then
          write(114, '(F15.7)', advance='no') TAB_CONTACT(i)%rn
        else 
          write(114, '(F15.7)') TAB_CONTACT(i)%rn
          k=0
        end if 
      end do
    else 
      do j=1, 6
        k=k+1
        if(k .le. 9) then
          write(114, '(F15.7)', advance='no') TAB_CONTACT(i)%rn
        else 
          write(114, '(F15.7)') TAB_CONTACT(i)%rn
          k=0
        end if 
      end do
    end if 
  end do
  ! And jump
  write(114, '(A)') ' '
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Tangential force (Rs+Rt)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! Naming the field, the dimension of the data, and number of lines (as before FIELD), and data type
  write(114,'(A,I1,I8,A)') 'RT ', 1, (f_counter+f_counter_periodic*2)*6, ' float'
  k = 0
  do i=1, nb_ligneCONTACT
    if (TAB_CONTACT(i)%nature /= 'SPSPx') cycle
    l_cdt = TAB_CONTACT(i)%icdent
    l_ant = TAB_CONTACT(i)%ianent

    if (l_cdt .gt. n_particles .or. l_ant .gt. n_particles) cycle

    ! Far from the walls 
    if (TAB_SPHERE(l_cdt)%center(3) .lt. pos_z_w1 + ave_rad*disp_wall_times) cycle 
    if (TAB_SPHERE(l_cdt)%center(3) .gt. pos_z_w2 - ave_rad*disp_wall_times) cycle 

    if (TAB_SPHERE(l_ant)%center(3) .lt. pos_z_w1 + ave_rad*disp_wall_times) cycle 
    if (TAB_SPHERE(l_ant)%center(3) .gt. pos_z_w2 - ave_rad*disp_wall_times) cycle 

    Lik(:) = TAB_SPHERE(l_cdt)%center(:)-TAB_SPHERE(l_ant)%center(:)

    if (TAB_CONTACT(i)%rn .le. 0) cycle
    
    ! Being careful with the periodic case
    if (sqrt(Lik(1)**2 + Lik(2)**2 + Lik(3)**2) .gt. (TAB_SPHERE(l_cdt)%Rmax + TAB_SPHERE(l_ant)%Rmax)) then 
      do j=1, 12
        k=k+1
        if(k .le. 9) then
          write(114, '(F15.7)', advance='no') sqrt(TAB_CONTACT(i)%rt**2 + TAB_CONTACT(i)%rs**2)
        else 
          write(114, '(F15.7)') sqrt(TAB_CONTACT(i)%rt**2 + TAB_CONTACT(i)%rs**2)
          k=0
        end if 
      end do
    else 
      do j=1, 6
        k=k+1
        if(k .le. 9) then
          write(114, '(F15.7)', advance='no') sqrt(TAB_CONTACT(i)%rt**2 + TAB_CONTACT(i)%rs**2)
        else 
          write(114, '(F15.7)') sqrt(TAB_CONTACT(i)%rt**2 + TAB_CONTACT(i)%rs**2)
          k=0
        end if 
      end do
    end if
  end do
  ! And jump
  write(114, '(A)') ' '


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Friction mobilization (Rs+Rt)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! Naming the field, the dimension of the data, and number of lines (as before FIELD), and data type
  write(114,'(A,I1,I8,A)') 'Fric_M ', 1, (f_counter+f_counter_periodic*2)*6, ' float'
  k = 0
  do i=1, nb_ligneCONTACT
    if (TAB_CONTACT(i)%nature /= 'SPSPx') cycle

    l_cdt = TAB_CONTACT(i)%icdent
    l_ant = TAB_CONTACT(i)%ianent

    if (l_cdt .gt. n_particles .or. l_ant .gt. n_particles) cycle

    ! Far from the walls 
    if (TAB_SPHERE(l_cdt)%center(3) .lt. pos_z_w1 + ave_rad*disp_wall_times) cycle 
    if (TAB_SPHERE(l_cdt)%center(3) .gt. pos_z_w2 - ave_rad*disp_wall_times) cycle 

    if (TAB_SPHERE(l_ant)%center(3) .lt. pos_z_w1 + ave_rad*disp_wall_times) cycle 
    if (TAB_SPHERE(l_ant)%center(3) .gt. pos_z_w2 - ave_rad*disp_wall_times) cycle 

    Lik(:) = TAB_SPHERE(l_cdt)%center(:)-TAB_SPHERE(l_ant)%center(:)

    if (TAB_CONTACT(i)%rn .le. 0) cycle 
    
    ! Attention.. periodic case
    if (sqrt(Lik(1)**2 + Lik(2)**2 + Lik(3)**2) .gt. (TAB_SPHERE(l_cdt)%Rmax + TAB_SPHERE(l_ant)%Rmax)) then 
      do j=1, 12
        k=k+1
        if(k .le. 9) then
          ! Careful--> I'm writing the coeficient of friction explicitly 
          write(114, '(F15.7)', advance='no') sqrt(TAB_CONTACT(i)%rt**2 + TAB_CONTACT(i)%rs**2)/(0.4*TAB_CONTACT(i)%rn)
        else 
          write(114, '(F15.7)') sqrt(TAB_CONTACT(i)%rt**2 + TAB_CONTACT(i)%rs**2)/(0.4*TAB_CONTACT(i)%rn)
          k=0
        end if 
      end do
    else 
      do j=1, 6
        k=k+1
        if(k .le. 9) then
          ! Careful--> I'm writing the coeficient of friction explicitly 
          write(114, '(F15.7)', advance='no') sqrt(TAB_CONTACT(i)%rt**2 + TAB_CONTACT(i)%rs**2)/(0.4*TAB_CONTACT(i)%rn)
        else 
          write(114, '(F15.7)') sqrt(TAB_CONTACT(i)%rt**2 + TAB_CONTACT(i)%rs**2)/(0.4*TAB_CONTACT(i)%rn)
          k=0
        end if 
      end do
    end if 
  end do
  ! And jump
  write(114, '(A)') ' '

  close(114)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!  FRICTION   !!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!  mobilization !!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  vtk_friction = './POSTPRO/VTK/fricMo_' // vtk_counter // '.vtk'
  
  open(unit=114, file=vtk_friction, status='replace')
  write(114,'(A)') '# vtk DataFile Version 3.0'
  write(114,'(A,I6)') 'FRICTION ', compteur_clout
  write(114,'(A)') 'ASCII'
  write(114,'(A)') 'DATASET POLYDATA'

  ! Writing the number of vertices (Each sphere ---> 100 vertices)
  sphere_n_vertices = 74
  sphere_n_faces = 144

  ! Finding the number of real contacts to draw
  real_ctc = 0 
  do i=1, nb_ligneCONTACT
    if (TAB_CONTACT(i)%nature /= 'SPSPx') cycle
    l_cdt = TAB_CONTACT(i)%icdent
    l_ant = TAB_CONTACT(i)%ianent

    if (l_cdt .gt. n_particles .or. l_ant .gt. n_particles) cycle

    ! Far from the walls 
    if (TAB_SPHERE(l_cdt)%center(3) .lt. pos_z_w1 + ave_rad*disp_wall_times) cycle 
    if (TAB_SPHERE(l_cdt)%center(3) .gt. pos_z_w2 - ave_rad*disp_wall_times) cycle 

    if (TAB_SPHERE(l_ant)%center(3) .lt. pos_z_w1 + ave_rad*disp_wall_times) cycle 
    if (TAB_SPHERE(l_ant)%center(3) .gt. pos_z_w2 - ave_rad*disp_wall_times) cycle 

    if (TAB_CONTACT(i)%rn .le. 0) cycle

    Lik(:) = TAB_SPHERE(l_cdt)%center(:)-TAB_SPHERE(l_ant)%center(:)
    ! Not joining contacts in the periodic zones
    !if (sqrt(Lik(1)**2 + Lik(2)**2 + Lik(3)**2) .gt. (TAB_SPHERE(l_cdt)%Rmax + TAB_SPHERE(l_ant)%Rmax)) cycle

    real_ctc = real_ctc + 1
  end do 

  n_l_vertices = real_ctc * sphere_n_vertices
  
  write(v_n_vertices, '(I8)') n_l_vertices
  
  vtk_n_points = 'POINTS ' // v_n_vertices // ' float'
  write(114,'(A)') vtk_n_points


  ! Writing the coordinates of the vertices of particles and walls
  ! 3 vertices per row
  k=0
  do i=1, nb_ligneCONTACT
    ! The current vertices of this sphere with its corresponding size in global coordinates 
    if (allocated(local_sphere)) deallocate(local_sphere)
    allocate(local_sphere(sphere_n_vertices,3))

    if (TAB_CONTACT(i)%nature /= 'SPSPx') cycle
    l_cdt = TAB_CONTACT(i)%icdent
    l_ant = TAB_CONTACT(i)%ianent

    if (l_cdt .gt. n_particles .or. l_ant .gt. n_particles) cycle

    ! Far from the walls 
    if (TAB_SPHERE(l_cdt)%center(3) .lt. pos_z_w1 + ave_rad*disp_wall_times) cycle 
    if (TAB_SPHERE(l_cdt)%center(3) .gt. pos_z_w2 - ave_rad*disp_wall_times) cycle 

    if (TAB_SPHERE(l_ant)%center(3) .lt. pos_z_w1 + ave_rad*disp_wall_times) cycle 
    if (TAB_SPHERE(l_ant)%center(3) .gt. pos_z_w2 - ave_rad*disp_wall_times) cycle 

    if (TAB_CONTACT(i)%rn .le. 0) cycle

    Lik(:) = TAB_SPHERE(l_cdt)%center(:)-TAB_SPHERE(l_ant)%center(:)
    ! Not joining contacts in the periodic zones
    !if (sqrt(Lik(1)**2 + Lik(2)**2 + Lik(3)**2) .gt. (TAB_SPHERE(l_cdt)%Rmax + TAB_SPHERE(l_ant)%Rmax)) cycle

    ! Size depends on the mobilization --> Max size at the beginning 
    do j=1, sphere_n_vertices
      local_sphere(j,1:3) = (ave_rad*0.2*sphere_vertices(j,1:3)*(sqrt(TAB_CONTACT(i)%rs**2 + TAB_CONTACT(i)%rt**2) / & 
                              (0.4*TAB_CONTACT(i)%rn))) + TAB_CONTACT(i)%coor_ctc(1:3)
    end do 

    do j=1, sphere_n_vertices
      k = k + 1
      if(k .lt. 3) then 
        if (local_sphere(j,1) .lt. 0) then
          write(114,'(F11.8,A)', advance = 'no') local_sphere(j,1), ' '
        else 
          write(114,'(F9.7,A)', advance = 'no') local_sphere(j,1), ' '
        end if
        
        if (local_sphere(j,2) .lt. 0) then
          write(114,'(F11.8,A)', advance = 'no') local_sphere(j,2), ' '
        else 
          write(114,'(F9.7,A)', advance = 'no') local_sphere(j,2), ' '
        end if
      
        if (local_sphere(j,3) .lt. 0) then
          write(114,'(F11.8,A)', advance = 'no') local_sphere(j,3), ' '
        else 
          write(114,'(F9.7,A)', advance = 'no') local_sphere(j,3), ' '
        end if
      else
        if (local_sphere(j,1) .lt. 0) then
          write(114,'(F11.8,A)', advance = 'no') local_sphere(j,1), ' '
        else 
          write(114,'(F9.7,A)', advance = 'no') local_sphere(j,1), ' '
        end if
        
        if (local_sphere(j,2) .lt. 0) then
          write(114,'(F11.8,A)', advance = 'no') local_sphere(j,2), ' '
        else 
          write(114,'(F9.7,A)', advance = 'no') local_sphere(j,2), ' '
        end if
        
        if (local_sphere(j,3) .lt. 0) then
          write(114,'(F11.8,A)', advance = 'no') local_sphere(j,3), ' '
        else 
          write(114,'(F9.7,A)') local_sphere(j,3), ' '
        end if
        k = 0
      end if 
    end do
  end do

  write(114, '(A)') ' '
  
  ! Finding the total number of faces
  n_l_faces = real_ctc * sphere_n_faces

  ! Writing the conectivity between vertices
  write(114,'(A)', advance='no') 'POLYGONS '
  write(114, '(2(I8,A))') n_l_faces, ' ' , (n_l_faces*4)
  
  curr_l_faces = 0

  do i=1, nb_ligneCONTACT
    if (TAB_CONTACT(i)%nature /= 'SPSPx') cycle
    l_cdt = TAB_CONTACT(i)%icdent
    l_ant = TAB_CONTACT(i)%ianent

    if (l_cdt .gt. n_particles .or. l_ant .gt. n_particles) cycle

    ! Far from the walls 
    if (TAB_SPHERE(l_cdt)%center(3) .lt. pos_z_w1 + ave_rad*disp_wall_times) cycle 
    if (TAB_SPHERE(l_cdt)%center(3) .gt. pos_z_w2 - ave_rad*disp_wall_times) cycle 

    if (TAB_SPHERE(l_ant)%center(3) .lt. pos_z_w1 + ave_rad*disp_wall_times) cycle 
    if (TAB_SPHERE(l_ant)%center(3) .gt. pos_z_w2 - ave_rad*disp_wall_times) cycle 

    if (TAB_CONTACT(i)%rn .le. 0) cycle

    Lik(:) = TAB_SPHERE(l_cdt)%center(:)-TAB_SPHERE(l_ant)%center(:)
    ! Not joining contacts in the periodic zones
    !if (sqrt(Lik(1)**2 + Lik(2)**2 + Lik(3)**2) .gt. (TAB_SPHERE(l_cdt)%Rmax + TAB_SPHERE(l_ant)%Rmax)) cycle

    do j=1, sphere_n_faces
      ! We have always triangles
      write(114, '(I1,A)', advance= 'no') 3, ' '
      
      ! First face... writing precisely 
      if (sphere_connect(j,1) + curr_l_faces -1 .lt. 10) then
        write(114, '(I1)', advance = 'no') sphere_connect(j,1) + curr_l_faces -1
      else if (sphere_connect(j,1) + curr_l_faces -1 .lt. 100) then
        write(114, '(I2)', advance = 'no') sphere_connect(j,1) + curr_l_faces -1
      else if (sphere_connect(j,1) + curr_l_faces -1 .lt. 1000) then
        write(114, '(I3)', advance = 'no') sphere_connect(j,1) + curr_l_faces -1
      else if (sphere_connect(j,1) + curr_l_faces -1 .lt. 10000) then
        write(114, '(I4)', advance = 'no') sphere_connect(j,1) + curr_l_faces -1
      else if (sphere_connect(j,1) + curr_l_faces -1 .lt. 100000) then
        write(114, '(I5)', advance = 'no') sphere_connect(j,1) + curr_l_faces -1
      else if (sphere_connect(j,1) + curr_l_faces -1 .lt. 1000000) then
        write(114, '(I6)', advance = 'no') sphere_connect(j,1) + curr_l_faces -1
      else if (sphere_connect(j,1) + curr_l_faces -1 .lt. 10000000) then
        write(114, '(I7)', advance = 'no') sphere_connect(j,1) + curr_l_faces -1
      end if
      
      write(114, '(A)', advance='no') ' '
      
      ! Second face
      if (sphere_connect(j,2) + curr_l_faces -1 .lt. 10) then
        write(114, '(I1)', advance = 'no') sphere_connect(j,2) + curr_l_faces -1
      else if (sphere_connect(j,2) + curr_l_faces -1 .lt. 100) then
        write(114, '(I2)', advance = 'no') sphere_connect(j,2) + curr_l_faces -1
      else if (sphere_connect(j,2) + curr_l_faces -1 .lt. 1000) then
        write(114, '(I3)', advance = 'no') sphere_connect(j,2) + curr_l_faces -1
      else if (sphere_connect(j,2) + curr_l_faces -1 .lt. 10000) then
        write(114, '(I4)', advance = 'no') sphere_connect(j,2) + curr_l_faces -1
      else if (sphere_connect(j,2) + curr_l_faces -1 .lt. 100000) then
        write(114, '(I5)', advance = 'no') sphere_connect(j,2) + curr_l_faces -1
      else if (sphere_connect(j,2) + curr_l_faces -1 .lt. 1000000) then
        write(114, '(I6)', advance = 'no') sphere_connect(j,2) + curr_l_faces -1
      else if (sphere_connect(j,2) + curr_l_faces -1 .lt. 10000000) then
        write(114, '(I7)', advance = 'no') sphere_connect(j,2) + curr_l_faces -1
      end if
      
      ! Third face
      write(114, '(A)', advance='no') ' '
      
      if (sphere_connect(j,3) + curr_l_faces -1 .lt. 10) then
        write(114, '(I1)') sphere_connect(j,3) + curr_l_faces -1
      else if (sphere_connect(j,3) + curr_l_faces -1 .lt. 100) then
        write(114, '(I2)') sphere_connect(j,3) + curr_l_faces -1
      else if (sphere_connect(j,3) + curr_l_faces -1 .lt. 1000) then
        write(114, '(I3)') sphere_connect(j,3) + curr_l_faces -1
      else if (sphere_connect(j,3) + curr_l_faces -1 .lt. 10000) then
        write(114, '(I4)') sphere_connect(j,3) + curr_l_faces -1
      else if (sphere_connect(j,3) + curr_l_faces -1 .lt. 100000) then
        write(114, '(I5)') sphere_connect(j,3) + curr_l_faces -1
      else if (sphere_connect(j,3) + curr_l_faces -1 .lt. 1000000) then
        write(114, '(I6)') sphere_connect(j,3) + curr_l_faces -1
      else if (sphere_connect(j,3) + curr_l_faces -1 .lt. 10000000) then
        write(114, '(I7)') sphere_connect(j,3) + curr_l_faces -1
      end if
      
    end do
    curr_l_faces = curr_l_faces + sphere_n_vertices
  end do

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!! FIELDS !!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! How many fields there will be?
  n_l_fields = 1                          ! Friction Mobilization
                                          ! ...
  
  ! A blank space
  write(114,'(A)') ''
  ! The cells begin
  write(114,'(A)', advance='no') 'CELL_DATA'
  
  ! Writing the number of data by field. It corresponds to the same number of faces
  write(114, '(2(I8,A))') n_l_faces
  write(114, '(A,I4)')  'FIELD FieldData', n_l_fields


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Friction mobilization (Rs+Rt)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! Naming the field, the dimension of the data, and number of lines (as before FIELD), and data type
  write(114,'(A,I1,I8,A)') 'Fric_M ', 1, n_l_faces, ' float'
  k = 0
  do i=1, nb_ligneCONTACT
    if (TAB_CONTACT(i)%nature /= 'SPSPx') cycle

    l_cdt = TAB_CONTACT(i)%icdent
    l_ant = TAB_CONTACT(i)%ianent

    if (l_cdt .gt. n_particles .or. l_ant .gt. n_particles) cycle

    ! Far from the walls 
    if (TAB_SPHERE(l_cdt)%center(3) .lt. pos_z_w1 + ave_rad*disp_wall_times) cycle 
    if (TAB_SPHERE(l_cdt)%center(3) .gt. pos_z_w2 - ave_rad*disp_wall_times) cycle 

    if (TAB_SPHERE(l_ant)%center(3) .lt. pos_z_w1 + ave_rad*disp_wall_times) cycle 
    if (TAB_SPHERE(l_ant)%center(3) .gt. pos_z_w2 - ave_rad*disp_wall_times) cycle 

    Lik(:) = TAB_SPHERE(l_cdt)%center(:)-TAB_SPHERE(l_ant)%center(:)

    if (TAB_CONTACT(i)%rn .le. 0) cycle 
    
    ! Not joining contacts in the periodic zones
    !if (sqrt(Lik(1)**2 + Lik(2)**2 + Lik(3)**2) .gt. (TAB_SPHERE(l_cdt)%Rmax + TAB_SPHERE(l_ant)%Rmax)) cycle

    do j=1, sphere_n_faces
      k=k+1
      if(k .le. 9) then
        ! Careful--> I'm writing the coeficient of friction explicitly 
        write(114, '(F15.7)', advance='no') sqrt(TAB_CONTACT(i)%rt**2 + TAB_CONTACT(i)%rs**2)/(0.4*TAB_CONTACT(i)%rn)
      else 
        write(114, '(F15.7)') sqrt(TAB_CONTACT(i)%rt**2 + TAB_CONTACT(i)%rs**2)/(0.4*TAB_CONTACT(i)%rn)
        k=0
      end if 
    end do
  end do
  ! And jump
  write(114, '(A)') ' '

  close(114)

  print*, 'Drawing in vtk  ---> Ok!'

end subroutine draw


!==============================================================================
! List of normal forces 
!==============================================================================
subroutine list_n_forces
  
  implicit none

  integer                                  :: i, cd, an
  real*8, dimension(3)                     :: Lik 
  character(len=24)                        :: l_file

  real*8                                   :: pos_z_w1, pos_z_w2
  
  pos_z_w1 = TAB_PLAN(1)%center(3)
  pos_z_w2 = TAB_PLAN(2)%center(3)

  ! The name of the file
  l_file = './POSTPRO/LIST__FN.     '

  ! Preparing the corresponding number of the file
  if (compteur_clout<10) then
    WRITE(l_file(20:21),'(I1)')   compteur_clout
  else if ( (compteur_clout>=10) .and. (compteur_clout<100) ) then
    WRITE(l_file(20:22),'(I2)')   compteur_clout
  else if ( (compteur_clout>=100).and. (compteur_clout<1000) ) then
    WRITE(l_file(20:23),'(I3)')   compteur_clout
  else if ( (compteur_clout>=1000).and. (compteur_clout<10000) ) then
    WRITE(l_file(20:24),'(I4)')   compteur_clout
  else
    print*, "Not enough file names :: List of normal forces"
    stop
  end if 

  ! Opening the file
  open(unit=115,file=l_file,status='replace')

  write(115,*) '#   FN     '

  do i=1,nb_ligneCONTACT
    ! Contacts between particles
    if (TAB_CONTACT(i)%nature /= 'SPSPx') cycle

    ! Active contacts
    if (TAB_CONTACT(i)%rn .le. 0.0) cycle
    
    cd   = TAB_CONTACT(i)%icdent
    an   = TAB_CONTACT(i)%ianent

    if ((cd .gt. n_particles) .or. (an .gt. n_particles)) cycle

    ! Counting a smaller box letting a space to the walls
    if (TAB_SPHERE(cd)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
    if (TAB_SPHERE(cd)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle

    if (TAB_SPHERE(an)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
    if (TAB_SPHERE(an)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle

    Lik(:) = TAB_SPHERE(cd)%center(:)-TAB_SPHERE(an)%center(:)

    ! Not joining contacts in the periodic zones
    if (sqrt(Lik(1)**2 + Lik(2)**2 + Lik(3)**2) .gt. (TAB_SPHERE(cd)%Rmax + TAB_SPHERE(an)%Rmax)) cycle

    write(115, '(F15.7)') TAB_CONTACT(i)%rn
  end do

  close(115)
  
  print*, 'Write List FN     ---> Ok!'

end subroutine list_n_forces


!==============================================================================
! List of tangential forces 
!==============================================================================
subroutine list_t_forces
  
  implicit none

  integer                                  :: i, cd, an
  real*8, dimension(3)                     :: Lik 
  character(len=24)                        :: l_file

  real*8                                   :: pos_z_w1, pos_z_w2
  
  pos_z_w1 = TAB_PLAN(1)%center(3)
  pos_z_w2 = TAB_PLAN(2)%center(3)

  ! The name of the file
  l_file = './POSTPRO/LIST__FT.     '

  ! Preparing the corresponding number of the file
  if (compteur_clout<10) then
    WRITE(l_file(20:21),'(I1)')   compteur_clout
  else if ( (compteur_clout>=10) .and. (compteur_clout<100) ) then
    WRITE(l_file(20:22),'(I2)')   compteur_clout
  else if ( (compteur_clout>=100).and. (compteur_clout<1000) ) then
    WRITE(l_file(20:23),'(I3)')   compteur_clout
  else if ( (compteur_clout>=1000).and. (compteur_clout<10000) ) then
    WRITE(l_file(20:24),'(I4)')   compteur_clout
  else
    print*, "Not enough file names :: List of normal forces"
    stop
  end if 

  ! Opening the file
  open(unit=125,file=l_file,status='replace')

  write(125,*) '#   FT     '

  do i=1,nb_ligneCONTACT
    ! Contacts between particles
    if (TAB_CONTACT(i)%nature /= 'SPSPx') cycle

    ! Active contacts
    if (TAB_CONTACT(i)%rn .le. 0.0) cycle
    
    cd   = TAB_CONTACT(i)%icdent
    an   = TAB_CONTACT(i)%ianent

    if ((cd .gt. n_particles) .or. (an .gt. n_particles)) cycle

    ! Counting a smaller box letting a space to the walls
    if (TAB_SPHERE(cd)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
    if (TAB_SPHERE(cd)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle

    if (TAB_SPHERE(an)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
    if (TAB_SPHERE(an)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle

    Lik(:) = TAB_SPHERE(cd)%center(:)-TAB_SPHERE(an)%center(:)

    ! Not joining contacts in the periodic zones
    if (sqrt(Lik(1)**2 + Lik(2)**2 + Lik(3)**2) .gt. (TAB_SPHERE(cd)%Rmax + TAB_SPHERE(an)%Rmax)) cycle

    write(125, '(F15.7)') TAB_CONTACT(i)%rt
  end do

  close(125)
  
  print*, 'Write List FT     ---> Ok!'

end subroutine list_t_forces


!==============================================================================
! List of normal forces 
!==============================================================================
subroutine list_branches
  
  implicit none

  integer                                  :: i, cd, an
  real*8, dimension(3)                     :: Lik
  character(len=24)                        :: l_file

  real*8                                   :: pos_z_w1, pos_z_w2
  
  pos_z_w1 = TAB_PLAN(1)%center(3)
  pos_z_w2 = TAB_PLAN(2)%center(3)

  ! The name of the file
  l_file = './POSTPRO/LIST_Lik.     '

  ! Preparing the corresponding number of the file
  if (compteur_clout<10) then
    WRITE(l_file(20:21),'(I1)')   compteur_clout
  else if ( (compteur_clout>=10) .and. (compteur_clout<100) ) then
    WRITE(l_file(20:22),'(I2)')   compteur_clout
  else if ( (compteur_clout>=100).and. (compteur_clout<1000) ) then
    WRITE(l_file(20:23),'(I3)')   compteur_clout
  else if ( (compteur_clout>=1000).and. (compteur_clout<10000) ) then
    WRITE(l_file(20:24),'(I4)')   compteur_clout
  else
    print*, "Not enough file names :: List of branches"
    stop
  end if 

  ! Opening the file
  open(unit=116,file=l_file,status='replace')

  write(116,*) '#  Branch    '

  do i=1,nb_ligneCONTACT
    ! Contacts between particles
    if (TAB_CONTACT(i)%nature /= 'SPSPx') cycle

    ! Active contacts
    if (TAB_CONTACT(i)%rn .le. 0.0) cycle
    
    cd   = TAB_CONTACT(i)%icdent
    an   = TAB_CONTACT(i)%ianent

    if ((cd .gt. n_particles) .or. (an .gt. n_particles)) cycle

    ! Counting a smaller box letting a space to the walls
    if (TAB_SPHERE(cd)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
    if (TAB_SPHERE(cd)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle

    if (TAB_SPHERE(an)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
    if (TAB_SPHERE(an)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle

    Lik(:) = TAB_SPHERE(cd)%center(:)-TAB_SPHERE(an)%center(:)

    ! Not joining contacts in the periodic zones
    if (sqrt(Lik(1)**2 + Lik(2)**2 + Lik(3)**2) .gt. (TAB_SPHERE(cd)%Rmax + TAB_SPHERE(an)%Rmax)) cycle

    write(116, '(F15.7)') sqrt(Lik(1)**2 + Lik(2)**2 + Lik(3)**2)
  end do
  
  close(116)

  print*, 'Write List Lik     ---> Ok!'

end subroutine list_branches


!==============================================================================
! List of friction mobilization
!==============================================================================
subroutine list_mobilization
  
  implicit none

  integer                                  :: i, cd, an
  real*8, dimension(3)                     :: Lik
  character(len=24)                        :: l_file

  real*8                                   :: pos_z_w1, pos_z_w2
  
  pos_z_w1 = TAB_PLAN(1)%center(3)
  pos_z_w2 = TAB_PLAN(2)%center(3)

  ! The name of the file
  l_file = './POSTPRO/LIST_MFr.     '

  ! Preparing the corresponding number of the file
  if (compteur_clout<10) then
    WRITE(l_file(20:21),'(I1)')   compteur_clout
  else if ( (compteur_clout>=10) .and. (compteur_clout<100) ) then
    WRITE(l_file(20:22),'(I2)')   compteur_clout
  else if ( (compteur_clout>=100).and. (compteur_clout<1000) ) then
    WRITE(l_file(20:23),'(I3)')   compteur_clout
  else if ( (compteur_clout>=1000).and. (compteur_clout<10000) ) then
    WRITE(l_file(20:24),'(I4)')   compteur_clout
  else
    print*, "Not enough file names :: List of friction mobilization"
    stop
  end if 

  ! Opening the file
  open(unit=117,file=l_file,status='replace')

  write(117,*) '#  Branch    '

  do i=1,nb_ligneCONTACT
    ! Contacts between particles
    if (TAB_CONTACT(i)%nature /= 'SPSPx') cycle

    ! Active contacts
    if (TAB_CONTACT(i)%rn .le. 0.0) cycle
    
    cd   = TAB_CONTACT(i)%icdent
    an   = TAB_CONTACT(i)%ianent

    if ((cd .gt. n_particles) .or. (an .gt. n_particles)) cycle

    ! Counting a smaller box letting a space to the walls
    if (TAB_SPHERE(cd)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
    if (TAB_SPHERE(cd)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle

    if (TAB_SPHERE(an)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
    if (TAB_SPHERE(an)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle

    Lik(:) = TAB_SPHERE(cd)%center(:)-TAB_SPHERE(an)%center(:)

    ! Not joining contacts in the periodic zones
    if (sqrt(Lik(1)**2 + Lik(2)**2 + Lik(3)**2) .gt. (TAB_SPHERE(cd)%Rmax + TAB_SPHERE(an)%Rmax)) cycle

    ! Attention ---> coeficient of friction explicitly written 
    write(117, '(F15.7)') sqrt(TAB_CONTACT(i)%rs**2 + TAB_CONTACT(i)%rt**2)/(0.4*TAB_CONTACT(i)%rn)
  end do

  close(117)
  
  print*, 'Write List mobilization ---> Ok!'

end subroutine list_mobilization

!==============================================================================
! Computing the normal force by class 
!==============================================================================
subroutine fn_class 
  
  implicit none
  
  integer                                  :: i, j, cd, an, n_bin
  real*8                                   :: r_max, r_min, interval
  real*8, dimension(:,:), allocatable      :: bins
  character(len=24)                        :: coor_c_file

  real*8                                   :: pos_z_w1, pos_z_w2

  pos_z_w1 = TAB_PLAN(1)%center(3)
  pos_z_w2 = TAB_PLAN(2)%center(3)  

  ! The name of the file
  coor_c_file = './POSTPRO/FN_CLASS.     '

  ! Preparing the corresponding number of the file
  if (compteur_clout<10) then
    WRITE(coor_c_file(20:21),'(I1)')   compteur_clout
  else if ( (compteur_clout>=10) .and. (compteur_clout<100) ) then
    WRITE(coor_c_file(20:22),'(I2)')   compteur_clout
  else if ( (compteur_clout>=100).and. (compteur_clout<1000) ) then
    WRITE(coor_c_file(20:23),'(I3)')   compteur_clout
  else if ( (compteur_clout>=1000).and. (compteur_clout<10000) ) then
    WRITE(coor_c_file(20:24),'(I4)')   compteur_clout
  else
    print*, "Not enough file names :: Coordination per class"
    stop
  end if 

  ! Opening the file
  open(unit=118,file=coor_c_file,status='replace')

  ! Initalizing
  r_max = -9999.
  r_min =  9999.

  ! Finding the max and min radius
  do i=1, n_particles
    ! Particles far from the walls
    if (TAB_SPHERE(i)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle 
    if (TAB_SPHERE(i)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle 

    r_min = min(r_min, TAB_SPHERE(i)%Rmax)
    r_max = max(r_max, TAB_SPHERE(i)%Rmax)
  end do 

  ! The number of bins - > Intervals in the size distribution
  n_bin = 36

  ! The size of the interval 
  interval = (r_max - r_min)/n_bin

  ! Allocating the bins 
  if (allocated(bins)) deallocate(bins)
  allocate(bins(n_bin,3))
  ! initializing
  bins(:,:) = 0

  ! Adding particles to bins as function of their size. First column: average size. 
  ! Second: number of contacts. Third: Total number of particles in this interval 

  ! Caracteristic size of each bin 
  do i=1, n_bin
    bins(i,1) = r_min + interval*(i - 0.5)
  end do

  do i=1, n_bin
    do j=1, n_particles
      if((TAB_SPHERE(j)%Rmax .ge. (r_min + (i-1)*interval)) .and. (TAB_SPHERE(j)%Rmax*0.999 .lt. (interval*i + r_min))) then 

        if (TAB_SPHERE(j)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle 
        if (TAB_SPHERE(j)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle 

        bins(i,2) = bins(i,2) + TAB_SPHERE(j)%av_fn
        bins(i,3) = bins(i,3) + 1
      end if
    end do
  end do

  ! Writing
  write(118,*) '#   Class     ', '    <fn>     '

  !print*, bins  
  
  do i=1, n_bin
    if (bins(1,3) .gt. 0.) then 
      write(118,'(2(1X,E12.5))') bins(i,1), real(bins(i,2)/bins(i,3))
    end if 
  end do 

  close(118)
                                              
  print*, 'Normal force Class  ---> Ok!'
  
end subroutine fn_class

!==============================================================================
! Computing the coordination by class 
!==============================================================================
subroutine mobilization_class 
  
  implicit none
  
  integer                                  :: i, j, cd, an, n_bin
  real*8                                   :: r_max, r_min, interval
  real*8, dimension(:,:), allocatable      :: bins
  character(len=24)                        :: coor_c_file

  real*8                                   :: pos_z_w1, pos_z_w2

  pos_z_w1 = TAB_PLAN(1)%center(3)
  pos_z_w2 = TAB_PLAN(2)%center(3)  

  ! The name of the file
  coor_c_file = './POSTPRO/MF_CLASS.     '

  ! Preparing the corresponding number of the file
  if (compteur_clout<10) then
    WRITE(coor_c_file(20:21),'(I1)')   compteur_clout
  else if ( (compteur_clout>=10) .and. (compteur_clout<100) ) then
    WRITE(coor_c_file(20:22),'(I2)')   compteur_clout
  else if ( (compteur_clout>=100).and. (compteur_clout<1000) ) then
    WRITE(coor_c_file(20:23),'(I3)')   compteur_clout
  else if ( (compteur_clout>=1000).and. (compteur_clout<10000) ) then
    WRITE(coor_c_file(20:24),'(I4)')   compteur_clout
  else
    print*, "Not enough file names :: Coordination per class"
    stop
  end if 

  ! Opening the file
  open(unit=119,file=coor_c_file,status='replace')

  ! Initalizing
  r_max = -9999.
  r_min =  9999.

  ! Finding the max and min radius
  do i=1, n_particles
    ! Particles far from the walls
    if (TAB_SPHERE(i)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle 
    if (TAB_SPHERE(i)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle 

    r_min = min(r_min, TAB_SPHERE(i)%Rmax)
    r_max = max(r_max, TAB_SPHERE(i)%Rmax)
  end do 

  ! The number of bins - > Intervals in the size distribution
  n_bin = 36

  ! The size of the interval 
  interval = (r_max - r_min)/n_bin

  ! Allocating the bins 
  if (allocated(bins)) deallocate(bins)
  allocate(bins(n_bin,3))
  ! initializing
  bins(:,:) = 0

  ! Adding particles to bins as function of their size. First column: average size. 
  ! Second: number of contacts. Third: Total number of particles in this interval 

  ! Caracteristic size of each bin 
  do i=1, n_bin
    bins(i,1) = r_min + interval*(i - 0.5)
  end do

  do i=1, n_bin
    do j=1, n_particles
      if((TAB_SPHERE(j)%Rmax .ge. (r_min + (i-1)*interval)) .and. (TAB_SPHERE(j)%Rmax*0.999 .lt. (interval*i + r_min))) then 

        if (TAB_SPHERE(j)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle 
        if (TAB_SPHERE(j)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle 

        bins(i,2) = bins(i,2) + TAB_SPHERE(j)%av_mob
        bins(i,3) = bins(i,3) + 1
      end if
    end do
  end do

  ! Writing
  write(119,*) '#   Class     ', ' <|ft|/mu fn>   '

  !print*, bins

  do i=1, n_bin
    if (bins(1,3) .gt. 0.) then
      write(119,'(2(1X,E12.5))') bins(i,1), real(bins(i,2)/bins(i,3))
    end if
  end do

  close(119)

  print*, 'Mobilization Class  ---> Ok!'

end subroutine mobilization_class


!==============================================================================
! Computing anisotropies per class of particles by branch
!==============================================================================
subroutine aniso_class

  implicit none

  integer                                  :: i,j,cd,an, n_bin
  integer                                  :: ierror,matz,lda, valid_c
  real*8                                   :: S1, S2, S3, cpt, ac, interval, norm_FN, norm_FT
  real*8                                   :: b_max, b_min, Rsik, Rnik, Rtik, norm_lik
  real*8                                   :: X1f, X1n, X2f, X2n, X3f, X3n
  real*8,dimension(3)                      :: nik, Fik, sik, tik, FikT, dirFT, Lik
  real*8,dimension(3)                      :: wr,wi, wrf, wrn, win, wif
  real*8,dimension(3,3)                    :: localframe, localframeN, localframeF
  real*8,dimension(:,:), allocatable       :: bins
  real*8,dimension(3,3)                    :: Fabric, Fabric_N, Fabric_T
  real*8,dimension(:,:), allocatable       :: all_fabric, all_forces_n, all_forces_t
  character(len=24)                        :: aniso_c_file

  real*8                                   :: pos_z_w1, pos_z_w2

  pos_z_w1 = TAB_PLAN(1)%center(3)
  pos_z_w2 = TAB_PLAN(2)%center(3)

  ! The name of the file
  aniso_c_file = './POSTPRO/AN_CLASS.     '

  ! Preparing the corresponding number of the file
  if (compteur_clout<10) then
    WRITE(aniso_c_file(20:21),'(I1)')   compteur_clout
  else if ( (compteur_clout>=10) .and. (compteur_clout<100) ) then
    WRITE(aniso_c_file(20:22),'(I2)')   compteur_clout
  else if ( (compteur_clout>=100).and. (compteur_clout<1000) ) then
    WRITE(aniso_c_file(20:23),'(I3)')   compteur_clout
  else if ( (compteur_clout>=1000).and. (compteur_clout<10000) ) then
    WRITE(aniso_c_file(20:24),'(I4)')   compteur_clout
  else
    print*, "Not enough file names :: Anisotropies per class"
    stop
  end if

  ! Opening the file
  open(unit=123,file=aniso_c_file,status='replace')

  ! Initializing
  b_max = -9999.
  b_min =  9999.

  ! Finding the max and min branch
  do i=1, nb_ligneCONTACT
    ! Active contacts
    if (TAB_CONTACT(i)%rn .le. 0.0) cycle

    cd = TAB_CONTACT(i)%icdent
    an = TAB_CONTACT(i)%ianent

    if ((cd .gt. n_particles) .or. (an .gt. n_particles)) cycle

    ! Counting a smaller box letting a space to the walls
    if (TAB_SPHERE(cd)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
    if (TAB_SPHERE(cd)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle

    if (TAB_SPHERE(an)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
    if (TAB_SPHERE(an)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle

    nik  = TAB_CONTACT(i)%n
    tik  = TAB_CONTACT(i)%t
    sik  = TAB_CONTACT(i)%s
    Rtik = TAB_CONTACT(i)%rt
    Rnik = TAB_CONTACT(i)%rn
    Rsik = TAB_CONTACT(i)%rs

    ! The branch
    Lik(1:3) = TAB_SPHERE(cd)%center(1:3)-TAB_SPHERE(an)%center(1:3)

    if (sqrt(Lik(1)**2 + Lik(2)**2 + Lik(3)**2) > 3*(TAB_SPHERE(cd)%Rmax + TAB_SPHERE(an)%Rmax)) then
      !Rayon_max = sqrt((Lik(1)**2 + Lik(2)**2 + Lik(3)**2))
      !Lik(:) = Lik(:) / Rayon_max
      Lik(:) = abs(TAB_SPHERE(cd)%Rmax+TAB_SPHERE(an)%Rmax)*nik(:)
    end if

    b_min = min(b_min, (Lik(1)**2 + Lik(2)**2 + Lik(3)**2)**0.5)
    b_max = max(b_max, (Lik(1)**2 + Lik(2)**2 + Lik(3)**2)**0.5)
  end do

  ! The number of bins - > Intervals in the size distribution
  n_bin = 36

  ! The size of the interval
  interval = (b_max - b_min)/n_bin

  ! Allocating the bins
  if (allocated(bins)) deallocate(bins)
  allocate(bins(n_bin,29))
  ! initializing
  bins(:,:) = 0

  ! The array of fabrics
  !if (allocated(all_fabric)) deallocate(all_fabric)
  !allocate(all_fabric(n_particles,9))
  !! initializing
  !all_fabric(:,:) = 0

  ! The array of forces tensor normal
  !if (allocated(all_forces_n)) deallocate(all_forces_n)
  !allocate(all_forces_n(n_particles,9))
  !! initializing
  !all_forces_n(:,:) = 0

  ! The array of forces tensor tangential
  !if (allocated(all_forces_t)) deallocate(all_forces_t)
  !allocate(all_forces_t(n_particles,9))
  !! initializing
  !all_forces_t(:,:) = 0

  !! Computing the fabric individually for all particles
  !do i=1, nb_ligneCONTACT
  !  if(TAB_CONTACT(i)%nature == 'SPPLx') cycle
  !  cd = TAB_CONTACT(i)%icdent
  !  an = TAB_CONTACT(i)%ianent
  !  if ((cd .gt. n_particles) .or. (an .gt. n_particles)) cycle
!
  !  ! Counting a smaller box letting a space to the walls
  !  if (TAB_SPHERE(cd)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
  !  if (TAB_SPHERE(cd)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle
!
  !  if (TAB_SPHERE(an)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
  !  if (TAB_SPHERE(an)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle
!
  !  Rnik = TAB_CONTACT(i)%rn
!
  !  ! Checking if it is an active contact
  !  if (Rnik .le. 0.D0) cycle
!
  !  nik  = TAB_CONTACT(i)%n
  !  tik  = TAB_CONTACT(i)%t
  !  sik  = TAB_CONTACT(i)%s
!
  !  Rtik = TAB_CONTACT(i)%rt
  !  Rsik = TAB_CONTACT(i)%rs
!
  !  ! The force vector
  !  Fik(:) = Rnik*nik(:) + Rtik*tik(:) + Rsik*sik(:)
!
  !  ! The norm of the normal force
  !  norm_FN = Fik(1)*nik(1) + Fik(2)*nik(2) + Fik(3)*nik(3)
!
  !  ! The tangential force vector
  !  FikT(1) = Fik(1) - Rnik*nik(1)
  !  FikT(2) = Fik(2) - Rnik*nik(2)
  !  FikT(3) = Fik(3) - Rnik*nik(3)
!
  !  ! The norm of FikT
  !  norm_FT = sqrt(FikT(1)**2 + FikT(2)**2 + FikT(3)**2)
!
  !  ! The direction of FikT
  !  dirFT = FikT/norm_FT
!
  !  all_fabric(cd,1) = all_fabric(cd,1) + nik(1)*nik(1)
  !  all_fabric(cd,2) = all_fabric(cd,2) + nik(1)*nik(2)
  !  all_fabric(cd,3) = all_fabric(cd,3) + nik(1)*nik(3)
  !  all_fabric(cd,4) = all_fabric(cd,4) + nik(2)*nik(1)
  !  all_fabric(cd,5) = all_fabric(cd,5) + nik(2)*nik(2)
  !  all_fabric(cd,6) = all_fabric(cd,6) + nik(2)*nik(3)
  !  all_fabric(cd,7) = all_fabric(cd,7) + nik(3)*nik(1)
  !  all_fabric(cd,8) = all_fabric(cd,8) + nik(3)*nik(2)
  !  all_fabric(cd,9) = all_fabric(cd,9) + nik(3)*nik(3)
!
  !  all_fabric(an,1) = all_fabric(an,1) + nik(1)*nik(1)
  !  all_fabric(an,2) = all_fabric(an,2) + nik(1)*nik(2)
  !  all_fabric(an,3) = all_fabric(an,3) + nik(1)*nik(3)
  !  all_fabric(an,4) = all_fabric(an,4) + nik(2)*nik(1)
  !  all_fabric(an,5) = all_fabric(an,5) + nik(2)*nik(2)
  !  all_fabric(an,6) = all_fabric(an,6) + nik(2)*nik(3)
  !  all_fabric(an,7) = all_fabric(an,7) + nik(3)*nik(1)
  !  all_fabric(an,8) = all_fabric(an,8) + nik(3)*nik(2)
  !  all_fabric(an,9) = all_fabric(an,9) + nik(3)*nik(3)
!
!
  !  all_forces_n(cd,1) = all_forces_n(cd,1) + norm_FN*nik(1)*nik(1)
  !  all_forces_n(cd,2) = all_forces_n(cd,2) + norm_FN*nik(1)*nik(2)
  !  all_forces_n(cd,3) = all_forces_n(cd,3) + norm_FN*nik(1)*nik(3)
  !  all_forces_n(cd,4) = all_forces_n(cd,4) + norm_FN*nik(2)*nik(1)
  !  all_forces_n(cd,5) = all_forces_n(cd,5) + norm_FN*nik(2)*nik(2)
  !  all_forces_n(cd,6) = all_forces_n(cd,6) + norm_FN*nik(2)*nik(3)
  !  all_forces_n(cd,7) = all_forces_n(cd,7) + norm_FN*nik(3)*nik(1)
  !  all_forces_n(cd,8) = all_forces_n(cd,8) + norm_FN*nik(3)*nik(2)
  !  all_forces_n(cd,9) = all_forces_n(cd,9) + norm_FN*nik(3)*nik(3)
!
  !  all_forces_n(an,1) = all_forces_n(an,1) + norm_FN*nik(1)*nik(1)
  !  all_forces_n(an,2) = all_forces_n(an,2) + norm_FN*nik(1)*nik(2)
  !  all_forces_n(an,3) = all_forces_n(an,3) + norm_FN*nik(1)*nik(3)
  !  all_forces_n(an,4) = all_forces_n(an,4) + norm_FN*nik(2)*nik(1)
  !  all_forces_n(an,5) = all_forces_n(an,5) + norm_FN*nik(2)*nik(2)
  !  all_forces_n(an,6) = all_forces_n(an,6) + norm_FN*nik(2)*nik(3)
  !  all_forces_n(an,7) = all_forces_n(an,7) + norm_FN*nik(3)*nik(1)
  !  all_forces_n(an,8) = all_forces_n(an,8) + norm_FN*nik(3)*nik(2)
  !  all_forces_n(an,9) = all_forces_n(an,9) + norm_FN*nik(3)*nik(3)
!
  !  all_forces_t(cd,1) = all_forces_t(cd,1) + norm_FT*nik(1)*dirFT(1)
  !  all_forces_t(cd,2) = all_forces_t(cd,2) + norm_FT*nik(1)*dirFT(2)
  !  all_forces_t(cd,3) = all_forces_t(cd,3) + norm_FT*nik(1)*dirFT(3)
  !  all_forces_t(cd,4) = all_forces_t(cd,4) + norm_FT*nik(2)*dirFT(1)
  !  all_forces_t(cd,5) = all_forces_t(cd,5) + norm_FT*nik(2)*dirFT(2)
  !  all_forces_t(cd,6) = all_forces_t(cd,6) + norm_FT*nik(2)*dirFT(3)
  !  all_forces_t(cd,7) = all_forces_t(cd,7) + norm_FT*nik(3)*dirFT(1)
  !  all_forces_t(cd,8) = all_forces_t(cd,8) + norm_FT*nik(3)*dirFT(2)
  !  all_forces_t(cd,9) = all_forces_t(cd,9) + norm_FT*nik(3)*dirFT(3)
!
  !  all_forces_t(an,1) = all_forces_t(an,1) + norm_FT*nik(1)*dirFT(1)
  !  all_forces_t(an,2) = all_forces_t(an,2) + norm_FT*nik(1)*dirFT(2)
  !  all_forces_t(an,3) = all_forces_t(an,3) + norm_FT*nik(1)*dirFT(3)
  !  all_forces_t(an,4) = all_forces_t(an,4) + norm_FT*nik(2)*dirFT(1)
  !  all_forces_t(an,5) = all_forces_t(an,5) + norm_FT*nik(2)*dirFT(2)
  !  all_forces_t(an,6) = all_forces_t(an,6) + norm_FT*nik(2)*dirFT(3)
  !  all_forces_t(an,7) = all_forces_t(an,7) + norm_FT*nik(3)*dirFT(1)
  !  all_forces_t(an,8) = all_forces_t(an,8) + norm_FT*nik(3)*dirFT(2)
  !  all_forces_t(an,9) = all_forces_t(an,9) + norm_FT*nik(3)*dirFT(3)
  !end do

  ! Filling the branch ranges

  do i=1, n_bin
    bins(i,1) = b_min + interval*(i - 0.5)
  end do

  ! Valid number of contacts
  valid_c = 0

 ! Computing the branch and putting the corresponding stress tensor in the corresponding bin
  do i=1, n_bin
    do j=1, nb_ligneCONTACT

      ! Active contacts
      if (TAB_CONTACT(j)%rn .le. 0.0) cycle

      cd = TAB_CONTACT(j)%icdent
      an = TAB_CONTACT(j)%ianent

      if ((cd .gt. n_particles) .or. (an .gt. n_particles)) cycle

      ! Counting a smaller box letting a space to the walls
      if (TAB_SPHERE(cd)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
      if (TAB_SPHERE(cd)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle

      if (TAB_SPHERE(an)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
      if (TAB_SPHERE(an)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle

      nik  = TAB_CONTACT(j)%n
      tik  = TAB_CONTACT(j)%t
      sik  = TAB_CONTACT(j)%s
      Rtik = TAB_CONTACT(j)%rt
      Rnik = TAB_CONTACT(j)%rn
      Rsik = TAB_CONTACT(j)%rs

      ! The branch
      Lik(1:3) = TAB_SPHERE(cd)%center(1:3)-TAB_SPHERE(an)%center(1:3)

      if (sqrt(Lik(1)**2 + Lik(2)**2 + Lik(3)**2) > 3*(TAB_SPHERE(cd)%Rmax + TAB_SPHERE(an)%Rmax)) then
        !Rayon_max = sqrt((Lik(1)**2 + Lik(2)**2 + Lik(3)**2))
        !Lik(:) = Lik(:) / Rayon_max
        Lik(:) = abs(TAB_SPHERE(cd)%Rmax+TAB_SPHERE(an)%Rmax)*nik(:)
      end if

      norm_lik = (Lik(1)**2 + Lik(2)**2 + Lik(3)**2)**0.5

      ! Checking the range of branch
      if((norm_lik .ge. (b_min + (i-1)*interval)) .and. (norm_lik*0.999 .lt. (interval*i + b_min))) then

        valid_c = valid_c + 1
        ! The counter of contacts
        bins(i,2)  = bins(i,2)  + 1

        ! The fabric tensor
        bins(i,3)  = bins(i,3)  + nik(1)*nik(1)
        bins(i,4)  = bins(i,4)  + nik(1)*nik(2)
        bins(i,5)  = bins(i,5)  + nik(1)*nik(3)
        bins(i,6)  = bins(i,6)  + nik(2)*nik(1)
        bins(i,7)  = bins(i,7)  + nik(2)*nik(2)
        bins(i,8)  = bins(i,8)  + nik(2)*nik(3)
        bins(i,9)  = bins(i,9)  + nik(3)*nik(1)
        bins(i,10) = bins(i,10) + nik(3)*nik(2)
        bins(i,11) = bins(i,11) + nik(3)*nik(3)

        ! For the forces tensors
        ! The force vector
        Fik(:) = Rnik*nik(:) + Rtik*tik(:) + Rsik*sik(:)

        ! The norm of the normal force
        norm_FN = Fik(1)*nik(1) + Fik(2)*nik(2) + Fik(3)*nik(3)

        ! The tangential force vector
        FikT(1) = Fik(1) - Rnik*nik(1)
        FikT(2) = Fik(2) - Rnik*nik(2)
        FikT(3) = Fik(3) - Rnik*nik(3)

        ! The norm of FikT
        norm_FT = sqrt(FikT(1)**2 + FikT(2)**2 + FikT(3)**2)

        ! The chi tensor of normal forces
        bins(i,12)  = bins(i,12)  + norm_FN*nik(1)*nik(1)
        bins(i,13)  = bins(i,13)  + norm_FN*nik(1)*nik(2)
        bins(i,14)  = bins(i,14)  + norm_FN*nik(1)*nik(3)
        bins(i,15)  = bins(i,15)  + norm_FN*nik(2)*nik(1)
        bins(i,16)  = bins(i,16)  + norm_FN*nik(2)*nik(2)
        bins(i,17)  = bins(i,17)  + norm_FN*nik(2)*nik(3)
        bins(i,18)  = bins(i,18)  + norm_FN*nik(3)*nik(1)
        bins(i,19)  = bins(i,19)  + norm_FN*nik(3)*nik(2)
        bins(i,20)  = bins(i,20)  + norm_FN*nik(3)*nik(3)

        if (norm_FT .gt. 0.0) then

          ! The direction of FikT
          dirFT = FikT/norm_FT

          ! The chi tensor of tangential forces
          bins(i,21)  = bins(i,21)  + norm_FT*nik(1)*dirFT(1)
          bins(i,22)  = bins(i,22)  + norm_FT*nik(1)*dirFT(2)
          bins(i,23)  = bins(i,23)  + norm_FT*nik(1)*dirFT(3)
          bins(i,24)  = bins(i,24)  + norm_FT*nik(2)*dirFT(1)
          bins(i,25)  = bins(i,25)  + norm_FT*nik(2)*dirFT(2)
          bins(i,26)  = bins(i,26)  + norm_FT*nik(2)*dirFT(3)
          bins(i,27)  = bins(i,27)  + norm_FT*nik(3)*dirFT(1)
          bins(i,28)  = bins(i,28)  + norm_FT*nik(3)*dirFT(2)
          bins(i,29)  = bins(i,29)  + norm_FT*nik(3)*dirFT(3)
        end if
      end if
    end do
  end do

  ! Normalizing by the number of contacts
  do i=1,n_bin
    bins(i,3:29) = bins(i,3:29) / valid_c!bins(i,2)
  end do

  ! The full force tensor
  bins(:,21) = bins(:,21) + bins(:,12)
  bins(:,22) = bins(:,22) + bins(:,13)
  bins(:,23) = bins(:,23) + bins(:,14)
  bins(:,24) = bins(:,24) + bins(:,15)
  bins(:,25) = bins(:,25) + bins(:,16)
  bins(:,26) = bins(:,26) + bins(:,17)
  bins(:,27) = bins(:,27) + bins(:,18)
  bins(:,28) = bins(:,28) + bins(:,19)
  bins(:,29) = bins(:,29) + bins(:,20)

  print*, bins

  ! Writing the heading
  write(123,*) '#   Class     ', '      ac      ', ' 2(X3-1)/(X1+3)n', ' 2(X3-1)/(X1+3)f'

  ! Finding the eigen values for each bin
  do i=1, n_bin
    ! The fabric
    Fabric(1,1) = bins(i,3)
    Fabric(1,2) = bins(i,4)
    Fabric(1,3) = bins(i,5)
    Fabric(2,1) = bins(i,6)
    Fabric(2,2) = bins(i,7)
    Fabric(2,3) = bins(i,8)
    Fabric(3,1) = bins(i,9)
    Fabric(3,2) = bins(i,10)
    Fabric(3,3) = bins(i,11)

    lda  = 3
    matz = 1

    call rg (lda, 3, Fabric, wr, wi, matz, localframe, ierror)

    if ( wr(1)==wr(2) .and. (wr(2)==wr(3)) ) then
      S1=wr(1)
      S2=S1
      S3=S1
    else
      S3 = max( wr(1),max(wr(2),wr(3)))
      S1 = min( wr(1),min(wr(2),wr(3)))
      if (wr(1)==wr(2)) then
        if (S1==wr(1)) S2=S1
        if (S3==wr(1)) S2=S3
      else if (wr(3)==wr(2)) then
        if (S1==wr(2)) S2=S1
        if (S3==wr(3)) S2=S3
      else
        do j=1,3
          if ((wr(j)<S3) .and. (wr(j)>S1)) S2=wr(j)
        end do
      end if
    end if

    ! the contact anisotropy
    ac = 2*(S3-S1)/(S3+S1)

    ! Computing the fabric forces tensors
    Fabric_N(1,1) = bins(i,12)
    Fabric_N(1,2) = bins(i,13)
    Fabric_N(1,3) = bins(i,14)
    Fabric_N(2,1) = bins(i,15)
    Fabric_N(2,2) = bins(i,16)
    Fabric_N(2,3) = bins(i,17)
    Fabric_N(3,1) = bins(i,18)
    Fabric_N(3,2) = bins(i,19)
    Fabric_N(3,3) = bins(i,20)

    Fabric_T(1,1) = bins(i,21)
    Fabric_T(1,2) = bins(i,22)
    Fabric_T(1,3) = bins(i,23)
    Fabric_T(2,1) = bins(i,24)
    Fabric_T(2,2) = bins(i,25)
    Fabric_T(2,3) = bins(i,26)
    Fabric_T(3,1) = bins(i,27)
    Fabric_T(3,2) = bins(i,28)
    Fabric_T(3,3) = bins(i,29)


    ! Finding the eigen values of the fabric normal forces tensor
    lda = 3
    matz = 1

    call rg(lda, 3, Fabric_N, wrn, win, matz, localframeN, ierror)

    if ( wrn(1)==wrn(2) .and. (wrn(2)==wrn(3)) ) then
      X1n=wrn(1)
      X2n=X1n
      X3n=X1n
    else
      X3n = max( wrn(1),max(wrn(2),wrn(3)))
      X1n = min( wrn(1),min(wrn(2),wrn(3)))
      if (wrn(1)==wrn(2)) then
        if (X1n==wrn(1)) X2n=X1n
        if (X3n==wrn(1)) X2n=X3n
      else if (wrn(3)==wrn(2)) then
        if (X1n==wrn(2)) X2n=X1n
        if (X3n==wrn(3)) X2n=X3n
      else
        do j=1,3
          if ((wrn(j)<X3n) .and. (wrn(j)>X1n)) X2n=wrn(j)
        end do
      end if
    end if

    ! Finding the eigen values of the full fabric forces tensor
    call rg(lda, 3, Fabric_T, wrf, wif, matz, localframeF, ierror)

    if ( wrf(1)==wrf(2) .and. (wrf(2)==wrf(3)) ) then
      X1f=wrf(1)
      X2f=X1f
      X3f=X1f
    else
      X3f = max(wrf(1),max(wrf(2),wrf(3)))
      X1f = min(wrf(1),min(wrf(2),wrf(3)))
      if (wrf(1)==wrf(2)) then
        if (X1f==wrf(1)) X2f=X1f
        if (X3f==wrf(1)) X2f=X3f
      else if (wrf(3)==wrf(2)) then
        if (X1f==wrf(2)) X2f=X1f
        if (X3f==wrf(3)) X2f=X3f
      else
        do j=1,3
          if ((wrf(j)<X3f) .and. (wrf(j)>X1f)) X2f=wrf(j)
        end do
      end if
    end if

    write(123,'(4(1X,E12.5))') bins(i,1), ac, 2*(X3n-X1n)/(X1n+X3n), 2*(X3f-X1f)/(X1f+X3f)
  end do

  close(123)

  print*, 'Write Anisotropy class    ---> Ok!'

end subroutine aniso_class


!==============================================================================
! Computing anisotropies per class of particles by branch
!==============================================================================
subroutine aniso_part

  implicit none

  integer                                  :: i,j,cd,an, n_bin
  integer                                  :: ierror,matz,lda, valid_c
  real*8                                   :: S1, S2, S3, cpt, ac, interval, norm_FN, norm_FT
  real*8                                   :: b_max, b_min, Rsik, Rnik, Rtik, norm_lik
  real*8                                   :: X1f, X1n, X2f, X2n, X3f, X3n
  real*8,dimension(3)                      :: nik, Fik, sik, tik, FikT, dirFT, Lik
  real*8,dimension(3)                      :: wr,wi, wrf, wrn, win, wif
  real*8,dimension(3,3)                    :: localframe, localframeN, localframeF
  real*8,dimension(:,:), allocatable       :: bins
  real*8,dimension(3,3)                    :: Fabric, Fabric_N, Fabric_T
  real*8,dimension(:,:), allocatable       :: all_fabric, all_forces_n, all_forces_t
  character(len=24)                        :: aniso_c_file

  real*8                                   :: pos_z_w1, pos_z_w2

  pos_z_w1 = TAB_PLAN(1)%center(3)
  pos_z_w2 = TAB_PLAN(2)%center(3)

  ! The name of the file
  aniso_c_file = './POSTPRO/AN_CLA_P.     '

  ! Preparing the corresponding number of the file
  if (compteur_clout<10) then
    WRITE(aniso_c_file(20:21),'(I1)')   compteur_clout
  else if ( (compteur_clout>=10) .and. (compteur_clout<100) ) then
    WRITE(aniso_c_file(20:22),'(I2)')   compteur_clout
  else if ( (compteur_clout>=100).and. (compteur_clout<1000) ) then
    WRITE(aniso_c_file(20:23),'(I3)')   compteur_clout
  else if ( (compteur_clout>=1000).and. (compteur_clout<10000) ) then
    WRITE(aniso_c_file(20:24),'(I4)')   compteur_clout
  else
    print*, "Not enough file names :: Anisotropies per class"
    stop
  end if

  ! Opening the file
  open(unit=124,file=aniso_c_file,status='replace')

  ! Initializing
  b_max = -9999.
  b_min =  9999.

  ! Finding the max and min particle size
  do i=1, n_particles
    ! Counting a smaller box letting a space to the walls
    if (TAB_SPHERE(i)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
    if (TAB_SPHERE(i)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle

    b_min = min(b_min, TAB_SPHERE(i)%Rmax)
    b_max = max(b_max, TAB_SPHERE(i)%Rmax)
  end do

  ! The number of bins - > Intervals in the size distribution
  n_bin = 36

  ! The size of the interval
  interval = (b_max - b_min)/n_bin

  ! Allocating the bins
  if (allocated(bins)) deallocate(bins)
  allocate(bins(n_bin,29))
  ! initializing
  bins(:,:) = 0

  ! The array of fabrics
  if (allocated(all_fabric)) deallocate(all_fabric)
  allocate(all_fabric(n_particles,9))
  ! initializing
  all_fabric(:,:) = 0

  ! The array of forces tensor normal
  if (allocated(all_forces_n)) deallocate(all_forces_n)
  allocate(all_forces_n(n_particles,9))
  ! initializing
  all_forces_n(:,:) = 0

  ! The array of forces tensor tangential
  if (allocated(all_forces_t)) deallocate(all_forces_t)
  allocate(all_forces_t(n_particles,9))
  ! initializing
  all_forces_t(:,:) = 0

  ! Valid number of contacts
  valid_c = 0

  ! Computing the fabric individually for all particles
  do i=1, nb_ligneCONTACT
    if(TAB_CONTACT(i)%nature == 'SPPLx') cycle
    cd = TAB_CONTACT(i)%icdent
    an = TAB_CONTACT(i)%ianent
    if ((cd .gt. n_particles) .or. (an .gt. n_particles)) cycle

    ! Counting a smaller box letting a space to the walls
    if (TAB_SPHERE(cd)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
    if (TAB_SPHERE(cd)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle

    if (TAB_SPHERE(an)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
    if (TAB_SPHERE(an)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle

    Rnik = TAB_CONTACT(i)%rn

    ! Checking if it is an active contact
    if (Rnik .le. 0.D0) cycle

    valid_c = valid_c + 1

    nik  = TAB_CONTACT(i)%n
    tik  = TAB_CONTACT(i)%t
    sik  = TAB_CONTACT(i)%s

    Rtik = TAB_CONTACT(i)%rt
    Rsik = TAB_CONTACT(i)%rs

    ! The force vector
    Fik(:) = Rnik*nik(:) + Rtik*tik(:) + Rsik*sik(:)

    ! The norm of the normal force
    norm_FN = Fik(1)*nik(1) + Fik(2)*nik(2) + Fik(3)*nik(3)

    ! The tangential force vector
    FikT(1) = Fik(1) - Rnik*nik(1)
    FikT(2) = Fik(2) - Rnik*nik(2)
    FikT(3) = Fik(3) - Rnik*nik(3)

    ! The norm of FikT
    norm_FT = sqrt(FikT(1)**2 + FikT(2)**2 + FikT(3)**2)

    all_fabric(cd,1) = all_fabric(cd,1) + nik(1)*nik(1)
    all_fabric(cd,2) = all_fabric(cd,2) + nik(1)*nik(2)
    all_fabric(cd,3) = all_fabric(cd,3) + nik(1)*nik(3)
    all_fabric(cd,4) = all_fabric(cd,4) + nik(2)*nik(1)
    all_fabric(cd,5) = all_fabric(cd,5) + nik(2)*nik(2)
    all_fabric(cd,6) = all_fabric(cd,6) + nik(2)*nik(3)
    all_fabric(cd,7) = all_fabric(cd,7) + nik(3)*nik(1)
    all_fabric(cd,8) = all_fabric(cd,8) + nik(3)*nik(2)
    all_fabric(cd,9) = all_fabric(cd,9) + nik(3)*nik(3)

    all_fabric(an,1) = all_fabric(an,1) + nik(1)*nik(1)
    all_fabric(an,2) = all_fabric(an,2) + nik(1)*nik(2)
    all_fabric(an,3) = all_fabric(an,3) + nik(1)*nik(3)
    all_fabric(an,4) = all_fabric(an,4) + nik(2)*nik(1)
    all_fabric(an,5) = all_fabric(an,5) + nik(2)*nik(2)
    all_fabric(an,6) = all_fabric(an,6) + nik(2)*nik(3)
    all_fabric(an,7) = all_fabric(an,7) + nik(3)*nik(1)
    all_fabric(an,8) = all_fabric(an,8) + nik(3)*nik(2)
    all_fabric(an,9) = all_fabric(an,9) + nik(3)*nik(3)


    all_forces_n(cd,1) = all_forces_n(cd,1) + norm_FN*nik(1)*nik(1)
    all_forces_n(cd,2) = all_forces_n(cd,2) + norm_FN*nik(1)*nik(2)
    all_forces_n(cd,3) = all_forces_n(cd,3) + norm_FN*nik(1)*nik(3)
    all_forces_n(cd,4) = all_forces_n(cd,4) + norm_FN*nik(2)*nik(1)
    all_forces_n(cd,5) = all_forces_n(cd,5) + norm_FN*nik(2)*nik(2)
    all_forces_n(cd,6) = all_forces_n(cd,6) + norm_FN*nik(2)*nik(3)
    all_forces_n(cd,7) = all_forces_n(cd,7) + norm_FN*nik(3)*nik(1)
    all_forces_n(cd,8) = all_forces_n(cd,8) + norm_FN*nik(3)*nik(2)
    all_forces_n(cd,9) = all_forces_n(cd,9) + norm_FN*nik(3)*nik(3)

    all_forces_n(an,1) = all_forces_n(an,1) + norm_FN*nik(1)*nik(1)
    all_forces_n(an,2) = all_forces_n(an,2) + norm_FN*nik(1)*nik(2)
    all_forces_n(an,3) = all_forces_n(an,3) + norm_FN*nik(1)*nik(3)
    all_forces_n(an,4) = all_forces_n(an,4) + norm_FN*nik(2)*nik(1)
    all_forces_n(an,5) = all_forces_n(an,5) + norm_FN*nik(2)*nik(2)
    all_forces_n(an,6) = all_forces_n(an,6) + norm_FN*nik(2)*nik(3)
    all_forces_n(an,7) = all_forces_n(an,7) + norm_FN*nik(3)*nik(1)
    all_forces_n(an,8) = all_forces_n(an,8) + norm_FN*nik(3)*nik(2)
    all_forces_n(an,9) = all_forces_n(an,9) + norm_FN*nik(3)*nik(3)

    if (norm_FT .gt. 0.0) then

      ! The direction of FikT
      dirFT = FikT/norm_FT

      all_forces_t(cd,1) = all_forces_t(cd,1) + norm_FT*nik(1)*dirFT(1)
      all_forces_t(cd,2) = all_forces_t(cd,2) + norm_FT*nik(1)*dirFT(2)
      all_forces_t(cd,3) = all_forces_t(cd,3) + norm_FT*nik(1)*dirFT(3)
      all_forces_t(cd,4) = all_forces_t(cd,4) + norm_FT*nik(2)*dirFT(1)
      all_forces_t(cd,5) = all_forces_t(cd,5) + norm_FT*nik(2)*dirFT(2)
      all_forces_t(cd,6) = all_forces_t(cd,6) + norm_FT*nik(2)*dirFT(3)
      all_forces_t(cd,7) = all_forces_t(cd,7) + norm_FT*nik(3)*dirFT(1)
      all_forces_t(cd,8) = all_forces_t(cd,8) + norm_FT*nik(3)*dirFT(2)
      all_forces_t(cd,9) = all_forces_t(cd,9) + norm_FT*nik(3)*dirFT(3)

      all_forces_t(an,1) = all_forces_t(an,1) + norm_FT*nik(1)*dirFT(1)
      all_forces_t(an,2) = all_forces_t(an,2) + norm_FT*nik(1)*dirFT(2)
      all_forces_t(an,3) = all_forces_t(an,3) + norm_FT*nik(1)*dirFT(3)
      all_forces_t(an,4) = all_forces_t(an,4) + norm_FT*nik(2)*dirFT(1)
      all_forces_t(an,5) = all_forces_t(an,5) + norm_FT*nik(2)*dirFT(2)
      all_forces_t(an,6) = all_forces_t(an,6) + norm_FT*nik(2)*dirFT(3)
      all_forces_t(an,7) = all_forces_t(an,7) + norm_FT*nik(3)*dirFT(1)
      all_forces_t(an,8) = all_forces_t(an,8) + norm_FT*nik(3)*dirFT(2)
      all_forces_t(an,9) = all_forces_t(an,9) + norm_FT*nik(3)*dirFT(3)
    end if
  end do

  do i=1, n_bin
    bins(i,1) = b_min + interval*(i - 0.5)
  end do

 ! Computing the branch and putting the corresponding stress tensor in the corresponding bin
  do i=1, n_bin
    do j=1, n_particles

      ! Counting a smaller box letting a space to the walls
      if (TAB_SPHERE(j)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
      if (TAB_SPHERE(j)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle

      if((TAB_SPHERE(j)%Rmax .ge. (b_min + (i-1)*interval)) .and. (TAB_SPHERE(j)%Rmax*0.999 .lt. (interval*i + b_min))) then

        ! The counter of particles
        bins(i,2)  = bins(i,2)  + 1

        ! The fabric tensor
        bins(i,3)  = bins(i,3)  + all_fabric(j,1)
        bins(i,4)  = bins(i,4)  + all_fabric(j,2)
        bins(i,5)  = bins(i,5)  + all_fabric(j,3)
        bins(i,6)  = bins(i,6)  + all_fabric(j,4)
        bins(i,7)  = bins(i,7)  + all_fabric(j,5)
        bins(i,8)  = bins(i,8)  + all_fabric(j,6)
        bins(i,9)  = bins(i,9)  + all_fabric(j,7)
        bins(i,10) = bins(i,10) + all_fabric(j,8)
        bins(i,11) = bins(i,11) + all_fabric(j,9)

        ! The chi tensor of normal forces
        bins(i,12)  = bins(i,12)  + all_forces_n(j,1)
        bins(i,13)  = bins(i,13)  + all_forces_n(j,2)
        bins(i,14)  = bins(i,14)  + all_forces_n(j,3)
        bins(i,15)  = bins(i,15)  + all_forces_n(j,4)
        bins(i,16)  = bins(i,16)  + all_forces_n(j,5)
        bins(i,17)  = bins(i,17)  + all_forces_n(j,6)
        bins(i,18)  = bins(i,18)  + all_forces_n(j,7)
        bins(i,19)  = bins(i,19)  + all_forces_n(j,8)
        bins(i,20)  = bins(i,20)  + all_forces_n(j,9)

        if (norm_FT .gt. 0.0) then

          ! The direction of FikT
          dirFT = FikT/norm_FT

          ! The chi tensor of tangential forces
          bins(i,21)  = bins(i,21)  + all_forces_t(j,1)
          bins(i,22)  = bins(i,22)  + all_forces_t(j,2)
          bins(i,23)  = bins(i,23)  + all_forces_t(j,3)
          bins(i,24)  = bins(i,24)  + all_forces_t(j,4)
          bins(i,25)  = bins(i,25)  + all_forces_t(j,5)
          bins(i,26)  = bins(i,26)  + all_forces_t(j,6)
          bins(i,27)  = bins(i,27)  + all_forces_t(j,7)
          bins(i,28)  = bins(i,28)  + all_forces_t(j,8)
          bins(i,29)  = bins(i,29)  + all_forces_t(j,9)
        end if
      end if
    end do
  end do

  ! Normalizing by the total number of valid contacts
  do i=1,n_bin
    bins(i,3:29) = bins(i,3:29) / (2*valid_c) !bins(i,2)
  end do

  ! The full force tensor
  bins(:,21) = bins(:,21) + bins(:,12)
  bins(:,22) = bins(:,22) + bins(:,13)
  bins(:,23) = bins(:,23) + bins(:,14)
  bins(:,24) = bins(:,24) + bins(:,15)
  bins(:,25) = bins(:,25) + bins(:,16)
  bins(:,26) = bins(:,26) + bins(:,17)
  bins(:,27) = bins(:,27) + bins(:,18)
  bins(:,28) = bins(:,28) + bins(:,19)
  bins(:,29) = bins(:,29) + bins(:,20)

  print*, bins

  ! Writing the heading
  write(124,*) '#   Class     ', '      ac      ', ' 2(X3-1)/(X1+3)n', ' 2(X3-1)/(X1+3)f'

  ! Finding the eigen values for each bin
  do i=1, n_bin
    ! The fabric
    Fabric(1,1) = bins(i,3)
    Fabric(1,2) = bins(i,4)
    Fabric(1,3) = bins(i,5)
    Fabric(2,1) = bins(i,6)
    Fabric(2,2) = bins(i,7)
    Fabric(2,3) = bins(i,8)
    Fabric(3,1) = bins(i,9)
    Fabric(3,2) = bins(i,10)
    Fabric(3,3) = bins(i,11)

    ! Computing the fabric forces tensors
    Fabric_N(1,1) = bins(i,12)
    Fabric_N(1,2) = bins(i,13)
    Fabric_N(1,3) = bins(i,14)
    Fabric_N(2,1) = bins(i,15)
    Fabric_N(2,2) = bins(i,16)
    Fabric_N(2,3) = bins(i,17)
    Fabric_N(3,1) = bins(i,18)
    Fabric_N(3,2) = bins(i,19)
    Fabric_N(3,3) = bins(i,20)

    Fabric_T(1,1) = bins(i,21)
    Fabric_T(1,2) = bins(i,22)
    Fabric_T(1,3) = bins(i,23)
    Fabric_T(2,1) = bins(i,24)
    Fabric_T(2,2) = bins(i,25)
    Fabric_T(2,3) = bins(i,26)
    Fabric_T(3,1) = bins(i,27)
    Fabric_T(3,2) = bins(i,28)
    Fabric_T(3,3) = bins(i,29)

    lda  = 3
    matz = 1

    call rg (lda, 3, Fabric, wr, wi, matz, localframe, ierror)

    if ( wr(1)==wr(2) .and. (wr(2)==wr(3)) ) then
      S1=wr(1)
      S2=S1
      S3=S1
    else
      S3 = max( wr(1),max(wr(2),wr(3)))
      S1 = min( wr(1),min(wr(2),wr(3)))
      if (wr(1)==wr(2)) then
        if (S1==wr(1)) S2=S1
        if (S3==wr(1)) S2=S3
      else if (wr(3)==wr(2)) then
        if (S1==wr(2)) S2=S1
        if (S3==wr(3)) S2=S3
      else
        do j=1,3
          if ((wr(j)<S3) .and. (wr(j)>S1)) S2=wr(j)
        end do
      end if
    end if

    ! the contact anisotropy
    ac = 2*(S3-S1)/(S3+S1)


    ! Finding the eigen values of the fabric normal forces tensor
    lda = 3
    matz = 1

    call rg(lda, 3, Fabric_N, wrn, win, matz, localframeN, ierror)

    if ( wrn(1)==wrn(2) .and. (wrn(2)==wrn(3)) ) then
      X1n=wrn(1)
      X2n=X1n
      X3n=X1n
    else
      X3n = max( wrn(1),max(wrn(2),wrn(3)))
      X1n = min( wrn(1),min(wrn(2),wrn(3)))
      if (wrn(1)==wrn(2)) then
        if (X1n==wrn(1)) X2n=X1n
        if (X3n==wrn(1)) X2n=X3n
      else if (wrn(3)==wrn(2)) then
        if (X1n==wrn(2)) X2n=X1n
        if (X3n==wrn(3)) X2n=X3n
      else
        do j=1,3
          if ((wrn(j)<X3n) .and. (wrn(j)>X1n)) X2n=wrn(j)
        end do
      end if
    end if

    ! Finding the eigen values of the full fabric forces tensor
    call rg(lda, 3, Fabric_T, wrf, wif, matz, localframeF, ierror)

    if ( wrf(1)==wrf(2) .and. (wrf(2)==wrf(3)) ) then
      X1f=wrf(1)
      X2f=X1f
      X3f=X1f
    else
      X3f = max(wrf(1),max(wrf(2),wrf(3)))
      X1f = min(wrf(1),min(wrf(2),wrf(3)))
      if (wrf(1)==wrf(2)) then
        if (X1f==wrf(1)) X2f=X1f
        if (X3f==wrf(1)) X2f=X3f
      else if (wrf(3)==wrf(2)) then
        if (X1f==wrf(2)) X2f=X1f
        if (X3f==wrf(3)) X2f=X3f
      else
        do j=1,3
          if ((wrf(j)<X3f) .and. (wrf(j)>X1f)) X2f=wrf(j)
        end do
      end if
    end if

    write(124,'(4(1X,E12.5))') bins(i,1), ac, 2*(X3n-X1n)/(X1n+X3n), 2*(X3f-X1f)/(X1f+X3f)
  end do

  close(124)

  print*, 'Write Anisotropy class part ---> Ok!'

end subroutine aniso_part


!==============================================================================
! Contacts orientation
!==============================================================================
subroutine ctc_dir

  implicit none

  real(kind=8)                                 :: interval
  real(kind=8)                                 :: angcurr, maxint, minint
  integer                                      :: i, j, ninter, contotal, cd, an
  real(kind=8), dimension(:,:), allocatable    :: theta_interval
  real(kind=8), dimension(3)                   :: nik
  character(len=26)                            :: file_name
  logical                                      :: dir_ctcdir

  real*8                                       :: pos_z_w1, pos_z_w2

  pos_z_w1 = TAB_PLAN(1)%center(3)
  pos_z_w2 = TAB_PLAN(2)%center(3)

  ! Cleaning or creating the folder if necessary
  if (first_over_all) then
    ! Asking if the file already exists
    inquire(file='./POSTPRO/CTCDIR', exist=dir_ctcdir)
    if(dir_ctcdir) then
      ! Cleaning
      call system('rm ./POSTPRO/CTCDIR/*')
    else
      ! Creating
      call system('mkdir ./POSTPRO/CTCDIR')
    end if
  end if

  ! The file name
  file_name       =  './POSTPRO/CTCDIR/CDIR.0000'

  if (compteur_clout<10) then
    WRITE(file_name(25:26),'(I1)')   compteur_clout
  else if ( (compteur_clout>=10) .and. (compteur_clout<100) ) then
    WRITE(file_name(24:26),'(I2)')   compteur_clout
  else if ( (compteur_clout>=100).and. (compteur_clout<1000) ) then
    WRITE(file_name(23:26),'(I3)')   compteur_clout
  end if
  ! Opening the file
  open(unit=120,file=file_name,status='replace')

  ! Number of intervals over pi
  ninter = 36

  ! Initializing variables
  interval = pi/ninter

  ! Allocating the vector that will contain the contact direction and frequency
  if (allocated(theta_interval)) deallocate(theta_interval)
  allocate(theta_interval(ninter,3))

  ! Computing the contact direction and initializing the frequency
  do i=1, ninter
    theta_interval(i,1) = (interval/2)*(2*i-1)
    theta_interval(i,2) = 0
  end do

  ! Computing the frequency for each interval
  do i=1, ninter
    ! The width of the interval
    minint = interval*(i - 1)
    maxint = interval*i
    do j=1,nb_ligneCONTACT
      cd = TAB_CONTACT(j)%icdent
      an = TAB_CONTACT(j)%ianent
      ! Only between discs
      if(TAB_CONTACT(j)%nature == 'SPPLx') cycle
      ! Only active contacts
      if(TAB_CONTACT(j)%rn .le. 0.D0) cycle
      if (TAB_SPHERE(cd)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
      if (TAB_SPHERE(cd)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle
      if (TAB_SPHERE(an)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
      if (TAB_SPHERE(an)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle

      ! The normal vector
      nik(1)  = TAB_CONTACT(j)%n(1)
      nik(2)  = TAB_CONTACT(j)%n(2)
      nik(3)  = TAB_CONTACT(j)%n(3)

      ! In the plane xy!
      ! The 1 and 2 quadrants only
      ! Changing the normal vector direction if necesary
      if(nik(3) .lt. 0) then
        nik(1)  = -nik(1)
        nik(3)  = -nik(3)
      end if

      ! Finding the current angle. Note the vector nik is already unitary!
      angcurr = acos(nik(1)/sqrt(nik(1)**2 + nik(3)**2))

      ! Checking the interval
      if(angcurr .lt. maxint .and. angcurr .ge. minint) then
        theta_interval(i,2) = theta_interval(i,2) + 1
      end if
    end do
  end do

  contotal = 0

  ! Computing the number of active contacts
  do i=1, ninter
    contotal = contotal + theta_interval(i,2)
  end do

  ! Normalizing
  theta_interval(:,3) = theta_interval(:,2)/(2*contotal)!*interval)????

  ! Writing the heading
  write(120,*) '# Theta(Rad) ', '  Freq_norm  '

  ! Writing
  do i=1, ninter
    write(120,'(3(1X,E12.5))') theta_interval(i,1), theta_interval(i,2), theta_interval(i,3)
  end do

  ! Writing the 3rd and 4th quadrants
  do i=1, ninter
    write(120,'(3(1X,E12.5))') theta_interval(i,1)+pi, theta_interval(i,2), theta_interval(i,3)
  end do

  close(120)

  print*, 'Write Contacts dir    ---> Ok!'

end subroutine ctc_dir

!==============================================================================
! Branch orientation
!==============================================================================
subroutine branch_dir

  implicit none

  integer                                      :: i, j, ninter, contotal, cd, an
  real(kind=8)                                 :: interval, pro_n, pro_t
  real(kind=8)                                 :: angcurr, maxint, minint
  real(kind=8), dimension(:,:), allocatable    :: theta_interval
  real(kind=8), dimension(3)                   :: nik, Lik, tik, sik, tanik
  character(len=26)                            :: file_name
  logical                                      :: dir_ctcdir

  real*8                                       :: pos_z_w1, pos_z_w2

  pos_z_w1 = TAB_PLAN(1)%center(3)
  pos_z_w2 = TAB_PLAN(2)%center(3)

  ! Cleaning or creating the folder if necessary
  if (first_over_all) then
    ! Asking if the file already exists
    inquire(file='./POSTPRO/BRCDIR', exist=dir_ctcdir)
    if(dir_ctcdir) then
      ! Cleaning
      call system('rm ./POSTPRO/BRCDIR/*')
    else
      ! Creating
      call system('mkdir ./POSTPRO/BRCDIR')
    end if
  end if

  ! The file name
  file_name       =  './POSTPRO/BRCDIR/BDIR.0000'

  if (compteur_clout<10) then
    WRITE(file_name(25:26),'(I1)')   compteur_clout
  else if ( (compteur_clout>=10) .and. (compteur_clout<100) ) then
    WRITE(file_name(24:26),'(I2)')   compteur_clout
  else if ( (compteur_clout>=100).and. (compteur_clout<1000) ) then
    WRITE(file_name(23:26),'(I3)')   compteur_clout
  end if
  ! Opening the file
  open(unit=121,file=file_name,status='replace')

  ! Number of intervals over pi
  ninter = 36

  ! Initializing variables
  contotal = 0
  interval = pi/ninter

  ! Allocating the vector that will contain the contact direction and frequency
  if (allocated(theta_interval)) deallocate(theta_interval)
  allocate(theta_interval(ninter,4))

  ! Computing the contact direction, initializing the frequency, and the total branch sum
  do i=1, ninter
    theta_interval(i,1) = (interval/2)*(2*i-1)
    theta_interval(i,2) = 0
    theta_interval(i,3) = 0
    theta_interval(i,4) = 0
  end do

  ! Computing the frequency for each interval
  do i=1, ninter
    ! The width of the interval
    minint = interval*(i - 1)
    maxint = interval*i
    do j=1,nb_ligneCONTACT
      cd = TAB_CONTACT(j)%icdent
      an = TAB_CONTACT(j)%ianent
      ! Only between discs
      if(TAB_CONTACT(j)%nature == 'SPPLx') cycle
      ! Only active contacts
      if(TAB_CONTACT(j)%rn .le. 0.D0) cycle
      if (TAB_SPHERE(cd)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
      if (TAB_SPHERE(cd)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle
      if (TAB_SPHERE(an)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
      if (TAB_SPHERE(an)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle

      ! The branch vector
      Lik(:) = TAB_SPHERE(cd)%center(:)-TAB_SPHERE(an)%center(:)

      ! Not joining contacts in the periodic zones
      if (sqrt(Lik(1)**2 + Lik(2)**2 + Lik(3)**2) .gt. (TAB_SPHERE(cd)%Rmax + TAB_SPHERE(an)%Rmax)) cycle

      ! The normal vector
      nik(1)  = TAB_CONTACT(j)%n(1)
      nik(2)  = TAB_CONTACT(j)%n(2)
      nik(3)  = TAB_CONTACT(j)%n(3)

      ! The tangential vector - xz plane!
      !tik = TAB_CONTACT(j)%t
      !sik = TAB_CONTACT(j)%s
      !tanik = sik*TAB_CONTACT(j)%rs + tik*TAB_CONTACT(j)%rt
      tanik(1) = nik(3)
      tanik(2) = 0.
      tanik(3) = -nik(1)

      ! unitary
      !tanik = tanik/sqrt(tanik(1)**2 + tanik(2)**2 + tanik(3)**2)

      ! Normal projection
      pro_n = abs(Lik(1)*nik(1) + Lik(2)*nik(2) + Lik(3)*nik(3))

      ! Tangential projection
      pro_t = abs(Lik(1)*tanik(1) + Lik(2)*tanik(2) + Lik(3)*tanik(3))

      ! The 1 and 2 quadrants only
      ! Changing the vector direction if necessary
      if(nik(3) .lt. 0) then
        nik(1) = -nik(1)
        nik(3) = -nik(3)
      end if

      ! Finding the current angle. Note the vector nik is already unitary!
      angcurr = acos(nik(1)/sqrt(nik(1)**2 + nik(3)**2))

      ! Checking the interval
      if(angcurr .lt. maxint .and. angcurr .ge. minint) then
        theta_interval(i,2) = theta_interval(i,2) + 1
        theta_interval(i,3) = theta_interval(i,3) + pro_n
        theta_interval(i,4) = theta_interval(i,4) + pro_t
      end if
    end do
  end do

  ! Computing the number of active contacts
  do i=1, ninter
    contotal = contotal + theta_interval(i,2)
  end do

  ! The averages
  theta_interval(:,3) = theta_interval(:,3)/(theta_interval(:,2)*2)
  theta_interval(:,4) = theta_interval(:,4)/(theta_interval(:,2)*2)

  ! Writing the heading
  write(121,*) '# Theta(Rad) ', '    <ln>    ', '    <lt>    '

  ! Writing
  do i=1, ninter
    write(121,'(3(1X,E12.5))') theta_interval(i,1), theta_interval(i,3), theta_interval(i,4)
  end do

  ! Writing the 3rd and 4th quadrants
  do i=1, ninter
    write(121,'(3(1X,E12.5))') theta_interval(i,1)+pi, theta_interval(i,3), theta_interval(i,4)
  end do

  close(121)

  print*, 'Write Branch dir    ---> Ok!'

end subroutine branch_dir


!==============================================================================
! Contacts orientation
!==============================================================================
subroutine forces_dir

  implicit none

  integer                                      :: i, j, ninter, contotal, cd, an
  real(kind=8)                                 :: interval, pro_n, pro_t, Rnik, Rtik, Rsik
  real(kind=8)                                 :: angcurr, maxint, minint
  real(kind=8), dimension(:,:), allocatable    :: theta_interval
  real(kind=8), dimension(3)                   :: nik, Fik, tik, sik, tanik, tanik_plane
  character(len=26)                            :: file_name
  logical                                      :: dir_ctcdir

  real*8                                       :: pos_z_w1, pos_z_w2

  pos_z_w1 = TAB_PLAN(1)%center(3)
  pos_z_w2 = TAB_PLAN(2)%center(3)

  ! Cleaning or creating the folder if necessary
  if (first_over_all) then
    ! Asking if the file already exists
    inquire(file='./POSTPRO/FRCDIR', exist=dir_ctcdir)
    if(dir_ctcdir) then
      ! Cleaning
      call system('rm ./POSTPRO/FRCDIR/*')
    else
      ! Creating
      call system('mkdir ./POSTPRO/FRCDIR')
    end if
  end if

  ! The file name
  file_name       =  './POSTPRO/FRCDIR/FDIR.0000'

  if (compteur_clout<10) then
    WRITE(file_name(25:26),'(I1)')   compteur_clout
  else if ( (compteur_clout>=10) .and. (compteur_clout<100) ) then
    WRITE(file_name(24:26),'(I2)')   compteur_clout
  else if ( (compteur_clout>=100).and. (compteur_clout<1000) ) then
    WRITE(file_name(23:26),'(I3)')   compteur_clout
  end if
  ! Opening the file
  open(unit=122,file=file_name,status='replace')

  ! Number of intervals over pi
  ninter = 36

  ! Initializing variables
  contotal = 0
  interval = pi/ninter

  ! Allocating the vector that will contain the contact direction and frequency
  if (allocated(theta_interval)) deallocate(theta_interval)
  allocate(theta_interval(ninter,4))

  ! Computing the contact direction, initializing the frequency, and the total force sum
  do i=1, ninter
    theta_interval(i,1) = (interval/2)*(2*i-1)
    theta_interval(i,2) = 0
    theta_interval(i,3) = 0
    theta_interval(i,4) = 0
  end do

  ! Computing the frequency for each interval
  do i=1, ninter
    ! The width of the interval
    minint = interval*(i - 1)
    maxint = interval*i
    do j=1,nb_ligneCONTACT
      cd = TAB_CONTACT(j)%icdent
      an = TAB_CONTACT(j)%ianent
      ! Only between discs
      if(TAB_CONTACT(j)%nature == 'SPPLx') cycle
      ! Only active contacts
      if(TAB_CONTACT(j)%rn .le. 0.D0) cycle
      if (TAB_SPHERE(cd)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
      if (TAB_SPHERE(cd)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle
      if (TAB_SPHERE(an)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
      if (TAB_SPHERE(an)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle

      Rnik = TAB_CONTACT(j)%rn
      Rtik = TAB_CONTACT(j)%rt
      Rsik = TAB_CONTACT(j)%rs

      ! The normal vector
      nik = TAB_CONTACT(j)%n

      ! The tangential vector
      tik = TAB_CONTACT(j)%t
      sik = TAB_CONTACT(j)%s
      ! The real tangential vector - in xz plane!
      !tanik(1:3) = Rsik*sik(1:3) + Rtik*tik(1:3)
      tanik(1) = nik(3)
      tanik(2) = 0.
      tanik(3) = -nik(1)

      ! The force vector
      Fik(1:3) = (Rnik*nik(1:3)+Rtik*tik(1:3)+Rsik*sik(1:3))

      ! unitary!
      !tanik = tanik/sqrt(tanik(1)**2 + tanik(2)**2 + tanik(3)**2)

      ! Normal projection
      pro_n = abs(Fik(1)*nik(1) + Fik(2)*nik(2) + Fik(3)*nik(3))

      !print*, 'ini'
      !print*, pro_n
      !print*, TAB_CONTACT(j)%rn
      !print*, TAB_CONTACT(j)%n
      !print*, Fik
      !print*, 'asdasd'
      !tanik_plane(1) = Rsik*sik(1) + Rtik*tik(1)
      !tanik_plane(2) = 0
      !tanik_plane(3) = Rsik*sik(3) + Rtik*tik(3)

      !tanik_plane = tanik_plane/sqrt(tanik_plane(1)**2+tanik_plane(3)**2)

      ! Tangential projection in the xz plane!!!
      !pro_t = abs(Fik(1)*tanik(1) + Fik(2)*tanik(2) + Fik(3)*tanik(3))
      pro_t = Fik(1)*tanik(1) + Fik(2)*tanik(2) + Fik(3)*tanik(3)
      !print*, 'ini'
      !print*, pro_t
      !print*, Rtik
      !print*, Rsik
      !print*, sqrt(Rsik**2+Rtik**2)
      !print*, 'asdasd'

      ! The 1 and 2 quadrants only
      ! Changing the vector direction if necessary
      if(nik(3) .lt. 0) then
        nik(1) = -nik(1)
        nik(3) = -nik(3)
      end if

      ! Finding the current angle. Note the vector nik is already unitary!
      angcurr = acos(nik(1)/sqrt(nik(1)**2 + nik(3)**2))

      ! Checking the interval
      if(angcurr .lt. maxint .and. angcurr .ge. minint) then
        theta_interval(i,2) = theta_interval(i,2) + 1
        theta_interval(i,3) = theta_interval(i,3) + pro_n
        theta_interval(i,4) = theta_interval(i,4) + pro_t
      end if
    end do
  end do

  ! Computing the number of active contacts
  !do i=1, ninter
  !  contotal = contotal + theta_interval(i,2)
  !end do

  ! The averages
  theta_interval(:,3) = theta_interval(:,3)/(theta_interval(:,2)*2)
  theta_interval(:,4) = theta_interval(:,4)/(theta_interval(:,2)*2)

  ! Normalizing by the normal contact force
  !theta_interval(:,4) = theta_interval(:,4)/(theta_interval(:,3))

  ! Writing the heading
  write(122,*) '# Theta(Rad) ', '    <fn>    ', '    <ft>    '

  ! Writing
  do i=1, ninter
    write(122,'(3(1X,E12.5))') theta_interval(i,1), theta_interval(i,3), theta_interval(i,4)
  end do

  ! Writing the 3rd and 4th quadrants
  do i=1, ninter
    write(122,'(3(1X,E12.5))') theta_interval(i,1)+pi, theta_interval(i,3), theta_interval(i,4)
  end do

  close(122)

  print*, 'Write Forces dir    ---> Ok!'

end subroutine forces_dir

!==============================================================================
! Computing the friction mobilization
!==============================================================================
subroutine mobilization

  implicit none

  integer                                  :: i, j, cd, an, n_bin, disp_wall_times
  real*8                                   :: mob_ave, act_ctc
  character(len=24)                        :: coor_c_file

  real*8                                   :: pos_z_w1, pos_z_w2

  pos_z_w1 = TAB_PLAN(1)%center(3)
  pos_z_w2 = TAB_PLAN(2)%center(3)

  mob_ave = 0.
  act_ctc = 0.
  disp_wall_times = 5

  do i=1, nb_ligneCONTACT
    if (TAB_CONTACT(i)%nature /= 'SPSPx') cycle
    cd = TAB_CONTACT(i)%icdent
    an = TAB_CONTACT(i)%ianent

    if (cd .gt. n_particles .or. an .gt. n_particles) cycle

    ! Far from the walls
    if (TAB_SPHERE(cd)%center(3) .lt. pos_z_w1 + 0.01*disp_wall_times) cycle
    if (TAB_SPHERE(cd)%center(3) .gt. pos_z_w2 - 0.01*disp_wall_times) cycle

    if (TAB_SPHERE(an)%center(3) .lt. pos_z_w1 + 0.01*disp_wall_times) cycle
    if (TAB_SPHERE(an)%center(3) .gt. pos_z_w2 - 0.01*disp_wall_times) cycle

    if (TAB_CONTACT(i)%rn .le. 0) cycle

    mob_ave = mob_ave + sqrt(TAB_CONTACT(i)%rs**2 + TAB_CONTACT(i)%rt**2)/(0.4*TAB_CONTACT(i)%rn)
    act_ctc = act_ctc + 1
  end do

  ! The average
  mob_ave = mob_ave/act_ctc

  if (first_over_all) then
    ! Writing
    write(126,*) '#   Time    ', ' <|ft|/mu fn>   '
  end if

  write(126,'(2(1X,E12.5))') time, mob_ave

  print*, 'Write mobilization     ---> Ok!'

end subroutine mobilization

!==============================================================================
! Computing the radial dist 
!==============================================================================
subroutine radial_dist

  implicit none

  integer                                  :: i, j, k, n_bins
  real*8                                   :: mean_rad, dist_part
  real*8, dimension(:,:), allocatable      :: array_dist
  character(len=24)                        :: coor_c_file


  real*8                                   :: pos_z_w1, pos_z_w2

  pos_z_w1 = TAB_PLAN(1)%center(3)
  pos_z_w2 = TAB_PLAN(2)%center(3)

  ! The name of the file
  coor_c_file = './POSTPRO/RAD_DIST.     '

  ! Preparing the corresponding number of the file
  if (compteur_clout<10) then
    WRITE(coor_c_file(20:21),'(I1)')   compteur_clout
  else if ( (compteur_clout>=10) .and. (compteur_clout<100) ) then
    WRITE(coor_c_file(20:22),'(I2)')   compteur_clout
  else if ( (compteur_clout>=100).and. (compteur_clout<1000) ) then
    WRITE(coor_c_file(20:23),'(I3)')   compteur_clout
  else if ( (compteur_clout>=1000).and. (compteur_clout<10000) ) then
    WRITE(coor_c_file(20:24),'(I4)')   compteur_clout
  else
    print*, "Not enough file names :: Radial dist"
    stop
  end if

  ! Opening the file
  open(unit=127,file=coor_c_file,status='replace')

  n_bins = 100
  mean_rad = 0.01

  if (allocated(array_dist)) deallocate(array_dist)
  allocate(array_dist(n_bins,2))

  array_dist(:,:) = 0.0

  ! Filling the distance columns
  do i=1, n_bins
    array_dist(i,1) = (10.*mean_rad/(n_bins))*i
  end do

  ! Filling the function g
  do k=1, n_bins
    print*, k
    do i=1, n_particles-1
      do j=i+1, n_particles
        dist_part = (TAB_SPHERE(i)%center(1) - TAB_SPHERE(j)%center(1))**2 + &
                    (TAB_SPHERE(i)%center(2) - TAB_SPHERE(j)%center(2))**2 + &
                    (TAB_SPHERE(i)%center(3) - TAB_SPHERE(j)%center(3))**2
        dist_part = dist_part**0.5
        if (dist_part .gt. 10.*mean_rad) cycle
        if (dist_part .lt. array_dist(k+1,1) .and. dist_part .ge. array_dist(k,1)) then
          array_dist(k,2) = array_dist(k,2) + 1.
        end if
      end do
    end do
  end do

  write(127,*) '#      r     ', '     N(r)    '

  do i=1, n_bins
    if (array_dist(i,2) .gt. 0) then
      write(127,'(2(1X,E12.5))') array_dist(i,1), array_dist(i,2)
    end if
  end do

  close(127)

  print*, 'Write radial dist     ---> Ok!'

end subroutine radial_dist

!==============================================================================
! Computing a segregation index
!==============================================================================

subroutine segregation_index

  implicit none

  integer                               :: i, nbl_ctc, counter_l
  real*8                                :: vl_total
  real*8                                :: x_mean, y_mean, z_mean
  real*8                                :: x_mean_vol, y_mean_vol, z_mean_vol
  real*8                                :: x_sdev, y_sdev, z_sdev
  real*8                                :: x_sdev_vol, y_sdev_vol, z_sdev_vol
  real*8                                :: xl_min, xl_max, yl_min, yl_max, zl_min, zl_max
  real*8, dimension(3)                  :: geo_center
  real*8, dimension(:,:), allocatable   :: correla_radius
  real*8                                :: r1_mean, r2_mean, var_corre_1, var_corre_2, var_corre_3
  real*8                                :: crr

  ! Initializing
  x_mean = 0.
  y_mean = 0.
  z_mean = 0.

  x_mean_vol = 0.
  y_mean_vol = 0.
  z_mean_vol = 0.

  ! Finding the geometrical center of the box
  xl_max =-10000000
  xl_min = 10000000
  yl_max =-10000000
  yl_min = 10000000
  zl_max =-10000000
  zl_min = 10000000

  do i=1,n_particles

    xl_max = max(xl_max, TAB_SPHERE(i)%center(1) + TAB_SPHERE(i)%Rmax)
    xl_min = min(xl_min, TAB_SPHERE(i)%center(1) - TAB_SPHERE(i)%Rmax)

    yl_max = max(yl_max, TAB_SPHERE(i)%center(2) + TAB_SPHERE(i)%Rmax)
    yl_min = min(yl_min, TAB_SPHERE(i)%center(2) - TAB_SPHERE(i)%Rmax)

    zl_max = max(zl_max, TAB_SPHERE(i)%center(3) + TAB_SPHERE(i)%Rmax)
    zl_min = min(zl_min, TAB_SPHERE(i)%center(3) - TAB_SPHERE(i)%Rmax)
  end do

  geo_center(1) = (xl_max+xl_min)/2
  geo_center(2) = (yl_max+yl_min)/2
  geo_center(3) = (zl_max+zl_min)/2

  ! Computing the total volume of grains
  vl_total = 0
  do i=1, n_particles
    vl_total = vl_total + TAB_SPHERE(i)%volume
  end do

  ! Computing a simple mean of coordinates
  do i=1, n_particles
    x_mean = x_mean + TAB_SPHERE(i)%center(1)
    y_mean = y_mean + TAB_SPHERE(i)%center(2)
    z_mean = z_mean + TAB_SPHERE(i)%center(3)
  end do

  x_mean = x_mean/n_particles
  y_mean = y_mean/n_particles
  z_mean = z_mean/n_particles

  ! Computing the standard deviation
  x_sdev = 0.
  y_sdev = 0.
  z_sdev = 0.

  do i=1, n_particles
    x_sdev = x_sdev + (TAB_SPHERE(i)%center(1) - x_mean)**2
    y_sdev = y_sdev + (TAB_SPHERE(i)%center(2) - y_mean)**2
    z_sdev = z_sdev + (TAB_SPHERE(i)%center(3) - z_mean)**2
  end do

  x_sdev = (x_sdev/(n_particles-1))**0.5
  y_sdev = (y_sdev/(n_particles-1))**0.5
  z_sdev = (z_sdev/(n_particles-1))**0.5

  ! Computing an average coordinate weigthed by volume
  do i=1, n_particles
    x_mean_vol = x_mean_vol + TAB_SPHERE(i)%center(1)*TAB_SPHERE(i)%volume
    y_mean_vol = y_mean_vol + TAB_SPHERE(i)%center(2)*TAB_SPHERE(i)%volume
    z_mean_vol = z_mean_vol + TAB_SPHERE(i)%center(3)*TAB_SPHERE(i)%volume
  end do

  x_mean_vol = x_mean_vol/(vl_total)
  y_mean_vol = y_mean_vol/(vl_total)
  z_mean_vol = z_mean_vol/(vl_total)

  ! Computing the standard deviation between geometrical center and the weigthed coordinates
  x_sdev_vol= 0.
  y_sdev_vol= 0.
  z_sdev_vol= 0.

  do i=1, n_particles
    x_sdev_vol = x_sdev_vol + (TAB_SPHERE(i)%center(1)*TAB_SPHERE(i)%volume - x_mean)**2
    y_sdev_vol = y_sdev_vol + (TAB_SPHERE(i)%center(2)*TAB_SPHERE(i)%volume - y_mean)**2
    z_sdev_vol = z_sdev_vol + (TAB_SPHERE(i)%center(3)*TAB_SPHERE(i)%volume - z_mean)**2
  end do

  x_sdev_vol = (x_sdev_vol/((n_particles-1)*vl_total))**0.5
  y_sdev_vol = (y_sdev_vol/((n_particles-1)*vl_total))**0.5
  z_sdev_vol = (z_sdev_vol/((n_particles-1)*vl_total))**0.5

  ! Looking for valid contact
  nbl_ctc = 0
  do i=1, nb_ligneCONTACT
    if (TAB_CONTACT(i)%nature /= 'SPSPx') cycle
    if (TAB_CONTACT(i)%rn .gt. 0) then
      nbl_ctc = nbl_ctc + 1
    end if
  end do

  ! Computing the correlation between particle sizes at contact
  if (allocated(correla_radius)) deallocate(correla_radius)
  allocate(correla_radius(nbl_ctc,2))

  counter_l = 0
  ! Collecting the radius
  do i=1, nb_ligneCONTACT
    if (TAB_CONTACT(i)%nature /= 'SPSPx') cycle
    if (TAB_CONTACT(i)%rn .gt. 0) then
        counter_l = counter_l + 1
      correla_radius(counter_l,1) = TAB_SPHERE(TAB_CONTACT(i)%icdent)%Rmax
      correla_radius(counter_l,2) = TAB_SPHERE(TAB_CONTACT(i)%ianent)%Rmax
    end if
  end do

  r1_mean = 0.
  r2_mean = 0.

  ! The means
  do i=1, nbl_ctc
    r1_mean = r1_mean + correla_radius(i,1)
    r2_mean = r2_mean + correla_radius(i,2)
  end do

  r1_mean = r1_mean/nbl_ctc
  r2_mean = r2_mean/nbl_ctc

  ! Computing the terms of the correlation
  var_corre_1 = 0.
  var_corre_3 = 0.
  var_corre_2 = 0.

  do i=1, nbl_ctc
    var_corre_1 = var_corre_1 + (correla_radius(i,1) - r1_mean)*(correla_radius(i,2) - r2_mean)
    var_corre_2 = var_corre_2 + (correla_radius(i,1) - r1_mean)**2
    var_corre_3 = var_corre_3 + (correla_radius(i,2) - r2_mean)**2
  end do

  ! Correlation index
  crr = var_corre_1/((var_corre_2)**0.5 * (var_corre_3)**0.5)

  if (first_over_all) then
    write(128,*) '#    Time    ', '     <x>     ', '     <y>     ', '     <z>     ', &
                                  '     S(x)    ', '     S(y)    ', '     S(z)    ', &
                                  '    <x>_v    ', '    <y>_v    ', '    <z>_v    ', &
                                  '   S(x)_v    ', '   S(y)_v    ', '   S(z)_v    ', &
                                  '     Crr     '
  end if

  write(128,'(14(1X,E12.5))') time, x_mean,     y_mean,     z_mean,     &
                                    x_sdev,     y_sdev,     z_sdev,     &
                                    x_mean_vol, y_mean_vol, z_mean_vol, &
                                    x_sdev_vol, y_sdev_vol, z_sdev_vol, crr

  print*, 'Write segregation index ---> Ok!'

end subroutine segregation_index

!==============================================================================
! Computing the variance of particle positions weigthed by mass
! for lacey index
!==============================================================================

subroutine for_lacey

  implicit none

  integer                         ::  i
  real*8                          ::  l_variance, total_vol
  real*8, dimension(3)            ::  average_coor, curr_dif_vec


  average_coor(:) = 0.0
  total_vol = 0.0
  ! Computing the average particle positions weigthed by volume
  do i=1, n_particles
    average_coor(:) = average_coor(:) + TAB_SPHERE(i)%volume*TAB_SPHERE(i)%center(:)
    total_vol = total_vol + TAB_SPHERE(i)%volume
  end do

  average_coor(:) = average_coor(:)/total_vol

  ! Finding the variance of the distances of particles centers
  ! to the average particle center
  do i=1, n_particles
    curr_dif_vec(:) = TAB_SPHERE(i)%volume*TAB_SPHERE(i)%center(:) - average_coor(:)
    l_variance = l_variance + (curr_dif_vec(1)**2 + curr_dif_vec(2)**2+curr_dif_vec(3)**2)
  end do

  l_variance = l_variance/n_particles

  if (first_over_all) then
    write(129,*) '#    Time    ', '      S      '
  end if

  write(129,'(2(1X,E12.5))') time, l_variance

  print*, 'Write for Lacey index ---> Ok!'

end subroutine for_lacey


!==============================================================================
! List of active particles
!==============================================================================
subroutine list_active_par

  implicit none

  integer                                  :: i, cd, an
  character(len=24)                        :: l_file

  real*8                                   :: pos_z_w1, pos_z_w2

  pos_z_w1 = TAB_PLAN(1)%center(3)
  pos_z_w2 = TAB_PLAN(2)%center(3)

  ! The name of the file
  l_file = './POSTPRO/LIST_ACT.     '

  ! Preparing the corresponding number of the file
  if (compteur_clout<10) then
    WRITE(l_file(20:21),'(I1)')   compteur_clout
  else if ( (compteur_clout>=10) .and. (compteur_clout<100) ) then
    WRITE(l_file(20:22),'(I2)')   compteur_clout
  else if ( (compteur_clout>=100).and. (compteur_clout<1000) ) then
    WRITE(l_file(20:23),'(I3)')   compteur_clout
  else if ( (compteur_clout>=1000).and. (compteur_clout<10000) ) then
    WRITE(l_file(20:24),'(I4)')   compteur_clout
  else
    print*, "Not enough file names :: List of active particles"
    stop
  end if

  ! Opening the file
  open(unit=130,file=l_file,status='replace')

  write(130,*) '# Active diameters '

  do i=1,n_particles
    if (TAB_SPHERE(i)%nctc .ge. 2) then
      write(130, '(E12.5)') TAB_SPHERE(i)%Rmax
    end if
  end do

  close(130)

  print*, 'Write Active Part     ---> Ok!'

end subroutine list_active_par


!==============================================================================
! Computing the proportion of floating particles per class
!==============================================================================
subroutine float_class

  implicit none

  integer                                  :: i, j, cd, an, n_bin, active_par
  real*8                                   :: r_max, r_min, interval
  real*8, dimension(:,:), allocatable      :: bins
  integer, dimension(:), allocatable       :: valid_part
  character(len=27)                        :: coor_c_file

  real*8                                   :: pos_z_w1, pos_z_w2

  pos_z_w1 = TAB_PLAN(1)%center(3)
  pos_z_w2 = TAB_PLAN(2)%center(3)

  ! The name of the file
  coor_c_file = './POSTPRO/FLOAT_CLASS.     '

  ! Preparing the corresponding number of the file
  if (compteur_clout<10) then
    WRITE(coor_c_file(23:24),'(I1)')   compteur_clout
  else if ( (compteur_clout>=10) .and. (compteur_clout<100) ) then
    WRITE(coor_c_file(23:25),'(I2)')   compteur_clout
  else if ( (compteur_clout>=100).and. (compteur_clout<1000) ) then
    WRITE(coor_c_file(23:26),'(I3)')   compteur_clout
  else if ( (compteur_clout>=1000).and. (compteur_clout<10000) ) then
    WRITE(coor_c_file(23:27),'(I4)')   compteur_clout
  else
    print*, "Not enough file names :: Floating per class"
    stop
  end if

  ! Opening the file
  open(unit=131,file=coor_c_file,status='replace')

  ! Initalizing
  r_max = -9999.
  r_min =  9999.
  active_par = 0

  ! Finding the max and min radius
  do i=1, n_particles
    ! Particles far from the walls
    if (TAB_SPHERE(i)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
    if (TAB_SPHERE(i)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle

    active_par = active_par + 1

    r_min = min(r_min, TAB_SPHERE(i)%Rmax)
    r_max = max(r_max, TAB_SPHERE(i)%Rmax)
  end do

  ! The number of bins - > Intervals in the size distribution
  n_bin = 36

  ! The size of the interval
  interval = (r_max - r_min)/n_bin

  ! Allocating the bins
  if (allocated(bins)) deallocate(bins)
  allocate(bins(n_bin,3))
  ! initializing
  bins(:,:) = 0

  ! Adding particles to bins as function of their size. First column: average size.
  ! Second: number of particles

  ! Caracteristic size of each bin
  do i=1, n_bin
    bins(i,1) = r_min + interval*(i - 0.5)
  end do

  do i=1, n_bin
    do j=1, n_particles
      if (TAB_SPHERE(j)%nctc .gt. 1) cycle

      if((TAB_SPHERE(j)%Rmax .ge. (r_min + (i-1)*interval)) .and. (TAB_SPHERE(j)%Rmax*0.999 .lt. (interval*i + r_min))) then 

        if (TAB_SPHERE(j)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
        if (TAB_SPHERE(j)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle

        bins(i,2) = bins(i,2) + 1.0
      end if
    end do
  end do

  bins(:,3) = bins(:,2)/real(active_par)

  ! Writing
  write(131,*) '#   Class     ', '    n_float    ', '    prop_float   '

  do i=1, n_bin
    if (bins(i,2) .gt. 0.) then
      write(131,'(3(1X,E12.5))') bins(i,1), bins(i,2), bins(i,3)
    end if
  end do

  close(131)

  print*, 'Float Class           ---> Ok!'

end subroutine float_class


!==============================================================================
! Computing the proportion of particles with a given number of contacts
!==============================================================================
subroutine prop_nctc

  implicit none

  integer                                  :: i
  real*8                                   :: l_p0, l_p1, l_p2, l_p3, l_p4, &
                                              l_p5, l_p6, l_p7, l_p8, l_p9, &
                                              l_p10, l_p11, l_p12, l_p13, l_p14, &
                                              l_p15, l_p16

  l_p0 = 0.
  l_p1 = 0.
  l_p2 = 0.
  l_p3 = 0.
  l_p4 = 0.
  l_p5 = 0.
  l_p6 = 0.
  l_p7 = 0.
  l_p8 = 0.
  l_p9 = 0.
  l_p10 = 0.
  l_p11 = 0.
  l_p12 = 0.
  l_p13 = 0.
  l_p14 = 0.
  l_p15 = 0.
  l_p16 = 0.

  do i=1, n_particles
    if (TAB_SPHERE(i)%nctc == 0 ) then
      l_p0 = l_p0 + 1.
    else if (TAB_SPHERE(i)%nctc == 1) then
      l_p1 = l_p1 + 1.
    else if (TAB_SPHERE(i)%nctc == 2) then
      l_p2 = l_p2 + 1.
    else if (TAB_SPHERE(i)%nctc == 3) then
      l_p3 = l_p3 + 1.
    else if (TAB_SPHERE(i)%nctc == 4) then
      l_p4 = l_p4 + 1.
    else if (TAB_SPHERE(i)%nctc == 5) then
      l_p5 = l_p5 + 1.
    else if (TAB_SPHERE(i)%nctc == 6) then
      l_p6 = l_p6 + 1.
    else if (TAB_SPHERE(i)%nctc == 7) then
      l_p7 = l_p7 + 1.
    else if (TAB_SPHERE(i)%nctc == 8) then
      l_p8 = l_p8 + 1.
    else if (TAB_SPHERE(i)%nctc == 9) then
      l_p9 = l_p9 + 1.
    else if (TAB_SPHERE(i)%nctc == 10) then
      l_p10 = l_p10 + 1.
    else if (TAB_SPHERE(i)%nctc == 11) then
      l_p11 = l_p11 + 1.
    else if (TAB_SPHERE(i)%nctc == 12) then
      l_p12 = l_p12 + 1.
    else if (TAB_SPHERE(i)%nctc == 13) then
      l_p13 = l_p13 + 1.
    else if (TAB_SPHERE(i)%nctc == 14) then
      l_p14 = l_p14 + 1.
    else if (TAB_SPHERE(i)%nctc == 15) then
      l_p15 = l_p15 + 1.
    else if (TAB_SPHERE(i)%nctc == 16) then
      l_p16 = l_p16 + 1.
    end if
  end do

  l_p0 = l_p0/real(n_particles)
  l_p1 = l_p1/real(n_particles)
  l_p2 = l_p2/real(n_particles)
  l_p3 = l_p3/real(n_particles)
  l_p4 = l_p4/real(n_particles)
  l_p5 = l_p5/real(n_particles)
  l_p6 = l_p6/real(n_particles)
  l_p7 = l_p7/real(n_particles)
  l_p8 = l_p8/real(n_particles)
  l_p9 = l_p9/real(n_particles)
  l_p10 = l_p10/real(n_particles)
  l_p11 = l_p11/real(n_particles)
  l_p12 = l_p12/real(n_particles)
  l_p13 = l_p13/real(n_particles)
  l_p14 = l_p14/real(n_particles)
  l_p15 = l_p15/real(n_particles)
  l_p16 = l_p16/real(n_particles)

  if (first_over_all) then
    ! Writing
    write(132,*) '#   Time    ', '     p0    ','     p1    ','     p2    ', &
                                 '     p3    ','     p4    ','     p5    ', &
                                 '     p6    ','     p7    ','     p8    ', &
                                 '     p9    ','    p10    ','     11    ', &
                                 '    p12    ','    p13    ','    p14    ', &
                                 '    p15    ','    p16    '
  end if

  write(132,'(18(1X,E12.5))') time, l_p0, l_p1, l_p2, l_p3, l_p4, l_p5, l_p6, l_p7, &
                                    l_p8, l_p9, l_p10, l_p11, l_p12, l_p13, l_p14, l_p15, l_p16

  print*, 'Write prop n_contacts ---> Ok!'

end subroutine prop_nctc


!================================================
! Computing the mean compressive force per particle class - Mstick
!================================================
subroutine qoverp_class

  implicit none

  integer                                        :: i, j, k, cd, an, n_bin
  real(kind=8)                                   :: Np_total, interval
  real(kind=8)                                   :: max_rad, min_rad
  real(kind=8)                                   :: Rtik,Rnik,Rsik
  real(kind=8)                                   :: S1,S2,S3
  real(kind=8), dimension(3)                     :: nik,tik,sik,Lik,Fik, ctc_point, jc_vec, ic_vec
  real(kind=8), dimension(3,3)                   :: M
  real(kind=8), dimension(:,:), allocatable      :: bins, tensor_part, pandq_bins
  logical                                        :: check_interface
  character(len=24)                              :: stress_c_file
  ! For the computation of eigen values and vectors
  real(kind=8),dimension(3,3)                    :: localframe
  integer                                        :: ierror,matz,lda
  real(kind=8),dimension(3)                      :: wr,wi

  real*8                                   :: pos_z_w1, pos_z_w2

  pos_z_w1 = TAB_PLAN(1)%center(3)
  pos_z_w2 = TAB_PLAN(2)%center(3)

  ! The name of the file
  stress_c_file = './POSTPRO/QoPCLASS.     '

  ! Preparing the corresponding number of the file
  if (compteur_clout<10) then
    WRITE(stress_c_file(20:21),'(I1)')   compteur_clout
  else if ( (compteur_clout>=10) .and. (compteur_clout<100) ) then
    WRITE(stress_c_file(20:22),'(I2)')   compteur_clout
  else if ( (compteur_clout>=100).and. (compteur_clout<1000) ) then
    WRITE(stress_c_file(20:23),'(I3)')   compteur_clout
  else if ( (compteur_clout>=1000).and. (compteur_clout<10000) ) then
    WRITE(stress_c_file(20:24),'(I4)')   compteur_clout
  else
    print*, "Not enough file names :: Stresses per class"
    stop
  end if

  ! Opening the file
  open(unit=133,file=stress_c_file,status='replace')

  ! Initializing variables

  ! Creating array to store the number of non cohesive contacts per group and their volume
  if (allocated(tensor_part)) deallocate(tensor_part)
  allocate(tensor_part(n_particles,9))

  ! Initializing
  tensor_part(:,:) = 0.0

  ! For all contacts -> Building the average stress tensor per group
  do i=1, nb_ligneCONTACT

    cd   = TAB_CONTACT(i)%icdent
    an   = TAB_CONTACT(i)%ianent

    ! Only contacts with betwen particles
    if (TAB_CONTACT(i)%nature=='SPPLx') cycle

    if (TAB_SPHERE(cd)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
    if (TAB_SPHERE(cd)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle

    if (TAB_SPHERE(an)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
    if (TAB_SPHERE(an)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle

    ! Active contacts
    if (TAB_CONTACT(i)%rn .gt. 0.0) then
      ! However they cannot belong to the same group!!
      nik  = TAB_CONTACT(i)%n
      tik  = TAB_CONTACT(i)%t
      sik  = TAB_CONTACT(i)%s
      Rtik = TAB_CONTACT(i)%rt
      Rnik = TAB_CONTACT(i)%rn
      Rsik = TAB_CONTACT(i)%rs
      ctc_point = TAB_CONTACT(i)%coor_ctc

      ! The branch vector
      Lik(1:3) = TAB_SPHERE(cd)%center(1:3)-TAB_SPHERE(an)%center(1:3)

      if (sqrt(Lik(1)**2 + Lik(2)**2 + Lik(3)**2) > 3*(TAB_SPHERE(cd)%Rmax + TAB_SPHERE(an)%Rmax)) then
        !Rayon_max = sqrt((Lik(1)**2 + Lik(2)**2 + Lik(3)**2))
        !Lik(:) = Lik(:) / Rayon_max
        Lik(:) = abs(TAB_SPHERE(cd)%Rmax+TAB_SPHERE(an)%Rmax)*nik(:)
      end if

      ! The vector from the center of cd to the contact point
      !ic_vec(1:3) = ctc_point(1:3) - TAB_SPHERE(cd)%center(1:3)
      ic_vec(1:3) = -nik(:)*TAB_SPHERE(cd)%Rmax

      ! The vector from the center of an to the contact point
      !jc_vec(1:3) = ctc_point(1:3) - TAB_SPHERE(an)%center(1:3)
      jc_vec(1:3) = nik(:)*TAB_SPHERE(an)%Rmax

      ! The force vector
      Fik(1:3) = (Rnik*nik(1:3)+Rtik*tik(1:3)+Rsik*sik(1:3))

      ! Building the tensor for the candidate
      tensor_part(cd, 1) = tensor_part(cd,1) - Fik(1)*ic_vec(1)
      tensor_part(cd, 2) = tensor_part(cd,2) - Fik(1)*ic_vec(2)
      tensor_part(cd, 3) = tensor_part(cd,3) - Fik(1)*ic_vec(3)
      tensor_part(cd, 4) = tensor_part(cd,4) - Fik(2)*ic_vec(1)
      tensor_part(cd, 5) = tensor_part(cd,5) - Fik(2)*ic_vec(2)
      tensor_part(cd, 6) = tensor_part(cd,6) - Fik(2)*ic_vec(3)
      tensor_part(cd, 7) = tensor_part(cd,7) - Fik(3)*ic_vec(1)
      tensor_part(cd, 8) = tensor_part(cd,8) - Fik(3)*ic_vec(2)
      tensor_part(cd, 9) = tensor_part(cd,9) - Fik(3)*ic_vec(3)

      ! Building the tensor for the antagonist
      tensor_part(an, 1) = tensor_part(an,1) + Fik(1)*jc_vec(1)
      tensor_part(an, 2) = tensor_part(an,2) + Fik(1)*jc_vec(2)
      tensor_part(an, 3) = tensor_part(an,3) + Fik(1)*jc_vec(3)
      tensor_part(an, 4) = tensor_part(an,4) + Fik(2)*jc_vec(1)
      tensor_part(an, 5) = tensor_part(an,5) + Fik(2)*jc_vec(2)
      tensor_part(an, 6) = tensor_part(an,6) + Fik(2)*jc_vec(3)
      tensor_part(an, 7) = tensor_part(an,7) + Fik(3)*jc_vec(1)
      tensor_part(an, 8) = tensor_part(an,8) + Fik(3)*jc_vec(2)
      tensor_part(an, 9) = tensor_part(an,9) + Fik(3)*jc_vec(3)
    end if
  end do

  ! Dividing the moment tensor to get the stress tensor
  do i=1, n_particles
    tensor_part(i,:) = tensor_part(i,:)/TAB_SPHERE(i)%volume
  end do

  ! Finding the max and min sizes
  min_rad =  9999.
  max_rad = -9999.

  ! Changing from volume to equivalent radius and finding max and min
  do i=1, n_particles
    min_rad = min(min_rad, TAB_SPHERE(i)%Rmax)
    max_rad = max(max_rad, TAB_SPHERE(i)%Rmax)
  end do

  !print*, min_rad
  !print*, max_rad

  ! The number of bins - > Intervals in the size distribution
  n_bin = 36

  ! The size of the interval
  interval = (max_rad - min_rad)/n_bin

  ! Allocating the bins
  if (allocated(bins)) deallocate(bins)
  allocate(bins(n_bin,11))
  ! initializing
  bins(:,:) = 0

  ! Characteristic size of each bin
  do i=1, n_bin
    bins(i,1) = min_rad + interval*(i - 0.5)
  end do

  if (allocated(pandq_bins)) deallocate(pandq_bins)
  allocate(pandq_bins(n_bin,2))

  pandq_bins(:,:) = 0.0

  ! The number of particles per bin and the tensor
  do i=1, n_bin
    do j=1, n_particles
      if((TAB_SPHERE(j)%Rmax .ge. (min_rad + (i-1)*interval)) .and. (TAB_SPHERE(j)%Rmax .lt. (interval*i + min_rad))) then

        if (tensor_part(j,1) .gt. 0 .and. tensor_part(j,5) .gt. 0 .and. tensor_part(j,9) .gt. 0) then
          bins(i,2) = bins(i,2) + 1

          bins(i,3)  = bins(i,3)  + tensor_part(j,1)
          bins(i,4)  = bins(i,4)  + tensor_part(j,2)
          bins(i,5)  = bins(i,5)  + tensor_part(j,3)
          bins(i,6)  = bins(i,6)  + tensor_part(j,4)
          bins(i,7)  = bins(i,7)  + tensor_part(j,5)
          bins(i,8)  = bins(i,8)  + tensor_part(j,6)
          bins(i,9)  = bins(i,9)  + tensor_part(j,7)
          bins(i,10) = bins(i,10) + tensor_part(j,8)
          bins(i,11) = bins(i,11) + tensor_part(j,9)
        end if
      end if
    end do
  end do

  ! Finding P and Q for each bin
  do i=1, n_bin
    !print*, 'yooooooo'
    !print*, i
    lda  = 3
    matz = 1

    ! Building a 3x3 matrix for the current stress tensor
    M(1,1) = bins(i,3)
    M(1,2) = bins(i,4)
    M(1,3) = bins(i,5)
    M(2,1) = bins(i,6)
    M(2,2) = bins(i,7)
    M(2,3) = bins(i,8)
    M(3,1) = bins(i,9)
    M(3,2) = bins(i,10)
    M(3,3) = bins(i,11)

    call rg (lda, 3, M, wr, wi, matz, localframe, ierror)

    if ( wr(1)==wr(2) .and. (wr(2)==wr(3)) ) then
      S1=wr(1)
      S2=S1
      S3=S1
    else
      S3 = max( wr(1),max(wr(2),wr(3)))
      S1 = min( wr(1),min(wr(2),wr(3)))
      if (wr(1)==wr(2)) then
        if (S1==wr(1)) S2=S1
        if (S3==wr(1)) S2=S3
      else if (wr(3)==wr(2)) then
        if (S1==wr(2)) S2=S1
        if (S3==wr(3)) S2=S3
      else
        do j=1,3
          if ((wr(j)<S3) .and. (wr(j)>S1)) S2=wr(j)
        end do
      end if
    end if

    ! Computing P
    pandq_bins(i,1) = (S1+S2+S3)/3

    ! Computing Q
    pandq_bins(i,2) = (S3-S1)/2

    ! Writing
    if (i==1) then
      write(133,*) '#   Avg size  ', '   N-parti  ', '      p     ', '      q     ', '    q/p     '
    end if
    if (bins(i,2) .gt. 0.0 ) then
      write(133,'(5(E14.7,1X))') bins(i,1), bins(i,2), pandq_bins(i,1), pandq_bins(i,2), pandq_bins(i,2)/pandq_bins(i,1)
    end if
  end do

  ! Closing the file
  close(133)

  print*, 'Write qoverp per class  ---> Ok!'

end subroutine qoverp_class


!==============================================================================
! stresses
!==============================================================================
subroutine stresses

  implicit none

  integer                                  :: i, cd, an
  integer                                  :: ierror, matz, lda
  real*8                                   :: Rtik, Rnik, Rsik
  real*8                                   :: S1,S2,S3,Rayon_max,q,p,dmean
  real*8,dimension(3)                      :: nik,tik,sik,Lik,Fik
  real*8,dimension(3)                      :: wr,wi
  real*8,dimension(3,3)                    :: Moment,M
  real*8,dimension(3,3)                    :: localframe

  real*8                                   :: pos_z_w1, pos_z_w2


  ! Initalizing variables
  Moment(:,:)   = 0.D0
  Lik = 0.0
  Fik = 0.0

  pos_z_w1 = TAB_PLAN(1)%center(3)
  pos_z_w2 = TAB_PLAN(2)%center(3)

  do i=1,nb_ligneCONTACT
    ! Contacts between particles
    if (TAB_CONTACT(i)%nature == 'SPPLx') cycle

    ! Active contacts
    if (TAB_CONTACT(i)%rn .le. 0.0) cycle

    cd   = TAB_CONTACT(i)%icdent
    an   = TAB_CONTACT(i)%ianent

    if ((cd .gt. n_particles) .or. (an .gt. n_particles)) cycle

    ! Counting a smaller box letting a space to the walls
    if (TAB_SPHERE(cd)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
    if (TAB_SPHERE(cd)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle

    if (TAB_SPHERE(an)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
    if (TAB_SPHERE(an)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle

    nik  = TAB_CONTACT(i)%n
    tik  = TAB_CONTACT(i)%t
    sik  = TAB_CONTACT(i)%s
    Rtik = TAB_CONTACT(i)%rt
    Rnik = TAB_CONTACT(i)%rn
    Rsik = TAB_CONTACT(i)%rs

    ! The branch
    Lik(1:3) = TAB_SPHERE(cd)%center(1:3)-TAB_SPHERE(an)%center(1:3)

    if (sqrt(Lik(1)**2 + Lik(2)**2 + Lik(3)**2) > 3*(TAB_SPHERE(cd)%Rmax + TAB_SPHERE(an)%Rmax)) then
      !Rayon_max = sqrt((Lik(1)**2 + Lik(2)**2 + Lik(3)**2))
      !Lik(:) = Lik(:) / Rayon_max
      Lik(:) = abs(TAB_SPHERE(cd)%Rmax+TAB_SPHERE(an)%Rmax)*nik(:)
    end if

    ! The force
    Fik(1:3) = (Rnik*nik(1:3)+Rtik*tik(1:3)+Rsik*sik(1:3))

    Moment(1,1:3) = Fik(1)*Lik(1:3) + Moment(1,1:3)
    Moment(2,1:3) = Fik(2)*Lik(1:3) + Moment(2,1:3)
    Moment(3,1:3) = Fik(3)*Lik(1:3) + Moment(3,1:3)
  end do

  ! The stress tensor
  Moment = Moment / (large * width * height)

  if (first_over_all) then
    write(134,*) '#   time     ', '   S( 1,1 )  ', '   S( 1,2 )  ', '   S( 1,3 )  ', &
                                  '   S( 2,1 )  ', '   S( 2,2 )  ', '   S( 2,3 )  ', &
                                  '   S( 3,1 )  ', '   S( 3,2 )  ', '   S( 3,3 )  '
  end if

  write(134,'(27(1X,E12.5))') time, Moment(1,1), Moment(1,2), Moment(1,3), &
                                    Moment(2,1), Moment(2,2), Moment(2,3), &
                                    Moment(3,1), Moment(3,2), Moment(3,3)

  print*, 'Write stresses compon ---> Ok!'

end subroutine stresses


!==============================================================================
! Computing the velocity profile
!==============================================================================
subroutine vloc_profile

  implicit none

  integer                                  :: i, j
  real*8                                   :: l_zmin, l_zmax, l_height, z1, z2
  real*8,dimension(:,:),allocatable        :: c_profile
  real*8,dimension(n_particles)            :: counted_part
  character(len=27)                        :: vp_c_file

  real*8                                   :: pos_z_w1, pos_z_w2

  pos_z_w1 = TAB_PLAN(1)%center(3)
  pos_z_w2 = TAB_PLAN(2)%center(3)

  ! The name of the file
  vp_c_file = './POSTPRO/VEL_PROFILE.     '

  ! Preparing the corresponding number of the file
  if (compteur_clout<10) then
    WRITE(vp_c_file(23:24),'(I1)')   compteur_clout
  else if ( (compteur_clout>=10) .and. (compteur_clout<100) ) then
    WRITE(vp_c_file(23:25),'(I2)')   compteur_clout
  else if ( (compteur_clout>=100).and. (compteur_clout<1000) ) then
    WRITE(vp_c_file(23:26),'(I3)')   compteur_clout
  else if ( (compteur_clout>=1000).and. (compteur_clout<10000) ) then
    WRITE(vp_c_file(23:27),'(I4)')   compteur_clout
  else
    print*, "Not enough file names :: Velocity profile"
    stop
  end if

  ! Opening the file
  open(unit=135,file=vp_c_file,status='replace')

  ! Finding the limits of the box to probe
  l_zmin =  9999.9
  l_zmax = -9999.9

  do i=1, n_particles
    if (TAB_SPHERE(i)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
    if (TAB_SPHERE(i)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle

    l_zmin = min(l_zmin, TAB_SPHERE(i)%center(3))
    l_zmax = max(l_zmax, TAB_SPHERE(i)%center(3))
  end do

  l_height = l_zmax - l_zmin

  ! Initalizing
  if (allocated(c_profile)) deallocate(c_profile)
  allocate(c_profile(n_layers,3))

  c_profile(:,:) = 0.0
  counted_part(:) = -1.0

  ! Filing the first colum of the array
  ! i.e., the average height of each layer
  do i=1, n_layers
    c_profile(i,1) = (l_height/n_layers)*(i-1) + (l_height/(2*n_layers))
  end do

  !print*, c_profile(:,1)

  ! Finding the average velocity in the shear plane by layers
  ! and the number of particles by layer
  do i=1, n_layers
    !print*, 'laaaaaa'
    !print*, i
    ! Current limits
    z1 = (l_height/n_layers)*(i-1) + l_zmin
    z2 = (l_height/n_layers)*i + l_zmin

    do j=1, n_particles
      ! Particles far from the walls
      if (counted_part(j) .gt. 0) cycle
      if (TAB_SPHERE(j)%center(3) .lt. pos_z_w1 + 0.01*5.0) then
        counted_part(j)=1
        cycle
      end if
      if (TAB_SPHERE(j)%center(3) .gt. pos_z_w2 - 0.01*5.0) then
        counted_part(j)=1
        cycle
      end if

      if (TAB_SPHERE(j)%center(3) .ge. z1 .and. TAB_SPHERE(j)%center(3) .lt. z2) then
        c_profile(i,2) = c_profile(i,2) + TAB_SPHERE(j)%Vx
        c_profile(i,3) = c_profile(i,3) + 1
        counted_part(j) = 1
      end if
    end do
  end do

  ! Computing the average velocity
  c_profile(:,2) = c_profile(:,2)/c_profile(:,3)

  ! Writing
  write(135,*) '#  Avg height  ', '#  Avg veloc   ', ' # particles  '

  do i=1, n_layers
    if(c_profile(i,3) .lt. 1) cycle
    write(135,'(3(1X,E12.5))') c_profile(i,1), c_profile(i,2),c_profile(i,3)
  end do

  close(135)

  print*, 'Velocity profile      ---> Ok!'

end subroutine vloc_profile

!==============================================================================
! PDF by coordination number
!==============================================================================
subroutine pdf_z

  implicit none

  integer                                  :: i, j, z_max
  character(len=23), dimension(20)         :: vp_c_file

  real*8                                   :: pos_z_w1, pos_z_w2

  pos_z_w1 = TAB_PLAN(1)%center(3)
  pos_z_w2 = TAB_PLAN(2)%center(3)

  ! The name of the file
  vp_c_file(1)  = './POSTPRO/PDF__Z1.     '
  vp_c_file(2)  = './POSTPRO/PDF__Z2.     '
  vp_c_file(3)  = './POSTPRO/PDF__Z3.     '
  vp_c_file(4)  = './POSTPRO/PDF__Z4.     '
  vp_c_file(5)  = './POSTPRO/PDF__Z5.     '
  vp_c_file(6)  = './POSTPRO/PDF__Z6.     '
  vp_c_file(7)  = './POSTPRO/PDF__Z7.     '
  vp_c_file(8)  = './POSTPRO/PDF__Z8.     '
  vp_c_file(9)  = './POSTPRO/PDF__Z9.     '
  vp_c_file(10) = './POSTPRO/PDF_Z10.     '
  vp_c_file(11) = './POSTPRO/PDF_Z11.     '
  vp_c_file(12) = './POSTPRO/PDF_Z12.     '
  vp_c_file(13) = './POSTPRO/PDF_Z13.     '
  vp_c_file(14) = './POSTPRO/PDF_Z14.     '
  vp_c_file(15) = './POSTPRO/PDF_Z15.     '
  vp_c_file(16) = './POSTPRO/PDF_Z16.     '
  vp_c_file(17) = './POSTPRO/PDF_Z17.     '
  vp_c_file(18) = './POSTPRO/PDF_Z18.     '
  vp_c_file(19) = './POSTPRO/PDF_Z19.     '
  vp_c_file(20) = './POSTPRO/PDF_Z20.     '

  z_max = 0

  do i=1, n_particles
    z_max = int(max(TAB_SPHERE(i)%nctc,real(z_max)))
  end do

  ! Max 20 files
  z_max = min(z_max,20)

  do i=2, z_max
    if (compteur_clout<10) then
      WRITE(vp_c_file(i)(19:20),'(I1)')   compteur_clout
    else if ( (compteur_clout>=10) .and. (compteur_clout<100) ) then
      WRITE(vp_c_file(i)(19:21),'(I2)')   compteur_clout
    else if ( (compteur_clout>=100).and. (compteur_clout<1000) ) then
      WRITE(vp_c_file(i)(19:22),'(I3)')   compteur_clout
    else if ( (compteur_clout>=1000).and. (compteur_clout<10000) ) then
      WRITE(vp_c_file(i)(19:23),'(I4)')   compteur_clout
    else
      print*, "Not enough file names :: Velocity profile"
      stop
    end if
    ! Opening the file
    open(unit=136,file=vp_c_file(i),status='replace')
    do j=1, n_particles
      if (i==TAB_SPHERE(j)%nctc) then
        if (TAB_SPHERE(j)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
        if (TAB_SPHERE(j)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle
        write(136,'(1(1X,E12.5))') TAB_SPHERE(j)%av_fn
      end if
    end do
    close(136)
  end do

  print*, 'Pdf coordinatio       ---> Ok!'

end subroutine pdf_z


!==================================================================================================
!Procedure pour calcul des contacts persistant ...
!==================================================================================================

subroutine persiste_contact
  implicit none
  integer                           :: i,j,cd,an
  real(kind=8)                      :: cpt,cpt_total,el_diameter
  logical                           :: contact_trouve

  if ( (compteur_clout >= debut_persistance) ) then
    print*,' --> Persistance des contacts'

    do i=1,n_particles
      if (TAB_SPHERE(i)%behav/='PLEXx') cycle
      cpt = cpt +1
      el_diameter = el_diameter + 2*TAB_SPHERE(i)%Rmax
    end do

    cpt        = 0
    if (first_list) then
      nb_ligneCONTACT_first = nb_ligneCONTACT
      allocate(TAB_CONTACT_INITIAL(nb_ligneCONTACT))
      do i=1,nb_ligneCONTACT
        TAB_CONTACT_INITIAL(i)%icdent      = TAB_CONTACT(i)%icdent
        TAB_CONTACT_INITIAL(i)%ianent      = TAB_CONTACT(i)%ianent
        TAB_CONTACT_INITIAL(i)%type        = TAB_CONTACT(i)%type
        TAB_CONTACT_INITIAL(i)%n(1)        = TAB_CONTACT(i)%n(1)
        TAB_CONTACT_INITIAL(i)%n(2)        = TAB_CONTACT(i)%n(2)
        TAB_CONTACT_INITIAL(i)%n(3)        = TAB_CONTACT(i)%n(3)
        TAB_CONTACT_INITIAL(i)%coor_ctc(1) = TAB_CONTACT(i)%coor_ctc(1)
        TAB_CONTACT_INITIAL(i)%coor_ctc(2) = TAB_CONTACT(i)%coor_ctc(2)
        TAB_CONTACT_INITIAL(i)%coor_ctc(3) = TAB_CONTACT(i)%coor_ctc(3)
        TAB_CONTACT_INITIAL(i)%rn          = TAB_CONTACT(i)%rn
        TAB_CONTACT_INITIAL(i)%rt          = TAB_CONTACT(i)%rt
        TAB_CONTACT_INITIAL(i)%rs          = TAB_CONTACT(i)%rs
        TAB_CONTACT_INITIAL(i)%contact_perdu = .false.
      end do
      first_list = .false.
    end if
 
    do i=1,nb_ligneCONTACT_first
      cd   = TAB_CONTACT_INITIAL(i)%icdent
      an   = TAB_CONTACT_INITIAL(i)%ianent
      if (TAB_CONTACT_INITIAL(i)%rn<0.000000000001*Mean_total_Normal_force) cycle
      if ((TAB_SPHERE(cd)%behav/='PLEXx').or.(TAB_SPHERE(an)%behav/='PLEXx')) cycle
      cpt_total = cpt_total + 1.0
    end do
     
    do i=1,nb_ligneCONTACT_first
      cd   = TAB_CONTACT_INITIAL(i)%icdent
      an   = TAB_CONTACT_INITIAL(i)%ianent
      if (TAB_CONTACT_INITIAL(i)%rn<0.000000000001*Mean_total_Normal_force) cycle
      if ((TAB_SPHERE(cd)%behav/='PLEXx').or.(TAB_SPHERE(an)%behav/='PLEXx')) cycle
      if (TAB_CONTACT_INITIAL(i)%contact_perdu) cycle
      do j=1,nb_ligneCONTACT
        cd   = TAB_CONTACT(j)%icdent
        an   = TAB_CONTACT(j)%ianent
        if (TAB_CONTACT(j)%rn<0.000000000001*Mean_total_Normal_force) cycle
        if ((TAB_SPHERE(cd)%behav/='PLEXx').or.(TAB_SPHERE(an)%behav/='PLEXx')) cycle

        if ((TAB_CONTACT_INITIAL(i)%icdent == TAB_CONTACT(j)%icdent).and.&
           (TAB_CONTACT_INITIAL(i)%ianent == TAB_CONTACT(j)%ianent)) then
          cpt = cpt + 1 
          contact_trouve = .true.
          exit
        else
          contact_trouve = .false.
        end if

      end do
      if (.NOT.contact_trouve) TAB_CONTACT_INITIAL(i)%contact_perdu = .true.
    end do
    !write(134,'(6(1X,D14.7))') time,Haut,time/(el_diameter * sqrt(1000/pressure_in_sample)),cpt/cpt_total
    !print*, 0.05/Haut,cpt/cpt_total
  end if
end subroutine persiste_contact


!==============================================================================
! Calcul du profil moyen de vitesse
!==============================================================================
subroutine profils
implicit none
integer                       :: i,k,j
integer                       :: nb_ligne_tab
real(kind=8)                  :: dmoyen,ep,cpt_v,cpt_c,Rayon_max,vboite
!! POUR LE CALCUL DU PROFIL DE VITESSE !!
real(kind=8)                  :: Vmoyen,Vymoyen,Vzmoyen,Vxmoyen
real(kind=8)                  :: Vrmoyen,Vrxmoyen,Vrymoyen,Vrzmoyen
real(kind=8),dimension(:,:),allocatable :: Moyenne_vitesse_moyenne
real(kind=8),dimension(:,:),allocatable :: Moyenne_vitesse_moyenne_r
!! POUR LE CALCUL DU PROFIL DE CONTRAINTE !!
type T_moyenne
  real(kind=8),dimension(3,3)              :: sigma
end type
type(T_moyenne),dimension(:),allocatable   ::  Moyenne_contrainte_moyenne
real(kind=8)                               ::  Rtik,Rnik,Rsik
real(kind=8),dimension(3)                  ::  nik,sik,tik,Lik,Fik
integer                                    ::  ierror_sigma,matz_sigma,lda_sigma
integer                                    ::  cd,an
real(kind=8),dimension(3)                  ::  wr_sigma,wi_sigma
real(kind=8)                               ::  S1_sigma,S2_sigma,S3_sigma,S_P
real(kind=8),dimension(3,3)                ::  localframe_sigma
real(kind=8),dimension(3,3)                ::  M_sigma
!! POUR LE CALCUL DU PROFIL DE I !!
real(kind=8),dimension(:,:),allocatable    :: Moyenne_I_moyen
real(kind=8)                               :: Mean_pressure
!! POUR LE CALCUL DU PROFIL DE COMPACITY !!
logical                                 :: sphere_in_band=.false.,point_out=.false.
real(kind=8)                            :: center_band,xmax_box,xmin_box,ymax_box,ymin_box,zmax_box,zmin_box,&
                                           epsilonx,epsilony,epsilonz,test_x,test_y,test_z,Lx,Ly,Lz,&
                                           nb_point_boite,nb_point_sphere,&
                                           volume_sphere_in_band,Vol_sphere,x_echantillon_max,x_echantillon_min,&
                                           nb_sphere_in_band,which_case,y_echantillon_max,y_echantillon_min
real(kind=8),dimension(:,:),allocatable :: Moyenne_compacity_moyen
integer                                 :: ix,jy,kz

if (compteur_clout>=debut_profils) then

  print*,'.... Profils ....'

  dmoyen = 0  
  x_echantillon_max = 0
  x_echantillon_min = 100000000.
  y_echantillon_max = 0
  y_echantillon_min = 100000000.
  do i=1,n_particles
    dmoyen = TAB_SPHERE(i)%Rmax + dmoyen
  end do
  dmoyen = 2*dmoyen / real(n_particles,8)
  print*,'Mean diamater =',dmoyen
  ep = dmoyen

!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  if (the_first_time_profils) then
    nb_ligne_tab = 0
    do
      nb_ligne_tab  =   nb_ligne_tab + 1
      ep = ep + dmoyen/10
      !if (ep>Haut+dmoyen) exit
    end do
!! POUR LE CALCUL DU PROFIL DE VITESSE !!
    allocate(TAB_vitesse(fin-debut_profils+1))
    do i=1,fin-debut_profils+1
      allocate(TAB_vitesse(i)%tab_vmoyen(nb_ligne_tab+100,4))
      allocate(TAB_vitesse(i)%tab_vrmoyen(nb_ligne_tab+100,4))
      TAB_vitesse(i)%tab_vmoyen = 0
      TAB_vitesse(i)%tab_vrmoyen = 0
    end do
!! POUR LE CALCUL DU PROFIL DE CONTRAINTE !!
    allocate(TAB_contrainte(fin-debut_profils+1))
    do i=1,fin-debut_profils+1
      allocate(TAB_contrainte(i)%tab_sigma(nb_ligne_tab+100))
      do j=1,nb_ligne_tab
        TAB_contrainte(i)%tab_sigma(j)%sigma = 0
      end do
    end do
!! POUR LE CALCUL DU PROFIL DE COMPACITY !!
    allocate(TAB_compacity(fin-debut_profils+1))
    do i=1,fin-debut_profils+1
      allocate(TAB_compacity(i)%tab_compacity_z(nb_ligne_tab+100,3))
      TAB_compacity(i)%tab_compacity_z = 0
    end do
    the_first_time_profils = .false.
  end if
!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  ep = 0
  nb_ligne_tab = 0

  do
!! POUR LE CALCUL DU PROFIL DE VITESSE !!
    Vmoyen   = 0
    Vymoyen  = 0
    Vzmoyen  = 0
    Vxmoyen  = 0
    Vrmoyen   = 0
    Vrymoyen  = 0
    Vrzmoyen  = 0
    Vrxmoyen  = 0
!! POUR LE CALCUL DU PROFIL DE CONTRAINTE !!
!! non, la il n y a rien a initialiser

!! POUR LE CALCUL DU PROFIL DE I !!
!! non, la il n y a rien a initialiser

!! POUR LE CALCUL DU PROFIL DE COMPACITY !!
    volume_sphere_in_band  = 0
    nb_sphere_in_band      = 0
    nb_point_boite         = 0.
    nb_point_sphere        = 0.
    
    cpt_v    = 0
    cpt_c    = 0
    ep = ep+dmoyen/10
    nb_ligne_tab = nb_ligne_tab + 1

!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!! Boucle sur les particules
    do i=1,n_particles
      if (TAB_SPHERE(i)%behav/='PLEXx') cycle
      
      if ((TAB_SPHERE(i)%center(3)<ep+dmoyen/10).and.(TAB_SPHERE(i)%center(3)>ep)) then
!! POUR LE CALCUL DU PROFIL DE VITESSE !!
        Vmoyen = Vmoyen + sqrt(TAB_SPHERE(i)%Vx**2 + TAB_SPHERE(i)%Vy**2 + TAB_SPHERE(i)%Vz**2)
        Vxmoyen = Vxmoyen + TAB_SPHERE(i)%Vx 
        Vymoyen = Vymoyen + TAB_SPHERE(i)%Vy  
        Vzmoyen = Vzmoyen + TAB_SPHERE(i)%Vz  

        Vrmoyen = Vrmoyen + sqrt(TAB_SPHERE(i)%Wx**2 + TAB_SPHERE(i)%Wy**2 + TAB_SPHERE(i)%Wz**2)
        Vrxmoyen = Vrxmoyen + (TAB_SPHERE(i)%Wx)
        Vrymoyen = Vrymoyen + (TAB_SPHERE(i)%Wy) 
        Vrzmoyen = Vrzmoyen + (TAB_SPHERE(i)%Wz) 
    
        cpt_v = cpt_v + 1
      end if
!! POUR LE CALCUL DU PROFIL DE COMPACITY !! (3 cas à distinguer : 
      sphere_in_band = .false.
      if ( abs(TAB_SPHERE(i)%center(3) - ep)            <= TAB_SPHERE(i)%Rmax .and. &
           abs(TAB_SPHERE(i)%center(3) - (ep+dmoyen/10))<= TAB_SPHERE(i)%Rmax ) then
           sphere_in_band = .true. 
           zmax_box = ep+dmoyen/10
           zmin_box = ep
      end if
      if ( abs(TAB_SPHERE(i)%center(3) - ep)            >= TAB_SPHERE(i)%Rmax .and. &
           abs(TAB_SPHERE(i)%center(3) - (ep+dmoyen/10))<= TAB_SPHERE(i)%Rmax ) then
           sphere_in_band = .true. 
           zmax_box = ep+dmoyen/10
           zmin_box = TAB_SPHERE(i)%center(3) - TAB_SPHERE(i)%Rmax
      end if
           
      if ( abs(TAB_SPHERE(i)%center(3) - ep)            <= TAB_SPHERE(i)%Rmax .and. &
           abs(TAB_SPHERE(i)%center(3) - (ep+dmoyen/10))>= TAB_SPHERE(i)%Rmax ) then
           sphere_in_band = .true. 
           zmax_box = TAB_SPHERE(i)%center(3) + TAB_SPHERE(i)%Rmax
           zmin_box = ep
      end if
      xmax_box = TAB_SPHERE(i)%center(1)+TAB_SPHERE(i)%Rmax
      xmin_box = TAB_SPHERE(i)%center(1)-TAB_SPHERE(i)%Rmax
      ymax_box = TAB_SPHERE(i)%center(2)+TAB_SPHERE(i)%Rmax
      ymin_box = TAB_SPHERE(i)%center(2)-TAB_SPHERE(i)%Rmax
      
      if (sphere_in_band .eqv..false.) cycle !La sphere n'est pas dans la bande considérée.
      
      nb_sphere_in_band = nb_sphere_in_band + 1.
      
      !Maintenant on discretise la petite box
      Lx = xmax_box-xmin_box
      Ly = ymax_box-ymin_box
      Lz = zmax_box-zmin_box
      Epsilonx = Lx/50
      Epsilony = Ly/50
      Epsilonz = Lz/50
      test_x = xmin_box
      test_y = ymin_box
      test_z = zmin_box
      ! On remet a zero les compteurs
      nb_point_boite         = 0.
      nb_point_sphere        = 0.
      ! On cherche combien de point il y a dans cette box
      do ix=1,int(Lx/Epsilonx)
        test_x = test_x + Epsilonx
        test_y = ymin_box
        test_z = zmin_box
        do jy=1,int(Ly/Epsilony)
          test_y = test_y + Epsilony
          test_z = zmin_box
          do kz=1,int(Lz/Epsilonz)
            test_z = test_z + Epsilonz
            nb_point_boite = nb_point_boite + 1.D0
          end do   ! sur z
        end do     ! sur y
      end do       ! sur x
      ! On cherche combien de point il y a dans la sphere
      test_x = xmin_box
      test_y = ymin_box
      test_z = zmin_box
      do ix=1,int(Lx/Epsilonx)
        test_x = test_x + Epsilonx
        test_y = ymin_box
        test_z = zmin_box
        do jy=1,int(Ly/Epsilony)
          test_y = test_y + Epsilony
          test_z = zmin_box
          do kz=1,int(Lz/Epsilonz)
            test_z = test_z + Epsilonz
            if ( (test_x-TAB_SPHERE(i)%center(1))**2+&
                 (test_y-TAB_SPHERE(i)%center(2))**2+&
                 (test_z-TAB_SPHERE(i)%center(3))**2   <= TAB_SPHERE(i)%Rmax**2 )  nb_point_sphere = nb_point_sphere + 1.D0
          end do   ! sur z
        end do     ! sur y
      end do       ! sur x
      Vol_sphere = (nb_point_sphere/nb_point_boite) * Lx*Ly*Lz
      volume_sphere_in_band = volume_sphere_in_band + Vol_sphere
    end do
    
!! Boucle sur les contacts
    do i=1,nb_ligneCONTACT
      if ( (TAB_CONTACT(i)%coor_ctc(3) < ep+ dmoyen/10).and.(TAB_CONTACT(i)%coor_ctc(3) > ep) ) then
        if  (TAB_CONTACT(i)%nature == 'SPPLx') cycle
        cd   = TAB_CONTACT(i)%icdent
        an   = TAB_CONTACT(i)%ianent
        nik  = TAB_CONTACT(i)%n
        tik  = TAB_CONTACT(i)%t
        sik  = TAB_CONTACT(i)%s
        Rtik = TAB_CONTACT(i)%rt
        Rnik = TAB_CONTACT(i)%rn
        Rsik = TAB_CONTACT(i)%rs
        if (Rnik==0) cycle
        Lik(1:3) = TAB_SPHERE(cd)%center(1:3)-TAB_SPHERE(an)%center(1:3)
        Fik(1:3) = (Rnik*nik(1:3)+Rtik*tik(1:3)+Rsik*sik(1:3))
    
        if ( sqrt( Lik(1)**2 + Lik(2)**2 + Lik(3)**2) > 2*( TAB_SPHERE(cd)%Rmax + TAB_SPHERE(an)%Rmax ) ) then
          Rayon_max = sqrt( ( Lik(1)**2 + Lik(2)**2 + Lik(3)**2) ) 
          lik(:) = Lik(:) / Rayon_max
          Lik(:) = abs(TAB_SPHERE(cd)%Rmax+TAB_SPHERE(an)%Rmax)*Lik(:)
        end if

!! POUR LE CALCUL DU PROFIL DE CONTRAINTE !!
        TAB_contrainte(compteur_clout-debut_profils+1)%tab_sigma(nb_ligne_tab)%sigma(1,1:3) =&
        Fik(1)*Lik(1:3) + TAB_contrainte(compteur_clout-debut_profils+1)%tab_sigma(nb_ligne_tab)%sigma(1,1:3) 

        TAB_contrainte(compteur_clout-debut_profils+1)%tab_sigma(nb_ligne_tab)%sigma(2,1:3) =&
        Fik(2)*Lik(1:3) + TAB_contrainte(compteur_clout-debut_profils+1)%tab_sigma(nb_ligne_tab)%sigma(2,1:3)

        TAB_contrainte(compteur_clout-debut_profils+1)%tab_sigma(nb_ligne_tab)%sigma(3,1:3) =&
        Fik(3)*Lik(1:3) + TAB_contrainte(compteur_clout-debut_profils+1)%tab_sigma(nb_ligne_tab)%sigma(3,1:3)
 
        cpt_c = cpt_c + 1
      end if
    end do
    
    !vboite = Long*Larg*dmoyen/10

!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    if (cpt_v==0) then
!! POUR LE CALCUL DU PROFIL DE VITESSE !!
      Vmoyen  = 0.D0
      Vxmoyen = 0.D0
      Vymoyen = 0.D0
      Vzmoyen = 0.D0

      Vrmoyen  = 0.D0
      Vrxmoyen = 0.D0
      Vrymoyen = 0.D0
      Vrzmoyen = 0.D0
    else if (cpt_v>0)  then
!! POUR LE CALCUL DU PROFIL DE VITESSE !!
      Vmoyen  = Vmoyen / real(cpt_v,8)
      Vxmoyen = Vxmoyen / real(cpt_v,8)
      Vymoyen = Vymoyen / real(cpt_v,8)
      Vzmoyen = Vzmoyen / real(cpt_v,8)

      Vrmoyen  = Vrmoyen / real(cpt_v,8)
      Vrxmoyen = Vrxmoyen / real(cpt_v,8)
      Vrymoyen = Vrymoyen / real(cpt_v,8)
      Vrzmoyen = Vrzmoyen / real(cpt_v,8)
    end if

!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ fin-debut_profils+1
    
!! POUR LE CALCUL DU PROFIL DE VITESSE !!

    TAB_vitesse(compteur_clout-debut_profils+1)%tab_vmoyen(nb_ligne_tab,1) = ep / dmoyen
    TAB_vitesse(compteur_clout-debut_profils+1)%tab_vmoyen(nb_ligne_tab,2) = Vxmoyen
    TAB_vitesse(compteur_clout-debut_profils+1)%tab_vmoyen(nb_ligne_tab,3) = Vymoyen
    TAB_vitesse(compteur_clout-debut_profils+1)%tab_vmoyen(nb_ligne_tab,4) = Vzmoyen

    TAB_vitesse(compteur_clout-debut_profils+1)%tab_vrmoyen(nb_ligne_tab,1) = ep
    TAB_vitesse(compteur_clout-debut_profils+1)%tab_vrmoyen(nb_ligne_tab,2) = Vrxmoyen
    TAB_vitesse(compteur_clout-debut_profils+1)%tab_vrmoyen(nb_ligne_tab,3) = Vrymoyen
    TAB_vitesse(compteur_clout-debut_profils+1)%tab_vrmoyen(nb_ligne_tab,4) = Vrzmoyen

!! POUR LE CALCUL DU PROFIL DE CONTRAINTE !!
    TAB_contrainte(compteur_clout-debut_profils+1)%tab_sigma(nb_ligne_tab)%sigma = &  
    TAB_contrainte(compteur_clout-debut_profils+1)%tab_sigma(nb_ligne_tab)%sigma / vboite

!! POUR LE CALCUL DU PROFIL DE COMPACITY !!

    TAB_compacity(compteur_clout-debut_profils+1)%tab_compacity_z(nb_ligne_tab,1) = ep / dmoyen
    TAB_compacity(compteur_clout-debut_profils+1)%tab_compacity_z(nb_ligne_tab,2) = volume_sphere_in_band/vboite
                           
    !if (ep>Haut+dmoyen) exit
  end do

!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  if (compteur_clout == fin) then
!! POUR LE CALCUL DU PROFIL DE VITESSE !!
    allocate (Moyenne_vitesse_moyenne(nb_ligne_tab+100,4))
    allocate (Moyenne_vitesse_moyenne_r(nb_ligne_tab+100,5))
    Moyenne_vitesse_moyenne(:,:) = 0.D0
    Moyenne_vitesse_moyenne_r(:,:) = 0.D0
!! POUR LE CALCUL DU PROFIL DE CONTRAINTE !!
    allocate (Moyenne_contrainte_moyenne(nb_ligne_tab+100))
    do i=1,nb_ligne_tab+100
      Moyenne_contrainte_moyenne(i)%sigma = 0.D0
    end do
!! POUR LE CALCUL DU PROFIL DE I !!
    allocate (Moyenne_I_moyen(nb_ligne_tab+100,4))
!! POUR LE CALCUL DU PROFIL DE COMPACITY !!
    allocate (Moyenne_compacity_moyen(nb_ligne_tab+200,3))

    do i=1,nb_ligne_tab+100
      do k= 1,fin-debut_profils+1
!! POUR LE CALCUL DU PROFIL DE VITESSE !!
        Moyenne_vitesse_moyenne(i,1) = Moyenne_vitesse_moyenne(i,1) + TAB_vitesse(k)%tab_vmoyen(i,1)
        Moyenne_vitesse_moyenne(i,2) = Moyenne_vitesse_moyenne(i,2) + TAB_vitesse(k)%tab_vmoyen(i,2)
        Moyenne_vitesse_moyenne(i,3) = Moyenne_vitesse_moyenne(i,3) + TAB_vitesse(k)%tab_vmoyen(i,3)
        Moyenne_vitesse_moyenne(i,4) = Moyenne_vitesse_moyenne(i,4) + TAB_vitesse(k)%tab_vmoyen(i,4)

        Moyenne_vitesse_moyenne_r(i,1) = Moyenne_vitesse_moyenne_r(i,1) + TAB_vitesse(k)%tab_vrmoyen(i,1)
        Moyenne_vitesse_moyenne_r(i,2) = Moyenne_vitesse_moyenne_r(i,2) + TAB_vitesse(k)%tab_vrmoyen(i,2)
        Moyenne_vitesse_moyenne_r(i,3) = Moyenne_vitesse_moyenne_r(i,3) + TAB_vitesse(k)%tab_vrmoyen(i,3)
        Moyenne_vitesse_moyenne_r(i,4) = Moyenne_vitesse_moyenne_r(i,4) + TAB_vitesse(k)%tab_vrmoyen(i,4)

!! POUR LE CALCUL DU PROFIL DE CONTRAINTE !!
        Moyenne_contrainte_moyenne(i)%sigma = &
        Moyenne_contrainte_moyenne(i)%sigma + TAB_contrainte(k)%tab_sigma(i)%sigma

!! POUR LE CALCUL DU PROFIL DE COMPACITY !!
        Moyenne_compacity_moyen(i,1) = Moyenne_compacity_moyen(i,1) + TAB_compacity(k)%tab_compacity_z(i,1)
        Moyenne_compacity_moyen(i,2) = Moyenne_compacity_moyen(i,2) + TAB_compacity(k)%tab_compacity_z(i,2)

      end do
    end do
!! POUR LE CALCUL DU PROFIL DE VITESSE !!
    Moyenne_vitesse_moyenne(:,:)   = Moyenne_vitesse_moyenne(:,:) / real(fin-debut_profils+1,8)
    Moyenne_vitesse_moyenne_r(:,:) = Moyenne_vitesse_moyenne_r(:,:) / real(fin-debut_profils+1,8)

!! POUR LE CALCUL DU PROFIL DE CONTRAINTE !!
    do i=1,nb_ligne_tab+100
      Moyenne_contrainte_moyenne(i)%sigma(:,:)= Moyenne_contrainte_moyenne(i)%sigma(:,:) / real(fin-debut_profils+1,8)
    end do

!! POUR LE CALCUL DU PROFIL DE I !!
    do i=1,nb_ligne_tab+100
      Mean_pressure = (Moyenne_contrainte_moyenne(i)%sigma(1,1)+&
                       Moyenne_contrainte_moyenne(i)%sigma(2,2)+&
                       Moyenne_contrainte_moyenne(i)%sigma(3,3))/3
                       
      Moyenne_I_moyen(i,1) =  Moyenne_vitesse_moyenne(i,1)
      Moyenne_I_moyen(i,2) = ( Moyenne_vitesse_moyenne(i+1,2)-Moyenne_vitesse_moyenne(i,2) ) / &      ! ATTENTION IL FAUT SE RAPPELER QU ON
                             ( dmoyen*(Moyenne_vitesse_moyenne(i+1,1)-Moyenne_vitesse_moyenne(i,1)) )  
      Moyenne_I_moyen(i,3) =  Moyenne_I_moyen(i,2) * dmoyen * sqrt(1000/Mean_pressure)
      Moyenne_I_moyen(i,4) =  abs(2*Moyenne_vitesse_moyenne_r(i,3)) * dmoyen * sqrt(1000/Mean_pressure)
    end do    

!! POUR LE CALCUL DU PROFIL DE COMPACITY !!
    Moyenne_compacity_moyen(:,:)   = Moyenne_compacity_moyen(:,:) / real(fin-debut_profils+1,8)

!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!! POUR LE CALCUL DU PROFIL DE VITESSE !!
    do i=1,nb_ligne_tab+20
      write(131,'(11(1X,D14.7))') Moyenne_vitesse_moyenne(i,1),Moyenne_vitesse_moyenne(i,2),&
                                  Moyenne_vitesse_moyenne(i,3),Moyenne_vitesse_moyenne(i,4),&
                                  Moyenne_vitesse_moyenne_r(i,2),Moyenne_vitesse_moyenne_r(i,3),&
                                  Moyenne_vitesse_moyenne_r(i,4)
!! POUR LE CALCUL DU PROFIL DE COMPACITY !!
      write(134,'(11(1X,D14.7))') Moyenne_compacity_moyen(i,1),Moyenne_compacity_moyen(i,2)
    end do
    
    
!! POUR LE CALCUL DU PROFIL DE CONTRAINTE !!
    do i=1,nb_ligne_tab+20
      lda_sigma        = 3
      matz_sigma       = 1
      M_sigma          = Moyenne_contrainte_moyenne(i)%sigma
      wr_sigma         = 0.D0
      wi_sigma         = 0.D0
      localframe_sigma = 0.D0
      call rg ( lda_sigma, 3, M_sigma, wr_sigma, wi_sigma, matz_sigma, localframe_sigma, ierror_sigma )
      if ( wr_sigma(1)==wr_sigma(2) .and. (wr_sigma(2)==wr_sigma(3)) ) then
        S1_sigma = wr_sigma(1)
        S2_sigma = S1_sigma
        S3_sigma = S1_sigma
      else
        S3_sigma = max( wr_sigma(1),max(wr_sigma(2),wr_sigma(3)))
        S1_sigma = min( wr_sigma(1),min(wr_sigma(2),wr_sigma(3)))
        if (wr_sigma(1)==wr_sigma(2)) S2_sigma=S1_sigma
        if (wr_sigma(3)==wr_sigma(2)) S2_sigma=S3_sigma
        do j=1,3
          if ((wr_sigma(j)<S3_sigma) .and. (wr_sigma(j)>S1_sigma)) S2_sigma=wr_sigma(j)
        end do
      end if

      if (Moyenne_contrainte_moyenne(i)%sigma(3,3) == 0 ) then
        S_P = 0.D0
      else 
        S_P = abs( Moyenne_contrainte_moyenne(i)%sigma(2,3) / Moyenne_contrainte_moyenne(i)%sigma(3,3) )
      end if

      write(132,'(11(1X,D14.7))') Moyenne_vitesse_moyenne(i,1), &
                                  (S3_sigma-S1_sigma)/(S3_sigma+S2_sigma+S1_sigma),&
                                  S_P
      
    end do

!! POUR LE CALCUL DU PROFIL DE I !!
    do i=1,nb_ligne_tab+20
      write(133,'(11(1X,D14.7))') Moyenne_vitesse_moyenne(i,1),Moyenne_I_moyen(i,2),Moyenne_I_moyen(i,3),Moyenne_I_moyen(i,4)
      !                           z, \dot \gamma, I, I avec les rotations
    end do

!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  end if
end if  
end subroutine profils


!==============================================================================
! Calcul du PDF
!==============================================================================
subroutine pdf
implicit none
real(kind=8)                             :: increment
real(kind=8),dimension(3,3)              :: Fabric,Fabric_LN,Fabric_LT,Fabric_L,localframe
real(kind=8),dimension(3,3)              :: Fabric_FN,Fabric_FT,Fabric_F
real(kind=8),dimension(3)                :: nik,tik,sik,Lik,Fik
real(kind=8)                             :: Rtik,Rnik,Rsik,norm,FN,FT,LN,LT,Norm_F,Norm_L,cpt,&
                                            ac1,ac2,acS,acS2,Sa,Sa_NS,Theta_c,Theta_fabric,&
                                            Theta_c_s,Theta_fabric_s,S1,S2,S3,theta2
integer                                  :: i,cd,an,j
real(kind=8),dimension(3,3)              :: a,aln,alt,afn,aft
real(kind=8),dimension(3)                :: wr,wi
integer                                  :: ierror,matz,lda

print*, compteur_clout,step_pdf_start,step_pdf_stop

do i=1,n_particles
  TAB_SPHERE(i)%Fabric = 0.D0
  TAB_SPHERE(i)%Fabric_N = 0.D0
  TAB_SPHERE(i)%Fabric_T = 0.D0
  TAB_SPHERE(i)%nb_ctc = 0
  TAB_SPHERE(i)%Norm_F = 0.D0
end do

if ( (compteur_clout >= step_pdf_start).and.(compteur_clout <= step_pdf_stop) ) then

  print*, '---> PDF'
  do i = 1, nb_ligneCONTACT
    cd   = TAB_CONTACT(i)%icdent
    an   = TAB_CONTACT(i)%ianent
    if (TAB_CONTACT(i)%rn<0.000000000001*Mean_total_Normal_force) cycle
    if ((TAB_SPHERE(cd)%behav/='PLEXx').or.(TAB_SPHERE(an)%behav/='PLEXx')) cycle
    write(106,'(4(1X,D12.5))') TAB_CONTACT(i)%rn,&
                                sqrt(TAB_CONTACT(i)%rs**2+TAB_CONTACT(i)%rt**2)
  end do
  
  do i=1,nb_ligneCONTACT
    if  (TAB_CONTACT(i)%nature == 'SPPLx') cycle 
    if  (TAB_CONTACT(i)%rn==0.D0) cycle
    cd   = TAB_CONTACT(i)%icdent
    an   = TAB_CONTACT(i)%ianent
    nik  = TAB_CONTACT(i)%n
    tik  = TAB_CONTACT(i)%t
    sik  = TAB_CONTACT(i)%s
    Rtik = TAB_CONTACT(i)%rt
    Rnik = TAB_CONTACT(i)%rn
    Rsik = TAB_CONTACT(i)%rs

    TAB_SPHERE(cd)%nb_ctc = TAB_SPHERE(cd)%nb_ctc + 1
    TAB_SPHERE(an)%nb_ctc = TAB_SPHERE(an)%nb_ctc + 1

    Fik(1:3) = (Rnik*nik(1:3)+Rtik*tik(1:3)+Rsik*sik(1:3))
    Norm_F   = Norm_F + sqrt( Fik(1)**2 + Fik(2)**2 + Fik(3)**2 )
  
    Lik(1:3) = TAB_SPHERE(cd)%center(1:3)-TAB_SPHERE(an)%center(1:3)
    Norm_L   = Norm_L + sqrt( Lik(1)**2 + Lik(2)**2 + Lik(3)**2 )

    LN       =  sqrt( Lik(1)**2 + Lik(2)**2 + Lik(3)**2 ) 
    FN       =  Rnik
   
    tik(1:3) =  Fik(1:3) - FN*nik(1:3)
    FT       =  sqrt( tik(1)**2 + tik(2)**2 + tik(3)**2 )

    if ( FT  > 0.00000000001) tik(1:3) = tik(1:3) / FT
    if ( FT <= 0.00000000001) tik(1:3) = 0.D0

    FT       =  Fik(1)*tik(1) + Fik(2)*tik(2) + Fik(3)*tik(3)  

    TAB_SPHERE(cd)%Norm_F = TAB_SPHERE(cd)%Norm_F + sqrt( FN**2 + FT**2 ) 
    TAB_SPHERE(an)%Norm_F = TAB_SPHERE(an)%Norm_F + sqrt( FN**2 + FT**2 ) 

    Fabric(1,1:3) = Fabric(1,1:3) + nik(1)*nik(1:3)
    Fabric(2,1:3) = Fabric(2,1:3) + nik(2)*nik(1:3)
    Fabric(3,1:3) = Fabric(3,1:3) + nik(3)*nik(1:3)

    TAB_SPHERE(cd)%Fabric(1,1:3) = TAB_SPHERE(cd)%Fabric(1,1:3) + nik(1)*nik(1:3)
    TAB_SPHERE(cd)%Fabric(2,1:3) = TAB_SPHERE(cd)%Fabric(2,1:3) + nik(2)*nik(1:3)
    TAB_SPHERE(cd)%Fabric(3,1:3) = TAB_SPHERE(cd)%Fabric(3,1:3) + nik(3)*nik(1:3) 
    TAB_SPHERE(an)%Fabric(1,1:3) = TAB_SPHERE(an)%Fabric(1,1:3) + nik(1)*nik(1:3)
    TAB_SPHERE(an)%Fabric(2,1:3) = TAB_SPHERE(an)%Fabric(2,1:3) + nik(2)*nik(1:3)
    TAB_SPHERE(an)%Fabric(3,1:3) = TAB_SPHERE(an)%Fabric(3,1:3) + nik(3)*nik(1:3)

    cpt = cpt + 1
  end do

  do i=1,n_particles
    if (TAB_SPHERE(i)%nb_ctc<2) cycle

   TAB_SPHERE(i)%Fabric(:,:) = TAB_SPHERE(i)%Fabric(:,:) / TAB_SPHERE(i)%nb_ctc

    acS = ( 5*(TAB_SPHERE(i)%Fabric(3,3)-TAB_SPHERE(i)%Fabric(1,1))/2 ) 
    acS2= ( 5*(TAB_SPHERE(i)%Fabric(3,3)-TAB_SPHERE(i)%Fabric(2,2))/2 )
     
    lda  = 3
    matz = 1
    call rg ( lda, 3, TAB_SPHERE(i)%Fabric, wr, wi, matz, localframe, ierror )
    if ( wr(1)==wr(2) .and. (wr(2)==wr(3)) ) then
      S1=wr(1)
      S2=S1
      S3=S1
    else
      S3 = max( wr(1),max(wr(2),wr(3)))
      S1 = min( wr(1),min(wr(2),wr(3)))
      if (wr(1)==wr(2)) S2=S1
      if (wr(3)==wr(2)) S2=S3
      do j=1,3
        if ((wr(j)<S3) .and. (wr(j)>S1)) S2=wr(j)
      end do
    end if



    ac1 = 5*(S3-S1)/2
    ac2 = 5*(S3-S2)/2

    theta2 = acS2/ac2

    if (abs(theta2) >1) theta2=0 

    write(108,'(10(1X,D12.5))') acS,acS2,ac1,ac2,acS/ac1,theta2
    
    if (TAB_SPHERE(i)%nb_ctc==2) write(109,'(10(1X,D12.5))') acS,acS2,ac1,ac2,acS/ac1,theta2
    if (TAB_SPHERE(i)%nb_ctc==3) write(110,'(10(1X,D12.5))') acS,acS2,ac1,ac2,acS/ac1,theta2
    if (TAB_SPHERE(i)%nb_ctc==4) write(111,'(10(1X,D12.5))') acS,acS2,ac1,ac2,acS/ac1,theta2
    if (TAB_SPHERE(i)%nb_ctc==5) write(112,'(10(1X,D12.5))') acS,acS2,ac1,ac2,acS/ac1,theta2
    if (TAB_SPHERE(i)%nb_ctc==6) write(113,'(10(1X,D12.5))') acS,acS2,ac1,ac2,acS/ac1,theta2
    if (TAB_SPHERE(i)%nb_ctc==7) write(114,'(10(1X,D12.5))') acS,acS2,ac1,ac2,acS/ac1,theta2
    if (TAB_SPHERE(i)%nb_ctc==8) write(115,'(10(1X,D12.5))') acS,acS2,ac1,ac2,acS/ac1,theta2
    if (TAB_SPHERE(i)%nb_ctc==9) write(116,'(10(1X,D12.5))') acS,acS2,ac1,ac2,acS/ac1,theta2
    if (TAB_SPHERE(i)%nb_ctc==10) write(117,'(10(1X,D12.5))') acS,acS2,ac1,ac2,acS/ac1,theta2
    if (TAB_SPHERE(i)%nb_ctc==11) write(118,'(10(1X,D12.5))') acS,acS2,ac1,ac2,acS/ac1,theta2
    if (TAB_SPHERE(i)%nb_ctc==12) write(119,'(10(1X,D12.5))') acS,acS2,ac1,ac2,acS/ac1,theta2
    
  end do

  
else if ( (compteur_clout >= step_pdf_stop ) ) then
  stop
end if 

end subroutine pdf


!==============================================================================
! Vitesse moyenne
!==============================================================================
subroutine vitesse_moyenne
implicit none
real(kind=8)                             :: V,Vx,Vy,Vz,Vrx,cpt,el_diameter,zc,viscosity,&
                                            I_Calculated,Ec,V_center_mass,&
                                            Nombre_Inertie,mean_mass,applied_pressure
integer                                  :: i

print*,'.... Vitesse moyenne ....'

V   = 0.D0
Vx  = 0.D0
Vy  = 0.D0
Vz  = 0.D0
cpt = 0.D0
Vrx = 0.D0
zc  = 0.D0
Ec  = 0
V_center_mass = 0
Nombre_Inertie = 0
mean_mass = 0

cpt = 0
do i=1,n_particles
  if (TAB_SPHERE(i)%behav/='PLEXx') cycle
  V  = V  + TAB_SPHERE(i)%Norme_V
  Vx = Vx + TAB_SPHERE(i)%Vx
  Vy = Vy + TAB_SPHERE(i)%Vy
  Vz = Vz + TAB_SPHERE(i)%Vz
  Vrx = Vrx + TAB_SPHERE(i)%Wx
  cpt = cpt + 1
  el_diameter = el_diameter + 2*TAB_SPHERE(i)%Rmax
  mean_mass   = mean_mass + 1000 * TAB_SPHERE(i)%volume
end do

V  = V/cpt
Vx = Vx/cpt
Vy = Vy/cpt
Vz = Vz/cpt
Vrx = Vrx/cpt
el_diameter = el_diameter / cpt
applied_pressure = 6.5790026 !/ (Long*Larg)

!I_Calculated   = (0.047472569/Haut) * el_diameter * sqrt( 1000./ pressure_in_sample )
!Nombre_Inertie = (0.047472569/Haut) * el_diameter * sqrt( 1000./ applied_pressure   )

!write(101,'(25(1X,D12.5))') time,Def,Haut,&
                            !V,Vx,Vy,Vz,Vrx,&
                            !Nombre_Inertie,I_Calculated

 end subroutine vitesse_moyenne


!==================================================================================================
!POUR UN SEUL PAS DE TEMPS 
!=================================================================================================




!==============================================================================
! Orientation des contacts
!==============================================================================
subroutine contact_orientation
implicit none
integer                                             ::  i,cd,an,Nsect,j,cpt_theta,cpt_phi,cpt_zero,k,cpt
real(kind=8)                                        ::  sect,sect2,theta,phi,fmoyen,norm,SinTheta,ThetaMax,ThetaMin,ThetaMean,FT,&
                                                        Rnik,Rtik,Rsik,fnorm
real(kind=8),dimension(3)                           ::  Lik,tik,nik,sik,Fik
real(kind=8),dimension(:),allocatable               ::  tab_alpha,tab_alpha2
real(kind=8),dimension(:),allocatable               ::  tab_number_n_theta,tab_number_n_xy,tab_number_v2_moyen
real(kind=8),dimension(3)                           ::  V
real(kind=8),dimension(3,3)                         ::  Fabric
print*,'.... Orientation des contacts ...'

Nsect = 40

allocate(tab_alpha2(Nsect))
tab_alpha2 = 0
sect2 = PI / real(Nsect,8)
tab_alpha2(1)=sect2
do i=2,Nsect
  tab_alpha2(i) = tab_alpha2(i-1) + sect2
end do
allocate(tab_number_n_xy(Nsect))
tab_number_n_xy = 0


cpt=0
Fabric(:,:) = 0.D0
do i=1,nb_ligneCONTACT
  if  (TAB_CONTACT(i)%nature == 'SPPLx') cycle
  if (TAB_CONTACT(i)%rn<0.000000000001*Mean_total_Normal_force) cycle

  nik  = TAB_CONTACT(i)%n

  Fabric(1,1:3) = nik(1)*nik(1:3) + Fabric(1,1:3)
  Fabric(2,1:3) = nik(2)*nik(1:3) + Fabric(2,1:3)
  Fabric(3,1:3) = nik(3)*nik(1:3) + Fabric(3,1:3)
  cpt=cpt+1
end do
Fabric(:,:) = Fabric(:,:) / real(cpt,8)
print*,Fabric(1,1),Fabric(2,2),Fabric(3,3)
cpt = 0

!----------------ORIENTATION DES NORMALES DANS LE PLAN XY ------------------
cpt_phi   = 0
do i=1,nb_ligneCONTACT
  if (TAB_CONTACT(i)%nature == 'SPPLx') cycle
  if (TAB_CONTACT(i)%rn<0.000000000001*Mean_total_Normal_force) cycle
  cd  = TAB_CONTACT(i)%icdent
  an  = TAB_CONTACT(i)%ianent
  Lik(1:3) = TAB_SPHERE(cd)%center(1:3)-TAB_SPHERE(an)%center(1:3)
  norm = sqrt( Lik(1)**2 + Lik(2)**2 + Lik(3)**2 )
  Lik(1:3) = Lik(1:3) / norm
  if (Lik(1)**2+Lik(2)**2 ==0.D0 ) cycle
  cpt_phi = cpt_phi + 1
end do

do i=1,Nsect
  do j=1,nb_ligneCONTACT
    cd  = TAB_CONTACT(j)%icdent
    an  = TAB_CONTACT(j)%ianent
    Lik(1:3) = TAB_SPHERE(cd)%center(1:3)-TAB_SPHERE(an)%center(1:3)
    norm = sqrt( Lik(1)**2 + Lik(2)**2 + Lik(3)**2 )
    Lik(1:3) = Lik(1:3) / norm

    if (Lik(1)**2+Lik(2)**2 ==0.D0 ) cycle

    phi = acos( Lik(1) / sqrt(Lik(1)**2+Lik(2)**2) )
        
    if (i==1) then
      if (phi<=tab_alpha2(i)) then
        tab_number_n_xy(i) = tab_number_n_xy(i) + 1
      end if
    end if

    if ( (i>1).and.(i<Nsect)) then
      if ((phi>tab_alpha2(i-1)).and.(phi<=tab_alpha2(i))) then
        tab_number_n_xy(i) = tab_number_n_xy(i) + 1
      end if
    end if

    if (i==Nsect) then
      if ((phi<=tab_alpha2(i)).and.(phi>tab_alpha2(i-1))) then
        tab_number_n_xy(i) = tab_number_n_xy(i) + 1
      end if
    end if
  end do
end do
tab_number_n_xy(:) = tab_number_n_xy(:) / real(cpt_phi,8)


do i=1,Nsect
  write(126,'(8(1X,D14.7))') sect2*(i-1) , tab_number_n_xy(i)
  write(126,'(8(1X,D14.7))') sect2*(i)   , tab_number_n_xy(i)
end do
do i=1,Nsect
  write(126,'(8(1X,D14.7))') sect2*(i-1)+PI , tab_number_n_xy(i)
  write(126,'(8(1X,D14.7))') sect2*(i)+PI   , tab_number_n_xy(i)
end do

!----------------ORIENTATION DES NORMALES DANS LE PLAN YZ ------------------

Nsect = 40
allocate(tab_alpha(Nsect))
tab_alpha = 0
sect =   Pi / real(Nsect,8)
tab_alpha(1)=0.D0
do i=2,Nsect
  tab_alpha(i) = tab_alpha(i-1) + sect
end do
allocate(tab_number_n_theta(Nsect))
allocate(tab_number_v2_moyen(Nsect))

tab_number_n_theta  = 0
tab_number_v2_moyen = 0
cpt_theta   = 0
do i=1,nb_ligneCONTACT
  if (TAB_CONTACT(i)%nature == 'SPPLx') cycle
  if (TAB_CONTACT(i)%rn<0.000000000001*Mean_total_Normal_force) cycle
  cd  = TAB_CONTACT(i)%icdent
  an  = TAB_CONTACT(i)%ianent
  
  if ((TAB_SPHERE(cd)%behav/='PLEXx').or.(TAB_SPHERE(an)%behav/='PLEXx')) cycle
  
!  Lik(1:3) = TAB_SPHERE(cd)%center(1:3)-TAB_SPHERE(an)%center(1:3)
  Lik(1:3) = TAB_CONTACT(i)%n
  norm = sqrt( Lik(1)**2 + Lik(2)**2 + Lik(3)**2 )
  Lik(1:3) = Lik(1:3) / norm
  phi = acos( (Lik(1)) / sqrt(Lik(1)**2+Lik(2)**2) )
  if ( (phi> Pi/2 + 2*sect2).and.(phi< 3*Pi/2 - 2*sect2)  ) cycle
  if ( (phi< Pi/2 - 2*sect2).and.(phi> 3*Pi/2 + 2*sect2)  ) cycle

  cpt_theta = cpt_theta + 1
end do

do i=1,Nsect
  do j=1,nb_ligneCONTACT
    cd  = TAB_CONTACT(j)%icdent
    an  = TAB_CONTACT(j)%ianent
    if ((TAB_SPHERE(cd)%behav/='PLEXx').or.(TAB_SPHERE(an)%behav/='PLEXx')) cycle
    if (TAB_CONTACT(j)%rn<0.0000001*Mean_total_Normal_force) cycle

    nik  = TAB_CONTACT(j)%n
    tik  = TAB_CONTACT(j)%t
    sik  = TAB_CONTACT(j)%s
    
    Rtik = TAB_CONTACT(j)%rt
    Rnik = TAB_CONTACT(j)%rn
    Rsik = TAB_CONTACT(j)%rs
   
    Fik(1:3) = (Rnik*nik(1:3)+Rtik*tik(1:3)+Rsik*sik(1:3))
    fnorm    = sqrt( Fik(1)**2 + Fik(2)**2 + Fik(3)**2)
        
!    tik(1:3) =  Fik(1:3) - TAB_CONTACT(j)%rn*nik(1:3)
!    FT =  sqrt( tik(1)**2 + tik(2)**2 + tik(3)**2 )
    
!    tik = tik / FT

    tik(1) = 0.
    tik(2) = TAB_CONTACT(j)%n(3)
    tik(3) =-TAB_CONTACT(j)%n(2)
    
    FT = Fik(1)*tik(1) + Fik(2)*tik(2) + Fik(3)*tik(3) 
    
!    if ( nik(1)*tik(3)-nik(3)*tik(1) < 0 ) FT =  -FT 

    Lik(1:3) = TAB_CONTACT(j)%n
    if ( Lik(3) < 0.D0 )  Lik(:)=-Lik(:)
    
    phi = acos( (Lik(1)) / sqrt(Lik(1)**2+Lik(2)**2) )

    if ( (phi> Pi/2 + 2*sect2).and.(phi< 3*Pi/2 - 2*sect2)  ) cycle
    if ( (phi< Pi/2 - 2*sect2).and.(phi> 3*Pi/2 + 2*sect2)  ) cycle

    theta = acos( Lik(2) / sqrt( Lik(2)**2+Lik(3)**2 ) )

    if (i==1) then
      if ( (theta)<=(tab_alpha(i)) ) then
!        tab_number_n_theta(i) = tab_number_n_theta(i) + fnorm  !/Mean_total_Normal_force
        tab_number_n_theta(i) = tab_number_n_theta(i) + FT  /Mean_total_Normal_force
!        tab_number_n_theta(i) = tab_number_n_theta(i) + Rnik/Mean_total_Normal_force
!       tab_number_n_theta(i) = tab_number_n_theta(i) + 1 !sqrt(1-Lik(3)**2)
      end if
    end if

    if ( (i>1).and.(i<Nsect)) then
      if ( ((theta)>=(tab_alpha(i-1)) ) .and. ( (theta)<=(tab_alpha(i)) ) ) then
!        tab_number_n_theta(i) = tab_number_n_theta(i) + fnorm  !/Mean_total_Normal_force
        tab_number_n_theta(i) = tab_number_n_theta(i) + FT  /Mean_total_Normal_force
!        tab_number_n_theta(i) = tab_number_n_theta(i) + Rnik/Mean_total_Normal_force
!       tab_number_n_theta(i) = tab_number_n_theta(i) + 1 !sqrt(1-Lik(3)**2)
      end if
    end if

    if (i==Nsect) then
      if ( ((theta)<=(tab_alpha(i))) .and. ((theta)>=(tab_alpha(i-1))) ) then
!        tab_number_n_theta(i) = tab_number_n_theta(i) + fnorm  !/Mean_total_Normal_force
        tab_number_n_theta(i) = tab_number_n_theta(i) + FT  /Mean_total_Normal_force
!        tab_number_n_theta(i) = tab_number_n_theta(i) + Rnik/Mean_total_Normal_force
!       tab_number_n_theta(i) = tab_number_n_theta(i) + 1 !sqrt(1-Lik(3)**2)
      end if
    end if
  end do
end do
tab_number_n_theta(:)  = tab_number_n_theta(:) / real(cpt_theta,8)

do i=1,Nsect
  write(124,'(8(1X,D14.7))') sect*(i-1) , tab_number_n_theta(i)
  write(124,'(8(1X,D14.7))') sect*(i)   , tab_number_n_theta(i)
end do
do i=1,Nsect
  write(124,'(8(1X,D14.7))') sect*(i-1)+PI , tab_number_n_theta(i)
  write(124,'(8(1X,D14.7))') sect*(i)+PI   , tab_number_n_theta(i)
end do
!do i=1,Nsect
!  write(124,'(8(1X,D14.7))') -1.D0 + sect*(i-1) , tab_number_n_theta(i)
!  write(124,'(8(1X,D14.7))') -1.D0 + sect*(i)   , tab_number_n_theta(i)
!end do

stop
end subroutine contact_orientation


!==============================================================================
! Orientation des contacts EN 3D TENTATIVE !!!!!
!==============================================================================
subroutine contact_orientation_3D
implicit none
integer                                             ::  i,cd,an,Nsect,j,cpt_theta,cpt_phi,k
real(kind=8)                                        ::  sect,theta,phi,norm
real(kind=8),dimension(3)                           ::  nik,Lik
real(kind=8),dimension(:),allocatable               ::  tab_alpha_theta,tab_alpha_phi
real(kind=8),dimension(:,:),allocatable             ::  tab_number_n

print*,'.... Orientation des contacts ...'

Nsect = 40
allocate(tab_alpha_theta(Nsect))
allocate(tab_alpha_phi(Nsect))
tab_alpha_theta = 0
tab_alpha_phi = 0
sect = PI / real(Nsect,8)
tab_alpha_theta(1)  =  0.D0
do i=2,Nsect
  tab_alpha_theta(i) = tab_alpha_theta(i-1) + sect
end do
tab_alpha_phi(1)    =  0.D0
do i=2,Nsect
  tab_alpha_phi(i) = tab_alpha_phi(i-1) + 2*sect
end do
allocate(tab_number_n(Nsect,Nsect))

tab_number_n = 0

cpt_phi   = 0
cpt_theta = 0
do i=1,nb_ligneCONTACT
  if (TAB_CONTACT(i)%nature == 'SPPLx') cycle
  if (TAB_CONTACT(i)%rn==0.D0)          cycle
  cd  = TAB_CONTACT(i)%icdent
  an  = TAB_CONTACT(i)%ianent
  Lik(1:3) = TAB_SPHERE(cd)%center(1:3)-TAB_SPHERE(an)%center(1:3)
  norm = sqrt( Lik(1)**2 + Lik(2)**2 + Lik(3)**2 )
  Lik(1:3) = Lik(1:3) / norm
  if (Lik(1)**2+Lik(2)**2 == 0.D0 ) cycle
  cpt_phi = cpt_phi + 1
end do

do i=1,Nsect    !  Boucle sur phi

  do j=1,Nsect  !  Boucle sur theta

    do k=1,nb_ligneCONTACT   !  Boucle sur les contacts
      cd        =   TAB_CONTACT(k)%icdent
      an        =   TAB_CONTACT(j)%ianent
      Lik(1:3)  =   TAB_SPHERE(cd)%center(1:3)-TAB_SPHERE(an)%center(1:3)
      norm      =   sqrt( Lik(1)**2 + Lik(2)**2 + Lik(3)**2 )
      Lik(1:3)  =   Lik(1:3) / norm
      
      phi       = acos( Lik(1) / sqrt(Lik(1)**2+Lik(2)**2) )
      if (Lik(1)**2+Lik(2)**2 == 0.D0 ) cycle
      
      
      
    end do                   !  Boucle sur les contacts
    
  end do        !  Boucle sur theta

end do          !  Boucle sur phi

end subroutine contact_orientation_3D


!==================================================================================================
!On ferme tout les fichiers
!=================================================================================================
subroutine close_all
implicit none
integer           :: i

print*,'---> fermeture des fichiers'

if (calcul_coordination                   == 1) close (100)
if (calcul_vitesse_moyenne                == 1) close (101)
if (calcul_compacity                      == 1) close (102)
if (calcul_qoverp                         == 1) close (103)
if (calcul_anisotropy_contact             == 1) close (104)
if (calcul_anisotropy_force               == 1) close (105)
if (calcul_anisotropy_branch              == 1) close (106)
if (calcul_walls_positions                == 1) close (108)
if (calcul_walls_forces                   == 1) close (109)
if (c_float_particles                     == 1) close (112)
if (c_mobilization                        == 1) close (126)
if (c_segregation_index                   == 1) close (128)
if (c_lacey_index                         == 1) close (129)
if (c_prop_nctc                           == 1) close (132)
if (c_stresses                            == 1) close (134)

deallocate(TAB_SPHERE)

print*,'--> Fin du postraitement <--'
stop

end subroutine close_all


!==============================================================================
! Reconstruction d'un BODIES.DAT, DOF.INI, DRV_DOF.DAT 
!==============================================================================



subroutine construct_bodies
implicit none
real(kind=8)                             :: dmean,H_max
integer                                  :: i,cpt,nbPart,nb_wallh,nb_wallb

if  (compteur_clout == quel_pas_de_temps) then
  dmean = 0.D0
  H_max = 0.D0
  do i=1,n_particles
    dmean = dmean + 2*TAB_SPHERE(i)%Rmax
    H_max = max(H_max,TAB_SPHERE(i)%center(3))
  end do
  dmean = dmean / real(n_particles,8)
  print*,dmean
  write(128,'(A72)')'! File BODIES                                                           ' 
  write(128,'(A72)')'!                                                                       ' 
  write(128,'(A72)')'! The symbol    $       preceeds a keyword used in scanning files.      ' 
  write(128,'(A72)')'!                                                                       '  
  write(128,'(A72)')'! The symbol    bdyty   stands for  body type data.                     ' 
  write(128,'(A72)')'! These data are distributed according to some species.                 ' 
  write(128,'(A72)')'!                                                                       '  
  write(128,'(A72)')'! the specy     blmty   stands for  bulk element type data ,            ' 
  write(128,'(A72)')'! i.e. part or total bulk geometric description,                        ' 
  write(128,'(A72)')'! and bulk behaviour laws;                                              ' 
  write(128,'(A72)')'!                                                                       ' 
  write(128,'(A72)')'! the specy     nodty   stands for  node type data ,                    ' 
  write(128,'(A72)')'! i.e. degrees of freedom data;                                         ' 
  write(128,'(A72)')'!                                                                       ' 
  write(128,'(A72)')'! the specy     tacty   stands for  contactor type data ;               ' 
  write(128,'(A72)')'!                                                                       ' 
  write(128,'(A72)')'! the keyword   $$$$$$  ends a body record.                             ' 
  write(128,'(A72)')'!                                                                       ' 
  write(128,'(A6)')'      ' 

  do i=1,n_particles
    TAB_SPHERE(i)%behav = 'PLEXx'
    if ( TAB_SPHERE(i)%center(3) < 2*dmean ) TAB_SPHERE(i)%behav='wallb'
    if ( TAB_SPHERE(i)%center(3) > H_max - 2*dmean ) TAB_SPHERE(i)%behav='wallh'
  end do

  cpt = 0
  do i=1,n_particles
    if (TAB_SPHERE(i)%behav /= 'PLEXx')  cycle
    cpt = cpt + 1
    
    write(128,'(A6)')                                        '$bdyty'
    write(128,'(1X,A5,1X,i6)')                               'RBDY3',cpt
    write(128,'(A6)')                                        '$blmty'
    write(128,'(1X,A5,5X,i2,2X,A5,2X,A5,1X,A6,D14.7)')       'PLAIN',1,'behav','PLEXx','avrd=',0.D0
    write(128,'(29X,3(A5,D14.7,2X))')  'I1  =',TAB_SPHERE(i)%I1,'I2  =',TAB_SPHERE(i)%I2,'I3  =',TAB_SPHERE(i)%I3
    write(128,'(A6)')                                        '$nodty'
    write(128,'(1X,A5,5X,i2,16X,3(A5,D14.7,2X))')            'NO6xx',1,'coo1=',TAB_SPHERE(i)%center(1),&
                                                                       'coo2=',TAB_SPHERE(i)%center(2),&
                                                                       'coo3=',TAB_SPHERE(i)%center(3)
    write(128,'(29X,3(A5,D14.7,2X))') 'coo4=',0.D0,'coo5=',0.D0,'coo6=',0.D0
    write(128,'(A6)')                                        '$tacty'
    write(128,'(1X,A5,5X,i2,2X,A5,2X,A5,2X,A5,D14.7)')'SPHER',1,'color',TAB_SPHERE(i)%color,'byrd=',TAB_SPHERE(i)%Rmax
    write(128,'(A6)')                                        '$$$$$$'
  end do

  nbPart = cpt
  nb_wallh = 0
  nb_wallb = 0
  
  do i=1,n_particles
    if (TAB_SPHERE(i)%behav /= 'wallh')  cycle
    cpt = cpt + 1
    nb_wallh = nb_wallh + 1
    
    write(128,'(A6)')                                        '$bdyty'
    write(128,'(1X,A5,1X,i6)')                               'RBDY3',cpt
    write(128,'(A6)')                                        '$blmty'
    write(128,'(1X,A5,5X,i2,2X,A5,2X,A5,1X,A6,D14.7)')       'PLAIN',1,'behav',TAB_SPHERE(i)%behav,'avrd=',0.D0
    write(128,'(29X,3(A5,D14.7,2X))')  'I1  =',TAB_SPHERE(i)%I1,'I2  =',TAB_SPHERE(i)%I2,'I3  =',TAB_SPHERE(i)%I3
    write(128,'(A6)')                                        '$nodty'
    write(128,'(1X,A5,5X,i2,16X,3(A5,D14.7,2X))')            'NO6xx',1,'coo1=',TAB_SPHERE(i)%center(1),&
                                                                       'coo2=',TAB_SPHERE(i)%center(2),&
                                                                       'coo3=',TAB_SPHERE(i)%center(3)
    write(128,'(29X,3(A5,D14.7,2X))') 'coo4=',0.D0,'coo5=',0.D0,'coo6=',0.D0
    write(128,'(A6)')                                        '$tacty'
    write(128,'(1X,A5,5X,i2,2X,A5,2X,A5,2X,A5,D14.7)')'SPHER',1,'color','TATAh','byrd=',TAB_SPHERE(i)%Rmax
    write(128,'(A6)')                                        '$$$$$$'
  end do

  cpt = cpt + 1
  write(128,'(A6)')                                        '$bdyty'
  write(128,'(1X,A5,1X,i6)')                               'RBDY3',cpt
  write(128,'(A6)')                                        '$blmty'
  write(128,'(1X,A5,2X,I5,2X,A5,2X,A5,2(2X,A5,D14.7))') 'PLAIN',1,'behav','wallb','avrd=',0.D0,'gyrd=',0.D0
  write(128,'(A6)')                                        '$nodty'
  write(128,'(1X,A5,5X,i2,16X,3(A5,D14.7,2X))')            'NO6xx',1,'coo1=',0.5/2,&
                                                                     'coo2=',0.6/2,&
                                                                     'coo3=',0.D0
  write(128,'(29X,3(A5,D14.7,2X))') 'coo4=',0.D0,'coo5=',0.D0,'coo6=',0.D0
  write(128,'(A6)')                                        '$tacty'
  do i=1,n_particles
    if (TAB_SPHERE(i)%behav /= 'wallb')  cycle  
    nb_wallb = nb_wallb + 1
    write(128,'(1X,A5,5X,i2,2X,A5,2X,A5,2X,A5,D14.7)')'SPHEb',1,'color','TATAb','byrd=',TAB_SPHERE(i)%Rmax
    write(128,'(29X,3(A5,D14.7,2X))') 'coo1=',TAB_SPHERE(i)%center(1)-0.5/2,&
                                      'coo2=',TAB_SPHERE(i)%center(2)-0.6/2,&
                                      'coo3=',TAB_SPHERE(i)%center(3)-0.D0
  enddo
  write(128,'(A6)')                                        '$$$$$$'


  print*,'nbPart=', nbPart
  print*,'nb_wallh=',nb_wallh
  print*,'nb_wallb=',nb_wallb



  write(130,'(A6)')  '! DOF'
  write(130,'(A1)')  ' '
                   !12345678901234567890123456789012345678901234567
  write(130,'(A47)') '$steps      0                time= 0.0000000D+00'
  write(130,'(A1)')  ' '
  write(130,'(A72)')'!-----------------------------------------------------------------------'



  do i=1,nb_wallh
    write(129,'(A72)') '!-----------------------------------------------------------------------' 
    write(129,'(A1)')  ' '
    write(129,'(A6)')  '$bdyty'
    write(129,'(1X,A5,2X,I5)') 'RBDY3',i+nbPart
    write(129,'(A6)')  '$nodty'
    write(129,'(1X,A5,2X,I5)') 'NO6xx',1
    write(129,'(A103)')  &
           '$dofty        [CT......+......AMP..*..cos.(..OMEGA.*.time.+.PHI..)]...*...[RAMPI.....+.....RAMP.*.time]'

    write(129,'(1X,A5,2X,I5,6(1X,D14.7))') 'vlocy',1,0.  ,0.,0.,0.,1.,0.
    write(129,'(1X,A5,2X,I5,6(1X,D14.7))') 'vlocy',2,0.2 ,0.,0.,0.,1.,0.
    write(129,'(1X,A5,2X,I5,6(1X,D14.7))') 'vlocy',3,0.  ,0.,0.,0.,1.,0.
    write(129,'(1X,A5,2X,I5,6(1X,D14.7))') 'vlocy',4,0.  ,0.,0.,0.,1.,0.
    write(129,'(1X,A5,2X,I5,6(1X,D14.7))') 'vlocy',5,0.  ,0.,0.,0.,1.,0.
    write(129,'(1X,A5,2X,I5,6(1X,D14.7))') 'vlocy',6,0.  ,0.,0.,0.,1.,0.
    write(129,'(A6)')'$$$$$$'
  end do

  write(129,'(A72)') '!-----------------------------------------------------------------------' 
  write(129,'(A1)')  ' '
  write(129,'(A6)')  '$bdyty'
  write(129,'(1X,A5,2X,I5)') 'RBDY3',nb_wallh+nbPart + 1
  write(129,'(A6)')  '$nodty'
  write(129,'(1X,A5,2X,I5)') 'NO6xx',1
  write(129,'(A103)')  &
         '$dofty        [CT......+......AMP..*..cos.(..OMEGA.*.time.+.PHI..)]...*...[RAMPI.....+.....RAMP.*.time]'
  write(129,'(1X,A5,2X,I5,6(1X,D14.7))') 'vlocy',1,0.,0.,0.,0.,1.,0.
  write(129,'(1X,A5,2X,I5,6(1X,D14.7))') 'vlocy',2,0.,0.,0.,0.,1.,0.
  write(129,'(1X,A5,2X,I5,6(1X,D14.7))') 'force',3,0.,0.,0.,0.,1.,0.
  write(129,'(1X,A5,2X,I5,6(1X,D14.7))') 'vlocy',4,0.,0.,0.,0.,1.,0.
  write(129,'(1X,A5,2X,I5,6(1X,D14.7))') 'vlocy',5,0.,0.,0.,0.,1.,0.
  write(129,'(1X,A5,2X,I5,6(1X,D14.7))') 'vlocy',6,0.,0.,0.,0.,1.,0.
  write(129,'(A6)')'$$$$$$'
  

end if

end subroutine construct_bodies









!---------------------------------------------------
! subroutine pour calculer les vecteurs et valeurs
! propre d'une matrice N*N
!---------------------------------------------------
subroutine rg ( lda, n, a, wr, wi, matz, z, ierror )
!
!*******************************************************************************
!
!! RG finds the eigenvalues and eigenvectors of a real(kind=8) general matrix.
!
!
!  Modified:
!
!    01 May 2000
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of A and Z.
!    LDA must be at least N.
!
!    Input, integer N, the number of rows and columns of A.
!
!    Input/output, real(kind=8) A(LDA,N).
!
!    On input, A contains the N by N matrix whose eigenvalues and
!    eigenvectors are desired.
!
!    On output, A has been overwritten by other information.
!
!    Output, real(kind=8) WR(N), WI(N), contain the real(kind=8) and imaginary parts,
!    respectively, of the eigenvalues.  Complex conjugate
!    pairs of eigenvalues appear consecutively with the
!    eigenvalue having the positive imaginary part first.
!
!    Input, integer MATZ, zero if only eigenvalues are desired.
!    non-zero, for both eigenvalues and eigenvectors.
!
!    Output, real(kind=8) Z(LDA,N), contains the real(kind=8) and imaginary parts of
!    the eigenvectors if MATZ is not zero.  If the J-th eigenvalue
!    is real(kind=8), the J-th column of Z contains its eigenvector.  If the
!    J-th eigenvalue is complex with positive imaginary part, the
!    J-th and (J+1)-th columns of Z contain the real(kind=8) and
!    imaginary parts of its eigenvector.  The conjugate of this
!    vector is the eigenvector for the conjugate eigenvalue.
!
!    Output, integer IERROR, error flag.
!    0, no error.
!    nonzero, an error occurred.
!
implicit none
!
integer lda
integer n
integer nm
!
real(kind=8) a(lda,n)
real(kind=8) fv1(n)
integer ierror
integer is1
integer is2
integer iv1(n)
integer matz
real(kind=8) wi(n)
real(kind=8) wr(n)
real(kind=8) z(lda,n)
!
ierror = 0

if ( n > lda ) then
ierror = 10 * n
return
end if
!
!  Balance the matrix.
!
call balanc ( lda, n, a, is1, is2, fv1 )
!
!  Put the matrix into upper Hessenberg form.
!
call elmhes ( lda, n, is1, is2, a, iv1 )

if ( matz == 0 ) then

call hqr ( lda, n, is1, is2, a, wr, wi, ierror )

if ( ierror /= 0 ) then
  return
end if

else

call eltran ( lda, n, is1, is2, a, iv1, z )

call hqr2 ( lda, n, is1, is2, a, wr, wi, z, ierror )

if ( ierror /= 0 ) then
  return
end if

call balbak ( lda, n, is1, is2, fv1, n, z )

end if

end subroutine rg

subroutine balanc ( nm, n, a, low, igh, scale )
!
!*******************************************************************************
!
!! BALANC balances a real(kind=8) matrix before eigenvalue calculations.
!
!
!  Discussion:
!
!    This subroutine balances a real(kind=8) matrix and isolates eigenvalues
!    whenever possible.
!
!    Suppose that the principal submatrix in rows LOW through IGH
!    has been balanced, that P(J) denotes the index interchanged
!    with J during the permutation step, and that the elements
!    of the diagonal matrix used are denoted by D(I,J).  Then
!
!      SCALE(J) = P(J),    J = 1,...,LOW-1,
!               = D(J,J),  J = LOW,...,IGH,
!               = P(J)     J = IGH+1,...,N.
!
!    The order in which the interchanges are made is N to IGH+1,
!    then 1 to LOW-1.
!
!    Note that 1 is returned for LOW if IGH is zero formally.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Modified:
!
!    15 March 2001
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of A, which must
!    be at least N.
!
!    Input, integer N, the order of the matrix.
!
!    Input/output, real(kind=8) A(NM,N), the N by N matrix.  On output,
!    the matrix has been balanced.
!
!    Output, integer LOW, IGH, indicate that A(I,J) is equal to zero if
!    (1) I is greater than J and
!    (2) J=1,...,LOW-1 or I=IGH+1,...,N.
!
!    Output, real(kind=8) SCALE(N), contains information determining the
!    permutations and scaling factors used.
!
implicit none
!
integer nm
integer n
!
real(kind=8) a(nm,n)
real(kind=8) b2
real(kind=8) c
real(kind=8) f
real(kind=8) g
integer i
integer iexc
integer igh
integer j
integer k
integer l
integer low
integer m
logical noconv
real(kind=8) r
real(kind=8), parameter :: radix = 16.0D+00
real(kind=8) s
real(kind=8) scale(n)
!
iexc = 0
j = 0
m = 0

b2 = radix**2
k = 1
l = n
go to 100

20 continue

scale(m) = j

if ( j /= m ) then

do i = 1, l
  call r_swap ( a(i,j), a(i,m) )
end do

do i = k, n
  call r_swap ( a(j,i), a(m,i) )
end do

end if

continue

if ( iexc == 2 ) go to 130
!
!  Search for rows isolating an eigenvalue and push them down.
!
continue

if ( l == 1 ) then
low = k
igh = l
return
end if

l = l - 1

100 continue

do j = l, 1, -1

 do i = 1, l
   if ( i /= j ) then
     if ( a(j,i) /= 0.0D+00 ) then
       go to 120
     end if
   end if
 end do

 m = l
 iexc = 1
 go to 20

120  continue

end do

go to 140
!
!  Search for columns isolating an eigenvalue and push them left.
!
130 continue

k = k + 1

140 continue

do j = k, l

do i = k, l
  if ( i /= j ) then
    if ( a(i,j) /= 0.0D+00 ) then
      go to 170
    end if
  end if
end do

m = k
iexc = 2
go to 20

170 continue

end do
!
!  Balance the submatrix in rows K to L.
!
scale(k:l) = 1.0D+00
!
!  Iterative loop for norm reduction.
!
noconv = .true.

do while ( noconv )

noconv = .false.

do i = k, l

  c = 0.0D+00
  r = 0.0D+00

  do j = k, l
    if ( j /= i ) then
      c = c + abs ( a(j,i) )
      r = r + abs ( a(i,j) )
    end if
  end do
!
!  Guard against zero C or R due to underflow.
!
  if ( c /= 0.0D+00 .and. r /= 0.0D+00 ) then

    g = r / radix
    f = 1.0D+00
    s = c + r

    do while ( c < g )
      f = f * radix
      c = c * b2
    end do

    g = r * radix

    do while ( c >= g )
      f = f / radix
      c = c / b2
    end do
!
!  Balance.
!
    if ( ( c + r ) / f < 0.95D+00 * s ) then

      g = 1.0D+00 / f
      scale(i) = scale(i) * f
      noconv = .true.

      a(i,k:n) = a(i,k:n) * g
      a(1:l,i) = a(1:l,i) * f

    end if

  end if

end do

end do

low = k
igh = l


end subroutine balanc
subroutine balbak ( lda, n, low, igh, scale, m, z )
!
!*******************************************************************************
!
!! BALBAK back transforms eigenvectors to undo the effect of BALANC.
!
!
!  Discussion:
!
!    This subroutine forms the eigenvectors of a real(kind=8) general
!    matrix by back transforming those of the corresponding
!    balanced matrix determined by BALANC.
!
!  Reference:
!
!    Parlett and Reinsch,
!    Numerische Mathematik,
!    Volume 13, pages 293-304, 1969.
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Modified:
!
!    18 February 2001
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of Z.
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer LOW, IGH, column indices determined by BALANC.
!
!    Input, real(kind=8) SCALE(N), contains information determining
!    the permutations and scaling factors used by BALANC.
!
!    Input, integer M, the number of columns of Z to be back-transformed.
!
!    Input/output, real(kind=8) Z(LDA,M), contains the real(kind=8) and imaginary parts
!    of the eigenvectors, which, on return, have been back-transformed.
!
implicit none
!
integer lda
integer m
integer n
!
integer i
integer igh
integer ii
integer j
integer k
integer low
real(kind=8) scale(n)
real(kind=8) z(lda,m)
!
if ( m <= 0 ) then
return
end if

if ( igh /= low ) then
do i = low, igh
  z(i,1:m) = scale(i) * z(i,1:m)
end do
end if

do ii = 1, n

i = ii

if ( i < low .or. i > igh ) then

  if ( i < low ) then
    i = low - ii
  end if

  k = int ( scale(i) )

  if ( k /= i ) then

    do j = 1, m
      call r_swap ( z(i,j), z(k,j) )
    end do

  end if

end if

end do


end subroutine balbak

subroutine elmhes ( lda, n, low, igh, a, jint )
!
!*******************************************************************************
!
!! ELMHES reduces all, or a portion of a matrix, to upper Hessenberg form.
!
!
!  Discussion:
!
!    The routine uses stabilized elementary similarity transformations.
!
!  Modified:
!
!    15 March 2001
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of A.
!
!    Input, integer N, the number of rows and columns in A.
!
!    Input, integer LOW, IGH.  If subroutine BALANC was called
!    before ELMHES, then it will set LOW and IGH.  If BALANC
!    was not called, then the user should set LOW = 1 and IGH=N.
!
!    Input/output, real(kind=8) A(LDA,N).
!
!    On input, A contains the matrix to be transformed.
!
!    On output, A contains the Hessenberg matrix.  The multipliers
!    which were used in the reduction are stored in the
!    remaining triangle under the Hessenberg matrix.
!
!    Output, integer JINT(IGH), contains information on the rows
!    and columns interchanged in the reduction.
!    Only elements LOW through IGH are used.
!
implicit none
!
integer igh
integer lda
integer n
!
real(kind=8) a(lda,n)
integer i
integer j
integer jint(igh)
integer low
integer m
real(kind=8) x
real(kind=8) y
!
do m = low+1, igh-1
!
!  Look for the largest element in the column A(J,M-1), where
!  J goes from M to IGH.   Store the row number as I.
!
x = 0.0D+00
i = m

do j = m, igh

  if ( abs ( a(j,m-1) ) > abs ( x ) ) then
    x = a(j,m-1)
    i = j
  end if

end do

jint(m) = i
!
!  If I is not M, interchange rows and columns I and M of A.
!
if ( i /= m ) then

  do j = m-1, n
    call r_swap ( a(i,j), a(m,j) )
  end do

  do j = 1, igh
    call r_swap ( a(j,i), a(j,m) )
  end do

end if

if ( x /= 0.0D+00 ) then

  do i = m+1, igh
    y = a(i,m-1)

    if ( y /= 0.0D+00 ) then

      y = y / x
      a(i,m-1) = y
      a(i,m:n) = a(i,m:n) - y * a(m,m:n)

      do j = 1, igh
        a(j,m) = a(j,m) + y * a(j,i)
      end do

    end if

  end do

end if

end do


end subroutine elmhes
subroutine eltran ( lda, n, low, igh, a, jint, z )
!
!*******************************************************************************
!
!! ELTRAN accumulates transformations used by ELMHES.
!
!
!  Modified:
!
!    15 March 2001
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of A and Z.
!
!    Input, integer N, the number of rows and columns in A.
!
!    Input, integer LOW, IGH.  If BALANC was called
!    before ELMHES, then it will set LOW and IGH.  If BALANC
!    was not called, then the user should set LOW = 1 and IGH=N.
!
!    Input, integer A(LDA,N), contains, in its lower triangle,
!    the multipliers used by ELMHES in the reduction.
!
!    Input, integer JINT(IGH), contains information on the rows
!    and columns interchanged in the reduction.
!
!    Output, real(kind=8) Z(LDA,N), contains the transformation matrix
!    produced in the reduction by ELMHES.
!
implicit none
!
integer igh
integer lda
integer n
!
real(kind=8) a(lda,igh)
integer i
integer j
integer jint(igh)
integer low
integer mm
integer mp
real(kind=8) z(lda,n)
!
!  Initialize Z to the identity matrix.
!
call rmat_identity ( lda, n, z )

do mm = 1, igh-low-1

mp = igh - mm

do i = mp+1, igh
  z(i,mp) = a(i,mp-1)
end do

i = jint(mp)

if ( i /= mp ) then

  do j = mp, igh
    z(mp,j) = z(i,j)
    z(i,j) = 0.0D+00
  end do

  z(i,mp) = 1.0D+00

end if

end do


end subroutine eltran

subroutine hqr ( lda, n, low, igh, h, wr, wi, ierror )
!
!*******************************************************************************
!
!! HQR finds the eigenvalues of a real(kind=8) upper Hessenberg matrix by the QR method.
!
!
!  Modified:
!
!    15 March 2001
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of H.  LDA must be
!    at least N.
!
!    Input, integer N, the number of rows and columns in the matrix.
!
!    Input, integer LOW, IGH, column indices determined by
!    BALANC.  If BALANC is not used, set LOW = 1, IGH=N.
!
!    Input/output, real(kind=8) H(LDA,N).
!
!    On input, H contains the upper Hessenberg matrix.  Information
!    about the transformations used in the reduction to Hessenberg
!    form by ELMHES or ORTHES, if performed, is stored
!    in the remaining triangle under the Hessenberg matrix.
!
!    On output, the information that was in H has been destroyed.
!
!    Output, real(kind=8) WR(N), WI(N), contain the real(kind=8) and imaginary parts,
!    respectively, of the eigenvalues.  The eigenvalues
!    are unordered except that complex conjugate pairs
!    of values appear consecutively with the eigenvalue
!    having the positive imaginary part first.  If an
!    error exit is made, the eigenvalues should be correct
!    for indices IERROR+1 through N.
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    J, if the limit of 30*N iterations is exhausted
!    while the J-th eigenvalue is being sought.
!
implicit none
!
integer lda
integer n
!
integer en
integer enm2
real(kind=8) h(lda,n)
real(kind=8) hnorm
integer i
integer ierror
integer igh
integer itn
integer its
integer j
integer k
integer l
integer ll
integer low
integer m
integer mm
integer na
logical notlas
real(kind=8) p
real(kind=8) q
real(kind=8) r
real(kind=8) s
real(kind=8) t
real(kind=8) tst1
real(kind=8) tst2
real(kind=8) w
real(kind=8) wi(n)
real(kind=8) wr(n)
real(kind=8) x
real(kind=8) y
real(kind=8) zz
!
ierror = 0
!
!  Compute the norm of the upper Hessenberg matrix.
!
hnorm = 0.0D+00
do i = 1, n
do j = max ( i-1, 1 ), n
  hnorm = hnorm + abs ( h(i,j) )
end do
end do
!
!  Store roots isolated by BALANC.
!
do i = 1, n

if (i < low .or. i > igh ) then
  wr(i) = h(i,i)
  wi(i) = 0.0D+00
end if

end do

en = igh
t = 0.0D+00
itn = 60 * n
!
!  Search for next eigenvalues.
!
60 continue

if ( en < low ) then
return
end if

its = 0
na = en - 1
enm2 = na - 1
!
!  Look for single small sub-diagonal element.
!
70 continue

do ll = low, en

l = en + low - ll
if ( l == low ) then
  exit
end if

s = abs ( h(l-1,l-1) ) + abs ( h(l,l) )
if ( s == 0.0D+00 ) then
  s = hnorm
end if

tst1 = s
tst2 = tst1 + abs ( h(l,l-1) )
if ( tst2 == tst1 ) then
  exit
end if

end do
!
!  Form shift.
!
x = h(en,en)

if ( l == en ) then
wr(en) = x + t
wi(en) = 0.0D+00
en = na
go to 60
end if

y = h(na,na)
w = h(en,na) * h(na,en)
if ( l == na) then
go to 280
end if

if ( itn == 0 ) then
ierror = en
return
end if
!
!  Form exceptional shift.
!
if ( its == 10 .or. its == 20 ) then

t = t + x

do i = low, en
  h(i,i) = h(i,i) - x
end do

s = abs ( h(en,na) ) + abs ( h(na,enm2) )
x = 0.75D+00 * s
y = x
w = -0.4375D+00 * s * s

end if

its = its + 1
itn = itn - 1
!
!  Look for two consecutive small sub-diagonal elements.
!
do mm = l, enm2

m = enm2 + l - mm
zz = h(m,m)
r = x - zz
s = y - zz
p = (r * s - w) / h(m+1,m) + h(m,m+1)
q = h(m+1,m+1) - zz - r - s
r = h(m+2,m+1)
s = abs ( p ) + abs ( q ) + abs ( r )
p = p / s
q = q / s
r = r / s
if ( m == l ) then
  exit
end if

tst1 = abs ( p ) * ( abs ( h(m-1,m-1) ) + abs ( zz ) + abs ( h(m+1,m+1) ) )
tst2 = tst1 + abs ( h(m,m-1) ) * ( abs ( q ) + abs ( r ) )
if ( tst2 == tst1 ) then
  exit
end if

end do

do i = m+2, en
h(i,i-2) = 0.0D+00
if ( i /= m+2 ) then
  h(i,i-3) = 0.0D+00
end if
end do
!
!  Double QR step involving rows l to en and columns m to en.
!
do k = m, na

notlas = k /= na

if ( k /= m ) then

  p = h(k,k-1)
  q = h(k+1,k-1)
  r = 0.0D+00
  if ( notlas ) then
    r = h(k+2,k-1)
  end if

  x = abs ( p ) + abs ( q ) + abs ( r )
  if ( x == 0.0D+00 ) then
    cycle
  end if

  p = p / x
  q = q / x
  r = r / x
end if

s = sign ( sqrt ( p*p+q*q+r*r), p )

if ( k /= m ) then
  h(k,k-1) = -s * x
else
  if ( l /= m ) then
    h(k,k-1) = -h(k,k-1)
  end if
end if

p = p + s
x = p / s
y = q / s
zz = r / s
q = q / p
r = r / p

if ( .not. notlas ) then
!
!  Row modification.
!
  do j = k, n
    p = h(k,j) + q * h(k+1,j)
    h(k,j) = h(k,j) - p * x
    h(k+1,j) = h(k+1,j) - p * y
  end do

  j = min ( en, k+3 )
!
!  Column modification.
!
  do i = 1, j
    p = x * h(i,k) + y * h(i,k+1)
    h(i,k) = h(i,k) - p
    h(i,k+1) = h(i,k+1) - p * q
  end do

else
!
!  Row modification.
!
  do j = k, n
    p = h(k,j) + q * h(k+1,j) + r * h(k+2,j)
    h(k,j) = h(k,j) - p * x
    h(k+1,j) = h(k+1,j) - p * y
    h(k+2,j) = h(k+2,j) - p * zz
  end do

  j = min(en,k+3)
!
!  Column modification.
!
  do i = 1, j
    p = x * h(i,k) + y * h(i,k+1) + zz * h(i,k+2)
    h(i,k) = h(i,k) - p
    h(i,k+1) = h(i,k+1) - p * q
    h(i,k+2) = h(i,k+2) - p * r
  end do

end if

end do

go to 70
!
!  Two roots found.
!
280 continue

p = ( y - x ) / 2.0D+00
q = p * p + w
zz = sqrt ( abs ( q ) )
x = x + t
!
!  Real(Kind=8) pair.
!
if (q >= 0.0D+00 ) then
zz = p + sign(zz,p)
wr(na) = x + zz
if ( zz == 0.0D+00 ) then
  wr(en) = wr(na)
else
  wr(en) = x - w / zz
end if
wi(na) = 0.0D+00
wi(en) = 0.0D+00
!
!  Complex pair.
!
else
wr(na) = x + p
wr(en) = x + p
wi(na) = zz
wi(en) = -zz
end if
!
!  Deduct the two eigenvalues we have found from the total to
!  be found, and proceed.
!
en = enm2

go to 60
end subroutine hqr
subroutine hqr2 ( lda, n, low, igh, h, wr, wi, z, ierror )
!
!*******************************************************************************
!
!! HQR2 finds the eigenvalues and eigenvectors of a real(kind=8) upper Hessenberg matrix by the QR method.
!
!
!  Modified:
!
!    15 March 2001
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of H and Z.
!    LDA must be at least N.
!
!    Input, integer N, the number of rows and columns in the
!    matrix H.
!
!    Input, integer LOW, IGH, column indices determined by
!    BALANC.  If BALANC is not used, set LOW = 1, IGH=N.
!
!    Input/output, real(kind=8) H(LDA,N).
!    On input, H contains the upper Hessenberg matrix.
!    On output, the information that was in H has been destroyed.
!
!    Output, real(kind=8) WR(N), WI(N), contain the real(kind=8) and imaginary parts,
!    respectively, of the eigenvalues.  The eigenvalues
!    are unordered except that complex conjugate pairs
!    of values appear consecutively with the eigenvalue
!    having the positive imaginary part first.  If an
!    error exit is made, the eigenvalues should be correct
!    for indices IERROR+1 through N.
!
!    Input/output, real(kind=8) Z(LDA,N).
!
!    On input, Z contains the transformation matrix produced by
!    ELTRAN after the reduction by ELMHES if performed.  If the
!    eigenvectors of the Hessenberg matrix are desired, Z must
!    contain the identity matrix.
!
!    On output, Z contains the real(kind=8) and imaginary parts of the
!    eigenvectors.  If the I-th eigenvalue is real(kind=8), the I-th column
!    of Z contains its eigenvector.  If the I-th eigenvalue is complex
!    with positive imaginary part, the I-th and (I+1)-th
!    columns of Z contain the real(kind=8) and imaginary parts of its
!    eigenvector.  The eigenvectors are unnormalized.  If an
!    error exit is made, none of the eigenvectors has been found.
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    J, if the limit of 30*N iterations is exhausted
!    while the J-th eigenvalue is being sought.
!
implicit none
!
integer lda
integer n
!
integer en
integer enm2
real(kind=8) h(lda,n)
real(kind=8) hnorm
integer i
integer ierror
integer igh
integer ii
integer itn
integer its
integer j
integer jj
integer k
integer l
integer ll
integer low
integer m
integer mm
integer na
integer nn
logical notlas
real(kind=8) p
real(kind=8) q
real(kind=8) r
real(kind=8) ra
real(kind=8) s
real(kind=8) sa
real(kind=8) t
real(kind=8) temp
real(kind=8) tst1
real(kind=8) tst2
real(kind=8) vi
real(kind=8) vr
real(kind=8) w
real(kind=8) wi(n)
real(kind=8) wr(n)
real(kind=8) x
real(kind=8) y
real(kind=8) z(lda,n)
real(kind=8) zz
!
ierror = 0
!
!  Compute the norm of the upper Hessenberg matrix.
!
hnorm = 0.0D+00
do i = 1, n
do j = max ( i-1, 1 ), n
  hnorm = hnorm + abs ( h(i,j) )
end do
end do
!
!  Store roots isolated by BALANC.
!
do i = 1, n

if ( i < low .or. i > igh ) then
  wr(i) = h(i,i)
  wi(i) = 0.0D+00
end if

end do

en = igh
t = 0.0D+00
itn = 60 * n
!
!  Search for next eigenvalues
!
60 continue

if ( en < low ) then
go to 340
end if

its = 0
na = int(en-1)
enm2 = na - 1
!
!  Look for single small sub-diagonal element.
!
70 continue

do ll = low, en

l = en + low - ll
if ( l == low ) then
  exit
end if

s = abs ( h(l-1,l-1) ) + abs ( h(l,l) )
if ( s == 0.0D+00 ) then
  s = hnorm
end if

tst1 = s
tst2 = tst1 + abs ( h(l,l-1) )
if ( tst2 == tst1 ) then
  exit
end if

end do
!
!  Form shift.
!
x = h(en,en)
if ( l == en ) then
go to 270
end if

y = h(na,na)
w = h(en,na) * h(na,en)
if ( l == na ) then
go to 280
end if

if ( itn == 0 ) then
ierror = en
return
end if
!
!  Form exceptional shift.
!
if ( its == 10 .or. its == 20 ) then

t = t + x

do i = low, en
  h(i,i) = h(i,i) - x
end do

s = abs ( h(en,na) ) + abs ( h(na,enm2) )
x = 0.75D+00 * s
y = x
w = -0.4375D+00 * s * s

end if

its = its + 1
itn = itn - 1
!
!  Look for two consecutive small sub-diagonal elements.
!
do mm = l, enm2

m = enm2 + l - mm
zz = h(m,m)
r = x - zz
s = y - zz
p = (r * s - w) / h(m+1,m) + h(m,m+1)
q = h(m+1,m+1) - zz - r - s
r = h(m+2,m+1)
s = abs ( p ) + abs ( q ) + abs ( r )
p = p / s
q = q / s
r = r / s

if ( m == l) then
  exit
end if

tst1 = abs ( p ) * ( abs ( h(m-1,m-1) ) + abs ( zz ) + abs ( h(m+1,m+1)))
tst2 = tst1 + abs ( h(m,m-1) ) * ( abs ( q ) + abs ( r ) )
if ( tst2 == tst1 ) then
  exit
end if

end do

do i = m+2, en
h(i,i-2) = 0.0D+00
if ( i /= m+2 ) then
  h(i,i-3) = 0.0D+00
end if
end do
!
!  Double QR step involving rows L to EN and columns M to EN.
!
do k = m, na

 notlas = k /= na

 if ( k /= m ) then
   p = h(k,k-1)
   q = h(k+1,k-1)
   r = 0.0D+00
   if ( notlas ) r = h(k+2,k-1)
   x = abs ( p ) + abs ( q ) + abs ( r )
   if ( x == 0.0D+00 ) then
     cycle
   end if
   p = p / x
   q = q / x
   r = r / x
 end if

 s = sign ( sqrt ( p*p + q*q + r*r ), p )

 if ( k /= m ) then
   h(k,k-1) = -s * x
 else
   if ( l /= m ) h(k,k-1) = -h(k,k-1)
 end if

 p = p + s
 x = p / s
 y = q / s
 zz = r / s
 q = q / p
 r = r / p

 if ( .not. notlas ) then
!
!  Row modification.
!
 do j = k, n
    p = h(k,j) + q * h(k+1,j)
    h(k,j) = h(k,j) - p * x
    h(k+1,j) = h(k+1,j) - p * y
 end do

 j = min ( en, k + 3 )
!
!  Column modification.
!
 do i = 1, j
    p = x * h(i,k) + y * h(i,k+1)
    h(i,k) = h(i,k) - p
    h(i,k+1) = h(i,k+1) - p * q
 end do
!
!  Accumulate transformations.
!
 do i = low, igh
   p = x * z(i,k) + y * z(i,k+1)
   z(i,k) = z(i,k) - p
   z(i,k+1) = z(i,k+1) - p * q
 end do

 else
!
!  Row modification.
!
 do j = k, n
    p = h(k,j) + q * h(k+1,j) + r * h(k+2,j)
    h(k,j) = h(k,j) - p * x
    h(k+1,j) = h(k+1,j) - p * y
    h(k+2,j) = h(k+2,j) - p * zz
 end do

 j=min(en,k+3)
!
!  Column modification.
!
 do i = 1, j
   p = x * h(i,k) + y * h(i,k+1) + zz * h(i,k+2)
   h(i,k) = h(i,k) - p
   h(i,k+1) = h(i,k+1) - p * q
   h(i,k+2) = h(i,k+2) - p * r
 end do
!
!  Accumulate transformations.
!
 do i = low, igh
    p = x * z(i,k) + y * z(i,k+1) + zz * z(i,k+2)
    z(i,k) = z(i,k) - p
    z(i,k+1) = z(i,k+1) - p * q
    z(i,k+2) = z(i,k+2) - p * r
 end do

end if

end do

go to 70
!
!  One root found.
!
270 continue

h(en,en) = x + t
wr(en) = h(en,en)
wi(en) = 0.0D+00
en = na
go to 60
!
!  Two roots found.
!
280 p = ( y - x ) / 2.0D+00
q = p * p + w
zz = sqrt ( abs ( q ) )
h(en,en) = x + t
x = h(en,en)
h(na,na) = y + t

if ( q < 0.0D+00 ) then
go to 320
end if
!
!  Real(Kind=8) pair.
!
zz = p + sign(zz,p)
wr(na) = x + zz
wr(en)=wr(na)
if ( zz /= 0.0D+00 ) wr(en) = x - w / zz
wi(na) = 0.0D+00
wi(en) = 0.0D+00
x = h(en,na)
s = abs ( x ) + abs ( zz )
p = x / s
q = zz / s
r = sqrt(p*p+q*q)
p = p / r
q = q / r
!
!  Row modification.
!
do j = na, n
 zz = h(na,j)
 h(na,j) = q * zz + p * h(en,j)
 h(en,j) = q * h(en,j) - p * zz
end do
!
!  Column modification.
!
do i = 1, en
 zz = h(i,na)
 h(i,na) = q * zz + p * h(i,en)
 h(i,en) = q * h(i,en) - p * zz
end do
!
!  Accumulate transformations.
!
do i = low, igh
 zz = z(i,na)
 z(i,na) = q * zz + p * z(i,en)
 z(i,en) = q * z(i,en) - p * zz
end do

go to 330
!
!  Complex pair.
!
320 continue

wr(na) = x + p
wr(en) = x + p
wi(na) = zz
wi(en) = -zz

330 continue
en = enm2
go to 60
!
!  All roots found.  Backsubstitute to find vectors of upper
!  triangular form.
!
340 continue

if ( hnorm == 0.0D+00 ) then
return
end if

do nn = 1, n
 en = n + 1 - nn
 p=wr(en)
 q=wi(en)
 na = en - 1

 if ( q < 0.0D+00 ) then
   go to 710
 end if

 if ( q > 0.0D+00 ) then
   go to 800
 end if
!
!  Real(Kind=8) vector.
!
 m = en
 h(en,en) = 1.0D+00
 if ( na == 0 ) then
   go to 800
 end if

 do ii = 1, na
    i = en - ii
    w = h(i,i) - p
    r = 0.0D+00

    do j = m, en
      r = r + h(i,j) * h(j,en)
    end do

    if ( wi(i) < 0.0D+00 ) then
      zz = w
      s = r
      go to 700
    end if

    m = i

    if ( wi(i) == 0.0D+00 ) then

      t = w

      if ( t == 0.0D+00 ) then
        tst1 = hnorm
        t = tst1
632       continue
        t = 0.01D+00 * t
        tst2 = hnorm + t
        if ( tst2 > tst1 ) then
          go to 632
        end if
      end if

      h(i,en) = -r / t
      go to 680

    end if
!
!  Solve real(kind=8) equations.
!
    x = h(i,i+1)
    y = h(i+1,i)
    q = (wr(i) - p) * (wr(i) - p) + wi(i) * wi(i)
    t = (x * s - zz * r) / q
    h(i,en) = t

    if ( abs ( x ) > abs ( zz ) ) then
      h(i+1,en) = (-r - w * t) / x
    else
      h(i+1,en) = (-s - y * t) / zz
    end if
!
!  Overflow control.
!
680       continue

    t = abs ( h(i,en) )
    if ( t == 0.0D+00 ) then
      go to 700
    end if
    tst1 = t
    tst2 = tst1 + 1.0D+00 / tst1

    if ( tst2 <= tst1 ) then

      do j = i, en
        h(j,en)=h(j,en)/t
      end do

    end if

700       continue
 end do
!
!  End real(kind=8) vector.
!
 go to 800
!
!  Complex vector.
!
710    continue

m = na
!
!  Last vector component chosen imaginary so that
!  eigenvector matrix is triangular.
!
 if ( abs ( h(en,na) ) > abs ( h(na,en) ) ) then
   h(na,na) = q / h(en,na)
   h(na,en) = -(h(en,en) - p) / h(en,na)
 else
   temp = 0.0D+00
   call cdiv(temp,-h(na,en),h(na,na)-p,q,h(na,na),h(na,en))
 end if

 h(en,na) = 0.0D+00
 h(en,en) = 1.0D+00
 enm2 = na - 1

 do ii = 1, enm2

    i = na - ii
    w = h(i,i) - p
    ra = 0.0D+00
    sa = 0.0D+00

    do j = m, en
       ra = ra + h(i,j) * h(j,na)
       sa = sa + h(i,j) * h(j,en)
    end do

    if ( wi(i) < 0.0D+00 ) then
      zz=w
      r = ra
      s = sa
      go to 795
    end if

    m = i

    if ( wi(i) == 0.0D+00 ) then
      call cdiv ( -ra, -sa, w, q, h(i,na), h(i,en) )
      go to 790
    end if
!
!  Solve complex equations.
!
    x = h(i,i+1)
    y = h(i+1,i)
    vr = (wr(i) - p) * (wr(i) - p) + wi(i) * wi(i) - q * q
    vi = (wr(i) - p) * 2.0D+00 * q

    if ( vr == 0.0D+00 .and. vi == 0.0D+00 ) then

      tst1 = hnorm * ( abs ( w ) + abs ( q ) + abs ( x ) + abs ( y ) &
        + abs ( zz ) )
      vr = tst1

783     continue

      vr = 0.01D+00 * vr
      tst2 = tst1 + vr
      if ( tst2 > tst1) then
        go to 783
      end if

    end if

    call cdiv(x*r-zz*ra+q*sa,x*s-zz*sa-q*ra,vr,vi,h(i,na),h(i,en))

    if ( abs ( x ) > abs ( zz ) + abs ( q ) ) then
      h(i+1,na) = (-ra - w * h(i,na) + q * h(i,en)) / x
      h(i+1,en) = (-sa - w * h(i,en) - q * h(i,na)) / x
    else
      call cdiv ( -r-y*h(i,na), -s-y*h(i,en), zz, q, h(i+1,na), h(i+1,en) )
    end if
!
!  Overflow control
!
790     continue

    t = max ( abs ( h(i,na) ), abs ( h(i,en) ) )

    if ( t /= 0.0D+00 ) then

      tst1 = t
      tst2 = tst1 + 1.0D+00 / tst1
      if ( tst2 <= tst1 ) then
        h(i:en,na) = h(i:en,na) / t
        h(i:en,en) = h(i:en,en) / t
      end if

    end if

795     continue

 end do
!
!  End complex vector
!
800    continue
end do
!
!  End back substitution.
!
!  Vectors of isolated roots
!
do i = 1, n

if ( i < low .or. i > igh ) then
  z(i,i:n) = h(i,i:n)
end if

end do
!
!  Multiply by transformation matrix to give
!  vectors of original full matrix.
!
do jj = low, n

j = n + low - jj
m = min ( j, igh )

do i = low, igh

  zz = 0.0D+00
  do  k = low, m
    zz = zz + z(i,k) * h(k,j)
  end do

  z(i,j) = zz

end do

end do


end subroutine hqr2


subroutine cdiv ( ar, ai, br, bi, cr, ci )
!
!*******************************************************************************
!
!! CDIV carries out complex division.
!
!
!  Discussion:
!
!    CDIV computes:
!
!      (CR,CI) = (AR,AI) / (BR,BI)
!
!    using real(kind=8) arithmetic.
!
!  Modified:
!
!    15 March 2001
!
!  Parameters:
!
!    Input, real(kind=8) AR, AI, the real(kind=8) and imaginary parts of the
!    number to be divided.
!
!    Input, real(kind=8) BR, BI, the real(kind=8) and imaginary parts of the divisor.
!
!    Output, real(kind=8) CR, CI, the real(kind=8) and imaginary parts of the resultant.
!
implicit none
!
real(kind=8) ai
real(kind=8) ais
real(kind=8) ar
real(kind=8) ars
real(kind=8) bi
real(kind=8) bis
real(kind=8) br
real(kind=8) brs
real(kind=8) ci
real(kind=8) cr
real(kind=8) s
!
s = abs ( br ) + abs ( bi )
ars = ar / s
ais = ai / s
brs = br / s
bis = bi / s
s = brs**2 + bis**2
cr = ( ars * brs + ais * bis ) / s
ci = ( ais * brs - ars * bis ) / s


end subroutine cdiv


subroutine r_swap ( x, y )
!
!*******************************************************************************
!
!! R_SWAP switches two real(kind=8) values.
!
!
!  Modified:
!
!    30 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real(kind=8) X, Y.  On output, the values of X and
!    Y have been interchanged.
!
implicit none
!
real(kind=8) x
real(kind=8) y
real(kind=8) z
!
z = x
x = y
y = z


end subroutine r_swap

subroutine rmat_identity ( lda, n, a )
!
!*******************************************************************************
!
!! RMAT_IDENTITY sets the square matrix A to the identity.
!
!
!  Modified:
!
!    24 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of A.
!
!    Input, integer N, the order of A.
!
!    Output, real(kind=8) A(LDA,N), the matrix which has been
!    set to the identity.
!
implicit none
!
integer lda
integer n
!
real(kind=8) a(lda,n)
integer i
integer j
!
do i = 1, n
do j = 1, n
  if ( i == j ) then
    a(i,j) = 1.0D+00
  else
    a(i,j) = 0.0D+00
  end if
end do
end do


end subroutine rmat_identity

end program post3D_spheres
