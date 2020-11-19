!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! POST-PROCESSOR !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! POLYGONS!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program post2D_polyg

implicit none

type  ::  T_CORPS
   integer                                  ::  num,lieu,nb_vertex,ctc, group
   real(kind=8)                             ::  radius,Area,ax1,ax2, s11, s12, s21, s22
   real(kind=8)                             ::  p, q, qoverp, crit, Area_temp
   real(kind=8),dimension(3)                ::  center,center_ref
   real(kind=8),dimension(3)                ::  V
   real(kind=8),allocatable,dimension(:,:)  ::  vertex,vertex_ref
   character(len=5)                         ::  behav
end type T_CORPS

type  ::  T_CONTACT
   integer                                  ::  icdent,ianent, cdver, sect, adj
   integer                                  ::  type, status_points
   real(kind=8),dimension(3)                ::  n
   real(kind=8),dimension(3)                ::  t
   real(kind=8),dimension(3)                ::  coor_ctc
   real(kind=8)                             ::  rn,rt,rs
   real(kind=8)                             ::  vn,vt,vs
   real(kind=8)                             ::  gap0, gap, cd_len, tang_disp
   character(len=5)                         ::  nature, status, law
end type T_CONTACT

type  ::  T_DOF
   integer                                  ::  num
   real(kind=8),dimension(3)                ::  X
   real(kind=8),dimension(3)                ::  V
end type T_DOF

type(T_CORPS),allocatable,dimension(:)      ::  TAB_POLYG,TAB_PLAN,TAB_POLYG_INI, TAB_POINT
type(T_CONTACT),allocatable,dimension(:)    ::  TAB_CONTACT,TAB_CONTACT_persiste
type(T_CONTACT),allocatable,dimension(:)    ::  TAB_CONTACT_POLYG


!DEFINITIONS Vloc-Dof-Bodies
!================================================
integer                                     ::  n_walls=0,fin=0,n_particles=0, n_r_grain=0, n_pt_points=0
integer                                     ::  compteur_clout=0,all_step=0
logical                                     ::  first_bodies=.true., fin_post=.false.
logical                                     ::  first_over_all = .true.
integer                                     ::  step,nb_ligneCONTACT,nb_ligneCONTACT_POLYG
real(kind=8)                                ::  time

!DEFINITIONS BOITE
!================================================
real(kind=8),dimension(2,4)                ::  tab_sommet_box
real(kind=8)                               ::  long,haut,height,width

!COMMANDES D APPELLE
!================================================
character(len=30)                          ::  command
logical                                    ::  first_list_contact=.true.
integer                                    ::  nb_ligneCONTACT_first=0
integer                                    ::  calcul_coordination=0,calcul_vitesse_moyenne=0,&
                                               calcul_compacity=0,calcul_qoverp=0,calcul_anisotropy_contact=0,&
                                               calcul_anisotropy_force=0,calcul_persiste_contact, &
                                               calcul_granulo=0, walls_position = 0, walls_f = 0, contact_dir = 0
integer                                    ::  calcul_construct_bodies=0, &
                                               c_draw = 0, calcul_anisotropy_branch=0, c_n_ctc_probability = 0, &
                                               c_f_list=0, option_cohe=0, c_fail_mode=0, c_list_interact=0, &
                                               c_brc_dir=0, c_frc_dir=0, c_sign_aniso=0, c_clean_tresca=0, &
                                               c_multi_part_frac=0, c_avg_shape_ratio=0, c_sign_aniso_l=0, &
                                               c_list_interact_l = 0

!Variables DavidC.
real(kind=8), allocatable, dimension(:,:)  ::  tab_temp_disp
real(kind=8), allocatable, dimension(:)    ::  tab_temp_rot
integer                                    ::  l, nb_part_invi, nb_part_visi, nb_vi_ini
integer                                    ::  N_PLPLx, nc_poly, n_ctc_cohe
integer                                    ::  N_PLJCx, nc_wall, N_PTPT2, nc_ptpt
logical                                    ::  copy_input=.false.
integer, dimension(:), allocatable         ::  t_failure
integer, dimension(:,:), allocatable       ::  cohe_couples
real(kind=8), dimension(:,:), allocatable  ::  l_forces_c

!PARAMETRE
!================================================
real(kind=8),parameter                     ::  pi=3.141593D0

!================================================
!Lecture du fichier de commande
!================================================
print*,'-------------------------------------------'
print*,'!                                         !'
print*,'!      Module de Post-traitement 2D       !'
print*,'!         from Emilien Azema Code         !'
print*,'!                 David C.                !'
print*,'-------------------------------------------'
print*,'Lecture des commandes'
open(unit=1,file='POST_INPUT.DAT',status='old')

do
  
  read(1,'(A30)') command
  print*,command

  if (command=='MULTIPART FRAG               :') then
    c_multi_part_frac = 1
    cycle
  end if
  
  if (command=='COPY IMPUT                   :') then
    copy_input = .true.
    cycle
  end if
  
  if (command=='VLOC_DOF_BODIES              :') then
    read(1,*) n_walls
    read(1,*) n_particles
    if (c_multi_part_frac == 1) then
      read(1,*) n_r_grain
      read(1,*) n_pt_points
    end if 
    cycle
  end if
  
  if (command=='COORDINATION                 :') then
    calcul_coordination=1
    open (unit=100,file='./POSTPRO/COORDINATION.DAT',status='replace')
    cycle
  end if
  
  if (command=='GRAIN SIZE DISTRIB           :') then
    calcul_granulo=1
    cycle
  end if
  
  if (command=='COMPACITY                    :') then
    calcul_compacity=1
    open (unit=102,file='./POSTPRO/COMPACITY.DAT',status='replace')
    cycle
  end if
  
  if (command=='QoverP                       :') then
    calcul_qoverp=1
    open (unit=103,file='./POSTPRO/QoverP.DAT',status='replace')
    cycle
  end if
  
  
  if (command=='VITESSE MOYENNE              :') then
    calcul_vitesse_moyenne=1
    open (unit=101,file='./POSTPRO/VITESSE_MOYENNE.DAT',status='replace')
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
  
  if (command=='PERSISTE CONTACT             :') then
    calcul_persiste_contact=1
    ! OJO ya estoy utilizando el puerto 106!!!!
    !open (unit=106,file='./POSTPRO/PERSISTE_CONTACT.DAT',status='replace')
    cycle
  end if
  
  if (command=='WALLS POSITION               :') then
    walls_position=1
    open (unit=108,file='./POSTPRO/WALLSPOS.DAT',status='replace')
    cycle
  end if
  
  if (command=='WALLS FORCES                 :') then
    walls_f=1
    open (unit=109,file='./POSTPRO/WALLSFORCE.DAT',status='replace')
    cycle
  end if
  
  if (command=='CONTACTS ORIENTATION         :') then
    contact_dir=1
    cycle
    ! Uses port 110
  end if
  
  if (command=='CONTACT NUMBER PROBABILITY   :') then
    c_n_ctc_probability = 1
    open (unit=111,file='./POSTPRO/CTC_PROBABILITY.DAT',status='replace')
    cycle
  end if

  if (command=='CONTACT FORCES LIST          :') then
    c_f_list=1
    cycle
    ! Uses port 112
  end if

  if (command=='FAILURE MODE                 :') then
    c_fail_mode=1
    ! Uses port 113
    open (unit=113,file='./POSTPRO/FAIL_MODE.DAT',status='replace')
    cycle
  end if

  if (command=='LIST INTERACTIONS            :') then
    c_list_interact=1
    ! Uses port 114
    cycle
  end if

  if (command=='BRANCH DISTRIBUTION          :') then
    c_brc_dir = 1
    ! Uses port 115
    cycle
  end if

  if (command=='FORCES DISTRIBUTION          :') then
    c_frc_dir = 1
    ! Uses port 116
    cycle
  end if

  if (command=='SIGNED ANISOTROPIES          :') then
    c_sign_aniso = 1
    ! Uses port 117
    open (unit=117,file='./POSTPRO/SGN_ANISO.DAT',status='replace')
    cycle
  end if

  if (command=='CLEAN TRESCA                 :') then
    c_clean_tresca = 1
    ! Uses port 118 and 123
    open (unit=123,file='./POSTPRO/Vloc_Rloc.INI',status='replace')
    cycle
  end if

  if (command=='SHAPE RATIO                  :') then
    c_avg_shape_ratio = 1
    ! Uses port 124
    open (unit=124,file='./POSTPRO/AVG_SHAPE.DAT',status='replace')
    cycle
  end if

  if (command=='SIGNED ANISOTROPIES L        :') then
    c_sign_aniso_l = 1
    ! Uses port 125
    open (unit=125,file='./POSTPRO/SGN_ANISO_L.DAT',status='replace')
    cycle
  end if

  if (command=='LIST INTERACTIONS L          :') then
    c_list_interact_l=1
    ! Uses port 126
    cycle
  end if

  if (command=='IS COHESIVE                  :') then
    option_cohe=1
    cycle
  end if

  if (command=='DRAW                         :') then
    c_draw = 1
    ! Uses port 121
    cycle
  end if
  
  if (command=='END                          :') exit
end do
close(1)

!==============================================================================
!Appel des differente ressources lues dans le fichier de commande
!==============================================================================
  do
    
    if (copy_input .and. first_over_all) then
      call system('cp DATBOX/DOF.INI OUTBOX/DOF.OUT.0')
      call system('cp DATBOX/Vloc_Rloc.INI OUTBOX/Vloc_Rloc.OUT.0')
    endif
    
    call read_Vloc_dof_bodies
    
    call test_position_contact
    
    ! Calling the Subroutines
    if (calcul_coordination         == 1)  call nb_coordination
    if (walls_position              == 1)  call wallsposition
    if (walls_f                     == 1)  call walls_force
    if (calcul_granulo              == 1)  call granulometry
    if (contact_dir                 == 1)  call contact_direction
    if (calcul_vitesse_moyenne      == 1)  call vitesse_moyenne
    if (calcul_compacity            == 1)  call compacity
    if (calcul_qoverp               == 1)  call qoverp
    if (calcul_anisotropy_contact   == 1)  call anisotropy_contact
    if (calcul_anisotropy_force     == 1)  call anisotropy_force
    if (calcul_anisotropy_branch    == 1)  call anisotropy_branch
    if (calcul_persiste_contact     == 1)  call persiste_contact
    if (c_n_ctc_probability         == 1)  call ctc_probability
    if (c_f_list                    == 1)  call f_list
    if (c_fail_mode                 == 1)  call failure_mode
    if (c_list_interact             == 1)  call list_interact
    if (c_draw                      == 1)  call draw
    if (c_brc_dir                   == 1)  call branch_dir
    if (c_frc_dir                   == 1)  call forces_dir
    if (c_sign_aniso                == 1)  call signed_anisotropy
    if (c_clean_tresca              == 1)  call clean_tresca
    if (c_sign_aniso_l              == 1)  call signed_anisotropy_l
    if (c_list_interact_l           == 1)  call list_interact_l

    ! One step subroutines 
    if (first_over_all) then
      if (c_avg_shape_ratio           == 1)  call avg_shape_ratio
    end if

    if (fin_post) then
      call close_all
      print*,'--> End of post-processing <--'
      exit
    end if
    first_over_all = .false.
  end do

!==============================================================================
! Subroutines
!==============================================================================

 contains

!==============================================================================
! Lecture des donnees
!==============================================================================
SUBROUTINE read_Vloc_dof_bodies
  
  implicit none
  
  integer                       ::  i,j,k,err,num_part,nb_vertex,gen, n_groups
  character(len=22)             ::  clout_DOF
  character(len=28)             ::  clout_Vloc
  character(len=20)             ::  clout_Bodies
  character(len=5)              ::  status, i_law
  character(len=6)              ::  text
  character(len=13)             ::  text2
  real(kind=8)                  ::  center1,center2,center3
  real(kind=8)                  ::  vertex_coor1,vertex_coor2,vertex_coor3
  real(kind=8),dimension(3)     ::  g_gravity
  real(kind=8)                  ::  coor1,coor2,coor3
  real(kind=8)                  ::  Vx,Vy,Vr
  real(kind=8)                  ::  radius,ax1,ax2
  integer                       ::  icdent,ianent,cpt_type, l_stat_p,l_ver, l_seg, l_adj
  real(kind=8)                  ::  n1,n2,n3
  real(kind=8)                  ::  rn,rt,rs,rn_temp,rt_temp,rs_temp, gap_temp
  real(kind=8)                  ::  vln,vlt,vls
  real(kind=8)                  ::  coor_ctc1,coor_ctc2,coor_ctc3
  real(kind=8)                  ::  coor_ctc1_temp,coor_ctc2_temp,coor_ctc3_temp
  real(kind=8)                  ::  gap_c, gap0, cd_len, tang_disp
  integer                       ::  n_simple,n_double
  logical                       ::  first_double=.true.
  
  clout_DOF    = './OUTBOX/DOF.OUT.    '
  clout_Vloc   = './OUTBOX/Vloc_Rloc.OUT.    '
  clout_Bodies = './OUTBOX/BODIES.OUT'
  
  ! Setting the index of the output files
  compteur_clout = compteur_clout+1
  
  if (copy_input .and. first_over_all) then
    compteur_clout = compteur_clout - 1
  endif
  
  if (compteur_clout<10) then
    WRITE(clout_DOF(18:19),'(I1)')  compteur_clout
    WRITE(clout_Vloc(24:25),'(I1)')  compteur_clout
  else if ((compteur_clout>=10) .and. (compteur_clout<100)) then
    WRITE(clout_DOF(18:20),'(I2)')   compteur_clout
    WRITE(clout_Vloc(24:26),'(I2)')  compteur_clout
  else if ((compteur_clout>=100) .and. (compteur_clout<1000)) then
    WRITE(clout_DOF(18:21),'(I3)')   compteur_clout
    WRITE(clout_Vloc(24:27),'(I3)')  compteur_clout
  else if ((compteur_clout>=1000) .and. (compteur_clout<10000)) then
    WRITE(clout_DOF(18:22),'(I4)')   compteur_clout
    WRITE(clout_Vloc(24:28),'(I4)')  compteur_clout
  end if
  
  ! Opening the BODIES.OUT file for its reading
  if (compteur_clout == 1) then
    open(unit=2,file=clout_Bodies, iostat=err ,status='old')
    if (err/=0) then
      fin_post=.true.
    endif
  elseif (compteur_clout == 0 .and. copy_input) then
    open(unit=2,file=clout_Bodies, iostat=err ,status='old')
    if (err/=0) then
      fin_post=.true.
    endif
  endif
  
  if (first_bodies) then
    allocate(TAB_POLYG(n_particles))
    allocate(TAB_PLAN(n_walls))
    allocate(TAB_POINT(n_pt_points))
  endif
  
  ! Counter
  i=0
  
  ! Reading the BODIES.OUT
  if((fin_post .neqv. .true.) .and. (first_bodies .eqv. .true.) ) then
    print*,'--->',clout_Bodies
    do
      read(2,'(A6)') text
      if (text == '$bdyty') then
        i = i +1
        if (i<n_particles+1) then
          TAB_POLYG(i)%num=i
          read(2,*) 
          read(2,*) 
          read(2,'(35X,D14.7)') radius
          read(2,*)
          
          TAB_POLYG(i)%behav = 'PLEXx'
          read(2,'(29X, 3(5x,D14.7,2X))') center1,center2,center3
          TAB_POLYG(i)%center_ref(1)=center1
          TAB_POLYG(i)%center_ref(2)=center2
          TAB_POLYG(i)%center_ref(3)=center3
          read(2,*)
          read(2,'(40X,i6,9X,D14.7)') nb_vertex
          TAB_POLYG(i)%radius=radius
          TAB_POLYG(i)%nb_vertex=nb_vertex
          
          ! Allocating the space for the vertices coordinates
          allocate(TAB_POLYG(i)%vertex_ref(3,nb_vertex))
          allocate(TAB_POLYG(i)%vertex(3,nb_vertex))
          
          ! Reading the vertices coordinates
          do j=1,nb_vertex
            read(2,'(29X, 2(5x,D14.7,2X))') vertex_coor1,vertex_coor2
            TAB_POLYG(i)%vertex_ref(1,j) = vertex_coor1
            TAB_POLYG(i)%vertex_ref(2,j) = vertex_coor2
            TAB_POLYG(i)%vertex_ref(3,j) = 0
            g_gravity(1) = g_gravity(1) + vertex_coor1
            g_gravity(2) = g_gravity(2) + vertex_coor2
            g_gravity(3) = g_gravity(3) + 0
          end do
          read(2,*)
          read(2,*)
          
          ! Computing the grain area
          TAB_POLYG(i)%Area = 0.D0
          do j=1,nb_vertex-2
            TAB_POLYG(i)%Area = TAB_POLYG(i)%Area + &
                       0.5d0*(((TAB_POLYG(i)%vertex_ref(1,j+1) - TAB_POLYG(i)%vertex_ref(1,1)) &
                              *(TAB_POLYG(i)%vertex_ref(2,j+2) - TAB_POLYG(i)%vertex_ref(2,1))) &
                             -((TAB_POLYG(i)%vertex_ref(2,j+1) - TAB_POLYG(i)%vertex_ref(2,1)) &
                              *(TAB_POLYG(i)%vertex_ref(1,j+2) - TAB_POLYG(i)%vertex_ref(1,1))))
          enddo
        else if (i<n_particles + n_walls + 1) then
          TAB_PLAN(i-n_particles)%num=i
          read(2,*)
          read(2,*)
          read(2,*)
          read(2,*)
          read(2,'(29X, 3(5x,D14.7,2X))') center1,center2,center3
          TAB_PLAN(i-n_particles)%center_ref(1)=center1
          TAB_PLAN(i-n_particles)%center_ref(2)=center2
          TAB_PLAN(i-n_particles)%center_ref(3)=center3
          read(2,*)
          read(2,'(29X, 2(5x,D14.7,2X))') ax1,ax2
          TAB_PLAN(i-n_particles)%ax1=ax1
          TAB_PLAN(i-n_particles)%ax2=ax2
          TAB_PLAN(i-n_particles)%center(:)=0.D0
          TAB_PLAN(i-n_particles)%behav = 'WALL_'
          
          ! Similar to poly, i'm going to create vertices (8 vertices in each wall)
          allocate(TAB_PLAN(i-n_particles)%vertex_ref(3,4))
          allocate(TAB_PLAN(i-n_particles)%vertex(3,4))
          
          ! Vertices
          TAB_PLAN(i-n_particles)%vertex_ref(1,1) =  ax1
          TAB_PLAN(i-n_particles)%vertex_ref(2,1) =  ax2
          TAB_PLAN(i-n_particles)%vertex_ref(3,1) =  0.D0
          
          TAB_PLAN(i-n_particles)%vertex_ref(1,2) = (- ax1)
          TAB_PLAN(i-n_particles)%vertex_ref(2,2) =  ax2
          TAB_PLAN(i-n_particles)%vertex_ref(3,2) =  0.D0
          
          TAB_PLAN(i-n_particles)%vertex_ref(1,3) = (- ax1)
          TAB_PLAN(i-n_particles)%vertex_ref(2,3) = (- ax2)
          TAB_PLAN(i-n_particles)%vertex_ref(3,3) =  0.D0
          
          TAB_PLAN(i-n_particles)%vertex_ref(1,4) =  ax1
          TAB_PLAN(i-n_particles)%vertex_ref(2,4) = (- ax2)
          TAB_PLAN(i-n_particles)%vertex_ref(3,4) =  0.D0
          
        else 
          TAB_POINT(i-n_particles-n_walls)%num=i
          read(2,*) 
          read(2,*) 
          read(2,'(35X,D14.7)') radius
          read(2,*)
          
          TAB_POINT(i-n_particles-n_walls)%behav = 'PLEXx'
          read(2,'(29X, 3(5x,D14.7,2X))') center1,center2,center3
          TAB_POINT(i-n_particles-n_walls)%center_ref(1)=center1
          TAB_POINT(i-n_particles-n_walls)%center_ref(2)=center2
          TAB_POINT(i-n_particles-n_walls)%center_ref(3)=center3
        end if
        
        if (i==n_particles+n_walls+n_pt_points) exit
      end if
      
      first_bodies = .false.
    end do
  end if
  close(2)
  
  ! Reading the DOF.OUT.#
  open(unit=5,file=clout_DOF,iostat=err,status='old')
  if (err/=0) then
    fin_post=.true.
  else
    fin_post=.false.
    i = 0
    read(5,*)
    read(5,*)
    read(5,*)
    read(5,'(7X,i8,19X,D14.7)') step,time
    print*,'--->',clout_DOF
    do
      read(5,'(A6)') text
      if (text == '$bdyty') then
        i = i+1
        read(5,'(6X,i9)') num_part
        
        if (i<n_particles+1) then
          read(5,*)
          read(5,'(29X, 3(5x,D14.7,2X))') coor1,coor2,coor3
          read(5,'(29X, 3(5x,D14.7,2X))') Vx,Vy,Vr
          
          ! Appling the relative displacement
          TAB_POLYG(i)%center(1)=coor1+TAB_POLYG(i)%center_ref(1)
          TAB_POLYG(i)%center(2)=coor2+TAB_POLYG(i)%center_ref(2)
          
          ! The third component is a rotation not a displacement
          !TAB_POLYG(i)%center(3)=coor3+TAB_POLYG(i)%center_ref(3)
          
          ! Applying the relative rotation 
          do j=1,TAB_POLYG(i)%nb_vertex
            TAB_POLYG(i)%vertex(1,j) = TAB_POLYG(i)%vertex_ref(1,j) * cos(coor3) - &
                                       TAB_POLYG(i)%vertex_ref(2,j) * sin(coor3) + &
                                       TAB_POLYG(i)%center(1)
            TAB_POLYG(i)%vertex(2,j) = TAB_POLYG(i)%vertex_ref(1,j) * sin(coor3) + &
                                       TAB_POLYG(i)%vertex_ref(2,j) * cos(coor3) + &
                                       TAB_POLYG(i)%center(2)
            ! Initializing the third component (it's necessary for vtk)
            TAB_POLYG(i)%vertex(3,j) = 0.D0
          end do
          
          TAB_POLYG(i)%V(1) = Vx
          TAB_POLYG(i)%V(2) = Vy
          TAB_POLYG(i)%V(3) = Vr
          
        else if (i<n_particles+n_walls+1) then
          read(5,*)
          read(5,'(29X, 3(5x,D14.7,2X))') coor1,coor2,coor3
          read(5,'(29X, 3(5x,D14.7,2X))') Vx,Vy,Vr
          TAB_PLAN(i-n_particles)%center(1)=coor1+TAB_PLAN(i-n_particles)%center_ref(1)
          TAB_PLAN(i-n_particles)%center(2)=coor2+TAB_PLAN(i-n_particles)%center_ref(2)
          
          ! The third componennt is a rotation not a displacement CHECK THIS
          !TAB_PLAN(i-n_particles)%center(3)=coor3+TAB_PLAN(i-n_particles)%center_ref(3)
          TAB_PLAN(i-n_particles)%V(1) = Vx
          TAB_PLAN(i-n_particles)%V(2) = Vy
          TAB_PLAN(i-n_particles)%V(3) = Vr
          
          ! Applying the relative rotation
          do j=1, 4
            TAB_PLAN(i-n_particles)%vertex(1,j) =  TAB_PLAN(i-n_particles)%vertex_ref(1,j) * cos(coor3) - &
                                                   TAB_PLAN(i-n_particles)%vertex_ref(2,j) * sin(coor3) + &
                                                   TAB_PLAN(i-n_particles)%center(1)
            TAB_PLAN(i-n_particles)%vertex(2,j) =  TAB_PLAN(i-n_particles)%vertex_ref(1,j) * sin(coor3) + &
                                                   TAB_PLAN(i-n_particles)%vertex_ref(2,j) * cos(coor3) + &
                                                   TAB_PLAN(i-n_particles)%center(2)
            ! Initializing the third component (it's necessary for vtk)
            TAB_PLAN(i-n_particles)%vertex(3,j) = 0.D0
          end do
        else
          read(5,*)
          read(5,'(29X, 3(5x,D14.7,2X))') coor1,coor2,coor3
          read(5,'(29X, 3(5x,D14.7,2X))') Vx,Vy,Vr
          TAB_POINT(i-n_particles-n_walls)%center(1)=coor1+TAB_POINT(i-n_particles-n_walls)%center_ref(1)
          TAB_POINT(i-n_particles-n_walls)%center(2)=coor2+TAB_POINT(i-n_particles-n_walls)%center_ref(2)
          TAB_POINT(i-n_particles-n_walls)%V(1) = Vx
          TAB_POINT(i-n_particles-n_walls)%V(2) = Vy
          TAB_POINT(i-n_particles-n_walls)%V(3) = Vr
        end if
      end if
      
      if (i==n_particles+n_walls+n_pt_points) exit
    end do
  end if
  
  N_PLJCx = 0
  N_PLPLx = 0
  N_PTPT2 = 0
  
  ! Reading the Vloc_Rloc.OUT.#
  ! This is a pre-reading for counting the number of contacts
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
      if (text2 == '$icdan  PLPLx') then
        N_PLPLx = N_PLPLx +1
      end if
      if (text2 == '$icdan  PLJCx') then
        N_PLJCx = N_PLJCx +1
      end if
      if (text2 == '$icdan  PTPT2') then
        N_PTPT2 = N_PTPT2 +1
      end if
    end do

    close(4)
    
    ! Allocating the corresponding space for the number of contacts
    if (allocated(TAB_CONTACT)) deallocate(TAB_CONTACT)
    allocate(TAB_CONTACT(N_PLPLx+N_PLJCx+N_PTPT2))
    
    ! Storing the contacts information
    ! In the array the particle-particle contacts are stored first
    nb_ligneCONTACT=0
    nc_poly = 0
    nc_wall = 0
    nc_ptpt = 0
    
    open(unit=4,file=clout_Vloc,status='old')
    read(4,*) 
    read(4,*) 
    read(4,*) 
    read(4,*) 
    read(4,*) 
    read(4,*) 
    
    do
      read(4,'(A13)',iostat=err) text2
      
      if (text2 == '$icdan  PLPLx') then
        nc_poly = nc_poly + 1
        read(4,*)
        read(4,'(8X,i5,23X,i5,2X,A5,9X,i5,23X,i5,2X,A5)') icdent, l_ver, i_law, ianent,l_seg, status
        read(4,'(29X, 3(5x,D14.7,2X))') rt,rn,rs
        read(4,'(29X, 3(5x,D14.7,2X))') vlt,vln,vls
        read(4,'(29X, 1(5x,D14.7,2X))') gap_c
        read(4,'(29X, 3(5x,D14.7,2X))') n1,n2,n3
        read(4,'(29X, 3(5x,D14.7,2X))') coor_ctc1,coor_ctc2,coor_ctc3

        !if (status(1:1) == "W") then
        read(4,*)
        read(4,'(3(1X,D14.7))') gap0, cd_len, tang_disp
        !end if

        nb_ligneCONTACT=nb_ligneCONTACT+1

        TAB_CONTACT(nc_poly)%icdent=icdent
        TAB_CONTACT(nc_poly)%ianent=ianent
        TAB_CONTACT(nc_poly)%status=status
        TAB_CONTACT(nc_poly)%sect  =l_seg
        TAB_CONTACT(nc_poly)%cdver =l_ver

        ! Giving some value to the status to average the TAB_CONTACT_POLYG somehow
        if (status == 'noctc') then
          TAB_CONTACT(nc_poly)%status_points = 3
        else if(status(1:3) == 'sli') then
          TAB_CONTACT(nc_poly)%status_points = 2
        else if(status == 'stick') then
          TAB_CONTACT(nc_poly)%status_points = 1
        else if(status == 'Wnctc') then
          TAB_CONTACT(nc_poly)%status_points = -1
        else if(status(1:3) == 'Wsl') then
          TAB_CONTACT(nc_poly)%status_points = -2
        else if(status == 'Wstck') then
          TAB_CONTACT(nc_poly)%status_points = -3
        end if

        !TAB_CONTACT(nc_poly)%law=i_law
        !if (status(1:1) == "W") then
        TAB_CONTACT(nc_poly)%gap0=gap0
        TAB_CONTACT(nc_poly)%gap = gap_c
        TAB_CONTACT(nc_poly)%cd_len=cd_len
        TAB_CONTACT(nc_poly)%tang_disp = tang_disp
        !else
        !  TAB_CONTACT(nc_poly)%gap0=-99.
        !  TAB_CONTACT(nc_poly)%gap = gap_c
        !  TAB_CONTACT(nc_poly)%cd_len=-99.
        !  TAB_CONTACT(nc_poly)%tang_disp = -99.
        !end if

        TAB_CONTACT(nc_poly)%n(1)=n1
        TAB_CONTACT(nc_poly)%n(2)=n2
        TAB_CONTACT(nc_poly)%n(3)=n3
        TAB_CONTACT(nc_poly)%t(1)=n2
        TAB_CONTACT(nc_poly)%t(2)=-n1
        TAB_CONTACT(nc_poly)%t(3)=0.D0
        TAB_CONTACT(nc_poly)%coor_ctc(1)=coor_ctc1
        TAB_CONTACT(nc_poly)%coor_ctc(2)=coor_ctc2
        TAB_CONTACT(nc_poly)%coor_ctc(3)=coor_ctc3
        TAB_CONTACT(nc_poly)%rn=rn
        TAB_CONTACT(nc_poly)%rt=rt
        TAB_CONTACT(nc_poly)%rs=rs
        TAB_CONTACT(nc_poly)%nature='PLPLx'
        TAB_CONTACT(nc_poly)%type = 0
        TAB_CONTACT(nc_poly)%vn=vln
        TAB_CONTACT(nc_poly)%vt=vlt
        TAB_CONTACT(nc_poly)%vs=vls
        TAB_CONTACT(nc_poly)%law=i_law
      end if
      
      if (text2 == '$icdan  PLJCx') then
        nc_wall = nc_wall + 1
        read(4,*)
        read(4,'(8X,i5,30X,A5,9X,i5,30X,A5,2X,i5)') icdent, i_law, ianent,status,l_adj
        read(4,'(29X, 3(5x,D14.7,2X))') rt,rn,rs
        read(4,'(29X, 3(5x,D14.7,2X))') vlt,vln,vls
        read(4,'(29X, 1(5x,D14.7,2X))') gap_c
        read(4,'(29X, 3(5x,D14.7,2X))') n1,n2,n3
        read(4,'(29X, 3(5x,D14.7,2X))') coor_ctc1,coor_ctc2,coor_ctc3
        
        nb_ligneCONTACT=nb_ligneCONTACT+1
        
        TAB_CONTACT(N_PLPLx+nc_wall)%icdent=icdent
        TAB_CONTACT(N_PLPLx+nc_wall)%ianent=ianent
        TAB_CONTACT(N_PLPLx+nc_wall)%status=status
        TAB_CONTACT(N_PLPLx+nc_wall)%status_points=99
        TAB_CONTACT(N_PLPLx+nc_wall)%n(1)=n1
        TAB_CONTACT(N_PLPLx+nc_wall)%n(2)=n2
        TAB_CONTACT(N_PLPLx+nc_wall)%n(3)=n3
        TAB_CONTACT(N_PLPLx+nc_wall)%t(1)=n2
        TAB_CONTACT(N_PLPLx+nc_wall)%t(2)=-n1
        TAB_CONTACT(N_PLPLx+nc_wall)%t(3)=0.D0
        TAB_CONTACT(N_PLPLx+nc_wall)%coor_ctc(1)=coor_ctc1
        TAB_CONTACT(N_PLPLx+nc_wall)%coor_ctc(2)=coor_ctc2
        TAB_CONTACT(N_PLPLx+nc_wall)%coor_ctc(3)=coor_ctc3
        TAB_CONTACT(N_PLPLx+nc_wall)%rn=rn
        TAB_CONTACT(N_PLPLx+nc_wall)%rt=rt
        TAB_CONTACT(N_PLPLx+nc_wall)%rs=rs
        TAB_CONTACT(N_PLPLx+nc_wall)%nature='PLJCx'
        TAB_CONTACT(N_PLPLx+nc_wall)%law=i_law
        TAB_CONTACT(N_PLPLx+nc_wall)%type  = 0
        TAB_CONTACT(N_PLPLx+nc_wall)%vn=vln
        TAB_CONTACT(N_PLPLx+nc_wall)%vt=vlt
        TAB_CONTACT(N_PLPLx+nc_wall)%vs=vls
        TAB_CONTACT(N_PLPLx+nc_wall)%adj=l_adj

        ! Default initialization
        TAB_CONTACT(N_PLPLx+nc_wall)%gap0=-99.
        TAB_CONTACT(N_PLPLx+nc_wall)%gap = gap_c
        TAB_CONTACT(N_PLPLx+nc_wall)%cd_len=-99.
        TAB_CONTACT(N_PLPLx+nc_wall)%tang_disp = -99.
      end if

      if (text2 == '$icdan  PTPT2') then
        nc_ptpt = nc_ptpt + 1
        read(4,*)
        read(4,'(8X,i5,16X,A5,9X,i5,16X,A5,2X,i5)') icdent, i_law, ianent,status,l_adj
        read(4,'(29X, 3(5x,D14.7,2X))') rt,rn,rs
        read(4,'(29X, 3(5x,D14.7,2X))') vlt,vln,vls
        read(4,'(29X, 1(5x,D14.7,2X))') gap_c
        read(4,'(29X, 3(5x,D14.7,2X))') n1,n2,n3
        read(4,'(29X, 3(5x,D14.7,2X))') coor_ctc1,coor_ctc2,coor_ctc3
        
        nb_ligneCONTACT=nb_ligneCONTACT+1
        
        TAB_CONTACT(N_PLPLx+N_PLJCx+nc_ptpt)%icdent=icdent
        TAB_CONTACT(N_PLPLx+N_PLJCx+nc_ptpt)%ianent=ianent
        TAB_CONTACT(N_PLPLx+N_PLJCx+nc_ptpt)%status=status
        TAB_CONTACT(N_PLPLx+N_PLJCx+nc_ptpt)%status_points=99
        TAB_CONTACT(N_PLPLx+N_PLJCx+nc_ptpt)%n(1)=n1
        TAB_CONTACT(N_PLPLx+N_PLJCx+nc_ptpt)%n(2)=n2
        TAB_CONTACT(N_PLPLx+N_PLJCx+nc_ptpt)%n(3)=n3
        TAB_CONTACT(N_PLPLx+N_PLJCx+nc_ptpt)%t(1)=n2
        TAB_CONTACT(N_PLPLx+N_PLJCx+nc_ptpt)%t(2)=-n1
        TAB_CONTACT(N_PLPLx+N_PLJCx+nc_ptpt)%t(3)=0.D0
        TAB_CONTACT(N_PLPLx+N_PLJCx+nc_ptpt)%coor_ctc(1)=coor_ctc1
        TAB_CONTACT(N_PLPLx+N_PLJCx+nc_ptpt)%coor_ctc(2)=coor_ctc2
        TAB_CONTACT(N_PLPLx+N_PLJCx+nc_ptpt)%coor_ctc(3)=coor_ctc3
        TAB_CONTACT(N_PLPLx+N_PLJCx+nc_ptpt)%rn=rn
        TAB_CONTACT(N_PLPLx+N_PLJCx+nc_ptpt)%rt=rt
        TAB_CONTACT(N_PLPLx+N_PLJCx+nc_ptpt)%rs=rs
        TAB_CONTACT(N_PLPLx+N_PLJCx+nc_ptpt)%nature='PTPT2'
        TAB_CONTACT(N_PLPLx+N_PLJCx+nc_ptpt)%law=i_law
        TAB_CONTACT(N_PLPLx+N_PLJCx+nc_ptpt)%type  = 0
        TAB_CONTACT(N_PLPLx+N_PLJCx+nc_ptpt)%gap = gap_c
        TAB_CONTACT(N_PLPLx+N_PLJCx+nc_ptpt)%vn=vln
        TAB_CONTACT(N_PLPLx+N_PLJCx+nc_ptpt)%vt=vlt
        TAB_CONTACT(N_PLPLx+N_PLJCx+nc_ptpt)%vs=vls
        TAB_CONTACT(N_PLPLx+N_PLJCx+nc_ptpt)%adj=l_adj

        read(4,*)
        read(4,'(1X,D14.7)') gap0

        TAB_CONTACT(N_PLPLx+nc_wall+nc_ptpt)%gap0=gap0
        
        ! Default initialization
        TAB_CONTACT(N_PLPLx+nc_wall+nc_ptpt)%cd_len=-99.
        TAB_CONTACT(N_PLPLx+nc_wall+nc_ptpt)%tang_disp = -99.
      end if

      if (nc_poly+nc_wall+nc_ptpt == N_PLPLx + N_PLJCx + N_PTPT2) then
        exit
      endif
    end do

    print*,'Step:', step
    print*,'Time:', time
  end if
  
  ! It's necessary to define the type of contact (i.e., simple or double)
  do i=1,nb_ligneCONTACT-1
    if (TAB_CONTACT(i)%type>0) cycle
    if ((TAB_CONTACT(i)%icdent == TAB_CONTACT(i+1)%icdent) .and. &
        (TAB_CONTACT(i)%ianent == TAB_CONTACT(i+1)%ianent)) then
      TAB_CONTACT(i)%type=2
      TAB_CONTACT(i+1)%type=2
    else
      TAB_CONTACT(i)%type=1
    end if
  end do
  
  ! Storing the contacts
  n_simple = 0
  n_double = 0
  first_double=.true.
  
  do i=1,nb_ligneCONTACT
    if (TAB_CONTACT(i)%type==1) n_simple = n_simple + 1
    if (TAB_CONTACT(i)%type==2) n_double = n_double + 1
  end do
  
  if (allocated(TAB_CONTACT_POLYG)) deallocate(TAB_CONTACT_POLYG)
  allocate(TAB_CONTACT_POLYG(n_simple+n_double))
  
  nb_ligneCONTACT_POLYG = 0

  print*,'---> Checking doble contacts'
  
  do i=1,nb_ligneCONTACT
    if (TAB_CONTACT(i)%type==1) then
      nb_ligneCONTACT_POLYG                             = nb_ligneCONTACT_POLYG +1
      TAB_CONTACT_POLYG(nb_ligneCONTACT_POLYG)%icdent   = TAB_CONTACT(i)%icdent
      TAB_CONTACT_POLYG(nb_ligneCONTACT_POLYG)%ianent   = TAB_CONTACT(i)%ianent
      TAB_CONTACT_POLYG(nb_ligneCONTACT_POLYG)%n        = TAB_CONTACT(i)%n
      TAB_CONTACT_POLYG(nb_ligneCONTACT_POLYG)%t        = TAB_CONTACT(i)%t
      TAB_CONTACT_POLYG(nb_ligneCONTACT_POLYG)%coor_ctc = TAB_CONTACT(i)%coor_ctc
      TAB_CONTACT_POLYG(nb_ligneCONTACT_POLYG)%rn       = TAB_CONTACT(i)%rn
      TAB_CONTACT_POLYG(nb_ligneCONTACT_POLYG)%rt       = TAB_CONTACT(i)%rt
      TAB_CONTACT_POLYG(nb_ligneCONTACT_POLYG)%rs       = TAB_CONTACT(i)%rs
      TAB_CONTACT_POLYG(nb_ligneCONTACT_POLYG)%nature   = TAB_CONTACT(i)%nature
      TAB_CONTACT_POLYG(nb_ligneCONTACT_POLYG)%type     = 1
      TAB_CONTACT_POLYG(nb_ligneCONTACT_POLYG)%status_points = TAB_CONTACT(i)%status_points
      TAB_CONTACT_POLYG(nb_ligneCONTACT_POLYG)%status   = TAB_CONTACT(i)%status
    end if
    if (TAB_CONTACT(i)%type==2) then
      if (first_double) then
        coor_ctc1_temp = TAB_CONTACT(i)%coor_ctc(1)
        coor_ctc2_temp = TAB_CONTACT(i)%coor_ctc(2)
        coor_ctc3_temp = TAB_CONTACT(i)%coor_ctc(3)
        rn_temp        = TAB_CONTACT(i)%rn
        rt_temp        = TAB_CONTACT(i)%rt
        rs_temp        = TAB_CONTACT(i)%rs
        first_double   = .false.
        l_stat_p       = TAB_CONTACT(i)%status_points
        gap_temp       = TAB_CONTACT(i)%gap
        cycle
      else if (.not. first_double) then
        nb_ligneCONTACT_POLYG                                 = nb_ligneCONTACT_POLYG+1
        TAB_CONTACT_POLYG(nb_ligneCONTACT_POLYG)%icdent       = TAB_CONTACT(i)%icdent
        TAB_CONTACT_POLYG(nb_ligneCONTACT_POLYG)%ianent       = TAB_CONTACT(i)%ianent
        TAB_CONTACT_POLYG(nb_ligneCONTACT_POLYG)%n            = TAB_CONTACT(i)%n
        TAB_CONTACT_POLYG(nb_ligneCONTACT_POLYG)%t            = TAB_CONTACT(i)%t
        TAB_CONTACT_POLYG(nb_ligneCONTACT_POLYG)%coor_ctc(1)  = (TAB_CONTACT(i)%coor_ctc(1) + coor_ctc1_temp) /2
        TAB_CONTACT_POLYG(nb_ligneCONTACT_POLYG)%coor_ctc(2)  = (TAB_CONTACT(i)%coor_ctc(2) + coor_ctc2_temp) /2
        TAB_CONTACT_POLYG(nb_ligneCONTACT_POLYG)%coor_ctc(3)  = (TAB_CONTACT(i)%coor_ctc(3) + coor_ctc3_temp) /2
        TAB_CONTACT_POLYG(nb_ligneCONTACT_POLYG)%gap          = (TAB_CONTACT(i)%gap + gap_temp) /2
        TAB_CONTACT_POLYG(nb_ligneCONTACT_POLYG)%rn           = TAB_CONTACT(i)%rn + rn_temp
        TAB_CONTACT_POLYG(nb_ligneCONTACT_POLYG)%rt           = TAB_CONTACT(i)%rt + rt_temp
        TAB_CONTACT_POLYG(nb_ligneCONTACT_POLYG)%rs           = TAB_CONTACT(i)%rs + rs_temp
        TAB_CONTACT_POLYG(nb_ligneCONTACT_POLYG)%nature       = TAB_CONTACT(i)%nature
        TAB_CONTACT_POLYG(nb_ligneCONTACT_POLYG)%type         = 2

        ! Choosing the most stable of the status points
        TAB_CONTACT_POLYG(nb_ligneCONTACT_POLYG)%status_points= min(l_stat_p, TAB_CONTACT(i)%status_points)
        first_double                                          = .true.
      end if
    end if
  end do

  ! If there are cohesive particles they should know they are in a group
  n_groups = 0

  if (c_draw == 1) then
    print*,'---> Cohesive groups'

    TAB_POLYG(:)%group = 0
    ! For all the particles
    do i=1, n_particles
      ! If this particle does not have a group
      if (TAB_POLYG(i)%group == 0) then
        ! Lets check if its contacts have
        do j=1, nb_ligneCONTACT
          icdent = TAB_CONTACT(j)%icdent
          ianent = TAB_CONTACT(j)%ianent
          ! is it cohesive?
          if (TAB_CONTACT(j)%status(1:1) /= 'W') cycle
          if (TAB_CONTACT(j)%gap > 0.0001) cycle
          ! if found as candidate
          if (icdent==i) then
            ! and the antagonist has a group
            if (TAB_POLYG(ianent)%group > 0) then
              ! Assign the same
              TAB_POLYG(icdent)%group = TAB_POLYG(ianent)%group
            end if
          ! if found as antagonist
          else if (ianent==i) then
            ! and the antagonist has a group
            if (TAB_POLYG(icdent)%group > 0) then
              ! Assign the same
              TAB_POLYG(ianent)%group = TAB_POLYG(icdent)%group
            end if
          end if
          ! We can break the do if found
          if (TAB_POLYG(i)%group > 0) then
            exit
          end if
        end do
      end if

      ! If was assigned. Then loop for non initialized neighbors
      if (TAB_POLYG(i)%group > 0) then
        ! And we loop looking for neighbors
        do j=1, nb_ligneCONTACT
          icdent = TAB_CONTACT(j)%icdent
          ianent = TAB_CONTACT(j)%ianent
          ! is it cohesive?
          if (TAB_CONTACT(j)%status(1:1) /= 'W') cycle
          if (TAB_CONTACT(j)%gap > 0.0001) cycle
          ! if found as candidate
          if (icdent==i) then
            ! we change the antagonist
            TAB_POLYG(ianent)%group = TAB_POLYG(icdent)%group
          ! if found as antagonist
          else if (ianent==i) then
            ! we change the candidate
            TAB_POLYG(icdent)%group = TAB_POLYG(ianent)%group
          end if
        end do
      end if

      ! if after all this research the group is still zero... so a new group begins
      if (TAB_POLYG(i)%group == 0) then
        n_groups = n_groups + 1
        TAB_POLYG(i)%group = n_groups
        ! And we loop again to change its neighbors
        do j=1, nb_ligneCONTACT
          icdent = TAB_CONTACT(j)%icdent
          ianent = TAB_CONTACT(j)%ianent
          ! is it cohesive?
          if (TAB_CONTACT(j)%status(1:1) /= 'W') cycle
          if (TAB_CONTACT(j)%gap > 0.0001) cycle
          ! if found as candidate
          if (icdent==i) then
            ! we change the antagonist
            TAB_POLYG(ianent)%group = TAB_POLYG(icdent)%group
            ! Look for neighbors of the neighbor
            do k=1, nb_ligneCONTACT
              if (TAB_CONTACT(k)%status(1:1) /= 'W') cycle
              if (TAB_CONTACT(k)%gap > 0.0001) cycle
              if (TAB_CONTACT(k)%icdent == ianent) then
                TAB_POLYG(TAB_CONTACT(k)%icdent)%group = TAB_POLYG(ianent)%group
              else if (TAB_CONTACT(j)%ianent == ianent) then
                TAB_POLYG(TAB_CONTACT(k)%ianent)%group = TAB_POLYG(ianent)%group
              end if
            end do
          ! if found as antagonist
          else if (ianent==i) then
            ! we change the candidate
            TAB_POLYG(icdent)%group = TAB_POLYG(ianent)%group
            ! Look for neighbors of the neighbor
            do k=1, nb_ligneCONTACT
              if (TAB_CONTACT(k)%status(1:1) /= 'W') cycle
              if (TAB_CONTACT(k)%gap > 0.0001) cycle
              if (TAB_CONTACT(k)%icdent == icdent) then
                TAB_POLYG(TAB_CONTACT(k)%icdent)%group = TAB_POLYG(icdent)%group
              else if (TAB_CONTACT(j)%ianent == icdent) then
                TAB_POLYG(TAB_CONTACT(k)%ianent)%group = TAB_POLYG(icdent)%group
              end if
            end do
          end if
        end do
      end if
    end do
  end if 

!      n_groups = n_groups + 1
!      TAB_POLYG(i)%group = n_groups
!      ! We start to look for this particle in the contact list 
!      do j=1, nb_ligneCONTACT
!        ! If we find it 
!        if(TAB_CONTACT(j)%icdent == i) then
!          ! we ask it is cohesive  
!          if (TAB_CONTACT(j)%status(1:1) == 'W') then
!            ! If its antagonist does not have already a group
!            if (TAB_POLYG(TAB_CONTACT(j)%ianent)%group == 0) then
!              ! We assign the group to the antagonist
!              TAB_POLYG(TAB_CONTACT(j)%ianent)%group = n_groups
!
!            ! If its antagonist already has a group
!            else
!              ! We modify the current group
!              n_groups = TAB_POLYG(TAB_CONTACT(j)%ianent)%group
!              ! We reassign the group of the current particle
!              TAB_POLYG(i)%group = n_groups
!            end if 
!            ! We look for other contacts including the antagonist 
!            do k=j+1, nb_ligneCONTACT
!              if(TAB_CONTACT(k)%icdent == TAB_CONTACT(j)%ianent .or. TAB_CONTACT(k)%ianent == TAB_CONTACT(j)%ianent) then
!                if (TAB_CONTACT(k)%status(1:1) == 'W') then
!                  TAB_POLYG(TAB_CONTACT(k)%icdent)%group = n_groups
!                  TAB_POLYG(TAB_CONTACT(k)%ianent)%group = n_groups
!                end if
!              end if
!            end do
!          end if
!        end if
!      enddo
!    ! If the particle has already a group
!    else
!      ! We look for other contacts with this particle 
!      do j=1, nb_ligneCONTACT
!        ! If there it exits as candidate in the contact list 
!        if(TAB_CONTACT(j)%icdent == i) then
!          ! If it is cohesive
!          if (TAB_CONTACT(j)%status(1:1) == 'W') then
!            ! We asign the group to the current antagonist
!            TAB_POLYG(TAB_CONTACT(j)%ianent)%group = TAB_POLYG(i)%group
!          end if
!        end if
!      end do
!    end if
!  end do

  print*,'---> Reading Done!'
end subroutine read_Vloc_dof_bodies


!==============================================================================
! Walls positions
!==============================================================================
subroutine wallsposition
  
  implicit none
  
  if (n_walls == 4) then
    if (first_over_all) then
      write(108,*) '#   time     ', '   wall 1-X  ', '   wall 1-Y  ', &
                                    '   wall 2-X  ', '   wall 2-Y  ', &
                                    '   wall 3-X  ', '   wall 3-Y  ', &
                                    '   wall 4-X  ', '   wall 4-Y  '
    endif
    write(108,'(9(1X,E12.5))') time, TAB_PLAN(1)%center(1), TAB_PLAN(1)%center(2)+TAB_PLAN(1)%ax2 &
                               , TAB_PLAN(2)%center(1), TAB_PLAN(2)%center(2)-TAB_PLAN(2)%ax2&
                               , TAB_PLAN(3)%center(1)+TAB_PLAN(3)%ax2, TAB_PLAN(3)%center(2) &
                               , TAB_PLAN(4)%center(1)-TAB_PLAN(4)%ax2, TAB_PLAN(4)%center(2)

  else if (n_walls == 2) then
    if (first_over_all) then
      write(108,*) '#   time     ', '   wall 1-X  ', '   wall 1-Y  ', &
                                    '   wall 2-X  ', '   wall 2-Y  '
    end if
    write(108,'(5(1X,E12.5))') time, TAB_PLAN(1)%center(1), TAB_PLAN(1)%center(2)+TAB_PLAN(1)%ax2 &
                                 , TAB_PLAN(2)%center(1), TAB_PLAN(2)%center(2)-TAB_PLAN(2)%ax2
  else
    print*, 'Invalid number of walls'
    stop
  end if
  print*, 'Write Walls position            ---> Ok!'
end subroutine wallsposition


!==============================================================================
! Walls forces
!==============================================================================
subroutine walls_force
  
  implicit none
  
  integer                            :: i, idw1, idw2, idw3, idw4
  real (kind=8)                      :: fn1, fn2, fn3, fn4
  real (kind=8)                      :: ft1, ft2, ft3, ft4
  
  fn1 = 0
  fn2 = 0
  fn3 = 0
  fn4 = 0 
  
  ft1 = 0
  ft2 = 0
  ft3 = 0
  ft4 = 0 
  
  idw1 = n_particles + 1 
  idw2 = n_particles + 2 
  idw3 = n_particles + 3
  idw4 = n_particles + 4
  
  
  do i=1, nb_ligneCONTACT_POLYG
    if (TAB_CONTACT_POLYG(i)%nature== 'PLJCx') then
      if(TAB_CONTACT_POLYG(i)%ianent == idw1) then
        fn1 = fn1 + abs(TAB_CONTACT_POLYG(i)%rn)
        ft1 = ft1 + abs(TAB_CONTACT_POLYG(i)%rt)
        
      elseif(TAB_CONTACT_POLYG(i)%ianent == idw2) then
        fn2 = fn2 + abs(TAB_CONTACT_POLYG(i)%rn)
        ft2 = ft2 + abs(TAB_CONTACT_POLYG(i)%rt)
      end if 
      if (n_walls .gt. 2) then 
        if(TAB_CONTACT_POLYG(i)%ianent == idw3) then
          fn3 = fn3 + abs(TAB_CONTACT_POLYG(i)%rn)
          ft3 = ft3 + abs(TAB_CONTACT_POLYG(i)%rt)

        elseif(TAB_CONTACT_POLYG(i)%ianent == idw4) then
          fn4 = fn4 + abs(TAB_CONTACT_POLYG(i)%rn)
          ft4 = ft4 + abs(TAB_CONTACT_POLYG(i)%rt)
        end if
      end if
    end if
  end do
  
  if (n_walls == 4) then 
    if (first_over_all) then
      write(109,*) '#   time     ', '     fn1     ', '     ft1     ', &
                   '     fn2     ', '     ft2     ', &
                   '     fn3     ', '     ft3     ', &
                   '     fn4     ', '     ft4    '
    endif
    write(109,'(9(1X,E12.5))') time, fn1, ft1, fn2, ft2, fn3, ft3, fn4, ft4
  else if (n_walls == 2) then
    if (first_over_all) then
      write(109,*) '#   time     ', '     fn1     ', '     ft1     ', &
                   '     fn2     ', '     ft2     '
    endif
    write(109,'(5(1X,E12.5))') time, fn1, ft1, fn2, ft2
  end if

  print*, 'Write Walls forces  ---> Ok!'

end subroutine walls_force


!==============================================================================
! Calculation of compacity
!==============================================================================
subroutine compacity
  
  implicit none
  
  real(kind=8)                             :: Area_grain,Area_total
  integer                                  :: i
  
  Area_grain = 0.D0
  
  do i=1,n_particles
    Area_grain = Area_grain + TAB_POLYG(i)%Area
  end do
  
  Area_total = width * height
  
  if(first_over_all) then
    write(102,*) '#   time     ', '   height    ', '   width    ', '    V_s/V    '
  endif
  
  write(102,'(6(1X,E12.5))') time, height, width, Area_grain/Area_total
  
  print*, 'Write Compacity                 ---> Ok!'
  
end subroutine compacity


!==============================================================================
! Calculation of q/p
!==============================================================================
subroutine qoverp
  
  implicit none
  
  integer                                  :: i,cd,an
  integer                                  :: ierror,matz,lda
  real(kind=8)                             :: Rtik,Rnik
  real(kind=8)                             :: S1,S2, ang_c
  real(kind=8),dimension(2)                :: nik,tik,Lik,Fik
  real(kind=8),dimension(2)                :: wr,wi
  real(kind=8),dimension(2,2)              :: Moment
  real(kind=8),dimension(2,2)              :: localframe
  
  Moment(:,:) = 0.0
  
  do i=1,nb_ligneCONTACT_POLYG
    if  (TAB_CONTACT_POLYG(i)%nature == 'PLJCx') cycle
    cd      = TAB_CONTACT_POLYG(i)%icdent
    an      = TAB_CONTACT_POLYG(i)%ianent
    nik(1)  = TAB_CONTACT_POLYG(i)%n(1)
    nik(2)  = TAB_CONTACT_POLYG(i)%n(2)
    tik(1)  = TAB_CONTACT_POLYG(i)%t(1)
    tik(2)  = TAB_CONTACT_POLYG(i)%t(2)
    Rtik    = TAB_CONTACT_POLYG(i)%rt
    Rnik    = TAB_CONTACT_POLYG(i)%rn
    
    ! Only active contacts
    if (abs(Rnik) .le. 1.D-9 .and. abs(Rnik) .le. 1.D-9) cycle
    
    Lik(1) = TAB_POLYG(cd)%center(1)-TAB_POLYG(an)%center(1)
    Lik(2) = TAB_POLYG(cd)%center(2)-TAB_POLYG(an)%center(2)
    Fik(1) = (Rnik*nik(1)+Rtik*tik(1))
    Fik(2) = (Rnik*nik(2)+Rtik*tik(2))
    
    Moment(1,1:2) = Fik(1)*Lik(1:2) + Moment(1,1:2)
    Moment(2,1:2) = Fik(2)*Lik(1:2) + Moment(2,1:2)
  end do
  
  Moment = Moment / (height*width)
  
  if (first_over_all) then
    write(103,*) '#   time     ', '     s11     ', '     s12     ', '     s22     ', &
                                  '      S1     ', '      S2     ', '(S1-S2)/(S1+S2)', '  theta_s  '
  endif

  write(103,'(4(1X,E12.5))', advance='no') time, Moment(1,1), Moment(1,2), Moment(2,2)

  lda  = 2
  matz = 1
  
  ! Finding the eigenvalues of the global stress tensor
  call rg(lda, 2, Moment, wr, wi, matz, localframe, ierror)
  S1 = max( wr(1),wr(2) )
  S2 = min( wr(1),wr(2) )

  ! Computing the orientation of the stress tensor 
  ang_c = abs(0.5*atan2(-2*Moment(1,2),(Moment(1,1)-Moment(2,2))))

  ang_c = ang_c*180/pi

  write(103,'(4(1X,E12.5))') S1,S2,(S1-S2)/(S1+S2), ang_c
  
  print*, 'Write QoverP                    ---> Ok!'
  
end subroutine qoverp


!==============================================================================
! Coordination number
!==============================================================================
subroutine nb_coordination
  
  implicit none
  
  real(kind=8)                             :: z,zc,zp
  integer                                  :: i
  
  ! Initializing variables
  z =0.
  zc=0.
  zp=0.
  
  ! Computing the number of particles in the box
  zp = n_particles
  
  ! Computing the number of contacts
  do i=1, nb_ligneCONTACT_POLYG
    ! Only between particles
    if (TAB_CONTACT_POLYG(i)%nature /= 'PLPLx') cycle
    ! Contacts in the box
    ! Computing the number of contacts if they are active
    if (abs(TAB_CONTACT_POLYG(i)%rn) .lt. 1e-9 .and. abs(TAB_CONTACT_POLYG(i)%rt) .lt. 1e-9) then
      cycle 
    else 
      zc=zc+1
    end if
  end do
  
  z = 2 * zc / zp
  
  if (first_over_all) then
    write(100,*) '#   time     ','      z      '
  endif
  
  write(100,'(2(1X,E12.5))') time, z
  
  print*, 'Write Coordination              ---> Ok!'
  
end subroutine nb_coordination


!==============================================================================
! Grain size distribution
!==============================================================================
  subroutine granulometry
    
    implicit none
    
    real(kind=8)                             :: rmin, rmax, interval, porcentNum
    real(kind=8)                             :: areaTotal, porcentArea, area
    integer                                  :: i, j, until, sumaInt, vipar
    integer, dimension(20000)                :: cant
    real(kind=8), dimension(20000)           :: tams
    character(len=23)                        :: nom
    
    nom       =  './POSTPRO/GDISTRIB.0000'
    
    if (compteur_clout<10) then
      WRITE(nom(22:23),'(I1)')   compteur_clout
    else if ( (compteur_clout>=10) .and. (compteur_clout<100) ) then
      WRITE(nom(21:23),'(I2)')   compteur_clout
    else if ( (compteur_clout>=100).and. (compteur_clout<1000) ) then
      WRITE(nom(20:23),'(I3)')   compteur_clout
    end if
    
    open(unit=107,file=nom,status='replace')
    
    rmin = 9999.
    rmax = -1.
    vipar = 0
    
    ! Finding the maximun and minimun radius
    do i=1, n_particles
      vipar = vipar + 1
      if (TAB_POLYG(i)%radius<rmin) then
        rmin = TAB_POLYG(i)%radius
      end if 
      if (TAB_POLYG(i)%radius>rmax) then
        rmax = TAB_POLYG(i)%radius
      end if
    end do
    
    until = 500  !here is specified the number of intervals for the constrution of the distribution
    
    interval = (rmax - rmin)/until            !size ranges
    
    do j=1, until
      do i=1, n_particles
        if(TAB_POLYG(i)%radius .ge. (rmin + (j-1)*interval) .and. TAB_POLYG(i)%radius .lt. (interval*j + rmin)) then
          cant(j) = cant(j) + 1
        end if
      end do
    end do
    
    do i=1, until
      tams(i) = rmin + interval*i - (interval/2)
    end do
    
    
    sumaInt = 0
    porcentNum = 0
    porcentArea = 0
    area = 0
    areaTotal = 0 
    
    do i=1, until
      areaTotal = areaTotal + (tams(i)**2)*pi*cant(i)
    enddo
    
    do j = 1, until  
      area = area + (tams(j)**2)*pi*cant(j)
      porcentArea = (area/areaTotal)*100
      
      write(107,'(4(1X,E12.5))') tams(j), porcentArea
    enddo
    
    close(107)
    
    print*, 'Write GDI                       ---> Ok!'
    
  end subroutine granulometry


!================================================
! Drawing in vtk
!================================================
subroutine draw
  
  implicit none
  
  character(:), allocatable                      ::  vtk_part, vtk_inter, vtk_counter
  character(:), allocatable                      ::  vtk_n_points
  logical                                        ::  dir_vtk
  character(len=6)                               ::  vtk_c_temp,  v_n_vertices
  integer                                        ::  i, j, k
  integer                                        ::  n_l_vertices, curr_l_faces
  integer                                        ::  n_l_fields, l_counter
  real(kind=8), dimension(3)                     ::  curr_l_vector
  integer                                        ::  l_cdt, l_ant
  
  ! Variables for the forces
  character(:), allocatable                      ::  vtk_forces
  character(:), allocatable                      ::  vtk_n_v_forces
  character(len=6)                               ::  v_n_v_forces
  real(kind=8)                                   ::  vtk_ave_force, l_force_scale, l_force
  real(kind=8)                                   ::  l_rn, l_rt, l_rs
  real(kind=8), dimension(3)                     ::  vtk_cd_center, vtk_an_center, v_l_normal, v_l_t
  integer                                        ::  f_counter
  
  ! Variables for the contact network
  character(:), allocatable                      ::  vtk_ctc_net
  character(:), allocatable                      ::  vtk_n_v_ctcnet
  integer                                        ::  ctc_counter
  character(len=6)                               ::  v_n_v_ctcnet
  real(kind=8)                                   ::  vtk_ave_rad, l_ctc_scale
  
  
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
  write(vtk_c_temp, '(I6)') compteur_clout
  
  
  if (compteur_clout<10) then
    vtk_counter = vtk_c_temp(6:6)
  else if (compteur_clout >= 10 .and. compteur_clout < 100) then
    vtk_counter = vtk_c_temp(5:6)
  else if (compteur_clout >= 100 .and. compteur_clout<1000) then
    vtk_counter = vtk_c_temp(4:6)
  else
    print*, 'Cannot draw more than 1000 files'
    stop
  end if
  
  vtk_part = './POSTPRO/VTK/rigid_' // vtk_counter // '.vtk'
  
  open(unit=121, file=vtk_part, status='replace')
  write(121,'(A)') '# vtk DataFile Version 3.0'
  write(121,'(A,I6)') 'RIGID ', compteur_clout
  write(121,'(A)') 'ASCII'
  write(121,'(A)') 'DATASET POLYDATA'
  
  ! Writing the number of vertices
  n_l_vertices = 0
  do i=1, n_particles
    n_l_vertices = n_l_vertices + TAB_POLYG(i)%nb_vertex 
  end do
  
  ! Adding the vertices of the walls
  n_l_vertices = n_l_vertices + n_walls*4
  
  write(v_n_vertices, '(I6)') n_l_vertices
  
  vtk_n_points = 'POINTS ' // v_n_vertices // ' float'
  write(121,'(A)') vtk_n_points
  
  ! Writing the coordinates of the vertices of particles and walls
  ! 3 vertices per row
  ! 
  k=0
  do i=1, n_particles + n_walls
    !write(121,*) 'Particle', i
    if (i .le. n_particles) then
      do j=1, TAB_POLYG(i)%nb_vertex
        k = k + 1
        if(k .lt. 3) then 
          if (TAB_POLYG(i)%vertex(1,j) .lt. 0) then
            write(121,'(F11.8,A)', advance = 'no') TAB_POLYG(i)%vertex(1,j), ' '
          else
            write(121,'(F9.7,A)', advance = 'no') TAB_POLYG(i)%vertex(1,j), ' '
          end if
          
          if (TAB_POLYG(i)%vertex(2,j) .lt. 0) then
            write(121,'(F11.8,A)', advance = 'no') TAB_POLYG(i)%vertex(2,j), ' '
          else 
            write(121,'(F9.7,A)', advance = 'no') TAB_POLYG(i)%vertex(2,j), ' '
          end if
          
          if (TAB_POLYG(i)%vertex(3,j) .lt. 0) then
            write(121,'(F11.8,A)', advance = 'no') TAB_POLYG(i)%vertex(3,j), ' '
          else 
            write(121,'(F9.7,A)', advance = 'no') TAB_POLYG(i)%vertex(3,j), ' '
          end if
        else
          if (TAB_POLYG(i)%vertex(1,j) .lt. 0) then
            write(121,'(F11.8,A)', advance = 'no') TAB_POLYG(i)%vertex(1,j), ' '
          else 
            write(121,'(F9.7,A)', advance = 'no') TAB_POLYG(i)%vertex(1,j), ' '
          end if
          
          if (TAB_POLYG(i)%vertex(2,j) .lt. 0) then
            write(121,'(F11.8,A)', advance = 'no') TAB_POLYG(i)%vertex(2,j), ' '
          else
            write(121,'(F9.7,A)', advance = 'no') TAB_POLYG(i)%vertex(2,j), ' '
          end if
          
          if (TAB_POLYG(i)%vertex(3,j) .lt. 0) then
            write(121,'(F11.8,A)') TAB_POLYG(i)%vertex(3,j), ' '
          else 
            write(121,'(F9.7,A)') TAB_POLYG(i)%vertex(3,j), ' '
          end if
          k = 0
        end if
      end do
    else
      ! Writing the vertices of the walls
      do j=1, 4 ! each wall has 4 vertices
        k = k +1
        if(k .lt. 3) then 
          if (TAB_PLAN(i-n_particles)%vertex(1,j) .lt. 0) then
            write(121,'(F11.8,A)', advance = 'no') TAB_PLAN(i-n_particles)%vertex(1,j), ' '
          else
            write(121,'(F9.7,A)', advance = 'no') TAB_PLAN(i-n_particles)%vertex(1,j), ' '
          end if
          
          if (TAB_PLAN(i-n_particles)%vertex(2,j) .lt. 0) then
            write(121,'(F11.8,A)', advance = 'no') TAB_PLAN(i-n_particles)%vertex(2,j), ' '
          else 
            write(121,'(F9.7,A)', advance = 'no') TAB_PLAN(i-n_particles)%vertex(2,j), ' '
          end if
          
          if (TAB_PLAN(i-n_particles)%vertex(3,j) .lt. 0) then
            write(121,'(F11.8,A)', advance = 'no') TAB_PLAN(i-n_particles)%vertex(3,j), ' '
          else 
            write(121,'(F9.7,A)', advance = 'no') TAB_PLAN(i-n_particles)%vertex(3,j), ' '
          end if
        else
          if (TAB_PLAN(i-n_particles)%vertex(1,j) .lt. 0) then
            write(121,'(F11.8,A)', advance = 'no') TAB_PLAN(i-n_particles)%vertex(1,j), ' '
          else 
            write(121,'(F9.7,A)', advance = 'no') TAB_PLAN(i-n_particles)%vertex(1,j), ' '
          end if
          
          if (TAB_PLAN(i-n_particles)%vertex(2,j) .lt. 0) then
            write(121,'(F11.8,A)', advance = 'no') TAB_PLAN(i-n_particles)%vertex(2,j), ' '
          else 
            write(121,'(F9.7,A)', advance = 'no') TAB_PLAN(i-n_particles)%vertex(2,j), ' '
          end if
          
          if (TAB_PLAN(i-n_particles)%vertex(3,j) .lt. 0) then
            write(121,'(F11.8,A)') TAB_PLAN(i-n_particles)%vertex(3,j), ' '
          else 
            write(121,'(F9.7,A)') TAB_PLAN(i-n_particles)%vertex(3,j), ' '
          end if
          
          k = 0
          
        end if
      end do
    end if
  end do
  
  write(121, '(A)') ' '
  
  ! Writing the conectivity between vertices
  ! Plus 4 walls
  write(121,'(A)', advance='no') 'POLYGONS '
  write(121, '(2(I6,A))') n_particles + n_walls, ' ' , (n_l_vertices + n_particles)+(n_walls)
  
  curr_l_faces = 0
  
  do i=1, n_particles
    !write(121,*) 'Particle', i
    ! Write the number of vertices
    write(121, '(I2,A)', advance= 'no') TAB_POLYG(i)%nb_vertex, ' '
    
    ! Write its consecutive conectivity
    do j=1, TAB_POLYG(i)%nb_vertex
      
      ! ... writing precisely 
      if (j .lt. TAB_POLYG(i)%nb_vertex) then
        if (curr_l_faces - 1 + j .lt. 10) then
          write(121, '(I1)', advance = 'no') curr_l_faces - 1 + j
        else if (curr_l_faces - 1 + j .lt. 100) then
          write(121, '(I2)', advance = 'no') curr_l_faces - 1 + j
        else if (curr_l_faces - 1 + j .lt. 1000) then
          write(121, '(I3)', advance = 'no') curr_l_faces - 1 + j
        else if (curr_l_faces - 1 + j .lt. 10000) then
          write(121, '(I4)', advance = 'no') curr_l_faces - 1 + j
        else if (curr_l_faces - 1 + j .lt. 100000) then
          write(121, '(I5)', advance = 'no') curr_l_faces - 1 + j
        end if
        
        write(121, '(A)', advance='no') ' '
        
      else
        if (curr_l_faces - 1 + j .lt. 10) then
          write(121, '(I1)') curr_l_faces - 1 + j
        else if (curr_l_faces - 1 + j .lt. 100) then
          write(121, '(I2)') curr_l_faces - 1 + j
        else if (curr_l_faces - 1 + j .lt. 1000) then
          write(121, '(I3)') curr_l_faces - 1 + j
        else if (curr_l_faces - 1 + j .lt. 10000) then
          write(121, '(I4)') curr_l_faces - 1 + j
        else if (curr_l_faces - 1 + j .lt. 100000) then
          write(121, '(I5)') curr_l_faces - 1 + j
        end if
      end if
      
    end do
    curr_l_faces = curr_l_faces + TAB_POLYG(i)%nb_vertex
  end do
  
  ! Writing the conectivity for the walls
  ! For all the walls 
  do i=1, n_walls
    
    !!!!! 1st (1-2-3-4)
    write(121, '(I1,A)', advance= 'no') 4, ' '
    
    if (curr_l_faces .lt. 10) then
      write(121, '(I1)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 100) then
      write(121, '(I2)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 1000) then
      write(121, '(I3)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 10000) then
      write(121, '(I4)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 100000) then
      write(121, '(I5)', advance = 'no') curr_l_faces
    end if
    
    write(121, '(A)', advance='no') ' '
    
    if (curr_l_faces +1 .lt. 10) then
      write(121, '(I1)', advance = 'no') curr_l_faces +1 
    else if (curr_l_faces +1 .lt. 100) then
      write(121, '(I2)', advance = 'no') curr_l_faces +1 
    else if (curr_l_faces +1 .lt. 1000) then
      write(121, '(I3)', advance = 'no') curr_l_faces +1 
    else if (curr_l_faces +1 .lt. 10000) then
      write(121, '(I4)', advance = 'no') curr_l_faces +1
    else if (curr_l_faces +1 .lt. 100000) then
      write(121, '(I5)', advance = 'no') curr_l_faces +1 
    end if
    
    write(121, '(A)', advance='no') ' '
    
    if (curr_l_faces +2 .lt. 10) then
      write(121, '(I1)', advance = 'no') curr_l_faces +2 
    else if (curr_l_faces +2 .lt. 100) then
      write(121, '(I2)', advance = 'no') curr_l_faces +2 
    else if (curr_l_faces +2 .lt. 1000) then
      write(121, '(I3)', advance = 'no') curr_l_faces +2 
    else if (curr_l_faces +2 .lt. 10000) then
      write(121, '(I4)', advance = 'no') curr_l_faces +2
    else if (curr_l_faces +2 .lt. 100000) then
      write(121, '(I5)', advance = 'no') curr_l_faces +2 
    end if
    
    write(121, '(A)', advance='no') ' '
    
    if (curr_l_faces +3 .lt. 10) then
      write(121, '(I1)') curr_l_faces +3
    else if (curr_l_faces +3 .lt. 100) then
      write(121, '(I2)') curr_l_faces +3 
    else if (curr_l_faces +3 .lt. 1000) then
      write(121, '(I3)') curr_l_faces +3 
    else if (curr_l_faces +3 .lt. 10000) then
      write(121, '(I4)') curr_l_faces +3
    else if (curr_l_faces +3 .lt. 100000) then
      write(121, '(I5)') curr_l_faces +3
    end if
    
    curr_l_faces = curr_l_faces + 4
    
  end do
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!! FIELDS !!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! How many fields there will be?
  n_l_fields = 7                          ! They are: Id, Material, Disp, Veloc, Spin, Z (coordination), group
                                          ! ...
  
  ! A blank space
  write(121,'(A)') ''
  ! The cells begin
  write(121,'(A)', advance='no') 'CELL_DATA'
  
  ! Writing the number of data by field. It corresponds to the same number of particles
  write(121, '(2(I6,A))') n_particles + n_walls
  write(121, '(A,I4)')  'FIELD FieldData', n_l_fields
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Id
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! Naming the field, the dimension of the data, and number of lines (as before FIELD), and data type
  write(121,'(A,I1,I6,A)') 'Id ', 1, (n_particles + n_walls), ' float'
  k = 0
  l_counter = 0
  do i=1, n_particles + n_walls
    l_counter = l_counter + 1
    
    write(121, '(I6)') l_counter
    
  end do
  ! And jump
  !write(121, '(A)') ' '
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Material
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(121,'(A,I1,I6,A)') 'Material ', 1, (n_particles + n_walls), ' float'
  k = 0
  l_counter = 1
  
  do i=1, n_particles + n_walls
    if (i .le. n_particles) then
      ! Material number 1
      write(121, '(I6)') l_counter
    else
      ! Material number 2
      if (i==n_particles+1) then
        l_counter = l_counter +1
      end if
      
      write(121, '(I6)') l_counter
    end if
  end do
  
  ! And jump
  !write(121, '(A)') ' '
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Displacement
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  curr_l_vector(:) = 0.D0
  
  write(121,'(A,I1,I6,A)') 'Disp ', 3, (n_particles + n_walls), ' float'
  k = 0
  do i=1, n_particles + n_walls
    ! For the particles
    if (i .le. n_particles) then
      curr_l_vector(:) = TAB_POLYG(i)%center(:) - TAB_POLYG(i)%center_ref
      write(121, '(3(F12.9,A))') curr_l_vector(1), ' ', curr_l_vector(2), ' ', curr_l_vector(3)
    else
    !For the walls 
      curr_l_vector(:) = TAB_PLAN(i-n_particles)%center(:) - TAB_PLAN(i-n_particles)%center_ref
      write(121, '(3(F12.9,A))') curr_l_vector(1), ' ', curr_l_vector(2), ' ', curr_l_vector(3)
    end if
  end do
  
  ! And jump
  !write(121, '(A)') ' '
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Velocity
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  curr_l_vector(:) = 0.D0
  
  write(121,'(A,I1,I6,A)') 'Veloc ', 3, (n_particles + n_walls), ' float'
  k = 0
  do i=1, n_particles + n_walls
    if (i .le. n_particles) then
      curr_l_vector(1) = TAB_POLYG(i)%V(1)
      curr_l_vector(2) = TAB_POLYG(i)%V(2)
      curr_l_vector(3) = TAB_POLYG(i)%V(3)
      
      write(121, '(3(F13.9,A))') curr_l_vector(1), ' ', curr_l_vector(2), ' ', curr_l_vector(3)
    else
      curr_l_vector(1) = TAB_PLAN(i-n_particles)%V(1)
      curr_l_vector(2) = TAB_PLAN(i-n_particles)%V(2)
      curr_l_vector(3) = TAB_PLAN(i-n_particles)%V(3)
      
      write(121, '(3(F13.9,A))') curr_l_vector(1), ' ', curr_l_vector(2), ' ', curr_l_vector(3)
    end if
  end do
  
  ! And jump
  !write(121, '(A)') ' '
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Coordination
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  l_counter = 0
  
  write(121,'(A,I1,I6,A)') 'Z ', 1, (n_particles + n_walls), ' float'
  k = 0
  do i=1, n_particles + n_walls
    l_counter = 0
    if (i .le. n_particles) then
      ! Counting the number of contacts
      do j=1, nb_ligneCONTACT_POLYG
        l_cdt = TAB_CONTACT_POLYG(j)%icdent
        l_ant = TAB_CONTACT_POLYG(j)%ianent
        
        if(TAB_POLYG(l_cdt)%behav /= 'PLEXx' .or. TAB_POLYG(l_ant)%behav /= 'PLEXx') cycle
        if(l_cdt == i .or. l_ant == i) then
          l_counter = l_counter + 1
        end if
      end do
      
      write(121, '(I4,A)') l_counter
      
    else
      write(121, '(I4,A)') 0
    end if
  end do
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Float
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  l_counter = 0
  
  write(121,'(A,I1,I6,A)') 'Float_stat ', 1, (n_particles + n_walls), ' float'
  k = 0
  do i=1, n_particles + n_walls
    l_counter = 0
    if (i .le. n_particles) then
      ! Counting the number of contacts
      do j=1, nb_ligneCONTACT_POLYG
        
        if(TAB_CONTACT_POLYG(j)%nature /= 'PLPLx') cycle
        
        l_cdt = TAB_CONTACT_POLYG(j)%icdent
        l_ant = TAB_CONTACT_POLYG(j)%ianent
        
        if(TAB_POLYG(l_cdt)%behav /= 'PLEXx' .or. TAB_POLYG(l_ant)%behav /= 'PLEXx') cycle
        if(l_cdt == i .or. l_ant == i) then
          ! Only active contacts
          if (abs(TAB_CONTACT_POLYG(j)%rn) .gt. 0.0) then
            l_counter = l_counter + 1
          end if
        end if
      end do
      
      if (l_counter == 0) then
        write(121, '(I4,A)') 1
      else 
        write(121, '(I4,A)') 0
      end if
    else
      write(121, '(I4,A)') -1
    end if
  end do
  
  ! And jump
  !write(121, '(A)') ' '

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Group
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  write(121,'(A,I1,I6,A)') 'Group ', 1, (n_particles + n_walls), ' float'
  k = 0
  l_counter = 1
  
  do i=1, n_particles + n_walls
    if (i .le. n_particles) then
      ! Material number 1
      write(121, '(I6)') TAB_POLYG(i)%group
    else      
      write(121, '(I6)') -1
    end if
  end do
  
  ! And jump
  !write(121, '(A)') ' '
  
  ! Closing the files of the particles
  close(121)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!  FORCES   !!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  vtk_forces = './POSTPRO/VTK/forces_' // vtk_counter // '.vtk'
  
  open(unit=121, file=vtk_forces, status='replace')
  write(121,'(A)') '# vtk DataFile Version 3.0'
  write(121,'(A,I6)') 'FORCES ', compteur_clout
  write(121,'(A)') 'ASCII'
  write(121,'(A)') 'DATASET POLYDATA'
  
  ! Writing the number of vertices for the forces. These are parallelepipeds joining the center of 
  ! particles in contact
  ! Four vertices for each parallelepiped
  f_counter = 0
  do i=1, nb_ligneCONTACT_POLYG
    if(TAB_CONTACT_POLYG(i)%nature /= 'PLPLx') cycle
    f_counter = f_counter + 1
  end do
  
  write(v_n_v_forces, '(I6)') f_counter * 4
  
  vtk_n_v_forces = 'POINTS' // v_n_v_forces // ' float'
  write(121,'(A)') vtk_n_v_forces
  
  ! Finding the average force
  vtk_ave_force = 0.D0
  
  do i = 1, nb_ligneCONTACT_POLYG
    if (TAB_CONTACT_POLYG(i)%nature /= 'PLPLx') cycle
    vtk_ave_force = vtk_ave_force + abs(TAB_CONTACT_POLYG(i)%rn)
  end do
  
  vtk_ave_force = vtk_ave_force/nb_ligneCONTACT_POLYG
  
  ! Force scale parameter
  l_force_scale = 0.D0
  do i=1, n_particles
    l_force_scale = l_force_scale + TAB_POLYG(i)%radius
  end do

  ! The scale is here!
  ! 5% of average radius
  l_force_scale = (l_force_scale/n_particles)*0.05
  
  k=0
  do i=1, nb_ligneCONTACT_POLYG
    
    if (TAB_CONTACT_POLYG(i)%nature /= 'PLPLx') cycle
    ! The particles ids
    l_cdt = TAB_CONTACT_POLYG(i)%icdent
    l_ant = TAB_CONTACT_POLYG(i)%ianent
    ! The center of each particle
    vtk_cd_center(:) = TAB_POLYG(l_cdt)%center(:)
    vtk_an_center(:) = TAB_POLYG(l_ant)%center(:)
    ! The normal and tangential vectors
    v_l_normal(:) = TAB_CONTACT_POLYG(i)%n(:)
    v_l_t(:) = TAB_CONTACT_POLYG(i)%t(:)
    
    ! The forces
    l_rn = TAB_CONTACT_POLYG(i)%rn
    l_rt = TAB_CONTACT_POLYG(i)%rt
    
    l_force = (l_rn/vtk_ave_force)*l_force_scale
    
    ! Write! ... 4 vertices for each line
    
    ! Candidate --- First vertex. 
    write(121, '(3(F15.7,A))', advance='no') vtk_cd_center(1)+l_force*v_l_t(1), ' ', &
                                             vtk_cd_center(2)+l_force*v_l_t(2), ' ', &
                                             vtk_cd_center(3), ' '
    ! Candidate --- Second vertex. 
    write(121, '(3(F15.7,A))', advance='no') vtk_cd_center(1)-l_force*v_l_t(1), ' ', &
                                             vtk_cd_center(2)-l_force*v_l_t(2), ' ', &
                                             vtk_cd_center(3), ' '
    ! Candidate --- Third vertex. 
    write(121, '(3(F15.7,A))', advance='no') vtk_an_center(1)-l_force*v_l_t(1), ' ', &
                                             vtk_an_center(2)-l_force*v_l_t(2), ' ', &
                                             vtk_an_center(3), ' '
    ! Candidate --- 4th vertex. 
    write(121, '(3(F15.7,A))') vtk_an_center(1)+l_force*v_l_t(1), ' ', &
                               vtk_an_center(2)+l_force*v_l_t(2), ' ', &
                               vtk_an_center(3), ' '
  end do
  
  ! Writing the conectivity
  write(121,'(A)', advance='no') 'POLYGONS '
  ! 6 faces for each parallepiped
  write(121, '(2(I6,A))') f_counter, ' ' , (f_counter*5)
  
  curr_l_faces = 0
  
  do i=1, nb_ligneCONTACT_POLYG
    ! Conecting the 4 vertices for each rectangle
    if (TAB_CONTACT_POLYG(i)%nature /= 'PLPLx') cycle
    
    !!!!! 1st (1-2-3-4)
    write(121, '(I1,A)', advance= 'no') 4, ' '
    
    if (curr_l_faces .lt. 10) then
      write(121, '(I1)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 100) then
      write(121, '(I2)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 1000) then
      write(121, '(I3)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 10000) then
      write(121, '(I4)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 100000) then
      write(121, '(I5)', advance = 'no') curr_l_faces
    end if
    
    write(121, '(A)', advance='no') ' '
    
    if (curr_l_faces +1 .lt. 10) then
      write(121, '(I1)', advance = 'no') curr_l_faces +1 
    else if (curr_l_faces +1 .lt. 100) then
      write(121, '(I2)', advance = 'no') curr_l_faces +1 
    else if (curr_l_faces +1 .lt. 1000) then
      write(121, '(I3)', advance = 'no') curr_l_faces +1 
    else if (curr_l_faces +1 .lt. 10000) then
      write(121, '(I4)', advance = 'no') curr_l_faces +1
    else if (curr_l_faces +1 .lt. 100000) then
      write(121, '(I5)', advance = 'no') curr_l_faces +1 
    end if
    
    write(121, '(A)', advance='no') ' '
    
    if (curr_l_faces +2 .lt. 10) then
      write(121, '(I1)', advance = 'no') curr_l_faces +2 
    else if (curr_l_faces +2 .lt. 100) then
      write(121, '(I2)', advance = 'no') curr_l_faces +2 
    else if (curr_l_faces +2 .lt. 1000) then
      write(121, '(I3)', advance = 'no') curr_l_faces +2 
    else if (curr_l_faces +2 .lt. 10000) then
      write(121, '(I4)', advance = 'no') curr_l_faces +2
    else if (curr_l_faces +2 .lt. 100000) then
      write(121, '(I5)', advance = 'no') curr_l_faces +2 
    end if
    
    write(121, '(A)', advance='no') ' '
    
    if (curr_l_faces +3 .lt. 10) then
      write(121, '(I1)') curr_l_faces +3
    else if (curr_l_faces +3 .lt. 100) then
      write(121, '(I2)') curr_l_faces +3 
    else if (curr_l_faces +3 .lt. 1000) then
      write(121, '(I3)') curr_l_faces +3 
    else if (curr_l_faces +3 .lt. 10000) then
      write(121, '(I4)') curr_l_faces +3
    else if (curr_l_faces +3 .lt. 100000) then
      write(121, '(I5)') curr_l_faces +3
    end if
    
    curr_l_faces = curr_l_faces + 4
    
  end do
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!! FORCES FIELDS !!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! How many fields there will be?
  n_l_fields = 4                          ! They are: Type force (Compression - tension), Rn, Rt, status
                                          ! ...
  
  ! A blank space
  write(121,'(A)') ''
  ! The cells begin
  write(121,'(A)', advance='no') 'CELL_DATA'
  
  ! Writing the number of data by field. It corresponds to the same number of forces
  write(121, '(2(I6,A))') f_counter
  write(121, '(A,I4)')  'FIELD FieldData', n_l_fields
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Type force (1 Compression, -1 Tension)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! Naming the field, the dimension of the data, and number of lines (as before FIELD), and data type
  write(121,'(A,I1,I6,A)') 'Type ', 1, f_counter, ' float'
  
  do i=1, nb_ligneCONTACT_POLYG
    if (TAB_CONTACT_POLYG(i)%nature /= 'PLPLx') cycle
    
    if (TAB_CONTACT_POLYG(i)%rn .lt. 0) then
      write(121, '(I3)') -1
    else
      write(121, '(I3)') 1
    end if
  end do
  
  ! And jump
  write(121, '(A)') ' '
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Normal Force
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! Naming the field, the dimension of the data, and number of lines (as before FIELD), and data type
  write(121,'(A,I1,I6,A)') 'RN ', 1, f_counter, ' float'
  do i=1, nb_ligneCONTACT_POLYG
    if (TAB_CONTACT_POLYG(i)%nature /= 'PLPLx') cycle
    write(121, '(F15.7)') TAB_CONTACT_POLYG(i)%rn
  end do
  
  ! And jump
  write(121, '(A)') ' '
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Tangential force (Rs)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! Naming the field, the dimension of the data, and number of lines (as before FIELD), and data type
  write(121,'(A,I1,I6,A)') 'RT ', 1, f_counter, ' float'
  do i=1, nb_ligneCONTACT_POLYG
    if (TAB_CONTACT_POLYG(i)%nature /= 'PLPLx') cycle
    write(121, '(F15.7)') abs(TAB_CONTACT_POLYG(i)%rt)
  end do
  
  ! And jump
  write(121, '(A)') ' '

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Status (3: noctc, 2 Sli, 1 Stick, -1 Wnctc, -2Wsli, -3 Wstck)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! Naming the field, the dimension of the data, and number of lines (as before FIELD), and data type
  write(121,'(A,I1,I6,A)') 'Status ', 1, f_counter, ' float'
  do i=1, nb_ligneCONTACT_POLYG
    if (TAB_CONTACT_POLYG(i)%nature /= 'PLPLx') cycle
    write(121, '(I2)') TAB_CONTACT_POLYG(i)%status_points
  end do
  
  ! And jump
  write(121, '(A)') ' '
  
  close(121)
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!  ONLY CONTACT NETWORK   !!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  vtk_ctc_net = './POSTPRO/VTK/ctc_nt_' // vtk_counter // '.vtk'
  
  open(unit=121, file=vtk_ctc_net, status='replace')
  write(121,'(A)') '# vtk DataFile Version 3.0'
  write(121,'(A,I6)') 'CONTACT NETWORK ', compteur_clout
  write(121,'(A)') 'ASCII'
  write(121,'(A)') 'DATASET POLYDATA'
  
  ! Writing the number of vertices for the contact network. These are parallelepipeds joining the center of 
  ! particles in contact
  ! Four vertices for each parallelepiped
  ctc_counter = 0
  do i=1, nb_ligneCONTACT_POLYG
    if(TAB_CONTACT_POLYG(i)%nature /= 'PLPLx' .or. TAB_CONTACT_POLYG(i)%gap > 0.1E-2) cycle
    if (TAB_CONTACT_POLYG(i)%type == 1) cycle
    ctc_counter = ctc_counter + 1
  end do
  
  write(v_n_v_ctcnet, '(I6)') ctc_counter * 4
  
  vtk_n_v_ctcnet = 'POINTS' // v_n_v_ctcnet // ' float'
  write(121,'(A)') vtk_n_v_ctcnet
  
  ! Finding the average radius
  vtk_ave_rad = 0.D0
  
  do i = 1, n_particles
    vtk_ave_rad = vtk_ave_rad + TAB_POLYG(i)%radius
  end do
  
  vtk_ave_rad = vtk_ave_rad/n_particles
  
  ! The scale is here!
  ! 5% of average radius
  l_ctc_scale = (vtk_ave_rad)*0.05
  
  k=0
  do i=1, nb_ligneCONTACT_POLYG
  
    if (TAB_CONTACT_POLYG(i)%nature /= 'PLPLx' .or. TAB_CONTACT_POLYG(i)%gap > 0.1E-2) cycle
    if (TAB_CONTACT_POLYG(i)%type == 1) cycle
    
    ! The particles ids
    l_cdt = TAB_CONTACT_POLYG(i)%icdent
    l_ant = TAB_CONTACT_POLYG(i)%ianent
    ! The center of each particle
    vtk_cd_center(:) = TAB_POLYG(l_cdt)%center(:)
    vtk_an_center(:) = TAB_POLYG(l_ant)%center(:)
    ! The normal and tangential vectors
    v_l_normal(:) = TAB_CONTACT_POLYG(i)%n(:)
    v_l_t(:) = TAB_CONTACT_POLYG(i)%t(:)
    
    ! Write! ... 4 vertices for each line
    
    ! Candidate --- First vertex. 
    write(121, '(3(F15.7,A))', advance='no') vtk_cd_center(1)+l_ctc_scale*v_l_t(1), ' ', &
                                             vtk_cd_center(2)+l_ctc_scale*v_l_t(2), ' ', &
                                             0.D0, ' '
    ! Candidate --- Second vertex. 
    write(121, '(3(F15.7,A))', advance='no') vtk_cd_center(1)-l_ctc_scale*v_l_t(1), ' ', &
                                             vtk_cd_center(2)-l_ctc_scale*v_l_t(2), ' ', &
                                             0.D0, ' '
    ! Candidate --- Third vertex. 
    write(121, '(3(F15.7,A))', advance='no') vtk_an_center(1)-l_ctc_scale*v_l_t(1), ' ', &
                                             vtk_an_center(2)-l_ctc_scale*v_l_t(2), ' ', &
                                             0.D0, ' '
    ! Candidate --- 4th vertex. 
    write(121, '(3(F15.7,A))') vtk_an_center(1)+l_ctc_scale*v_l_t(1), ' ', &
                               vtk_an_center(2)+l_ctc_scale*v_l_t(2), ' ', &
                               0.D0, ' '
  end do
  
  ! Writing the conectivity
  write(121,'(A)', advance='no') 'POLYGONS '
  ! 6 faces for each parallepiped
  write(121, '(2(I6,A))') ctc_counter, ' ' , (ctc_counter*5)
  
  curr_l_faces = 0
  
  do i=1, nb_ligneCONTACT_POLYG
    ! Conecting the 4 vertices for each rectangle
    if (TAB_CONTACT_POLYG(i)%nature /= 'PLPLx' .or. TAB_CONTACT_POLYG(i)%gap > 0.1E-2) cycle
    if (TAB_CONTACT_POLYG(i)%type == 1) cycle
    
    !!!!! 1st (1-2-3-4)
    write(121, '(I1,A)', advance= 'no') 4, ' '
    
    if (curr_l_faces .lt. 10) then
      write(121, '(I1)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 100) then
      write(121, '(I2)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 1000) then
      write(121, '(I3)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 10000) then
      write(121, '(I4)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 100000) then
      write(121, '(I5)', advance = 'no') curr_l_faces
    end if
    
    write(121, '(A)', advance='no') ' '
    
    if (curr_l_faces +1 .lt. 10) then
      write(121, '(I1)', advance = 'no') curr_l_faces +1 
    else if (curr_l_faces +1 .lt. 100) then
      write(121, '(I2)', advance = 'no') curr_l_faces +1 
    else if (curr_l_faces +1 .lt. 1000) then
      write(121, '(I3)', advance = 'no') curr_l_faces +1 
    else if (curr_l_faces +1 .lt. 10000) then
      write(121, '(I4)', advance = 'no') curr_l_faces +1
    else if (curr_l_faces +1 .lt. 100000) then
      write(121, '(I5)', advance = 'no') curr_l_faces +1 
    end if
    
    write(121, '(A)', advance='no') ' '
    
    if (curr_l_faces +2 .lt. 10) then
      write(121, '(I1)', advance = 'no') curr_l_faces +2 
    else if (curr_l_faces +2 .lt. 100) then
      write(121, '(I2)', advance = 'no') curr_l_faces +2 
    else if (curr_l_faces +2 .lt. 1000) then
      write(121, '(I3)', advance = 'no') curr_l_faces +2 
    else if (curr_l_faces +2 .lt. 10000) then
      write(121, '(I4)', advance = 'no') curr_l_faces +2
    else if (curr_l_faces +2 .lt. 100000) then
      write(121, '(I5)', advance = 'no') curr_l_faces +2 
    end if
    
    write(121, '(A)', advance='no') ' '
    
    if (curr_l_faces +3 .lt. 10) then
      write(121, '(I1)') curr_l_faces +3
    else if (curr_l_faces +3 .lt. 100) then
      write(121, '(I2)') curr_l_faces +3 
    else if (curr_l_faces +3 .lt. 1000) then
      write(121, '(I3)') curr_l_faces +3 
    else if (curr_l_faces +3 .lt. 10000) then
      write(121, '(I4)') curr_l_faces +3
    else if (curr_l_faces +3 .lt. 100000) then
      write(121, '(I5)') curr_l_faces +3
    end if
    
    curr_l_faces = curr_l_faces + 4
    
  end do
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!! CONTACT FIELDS !!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! How many fields there will be?
  n_l_fields = 1                          ! They are: Status
                                          ! ...
  
  ! A blank space
  write(121,'(A)') ''
  ! The cells begin
  write(121,'(A)', advance='no') 'CELL_DATA'
  
  ! Writing the number of data by field. It corresponds to the same number of forces
  write(121, '(2(I6,A))') ctc_counter
  write(121, '(A,I4)')  'FIELD FieldData', n_l_fields
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Status (3: noctc, 2 Sli, 1 Stick, -1 Wnctc, -2Wsli, -3 Wstck)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! Naming the field, the dimension of the data, and number of lines (as before FIELD), and data type
  write(121,'(A,I1,I6,A)') 'Status ', 1, ctc_counter, ' float'
  do i=1, nb_ligneCONTACT_POLYG
    if (TAB_CONTACT_POLYG(i)%nature /= 'PLPLx' .or. TAB_CONTACT_POLYG(i)%gap > 0.1E-2) cycle
    if (TAB_CONTACT_POLYG(i)%type == 1) cycle
    write(121, '(I2)') TAB_CONTACT_POLYG(i)%status_points
  end do
  
  ! And jump
  write(121, '(A)') ' '

  close(121)
  
  
  print*, 'Drawing in vtk       ---> Ok!'
  
end subroutine draw


!==============================================================================
! Contacts orientation
!==============================================================================
subroutine contact_direction
  
  implicit none
  
  real(kind=8)                                  :: interval
  real(kind=8)                                  :: angcurr, maxint, minint
  integer                                       :: i, j, ninter, contotal, contmean
  real(kind=8), dimension(:,:), allocatable     :: theta_interval
  real(kind=8), dimension(2)                    :: nik
  character(len=27)                             :: nom
  logical                                       :: dir_ctcdir
  
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
  nom       =  './POSTPRO/CTCDIR/CDIR.    '
  
  if (compteur_clout<10) then
    WRITE(nom(23:24),'(I1)')   compteur_clout
  else if ( (compteur_clout>=10) .and. (compteur_clout<100) ) then
    WRITE(nom(23:25),'(I2)')   compteur_clout
  else if ( (compteur_clout>=100).and. (compteur_clout<1000) ) then
    WRITE(nom(23:26),'(I3)')   compteur_clout
  else if ( (compteur_clout>=1000).and. (compteur_clout<10000) ) then
    WRITE(nom(23:27),'(I4)')   compteur_clout
  else 
    print*, "Out of numbers ---> Orientations"
    stop
  end if
  
  ! Number of intervals over pi
  ninter = 18
  
  ! Initializing variables
  contotal = 0
  interval = pi/ninter
  
  ! Allocating the vector that will contain the contact direction and frequency
  if (allocated(theta_interval)) deallocate(theta_interval)
  allocate(theta_interval(ninter,2))
  
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
    do j=1,nb_ligneCONTACT_POLYG
      ! Only between discs
      if(TAB_CONTACT_POLYG(j)%nature /='PLPLx') cycle
      ! Only active contacts
      if(abs(TAB_CONTACT_POLYG(j)%rn) .le. 1.D-8) cycle
      ! The normal vector
      nik(1)  = TAB_CONTACT_POLYG(j)%n(1)
      nik(2)  = TAB_CONTACT_POLYG(j)%n(2)
      
      ! The 1 and 2 quadrants only
      ! Changing the normal vector direction if necessary
      if(nik(2) .lt. 0) then
        nik(1)  = -nik(1)
        nik(2)  = -nik(2)
      end if
      
      ! Finding the current angle. Note the vector nik is already unitary!
      angcurr = acos(nik(1))
      
      ! Checking the interval
      if(angcurr .lt. maxint .and. angcurr .ge. minint) then
        theta_interval(i,2) = theta_interval(i,2) + 1
      end if
    end do
  end do
  
  ! Computing the number of active contacts
  do i=1, ninter
    contotal = contotal + theta_interval(i,2)
  end do
  
  ! Normalizing
  theta_interval(:,2) = theta_interval(:,2)/(contotal*interval)
  
  ! Opening the file
  open(unit=110,file=nom,status='replace')
  
  ! Writing the heading 
  write(110,*) '# Theta(Rad) ', '  Freq_norm  '
  
  ! Writing
  do i=1, ninter
    write(110,'(2(1X,E12.5))') theta_interval(i,1), theta_interval(i,2)
  end do
  
  ! Writing the 3rd and 4th quadrants
  do i=1, ninter
    write(110,'(2(1X,E12.5))') theta_interval(i,1)+pi, theta_interval(i,2)
  end do
  
  close(110)
  
  print*, 'Write Contacts dir              ---> Ok!'
  
end subroutine contact_direction

!==============================================================================
! Contact forces lists for PDF and more
!==============================================================================
subroutine f_list

  implicit none

  integer                                       :: i
  real*8                                        :: Rnik, Rtik
  character(len=26)                             :: file_c
  logical                                       :: dir_ctcdir

  ! Cleaning or creating the folder if necessary
  if (first_over_all) then
    ! Asking if the file already exists
    inquire(file='./POSTPRO/F_LIST', exist=dir_ctcdir)
    if(dir_ctcdir) then
      ! Cleaning
      call system('rm ./POSTPRO/F_LIST/*')
    else
      ! Creating
      call system('mkdir ./POSTPRO/F_LIST')
    end if
  end if

  ! The file name
  file_c       =  './POSTPRO/F_LIST/F_L.    '

  if (compteur_clout<10) then
    WRITE(file_c(22:23),'(I1)')   compteur_clout
  else if ( (compteur_clout>=10) .and. (compteur_clout<100) ) then
    WRITE(file_c(22:24),'(I2)')   compteur_clout
  else if ( (compteur_clout>=100).and. (compteur_clout<1000) ) then
    WRITE(file_c(22:25),'(I3)')   compteur_clout
  else if ( (compteur_clout>=1000).and. (compteur_clout<10000) ) then
    WRITE(file_c(22:26),'(I4)')   compteur_clout
  end if

  ! Opening the file
  open(unit=112,file=file_c,status='replace')

  ! Writing the heading
  write(112,*) '     FN     ', '     FT      '

  do i=1, nb_ligneCONTACT_POLYG
    ! Only contacts between particles
    if (TAB_CONTACT_POLYG(i)%nature == 'PLJCx') cycle
    Rnik = TAB_CONTACT_POLYG(i)%rt
    Rtik = TAB_CONTACT_POLYG(i)%rn
    ! Only forces above a given threshold
    if (option_cohe == 1) then
      if (abs(Rnik) .lt. 1e-10 .and. abs(Rtik) .lt. 1e-10) cycle
    else
      if (Rnik .lt. 1e-10) cycle
    end if

    write(112,'(2(1X,E12.5))') Rnik, Rtik
  end do

  close(112)

  print*, 'Write Contacts Forces           ---> Ok!'

end subroutine f_list

!================================================
! Estimating the failure type of contacts
!================================================
subroutine failure_mode

  implicit none

  integer                                  ::  i, j
  integer                                  ::  n_tension, n_shear
  real*8                                   ::  CNCT_ratio
  logical                                  ::  cur_broken

  n_tension = 0
  n_shear = 0

  ! Careful! Setting the CN/CT ratio in hard 
  CNCT_ratio = 1.

  !print*, 'A1'
  !Finding the number of initial Cohesive contacts
  if (first_over_all) then
    n_ctc_cohe = 0
    do i=1, nb_ligneCONTACT
      if (TAB_CONTACT(i)%status(1:1) == "W") then
        n_ctc_cohe = n_ctc_cohe + 1
      end if
    end do

    ! Allocating array of cohe couples
    if (allocated(cohe_couples)) deallocate(cohe_couples)
    allocate(cohe_couples(n_ctc_cohe,2))

    ! Allocating the array for the local forces
    if (allocated(l_forces_c)) deallocate(l_forces_c)
    allocate(l_forces_c(n_ctc_cohe,3))

    ! Allocating the array for the type of failure
    if (allocated(t_failure)) deallocate(t_failure)
    allocate(t_failure(n_ctc_cohe))

    ! Initializing all that stuff
    l_forces_c(:,:) = 0.D0
    t_failure(:) = 0

    j=0
    do i=1, nb_ligneCONTACT
      if (TAB_CONTACT(i)%status(1:1) == "W") then
        j = j+1
        cohe_couples(j,1) = TAB_CONTACT(i)%icdent
        cohe_couples(j,2) = TAB_CONTACT(i)%ianent

        ! Status, normal and tangential forces
        l_forces_c(j,1) = 1.D0
        l_forces_c(j,2) = TAB_CONTACT(i)%rn
        l_forces_c(j,3) = TAB_CONTACT(i)%rt
      end if
    end do
  end if

  ! Searching broken contacts
  do i=1, n_ctc_cohe
    ! Checking the remaining cohesive contacts
    if (t_failure(i)==0) then
      ! Assuming this contact has been broken
      cur_broken = .true.
      ! Search if the contact still exists in the list
      do j=1, nb_ligneCONTACT
        if (TAB_CONTACT(j)%icdent == cohe_couples(i,1) .and. TAB_CONTACT(j)%ianent == cohe_couples(i,2)) then
          if (TAB_CONTACT(j)%status(1:1) == "W") then
            ! The contact exists and it's Cohe
            cur_broken = .false.

            ! The contact remains and forces are updated
            l_forces_c(i,2) = TAB_CONTACT(j)%rn
            l_forces_c(i,3) = TAB_CONTACT(j)%rt
          end if
        end if
      end do

      ! Checking cur_broken
      ! If it was broken
      if (cur_broken) then
        l_forces_c(i,1) = 0.D0
        ! Checking the most probable mode it was broken
        if (l_forces_c(i,2) .gt. abs(l_forces_c(i,3))) then
          ! It was by tension
          t_failure(i) = 1
        else
          ! It was by shear
          t_failure(i) = 2
        end if
      end if
    end if
  end do

  ! Counting the type of contacts broken
  do i=1, n_ctc_cohe
    if (t_failure(i)==1) then
      n_tension = n_tension + 1
    else if (t_failure(i)==2) then
      n_shear = n_shear + 1
    end if
  end do

  if(first_over_all) then
    write(113,*) '#   time     ', ' n_tension  ', '   n_shear   '
  endif

  write(113,'(1(1X,E12.5),2(1X,I8))') time, n_tension, n_shear

  print*, 'Write Failure Mode              ---> Ok!'
end subroutine failure_mode


!================================================
! Listing the interactions
!================================================
subroutine list_interact

  implicit none

  integer                                  :: i, cd_loc, an_loc
  real*8                                   :: Lnik, Ltik
  real*8, dimension(2)                     :: Lik, nik, tik
  character(len=29)                        :: file_c
  logical                                  :: dir_ctcdir

  ! Initializing
  Lik(:) = 0
  Lnik = 0
  Ltik = 0

  ! Cleaning or creating the folder if necessary
  if (first_over_all) then
    ! Asking if the file already exists
    inquire(file='./POSTPRO/INTER_LIST', exist=dir_ctcdir)
    if(dir_ctcdir) then
      ! Cleaning
      call system('rm ./POSTPRO/INTER_LIST/*')
    else
      ! Creating
      call system('mkdir ./POSTPRO/INTER_LIST')
    end if
  end if

  ! The file name
  file_c       =  './POSTPRO/INTER_LIST/F_L.    '

  if (compteur_clout<10) then
    write(file_c(26:26),'(I1)')   compteur_clout
  else if ( (compteur_clout>=10) .and. (compteur_clout<100) ) then
    write(file_c(26:27),'(I2)')   compteur_clout
  else if ( (compteur_clout>=100).and. (compteur_clout<1000) ) then
    write(file_c(26:28),'(I3)')   compteur_clout
  else if ( (compteur_clout>=1000).and. (compteur_clout<10000) ) then
    write(file_c(26:29),'(I4)')   compteur_clout
  else
    print*, "The list interact cannot handle more files"
    stop
  end if

  ! Opening the file
  open(unit=114,file=file_c,status='replace')

  ! Writing the heading
  write(114,*) '   CD     ', '   AN     ', '   CTC_LEN   ', &
               '      FN     ', '      FT     ', '    GAP_0    ', '     GAP     ', &
               '    DISP_T   ', '      LN     ', '      LT     '!, '   STATUS    '
  do i=1, nb_ligneCONTACT_POLYG
    if (TAB_CONTACT_POLYG(i)%nature == 'PLJCx') cycle
    cd_loc  = TAB_CONTACT_POLYG(i)%icdent
    an_loc  = TAB_CONTACT_POLYG(i)%ianent
    nik(1)  = TAB_CONTACT_POLYG(i)%n(1)
    nik(2)  = TAB_CONTACT_POLYG(i)%n(2)
    tik(1)  = TAB_CONTACT_POLYG(i)%t(1)
    tik(2)  = TAB_CONTACT_POLYG(i)%t(2)

    if (TAB_CONTACT_POLYG(i)%rn .lt. 1e-9 .and. TAB_CONTACT_POLYG(i)%rt .lt. 1e-9) cycle

    ! The branch vector
    Lik(1) = TAB_POLYG(cd_loc)%center(1) - TAB_POLYG(an_loc)%center(1)
    Lik(2) = TAB_POLYG(cd_loc)%center(2) - TAB_POLYG(an_loc)%center(2)
    
    ! The components
    Lnik = Lik(1)*nik(1)+Lik(2)*nik(2)
    Ltik = Lik(1)*tik(1)+Lik(2)*tik(2)

    write(114,'(2(1X,I9),8(1X,E12.5))') TAB_CONTACT_POLYG(i)%icdent, TAB_CONTACT_POLYG(i)%ianent, &
                                                  TAB_CONTACT_POLYG(i)%cd_len, &
                                                  TAB_CONTACT_POLYG(i)%rn, TAB_CONTACT_POLYG(i)%rt, &
                                                  TAB_CONTACT_POLYG(i)%gap0, TAB_CONTACT_POLYG(i)%gap, &
                                                  TAB_CONTACT_POLYG(i)%tang_disp, Lnik, Ltik!, &
                                                  !TAB_CONTACT_POLYG(i)%law!, TAB_CONTACT_POLYG(i)%status
  end do

  close(114)

  print*, 'Write Interaction List          ---> Ok!'
end subroutine list_interact


!================================================
! Listing the interactions on frame L
!================================================
subroutine list_interact_l

  implicit none

  integer                                  :: i, cd_loc, an_loc
  real*8                                   :: Lnik, Ltik, Rnik, Rtik
  real*8                                   :: Force_N, Force_T, Norme_LN, Norme_LT
  real*8, dimension(2)                     :: Lik, nik, tik, Fik, nik_c, tik_c
  character(len=31)                        :: file_c
  logical                                  :: dir_ctcdir

  ! Initializing
  Lik(:) = 0
  Lnik = 0
  Ltik = 0

  ! Cleaning or creating the folder if necessary
  if (first_over_all) then
    ! Asking if the file already exists
    inquire(file='./POSTPRO/INTER_LIST_L', exist=dir_ctcdir)
    if(dir_ctcdir) then
      ! Cleaning
      call system('rm ./POSTPRO/INTER_LIST_L/*')
    else
      ! Creating
      call system('mkdir ./POSTPRO/INTER_LIST_L')
    end if
  end if

  ! The file name
  file_c       =  './POSTPRO/INTER_LIST_L/F_L.    '

  if (compteur_clout<10) then
    write(file_c(28:28),'(I1)')   compteur_clout
  else if ( (compteur_clout>=10) .and. (compteur_clout<100) ) then
    write(file_c(28:29),'(I2)')   compteur_clout
  else if ( (compteur_clout>=100).and. (compteur_clout<1000) ) then
    write(file_c(28:30),'(I3)')   compteur_clout
  else if ( (compteur_clout>=1000).and. (compteur_clout<10000) ) then
    write(file_c(28:31),'(I4)')   compteur_clout
  else
    print*, "The list interact cannot handle more files"
    stop
  end if

  ! Opening the file
  open(unit=126,file=file_c,status='replace')

  ! Writing the heading
  write(126,*) '   CD     ', '   AN     ', '   CTC_LEN   ', &
               '      FN     ', '      FT     ', '    GAP_0    ', '     GAP     ', &
               '    DISP_T   ', '      LN     ', '      LT     '!, '   STATUS    '
  do i=1, nb_ligneCONTACT_POLYG
    if (TAB_CONTACT_POLYG(i)%nature == 'PLJCx') cycle
    cd_loc  = TAB_CONTACT_POLYG(i)%icdent
    an_loc  = TAB_CONTACT_POLYG(i)%ianent
    nik_c(1)  = TAB_CONTACT_POLYG(i)%n(1)
    nik_c(2)  = TAB_CONTACT_POLYG(i)%n(2)
    tik_c(1)  = TAB_CONTACT_POLYG(i)%t(1)
    tik_c(2)  = TAB_CONTACT_POLYG(i)%t(2)
    Rnik  = TAB_CONTACT_POLYG(i)%rn
    Rtik  = TAB_CONTACT_POLYG(i)%rt

    if (TAB_CONTACT_POLYG(i)%rn .lt. 1e-9 .and. TAB_CONTACT_POLYG(i)%rt .lt. 1e-9) cycle

    ! The branch vector
    Lik(1) = TAB_POLYG(cd_loc)%center(1) - TAB_POLYG(an_loc)%center(1)
    Lik(2) = TAB_POLYG(cd_loc)%center(2) - TAB_POLYG(an_loc)%center(2)

        ! The new frame on L
    nik(1) = Lik(1)/(Lik(1)**2 + Lik(2)**2)**0.5
    nik(2) = Lik(2)/(Lik(1)**2 + Lik(2)**2)**0.5
    tik(1) = nik(2)
    tik(2) = -nik(1)

    ! Building the stress tensor (on the contact frame)
    Fik(1) = (Rnik*nik_c(1)+Rtik*tik_c(1))
    Fik(2) = (Rnik*nik_c(2)+Rtik*tik_c(2))

    ! Forces on the L frame 
    Force_N = Fik(1)*nik(1) + Fik(2)*nik(2)
    Force_T = Fik(1)*tik(1) + Fik(2)*tik(2)

    Norme_LN = (Lik(1)*nik(1) + Lik(2)*nik(2))
    Norme_LT = (Lik(1)*tik(1) + Lik(2)*tik(2))
    
    ! The components
    Lnik = Lik(1)*nik(1)+Lik(2)*nik(2)
    Ltik = Lik(1)*tik(1)+Lik(2)*tik(2)

    write(126,'(2(1X,I9),8(1X,E12.5))') TAB_CONTACT_POLYG(i)%icdent, TAB_CONTACT_POLYG(i)%ianent, &
                                                  TAB_CONTACT_POLYG(i)%cd_len, &
                                                  Force_N, Force_T, &
                                                  TAB_CONTACT_POLYG(i)%gap0, TAB_CONTACT_POLYG(i)%gap, &
                                                  TAB_CONTACT_POLYG(i)%tang_disp, Norme_LN, Norme_LT !, &
                                                  !TAB_CONTACT_POLYG(i)%law!, TAB_CONTACT_POLYG(i)%status
  end do

  close(126)

  print*, 'Write Interaction List L       ---> Ok!'
end subroutine list_interact_l

!==================================================================================================
!Procedure pour calcul des contacts persistant ...
!==================================================================================================
  subroutine persiste_contact
  implicit none
    integer                           :: i,j,cpt

     cpt        = 0
     if (first_list_contact) then
       nb_ligneCONTACT_first = nb_ligneCONTACT
       allocate(TAB_CONTACT_persiste(nb_ligneCONTACT_first))
       do i=1,nb_ligneCONTACT_first
         TAB_CONTACT_persiste(i)%icdent      = TAB_CONTACT(i)%icdent
         TAB_CONTACT_persiste(i)%ianent      = TAB_CONTACT(i)%ianent
         TAB_CONTACT_persiste(i)%n(1)        = TAB_CONTACT(i)%n(1)
         TAB_CONTACT_persiste(i)%n(2)        = TAB_CONTACT(i)%n(2)
         TAB_CONTACT_persiste(i)%n(3)        = TAB_CONTACT(i)%n(3)
         TAB_CONTACT_persiste(i)%coor_ctc(1) = TAB_CONTACT(i)%coor_ctc(1)
         TAB_CONTACT_persiste(i)%coor_ctc(2) = TAB_CONTACT(i)%coor_ctc(2)
         TAB_CONTACT_persiste(i)%coor_ctc(3) = TAB_CONTACT(i)%coor_ctc(3)
         TAB_CONTACT_persiste(i)%rn          = TAB_CONTACT(i)%rn
         TAB_CONTACT_persiste(i)%rt          = TAB_CONTACT(i)%rt
         TAB_CONTACT_persiste(i)%rs          = TAB_CONTACT(i)%rs
       end do
       first_list_contact = .false.
     end if
     
     do i=1,nb_ligneCONTACT_first
       if (TAB_CONTACT(i)%nature == 'PLJCx') cycle
       do j=1,nb_ligneCONTACT
         if ((TAB_CONTACT_persiste(i)%icdent == TAB_CONTACT(j)%icdent).and. &
             (TAB_CONTACT_persiste(i)%ianent == TAB_CONTACT(j)%ianent)) then
           cpt = cpt + 1 
         end if
       end do
     end do

     !if (type_box==1) write(106,'(3(1X,D12.5))') time, real(cpt,8)
     !if (type_box==0) write(106,'(4(1X,D12.5))') time,height,real(cpt,8)
   end subroutine persiste_contact


!==============================================================================
! Calcul de l'anisotropie des contacts
!==============================================================================
subroutine anisotropy_contact
  
  implicit none
  
  real(kind=8),dimension(2,2)              :: Fabric
  real(kind=8),dimension(2)                :: nik
  integer                                  :: i,cd,an
  integer                                  :: ierror,matz,lda
  real(kind=8),dimension(2)                :: wr,wi
  real(kind=8)                             :: S1,S2,cpt, Rnik
  real(kind=8),dimension(2,2)              :: localframe
  
  ! Initializing the number of active contacts
  cpt = 0
  
  ! Initializing the fabric tensor
  Fabric(:,:) = 0.
  
  ! Building the fabric tensor
  do i=1,nb_ligneCONTACT
    ! Checking it is a contact between two bodies
    if(TAB_CONTACT_POLYG(i)%nature == 'PLJCx') cycle
    cd      = TAB_CONTACT_POLYG(i)%icdent
    an      = TAB_CONTACT_POLYG(i)%ianent
    nik(1)  = TAB_CONTACT_POLYG(i)%n(1)
    nik(2)  = TAB_CONTACT_POLYG(i)%n(2)
    Rnik    = TAB_CONTACT_POLYG(i)%rn
    
    ! Checking the normal force is not zero
    if (abs(Rnik) .le. 1.D-8) cycle
    
    Fabric(1,1:2) = nik(1)*nik(1:2) + Fabric(1,1:2)
    Fabric(2,1:2) = nik(2)*nik(1:2) + Fabric(2,1:2)
    cpt = cpt + 1
  end do
  
  ! Normalizing by the number of contacts
  Fabric = Fabric / cpt
  
  ! Finding the eigen values
  lda  = 2
  matz = 1
  
  call rg ( lda, 2, Fabric, wr, wi, matz, localframe, ierror )
  S1 = max( wr(1),wr(2) )
  S2 = min( wr(1),wr(2) )
  
  if (first_over_all) then
    write(104,*) '#   time     ', '     ac      '
  end if
  
  write(104,'(2(1X,E12.5))') time, 2*(S1-S2)
  
  print*, 'Write Anisotropy Conctacts      ---> Ok!'
  
end subroutine anisotropy_contact


!==============================================================================
! Calcul de l'anisotropie des contacts
!==============================================================================
subroutine anisotropy_force
  
  implicit none
  
  real(kind=8), dimension(2,2)                 :: Fabric_N,Fabric_T, Fabric_F
  real(kind=8), dimension(2)                   :: nik,tik
  integer                                      :: i,cd,an,j
  integer                                      :: ierror,matz,lda
  real(kind=8), dimension(2)                   :: wrn, win, wrf, wif
  real(kind=8)                                 :: X1n,X2n,X1t,X2t, X1f,X2f,cpt,Rtik,Rnik,av_force
  real(kind=8), dimension(2,2)                 :: localframeN, localframeF
  real(kind=8), dimension(2)                   :: Fik
  
  ! Initializing variables
  Fabric_N(:,:) = 0.
  Fabric_T(:,:) = 0.
  Fabric_F(:,:) = 0.
  cpt = 0.
  av_force = 0.
  Fik(:) = 0.
  
  ! Building the fabric forces tensor
  do i=1,nb_ligneCONTACT_POLYG
    if  (TAB_CONTACT_POLYG(i)%nature == 'PLJCx') cycle
    cd      = TAB_CONTACT_POLYG(i)%icdent
    an      = TAB_CONTACT_POLYG(i)%ianent
    nik(1)  = TAB_CONTACT_POLYG(i)%n(1)
    nik(2)  = TAB_CONTACT_POLYG(i)%n(2)
    tik(1)  = TAB_CONTACT_POLYG(i)%t(1)
    tik(2)  = TAB_CONTACT_POLYG(i)%t(2)
    Rnik    = TAB_CONTACT_POLYG(i)%rn
    Rtik    = TAB_CONTACT_POLYG(i)%rt
    
    ! Only active contacts
    if (abs(Rnik) .le. 1.D-8) cycle
    
    ! Active contacts
    cpt     = cpt + 1
    
    ! The force vector
    Fik(:) = Rnik*nik(:) + Rtik*tik(:)
    
    ! Average normal force
    av_force = av_force + Fik(1)*nik(1) + Fik(2)*nik(2)
    
    ! Building the tensor of normal forces
    Fabric_N(1,1:2) = Rnik*nik(1)*nik(1:2) + Fabric_N(1,1:2)
    Fabric_N(2,1:2) = Rnik*nik(2)*nik(1:2) + Fabric_N(2,1:2)
    
    ! Building the tensor of tangential forces
    Fabric_T(1,1:2) = Rtik*nik(1)*tik(1:2) + Fabric_T(1,1:2)
    Fabric_T(2,1:2) = Rtik*nik(2)*tik(1:2) + Fabric_T(2,1:2)
    
  end do
  
  ! Computing the average normal force
  av_force = av_force / cpt
  
  ! Computing the fabric forces tensors
  Fabric_N = Fabric_N / (cpt*av_force)
  Fabric_T = Fabric_T / (cpt*av_force)
  
  ! Building the full fabric forces tensor (Built this always before using rg!!!)
  Fabric_F = Fabric_N + Fabric_T
  
  ! Finding the eigen values of the fabric normal forces tensor
  lda = 2
  matz = 1
  
  call rg(lda, 2, Fabric_N, wrn, win, matz, localframeN, ierror)
  
  X1n = max(wrn(1),wrn(2))
  X2n = min(wrn(1),wrn(2))
  
  ! Finding the eigen values of the full fabric forces tensor
  call rg(lda, 2, Fabric_F, wrf, wif, matz, localframeF, ierror)
  X1f = max(wrf(1),wrf(2))
  X2f = min(wrf(1),wrf(2))
  
  if (first_over_all) then
     write(105,*) '#   time     ', '2(X1-2)/(X1+2)n', ' 2(X1-2)/(X1+2)f '
  end if 
  
  write(105,'(3(1X,E12.5))') time, 2*(X1n-X2n)/(X1n+X2n),2*(X1f-X2f)/(X1f+X2f)
  
  print*, 'Write Anisotropy Force          ---> Ok!'
end subroutine anisotropy_force


!==============================================================================
! Calcul de l'anisotropie des branches
!==============================================================================
subroutine anisotropy_branch
  
  implicit none
  
  real(kind=8), dimension(2,2)                 :: Fabric_N,Fabric_T, Fabric_F
  real(kind=8), dimension(2)                   :: nik,tik
  integer                                      :: i,cd,an,j
  integer                                      :: ierror,matz,lda
  real(kind=8), dimension(2)                   :: wrn, win, wrf, wif
  real(kind=8)                                 :: X1n,X2n,X1t,X2t, X1f,X2f,cpt,av_length, Lnik, Ltik, Rnik
  real(kind=8), dimension(2,2)                 :: localframeN, localframeF
  real(kind=8), dimension(2)                   :: Lik
  
  ! Initializing variables
  Fabric_N(:,:) = 0.
  Fabric_T(:,:) = 0.
  Fabric_F(:,:) = 0.
  cpt = 0.
  av_length = 0.
  Lik(:) = 0.
  
  ! Building the chi tensor with the length of branches 
  do i=1,nb_ligneCONTACT_POLYG
    if  (TAB_CONTACT_POLYG(i)%nature == 'PLJCx') cycle
    cd      = TAB_CONTACT_POLYG(i)%icdent
    an      = TAB_CONTACT_POLYG(i)%ianent
    nik(1)  = TAB_CONTACT_POLYG(i)%n(1)
    nik(2)  = TAB_CONTACT_POLYG(i)%n(2)
    tik(1)  = TAB_CONTACT_POLYG(i)%t(1)
    tik(2)  = TAB_CONTACT_POLYG(i)%t(2)
    Rnik    = TAB_CONTACT_POLYG(i)%rn
    
    ! Only active contacts
    if (abs(Rnik) .le. 1.D-8) cycle
    
    ! Active contacts
    cpt     = cpt + 1
    
    ! The branch vector
    Lik(1) = TAB_POLYG(cd)%center(1) - TAB_POLYG(an)%center(1)
    Lik(2) = TAB_POLYG(cd)%center(2) - TAB_POLYG(an)%center(2)
    
    ! Average branch length
    av_length = av_length + (Lik(1)**2 + Lik(2)**2)**0.5
    
    Lnik = Lik(1)*nik(1)+Lik(2)*nik(2)
    Ltik = Lik(1)*tik(1)+Lik(2)*tik(2)
    
    ! Building the chi_n
    Fabric_N(1,1:2) = Lnik*nik(1)*nik(1:2) + Fabric_N(1,1:2)
    Fabric_N(2,1:2) = Lnik*nik(2)*nik(1:2) + Fabric_N(2,1:2)
    
    ! Building the tensor of tangential forces
    Fabric_T(1,1:2) = Ltik*nik(1)*tik(1:2) + Fabric_T(1,1:2)
    Fabric_T(2,1:2) = Ltik*nik(2)*tik(1:2) + Fabric_T(2,1:2)
  end do
  
  ! Computing the average normal force
  av_length = av_length / cpt
  
  ! Computing the fabric forces tensors
  Fabric_N = Fabric_N / (cpt*av_length)
  Fabric_T = Fabric_T / (cpt*av_length)
  
  ! Building the full fabric forces tensor (Built this always before using rg!!!)
  Fabric_F = Fabric_N + Fabric_T
  
  ! Finding the eigen values of the fabric normal forces tensor
  lda = 2
  matz = 1
  
  call rg(lda, 2, Fabric_N, wrn, win, matz, localframeN, ierror)
  
  X1n = max(wrn(1),wrn(2))
  X2n = min(wrn(1),wrn(2))
  
  ! Finding the eigen values of the full fabric forces tensor
  call rg(lda, 2, Fabric_F, wrf, wif, matz, localframeF, ierror)
  X1f = max(wrf(1),wrf(2))
  X2f = min(wrf(1),wrf(2))
  
  if (first_over_all) then
     write(106,*) '#   time     ', '2(X1-2)/(X1+2)n', ' 2(X1-2)/(X1+2)f '
  end if 
  
  write(106,'(3(1X,E12.5))') time, 2*(X1n-X2n)/(X1n+X2n),2*(X1f-X2f)/(X1f+X2f)
  
  print*, 'Write Anisotropy Branch          ---> Ok!'
end subroutine anisotropy_branch


!==============================================================================
! Branch orientation
!==============================================================================
subroutine branch_dir

  implicit none

  integer                                      :: i, j, ninter, contotal, cd, an
  real(kind=8)                                 :: interval, pro_n, pro_t
  real(kind=8)                                 :: angcurr, maxint, minint
  real(kind=8), dimension(3)                   :: nik, Lik, tik, sik, tanik
  real(kind=8), dimension(:,:), allocatable    :: theta_interval
  character(len=26)                            :: file_name
  logical                                      :: dir_ctcdir

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
  file_name       =  './POSTPRO/BRCDIR/BDIR.    '

  if (compteur_clout<10) then
    WRITE(file_name(23:23),'(I1)')   compteur_clout
  else if ( (compteur_clout>=10) .and. (compteur_clout<100) ) then
    WRITE(file_name(23:24),'(I2)')   compteur_clout
  else if ( (compteur_clout>=100).and. (compteur_clout<1000) ) then
    WRITE(file_name(23:25),'(I3)')   compteur_clout
  else if ( (compteur_clout>=1000).and. (compteur_clout<10000) ) then
    WRITE(file_name(23:26),'(I4)')   compteur_clout
  end if

  ! Opening the file
  open(unit=115,file=file_name,status='replace')

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
    do j=1,nb_ligneCONTACT_POLYG
      cd = TAB_CONTACT_POLYG(j)%icdent
      an = TAB_CONTACT_POLYG(j)%ianent
      ! Only between discs
      if(TAB_CONTACT_POLYG(j)%nature == 'PLJCx') cycle

      ! Only active contacts
      if(abs(TAB_CONTACT_POLYG(j)%rn) .lt. 1.D-9) cycle

      ! The branch vector
      Lik(:) = TAB_POLYG(cd)%center(:)-TAB_POLYG(an)%center(:)

      ! Not joining contacts in the periodic zones
      if (sqrt(Lik(1)**2 + Lik(2)**2) .gt. 3.*(TAB_POLYG(cd)%radius + TAB_POLYG(an)%radius)) cycle

      ! The normal vector
      nik  = TAB_CONTACT_POLYG(j)%n

      ! The tangential vector - xz plane!
      !tik = TAB_CONTACT(j)%t
      !sik = TAB_CONTACT(j)%s
      !tanik = sik*TAB_CONTACT(j)%rs + tik*TAB_CONTACT(j)%rt

      ! unitary
      !tanik = tanik/sqrt(tanik(1)**2 + tanik(2)**2 + tanik(3)**2)

      ! Normal projection
      pro_n = abs(Lik(1)*nik(1) + Lik(2)*nik(2))

      ! The 1 and 2 quadrants only
      ! Changing the vector direction if necessary
      if(nik(2) .lt. 0) then
        nik(1) = -nik(1)
        nik(2) = -nik(2)
      end if

      tanik(1) = nik(2)
      tanik(2) = -nik(1)

      ! Tangential projection
      pro_t = Lik(1)*tanik(1) + Lik(2)*tanik(2)

      !if (pro_n .lt. 0.6*pro_t) cycle

      ! Finding the current angle. Note the vector nik is unitary!
      angcurr = acos(nik(1)/sqrt(nik(1)**2 + nik(2)**2))

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
  write(115,*) '# Theta(Rad) ', '    <ln>    ', '    <lt>    '

  ! Writing
  do i=1, ninter
    write(115,'(3(1X,E12.5))') theta_interval(i,1), theta_interval(i,3), theta_interval(i,4)
  end do

  ! Writing the 3rd and 4th quadrants
  do i=1, ninter
    write(115,'(3(1X,E12.5))') theta_interval(i,1)+pi, theta_interval(i,3), theta_interval(i,4)
  end do

  close(115)

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
  file_name       =  './POSTPRO/FRCDIR/FDIR.    '

  if (compteur_clout<10) then
    WRITE(file_name(23:23),'(I1)')   compteur_clout
  else if ( (compteur_clout>=10) .and. (compteur_clout<100) ) then
    WRITE(file_name(23:24),'(I2)')   compteur_clout
  else if ( (compteur_clout>=100).and. (compteur_clout<1000) ) then
    WRITE(file_name(23:25),'(I3)')   compteur_clout
  else if ( (compteur_clout>=1000).and. (compteur_clout<10000) ) then
    WRITE(file_name(23:26),'(I4)')   compteur_clout
  end if

  ! Opening the file
  open(unit=116,file=file_name,status='replace')

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
    do j=1,nb_ligneCONTACT_POLYG
      cd = TAB_CONTACT_POLYG(j)%icdent
      an = TAB_CONTACT_POLYG(j)%ianent
      ! Only between discs
      if(TAB_CONTACT_POLYG(j)%nature == 'PLJCx') cycle

      ! Only active contacts
      if(abs(TAB_CONTACT_POLYG(j)%rn) .lt. 1.D-9) cycle

      Rnik = TAB_CONTACT_POLYG(j)%rn
      Rtik = TAB_CONTACT_POLYG(j)%rt

      ! The normal vector
      nik = TAB_CONTACT_POLYG(j)%n

      ! The tangential vector
      tik = TAB_CONTACT_POLYG(j)%t
      ! The real tangential vector - in xz plane!
      !tanik(1:3) = Rsik*sik(1:3) + Rtik*tik(1:3)
      tanik(1) = nik(2)
      tanik(2) = -nik(1)

      ! The force vector
      Fik(1:2) = (Rnik*nik(1:2)+Rtik*tik(1:2))

      ! unitary!
      !tanik = tanik/sqrt(tanik(1)**2 + tanik(2)**2 + tanik(3)**2)

      ! Normal projection
      pro_n = abs(Fik(1)*nik(1) + Fik(2)*nik(2))

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
      pro_t = Fik(1)*tanik(1) + Fik(2)*tanik(2)
      !print*, 'ini'
      !print*, pro_t
      !print*, Rtik
      !print*, Rsik
      !print*, sqrt(Rsik**2+Rtik**2)
      !print*, 'asdasd'

      ! The 1 and 2 quadrants only
      ! Changing the vector direction if necessary
      if(nik(2) .lt. 0) then
        nik(1) = -nik(1)
        nik(2) = -nik(2)
      end if

      ! Finding the current angle. Note the vector nik is already unitary!
      angcurr = acos(nik(1)/sqrt(nik(1)**2 + nik(2)**2))

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
  write(116,*) '# Theta(Rad) ', '    <fn>    ', '    <ft>    '

  ! Writing
  do i=1, ninter
    write(116,'(3(1X,E12.5))') theta_interval(i,1), theta_interval(i,3), theta_interval(i,4)
  end do

  ! Writing the 3rd and 4th quadrants
  do i=1, ninter
    write(116,'(3(1X,E12.5))') theta_interval(i,1)+pi, theta_interval(i,3), theta_interval(i,4)
  end do

  close(116)

  print*, 'Write Forces dir    ---> Ok!'

end subroutine forces_dir

!==============================================================================
! Computing the probability of N number of contacts 
!==============================================================================
subroutine ctc_probability
  
  implicit none 
  
  integer                                          :: i, j, temp_counter, an, cd
  real (kind=8), dimension(20)                     :: prob_n_ctc
  
  
  ! Initializing variables
  prob_n_ctc(:) = 0
  
  ! For all the particles
  do i=1, n_particles
    temp_counter = 0
    ! Going to find this particle in the contact list
    
    do j=1, nb_ligneCONTACT_POLYG
      
      ! For the active contacts
      if (TAB_CONTACT_POLYG(j)%rn .le. 0.D0) cycle
      if (TAB_CONTACT_POLYG(j)%nature == 'PLJCx') cycle
      
      cd = TAB_CONTACT_POLYG(j)%icdent
      an = TAB_CONTACT_POLYG(j)%ianent
      if (cd == i .or. an == i) then
        temp_counter = temp_counter + 1
      end if
    end do
    prob_n_ctc(temp_counter + 1) = prob_n_ctc(temp_counter + 1) + 1
  end do
  
  prob_n_ctc(:) = prob_n_ctc(:)/n_particles
  
  if (first_over_all) then
     write(111,*) '#   time     ', '    Pc(0)    ', '    Pc(1)    ', '    Pc(2)    ', '    Pc(3)    ', '    Pc(4)    ', &
                  '    Pc(5)    ', '    Pc(6)    ', '    Pc(7)    ', '    Pc(8)    ', '    Pc(9)    ', '    Pc(10)    '
  end if 
  
  write(111,'(12(1X,E12.5))') time, prob_n_ctc(1), prob_n_ctc(2), prob_n_ctc(3), prob_n_ctc(4), prob_n_ctc(5), prob_n_ctc(6),&
                                   prob_n_ctc(7), prob_n_ctc(8), prob_n_ctc(9), prob_n_ctc(10), prob_n_ctc(11) 
  
  print*, 'Write Contact Number Probability---> Ok!'
  
end subroutine ctc_probability


!==============================================================================
! Calcul des ansiotropies de contacts signees a partir des tenseurs H
!==============================================================================
subroutine signed_anisotropy

  implicit none

  integer                                             ::  i,j,cd,an,Nsect_orientation
  real(kind=8)                                        ::  ac,aln,alt,afn,aft,&
                                                          ac_signe,aln_signe,alt_signe,afn_signe,aft_signe,&
                                                          theta_c,theta_ln,theta_lt,theta_fn,theta_ft, theta_sigma, &
                                                          double_product,double_product_sign

  real(kind=8)                                        ::  Force_N,Force_T,Norme_LN,Norme_LT,alpha, Rtik, Rnik, &
                                                          sect,cpt,mean_FN,mean_LN,mean_LT,mean_FT

  real(kind=8),dimension(2)                           ::  nik,tik,Lik,nik_alpha,tik_alpha, Fik
  real(kind=8),dimension(2,2)                         ::  Fabric,HN,HT,LN,LT,Matrice, Moment
  real(kind=8),dimension(:),allocatable               ::  tab_alpha,tab_FN,tab_FT,tab_LN,tab_LT,tab_C,tab_nx

  real(kind=8)                             :: S1,S2
  real(kind=8),dimension(2,2)              :: localframe
  real(kind=8),dimension(2)                :: wr,wi
  integer                                  :: ierror,matz,lda


  Nsect_orientation = 18
  allocate(tab_alpha(Nsect_orientation))
  tab_alpha = 0
  sect = pi / real(Nsect_orientation,8)
  tab_alpha(1) = sect

  do i=2,Nsect_orientation
    tab_alpha(i)=tab_alpha(i-1)+sect
  end do

  if (allocated(tab_C)) deallocate(tab_C)
  allocate(tab_C(Nsect_orientation))

  if (allocated(tab_FN)) deallocate(tab_FN)
  allocate(tab_FN(Nsect_orientation))

  if (allocated(tab_FT)) deallocate(tab_FT)
  allocate(tab_FT(Nsect_orientation))

  if (allocated(tab_LN)) deallocate(tab_LN)
  allocate(tab_LN(Nsect_orientation))

  if (allocated(tab_LT)) deallocate(tab_LT)
  allocate(tab_LT(Nsect_orientation))

  if (allocated(tab_nx)) deallocate(tab_nx)
  allocate(tab_nx(Nsect_orientation))

  tab_C   =  0
  tab_FN  =  0
  tab_FT  =  0
  tab_LN  =  0
  tab_LT  =  0
  tab_nx  =  0

  mean_LN = 0
  mean_LT = 0
  mean_FN = 0
  mean_FT = 0

  cpt = 0

  do i=1,nb_ligneCONTACT_POLYG
    cd      = TAB_CONTACT_POLYG(i)%icdent
    an      = TAB_CONTACT_POLYG(i)%ianent

    if (TAB_CONTACT_POLYG(i)%nature == 'PLJCx') cycle
    if (abs(TAB_CONTACT_POLYG(i)%rn) .lt. 1.D-9  .and. abs(TAB_CONTACT_POLYG(i)%rt) .lt. 1.D-9) cycle

    nik(1) = TAB_CONTACT_POLYG(i)%n(1)
    nik(2) = TAB_CONTACT_POLYG(i)%n(2)
    tik(1) = TAB_CONTACT_POLYG(i)%t(1)
    tik(2) = TAB_CONTACT_POLYG(i)%t(2)
    Rtik   = TAB_CONTACT_POLYG(i)%rt
    Rnik   = TAB_CONTACT_POLYG(i)%rn

    Lik(1) = TAB_POLYG(cd)%center(1)-TAB_POLYG(an)%center(1)
    Lik(2) = TAB_POLYG(cd)%center(2)-TAB_POLYG(an)%center(2)

    ! Building the stress tensor
    Fik(1) = (Rnik*nik(1)+Rtik*tik(1))
    Fik(2) = (Rnik*nik(2)+Rtik*tik(2))
    
    Moment(1,1:2) = Fik(1)*Lik(1:2) + Moment(1,1:2)
    Moment(2,1:2) = Fik(2)*Lik(1:2) + Moment(2,1:2)
  
    Force_N = TAB_CONTACT_POLYG(i)%rn !TAB_CONTACT_POLYG(i)%F(1)*nik(1) + TAB_CONTACT_POLYG(i)%F(2)*nik(2)
    Force_T = TAB_CONTACT_POLYG(i)%rt !TAB_CONTACT_POLYG(i)%F(1)*tik(1) + TAB_CONTACT_POLYG(i)%F(2)*tik(2)

    Norme_LN = (Lik(1)*nik(1) + Lik(2)*nik(2))
    Norme_LT = (Lik(1)*tik(1) + Lik(2)*tik(2))

    mean_FN = mean_FN + Force_N
    mean_LN = mean_LN + Norme_LN
    mean_FT = mean_FT + Force_T
    mean_LT = mean_LT + Norme_LT

    cpt = cpt + 1.0
  end do

  ! Computing the stress tensor
  Moment = Moment / (height*width)
  
  ! Finding the orientation of the stress tensor
  theta_sigma = abs( 0.5* atan2( -2*Moment(1,2) , (Moment(1,1)-Moment(2,2)) ) )

  mean_FN = mean_FN/cpt
  mean_LN = mean_LN/cpt
  mean_LT = mean_LT/cpt
  mean_FT = mean_FT/cpt

  Force_N  = 0.
  Force_T  = 0.
  Norme_LN = 0.
  Norme_LT = 0.

  do j=1,nb_ligneCONTACT_POLYG
    cd      = TAB_CONTACT_POLYG(j)%icdent
    an      = TAB_CONTACT_POLYG(j)%ianent

    if (TAB_CONTACT_POLYG(j)%nature == 'PLJCx') cycle
    if (abs(TAB_CONTACT_POLYG(j)%rn) .lt. 1.D-9  .and. abs(TAB_CONTACT_POLYG(j)%rt) .lt. 1.D-9) cycle

    nik(1) = TAB_CONTACT_POLYG(j)%n(1)
    nik(2) = TAB_CONTACT_POLYG(j)%n(2)
    tik(1) = TAB_CONTACT_POLYG(j)%t(1)
    tik(2) = TAB_CONTACT_POLYG(j)%t(2)

    Lik(1) = TAB_POLYG(cd)%center(1)-TAB_POLYG(an)%center(1)
    Lik(2) = TAB_POLYG(cd)%center(2)-TAB_POLYG(an)%center(2)

    Force_N  = TAB_CONTACT_POLYG(j)%rn !TAB_CONTACT_POLYG(i)%F(1)*nik(1) + TAB_CONTACT_POLYG(i)%F(2)*nik(2)
    Force_T  = TAB_CONTACT_POLYG(j)%rt !TAB_CONTACT_POLYG(i)%F(1)*tik(1) + TAB_CONTACT_POLYG(i)%F(2)*tik(2)
    Norme_LN = ( Lik(1)*nik(1) + Lik(2)*nik(2) )
    Norme_LT = ( Lik(1)*tik(1) + Lik(2)*tik(2) )

    if (Norme_LT<0.0000001) Norme_LT = 0
    if ( nik(2)<0 ) then
      nik(1) = -nik(1)
      nik(2) = -nik(2)
    end If

    alpha = acos(nik(1))

    do i=1,Nsect_orientation
    if (i==1) then
      if (alpha<tab_alpha(i)) then
        tab_C(i)     = tab_C(i)  + 1.0
        tab_FN(i)    = tab_FN(i) + Force_N
        tab_FT(i)    = tab_FT(i) + Force_T
        tab_LN(i)    = tab_LN(i) + Norme_LN
        tab_LT(i)    = tab_LT(i) + Norme_LT
        tab_nx(i)    = tab_nx(i) + alpha
        exit
      end if
      else if ((i>1).and.(i<Nsect_orientation)) then
      if ((alpha>=tab_alpha(i-1)).and.(alpha<tab_alpha(i))) then
        tab_C(i)     = tab_C(i)  + 1.0
        tab_FN(i)    = tab_FN(i) + Force_N
        tab_FT(i)    = tab_FT(i) + Force_T
        tab_LN(i)    = tab_LN(i) + Norme_LN
        tab_LT(i)    = tab_LT(i) + Norme_LT
        tab_nx(i)    = tab_nx(i) + alpha
        exit
      end if
      else if (i==Nsect_orientation) then
        if ((alpha<=tab_alpha(i)).and.(alpha>=tab_alpha(i-1))) then
          tab_C(i)     = tab_C(i)  + 1.0
          tab_FN(i)    = tab_FN(i) + Force_N
          tab_FT(i)    = tab_FT(i) + Force_T
          tab_LN(i)    = tab_LN(i) + Norme_LN
          tab_LT(i)    = tab_LT(i) + Norme_LT
          tab_nx(i)    = tab_nx(i) + alpha
        end if
      end if
    end do
  end do

  do i=1,Nsect_orientation
    if (tab_C(i)==0) cycle
    tab_FN(i) = tab_FN(i)/tab_C(i)
    tab_FT(i) = tab_FT(i)/tab_C(i)
    tab_LN(i) = tab_LN(i)/tab_C(i)
    tab_LT(i) = tab_LT(i)/tab_C(i)
    tab_nx(i) = tab_nx(i)/tab_C(i)
  end do

  HN(:,:) = 0
  HT(:,:) = 0
  LN(:,:) = 0
  LT(:,:) = 0

  do i=1,Nsect_orientation
    nik_alpha(1) = cos(tab_nx(i))
    nik_alpha(2) = sin(tab_nx(i))
    tik_alpha(1) =  nik_alpha(2)
    tik_alpha(2) = -nik_alpha(1)

    !######
    HN(1,1:2) = HN(1,1:2) + tab_FN(i)*nik_alpha(1)*nik_alpha(1:2)
    HN(2,1:2) = HN(2,1:2) + tab_FN(i)*nik_alpha(2)*nik_alpha(1:2)
    !######
    HT(1,1:2) = HT(1,1:2) + tab_FT(i)*nik_alpha(1)*tik_alpha(1:2)
    HT(2,1:2) = HT(2,1:2) + tab_FT(i)*nik_alpha(2)*tik_alpha(1:2)
    !######
    LN(1,1:2) = LN(1,1:2) + tab_LN(i)*nik_alpha(1)*nik_alpha(1:2)
    LN(2,1:2) = LN(2,1:2) + tab_LN(i)*nik_alpha(2)*nik_alpha(1:2)
    !######
    LT(1,1:2) = LT(1,1:2) + tab_LT(i)*nik_alpha(1)*tik_alpha(1:2)
    LT(2,1:2) = LT(2,1:2) + tab_LT(i)*nik_alpha(2)*tik_alpha(1:2)
    !######
  end do

  HN = HN / ( mean_FN*real(Nsect_orientation,8) )
  HT = HT / ( mean_FN*real(Nsect_orientation,8) )
  LN = LN / ( mean_LN*real(Nsect_orientation,8) )
  LT = LT / ( mean_LN*real(Nsect_orientation,8) )

  cpt = 0.D0
  Fabric = 0.D0

  do i=1,nb_ligneCONTACT_POLYG
    cd = TAB_CONTACT_POLYG(i)%icdent
    an = TAB_CONTACT_POLYG(i)%ianent

    if (TAB_CONTACT_POLYG(i)%nature == 'PLJCx') cycle
    if (abs(TAB_CONTACT_POLYG(i)%rn) .lt. 1.D-9  .and. abs(TAB_CONTACT_POLYG(i)%rt) .lt. 1.D-9) cycle

    nik(1)  = TAB_CONTACT_POLYG(i)%n(1)
    nik(2)  = TAB_CONTACT_POLYG(i)%n(2)

    Fabric(1,1:2) = nik(1)*nik(1:2) + Fabric(1,1:2)
    Fabric(2,1:2) = nik(2)*nik(1:2) + Fabric(2,1:2)

    cpt = cpt + 1
  end do

  Fabric = Fabric / cpt

  !!!!! CONTACT !!!!!
  Matrice = Fabric
  ac = ( 2*sqrt(  (Matrice(1,1)-Matrice(2,2))**2 + 4*Matrice(1,2)**2))/(Matrice(1,1)+Matrice(2,2))
  theta_c = abs( 0.5* atan2( -2*Matrice(1,2) , (Matrice(1,1)-Matrice(2,2)) ) ) ! ATTENTION IL FAUT UTILISER LA FONCTION TAN2 DE FORTRAN !
  ac_signe = ac*cos(2*(theta_c-theta_sigma))

  !!!!! FORCE NORMALE !!!!!
  Matrice = HN
  afn = ( 2*sqrt(  (Matrice(1,1)-Matrice(2,2))**2 + 4*Matrice(1,2)**2   ) )/(Matrice(1,1)+Matrice(2,2))
  theta_fn = abs( 0.5* atan2( -2*Matrice(1,2) , (Matrice(1,1)-Matrice(2,2)) ) ) ! ATTENTION IL FAUT UTILISER LA FONCTION TAN2 DE FORTRAN !
  afn_signe = afn*cos(2*(theta_fn-theta_sigma))

  !!!!! FORCE TANGENTE !!!!!
  Matrice = HT
  if ( (Matrice(1,1)==0.0) .and. (Matrice(2,2) == 0.0) ) Then
    aft = 0
    aft_signe = 0
    theta_ft = 0
  else
    aft = ( 2*sqrt(  (Matrice(1,1)-Matrice(2,2))**2 + 4*Matrice(1,2)**2   ) )/(HN(1,1)+HN(2,2))
    theta_ft = abs( 0.5* atan2( -2*Matrice(1,2) , (Matrice(1,1)-Matrice(2,2)) ) ) ! ATTENTION IL FAUT UTILISER LA FONCTION TAN2 DE FORTRAN !
    aft_signe = aft*cos(2*(theta_ft-theta_sigma))
    if (aft .lt. 1D-8) then
      theta_ft=0.
      aft = 0.
      aft_signe =0.
    end if
  end if

  !!!!! BRANCHE NORMALE !!!!!
  Matrice = LN
  aln = ( 2*sqrt(  (Matrice(1,1)-Matrice(2,2))**2 + 4*Matrice(1,2)**2   ) )/(Matrice(1,1)+Matrice(2,2))
  theta_ln = abs( 0.5* atan2( -2*Matrice(1,2) , (Matrice(1,1)-Matrice(2,2)) ) ) ! ATTENTION IL FAUT UTILISER LA FONCTION TAN2 DE FORTRAN !
  aln_signe = aln*cos(2*(theta_ln-theta_sigma))

  !!!!! BRANCHE TANGENTE !!!!!
  Matrice = LT
  alt = ( 2*sqrt(  (Matrice(1,1)-Matrice(2,2))**2 + 4*Matrice(1,2)**2   ) )/(LN(1,1)+LN(2,2))
  Theta_lt = abs( 0.5* atan2( -2*Matrice(1,2) , (Matrice(1,1)-Matrice(2,2)) ) ) ! ATTENTION IL FAUT UTILISER LA FONCTION TAN2 DE FORTRAN !
  alt_signe = alt*cos(2*(theta_lt-theta_sigma))

  double_product = 0.5*(ac*aln + ac*afn + aln*afn + alt*aft)
  double_product_sign = 0.5*(ac_signe*aln_signe + ac_signe*afn_signe + aln_signe*afn_signe + alt_signe*aft_signe)

  !write(117,'(30(1X,D12.5))')   temps,epsilon1,epsilon_q,&
  !ac,aln,alt,afn,aft,&
  !0.5*( ac+afn+aft+aln+alt ), &
  !0.5*( ac+afn+aft+aln+alt ) / (1+double_product),&
  !ac_signe,aln_signe,alt_signe,afn_signe,aft_signe,&
  !0.5*( ac_signe+afn_signe+aft_signe+aln_signe+alt_signe ), &
  !0.5*( ac_signe+afn_signe+aft_signe+aln_signe+alt_signe )/ (1+double_product_sign), &
  !theta_c,theta_ln,Theta_lt,theta_fn,theta_ft,theta_sigma

  if (first_over_all) then
     write(117,*) '#   time     ', '      ac     ', '      aln    ', '      alt    ', '      afn    ', '      aft    ', &
                                   '     ac_s    ', '     aln_s   ', '     alt_s   ', '     afn_s   ', '     aft_s   ', &
                                   '     t_c     ', '     t_ln    ', '     t_lt    ', '     t_fn    ', '     t_ft    ', &
                                   '     t_si    ', '     FN    ', '     FT    ', '     LN    ', '     LT    '
  end if 
  write(117,'(21(1X,E12.5))')   time,ac,aln,alt,afn,aft,&
                                      ac_signe,aln_signe,alt_signe,afn_signe,aft_signe,&
                                      theta_c,theta_ln,theta_lt,theta_fn,theta_ft, &
                                      theta_sigma, mean_FN, mean_FT, mean_LN, mean_LT
  deallocate(tab_C)
  deallocate(tab_FN)
  deallocate(tab_FT)
  deallocate(tab_LN)
  deallocate(tab_LT)
  deallocate(tab_nx)

  print*, 'Signed anisotropies     ---> Ok!'

end subroutine signed_anisotropy


!==============================================================================
! Calcul des ansiotropies de contacts signees a partir des tenseurs H sur le repere L
!==============================================================================
subroutine signed_anisotropy_l

  implicit none

  integer                                             ::  i,j,cd,an,Nsect_orientation
  real(kind=8)                                        ::  ac,aln,alt,afn,aft,&
                                                          ac_signe,aln_signe,alt_signe,afn_signe,aft_signe,&
                                                          theta_c,theta_ln,theta_lt,theta_fn,theta_ft, theta_sigma, &
                                                          double_product,double_product_sign

  real(kind=8)                                        ::  Force_N,Force_T,Norme_LN,Norme_LT,alpha, Rtik, Rnik, &
                                                          sect,cpt,mean_FN,mean_LN,mean_LT,mean_FT

  real(kind=8),dimension(2)                           ::  nik,tik,Lik,nik_alpha,tik_alpha, Fik, nik_c,tik_c
  real(kind=8),dimension(2,2)                         ::  Fabric,HN,HT,LN,LT,Matrice, Moment
  real(kind=8),dimension(:),allocatable               ::  tab_alpha,tab_FN,tab_FT,tab_LN,tab_LT,tab_C,tab_nx

  real(kind=8)                             :: S1,S2
  real(kind=8),dimension(2,2)              :: localframe
  real(kind=8),dimension(2)                :: wr,wi
  integer                                  :: ierror,matz,lda


  Nsect_orientation = 18
  allocate(tab_alpha(Nsect_orientation))
  tab_alpha = 0
  sect = pi / real(Nsect_orientation,8)
  tab_alpha(1) = sect

  do i=2,Nsect_orientation
    tab_alpha(i)=tab_alpha(i-1)+sect
  end do

  if (allocated(tab_C)) deallocate(tab_C)
  allocate(tab_C(Nsect_orientation))

  if (allocated(tab_FN)) deallocate(tab_FN)
  allocate(tab_FN(Nsect_orientation))

  if (allocated(tab_FT)) deallocate(tab_FT)
  allocate(tab_FT(Nsect_orientation))

  if (allocated(tab_LN)) deallocate(tab_LN)
  allocate(tab_LN(Nsect_orientation))

  if (allocated(tab_LT)) deallocate(tab_LT)
  allocate(tab_LT(Nsect_orientation))

  if (allocated(tab_nx)) deallocate(tab_nx)
  allocate(tab_nx(Nsect_orientation))

  tab_C   =  0
  tab_FN  =  0
  tab_FT  =  0
  tab_LN  =  0
  tab_LT  =  0
  tab_nx  =  0

  mean_LN = 0
  mean_LT = 0
  mean_FN = 0
  mean_FT = 0

  cpt = 0

  do i=1,nb_ligneCONTACT_POLYG
    cd      = TAB_CONTACT_POLYG(i)%icdent
    an      = TAB_CONTACT_POLYG(i)%ianent

    if (TAB_CONTACT_POLYG(i)%nature == 'PLJCx') cycle
    if (abs(TAB_CONTACT_POLYG(i)%rn) .lt. 1.D-9  .and. abs(TAB_CONTACT_POLYG(i)%rt) .lt. 1.D-9) cycle

    nik_c(1) = TAB_CONTACT_POLYG(i)%n(1)
    nik_c(2) = TAB_CONTACT_POLYG(i)%n(2)
    tik_c(1) = TAB_CONTACT_POLYG(i)%t(1)
    tik_c(2) = TAB_CONTACT_POLYG(i)%t(2)
    Rtik     = TAB_CONTACT_POLYG(i)%rt
    Rnik     = TAB_CONTACT_POLYG(i)%rn

    Lik(1) = TAB_POLYG(cd)%center(1)-TAB_POLYG(an)%center(1)
    Lik(2) = TAB_POLYG(cd)%center(2)-TAB_POLYG(an)%center(2)

    ! The new frame on L
    nik(1) = Lik(1)/(Lik(1)**2 + Lik(2)**2)**0.5
    nik(2) = Lik(2)/(Lik(1)**2 + Lik(2)**2)**0.5
    tik(1) = nik(2)
    tik(2) = -nik(1)

    ! Building the stress tensor (on the contact frame)
    Fik(1) = (Rnik*nik_c(1)+Rtik*tik_c(1))
    Fik(2) = (Rnik*nik_c(2)+Rtik*tik_c(2))
    
    Moment(1,1:2) = Fik(1)*Lik(1:2) + Moment(1,1:2)
    Moment(2,1:2) = Fik(2)*Lik(1:2) + Moment(2,1:2)
  
    !Force_N = TAB_CONTACT_POLYG(i)%rn !TAB_CONTACT_POLYG(i)%F(1)*nik(1) + TAB_CONTACT_POLYG(i)%F(2)*nik(2)
    !Force_T = TAB_CONTACT_POLYG(i)%rt !TAB_CONTACT_POLYG(i)%F(1)*tik(1) + TAB_CONTACT_POLYG(i)%F(2)*tik(2)

    ! Forces on the L frame 
    Force_N = Fik(1)*nik(1) + Fik(2)*nik(2)
    Force_T = Fik(1)*tik(1) + Fik(2)*tik(2)

    Norme_LN = (Lik(1)*nik(1) + Lik(2)*nik(2))
    Norme_LT = (Lik(1)*tik(1) + Lik(2)*tik(2))

    mean_FN = mean_FN + Force_N
    mean_LN = mean_LN + Norme_LN
    mean_FT = mean_FT + Force_T
    mean_LT = mean_LT + Norme_LT

    cpt = cpt + 1.0
  end do

  ! Computing the stress tensor
  Moment = Moment / (height*width)
  
  ! Finding the orientation of the stress tensor
  theta_sigma = abs( 0.5* atan2( -2*Moment(1,2) , (Moment(1,1)-Moment(2,2)) ) )

  mean_FN = mean_FN/cpt
  mean_LN = mean_LN/cpt
  mean_LT = mean_LT/cpt
  mean_FT = mean_FT/cpt

  Force_N  = 0.
  Force_T  = 0.
  Norme_LN = 0.
  Norme_LT = 0.

  do j=1,nb_ligneCONTACT_POLYG
    cd      = TAB_CONTACT_POLYG(j)%icdent
    an      = TAB_CONTACT_POLYG(j)%ianent

    if (TAB_CONTACT_POLYG(j)%nature == 'PLJCx') cycle
    if (abs(TAB_CONTACT_POLYG(j)%rn) .lt. 1.D-9  .and. abs(TAB_CONTACT_POLYG(j)%rt) .lt. 1.D-9) cycle

    nik_c(1) = TAB_CONTACT_POLYG(j)%n(1)
    nik_c(2) = TAB_CONTACT_POLYG(j)%n(2)
    tik_c(1) = TAB_CONTACT_POLYG(j)%t(1)
    tik_c(2) = TAB_CONTACT_POLYG(j)%t(2)
    Rtik     = TAB_CONTACT_POLYG(j)%rt
    Rnik     = TAB_CONTACT_POLYG(j)%rn

    Lik(1) = TAB_POLYG(cd)%center(1)-TAB_POLYG(an)%center(1)
    Lik(2) = TAB_POLYG(cd)%center(2)-TAB_POLYG(an)%center(2)

    ! The new frame on L
    nik(1) = Lik(1)/(Lik(1)**2 + Lik(2)**2)**0.5
    nik(2) = Lik(2)/(Lik(1)**2 + Lik(2)**2)**0.5
    tik(1) = nik(2)
    tik(2) = -nik(1)

    ! Building the stress tensor (on the contact frame)
    Fik(1) = (Rnik*nik_c(1)+Rtik*tik_c(1))
    Fik(2) = (Rnik*nik_c(2)+Rtik*tik_c(2))

    !Force_N  = TAB_CONTACT_POLYG(j)%rn !TAB_CONTACT_POLYG(i)%F(1)*nik(1) + TAB_CONTACT_POLYG(i)%F(2)*nik(2)
    !Force_T  = TAB_CONTACT_POLYG(j)%rt !TAB_CONTACT_POLYG(i)%F(1)*tik(1) + TAB_CONTACT_POLYG(i)%F(2)*tik(2)

    ! Forces on the L frame 
    Force_N = Fik(1)*nik(1) + Fik(2)*nik(2)
    Force_T = Fik(1)*tik(1) + Fik(2)*tik(2)

    Norme_LN = ( Lik(1)*nik(1) + Lik(2)*nik(2) )
    Norme_LT = ( Lik(1)*tik(1) + Lik(2)*tik(2) )

    if (Norme_LT<0.0000001) Norme_LT = 0
    if ( nik(2)<0 ) then
      nik(1) = -nik(1)
      nik(2) = -nik(2)
    end If

    alpha = acos(nik(1))

    do i=1,Nsect_orientation
    if (i==1) then
      if (alpha<tab_alpha(i)) then
        tab_C(i)     = tab_C(i)  + 1.0
        tab_FN(i)    = tab_FN(i) + Force_N
        tab_FT(i)    = tab_FT(i) + Force_T
        tab_LN(i)    = tab_LN(i) + Norme_LN
        tab_LT(i)    = tab_LT(i) + Norme_LT
        tab_nx(i)    = tab_nx(i) + alpha
        exit
      end if
      else if ((i>1).and.(i<Nsect_orientation)) then
      if ((alpha>=tab_alpha(i-1)).and.(alpha<tab_alpha(i))) then
        tab_C(i)     = tab_C(i)  + 1.0
        tab_FN(i)    = tab_FN(i) + Force_N
        tab_FT(i)    = tab_FT(i) + Force_T
        tab_LN(i)    = tab_LN(i) + Norme_LN
        tab_LT(i)    = tab_LT(i) + Norme_LT
        tab_nx(i)    = tab_nx(i) + alpha
        exit
      end if
      else if (i==Nsect_orientation) then
        if ((alpha<=tab_alpha(i)).and.(alpha>=tab_alpha(i-1))) then
          tab_C(i)     = tab_C(i)  + 1.0
          tab_FN(i)    = tab_FN(i) + Force_N
          tab_FT(i)    = tab_FT(i) + Force_T
          tab_LN(i)    = tab_LN(i) + Norme_LN
          tab_LT(i)    = tab_LT(i) + Norme_LT
          tab_nx(i)    = tab_nx(i) + alpha
        end if
      end if
    end do
  end do

  do i=1,Nsect_orientation
    if (tab_C(i)==0) cycle
    tab_FN(i) = tab_FN(i)/tab_C(i)
    tab_FT(i) = tab_FT(i)/tab_C(i)
    tab_LN(i) = tab_LN(i)/tab_C(i)
    tab_LT(i) = tab_LT(i)/tab_C(i)
    tab_nx(i) = tab_nx(i)/tab_C(i)
  end do

  HN(:,:) = 0
  HT(:,:) = 0
  LN(:,:) = 0
  LT(:,:) = 0

  do i=1,Nsect_orientation
    nik_alpha(1) = cos(tab_nx(i))
    nik_alpha(2) = sin(tab_nx(i))
    tik_alpha(1) =  nik_alpha(2)
    tik_alpha(2) = -nik_alpha(1)

    !######
    HN(1,1:2) = HN(1,1:2) + tab_FN(i)*nik_alpha(1)*nik_alpha(1:2)
    HN(2,1:2) = HN(2,1:2) + tab_FN(i)*nik_alpha(2)*nik_alpha(1:2)
    !######
    HT(1,1:2) = HT(1,1:2) + tab_FT(i)*nik_alpha(1)*tik_alpha(1:2)
    HT(2,1:2) = HT(2,1:2) + tab_FT(i)*nik_alpha(2)*tik_alpha(1:2)
    !######
    LN(1,1:2) = LN(1,1:2) + tab_LN(i)*nik_alpha(1)*nik_alpha(1:2)
    LN(2,1:2) = LN(2,1:2) + tab_LN(i)*nik_alpha(2)*nik_alpha(1:2)
    !######
    LT(1,1:2) = LT(1,1:2) + tab_LT(i)*nik_alpha(1)*tik_alpha(1:2)
    LT(2,1:2) = LT(2,1:2) + tab_LT(i)*nik_alpha(2)*tik_alpha(1:2)
    !######
  end do

  HN = HN / ( mean_FN*real(Nsect_orientation,8) )
  HT = HT / ( mean_FN*real(Nsect_orientation,8) )
  LN = LN / ( mean_LN*real(Nsect_orientation,8) )
  LT = LT / ( mean_LN*real(Nsect_orientation,8) )

  cpt = 0.D0
  Fabric = 0.D0

  do i=1,nb_ligneCONTACT_POLYG
    cd = TAB_CONTACT_POLYG(i)%icdent
    an = TAB_CONTACT_POLYG(i)%ianent

    if (TAB_CONTACT_POLYG(i)%nature == 'PLJCx') cycle
    if (abs(TAB_CONTACT_POLYG(i)%rn) .lt. 1.D-9  .and. abs(TAB_CONTACT_POLYG(i)%rt) .lt. 1.D-9) cycle

    !nik_c(1) = TAB_CONTACT_POLYG(i)%n(1)
    !nik_c(2) = TAB_CONTACT_POLYG(i)%n(2)
    !tik_c(1) = TAB_CONTACT_POLYG(i)%t(1)
    !tik_c(2) = TAB_CONTACT_POLYG(i)%t(2)

    Lik(1) = TAB_POLYG(cd)%center(1)-TAB_POLYG(an)%center(1)
    Lik(2) = TAB_POLYG(cd)%center(2)-TAB_POLYG(an)%center(2)

    ! The new frame on L
    nik(1) = Lik(1)/(Lik(1)**2 + Lik(2)**2)**0.5
    nik(2) = Lik(2)/(Lik(1)**2 + Lik(2)**2)**0.5
    tik(1) = nik(2)
    tik(2) = -nik(1)

    !nik(1)  = TAB_CONTACT_POLYG(i)%n(1)
    !nik(2)  = TAB_CONTACT_POLYG(i)%n(2)

    Fabric(1,1:2) = nik(1)*nik(1:2) + Fabric(1,1:2)
    Fabric(2,1:2) = nik(2)*nik(1:2) + Fabric(2,1:2)

    cpt = cpt + 1
  end do

  Fabric = Fabric / cpt

  !!!!! CONTACT !!!!!
  Matrice = Fabric
  ac = ( 2*sqrt(  (Matrice(1,1)-Matrice(2,2))**2 + 4*Matrice(1,2)**2))/(Matrice(1,1)+Matrice(2,2))
  theta_c = abs( 0.5* atan2( -2*Matrice(1,2) , (Matrice(1,1)-Matrice(2,2)) ) ) ! ATTENTION IL FAUT UTILISER LA FONCTION TAN2 DE FORTRAN !
  ac_signe = ac*cos(2*(theta_c-theta_sigma))

  !!!!! FORCE NORMALE !!!!!
  Matrice = HN
  afn = ( 2*sqrt(  (Matrice(1,1)-Matrice(2,2))**2 + 4*Matrice(1,2)**2   ) )/(Matrice(1,1)+Matrice(2,2))
  theta_fn = abs( 0.5* atan2( -2*Matrice(1,2) , (Matrice(1,1)-Matrice(2,2)) ) ) ! ATTENTION IL FAUT UTILISER LA FONCTION TAN2 DE FORTRAN !
  afn_signe = afn*cos(2*(theta_fn-theta_sigma))

  !!!!! FORCE TANGENTE !!!!!
  Matrice = HT
  if ( (Matrice(1,1)==0.0) .and. (Matrice(2,2) == 0.0) ) Then
    aft = 0
    aft_signe = 0
    theta_ft = 0
  else
    aft = ( 2*sqrt(  (Matrice(1,1)-Matrice(2,2))**2 + 4*Matrice(1,2)**2   ) )/(HN(1,1)+HN(2,2))
    theta_ft = abs( 0.5* atan2( -2*Matrice(1,2) , (Matrice(1,1)-Matrice(2,2)) ) ) ! ATTENTION IL FAUT UTILISER LA FONCTION TAN2 DE FORTRAN !
    aft_signe = aft*cos(2*(theta_ft-theta_sigma))
    if (aft .lt. 1D-8) then
      theta_ft=0.
      aft = 0.
      aft_signe =0.
    end if
  end if

  !!!!! BRANCHE NORMALE !!!!!
  Matrice = LN
  aln = ( 2*sqrt(  (Matrice(1,1)-Matrice(2,2))**2 + 4*Matrice(1,2)**2   ) )/(Matrice(1,1)+Matrice(2,2))
  theta_ln = abs( 0.5* atan2( -2*Matrice(1,2) , (Matrice(1,1)-Matrice(2,2)) ) ) ! ATTENTION IL FAUT UTILISER LA FONCTION TAN2 DE FORTRAN !
  aln_signe = aln*cos(2*(theta_ln-theta_sigma))

  !!!!! BRANCHE TANGENTE !!!!!
  Matrice = LT
  alt = ( 2*sqrt(  (Matrice(1,1)-Matrice(2,2))**2 + 4*Matrice(1,2)**2   ) )/(LN(1,1)+LN(2,2))
  Theta_lt = abs( 0.5* atan2( -2*Matrice(1,2) , (Matrice(1,1)-Matrice(2,2)) ) ) ! ATTENTION IL FAUT UTILISER LA FONCTION TAN2 DE FORTRAN !
  alt_signe = alt*cos(2*(theta_lt-theta_sigma))

  double_product = 0.5*(ac*aln + ac*afn + aln*afn + alt*aft)
  double_product_sign = 0.5*(ac_signe*aln_signe + ac_signe*afn_signe + aln_signe*afn_signe + alt_signe*aft_signe)

  !write(117,'(30(1X,D12.5))')   temps,epsilon1,epsilon_q,&
  !ac,aln,alt,afn,aft,&
  !0.5*( ac+afn+aft+aln+alt ), &
  !0.5*( ac+afn+aft+aln+alt ) / (1+double_product),&
  !ac_signe,aln_signe,alt_signe,afn_signe,aft_signe,&
  !0.5*( ac_signe+afn_signe+aft_signe+aln_signe+alt_signe ), &
  !0.5*( ac_signe+afn_signe+aft_signe+aln_signe+alt_signe )/ (1+double_product_sign), &
  !theta_c,theta_ln,Theta_lt,theta_fn,theta_ft,theta_sigma

  if (first_over_all) then
     write(125,*) '#   time     ', '      ac     ', '      aln    ', '      alt    ', '      afn    ', '      aft    ', &
                                   '     ac_s    ', '     aln_s   ', '     alt_s   ', '     afn_s   ', '     aft_s   ', &
                                   '     t_c     ', '     t_ln    ', '     t_lt    ', '     t_fn    ', '     t_ft    ', &
                                   '     t_si    ', '     FN    ', '     FT    ', '     LN    ', '     LT    '
  end if 
  write(125,'(21(1X,E12.5))')   time, ac, aln,  alt,  afn,  aft,&
                                      ac_signe,aln_signe,alt_signe,afn_signe,aft_signe,&
                                      theta_c, theta_ln,theta_lt,theta_fn,theta_ft, & 
                                      theta_sigma, mean_FN, mean_FT, mean_LN, mean_LT

  deallocate(tab_C)
  deallocate(tab_FN)
  deallocate(tab_FT)
  deallocate(tab_LN)
  deallocate(tab_LT)
  deallocate(tab_nx)

  print*, 'Signed anisotropies L   ---> Ok!'

end subroutine signed_anisotropy_l

!==============================================================================
! Clean TRESCA from parasit aditional contacts
!==============================================================================
subroutine clean_tresca
  
  implicit none
  
  integer                            :: i, j, k, cd, an, l_cdver, curr_check, n_l_cell
  integer                            :: l_counter_plpl, l_counter_pljc
  integer,dimension(:),allocatable   :: n_c_pp
  
  ! Checks the contact list and changes the status of contacts from Mstck to stick to cells not belonging to the same grain
  curr_check = 0

  ! How many cells per grain?   
  if (fin_post) then
    open(unit=118,file='CELL_NUMBER.DAT',status='old')

    ! Allocate the array
    if (allocated(n_c_pp)) deallocate(n_c_pp)
    allocate(n_c_pp(n_r_grain))

    do i=1, n_r_grain
      read(118,*) n_l_cell
      n_c_pp(i) = n_l_cell
    end do

    ! Now look for the good neighbors
    ! For all the particles
    do i=1, n_r_grain
      do j=curr_check+1, curr_check+n_c_pp(i)
        do k=1, nb_ligneCONTACT
          !if (TAB_CONTACT(k)%icdent == 20611 .or. TAB_CONTACT(k)%ianent ==20611) then
          !  print*, TAB_CONTACT(k)%icdent, ' ', TAB_CONTACT(k)%ianent, ' ', TAB_CONTACT(k)%nature, ' ', &
          !  TAB_CONTACT(k)%status(1:1), ' ', TAB_CONTACT(k)%gap
          !end if
          ! Taking care of polyg polyg contacts only
          if (TAB_CONTACT(k)%nature /= 'PLPLx') cycle

          ! Only Checking cohesive contacts
          if (TAB_CONTACT(k)%status(1:1) /= 'W') cycle

          ! Index of particles in contact
          cd = TAB_CONTACT(k)%icdent
          an = TAB_CONTACT(k)%ianent
        
          if (cd == j) then
            ! Check if this is a false positive case
            if (an .lt. curr_check .or. an .gt. curr_check+n_c_pp(i)) then
              TAB_CONTACT(k)%status = 'noctc'
              TAB_CONTACT(k)%gap0   = 0.D0
              TAB_CONTACT(k)%cd_len = 0.D0
              cycle
            end if
          else if (an==j) then
            ! Check if this is a false positive case
            if (cd .lt. curr_check .or. cd .gt. curr_check+n_c_pp(i)) then
              TAB_CONTACT(k)%status = 'noctc'
              TAB_CONTACT(k)%gap0   = 0.D0
              TAB_CONTACT(k)%cd_len = 0.D0
              cycle
            end if
          end if

          ! There's the case that it is inside the particle but it is not adjacent!
          ! If they are too far away... we need to skip!!!
          ! By hand! CAREFULL. 10% of Wth!!!
          if (TAB_CONTACT(k)%gap > (0.2000000D-01/20)) then
            TAB_CONTACT(k)%status = 'noctc'
            TAB_CONTACT(k)%gap0   = 0.D0
            TAB_CONTACT(k)%cd_len = 0.D0
          end if
        end do
      end do
      curr_check = curr_check + n_c_pp(i)
    end do
    
    ! Now it writes again the Vloc_Rloc file!
    write(123,'(A5)' ) '     '
    write(123,'(A11)') '! Vloc_Rloc'
    write(123,'(A5)' ) '     '
    write(123,'(A48)') '$steps        0              time= 0.0000000D+00'
    write(123,'(A1)' ) ' '
    write(123,'(A72)') '!-----------------------------------------------------------------------'
    
    l_counter_plpl = 0
    l_counter_pljc = 0
    do i=1, nb_ligneCONTACT
      
      ! Index of particles in contact
      cd = TAB_CONTACT(i)%icdent
      an = TAB_CONTACT(i)%ianent
      l_cdver = TAB_CONTACT(i)%cdver
      
      if (TAB_CONTACT(i)%nature == 'PLPLx') then
        l_counter_plpl = l_counter_plpl + 1 
        write(123,'(A6,2X,A5,2X,I7)')'$icdan','PLPLx',i
        write(123,'(A90)') ' cdbdy  numbr  cdtac  numbr  vertx  numbr  behav  anbdy  numbr  antac  numbr  segmt  numbr  sttus  iadj'
        write(123,'(1X,A5,2X,I5,2X,A5,2X,I5,2X,A5,2X,I5,2X,A5,2X,A5,2X,I5,2X,A5,2X,I5,2X,A5,2X,I5,2X,A5,1X,I5)') &
        'RBDY2', cd,'POLYG', 1, 'CDVER', l_cdver, TAB_CONTACT(i)%law, 'RBDY2', an,'POLYG',1, 'ANSEG', TAB_CONTACT(i)%sect, &
        TAB_CONTACT(i)%status, an
      else if (TAB_CONTACT(i)%nature == 'PLJCx') then
        l_counter_pljc = l_counter_pljc + 1 
        write(123,'(A6,2X,A5,2X,I7)')'$icdan','PLJCx', (i - l_counter_plpl)
        write(123,'(A103)') &
        ' cdbdy  numbr  cdtac  numbr  vertx  numbr  behav  anbdy  numbr  antac  numbr                sttus iadj '
        write(123,'(1X,A5,2X,I5,2X,A5,2X,I5,2X,A5,2X,I5,2X,A5,2X,A5,2X,I5,2X,A5,2X,I5,16X,A5,1X,I5)')   &
        'RBDY2',cd,'POLYG', 1,'CDVER',l_cdver, TAB_CONTACT(i)%law, 'RBDY2', an,'JONCx', 1, TAB_CONTACT(i)%status, &
        TAB_CONTACT(i)%adj
      else if (TAB_CONTACT(i)%nature == 'PTPT2') then
        write(123,'(A6,2X,A5,2X,I7)')'$icdan','PTPT2', (i - l_counter_plpl - l_counter_pljc)
        WRITE(123,'(A76)')' cdbdy  numbr  cdtac  numbr  behav  anbdy  numbr  antac  numbr  sttus   iadj'      
        WRITE(123,'(1X,A5,2X,I5,2X,A5,2X,I5,2X,A5,2X,A5,2X,I5,2X,A5,2X,I5,2X,A5,2X,I5)')   &
        'RBDY2', cd,'PT2Dx', 2, TAB_CONTACT(i)%law, 'RBDY2', an,'PT2Dx', 2,  TAB_CONTACT(i)%status, TAB_CONTACT(i)%adj
      else
        print*, 'Queeeeeee!'
        stop
      end if
      
      write(123,104) 'rlt/H', TAB_CONTACT(i)%rt         ,'rln/H', TAB_CONTACT(i)%rn         ,'rls/H', 0.0000000D+00
      write(123,104) 'vlt =', TAB_CONTACT(i)%vt         ,'vln =', TAB_CONTACT(i)%vn         ,'vls =', 0.0000000D+00
      write(123,103) 'gapTT', TAB_CONTACT(i)%gap
      write(123,104) 'n(1)=', TAB_CONTACT(i)%n(1)       ,'n(2)=', TAB_CONTACT(i)%n(2)       ,'n(3)=',TAB_CONTACT(i)%n(3)
      write(123,104) 'coo1=', TAB_CONTACT(i)%coor_ctc(1),'coo2=', TAB_CONTACT(i)%coor_ctc(2),'coo3=',TAB_CONTACT(i)%coor_ctc(3)
      
      if (TAB_CONTACT(i)%nature == 'PLPLx') then
        write(123,'(A)') '#  gapTTbegin      cd_length       disp_t                                  '
        write(123,'(3(1X,D14.7))') TAB_CONTACT(i)%gap0, TAB_CONTACT(i)%cd_len, TAB_CONTACT(i)%tang_disp
      else if (TAB_CONTACT(i)%nature == 'PTPT2') then
        write(123,'(A)') '#gapREF                                                                    '
        write(123,'(1X,D14.7)') TAB_CONTACT(i)%gap0
      end if
      write(123,*) ' '
    end do
    
    104 FORMAT(1X,5X,2X,5X,2X,5X,2X,5X,3(2X,A5,D14.7))
    103 FORMAT(1X,5X,2X,5X,2X,5X,2X,5X,2X,A5,D14.7)
    
    print*, 'Cleaning Tresca          ---> Ok!'
  end if
end subroutine clean_tresca

!==============================================================================
! Computing the average shape ratio
!==============================================================================
subroutine avg_shape_ratio
  
  implicit none 
  
  integer                                   :: i, j, n_sides
  real*8                                    :: short_h, long_l, avg_ratio, min_l, max_l, leng_l
  real*8                                    :: avg_radius
  real*8, dimension(2)                      :: max_vec, test_vec, orth_vec_uni

  ! Initializing variables
  short_h = 0
  long_l  = 0
  avg_ratio = 0
  avg_radius = 0
  
  ! For all the particles
  do i=1, n_particles

    ! Storing the radius
    avg_radius = avg_radius + TAB_POLYG(i)%radius
    
    n_sides = TAB_POLYG(i)%nb_vertex
    min_l =-9999.
    max_l =-9999.

    ! Finding the long side
    do j=1, n_sides
      leng_l = (TAB_POLYG(i)%vertex_ref(1,j)**2 + TAB_POLYG(i)%vertex_ref(2,j)**2)**0.5
      if (leng_l > max_l) then
        max_l = leng_l
        max_vec(1) = TAB_POLYG(i)%vertex_ref(1,j)
        max_vec(2) = TAB_POLYG(i)%vertex_ref(2,j)
      end if
    end do

    orth_vec_uni(1) = -max_vec(2)
    orth_vec_uni(2) =  max_vec(1)

    orth_vec_uni = orth_vec_uni/((orth_vec_uni(1)**2 + orth_vec_uni(2)**2)**0.5)

    ! Finding the orthogonal projection of 
    do j=1, n_sides
      ! The projection on the orthogonal vector
      test_vec(1) = TAB_POLYG(i)%vertex_ref(1,j)
      test_vec(2) = TAB_POLYG(i)%vertex_ref(2,j)
      leng_l = abs(test_vec(1)*orth_vec_uni(1)+test_vec(2)*orth_vec_uni(2))

      if (leng_l > min_l) then
        min_l = leng_l
      end if
    end do

    short_h = short_h + min_l
    long_l  = long_l  + max_l

    avg_ratio = avg_ratio + min_l/max_l
  end do
  
  short_h = short_h/n_particles
  long_l  = long_l/n_particles
  avg_ratio = avg_ratio/n_particles
  avg_radius = avg_radius/n_particles
  
  if (first_over_all) then
     write(124,*) '#   time     ', '     <h>     ', '     <L>     ', '    <h/L>    ', '    <rad>    '
  end if
  
  write(124,'(5(1X,E12.5))') time, short_h, long_l, avg_ratio, avg_radius

  print*, 'Write Average aspect ratio   ---> Ok!'
  
end subroutine avg_shape_ratio

!==============================================================================
! Vitesse moyenne
!==============================================================================
  subroutine vitesse_moyenne
  implicit none
    real(kind=8)                             :: V,Vx,Vy,Vr,cpt
    integer                                  :: i

    V   = 0.D0
    Vx  = 0.D0
    Vy  = 0.D0
    Vr  = 0.D0
    cpt = 0.D0
    do i=1,n_particles
      V  = V  + sqrt( TAB_POLYG(i)%V(1)**2 + TAB_POLYG(i)%V(2)**2 )
      Vx = Vx +       TAB_POLYG(i)%V(1)
      Vy = Vy +       TAB_POLYG(i)%V(2)
      Vr = Vr +       TAB_POLYG(i)%V(3)
      cpt = cpt +1
    end do

    V  = V/cpt
    Vx = Vx/cpt
    Vy = Vy/cpt
    Vr = Vr/cpt
    
    !if (type_box==1) write(101,'(5(1X,D12.5))') time,V,Vx,Vy,Vr
    !if (type_box==0) write(101,'(6(1X,D12.5))') time,height,V,Vx,Vy,Vr
  end subroutine vitesse_moyenne


!==============================================================================
! Test contact
!==============================================================================
subroutine test_position_contact
  
  implicit none
  
  integer                            ::  i
  real(kind=8)                       ::  xmax,ymax
  real(kind=8)                       ::  xmin,ymin
  
  xmax = -1000000000.D0
  ymax = -1000000000.D0
  xmin = 1000000000.D0
  ymin = 1000000000.D0
  
  nb_part_visi = 0
  
  ! Attention
  ! Case with 4 walls
  if (n_walls == 4) then 
    do i=1,n_walls
      xmax = max(xmax,TAB_PLAN(i)%center(1)-TAB_PLAN(i)%ax2)
      ymax = max(ymax,TAB_PLAN(i)%center(2)-TAB_PLAN(i)%ax2)
      xmin = min(xmin,TAB_PLAN(i)%center(1)+TAB_PLAN(i)%ax2)
      ymin = min(ymin,TAB_PLAN(i)%center(2)+TAB_PLAN(i)%ax2)
    end do
  ! Upper and down walls only
  else if (n_walls == 2) then 
    do i=1, n_walls
      ymax = max(ymax,TAB_PLAN(i)%center(2)-TAB_PLAN(i)%ax2)
      ymin = min(ymin,TAB_PLAN(i)%center(2)+TAB_PLAN(i)%ax2)
    end do

    do i=1,n_particles
      xmax = max(xmax,TAB_POLYG(i)%center(1))
      xmin = min(xmin,TAB_POLYG(i)%center(1))
    end do
  end if  
  
  width = (xmax - xmin)
  height = (ymax - ymin)
  
end subroutine test_position_contact

!==============================================================================
!CONSTRUCTION D UN NOUVEAU BODIES DOF ET DRV_DOF...
!==============================================================================
subroutine construct_bodies
implicit none
real(kind=8)                      :: ymax,ymin,dmin,xmax,xmin
integer                           :: i,j,cpt=0,nbPart,nb_wallh,nb_wallb

ymax = 0.D0
xmax = 0.D0
xmin = 10000.D0
ymin = 10000.D0
dmin = 0.D0
do i=1,n_particles
  ymax = max( ymax, TAB_POLYG(i)%center(2) )
  ymin = min( ymin, TAB_POLYG(i)%center(2) )
  xmax = max( xmax, TAB_POLYG(i)%center(1) )
  xmin = min( xmin, TAB_POLYG(i)%center(1) )
  dmin = dmin + 2*TAB_POLYG(i)%radius
end do
dmin = dmin / real(n_particles,8)

write(1004,'(A72)')'! File BODIES                                                           '
write(1004,'(A72)')'!                                                                       '
write(1004,'(A72)')'! The symbol    $       preceeds a keyword used in scanning files.      '
write(1004,'(A72)')'!                                                                       ' 
write(1004,'(A72)')'! The symbol    bdyty   stands for  body type data.                     '
write(1004,'(A72)')'! These data are distributed according to some species.                 '
write(1004,'(A72)')'!                                                                       ' 
write(1004,'(A72)')'! the specy     blmty   stands for  bulk element type data ,            '
write(1004,'(A72)')'! i.e. part or total bulk geometric description,                        '
write(1004,'(A72)')'! and bulk behaviour laws;                                              '
write(1004,'(A72)')'!                                                                       '
write(1004,'(A72)')'! the specy     nodty   stands for  node type data ,                    '
write(1004,'(A72)')'! i.e. degrees of freedom data;                                         '
write(1004,'(A72)')'!                                                                       '
write(1004,'(A72)')'! the specy     tacty   stands for  contactor type data ;               '
write(1004,'(A72)')'!                                                                       '
write(1004,'(A72)')'! the keyword   $$$$$$  ends a body record.                             '
write(1004,'(A72)')'!                                                                       '
write(1004,'(A6)')'      '

do i=1,n_particles
  TAB_POLYG(i)%behav = 'plexx'
  TAB_POLYG(i)%radius = 0.D0
  if ( TAB_POLYG(i)%center(2) < 2*dmin ) TAB_POLYG(i)%behav='wallb'
  if ( TAB_POLYG(i)%center(2) > ymax - 3*dmin ) TAB_POLYG(i)%behav='wallh'
end do

cpt = 0
do i=1,n_particles
  if (TAB_POLYG(i)%behav /= 'plexx')  cycle
  
  cpt = cpt + 1
  write(1004,'(A6)') '$bdyty'
  write(1004,'(1X,A5,2X,I5)') 'RBDY2',cpt
  write(1004,'(A6)') '$blmty'
  write(1004,'(1X,A5,2X,I5,2X,A5,2X,A5,2(2X,A5,D14.7))') &
       'PLAIN',1,'behav',TAB_POLYG(i)%behav,'avrd=',TAB_POLYG(i)%radius,'gyrd=',TAB_POLYG(i)%radius/sqrt(2.D0)
  write(1004,'(A6)') '$nodty'
  write(1004,'(1X,A5,2X,I5,2X,5X,2X,5X,3(2X,A5,D14.7))') &
       'NO3xx',1,'coo1=',TAB_POLYG(i)%center(1),'coo2=',TAB_POLYG(i)%center(2),'coo3=',0.D0
  write(1004,'(A6)') '$tacty'
  write(1004,'(1X,A5,2X,I5,2X,A5,2X,A5,2X,A10,I5,6X,A5,D14.7)') &
       'POLYG',1,'color','BLEUx','nb_vertex=',TAB_POLYG(i)%nb_vertex,'byrd=',TAB_POLYG(i)%radius
!  do j=1,TAB_POLYG(i)%nb_vertex
!    write(1004,'(29X,A5,D14.7,2X,A5,D14.7)') &
!        'coo1=',TAB_POLYG(i)%vertex(1,j)-TAB_POLYG(i)%center(1), &
!         'coo2=',TAB_POLYG(i)%vertex(2,j)-TAB_POLYG(i)%center(2)
!  enddo
  write(1004,'(A6)')'$$$$$$'
  write(1004,'(A6)')'      '
enddo
nbPart = cpt
nb_wallh = 0
do i=1,n_particles
  if (TAB_POLYG(i)%behav /= 'wallh')  cycle
  
  cpt = cpt + 1
  nb_wallh = nb_wallh + 1
  write(1004,'(A6)') '$bdyty'
  write(1004,'(1X,A5,2X,I5)') 'RBDY2',cpt
  write(1004,'(A6)') '$blmty'
  write(1004,'(1X,A5,2X,I5,2X,A5,2X,A5,2(2X,A5,D14.7))') &
       'PLAIN',1,'behav',TAB_POLYG(i)%behav,'avrd=',TAB_POLYG(i)%radius,'gyrd=',TAB_POLYG(i)%radius/sqrt(2.D0)
  write(1004,'(A6)') '$nodty'
  write(1004,'(1X,A5,2X,I5,2X,5X,2X,5X,3(2X,A5,D14.7))') &
       'NO3xx',1,'coo1=',TAB_POLYG(i)%center(1),'coo2=',TAB_POLYG(i)%center(2),'coo3=',0.D0
  write(1004,'(A6)') '$tacty'
  write(1004,'(1X,A5,2X,I5,2X,A5,2X,A5,2X,A10,I5,6X,A5,D14.7)') &
       'POLYG',1,'color','REDxx','nb_vertex=',TAB_POLYG(i)%nb_vertex,'byrd=',TAB_POLYG(i)%radius
!  do j=1,TAB_POLYG(i)%nb_vertex
!    write(1004,'(29X,A5,D14.7,2X,A5,D14.7)') &
!         'coo1=',TAB_POLYG(i)%vertex(1,j)-TAB_POLYG(i)%center(1), &
!         'coo2=',TAB_POLYG(i)%vertex(2,j)-TAB_POLYG(i)%center(2)
!  enddo
  write(1004,'(A6)')'$$$$$$'
  write(1004,'(A6)')'      '
enddo
nb_wallb = 0
!do i=1,n_particles
!  if (TAB_POLYG(i)%behav /= 'wallb')  cycle
!  
!  cpt = cpt + 1
!  nb_wallb = nb_wallb + 1
!  write(1004,'(A6)') '$bdyty'
!  write(1004,'(1X,A5,2X,I5)') 'RBDY2',cpt
!  write(1004,'(A6)') '$blmty'
!  write(1004,'(1X,A5,2X,I5,2X,A5,2X,A5,2(2X,A5,D14.7))') &
!       'PLAIN',1,'behav',TAB_POLYG(i)%behav,'avrd=',TAB_POLYG(i)%radius,'gyrd=',TAB_POLYG(i)%radius/sqrt(2.D0)
!  write(1004,'(A6)') '$nodty'
!  write(1004,'(1X,A5,2X,I5,2X,5X,2X,5X,3(2X,A5,D14.7))') &
!       'NO3xx',1,'coo1=',TAB_POLYG(i)%center(1),'coo2=',TAB_POLYG(i)%center(2),'coo3=',TAB_POLYG(i)%center(3)
!  write(1004,'(A6)') '$tacty'
!  write(1004,'(1X,A5,2X,I5,2X,A5,2X,A5,2X,A10,I5,6X,A5,D14.7)') &
!       'POLYG',1,'color','REDxx','nb_vertex=',TAB_POLYG(i)%nb_vertex,'byrd=',TAB_POLYG(i)%radius
!  do j=1,TAB_POLYG(i)%nb_vertex
!    write(1004,'(29X,A5,D14.7,2X,A5,D14.7)') &
!         'coo1=',TAB_POLYG(i)%vertex(1,j)-TAB_POLYG(i)%center(1), &
!         'coo2=',TAB_POLYG(i)%vertex(2,j)-TAB_POLYG(i)%center(2)
!  enddo
!  write(1004,'(A6)')'$$$$$$'
!  write(1004,'(A6)')'      '
!enddo

cpt = cpt + 1
write(1004,'(A6)') '$bdyty'
write(1004,'(1X,A5,2X,I5)') 'RBDY2',cpt
write(1004,'(A6)') '$blmty'
write(1004,'(1X,A5,2X,I5,2X,A5,2X,A5,2(2X,A5,D14.7))') &
   'PLAIN',1,'behav','wallb','avrd=',0.D0,'gyrd=',0.D0
write(1004,'(A6)') '$nodty'
write(1004,'(1X,A5,2X,I5,2X,5X,2X,5X,3(2X,A5,D14.7))') &
   'NO3xx',1,'coo1=',0.4750000D+02,'coo2=',3*dmin/2,'coo3=',0.D0
write(1004,'(A6)') '$tacty'
do i=1,n_particles
  if (TAB_POLYG(i)%behav /= 'wallb')  cycle
  
  nb_wallb = nb_wallb + 1
  write(1004,'(1X,A5,2X,I5,2X,A5,2X,A5,2X,A10,I5,6X,A5,D14.7)') &
       'POLYG',1,'color','REDxx','nb_vertex=',TAB_POLYG(i)%nb_vertex,'byrd=',0.D0
!  do j=1,TAB_POLYG(i)%nb_vertex
!    write(1004,'(29X,A5,D14.7,2X,A5,D14.7)') &
!         'coo1=',TAB_POLYG(i)%vertex(1,j)-0.4750000D+02, &
!         'coo2=',TAB_POLYG(i)%vertex(2,j)-3*dmin/2
!  enddo
enddo

cpt = cpt+1
write(1004,'(A6)') '$bdyty'
write(1004,'(1X,A5,2X,I5)') 'RBDY2',cpt
write(1004,'(A6)') '$blmty'
write(1004,'(1X,A5,2X,I5,2X,A5,2X,A5,2(2X,A5,D14.7))') &
   'PLAIN',1,'behav','TDURx','avrd=',0.1200000D+01,'gyrd=',0.8485281D+00
write(1004,'(A6)') '$nodty'
write(1004,'(1X,A5,2X,I5,2X,5X,2X,5X,3(2X,A5,D14.7))') &
   'NO3xx',1,'coo1=',0.4750000D+02,'coo2=',-0.1200000D+01,'coo3=',0.0000000D+00
write(1004,'(A6)') '$tacty'
write(1004,'(1X,A5,2X,I5,2X,A5,2X,A5,2(2X,A5,D14.7))') &
   'JONCx',1,'color','WALLx','ax1=',0.5000000D+02,'ax2=',0.1200000D+01
write(1004,'(A6)') '$$$$$$'
write(1004,'(A6)') '      '

print*,'nbPart=', nbPart
print*,'nb_wallh=',nb_wallh
print*,'nb_wallb=',nb_wallb


write(1006,'(A6)')  '! DOF'
write(1006,'(A1)')  ' '
                   !12345678901234567890123456789012345678901234567
write(1006,'(A47)') '$steps      0                time= 0.0000000D+00'
write(1006,'(A1)')  ' '
write(1006,'(A72)')'!-----------------------------------------------------------------------'



do i=1,nb_wallh
  write(1005,'(A72)') '!-----------------------------------------------------------------------' 
  write(1005,'(A1)')  ' '
  write(1005,'(A6)')  '$bdyty'
  write(1005,'(1X,A5,2X,I5)') 'RBDY2',i+nbPart
  write(1005,'(A6)')  '$nodty'
  write(1005,'(1X,A5,2X,I5)') 'NO3xx',1
  write(1005,'(A103)')  &
         '$dofty        [CT......+......AMP..*..cos.(..OMEGA.*.time.+.PHI..)]...*...[RAMPI.....+.....RAMP.*.time]'

  write(1005,'(1X,A5,2X,I5,6(1X,D14.7))') 'vlocy',1,0.2,0.,0.,0.,1.,0.
  write(1005,'(1X,A5,2X,I5,6(1X,D14.7))') 'vlocy',2,0.,0.,0.,0.,1.,0.
  write(1005,'(1X,A5,2X,I5,6(1X,D14.7))') 'vlocy',3,0.,0.,0.,0.,1.,0.
  write(1005,'(A6)')'$$$$$$'
end do
!do i=1,nb_wallb
!  write(1005,'(A72)') '!-----------------------------------------------------------------------' 
!  write(1005,'(A1)')  ' '
!  write(1005,'(A6)')  '$bdyty'
!  write(1005,'(1X,A5,2X,I5)') 'RBDY2',i+nbPart+nb_wallh
!  write(1005,'(A6)')  '$nodty'
!  write(1005,'(1X,A5,2X,I5)') 'NO3xx',1
!  write(1005,'(A103)')  &
!         '$dofty        [CT......+......AMP..*..cos.(..OMEGA.*.time.+.PHI..)]...*...[RAMPI.....+.....RAMP.*.time]'!
!
!  write(1005,'(1X,A5,2X,I5,6(1X,D14.7))') 'vlocy',1,0.,0.,0.,0.,1.,0.
!  write(1005,'(1X,A5,2X,I5,6(1X,D14.7))') 'force',2,0.,0.,0.,0.,1.,0.
!  write(1005,'(1X,A5,2X,I5,6(1X,D14.7))') 'vlocy',3,0.,0.,0.,0.,1.,0.
!  write(1005,'(A6)')'$$$$$$'
!end do

  write(1005,'(A72)') '!-----------------------------------------------------------------------' 
  write(1005,'(A1)')  ' '
  write(1005,'(A6)')  '$bdyty'
  write(1005,'(1X,A5,2X,I5)') 'RBDY2',nb_wallh+nbPart + 1
  write(1005,'(A6)')  '$nodty'
  write(1005,'(1X,A5,2X,I5)') 'NO3xx',1
  write(1005,'(A103)')  &
         '$dofty        [CT......+......AMP..*..cos.(..OMEGA.*.time.+.PHI..)]...*...[RAMPI.....+.....RAMP.*.time]'
  write(1005,'(1X,A5,2X,I5,6(1X,D14.7))') 'vlocy',1,0.,0.,0.,0.,1.,0.
  write(1005,'(1X,A5,2X,I5,6(1X,D14.7))') 'force',2,0.,0.,0.,0.,1.,0.
  write(1005,'(1X,A5,2X,I5,6(1X,D14.7))') 'vlocy',3,0.,0.,0.,0.,1.,0.
  write(1005,'(A6)')'$$$$$$'
  
  write(1005,'(A72)') '!-----------------------------------------------------------------------' 
  write(1005,'(A1)')  ' '
  write(1005,'(A6)')  '$bdyty'
  write(1005,'(1X,A5,2X,I5)') 'RBDY2',nb_wallh+nbPart + 2
  write(1005,'(A6)')  '$nodty'
  write(1005,'(1X,A5,2X,I5)') 'NO3xx',1
  write(1005,'(A103)')  &
         '$dofty        [CT......+......AMP..*..cos.(..OMEGA.*.time.+.PHI..)]...*...[RAMPI.....+.....RAMP.*.time]'

  write(1005,'(1X,A5,2X,I5,6(1X,D14.7))') 'vlocy',1,0.,0.,0.,0.,1.,0.
  write(1005,'(1X,A5,2X,I5,6(1X,D14.7))') 'force',2,10000000*(xmax-xmin),0.,0.,0.,1.,0.
  write(1005,'(1X,A5,2X,I5,6(1X,D14.7))') 'vlocy',3,0.,0.,0.,0.,1.,0.
  write(1005,'(A6)')'$$$$$$'

end subroutine construct_bodies

!==============================================================================
! Fermeture des fichiers
!==============================================================================

 subroutine close_all
 implicit none
 
   if (calcul_coordination              == 1)  close (100)
   if (calcul_vitesse_moyenne           == 1)  close (101)
   if (calcul_compacity                 == 1)  close (102)
   if (calcul_qoverp                    == 1)  close (103)
   if (calcul_anisotropy_contact        == 1)  close (104)
   if (calcul_anisotropy_force          == 1)  close (105)
   if (calcul_anisotropy_branch         == 1)  close (106)
!   if (calcul_persiste_contact          == 1)  close (106)
   if (walls_position                   == 1)  close (108)
   if (walls_f                          == 1)  close (109)
   if (c_n_ctc_probability              == 1)  close (111)
   if (c_fail_mode                      == 1)  close (113)
   if (c_sign_aniso                     == 1)  close (117)
   if (c_avg_shape_ratio                == 1)  close (124)
   if (c_sign_aniso_l                     == 1)  close (125)
  end subroutine close_all

!==============================================================================
!==============================================================================
!==============================================================================
!==============================================================================

!==============================================================================
!==============================================================================
!==============================================================================
!==============================================================================


!---------------------------------------------------
! subroutine pour calculer les vecteurs et valeurs
! propre d'une matrice N*N
!---------------------------------------------------
subroutine rg(lda, n, a, wr, wi, matz, z, ierror )
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




!===========================================================================
!That's all folks!
!===========================================================================
end program post2D_polyg
