!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!  POST-PROCESSOR !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!! DEFORMABLE BODIES !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program post2D_defo

implicit none

type  ::  T_CORPS_RIG
    integer                                  ::  num,nb_vertex,nb_neighbor,nb_ctc
    real(kind=8)                             ::  rad,surf, ax1, ax2
    real(kind=8),dimension(3)                ::  center,center_ref
    real(kind=8),dimension(3)                ::  veloc
    real(kind=8),allocatable,dimension(:,:)  ::  vertex,vertex_ref
    character(len=5)                         ::  behav
    character(len=6)                         ::  name
end type T_CORPS_RIG

type :: T_MESH_CONECT
  integer                                        ::  type_mesh
  integer,dimension(4)                           ::  connectivity
  character(len=5)                               ::  behav
  real*8                                         ::  det_j, surf_element
  real*8, dimension(2,2)                         ::  jacobi, inv_jacobi, def_grad
  real*8, dimension(2,2)                         ::  cg_right, cg_left, pk_2nd
  real*8, dimension(2,2)                         ::  green_E, dif_disp, cauchy_s
  real*8, dimension(2,2)                         ::  strain_gvp
end type T_MESH_CONECT

type :: T_NODE_INFO
  integer                                        ::  common_to
  real, dimension(2,2)                           ::  cg_right, cg_left, pk_2nd, def_grad
  real, dimension(2,2)                           ::  green_E, cauchy_s
end type T_NODE_INFO

type  ::  T_CORPS_DEFO
  integer                                        ::  id
  integer                                        ::  nb_mesh,nb_node,nb_tacty
  integer                                        ::  nctc
  real(kind=8)                                   ::  surface
  real(kind=8),dimension(3)                      ::  center
  real(kind=8),allocatable,dimension(:,:)        ::  veloc
  real(kind=8),allocatable,dimension(:,:)        ::  node,node_ref,tacty
  character(len=6)                               ::  name
  type(T_MESH_CONECT),allocatable,dimension(:)   ::  mesh_connect
  type(T_NODE_INFO), allocatable, dimension(:)   ::  node_info
end type T_CORPS_DEFO

type  ::  T_CONTACT
  integer                                        ::  icdent,ianent,sect,iantact,ictacty
  integer                                        ::  type, loc
  real(kind=8),dimension(2)                      ::  n, t
  real(kind=8),dimension(2)                      ::  F_tot
  real(kind=8),dimension(2)                      ::  coor_ctc
  real(kind=8)                                   ::  rn,rt,rs
  real(kind=8)                                   ::  vln,vlt,vs,gap, c_length
  character(len=5)                               ::  nature
  character(len=5)                               ::  status
  logical                                        ::  counted
end type T_CONTACT

type(T_CORPS_RIG),allocatable,dimension(:)       ::  TAB_WALL
type(T_CORPS_DEFO),allocatable,dimension(:)      ::  TAB_MESH
type(T_CONTACT),allocatable,dimension(:)         ::  TAB_CONTACT, TAB_CONTACT_BODY

!DEFINITIONS Vloc-Dof-Bodies
!================================================================
integer                            ::  n_walls=0, n_def_particles=0
real*8                             ::  fric_part = 0
integer                            ::  compteur_clout=0
logical                            ::  first_over_all=.true., fin_post=.false.
integer                            ::  time_step,nb_ligneCONTACT, nb_ligneCONTACT_BODY
real(kind=8)                       ::  time

!DEFINITIONS BOX
!================================================================
real(kind=8)                       ::  height, length, height_0, length_0
real(kind=8)                       ::  epsilon1, epsilon2, epsilon_q

!COMMANDS
!================================================================
character(len=30)                  ::  command
integer                            ::  c_compacity=0, c_wall_pos = 0, c_wall_f = 0
integer                            ::  c_coordination=0, c_list_forces_b=0, c_draw=0, &
                                       c_check_force=0, c_list_ctc_length=0, c_list_radius=0, &
                                       c_contact_anisotropy=0, c_force_anisotropy=0, &
                                       c_branch_anisotropy=0, c_ctc_distribution=0, c_force_distribution=0, &
                                       c_branch_distribution=0, c_qoverp=0, c_mesh_stresses=0, &
                                       c_qp_walls_nfric=0, c_strain_grain=0, c_green_lag=0, &
                                       c_avg_c_dist=0

! CONSTANTS
!================================================================
real(kind=8),parameter             :: pi=3.141593D0

! Variables DAC
integer                            ::  N_CLALp,N_CLJCx
real*8                             ::  f_b_scale, poisson_ratio, E_young

!================================================================
! Reading commands
!================================================================
print*,'-------------------------------------------'
print*,'!                                         !'
print*,'!           Post processing 2D            !'
print*,'!        from Emilien Azema code          !'
print*,'!                                         !'
print*,'-------------------------------------------'
print*,'Reading input file'
open(unit=1,file='./POST_INPUT.DAT',status='old')
do
  read(1,'(A30)') command
  print*,command

  if (command=='INFO SAMPLE                  :') then
    read(1,*) n_walls
    read(1,*) n_def_particles
    read(1,*) fric_part
    cycle
  end if

  if (command=='COMPACITY                    :') then
    c_compacity=1
    ! Uses port 100
    open (unit=100,file='./POSTPRO/COMPACITY.DAT',status='replace')
    cycle
  end if

  if (command=='WALL POSITION                :') then
    c_wall_pos=1
    ! Uses port 101
    open (unit=101,file='./POSTPRO/WALLSPOS.DAT',status='replace')
    cycle
  end if

  if (command=='WALL FORCE                   :') then
    c_wall_f=1
    ! Uses port 102
    open (unit=102,file='./POSTPRO/WALLFORCE.DAT',status='replace')
    cycle
  end if

  if (command=='COORDINATION                 :') then
    c_coordination=1
    ! Uses port 103
    open (unit=103,file='./POSTPRO/COORDINATION.DAT',status='replace')
    cycle
  end if

  if (command=='LIST FORCES ALONG BRANCH     :') then
    c_list_forces_b = 1
    ! Uses port 104
    cycle
  end if

  if (command=='LIST CONTACT LENGTH          :') then
    c_list_ctc_length = 1
    ! Uses port 105
    cycle
  end if

  if (command=='LIST RADIUS                  :') then
    c_list_radius = 1
    ! Uses port 106
    ! It is done over one single time step
    cycle
  end if

  if (command=='CONTACT ANISOTROPY           :') then
    c_contact_anisotropy=1
    ! Uses port 107
    open (unit=107,file='./POSTPRO/ANISOTROPY_CONTACT.DAT',status='replace')
    cycle
  end if

  if (command=='FORCE ANISOTROPY             :') then
    c_force_anisotropy=1
    ! Uses port 108
    open (unit=108,file='./POSTPRO/ANISOTROPY_FORCE.DAT',status='replace')
    cycle
  end if

  if (command=='BRANCH ANISOTROPY            :') then
    c_branch_anisotropy=1
    ! Uses port 109
    open (unit=109,file='./POSTPRO/ANISOTROPY_BRANCH.DAT',status='replace')
    cycle
  end if

  if (command=='CONTACTS DISTRIBUTION        :') then
    c_ctc_distribution=1
    ! Uses port 110
    cycle
  end if

  if (command=='FORCES DISTRIBUTION          :') then
    c_force_distribution=1
    ! Uses port 111
    cycle
  end if

  if (command=='BRANCH DISTRIBUTION          :') then
    c_branch_distribution=1
    ! Uses port 112
    cycle
  end if

  if (command=='Q/P                          :') then
    c_qoverp=1
    ! Uses port 113
    open (unit=113,file='./POSTPRO/Q_OVER_P.DAT',status='replace')
    cycle
  end if

  if (command=='QP WALLS NO FRIC             :') then
    c_qp_walls_nfric=1
    ! Uses port 114
    open (unit=114,file='./POSTPRO/QP_WALLS.DAT',status='replace')
    cycle
  end if

  if (command=='DRAW                         :') then
    c_draw = 1
    ! Uses port 121
    cycle
  end if

  if (command=='AVERAGE STRAIN GRAIN         :') then
    c_strain_grain = 1
    ! Uses port 122
    cycle
  end if

  if (command=='MESH STRESSES                :') then
    c_mesh_stresses=1
    read(1,*) poisson_ratio
    read(1,*) E_young
    c_green_lag = 1
    ! Green_Lab uses port 123
    cycle
  end if

  if (command=='AVERAGE INTERPART DIST       :') then
    c_avg_c_dist=1
    ! Uses port 124
    open (unit=124,file='./POSTPRO/AVG_CTC_DIST.DAT',status='replace')
    cycle
  end if

  if (command=='FORCE BARS SCALE             :') then
    read(1,*) f_b_scale
    cycle
  end if

  if (command=='CHECK FORCE                  :') then
    c_check_force = 1
    cycle
  end if

  if (command=='END                          :') exit
end do
close(1)

!==============================================================================
! Calling the subroutines
!==============================================================================
do

  call read_Vloc_dof_bodies
  call compute_size_box
  if (c_mesh_stresses==1) then
    call mesh_stresses
  end if

  if (c_compacity                  == 1)  call compacity
  if (c_wall_pos                   == 1)  call wall_pos
  if (c_wall_f                     == 1)  call wall_f
  if (c_coordination               == 1)  call coordination
  if (c_list_forces_b              == 1)  call list_forces_b
  if (c_check_force                == 1)  call check_force
  if (c_list_ctc_length            == 1)  call list_ctc_length
  if (c_contact_anisotropy         == 1)  call contact_anisotropy
  if (c_force_anisotropy           == 1)  call force_anisotropy
  if (c_branch_anisotropy          == 1)  call branch_anisotropy
  if (c_ctc_distribution           == 1)  call ctc_distribution
  if (c_force_distribution         == 1)  call force_distribution
  if (c_branch_distribution        == 1)  call branch_distribution
  if (c_qoverp                     == 1)  call qoverp
  if (c_qp_walls_nfric             == 1)  call qp_walls_nf
  if (c_strain_grain               == 1)  call strain_grain
  if (c_green_lag                  == 1)  call green_lagrange
  if (c_avg_c_dist                 == 1)  call avg_c_dist

  ! Over a single time step
  if (first_over_all) then
    if (c_list_radius                == 1)  call list_radius
  end if

  if (c_draw                       == 1)  call draw

  if (first_over_all) first_over_all=.false.
  if (fin_post) then
    call close_all
    exit
  end if
  ! Some void spaces for readability
  print*, ' '
  print*, ' '
end do
!==============================================================================

contains

!==============================================================================
! Reading historic files
!==============================================================================
subroutine read_Vloc_dof_bodies

  implicit none

  integer                       ::  i,j,k,kk,err,num_part,nodea,nodeb,nodec,noded
  integer                       ::  icdent,ianent,icdtact,iantact, cd, an
  integer                       ::  l_min, l_pos, cpt, l_counted
  real(kind=8)                  ::  center1,center2, f_norm
  real(kind=8)                  ::  vertex_coor1,vertex_coor2
  real(kind=8)                  ::  coor1,coor2,coor3,surf_mesh
  real(kind=8)                  ::  Vx,Vy,Vr,gap
  real(kind=8)                  ::  radius,ax1,ax2
  real(kind=8)                  ::  n1,n2,n3
  real(kind=8)                  ::  rn,rt,rs
  real(kind=8)                  ::  vln,vlt,vls
  real(kind=8)                  ::  coor_ctc1,coor_ctc2,coor_ctc3
  real(kind=8)                  ::  f_xx, f_xy, f_yx, f_yy
  real(kind=8)                  ::  sig_xx, sig_xy, sig_yx, sig_yy
  real(kind=8),dimension(3)     ::  g_gravity
  real(kind=8),dimension(2)     ::  X1,X2
  character(len=22)             ::  clout_DOF, clout_GPV
  character(len=28)             ::  clout_Vloc
  character(len=6)              ::  text
  character(len=13)             ::  text2
  character(len=5)              ::  behav,statut
  type(T_CONTACT),dimension(1)  ::  curr_obj

  ! Variables for building a convex hull (It is used for the contact length computation)
  integer                       ::  n_cur_points, c_real_ctc_p, c_loc_1, c_loc_2
  integer, dimension(:), allocatable  ::  c_points_order, check_stable
  real*8                        ::  c_ctc_length
  real*8,dimension(:,:), allocatable  ::  array_temp, array_temp_2
  logical                       ::  c_repeated


  compteur_clout = compteur_clout+1
  clout_DOF  =  './OUTBOX/DOF.OUT.    '
  clout_Vloc =  './OUTBOX/Vloc_Rloc.OUT.    '
  clout_GPV  =  './OUTBOX/GPV.OUT.    '

  ! Naming files
  if (compteur_clout<10) then
    write(clout_DOF(18:19),'(I1)')   compteur_clout
    write(clout_Vloc(24:25),'(I1)')  compteur_clout
    write(clout_GPV(18:19),'(I1)')  compteur_clout
  else if ( (compteur_clout>=10) .and. (compteur_clout<100) ) then
    write(clout_DOF(18:20),'(I2)')   compteur_clout
    write(clout_Vloc(24:26),'(I2)')  compteur_clout
    write(clout_GPV(18:20),'(I2)')  compteur_clout
  else if ( (compteur_clout>=100) .and. (compteur_clout<1000) ) then
    write(clout_DOF(18:21),'(I3)')   compteur_clout
    write(clout_Vloc(24:27),'(I3)')  compteur_clout
    write(clout_GPV(18:21),'(I3)')  compteur_clout
  else if ( (compteur_clout>=1000) .and. (compteur_clout<10000) ) then
    write(clout_DOF(18:22),'(I4)')   compteur_clout
    write(clout_Vloc(24:28),'(I4)')  compteur_clout
    write(clout_GPV(18:22),'(I4)')  compteur_clout
  end if

  ! Reading the bodies file
  if (first_over_all) then
    open(unit=2,file='./OUTBOX/BODIES.OUT',status='old')

    ! Reading the meshes for deformable particles
    allocate(TAB_MESH(n_def_particles))
    TAB_MESH(:)%nb_node  = 0
    TAB_MESH(:)%nb_mesh  = 0
    TAB_MESH(:)%nb_tacty = 0

    i = 0
    do
      read(2,'(A6)') text
      if (text == '$bdyty') then
        read(2,'(A6,1X,i6)') text, num_part
        if (text == ' MAILx') then
          i = i + 1
          TAB_MESH(num_part)%id = num_part
          read(2,'(A6)') text
          if (text == '$blmty') then
            do
              read(2,'(A6)') text
              if ((text==' Q4xxx').or.(text==' T3xxx')) TAB_MESH(num_part)%nb_mesh = &
                                                      TAB_MESH(num_part)%nb_mesh + 1
              if (text=='$nodty') exit
            end do
            do
              read(2,'(A6)') text
              if (text==' NO2xx') TAB_MESH(num_part)%nb_node = &
                                  TAB_MESH(num_part)%nb_node + 1
              if (text=='$tacty') exit
            end do
            do
              read(2,'(A6)') text
              if ((text==' ALpxx').or.(text=='+ALpxx')) TAB_MESH(num_part)%nb_tacty = &
                                                      TAB_MESH(num_part)%nb_tacty + 1
              if (text=='$$$$$$') exit
            end do
          end if
        end if
        if (i==n_def_particles) exit
      end if
    end do

    do i=1,n_def_particles
      allocate(TAB_MESH(i)%mesh_connect(TAB_MESH(i)%nb_mesh))
      allocate(TAB_MESH(i)%tacty(TAB_MESH(i)%nb_tacty,2))
      allocate(TAB_MESH(i)%node_ref(TAB_MESH(i)%nb_node,3))
      ! Initalizing with 3 components for vtk 
      allocate(TAB_MESH(i)%veloc(TAB_MESH(i)%nb_node,3))
      allocate(TAB_MESH(i)%node(TAB_MESH(i)%nb_node,3))
    end do

    rewind(2)

    i = 0
    do
      read(2,'(A6)') text
      if (text == '$bdyty') then
        read(2,'(A6,1X,i6)') text,num_part
        if (text == ' MAILx') then
          i = i + 1
          read(2,'(A6)') text
          if (text == '$blmty') then
            do j=1,TAB_MESH(num_part)%nb_mesh
              read(2,'(A6)') text
              if (text==' Q4xxx') TAB_MESH(num_part)%mesh_connect(j)%type_mesh = 4
              if (text==' T3xxx') TAB_MESH(num_part)%mesh_connect(j)%type_mesh = 3
              read(2,*)
            end do
            read(2,*)
            do j=1,TAB_MESH(num_part)%nb_node
              read(2,'(34X,D14.7,7X,D14.7)') vertex_coor1,vertex_coor2
              TAB_MESH(num_part)%node_ref(j,1) = vertex_coor1
              TAB_MESH(num_part)%node_ref(j,2) = vertex_coor2
              TAB_MESH(num_part)%node_ref(j,3) = 0.0
            end do
            read(2,*)
            do j=1,TAB_MESH(num_part)%nb_tacty
              read(2,'(A6,28X,i5,7X,i5)') text,nodea,nodeb
              if ((text==' CLxxx').or.(text=='+CLxxx')) cycle
              TAB_MESH(num_part)%tacty(j,1) = nodea
              TAB_MESH(num_part)%tacty(j,2) = nodeb
            end do
          end if
        end if
      end if
      if (i == n_def_particles) exit
    end do

    rewind(2)

    i = 0
    do
      read(2,'(A6)') text
      if (text == '$bdyty') then
        read(2,'(A6,1X,i6)') text,num_part
        if (text == ' MAILx') then
          i = i +1
          read(2,'(A6)') text
          if (text == '$blmty') then
            do j=1,TAB_MESH(num_part)%nb_mesh
              TAB_MESH(num_part)%mesh_connect(j)%connectivity(:)=0
              if (TAB_MESH(num_part)%mesh_connect(j)%type_mesh==4) then
                read(2,'(20X,4(i7))') nodea,nodeb,nodec,noded
                read(2,'(A36)') behav
                TAB_MESH(num_part)%mesh_connect(j)%connectivity(1)=nodea
                TAB_MESH(num_part)%mesh_connect(j)%connectivity(2)=nodeb
                TAB_MESH(num_part)%mesh_connect(j)%connectivity(3)=nodec
                TAB_MESH(num_part)%mesh_connect(j)%connectivity(4)=noded
                TAB_MESH(num_part)%mesh_connect(j)%behav=behav
              endif
              if (TAB_MESH(num_part)%mesh_connect(j)%type_mesh==3) then
                read(2,'(20X,3(i7))') nodea,nodeb,nodec
                read(2,'(A36)') behav
                TAB_MESH(num_part)%mesh_connect(j)%connectivity(1)=nodea
                TAB_MESH(num_part)%mesh_connect(j)%connectivity(2)=nodeb
                TAB_MESH(num_part)%mesh_connect(j)%connectivity(3)=nodec
                TAB_MESH(num_part)%mesh_connect(j)%behav=behav
              endif
            end do
          end if
        end if
      end if
      if (i == n_def_particles) exit
    end do

    ! Allocating space for walls
    allocate(TAB_WALL(n_walls))

    ! Preparing
    rewind(2)

    ! Reading the rigid walls
    i = 0
    do
      read(2,'(A6)') text
      if (text == '$bdyty') then
        read(2,'(A6)') text
        if (text == ' MAILx') cycle
        if (text == ' RBDY2') then
          i = i + 1
          read(2,*)
          read(2,'(22X,A5,7x,D14.7)') behav,radius
          TAB_WALL(i)%behav = behav
          TAB_WALL(i)%rad = radius
          read(2,*)
          read(2,'(29X, 2(5x,D14.7,2X))') center1,center2
          TAB_WALL(i)%center_ref(1) = center1
          TAB_WALL(i)%center_ref(2) = center2
          TAB_WALL(i)%center_ref(3) = 0.0
          read(2,*)
          read(2,'(29X, 2(5x,D14.7,2X))') ax1, ax2
          TAB_WALL(i)%ax1 = ax1
          TAB_WALL(i)%ax2 = ax2
          TAB_WALL(i)%nb_vertex=4

          allocate(TAB_WALL(i)%vertex_ref(3,4))
          allocate(TAB_WALL(i)%vertex(3,4))

          TAB_WALL(i)%vertex_ref(1,1) =  ax1
          TAB_WALL(i)%vertex_ref(2,1) =  ax2
          TAB_WALL(i)%vertex_ref(3,1) =  0.0

          TAB_WALL(i)%vertex_ref(1,2) =  - ax1
          TAB_WALL(i)%vertex_ref(2,2) =  ax2
          TAB_WALL(i)%vertex_ref(3,2) =  0.0

          TAB_WALL(i)%vertex_ref(1,3) =  - ax1
          TAB_WALL(i)%vertex_ref(2,3) =  - ax2
          TAB_WALL(i)%vertex_ref(3,3) =  0.0

          TAB_WALL(i)%vertex_ref(1,4) =  ax1
          TAB_WALL(i)%vertex_ref(2,4) =  - ax2
          TAB_WALL(i)%vertex_ref(3,4) =  0.0
        end if
      end if
      if (i==n_walls) exit
    end do
    close(2)
  end if

  ! Reading DOF files
  open(unit=3,file=clout_DOF,iostat=err,status='old')
  if (err/=0) then
    fin_post=.true.
  else
    i = 0
    k = 0
    kk= 0
    read(3,*)
    read(3,*)
    read(3,*)
    read(3,'(7X,i8,19X,D14.7)') time_step,time
    print*,'time_step:', time_step, time
    print*,'--->',clout_DOF
    do
      read(3,'(A6)') text
      if (text == '$bdyty') then
        i = i+1
        read(3,'(A6)') text
        if (text == ' MAILx') then
          kk = kk + 1
          read(3,*)
          read(3,*)
          read(3,*)
          ! Computing the center of gravity for this element
          g_gravity(:) = 0.D0
          do j=1,TAB_MESH(kk)%nb_node
            read(3,'(29X, 2(5x,D14.7,2X))') coor1,coor2
            read(3,'(29X, 2(5x,D14.7,2X))') Vx,Vy
            TAB_MESH(kk)%node(j,1) = TAB_MESH(kk)%node_ref(j,1) + coor1
            TAB_MESH(kk)%node(j,2) = TAB_MESH(kk)%node_ref(j,2) + coor2
            ! This is necessary for vtk
            TAB_MESH(kk)%node(j,3) = 0.0

            TAB_MESH(kk)%veloc(j,1)    = Vx
            TAB_MESH(kk)%veloc(j,2)    = Vy

            TAB_MESH(kk)%veloc(j,3) = 0.

            g_gravity(1) = g_gravity(1) + TAB_MESH(kk)%node(j,1)
            g_gravity(2) = g_gravity(2) + TAB_MESH(kk)%node(j,2)
            g_gravity(3) = 0.0
          end do
          g_gravity(:) = g_gravity(:) / real(TAB_MESH(kk)%nb_node,8)
          TAB_MESH(kk)%center(:) = g_gravity(:)
          ! Computing the surface for this element
          TAB_MESH(i)%surface = 0.
          do j=1,TAB_MESH(kk)%nb_mesh
            surf_mesh = 0.
            if (TAB_MESH(i)%mesh_connect(j)%type_mesh==4) then
              nodea = TAB_MESH(kk)%mesh_connect(j)%connectivity(1)
              nodeb = TAB_MESH(kk)%mesh_connect(j)%connectivity(2)
              nodec = TAB_MESH(kk)%mesh_connect(j)%connectivity(3)
              noded = TAB_MESH(kk)%mesh_connect(j)%connectivity(4)

              X1(:) = TAB_MESH(kk)%node(nodea,:)-TAB_MESH(kk)%node(nodeb,:)
              X2(:) = TAB_MESH(kk)%node(nodea,:)-TAB_MESH(kk)%node(nodec,:)

              surf_mesh = 0.5*abs( X1(1)*X2(2)-X1(2)*X2(1) )

              X1(:) = TAB_MESH(kk)%node(nodea,:)-TAB_MESH(kk)%node(nodec,:)
              X2(:) = TAB_MESH(kk)%node(nodea,:)-TAB_MESH(kk)%node(noded,:)

              surf_mesh = surf_mesh +0.5*abs( X1(1)*X2(2)-X1(2)*X2(1) )
            end if
            if (TAB_MESH(kk)%mesh_connect(j)%type_mesh==3) then
              nodea = TAB_MESH(kk)%mesh_connect(j)%connectivity(1)
              nodeb = TAB_MESH(kk)%mesh_connect(j)%connectivity(2)
              nodec = TAB_MESH(kk)%mesh_connect(j)%connectivity(3)

              X1(:) = TAB_MESH(kk)%node(nodea,:)-TAB_MESH(kk)%node(nodeb,:)
              X2(:) = TAB_MESH(kk)%node(nodea,:)-TAB_MESH(kk)%node(nodec,:)

              surf_mesh = 0.5*abs(X1(1)*X2(2)-X1(2)*X2(1))
            end if
            TAB_MESH(kk)%surface = TAB_MESH(kk)%surface + surf_mesh
          end do
        else if (text == ' RBDY2') then
          k = k + 1
          read(3,*)
          read(3,'(29X, 3(5x,D14.7,2X))') coor1,coor2,coor3
          read(3,'(29X, 3(5x,D14.7,2X))') Vx,Vy,Vr
          TAB_WALL(k)%veloc(1) = Vx
          TAB_WALL(k)%veloc(2) = Vy
          TAB_WALL(k)%veloc(3) = 0.0

          TAB_WALL(k)%center(1)=coor1+TAB_WALL(k)%center_ref(1)
          TAB_WALL(k)%center(2)=coor2+TAB_WALL(k)%center_ref(2)
          TAB_WALL(k)%center(3)=0.0

          do j=1,TAB_WALL(k)%nb_vertex
            TAB_WALL(k)%vertex(1,j) = TAB_WALL(k)%vertex_ref(1,j) * cos(coor3) - &
                                      TAB_WALL(k)%vertex_ref(2,j) * sin(coor3) + &
                                      TAB_WALL(k)%center(1)
            TAB_WALL(k)%vertex(2,j) = TAB_WALL(k)%vertex_ref(1,j) * sin(coor3) + &
                                      TAB_WALL(k)%vertex_ref(2,j) * cos(coor3) + &
                                      TAB_WALL(k)%center(2)
            TAB_WALL(k)%vertex(3,J) = 0.0
          end do
        end if
      end if
      if(i == n_def_particles+n_walls) exit
    end do
  end if

  ! Closing port
  close(3)

  N_CLALp = 0
  N_CLJCx = 0

  ! Reading the Vloc_Rloc files
  open(unit=4,file=clout_Vloc,iostat=err,status='old')
  if (err/=0) then
    print*, 'Error reading Vloc_Rloc file'
    stop
  else
    print*,'--->',clout_Vloc

    ! Reading the number of contacts
    do
      read(4,'(A13)',iostat=err) text2
      if (err/=0) then
        ! It should enter here only in case of end of file
        close(4)
        exit
      else
        if (text2 == '!-------------') cycle

        if (text2 == '$icdan  CLALp') then
          N_CLALp = N_CLALp +1
        elseif(text2 == '$icdan  CLJCx') then
          N_CLJCx = N_CLJCx +1
        end if
      end if
    end do
  end if

  open(unit=4,file=clout_Vloc,iostat=err,status='old')

  ! Allocating the tab for all contacts
  if (allocated(TAB_CONTACT)) deallocate(TAB_CONTACT)
  allocate(TAB_CONTACT(N_CLALp + N_CLALp))

  i = 0
  j = 0

  nb_ligneCONTACT=0

  ! Storing contact information.
  ! The idea is to store contacts sequentially.
  ! First, interparticle contacts. Then, particle walls contacts

  do
    read(4,'(A13)',iostat=err) text2
    if (err/=0) then
      print*, 'Error reading Vloc_Rloc p2'
      close(4)
      exit
    end if

    if (text2 == '!-------------') cycle

    if (text2 == '$icdan  CLALp') then
      read(4,*)
      read(4,'(8X,i5,9X,i5,16X,i5,9X,i5,9X,A5)') icdent,icdtact,ianent,iantact,statut
      read(4,'(29X, 3(5x,D14.7,2X))') rt,rn,rs
      read(4,'(29X, 3(5x,D14.7,2X))') vlt,vln,vls
      read(4,'(29X, (5x,D14.7,2X))')  gap
      read(4,'(29X, 3(5x,D14.7,2X))') n1,n2,n3
      read(4,'(29X, 3(5x,D14.7,2X))') coor_ctc1,coor_ctc2,coor_ctc3

      nb_ligneCONTACT=nb_ligneCONTACT+1
      i = i + 1
      TAB_CONTACT(i)%icdent=icdent
      TAB_CONTACT(i)%ianent=ianent
      TAB_CONTACT(i)%ictacty=icdtact
      TAB_CONTACT(i)%iantact=iantact
      TAB_CONTACT(i)%n(1)=n1
      TAB_CONTACT(i)%n(2)=n2
      TAB_CONTACT(i)%t(1)=n2
      TAB_CONTACT(i)%t(2)=-n1

      TAB_CONTACT(i)%coor_ctc(1)=coor_ctc1
      TAB_CONTACT(i)%coor_ctc(2)=coor_ctc2

      TAB_CONTACT(i)%rn=rn
      TAB_CONTACT(i)%rt=rt
      TAB_CONTACT(i)%rs=rs
      TAB_CONTACT(i)%vlt=vlt
      TAB_CONTACT(i)%gap=gap
      TAB_CONTACT(i)%nature='CLALp'
      TAB_CONTACT(i)%type  = 1
      TAB_CONTACT(i)%status  = statut
      TAB_CONTACT(i)%counted = .false.

    else if (text2 == '$icdan  CLJCx') then
      read(4,*)
      read(4,'(8X,i5,9X,i5,16X,i5,9X,i5,2X,A5)') icdent,icdtact,ianent,iantact,statut
      read(4,'(29X, 3(5x,D14.7,2X))') rt,rn,rs
      read(4,'(29X, 3(5x,D14.7,2X))') vlt,vln,vls
      read(4,*)
      read(4,'(29X, 3(5x,D14.7,2X))') n1,n2,n3
      read(4,'(29X, 3(5x,D14.7,2X))') coor_ctc1,coor_ctc2,coor_ctc3

      nb_ligneCONTACT=nb_ligneCONTACT+1
      j = j + 1
      TAB_CONTACT(N_CLALp+j)%icdent=icdent
      TAB_CONTACT(N_CLALp+j)%ianent=ianent
      TAB_CONTACT(N_CLALp+j)%n(1)=n1
      TAB_CONTACT(N_CLALp+j)%n(2)=n2

      TAB_CONTACT(N_CLALp+j)%t(1)=n2
      TAB_CONTACT(N_CLALp+j)%t(2)=-n1

      TAB_CONTACT(N_CLALp+j)%coor_ctc(1)=coor_ctc1
      TAB_CONTACT(N_CLALp+j)%coor_ctc(2)=coor_ctc2

      TAB_CONTACT(N_CLALp+j)%rn=rn
      TAB_CONTACT(N_CLALp+j)%rt=rt
      TAB_CONTACT(N_CLALp+j)%rs=rs
      TAB_CONTACT(N_CLALp+j)%vlt=vlt
      TAB_CONTACT(N_CLALp+j)%nature='CLJCx'
      TAB_CONTACT(N_CLALp+j)%status=statut
      TAB_CONTACT(N_CLALp+j)%type  = 1
      TAB_CONTACT(N_CLALp+j)%counted = .false.
    end if
    if(i+j==N_CLALp+N_CLJCx) exit
  end do

  ! Reading DOF files
  open(unit=5,file=clout_DOF,iostat=err,status='old')

  !print*, '0 state'
  !do i=1, nb_ligneCONTACT
  !  print*, TAB_CONTACT(i)%icdent, TAB_CONTACT(i)%ianent, TAB_CONTACT(i)%nature
  !end do

  ! Checking that antagonists have always bigger index than candidates
  ! Otherwise, the index are changed
  do i=1, nb_ligneCONTACT
    if (TAB_CONTACT(i)%nature == 'CLJCx') cycle
    cd = TAB_CONTACT(i)%icdent
    an = TAB_CONTACT(i)%ianent

    ! Changing if necessary
    if (cd .gt. an) then
      TAB_CONTACT(i)%icdent = an
      TAB_CONTACT(i)%ianent = cd
      ! Then the orientation of the local frame should be inversed
      TAB_CONTACT(i)%n(:) = -TAB_CONTACT(i)%n(:)
      TAB_CONTACT(i)%t(:) = -TAB_CONTACT(i)%t(:)
    end if
  end do

  ! The contacts table should be reordered now at two levels.
  ! The candidate level and the antagonist level

  !print*, 'Antes'
  !do i=1, nb_ligneCONTACT
  !  print*, TAB_CONTACT(i)%icdent, TAB_CONTACT(i)%ianent, TAB_CONTACT(i)%nature
  !end do

  do i=1, nb_ligneCONTACT
    if (TAB_CONTACT(i)%nature == 'CLJCx') cycle
    l_pos = i
    l_min = TAB_CONTACT(i)%icdent
    do j=i+1, nb_ligneCONTACT
      if (TAB_CONTACT(j)%nature == 'CLJCx') cycle
      if (TAB_CONTACT(j)%icdent .lt. l_min) then
        l_pos = j
        curr_obj(1) = TAB_CONTACT(l_pos)
        l_min = TAB_CONTACT(l_pos)%icdent
      end if
    end do
    if(l_pos .gt. i) then
      TAB_CONTACT(l_pos) = TAB_CONTACT(i)
      TAB_CONTACT(i) = curr_obj(1)
    end if
  end do

  !print*, 'Despues'
  !do i=1, nb_ligneCONTACT
  !  print*, TAB_CONTACT(i)%icdent, TAB_CONTACT(i)%ianent, TAB_CONTACT(i)%nature
  !end do

  ! Now ordering by antagonist
  do i=1, nb_ligneCONTACT
    if (TAB_CONTACT(i)%nature == 'CLJCx') cycle
    l_pos = i
    l_min = TAB_CONTACT(i)%ianent
    do j=i+1, nb_ligneCONTACT
      if (TAB_CONTACT(j)%nature == 'CLJCx') cycle
      if (TAB_CONTACT(j)%icdent .gt. TAB_CONTACT(i)%icdent) exit
      if (TAB_CONTACT(j)%ianent .lt. l_min) then
        l_pos = j
        curr_obj(1) = TAB_CONTACT(l_pos)
        l_min = TAB_CONTACT(l_pos)%ianent
      end if
    end do
    if (l_pos .gt. i) then
      TAB_CONTACT(l_pos) = TAB_CONTACT(i)
      TAB_CONTACT(i) = curr_obj(1)
    end if
  end do

  !print*, 'Despues_despues'
  !do i=1, nb_ligneCONTACT
  !  print*, TAB_CONTACT(i)%icdent, TAB_CONTACT(i)%ianent, TAB_CONTACT(i)%nature
  !end do

  ! Building the summarized array of contacts between particles
  ! for having the effective number of contacts (Emilien's style)

  l_counted= 0
  do i=1,nb_ligneCONTACT
    if (TAB_CONTACT(i)%counted) cycle
    cd = TAB_CONTACT(i)%icdent
    an = TAB_CONTACT(i)%ianent

    do j=i+1,nb_ligneCONTACT
      if (j-i > 20) exit     ! I suppose that there will be no more than 20 contacts per particle!
      if ((cd == TAB_CONTACT(j)%icdent).and.(an==TAB_CONTACT(j)%ianent).or.&
          (cd == TAB_CONTACT(j)%ianent).and.(an==TAB_CONTACT(j)%icdent) ) then

        TAB_CONTACT(j)%counted = .true.
        l_counted = l_counted + 1

        TAB_CONTACT(i)%type = TAB_CONTACT(i)%type + 1
        TAB_CONTACT(j)%type = TAB_CONTACT(i)%type
      end if
    end do
  end do

  ! Real number of contacts
  nb_ligneCONTACT_BODY = nb_ligneCONTACT - l_counted

  ! Allocating the array for the new contact table
  if (allocated(TAB_CONTACT_BODY)) deallocate(TAB_CONTACT_BODY)
  allocate(TAB_CONTACT_BODY(nb_ligneCONTACT_BODY))

  cpt = 0

  do i=1,nb_ligneCONTACT
    if (TAB_CONTACT(i)%counted) cycle
    cpt = cpt + 1

    ! Your data
    TAB_CONTACT_BODY(cpt)%loc      = i
    TAB_CONTACT_BODY(cpt)%icdent   = TAB_CONTACT(i)%icdent
    TAB_CONTACT_BODY(cpt)%ianent   = TAB_CONTACT(i)%ianent
    ! For deformable bodies it does not have a sense the tangential or normal directions
    !TAB_CONTACT_BODY(cpt)%n(:)     = TAB_CONTACT(i)%n(:)
    !TAB_CONTACT_BODY(cpt)%t(:)     = TAB_CONTACT(i)%t(:)
    !TAB_CONTACT_BODY(cpt)%rn       = TAB_CONTACT(i)%rn
    !TAB_CONTACT_BODY(cpt)%rt       = TAB_CONTACT(i)%rt
    TAB_CONTACT_BODY(cpt)%coor_ctc = TAB_CONTACT(i)%coor_ctc
    TAB_CONTACT_BODY(cpt)%vln      = TAB_CONTACT(i)%vln
    TAB_CONTACT_BODY(cpt)%vlt      = TAB_CONTACT(i)%vlt
    TAB_CONTACT_BODY(cpt)%type     = TAB_CONTACT(i)%type
    TAB_CONTACT_BODY(cpt)%nature   = TAB_CONTACT(i)%nature

    ! Building the total force vector
    TAB_CONTACT_BODY(cpt)%F_tot(1:2) = TAB_CONTACT(i)%rn*TAB_CONTACT(i)%n(1:2) + &
                                       TAB_CONTACT(i)%rt*TAB_CONTACT(i)%t(1:2)

    cd = TAB_CONTACT_BODY(cpt)%icdent
    an = TAB_CONTACT_BODY(cpt)%ianent

    ! Data of your friends
    do j=i+1, nb_ligneCONTACT
      if (j-i > 20) exit
      if ((cd == TAB_CONTACT(j)%icdent).and.(an==TAB_CONTACT(j)%ianent) .or. &
          (cd == TAB_CONTACT(j)%ianent).and.(an==TAB_CONTACT(j)%icdent) ) then

        TAB_CONTACT_BODY(cpt)%coor_ctc = TAB_CONTACT_BODY(cpt)%coor_ctc + TAB_CONTACT(j)%coor_ctc
        !TAB_CONTACT_BODY(cpt)%n(:)     = TAB_CONTACT_BODY(cpt)%n(:)     + TAB_CONTACT(j)%n(:)
        !TAB_CONTACT_BODY(cpt)%t(:)     = TAB_CONTACT_BODY(cpt)%t(:)     + TAB_CONTACT(j)%t(:)
        !TAB_CONTACT_BODY(cpt)%rn       = TAB_CONTACT_BODY(cpt)%rn       + TAB_CONTACT(j)%rn
        !TAB_CONTACT_BODY(cpt)%rt       = TAB_CONTACT_BODY(cpt)%rt       + TAB_CONTACT(j)%rt
        TAB_CONTACT_BODY(cpt)%vln      = TAB_CONTACT_BODY(cpt)%vln      + TAB_CONTACT(j)%vln
        TAB_CONTACT_BODY(cpt)%vlt      = TAB_CONTACT_BODY(cpt)%vlt      + TAB_CONTACT(j)%vlt


        ! For the total vector force
        TAB_CONTACT_BODY(cpt)%F_tot(1:2) = TAB_CONTACT_BODY(cpt)%F_tot(1:2) + &
                                           TAB_CONTACT(j)%rn*TAB_CONTACT(j)%n(1:2) + &
                                           TAB_CONTACT(j)%rt*TAB_CONTACT(j)%t(1:2)
      end if
    end do

    ! Finally
    !TAB_CONTACT_BODY(cpt)%n(:)     = TAB_CONTACT_BODY(cpt)%n(:)/TAB_CONTACT(i)%type
    !TAB_CONTACT_BODY(cpt)%t(:)     = TAB_CONTACT_BODY(cpt)%t(:)/TAB_CONTACT(i)%type

    TAB_CONTACT_BODY(cpt)%coor_ctc = TAB_CONTACT_BODY(cpt)%coor_ctc/TAB_CONTACT(i)%type
    TAB_CONTACT_BODY(cpt)%vln      = TAB_CONTACT_BODY(cpt)%vln/TAB_CONTACT(i)%type
    TAB_CONTACT_BODY(cpt)%vlt      = TAB_CONTACT_BODY(cpt)%vlt/TAB_CONTACT(i)%type

    ! Defining the status of the contact (Approximative)
    if (TAB_CONTACT_BODY(cpt)%vln**2 + TAB_CONTACT_BODY(cpt)%vlt**2 .gt. 0.000000001) then
      TAB_CONTACT_BODY(cpt)%status='slide'
    else
      TAB_CONTACT_BODY(cpt)%status='stick'
    end if
  end do

  !!! Now we can compute a couple of interesting quantities :)

  ! Computing the number of active contacts per particle
  TAB_MESH(:)%nctc = 0.

  do i=1,nb_ligneCONTACT_BODY
    icdent = TAB_CONTACT_BODY(i)%icdent
    ianent = TAB_CONTACT_BODY(i)%ianent

    ! Contacts between particles
    if (TAB_CONTACT_BODY(i)%nature == 'CLJCx') cycle
    ! Only active contacts
    f_norm = (TAB_CONTACT_BODY(i)%F_tot(1)**2 + TAB_CONTACT_BODY(i)%F_tot(1)**2)**0.5

    if (f_norm .le. 1.0E-9) cycle

    TAB_MESH(icdent)%nctc = TAB_MESH(icdent)%nctc + 1
    TAB_MESH(ianent)%nctc = TAB_MESH(ianent)%nctc + 1
  end do

  ! Computing the contact length
  i=1
  do
    if (allocated(array_temp)) deallocate(array_temp)
    allocate(array_temp(TAB_CONTACT(i)%type,2))

    ! Not arbitrary initialization
    array_temp(:,:) = -1

    if (allocated(c_points_order)) deallocate(c_points_order)
    allocate(c_points_order(TAB_CONTACT(i)%type))

    if (allocated(check_stable)) deallocate(check_stable)
    allocate(check_stable(TAB_CONTACT(i)%type))

    ! Aaaaaaaah fortrannnnnnn
    ! We need to check if there are repeated contact points for making the
    ! convex hull algorithm to work fine
    c_real_ctc_p = 0
    do j=1, TAB_CONTACT(i)%type
      if (j==1) then
        array_temp(j,1) = TAB_CONTACT(i+j-1)%coor_ctc(1)
        array_temp(j,2) = TAB_CONTACT(i+j-1)%coor_ctc(2)
        c_real_ctc_p = c_real_ctc_p +1
      else
        c_repeated = .false.
        do k=1, j-1
          !print*, j
          !print*, abs(TAB_CONTACT(i+j-1)%coor_ctc(1) - array_temp(k,1))
          !print*, abs(TAB_CONTACT(i+j-1)%coor_ctc(2) - array_temp(k,2))
          if (abs(TAB_CONTACT(i+j-1)%coor_ctc(1) - array_temp(k,1)) .lt. 1E-9 .and. &
              abs(TAB_CONTACT(i+j-1)%coor_ctc(2) - array_temp(k,2)) .lt. 1E-9) then
            c_repeated = .true.
          end if
        end do
        if (.not. c_repeated) then
          array_temp(j,1) = TAB_CONTACT(i+j-1)%coor_ctc(1)
          array_temp(j,2) = TAB_CONTACT(i+j-1)%coor_ctc(2)
          c_real_ctc_p = c_real_ctc_p + 1
        end if
      end if
    end do

    !print*, c_real_ctc_p
    c_ctc_length = 0

    ! aaaaah shit cases! with shit initializations!
    ! ....But computing at the end the contact length
    if (c_real_ctc_p .lt. TAB_CONTACT(i)%type) then
      if (allocated(array_temp_2)) deallocate(array_temp_2)
      allocate(array_temp_2(c_real_ctc_p,2))
      k=0
      do j=1, TAB_CONTACT(i)%type
        if(array_temp(j,1) .lt. 0) cycle
        k=k+1
        array_temp_2(k,1) = array_temp(j,1)
        array_temp_2(k,2) = array_temp(j,2)
      end do
      ! Building a convex hull with the set of contact points
      c_points_order(:) = 9999. !--> necessary
      call envelope(array_temp_2(:,1), array_temp_2(:,2), c_real_ctc_p, &
                c_points_order, n_cur_points, check_stable)
      !print*, array_temp_2
      !print*, c_points_order
      !print*, n_cur_points
      !print*, check_stable
      !stop
      do j=1, n_cur_points-1
        c_loc_1 = minloc(c_points_order(:)-j,1)
        c_points_order(c_loc_1) = 99.
        c_loc_2 = minloc(c_points_order(:)-(j+1),1)
        !print*, c_loc_1
        !print*, c_loc_2
        ! Distance between two points
        c_ctc_length = c_ctc_length + ((array_temp_2(c_loc_1,1) - array_temp_2(c_loc_2,1))**2 + &
                                       (array_temp_2(c_loc_1,2) - array_temp_2(c_loc_2,2))**2)**0.5
      end do
    else
      ! Building a convex hull with the set of contact points
      c_points_order(:) = 9999. !--> necessary
      call envelope(array_temp(:,1), array_temp(:,2), TAB_CONTACT(i)%type, &
                c_points_order, n_cur_points, check_stable)
      do j=1, n_cur_points-1
        c_loc_1 = minloc(c_points_order(:)-j,1)
        c_points_order(c_loc_1) = 99.
        c_loc_2 = minloc(c_points_order(:)-(j+1),1)

        ! Distance between two points
        c_ctc_length = c_ctc_length + ((array_temp(c_loc_1,1) - array_temp(c_loc_2,1))**2 + &
                                       (array_temp(c_loc_1,2) - array_temp(c_loc_2,2))**2)**0.5
      end do
    end if

    !print*, array_temp
    !print*, c_points_order
    !print*, n_cur_points
    !print*, check_stable

    ! Checking the index order
    TAB_CONTACT(i)%c_length = c_ctc_length
    !print*, c_ctc_length
    i=i+TAB_CONTACT(i)%type
    if (i .ge. nb_ligneCONTACT) exit
  end do

  close(5)

  ! Reading the GPV files
  open(unit=8,file=clout_GPV,iostat=err,status='old')
  if (err/=0) then
    print*, 'Error reading GPV file'
    close(8)
    stop
  else
    read(8,*)
    read(8,*)
    read(8,*)
    read(8,*)
    read(8,*)
    read(8,*)
    do i=1, n_def_particles
      read(8,*) text
      !print*, text
      !print*, i
      if (text == '$bdyty') then
        read(8,*)
        do j=1, TAB_MESH(i)%nb_mesh
          read(8,*)
          read(8,*)
          read(8,*)
          read(8,*)
          read(8,'(4(D14.7))') sig_xx, sig_xy, sig_yx, sig_yy
          read(8,'(4(D14.7))') f_xx, f_xy, f_yx, f_yy

          read(8,*)
          ! In the GVP is the deformation gradient that is written
          ! Now, I want to compute the Green Lagrange strain tensor E
          ! E = 1/2(C-I)
          ! With C = F^t * F
          ! Then
          TAB_MESH(i)%mesh_connect(j)%strain_gvp(1,1) = (f_xx*f_xx+f_yx*f_yx) - 1.
          TAB_MESH(i)%mesh_connect(j)%strain_gvp(1,2) = f_xx*f_xy+f_yx*f_yy
          TAB_MESH(i)%mesh_connect(j)%strain_gvp(2,1) = f_xy*f_xx+f_yy*f_yx
          TAB_MESH(i)%mesh_connect(j)%strain_gvp(2,2) = (f_xy*f_xy+f_yy*f_yy) - 1.

          TAB_MESH(i)%mesh_connect(j)%strain_gvp = 0.5*TAB_MESH(i)%mesh_connect(j)%strain_gvp
        end do
        read(8,*)
      else
        print*, 'Error while reading GPV'
        stop
      end if
    end do
  end if

  ! Closing port
  close(8)

end subroutine read_Vloc_dof_bodies


!==============================================================================
! Compute size box
!==============================================================================
subroutine compute_size_box

  implicit none

  integer                       ::  i,j
  real(kind=8)                  ::  xmax,ymax
  real(kind=8)                  ::  xmin,ymin

  xmax = -1000000000.E0
  ymax = -1000000000.E0
  xmin = 1000000000.E0
  ymin = 1000000000.E0

  do i=1,n_def_particles
    do j=1,TAB_MESH(i)%nb_node
      xmax = max(xmax,TAB_MESH(i)%node(j,1))
      ymax = max(ymax,TAB_MESH(i)%node(j,2))
      xmin = min(xmin,TAB_MESH(i)%node(j,1))
      ymin = min(ymin,TAB_MESH(i)%node(j,2))
    end do
  end do

  if (first_over_all) then
    height_0 = abs(ymax-ymin)
    length_0 = abs(xmax-xmin)
  end if

  height = abs(ymax-ymin)
  length = abs(xmax-xmin)

  epsilon1 = log(abs(height-height)/height_0 + 1)
  epsilon2 = log(abs(length-length_0)/length_0 + 1)

  epsilon_q = max(epsilon1,epsilon2)-min(epsilon1,epsilon2)
end subroutine compute_size_box

!==============================================================================
! Computing the solid fraction
!==============================================================================
subroutine compacity

  implicit none

  integer                                  :: i
  real(kind=8)                             :: surf_defor,surf_box

  surf_defor = 0.D0

  do i=1,n_def_particles
    surf_defor = surf_defor + TAB_MESH(i)%surface
  end do

  surf_box = length * height

  if (first_over_all) then
    write(100,*) '#   time     ', '   height    ', '    large    ',' V_grain/V_box '
  end if

  write(100,'(4(1X,E12.5))') time,height,length,surf_defor/surf_box

  print*, 'Writing Compacity               ---> Ok!'
end subroutine compacity
!==============================================================================

!==============================================================================
! Computing walls positions
!==============================================================================
subroutine wall_pos

  implicit none

  integer                         :: l_up(1), l_bottom(1), l_left(1), l_right(1)

  l_up     = maxloc((/TAB_WALL(1)%center(2),TAB_WALL(2)%center(2), TAB_WALL(3)%center(2),TAB_WALL(4)%center(2)/))
  l_bottom = minloc((/TAB_WALL(1)%center(2),TAB_WALL(2)%center(2), TAB_WALL(3)%center(2),TAB_WALL(4)%center(2)/))
  l_right  = maxloc((/TAB_WALL(1)%center(1),TAB_WALL(2)%center(1), TAB_WALL(3)%center(1),TAB_WALL(4)%center(1)/))
  l_left   = minloc((/TAB_WALL(1)%center(1),TAB_WALL(2)%center(1), TAB_WALL(3)%center(1),TAB_WALL(4)%center(1)/))

  if (first_over_all) then
    write(101,*) '#   time     ', '   down   X  ', '   down  Y   ', &
                                  '    up   X   ', '    up   Y   ', &
                                  '   left   X  ', '   left  Y   ', &
                                  '  right   X  ', '  right  Y   '
  endif

  write(101,'(9(1X,E12.5))') time, TAB_WALL(l_bottom)%center(1), TAB_WALL(l_bottom)%center(2)+TAB_WALL(l_bottom)%ax2 &
                                 , TAB_WALL(l_up)%center(1)    , TAB_WALL(l_up)%center(2)    -TAB_WALL(l_up)%ax2 &
                                 , TAB_WALL(l_left)%center(1)  + TAB_WALL(l_left)%ax2,        TAB_WALL(l_left)%center(2) &
                                 , TAB_WALL(l_right)%center(1) - TAB_WALL(l_right)%ax2,       TAB_WALL(l_right)%center(2)

  print*, 'Write Walls positions           ---> Ok!'

end subroutine wall_pos

!==============================================================================
! Computing forces on walls
!==============================================================================
subroutine wall_f

  implicit none

  integer                         :: i
  integer                         :: l_up(1), l_bottom(1), l_left(1), l_right(1)
  real*8,dimension(2)             :: f_up, f_bottom, f_left, f_right

  ! Vectors containing the forces.
  ! First component: rn, Second component: rt
  f_up(:) = 0
  f_bottom(:) = 0
  f_left(:) = 0
  f_right(:) = 0

  l_up     = maxloc((/TAB_WALL(1)%center(2),TAB_WALL(2)%center(2), TAB_WALL(3)%center(2),TAB_WALL(4)%center(2)/))
  l_bottom = minloc((/TAB_WALL(1)%center(2),TAB_WALL(2)%center(2), TAB_WALL(3)%center(2),TAB_WALL(4)%center(2)/))
  l_right  = maxloc((/TAB_WALL(1)%center(1),TAB_WALL(2)%center(1), TAB_WALL(3)%center(1),TAB_WALL(4)%center(1)/))
  l_left   = minloc((/TAB_WALL(1)%center(1),TAB_WALL(2)%center(1), TAB_WALL(3)%center(1),TAB_WALL(4)%center(1)/))


  ! For all contacts
  do i=1, nb_ligneCONTACT
    ! Only contacts with walls
    if (TAB_CONTACT(i)%nature == 'CLALp') cycle
    !Only active contacts
    if (TAB_CONTACT(i)%rn .lt. 1.0E-9) cycle

    if (TAB_CONTACT(i)%ianent==l_up(1)) then
      f_up(1) = f_up(1) + TAB_CONTACT(i)%rn
      f_up(2) = f_up(2) + TAB_CONTACT(i)%rt
    elseif(TAB_CONTACT(i)%ianent==l_bottom(1)) then
      f_bottom(1) = f_bottom(1) + TAB_CONTACT(i)%rn
      f_bottom(2) = f_bottom(2) + TAB_CONTACT(i)%rt
    elseif(TAB_CONTACT(i)%ianent==l_right(1)) then
      f_right(1) = f_right(1) + TAB_CONTACT(i)%rn
      f_right(2) = f_right(2) + TAB_CONTACT(i)%rt
    elseif(TAB_CONTACT(i)%ianent==l_left(1)) then
      f_left(1) = f_left(1) + TAB_CONTACT(i)%rn
      f_left(2) = f_left(2) + TAB_CONTACT(i)%rt
    else
      print*, 'Unknown contact'
      stop
    end if
  end do

  if (first_over_all) then
    write(102,*) '#   time     ', '   down   N  ', '   down  T   ', &
                                  '    up   N   ', '    up   T   ', &
                                  '   left   N  ', '   left  T   ', &
                                  '  right   N  ', '  right  T   '
  endif

  write(102,'(9(1X,E12.5))') time, f_bottom(1), f_bottom(2), f_up(1)   , f_up(2), &
                                   f_left(1)  , f_left(2)  , f_right(1), f_right(2)

  print*, 'Write Walls forces              ---> Ok!'

end subroutine wall_f

!==============================================================================
! Computing forces on walls
!==============================================================================
subroutine coordination

  implicit none

  integer                                  :: i, cd, an, float_part, n_l_particles
  real*8                                   :: z, zc, zp, l_f_norm
  integer, dimension(:), allocatable       :: valid_part

  ! Initalizing variables
  z   = 0
  zc  = 0
  zp  = 0
  float_part = 0
  n_l_particles = 0

  ! This part supposes that walls are always declared at the end of the bodies list
  if (allocated(valid_part)) deallocate(valid_part)
  allocate(valid_part(n_def_particles))

  valid_part(:) = 1

  do i=1, nb_ligneCONTACT_BODY
    cd = TAB_CONTACT_BODY(i)%icdent
    an = TAB_CONTACT_BODY(i)%ianent

    if(TAB_CONTACT_BODY(i)%nature == 'CLJCx') valid_part(cd) = 0
    if(TAB_MESH(cd)%nctc .lt. 2) valid_part(cd) = 0
    !if(TAB_MESH(an)%nctc .lt. 2) valid_part(an) = 0
  end do

  ! Computing the number of contacts
  do i=1, nb_ligneCONTACT_BODY
    ! Only between particles
    if (TAB_CONTACT_BODY(i)%nature == 'CLJCx') cycle
    cd = TAB_CONTACT_BODY(i)%icdent
    an = TAB_CONTACT_BODY(i)%ianent

    ! Very important! Removing contact between particles that are in contact with walls
    if (valid_part(cd) .lt. 1 .and. valid_part(an) .lt. 1) cycle

    ! Computing the number of contacts if they are active
    l_f_norm = (TAB_CONTACT_BODY(i)%F_tot(1)**2 + TAB_CONTACT_BODY(i)%F_tot(2)**2)**0.5

    if (l_f_norm .gt. 1e-9) then
      zc = zc + 1
    end if
  end do

  do i=1, n_def_particles
    if (valid_part(i) == 1) then
      zp = zp + 1
    end if
  end do

 !print*, n_l_particles
 !print*, float_part
 !print*, zp
 !print*, zc

  z = 2*zc/zp

  !print*, z

  if (first_over_all) then
    write(103,*) '#   time     ', '      z      '
  end if

  write(103,'(2(1X,E12.5))') time, z

  print*, 'Write Coordination              ---> Ok!'

end subroutine coordination

!==============================================================================
! List of normal forces
!==============================================================================
subroutine list_forces_b

  implicit none

  integer                                  :: i, cd, an
  real*8, dimension(3)                     :: Lik, F_b
  character(len=32)                        :: l_file
  logical                                  :: dir_c

  ! Cleaning or creating the folder if necessary
  if (first_over_all) then
    ! Asking if the file already exists
    inquire(file='./POSTPRO/LIST_FB', exist=dir_c)
    if(dir_c) then
      ! Cleaning
      call system('rm ./POSTPRO/LIST_FB/*')
    else
      ! Creating
      call system('mkdir ./POSTPRO/LIST_FB')
    end if
  end if

  ! The name of the file
  l_file = './POSTPRO/LIST_FB/F_Branch.     '

  ! Preparing the corresponding number of the file
  if (compteur_clout<10) then
    WRITE(l_file(28:29),'(I1)')   compteur_clout
  else if ( (compteur_clout>=10) .and. (compteur_clout<100) ) then
    WRITE(l_file(28:30),'(I2)')   compteur_clout
  else if ( (compteur_clout>=100).and. (compteur_clout<1000) ) then
    WRITE(l_file(28:31),'(I3)')   compteur_clout
  else if ( (compteur_clout>=1000).and. (compteur_clout<10000) ) then
    WRITE(l_file(28:32),'(I4)')   compteur_clout
  else
    print*, "Not enough file names :: List of normal forces"
    stop
  end if

  ! Opening the file
  open(unit=104,file=l_file,status='replace')

  write(104,*) '#   FN     '

  do i=1,nb_ligneCONTACT_BODY
    ! Contacts between particles
    if (TAB_CONTACT_BODY(i)%nature == 'CLJCx') cycle
    ! Active contacts
    if ((TAB_CONTACT_BODY(i)%F_tot(1)**2 + TAB_CONTACT_BODY(i)%F_tot(2)**2)**0.5 &
        .le. 1.0E-9) cycle

    cd = TAB_CONTACT_BODY(i)%icdent
    an = TAB_CONTACT_BODY(i)%ianent

    ! Computing the branch vector
    Lik(:) = TAB_MESH(cd)%center(:)-TAB_MESH(an)%center(:)

    ! Normalizing the branch vector
    Lik(:) = Lik(:)/((Lik(1)**2 + Lik(2)**2)**0.5)

    ! Projection of the force vector on the branch
    F_b(1) = TAB_CONTACT_BODY(i)%F_tot(1)*Lik(1)
    F_b(2) = TAB_CONTACT_BODY(i)%F_tot(2)*Lik(2)
    F_b(3) = 0

    ! Writing the norm of this force vector
    write(104, '(E12.5)') (F_b(1)**2 + F_b(2)**2)**0.5
  end do

  close(104)

  print*, 'Write List Forces Branch        ---> Ok!'

end subroutine list_forces_b

!==============================================================================
! List of contact lengths
!==============================================================================
subroutine list_ctc_length

  implicit none

  integer                                  :: i
  character(len=32)                        :: l_file
  logical                                  :: dir_c

  ! Cleaning or creating the folder if necessary
  if (first_over_all) then
    ! Asking if the file already exists
    inquire(file='./POSTPRO/LIST_CL', exist=dir_c)
    if(dir_c) then
      ! Cleaning
      call system('rm ./POSTPRO/LIST_CL/*')
    else
      ! Creating
      call system('mkdir ./POSTPRO/LIST_CL')
    end if
  end if

  ! The name of the file
  l_file = './POSTPRO/LIST_CL/C_LENGTH.     '

  ! Preparing the corresponding number of the file
  if (compteur_clout<10) then
    WRITE(l_file(28:29),'(I1)')   compteur_clout
  else if ( (compteur_clout>=10) .and. (compteur_clout<100) ) then
    WRITE(l_file(28:30),'(I2)')   compteur_clout
  else if ( (compteur_clout>=100).and. (compteur_clout<1000) ) then
    WRITE(l_file(28:31),'(I3)')   compteur_clout
  else if ( (compteur_clout>=1000).and. (compteur_clout<10000) ) then
    WRITE(l_file(28:32),'(I4)')   compteur_clout
  else
    print*, "Not enough file names :: List of contact length"
    stop
  end if

  ! Opening the file
  open(unit=105,file=l_file,status='replace')

  write(105,*) '#  Length   '

  do i=1, nb_ligneCONTACT_BODY
    ! Contacts between particles
    if (TAB_CONTACT_BODY(i)%nature == 'CLJCx') cycle
    ! Active contacts
    if ((TAB_CONTACT_BODY(i)%F_tot(1)**2 + TAB_CONTACT_BODY(i)%F_tot(2)**2)**0.5 &
        .le. 1.0E-9) cycle

    ! Writing the norm of this force vector
    write(105, '(E12.5)') TAB_CONTACT(TAB_CONTACT_BODY(i)%loc)%c_length
  end do

  close(105)

  print*, 'Write List Contact Length       ---> Ok!'

end subroutine list_ctc_length


!==============================================================================
! List of particle radius
!==============================================================================
subroutine list_radius

  implicit none

  integer                                  :: i
  character(len=22)                        :: l_file

  ! The name of the file
  l_file = './POSTPRO/P_RADIUS.DAT'

  ! Opening the file
  open(unit=106,file=l_file,status='replace')

  write(106,*) '#  Radius   '

  do i=1, n_def_particles
    ! Writing the equivalent radius for each particle
    write(106, '(E13.6)') (TAB_MESH(i)%surface/pi)**0.5
  end do

  close(106)

  print*, 'Write List Radius               ---> Ok!'

end subroutine list_radius


!==============================================================================
! Computing the contact anisotropy
!==============================================================================
subroutine contact_anisotropy

  implicit none

  integer                                  :: i,cd,an
  integer                                  :: ierror,matz,lda
  real*8                                   :: F1, F2, cpt, Rnik, ac
  real*8,dimension(2,2)                    :: Fabric
  real*8,dimension(2)                      :: nik, Lik
  real*8,dimension(2)                      :: wr,wi
  real*8,dimension(2,2)                    :: localframe

  ! Initializing the number of active contacts
  cpt = 0

  ! Initializing the fabric tensor
  Fabric(:,:) = 0.0

  ! Building the fabric tensor
  do i=1,nb_ligneCONTACT_BODY
    ! Checking it is a contact between two bodies
    if(TAB_CONTACT_BODY(i)%nature == 'CLJCx') cycle
    cd   = TAB_CONTACT_BODY(i)%icdent
    an   = TAB_CONTACT_BODY(i)%ianent

    Rnik = (TAB_CONTACT_BODY(i)%F_tot(1)**2+TAB_CONTACT_BODY(i)%F_tot(2)**2)**0.5

    ! Checking if it is an active contact
    if (Rnik .le. 1.E-9) cycle

    ! Computing the branch vector
    Lik(1) = TAB_MESH(cd)%center(1)-TAB_MESH(an)%center(1)
    Lik(2) = TAB_MESH(cd)%center(2)-TAB_MESH(an)%center(2)

    nik  = Lik/(Lik(1)**2+Lik(2)**2)**0.5

    Fabric(1,1:2) = Fabric(1,1:2) + nik(1)*nik(1:2)
    Fabric(2,1:2) = Fabric(2,1:2) + nik(2)*nik(1:2)

    cpt = cpt + 1
  end do

  ! Normalizing by the number of contacts
  Fabric = Fabric / cpt

  ! Finding the eigen values
  lda  = 2
  matz = 1

  call rg (lda, 2, Fabric, wr, wi, matz, localframe, ierror)

  F1 = max(wr(1),wr(2))
  F2 = min(wr(1),wr(2))

  ! Computing ac from the principal values of F

  ac = 2*(F1-F2)

  if (first_over_all) then
    write(107,*) '#   time     ', '     S1      ', '     S2      ', '     ac      '
  end if

  write(107,'(5(1X,E12.5))') time, F1, F2, ac

  print*, 'Write Anisotropy Conctacts      ---> Ok!'

end subroutine contact_anisotropy


!==============================================================================
! Computing the force anisotropy
!==============================================================================
subroutine force_anisotropy

  implicit none

  integer                                  :: i, cd, an
  integer                                  :: ierror,matz,lda
  real*8                                   :: X1n, X2n, X1f, X2f
  real*8                                   :: cpt, Rtik, Rnik, av_force, norm_FN, norm_FT
  real*8, dimension(2)                     :: wrn, win, wrf, wif
  real*8, dimension(2)                     :: nik,tik, dirFT
  real*8, dimension(2)                     :: Fik, FikT, Lik
  real*8, dimension(2,2)                   :: Fabric_N,Fabric_T, Fabric_F
  real*8, dimension(2,2)                   :: localframeN, localframeF

  ! Initializing variables
  Fabric_N(:,:) = 0.
  Fabric_T(:,:) = 0.
  Fabric_F(:,:) = 0.

  cpt = 0.

  av_force = 0.
  Fik(:)   = 0.
  FikT(:)  = 0.

  ! Building the fabric forces tensor
  do i=1,nb_ligneCONTACT_BODY
    if  (TAB_CONTACT_BODY(i)%nature == 'CLJCx') cycle
    cd   = TAB_CONTACT_BODY(i)%icdent
    an   = TAB_CONTACT_BODY(i)%ianent

    Rnik = (TAB_CONTACT_BODY(i)%F_tot(1)**2+TAB_CONTACT_BODY(i)%F_tot(2)**2)**0.5

    ! Checking if it is an active contact
    if (Rnik .le. 1.E-9) cycle

    ! Computing the branch vector
    Lik(1) = TAB_MESH(cd)%center(1)-TAB_MESH(an)%center(1)
    Lik(2) = TAB_MESH(cd)%center(2)-TAB_MESH(an)%center(2)

    nik  = Lik/(Lik(1)**2+Lik(2)**2)**0.5

    ! Building the tangential direction... perpendicular to nik
    tik(1) =  nik(2)
    tik(2) = -nik(1)

    ! Finding the tangential force
    Rtik = TAB_CONTACT_BODY(i)%F_tot(1)*tik(1) + TAB_CONTACT_BODY(i)%F_tot(2)*tik(2)

    !print*, "First"
    !print*, Rnik

    ! Active contacts
    cpt = cpt + 1

    ! The force vector
    Fik(:) = Rnik*nik(:) + Rtik*tik(:)

    ! Average normal force
    !av_force = av_force + Fik(1)*nik(1) + Fik(2)*nik(2) + Fik(3)*nik(3)
    ! Average total force
    av_force = av_force + sqrt(Fik(1)**2 + Fik(2)**2)

    ! The norm of the normal force...
    norm_FN = Fik(1)*nik(1) + Fik(2)*nik(2)

    !print*, norm_FN

    ! Building the tensor of normal forces
    Fabric_N(1,1:2) = Fabric_N(1,1:2) + norm_FN*nik(1)*nik(1:2)
    Fabric_N(2,1:2) = Fabric_N(2,1:2) + norm_FN*nik(2)*nik(1:2)

    ! The tangential force vector
    FikT(1) = Fik(1) - Rnik*nik(1)
    FikT(2) = Fik(2) - Rnik*nik(2)

    ! The norm of FikT
    norm_FT = sqrt(FikT(1)**2 + FikT(2)**2)

    ! The direction of FikT
    dirFT = FikT/norm_FT

    ! Building the tensor of tangential forces
    Fabric_T(1,1:2) = Fabric_T(1,1:2) + norm_FT*nik(1)*dirFT(1:2)
    Fabric_T(2,1:2) = Fabric_T(2,1:2) + norm_FT*nik(2)*dirFT(1:2)
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

  ! Finding the eigenvalues of the fabric normal forces tensor
  lda = 2
  matz = 1

  call rg(lda, 2, Fabric_N, wrn, win, matz, localframeN, ierror)

  X1n = max(wrn(1),wrn(2))
  X2n = min(wrn(1),wrn(2))

  ! Finding the eigen values of the full fabric forces tensor
  call rg(lda, 3, Fabric_F, wrf, wif, matz, localframeF, ierror)

  X1f = max(wrf(1),wrf(2))
  X2f = min(wrf(1),wrf(2))

  if (first_over_all) then
     write(108,*) '#   time    ', ' 2(X2-1)/(X2+3)n', ' 2(X3-1)/(X1+3)f '
  end if

  !print*, X1f,  X2f, X3f

  write(108,'(3(1X,E12.5))') time, 2*(X1n-X2n)/(X1n+X2n), 2*(X1f-X2f)/(X1f+X2f)

  print*, 'Write Anisotropy Force          ---> Ok!'
end subroutine force_anisotropy


!==============================================================================
! Computing the branch anisotropy
!==============================================================================
subroutine branch_anisotropy

  implicit none

  integer                                  :: i,cd,an
  integer                                  :: ierror,matz,lda
  real*8                                   :: X1n, X2n, X1f, X2f
  real*8                                   :: cpt, av_length, av_total, Lnik, Ltik, Rnik
  real*8, dimension(2)                     :: nik, tik, ikt
  real*8, dimension(2)                     :: wrn, win, wrf, wif
  real*8, dimension(2)                     :: Lik
  real*8, dimension(2,2)                   :: Fabric_N,Fabric_T, Fabric_F
  real*8, dimension(2,2)                   :: localframeN, localframeF


  ! Initializing variables
  Fabric_N(:,:) = 0.
  Fabric_T(:,:) = 0.
  Fabric_F(:,:) = 0.

  cpt = 0.

  av_length = 0.
  av_total = 0.
  Lik(:) = 0.

  ! Building the chi tensor with the length of branches
  do i=1,nb_ligneCONTACT_BODY
    if  (TAB_CONTACT_BODY(i)%nature == 'CLJCx') cycle
    cd   = TAB_CONTACT_BODY(i)%icdent
    an   = TAB_CONTACT_BODY(i)%ianent

    Rnik = (TAB_CONTACT_BODY(i)%F_tot(1)**2+TAB_CONTACT_BODY(i)%F_tot(2)**2)**0.5

    ! Checking if it is an active contact
    if (Rnik .le. 1.E-9) cycle

    ! Computing the branch vector
    Lik(1) = TAB_MESH(cd)%center(1)-TAB_MESH(an)%center(1)
    Lik(2) = TAB_MESH(cd)%center(2)-TAB_MESH(an)%center(2)

    nik  = Lik/(Lik(1)**2+Lik(2)**2)**0.5

    ! Building the tangential direction... perpendicular to nik
    tik(1) =  nik(2)
    tik(2) = -nik(1)

    ! Active contacts
    cpt = cpt + 1

    ! Average branch length in the normal contact direction
    !av_length = av_length + Lik(1)*nik(1) + Lik(2)*nik(2) + Lik(3)*nik(3)

    ! Average total length
    !av_total = av_total + sqrt(Lik(1)**2 + Lik(2)**2 + Lik(3)**2)

    ! The norm of the branch in the normal contact direction
    Lnik = Lik(1)*nik(1) + Lik(2)*nik(2)

    ! The branch in the tangential direction
    Ltik = Lik(1)*tik(1) + Lik(2)*tik(2)

    ! Building the chi_n
    Fabric_N(1,1:2) = Lnik*nik(1)*nik(1:2) + Fabric_N(1,1:2)
    Fabric_N(2,1:2) = Lnik*nik(2)*nik(1:2) + Fabric_N(2,1:2)

    ! Building the tensor of tangential forces
    Fabric_T(1,1:2) = Ltik*nik(1)*ikt(1:2) + Fabric_T(1,1:2)
    Fabric_T(2,1:2) = Ltik*nik(2)*ikt(1:2) + Fabric_T(2,1:2)
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
     write(109,*) '#   time     ', ' 2(X3-1)/(X1+3)n ', ' 2(X3-1)/(X1+3)f '
  end if

  write(109,'(3(1X,E12.5))') time, 2*(X1n-X2n)/(X1n+X2n), 2*(X1f-X2f)/(X1f+X2f)

  print*, 'Write Anisotropy Branch         ---> Ok!'
end subroutine branch_anisotropy


!==============================================================================
! Contacts orientation distribution
!==============================================================================
subroutine ctc_distribution

  implicit none

  integer                                      :: i, j, ninter, contotal, cd, an
  real(kind=8)                                 :: interval
  real(kind=8)                                 :: angcurr, maxint, minint, Rnik
  real(kind=8), dimension(2)                   :: nik, Lik
  real(kind=8), dimension(:,:), allocatable    :: theta_interval
  character(len=27)                            :: file_name
  logical                                      :: dir_ctcdir

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
  file_name       =  './POSTPRO/CTCDIR/CDIR.     '

  if (compteur_clout<10) then
    WRITE(file_name(23:24),'(I1)')   compteur_clout
  else if ( (compteur_clout>=10) .and. (compteur_clout<100) ) then
    WRITE(file_name(23:25),'(I2)')   compteur_clout
  else if ( (compteur_clout>=100).and. (compteur_clout<1000) ) then
    WRITE(file_name(23:26),'(I3)')   compteur_clout
  else if ( (compteur_clout>=1000).and. (compteur_clout<10000) ) then
    WRITE(file_name(23:27),'(I4)')   compteur_clout
  end if
  ! Opening the file
  open(unit=110,file=file_name,status='replace')

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
    do j=1,nb_ligneCONTACT_BODY
      cd = TAB_CONTACT_BODY(j)%icdent
      an = TAB_CONTACT_BODY(j)%ianent

      ! Only between discs
      if(TAB_CONTACT_BODY(j)%nature == 'CLJCx') cycle

      Rnik = (TAB_CONTACT_BODY(i)%F_tot(1)**2+TAB_CONTACT_BODY(i)%F_tot(2)**2)**0.5

      ! Checking if it is an active contact
      if (Rnik .le. 1.E-9) cycle

      ! Computing the branch vector
      Lik(1) = TAB_MESH(cd)%center(1)-TAB_MESH(an)%center(1)
      Lik(2) = TAB_MESH(cd)%center(2)-TAB_MESH(an)%center(2)

      ! The normal contact vector
      nik  = Lik/(Lik(1)**2+Lik(2)**2)**0.5

      ! In the plane xy!
      ! The 1 and 2 quadrants only
      ! Changing the normal vector direction if necesary
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

  contotal = 0

  ! Computing the number of active contacts
  do i=1, ninter
    contotal = contotal + theta_interval(i,2)
  end do

  ! Normalizing
  theta_interval(:,3) = theta_interval(:,2)/(2*contotal)!*interval)????

  ! Writing the heading
  write(110,*) '# Theta(Rad) ', '  Freq_norm  '

  ! Writing
  do i=1, ninter
    write(110,'(3(1X,E12.5))') theta_interval(i,1), theta_interval(i,2), theta_interval(i,3)
  end do

  ! Writing the 3rd and 4th quadrants
  do i=1, ninter
    write(110,'(3(1X,E12.5))') theta_interval(i,1)+pi, theta_interval(i,2), theta_interval(i,3)
  end do

  close(110)

  print*, 'Write Contacts distribution     ---> Ok!'

end subroutine ctc_distribution

!==============================================================================
! Contacts orientation
!==============================================================================
subroutine force_distribution

  implicit none

  integer                                      :: i, j, ninter, contotal, cd, an
  real(kind=8)                                 :: interval, pro_n, pro_t, Rnik, Rtik
  real(kind=8)                                 :: angcurr, maxint, minint
  real(kind=8), dimension(:,:), allocatable    :: theta_interval
  real(kind=8), dimension(2)                   :: nik, Fik, tik, tanik, tanik_plane, Lik
  character(len=27)                            :: file_name
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
  file_name       =  './POSTPRO/FRCDIR/FDIR.     '

  if (compteur_clout<10) then
    WRITE(file_name(23:24),'(I1)')   compteur_clout
  else if ( (compteur_clout>=10) .and. (compteur_clout<100) ) then
    WRITE(file_name(23:25),'(I2)')   compteur_clout
  else if ( (compteur_clout>=100).and. (compteur_clout<1000) ) then
    WRITE(file_name(23:26),'(I3)')   compteur_clout
  else if ( (compteur_clout>=1000).and. (compteur_clout<10000) ) then
    WRITE(file_name(23:27),'(I4)')   compteur_clout
  end if

  ! Opening the file
  open(unit=111,file=file_name,status='replace')

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
    do j=1,nb_ligneCONTACT_BODY
      cd = TAB_CONTACT_BODY(j)%icdent
      an = TAB_CONTACT_BODY(j)%ianent
      ! Only between discs
      if(TAB_CONTACT_BODY(j)%nature == 'CLJCx') cycle

      Rnik = (TAB_CONTACT_BODY(i)%F_tot(1)**2+TAB_CONTACT_BODY(i)%F_tot(2)**2)**0.5

      ! Checking if it is an active contact
      if (Rnik .le. 1.E-9) cycle

      ! Computing the branch vector
      Lik(1) = TAB_MESH(cd)%center(1)-TAB_MESH(an)%center(1)
      Lik(2) = TAB_MESH(cd)%center(2)-TAB_MESH(an)%center(2)

      ! The normal contact vector
      nik  = Lik/(Lik(1)**2+Lik(2)**2)**0.5

      tik(1) = nik(2)
      tik(2) = -nik(1)

      Rtik = TAB_CONTACT_BODY(i)%F_tot(1)*tik(1) + TAB_CONTACT_BODY(i)%F_tot(2)*tik(2)

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
      pro_t = abs(Fik(1)*tik(1) + Fik(2)*tik(2))
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
      angcurr = acos(nik(1))

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
  write(111,*) '# Theta(Rad) ', '    <fn>    ', '    <ft>    '

  ! Writing
  do i=1, ninter
    write(111,'(3(1X,E12.5))') theta_interval(i,1), theta_interval(i,3), theta_interval(i,4)
  end do

  ! Writing the 3rd and 4th quadrants
  do i=1, ninter
    write(111,'(3(1X,E12.5))') theta_interval(i,1)+pi, theta_interval(i,3), theta_interval(i,4)
  end do

  close(111)

  print*, 'Write Forces distribution       ---> Ok!'

end subroutine force_distribution


!==============================================================================
! Branch orientation distribution
!==============================================================================
subroutine branch_distribution

  implicit none

  integer                                      :: i, j, ninter, contotal, cd, an
  real(kind=8)                                 :: interval, pro_n, pro_t
  real(kind=8)                                 :: angcurr, maxint, minint, Rnik
  real(kind=8), dimension(:,:), allocatable    :: theta_interval
  real(kind=8), dimension(2)                   :: nik, Lik, tik, tanik
  character(len=27)                            :: file_name
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
  file_name       =  './POSTPRO/BRCDIR/BDIR.     '

  if (compteur_clout<10) then
    WRITE(file_name(23:24),'(I1)')   compteur_clout
  else if ( (compteur_clout>=10) .and. (compteur_clout<100) ) then
    WRITE(file_name(23:25),'(I2)')   compteur_clout
  else if ( (compteur_clout>=100).and. (compteur_clout<1000) ) then
    WRITE(file_name(23:26),'(I3)')   compteur_clout
  else if ( (compteur_clout>=1000).and. (compteur_clout<10000) ) then
    WRITE(file_name(23:27),'(I4)')   compteur_clout
  end if
  ! Opening the file
  open(unit=112,file=file_name,status='replace')

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
    do j=1,nb_ligneCONTACT_BODY
      cd = TAB_CONTACT_BODY(j)%icdent
      an = TAB_CONTACT_BODY(j)%ianent
      ! Only between discs
      if(TAB_CONTACT_BODY(j)%nature == 'CLJCx') cycle

      Rnik = (TAB_CONTACT_BODY(i)%F_tot(1)**2+TAB_CONTACT_BODY(i)%F_tot(2)**2)**0.5

      ! Checking if it is an active contact
      if (Rnik .le. 1.E-9) cycle

      ! Computing the branch vector
      Lik(1) = TAB_MESH(cd)%center(1)-TAB_MESH(an)%center(1)
      Lik(2) = TAB_MESH(cd)%center(2)-TAB_MESH(an)%center(2)

      ! The normal contact vector
      nik  = Lik/(Lik(1)**2+Lik(2)**2)**0.5

      ! The tangential vector - xz plane!
      !tik = TAB_CONTACT(j)%t
      !sik = TAB_CONTACT(j)%s
      !tanik = sik*TAB_CONTACT(j)%rs + tik*TAB_CONTACT(j)%rt
      tanik(1) = nik(2)
      tanik(2) = -nik(1)

      ! unitary
      !tanik = tanik/sqrt(tanik(1)**2 + tanik(2)**2 + tanik(3)**2)

      ! Normal projection
      pro_n = abs(Lik(1)*nik(1) + Lik(2)*nik(2))

      ! Tangential projection
      pro_t = abs(Lik(1)*tanik(1) + Lik(2)*tanik(2))

      ! The 1 and 2 quadrants only
      ! Changing the vector direction if necessary
      if(nik(2) .lt. 0) then
        nik(1) = -nik(1)
        nik(2) = -nik(2)
      end if

      ! Finding the current angle. Note the vector nik is already unitary!
      angcurr = acos(nik(1))

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
  write(112,*) '# Theta(Rad) ', '    <ln>    ', '    <lt>    '

  ! Writing
  do i=1, ninter
    write(112,'(3(1X,E12.5))') theta_interval(i,1), theta_interval(i,3), theta_interval(i,4)
  end do

  ! Writing the 3rd and 4th quadrants
  do i=1, ninter
    write(112,'(3(1X,E12.5))') theta_interval(i,1)+pi, theta_interval(i,3), theta_interval(i,4)
  end do

  close(112)

  print*, 'Write Branch distribution       ---> Ok!'

end subroutine branch_distribution

!==============================================================================
! Computing q/p
!==============================================================================
subroutine qoverp

  implicit none

  real(kind=8),dimension(2,2)              :: Moment
  real(kind=8),dimension(2)                :: nik,tik,Lik,Fik
  real(kind=8)                             :: Rtik,Rnik
  integer                                  :: i,cd,an
  integer                                  :: ierror,matz,lda
  real(kind=8),dimension(2)                :: wr,wi
  real(kind=8)                             :: S1,S2
  real(kind=8),dimension(2,2)              :: localframe

  Moment(:,:) = 0.0

  do i=1,nb_ligneCONTACT_BODY
    cd = TAB_CONTACT_BODY(i)%icdent
    an = TAB_CONTACT_BODY(i)%ianent

    ! Contacts between particles
    if  (TAB_CONTACT_BODY(i)%nature == 'CLJCx') cycle

    ! Only active contacts
    Rnik = (TAB_CONTACT_BODY(i)%F_tot(1)**2+TAB_CONTACT_BODY(i)%F_tot(2)**2)**0.5

    if (Rnik .le. 1.0D-9) cycle

    ! Computing the branch vector
    Lik(1) = TAB_MESH(cd)%center(1)-TAB_MESH(an)%center(1)
    Lik(2) = TAB_MESH(cd)%center(2)-TAB_MESH(an)%center(2)

    ! The normal contact vector
    nik  = Lik/(Lik(1)**2+Lik(2)**2)**0.5

    tik(1) = nik(2)
    tik(2) = -nik(1)

    Rtik = TAB_CONTACT_BODY(i)%F_tot(1)*tik(1) + TAB_CONTACT_BODY(i)%F_tot(2)*tik(2)

    Fik(1) = (Rnik*nik(1)+Rtik*tik(1))
    Fik(2) = (Rnik*nik(2)+Rtik*tik(2))

    Moment(1,1:2) = Fik(1)*Lik(1:2) + Moment(1,1:2)
    Moment(2,1:2) = Fik(2)*Lik(1:2) + Moment(2,1:2)
  end do

  Moment = Moment / (height*length)

  lda  = 2
  matz = 1

  ! Finding the eigenvalues of the global stress tensor
  call rg(lda, 2, Moment, wr, wi, matz, localframe, ierror)
  S1 = max( wr(1),wr(2) )
  S2 = min( wr(1),wr(2) )

  if (first_over_all) then
    write(113,*) '#   time     ', '      S1     ', '      S2     ', '(S1-S2)/(S1+S2)'
  endif

  write(113,'(4(1X,E12.5))') time, S1,S2,(S1-S2)/(S1+S2)

  print*, 'Write QoverP                    ---> Ok!'
end subroutine qoverp

!==============================================================================
! Computing q/p with walls forces. Note that only works with frictionless walls
!==============================================================================
subroutine qp_walls_nf

  implicit none

  integer                                  :: i,cd,an
  real(kind=8)                             :: Rnik
  real(kind=8)                             :: S11,S22, S1, S2
  real(kind=8),dimension(4)                :: walls_fn

  ! Initializing
  walls_fn(:) = 0.0

  do i=1,nb_ligneCONTACT
    cd = TAB_CONTACT(i)%icdent
    an = TAB_CONTACT(i)%ianent

    ! Contacts between particles
    !print*, TAB_CONTACT(i)%nature
    !print*, TAB_CONTACT(i)%rn

    if  (TAB_CONTACT(i)%nature == 'CLJCx') then

      ! Only active contacts
      Rnik = abs(TAB_CONTACT(i)%rn)
      !(TAB_CONTACT_BODY(i)%F_tot(1)**2+TAB_CONTACT_BODY(i)%F_tot(2)**2)**0.5

      if (Rnik .le. 1.0D-9) cycle

      if (an == 1) then
        walls_fn(1) = walls_fn(1) + Rnik
      else if(an == 2) then
        walls_fn(2) = walls_fn(2) + Rnik
      else if(an == 3) then
        walls_fn(3) = walls_fn(3) + Rnik
      else if(an == 4) then
        walls_fn(4) = walls_fn(4) + Rnik
      end if
    end if
  end do

  ! Note that this is very specific! I know the order of walls in the bodies file

  ! X direction
  S11 = (walls_fn(3) + walls_fn(4))/2.
  S11 = S11/height

  ! Y direction
  S22 = (walls_fn(1) + walls_fn(2))/2.
  S22 = S22/length

  S1 = max(S11,S22)
  S2 = min(S11,S22)

  print*, S1
  print*, S2

  if (first_over_all) then
    write(114,*) '#   time     ', '      S11    ', '      S22    ', '(S1-S2)/(S1+S2)'
  endif

  if (S1 .gt. 1.0D-9 .and. S2 .gt. 1.0D-9) then
    write(114,'(4(1X,E12.5))') time, S11,S22,(S1-S2)/(S1+S2)
  end if

  print*, 'Write QP frictionless walls     ---> Ok!'
end subroutine qp_walls_nf


!==============================================================================
! Checking sum of contact forces on particles
!==============================================================================
subroutine check_force

  implicit none

  integer                                        :: i, cd, an
  real*8, dimension(:,:), allocatable            :: f_sum


  if (allocated(f_sum)) deallocate(f_sum)
  allocate(f_sum(n_def_particles+n_walls,2))

  f_sum(:,:) = 0

  do i=1, nb_ligneCONTACT_BODY
    cd = TAB_CONTACT_BODY(i)%icdent
    an = TAB_CONTACT_BODY(i)%ianent

    f_sum(cd,1) = f_sum(cd,1) + TAB_CONTACT_BODY(i)%F_tot(1)
    f_sum(cd,2) = f_sum(cd,2) + TAB_CONTACT_BODY(i)%F_tot(2)

    if (TAB_CONTACT_BODY(i)%nature == 'CLJCx') then
      f_sum(an+n_def_particles,1) = f_sum(an + n_def_particles,1) - TAB_CONTACT_BODY(i)%F_tot(1)
      f_sum(an+n_def_particles,2) = f_sum(an + n_def_particles,2) - TAB_CONTACT_BODY(i)%F_tot(2)
    else
      f_sum(an,1) = f_sum(an,1) - TAB_CONTACT_BODY(i)%F_tot(1)
      f_sum(an,2) = f_sum(an,2) - TAB_CONTACT_BODY(i)%F_tot(2)
    end if
  end do

  do i=1, n_def_particles
    print*, i, ' ', f_sum(i,1), ' ', f_sum(i,2)
  end do

end subroutine check_force


!==============================================================================
! Computing the average strain tensor per grain
!==============================================================================
subroutine strain_grain

  implicit none

  integer                                  :: i, j
  real*8, dimension(2,2)                   :: avg_part
  logical                                  :: dir_c
  character(len=32)                        :: l_file


! Cleaning or creating the folder if necessary
  if (first_over_all) then
    ! Asking if the file already exists
    inquire(file='./POSTPRO/LIST_ST', exist=dir_c)
    if(dir_c) then
      ! Cleaning
      call system('rm ./POSTPRO/LIST_ST/*')
    else
      ! Creating
      call system('mkdir ./POSTPRO/LIST_ST')
    end if
  end if

  ! The name of the file
  l_file = './POSTPRO/LIST_ST/F_Strain.     '

  ! Preparing the corresponding number of the file
  if (compteur_clout<10) then
    WRITE(l_file(28:29),'(I1)')   compteur_clout
  else if ( (compteur_clout>=10) .and. (compteur_clout<100) ) then
    WRITE(l_file(28:30),'(I2)')   compteur_clout
  else if ( (compteur_clout>=100).and. (compteur_clout<1000) ) then
    WRITE(l_file(28:31),'(I3)')   compteur_clout
  else if ( (compteur_clout>=1000).and. (compteur_clout<10000) ) then
    WRITE(l_file(28:32),'(I4)')   compteur_clout
  else
    print*, "Not enough file names :: List of grain strain"
    stop
  end if

  ! Opening the file
  open(unit=122,file=l_file,status='replace')

  write(122,*) '#    E_xx    ', '    E_xy    ', '    E_yx    ', '    E_yy    '

  do i=1, n_def_particles
    ! Initializing
    avg_part(:,:) = 0.0
    do j=1, TAB_MESH(i)%nb_mesh
      avg_part(1,1) = avg_part(1,1) + TAB_MESH(i)%mesh_connect(j)%strain_gvp(1,1)
      avg_part(1,2) = avg_part(1,2) + TAB_MESH(i)%mesh_connect(j)%strain_gvp(1,2)
      avg_part(2,1) = avg_part(2,1) + TAB_MESH(i)%mesh_connect(j)%strain_gvp(2,1)
      avg_part(2,2) = avg_part(2,2) + TAB_MESH(i)%mesh_connect(j)%strain_gvp(2,2)
    end do
    avg_part = avg_part/TAB_MESH(i)%nb_mesh

    ! And writing
    write(122, '(4(1X,E12.5))') avg_part(1,1), avg_part(1,2), avg_part(2,1), avg_part(2,2)
  end do

  ! Closing the port
  close(122)

  print*, 'Computing strains    ---> Ok!'

end subroutine strain_grain


!================================================
! Drawing in vtk
!================================================
subroutine draw

  implicit none

  integer                                        ::  i, j, k
  integer                                        ::  n_l_vertices, n_l_faces, curr_l_faces
  integer                                        ::  n_l_fields, l_counter
  integer                                        ::  l_cdt, l_ant
  real(kind=8), dimension(3)                     ::  curr_l_vector, Lik
  character(len=8)                               ::  vtk_c_temp,  v_n_vertices
  character(:), allocatable                      ::  vtk_part, vtk_inter, vtk_counter
  character(:), allocatable                      ::  vtk_n_points
  logical                                        ::  dir_vtk

  ! Variables for the forces
  integer                                        ::  f_counter
  real(kind=8)                                   ::  vtk_ave_force, l_force_scale, l_force
  real(kind=8)                                   ::  l_rn, l_rt, l_ave_rad
  real(kind=8), dimension(3)                     ::  vtk_cd_center, vtk_an_center, v_l_normal, v_l_t
  character(:), allocatable                      ::  vtk_forces
  character(:), allocatable                      ::  vtk_n_v_forces
  character(len=8)                               ::  v_n_v_forces

  ! Variables for the contacts points
  integer                                        ::  n_l_sides
  real*8                                         ::  cur_angle_rad, rad_ctci
  real*8, dimension(2)                           ::  cur_disc_vertex


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

  vtk_part = './POSTPRO/VTK/meshx_' // vtk_counter // '.vtk'

  open(unit=121, file=vtk_part, status='replace')
  write(121,'(A)') '# vtk DataFile Version 3.0'
  write(121,'(A,I6)') 'RIGID ', compteur_clout
  write(121,'(A)') 'ASCII'
  write(121,'(A)') 'DATASET POLYDATA'

  ! Writing the number of vertices and faces
  n_l_vertices = 0
  n_l_faces = 0
  do i=1, n_def_particles
    n_l_vertices = n_l_vertices + TAB_MESH(i)%nb_node
    n_l_faces = n_l_faces + TAB_MESH(i)%nb_mesh
  end do

  ! Adding the vertices of the walls
  n_l_vertices = n_l_vertices + n_walls*4
  n_l_faces = n_l_faces + n_walls

  write(v_n_vertices, '(I8)') n_l_vertices

  vtk_n_points = 'POINTS ' // v_n_vertices // ' float'
  write(121,'(A)') vtk_n_points

  ! Writing the coordinates of the vertices of particles and walls
  ! 3 vertices per row
  k=0
  do i=1, n_def_particles + n_walls
    !write(121,*) 'Particle', i
    if (i .le. n_def_particles) then
      do j=1, TAB_MESH(i)%nb_node
        k = k + 1
        if(k .lt. 3) then
          write(121,'(E13.6,A)', advance = 'no') TAB_MESH(i)%node(j,1), ' '
          write(121,'(E13.6,A)', advance = 'no') TAB_MESH(i)%node(j,2), ' '
          write(121,'(E13.6,A)', advance = 'no') TAB_MESH(i)%node(j,3), ' '
        else
          write(121,'(E13.6,A)', advance = 'no') TAB_MESH(i)%node(j,1), ' '
          write(121,'(E13.6,A)', advance = 'no') TAB_MESH(i)%node(j,2), ' '
          write(121,'(E13.6,A)')                 TAB_MESH(i)%node(j,3), ' '
          k = 0
        end if
      end do
    else
      ! Writing the vertices of the walls
      do j=1, 4
        k = k +1
        if(k .lt. 3) then
          write(121,'(E13.6,A)', advance = 'no') TAB_WALL(i-n_def_particles)%vertex(1,j), ' '
          write(121,'(E13.6,A)', advance = 'no') TAB_WALL(i-n_def_particles)%vertex(2,j), ' '
          write(121,'(E13.6,A)', advance = 'no') TAB_WALL(i-n_def_particles)%vertex(3,j), ' '
        else
          write(121,'(E13.6,A)', advance = 'no') TAB_WALL(i-n_def_particles)%vertex(1,j), ' '
          write(121,'(E13.6,A)', advance = 'no') TAB_WALL(i-n_def_particles)%vertex(2,j), ' '
          write(121,'(E13.6,A)')                 TAB_WALL(i-n_def_particles)%vertex(3,j), ' '
          k = 0
        end if
      end do
    end if
  end do

  write(121, '(A)') ' '

  ! Writing the conectivity between vertices
  ! Plus 4 faces each wall
  write(121,'(A)', advance='no') 'POLYGONS '
  write(121, '(2(I8,A))') (n_l_faces) , ' ' , (n_l_faces*4)+n_walls

  curr_l_faces = 0

  do i=1, n_def_particles
    !write(121,*) 'Particle', i
    do j=1, TAB_MESH(i)%nb_mesh
      ! We have always triangles
      write(121, '(I1,A)', advance= 'no') 3, ' '

      ! First ...
      write(121, '(I7,A)', advance = 'no') TAB_MESH(i)%mesh_connect(j)%connectivity(1) + curr_l_faces -1, ' '

      ! Second
      write(121, '(I7,A)', advance = 'no') TAB_MESH(i)%mesh_connect(j)%connectivity(2) + curr_l_faces -1, ' '

      ! Third
      write(121, '(I7,A)') TAB_MESH(i)%mesh_connect(j)%connectivity(3) + curr_l_faces -1, ' '
    end do
    curr_l_faces = curr_l_faces + TAB_MESH(i)%nb_node
  end do

  ! Writing the conectivity for the walls
  ! For all the walls
  do i=1, n_walls
    ! Each one with 6 faces

    !!!!! 1st (1-2-3-4)
    write(121, '(I1,A)', advance= 'no') 4, ' '
    write(121, '(I7,A)', advance = 'no') curr_l_faces, ' '
    write(121, '(I7,A)', advance = 'no') curr_l_faces +1, ' '
    write(121, '(I7,A)', advance = 'no') curr_l_faces +2, ' '
    write(121, '(I7,A)') curr_l_faces +3

    curr_l_faces = curr_l_faces + 4
  end do

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!! CELL DATA  !!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 ! How many fields there will be?
  n_l_fields = 3                          ! They are: Id, Material, Z (coordination)
                                          ! ...

  ! A blank space
  write(121,'(A)') ''
  ! The cells begin
  write(121,'(A)', advance='no') 'CELL_DATA '

  ! Writing the number of data by field. It corresponds to the same number of faces
  write(121, '(2(I6,A))') n_l_faces
  write(121, '(A,I4)')  'FIELD FieldData', n_l_fields

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Id
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Naming the field, the dimension of the data, and number of lines (as before FIELD), and data type
  write(121,'(A,I1,I8,A)') 'Id ', 1, n_l_faces, ' float'
  k = 0
  l_counter = 0
  ! The particles
  do i=1, n_def_particles
    l_counter = l_counter + 1

    do j=1, TAB_MESH(i)%nb_mesh
      k=k+1
      if(k .le. 9) then
        write(121, '(I6)', advance='no') l_counter
      else
        write(121, '(I6)') l_counter
        k=0
      end if
    end do
  end do
  ! The walls
  do i=1, n_walls
    l_counter = l_counter + 1
    k=k+1
    if(k .le. 9) then
      write(121, '(I6)', advance='no') l_counter
    else
      write(121, '(I6)') l_counter
      k=0
    end if
  end do
  ! And jump
  write(121, '(A)') ' '


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Material
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(121,'(A,I1,I8,A)') 'Material ', 1, n_l_faces, ' float'
  k = 0
  l_counter = 1
  ! The particles
  do i=1, n_def_particles
    ! Material number 1
    do j=1, TAB_MESH(i)%nb_mesh
      k=k+1
      if(k .le. 9) then
        write(121, '(I6)', advance='no') l_counter
      else
        write(121, '(I6)') l_counter
        k=0
      end if
    end do
  end do
  l_counter = l_counter + 1
  ! The walls
  do i=1, n_walls
  ! Material number 2
    k=k+1
    if(k .le. 9) then
      write(121, '(I6)', advance='no') l_counter
    else
      write(121, '(I6)') l_counter
      k=0
    end if
  end do

  ! And jump
  write(121, '(A)') ' '

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Coordination
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  l_counter = 0

  write(121,'(A,I1,I8,A)') 'Z ', 1, n_l_faces, ' float'
  k = 0
  ! The particles
  do i=1, n_def_particles
    do j=1, TAB_MESH(i)%nb_mesh
      k=k+1
      if(k .le. 9) then
        write(121, '(I4,A)', advance='no') TAB_MESH(i)%nctc , ' '
      else
        write(121, '(I4,A)') TAB_MESH(i)%nctc
        k=0
      end if
    end do
  end do
  ! The walls
  do j=1, n_walls
    k=k+1
    if(k .le. 9) then
      write(121, '(I4,A)', advance='no') 0, ' '
    else
      write(121, '(I4,A)') 0
      k=0
    end if
  end do

  ! And jump
  write(121, '(A)') ' '


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!! POINT DATA  !!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 ! How many fields there will be?
 if (c_mesh_stresses==1) then
    n_l_fields = 8                          ! They are: Displacement, velocity,
                                            ! Cauchy-Green_right, Second_Piola_Kirchhoff
                                            ! Cauchy-Green left, Deformation gradient
                                            ! Green_E, Cauchy_s
  else
    n_l_fields = 2
  end if

  ! A blank space
  write(121,'(A)') ''
  ! The cells begin
  write(121,'(A)', advance='no') 'POINT_DATA '

  ! Writing the number of data by field. It corresponds to the same number of faces
  write(121, '(2(I8,A))') n_l_vertices
  write(121, '(A,I4)')  'FIELD FieldData', n_l_fields

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Dispacement
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  curr_l_vector(:) = 0.D0

  write(121,'(A,I1,I8,A)') 'Disp ', 3, n_l_vertices, ' float'

  k = 0
  do i=1, n_def_particles + n_walls
    if (i .le. n_def_particles) then
      do j=1, TAB_MESH(i)%nb_node
        curr_l_vector(:) = TAB_MESH(i)%node(j,:) - TAB_MESH(i)%node_ref(j,:)
        k=k+3
        if(k .le. 9) then
          write(121, '(3(E13.6,A))', advance='no') curr_l_vector(1), ' ', curr_l_vector(2), ' ', curr_l_vector(3), ' '
        else
          write(121, '(3(E13.6,A))') curr_l_vector(1), ' ', curr_l_vector(2), ' ', curr_l_vector(3)
          k=0
        end if
      end do
    else
      curr_l_vector(:) = TAB_WALL(i-n_def_particles)%center(:) - TAB_WALL(i-n_def_particles)%center_ref(:)
      do j=1, 4
        k=k+3
        if(k .le. 9) then
          write(121, '(3(E13.6,A))', advance='no') curr_l_vector(1), ' ', curr_l_vector(2), ' ', curr_l_vector(3), ' '
        else
          write(121, '(3(F13.6,A))') curr_l_vector(1), ' ', curr_l_vector(2), ' ', curr_l_vector(3)
          k=0
        end if
      end do
    end if
  end do

  ! And jump
  write(121, '(A)') ' '

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Velocity
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  curr_l_vector(:) = 0.D0

  write(121,'(A,I1,I8,A)') 'Veloc ', 3, n_l_vertices, ' float'
  k = 0
  do i=1, n_def_particles + n_walls
    if (i .le. n_def_particles) then
      do j=1, TAB_MESH(i)%nb_node
        curr_l_vector(:) = TAB_MESH(i)%veloc(j,:)
        k=k+3
        if(k .le. 9) then
          write(121, '(3(E13.6,A))', advance='no') curr_l_vector(1), ' ', curr_l_vector(2), ' ', curr_l_vector(3), ' '
        else
          write(121, '(3(E13.6,A))') curr_l_vector(1), ' ', curr_l_vector(2), ' ', curr_l_vector(3)
          k=0
        end if
      end do
    else
      curr_l_vector(:) = TAB_WALL(i-n_def_particles)%veloc(:)
      do j=1, 4
        k=k+3
        if(k .le. 9) then
          write(121, '(3(E13.6,A))', advance='no') curr_l_vector(1), ' ', curr_l_vector(2), ' ', curr_l_vector(3), ' '
        else
          write(121, '(3(E13.6,A))') curr_l_vector(1), ' ', curr_l_vector(2), ' ', curr_l_vector(3)
          k=0
        end if
      end do
    end if
  end do

  ! And jump
  write(121, '(A)') ' '

  if (c_mesh_stresses == 1) then

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Second Piola-Kirchhoff stress tensor
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    write(121,'(A,I1,I8,A)') '2PK_stress ', 3, n_l_vertices, ' float'
    k = 0
    do i=1, n_def_particles + n_walls
      if (i .le. n_def_particles) then
        do j=1, TAB_MESH(i)%nb_node
          k=k+3
          if(k .le. 9) then
            write(121, '(3(E13.6,A))', advance='no') TAB_MESH(i)%node_info(j)%pk_2nd(1,1), ' ', &
                                                     TAB_MESH(i)%node_info(j)%pk_2nd(1,2), ' ', &
                                                     !TAB_MESH(i)%node_info(j)%pk_2nd(2,1), ' ', &
                                                     TAB_MESH(i)%node_info(j)%pk_2nd(2,2), ' '
          else
            write(121, '(3(E13.6,A))') TAB_MESH(i)%node_info(j)%pk_2nd(1,1), ' ', &
                                       TAB_MESH(i)%node_info(j)%pk_2nd(1,2), ' ', &
                                       !TAB_MESH(i)%node_info(j)%pk_2nd(2,1), ' ', &
                                       TAB_MESH(i)%node_info(j)%pk_2nd(2,2), ' '
            k=0
          end if
        end do
      else
        do j=1, 4
          k=k+3
          if(k .le. 9) then
            write(121, '(3(E13.6,A))', advance='no') 0.0, ' ', 0.0, ' ', 0.0, ' '!, 0.0, ' '
          else
            write(121, '(3(E13.6,A))') 0.0, ' ', 0.0, ' ', 0.0, ' '!, 0.0, ' '
            k=0
          end if
        end do
      end if
    end do

    ! And jump
    write(121, '(A)') ' '

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Right Cauchy-Green deformation tensor
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    write(121,'(A,I1,I8,A)') 'CG_right ', 3 , n_l_vertices, ' float'
    k = 0
    do i=1, n_def_particles + n_walls
      if (i .le. n_def_particles) then
        do j=1, TAB_MESH(i)%nb_node
          k=k+3
          if(k .le. 9) then
            write(121, '(3(E13.6,A))', advance='no') TAB_MESH(i)%node_info(j)%cg_right(1,1), ' ', &
                                                     TAB_MESH(i)%node_info(j)%cg_right(1,2), ' ', &
                                                     !TAB_MESH(i)%node_info(j)%cg_right(2,1), ' ', &
                                                     TAB_MESH(i)%node_info(j)%cg_right(2,2), ' '
          else
            write(121, '(3(E13.6,A))') TAB_MESH(i)%node_info(j)%cg_right(1,1), ' ', &
                                       TAB_MESH(i)%node_info(j)%cg_right(1,2), ' ', &
                                       !TAB_MESH(i)%node_info(j)%cg_right(2,1), ' ', &
                                       TAB_MESH(i)%node_info(j)%cg_right(2,2), ' '
            k=0
          end if
        end do
      else
        do j=1, 4
          k=k+3
          if(k .le. 9) then
            write(121, '(3(E13.6,A))', advance='no') 0.0, ' ', 0.0, ' ', 0.0, ' '!, 0.0, ' '
          else
            write(121, '(3(E13.6,A))') 0.0, ' ', 0.0, ' ', 0.0, ' '!, 0.0, ' '
            k=0
          end if
        end do
      end if
    end do

    ! And jump
    write(121, '(A)') ' '


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Left Cauchy-Green deformation tensor
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    write(121,'(A,I1,I8,A)') 'CG_left ', 3 , n_l_vertices, ' float'
    k = 0
    do i=1, n_def_particles + n_walls
      if (i .le. n_def_particles) then
        do j=1, TAB_MESH(i)%nb_node
          k=k+3
          if(k .le. 9) then
            write(121, '(3(E13.6E2,A))', advance='no') TAB_MESH(i)%node_info(j)%cg_left(1,1), ' ', &
                                                     TAB_MESH(i)%node_info(j)%cg_left(1,2), ' ', &
                                                     !TAB_MESH(i)%node_info(j)%cg_left(2,1), ' ', &
                                                     TAB_MESH(i)%node_info(j)%cg_left(2,2), ' '
          else
            write(121, '(3(E13.6E2,A))') TAB_MESH(i)%node_info(j)%cg_left(1,1), ' ', &
                                       TAB_MESH(i)%node_info(j)%cg_left(1,2), ' ', &
                                       !TAB_MESH(i)%node_info(j)%cg_left(2,1), ' ', &
                                       TAB_MESH(i)%node_info(j)%cg_left(2,2), ' '
            k=0
          end if
        end do
      else
        do j=1, 4
          k=k+3
          if(k .le. 9) then
            write(121, '(3(E13.6E2,A))', advance='no') 0.0, ' ', 0.0, ' ', 0.0, ' '!, 0.0, ' '
          else
            write(121, '(3(E13.6E2,A))') 0.0, ' ', 0.0, ' ', 0.0, ' '!, 0.0, ' '
            k=0
          end if
        end do
      end if
    end do

    ! And jump
    write(121, '(A)') ' '

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Deformation gradient
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    write(121,'(A,I1,I8,A)') 'F_ ', 3 , n_l_vertices, ' float'
    k = 0
    do i=1, n_def_particles + n_walls
      if (i .le. n_def_particles) then
        do j=1, TAB_MESH(i)%nb_node
          k=k+3
          if(k .le. 9) then
            write(121, '(3(E13.6E2,A))', advance='no') TAB_MESH(i)%node_info(j)%def_grad(1,1), ' ', &
                                                     TAB_MESH(i)%node_info(j)%def_grad(1,2), ' ', &
                                                     !TAB_MESH(i)%node_info(j)%def_grad(2,1), ' ', &
                                                     TAB_MESH(i)%node_info(j)%def_grad(2,2), ' '
          else
            write(121, '(3(E13.6E2,A))') TAB_MESH(i)%node_info(j)%def_grad(1,1), ' ', &
                                       TAB_MESH(i)%node_info(j)%def_grad(1,2), ' ', &
                                       !TAB_MESH(i)%node_info(j)%def_grad(2,1), ' ', &
                                       TAB_MESH(i)%node_info(j)%def_grad(2,2), ' '
            k=0
          end if
        end do
      else
        do j=1, 4
          k=k+3
          if(k .le. 9) then
            write(121, '(3(E13.6E2,A))', advance='no') 0.0, ' ', 0.0, ' ', 0.0, ' '!, 0.0, ' '
          else
            write(121, '(3(E13.6E2,A))') 0.0, ' ', 0.0, ' ', 0.0, ' '!, 0.0, ' '
            k=0
          end if
        end do
      end if
    end do

    ! And jump
    write(121, '(A)') ' '

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Green_E
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    write(121,'(A,I1,I8,A)') 'E ', 3 , n_l_vertices, ' float'
    k = 0
    do i=1, n_def_particles + n_walls
      if (i .le. n_def_particles) then
        do j=1, TAB_MESH(i)%nb_node
          k=k+3
          if(k .le. 9) then
            write(121, '(3(E13.6E2,A))', advance='no') TAB_MESH(i)%node_info(j)%green_E(1,1), ' ', &
                                                     TAB_MESH(i)%node_info(j)%green_E(1,2), ' ', &
                                                     !TAB_MESH(i)%node_info(j)%def_grad(2,1), ' ', &
                                                     TAB_MESH(i)%node_info(j)%green_E(2,2), ' '
          else
            write(121, '(3(E13.6E2,A))') TAB_MESH(i)%node_info(j)%green_E(1,1), ' ', &
                                       TAB_MESH(i)%node_info(j)%green_E(1,2), ' ', &
                                       !TAB_MESH(i)%node_info(j)%def_grad(2,1), ' ', &
                                       TAB_MESH(i)%node_info(j)%green_E(2,2), ' '
            k=0
          end if
        end do
      else
        do j=1, 4
          k=k+3
          if(k .le. 9) then
            write(121, '(3(E13.6E2,A))', advance='no') 0.0, ' ', 0.0, ' ', 0.0, ' '!, 0.0, ' '
          else
            write(121, '(3(E13.6E2,A))') 0.0, ' ', 0.0, ' ', 0.0, ' '!, 0.0, ' '
            k=0
          end if
        end do
      end if
    end do

    ! And jump
    write(121, '(A)') ' '

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Cauchy S
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    write(121,'(A,I1,I8,A)') 'Sigma ', 3 , n_l_vertices, ' float'
    k = 0
    do i=1, n_def_particles + n_walls
      if (i .le. n_def_particles) then
        do j=1, TAB_MESH(i)%nb_node
          k=k+3
          if(k .le. 9) then
            write(121, '(3(E13.6E2,A))', advance='no') TAB_MESH(i)%node_info(j)%cauchy_s(1,1), ' ', &
                                                     TAB_MESH(i)%node_info(j)%cauchy_s(1,2), ' ', &
                                                     !TAB_MESH(i)%node_info(j)%def_grad(2,1), ' ', &
                                                     TAB_MESH(i)%node_info(j)%cauchy_s(2,2), ' '
          else
            write(121, '(3(E13.6E2,A))') TAB_MESH(i)%node_info(j)%cauchy_s(1,1), ' ', &
                                       TAB_MESH(i)%node_info(j)%cauchy_s(1,2), ' ', &
                                       !TAB_MESH(i)%node_info(j)%def_grad(2,1), ' ', &
                                       TAB_MESH(i)%node_info(j)%cauchy_s(2,2), ' '
            k=0
          end if
        end do
      else
        do j=1, 4
          k=k+3
          if(k .le. 9) then
            write(121, '(3(E13.6E2,A))', advance='no') 0.0, ' ', 0.0, ' ', 0.0, ' '!, 0.0, ' '
          else
            write(121, '(3(E13.6E2,A))') 0.0, ' ', 0.0, ' ', 0.0, ' '!, 0.0, ' '
            k=0
          end if
        end do
      end if
    end do

    ! And jump
    write(121, '(A)') ' '
  end if



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

  ! Writing the number of vertices for the forces. These are rectangles joining the center of
  ! particles in contact
  ! Four vertices for each rectangle
  f_counter = 0
  do i=1, nb_ligneCONTACT_BODY
    if(TAB_CONTACT_BODY(i)%nature == 'CLJCx') cycle
    f_counter = f_counter + 1
  end do

  write(v_n_v_forces, '(I8)') f_counter * 4

  vtk_n_v_forces = 'POINTS' // v_n_v_forces // ' float'
  write(121,'(A)') vtk_n_v_forces

  ! Finding the average force
  vtk_ave_force = 0.D0

  do i = 1, nb_ligneCONTACT_BODY
    if (TAB_CONTACT_BODY(i)%nature == 'CLJCx') cycle
    l_cdt = TAB_CONTACT_BODY(i)%icdent
    l_ant = TAB_CONTACT_BODY(i)%ianent
    ! The branch
    Lik(:) = TAB_MESH(l_cdt)%center(:)-TAB_MESH(l_ant)%center(:)
    !Normalizing
    Lik(:) = Lik(:)/(Lik(1)**2 + Lik(2)**2)**0.5
    ! The projection on the branch
    vtk_ave_force = vtk_ave_force + abs(TAB_CONTACT_BODY(i)%F_tot(1)*Lik(1) + &
                                        TAB_CONTACT_BODY(i)%F_tot(2)*Lik(2))
  end do

  ! The average force
  vtk_ave_force = vtk_ave_force/nb_ligneCONTACT_BODY

  ! The average particle radius
  l_ave_rad = 0
  do i=1, n_def_particles
    l_ave_rad = l_ave_rad + (TAB_MESH(i)%surface/pi)**0.5
  end do

  l_ave_rad = l_ave_rad/n_def_particles

  ! Force scale parameter
  l_force_scale = l_ave_rad*f_b_scale

  k=0
  do i=1, nb_ligneCONTACT_BODY

    if (TAB_CONTACT_BODY(i)%nature == 'CLJCx') cycle

    ! The particles ids
    l_cdt = TAB_CONTACT_BODY(i)%icdent
    l_ant = TAB_CONTACT_BODY(i)%ianent

    ! The center of each particle
    vtk_cd_center(:) = TAB_MESH(l_cdt)%center(:)
    vtk_an_center(:) = TAB_MESH(l_ant)%center(:)

    ! The branch
    Lik(:) = TAB_MESH(l_cdt)%center(:)-TAB_MESH(l_ant)%center(:)

    !Normalizing
    Lik(:) = Lik(:)/(Lik(1)**2 + Lik(2)**2)**0.5

    ! The normal and tangential vectors... are replaced following the branch
    v_l_normal(:) = Lik(:)
    v_l_t(1) = Lik(2)
    v_l_t(2) = -Lik(1)

    ! The forces
    l_rn = TAB_CONTACT_BODY(i)%F_tot(1)*v_l_normal(1) + TAB_CONTACT_BODY(i)%F_tot(2)*v_l_normal(2)
    l_rt = abs(TAB_CONTACT_BODY(i)%F_tot(1)*v_l_t(1) + TAB_CONTACT_BODY(i)%F_tot(2)*v_l_t(2))

    l_force = (l_rn/vtk_ave_force)*l_force_scale

    ! Write! ... 4 vertices for each line

    !  --- First vertex.
    write(121, '(3(E13.6,A))', advance='no') vtk_cd_center(1)+l_force*v_l_t(1), ' ', &
                                             vtk_cd_center(2)+l_force*v_l_t(2), ' ', &
                                             vtk_cd_center(3), ' '
    !  --- Second vertex.
    write(121, '(3(E13.6,A))', advance='no') vtk_cd_center(1)-l_force*v_l_t(1), ' ', &
                                             vtk_cd_center(2)-l_force*v_l_t(2), ' ', &
                                             vtk_cd_center(3), ' '
    !  --- Third vertex.
    write(121, '(3(E13.6,A))', advance='no') vtk_an_center(1)-l_force*v_l_t(1), ' ', &
                                             vtk_an_center(2)-l_force*v_l_t(2), ' ', &
                                             vtk_an_center(3), ' '
    !  --- 4th vertex.
    write(121, '(3(E13.6,A))') vtk_an_center(1)+l_force*v_l_t(1), ' ', &
                               vtk_an_center(2)+l_force*v_l_t(2), ' ', &
                               vtk_an_center(3), ' '

  end do

  ! Writing the conectivity
  write(121,'(A)', advance='no') 'POLYGONS '
  ! 6 faces for each parallepiped
  write(121, '(2(I8,A))') f_counter, ' ' , (f_counter*5)

  curr_l_faces = 0

  do i=1, nb_ligneCONTACT_BODY
    ! Each one with 1 face
    if (TAB_CONTACT_BODY(i)%nature == 'CLJCx') cycle

    !!!!! 1st (1-2-3-4)
    write(121, '(I1,A)', advance= 'no') 4, ' '
    write(121, '(I7,A)', advance = 'no') curr_l_faces, ' '
    write(121, '(I7,A)', advance = 'no') curr_l_faces +1, ' '
    write(121, '(I7,A)', advance = 'no') curr_l_faces +2, ' '
    write(121, '(I7,A)') curr_l_faces +3, ' '

    curr_l_faces = curr_l_faces + 4
  end do


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!! FORCES FIELDS !!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! How many fields there will be?
  n_l_fields = 3                          ! They are: Type force (Compression - tension), RN, RT
                                          ! ...

  ! A blank space
  write(121,'(A)') ''
  ! The cells begin
  write(121,'(A)', advance='no') 'CELL_DATA'

  ! Writing the number of data by field. It corresponds to the same number of faces
  write(121, '(2(I8,A))') f_counter
  write(121, '(A,I4)')  'FIELD FieldData', n_l_fields

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Type force (1 Compression, -1 Tension)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Naming the field, the dimension of the data, and number of lines (as before FIELD), and data type
  write(121,'(A,I1,I8,A)') 'Type ', 1, f_counter, ' float'
  k = 0
  do i=1, nb_ligneCONTACT_BODY
    if (TAB_CONTACT_BODY(i)%nature == 'CLJCx') cycle
    ! The particles ids
    l_cdt = TAB_CONTACT_BODY(i)%icdent
    l_ant = TAB_CONTACT_BODY(i)%ianent

    ! The center of each particle
    vtk_cd_center(:) = TAB_MESH(l_cdt)%center(:)
    vtk_an_center(:) = TAB_MESH(l_ant)%center(:)

    ! The branch
    Lik(:) = TAB_MESH(l_cdt)%center(:)-TAB_MESH(l_ant)%center(:)

    !Normalizing
    Lik(:) = Lik(:)/(Lik(1)**2 + Lik(2)**2)**0.5

    ! The normal and tangential vectors... are replaced following the branch
    v_l_normal(:) = Lik(:)
    v_l_t(1) = Lik(2)
    v_l_t(2) = -Lik(1)

    l_rn = TAB_CONTACT_BODY(i)%F_tot(1)*v_l_normal(1) + TAB_CONTACT_BODY(i)%F_tot(2)*v_l_normal(2)
    l_rt = abs(TAB_CONTACT_BODY(i)%F_tot(1)*v_l_t(1) + TAB_CONTACT_BODY(i)%F_tot(2)*v_l_t(2))

    k=k+1
    if(k .le. 9) then
      if (l_rn .lt. 0) then
        write(121, '(I3)', advance='no') -1
      else
        write(121, '(I3)', advance='no') 1
      end if
    else
      if (l_rn .lt. 0) then
        write(121, '(I3)') -1
      else
        write(121, '(I3)') 1
      end if
      k=0
    end if
  end do
  ! And jump
  write(121, '(A)') ' '

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Normal Force
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Naming the field, the dimension of the data, and number of lines (as before FIELD), and data type
  write(121,'(A,I1,I8,A)') 'RN ', 1, f_counter, ' float'
  k = 0
  do i=1, nb_ligneCONTACT_BODY
    if (TAB_CONTACT_BODY(i)%nature == 'CLJCx') cycle
    ! The particles ids
    l_cdt = TAB_CONTACT_BODY(i)%icdent
    l_ant = TAB_CONTACT_BODY(i)%ianent

    ! The center of each particle
    vtk_cd_center(:) = TAB_MESH(l_cdt)%center(:)
    vtk_an_center(:) = TAB_MESH(l_ant)%center(:)

    ! The branch
    Lik(:) = TAB_MESH(l_cdt)%center(:)-TAB_MESH(l_ant)%center(:)

    !Normalizing
    Lik(:) = Lik(:)/(Lik(1)**2 + Lik(2)**2)**0.5

    ! The normal and tangential vectors... are replaced following the branch
    v_l_normal(:) = Lik(:)
    v_l_t(1) = Lik(2)
    v_l_t(2) = -Lik(1)

    l_rn = TAB_CONTACT_BODY(i)%F_tot(1)*v_l_normal(1) + TAB_CONTACT_BODY(i)%F_tot(2)*v_l_normal(2)
    l_rt = abs(TAB_CONTACT_BODY(i)%F_tot(1)*v_l_t(1) + TAB_CONTACT_BODY(i)%F_tot(2)*v_l_t(2))

    k=k+1

    if(k .le. 9) then
      write(121, '(E13.6)', advance='no') l_rn
    else
      write(121, '(E13.6)') l_rn
      k=0
    end if
  end do
  ! And jump
  write(121, '(A)') ' '

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Tangential force
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (along the perpendicular orientation of the branch)

  ! Naming the field, the dimension of the data, and number of lines (as before FIELD), and data type
  write(121,'(A,I1,I6,A)') 'RT ', 1, f_counter, ' float'
  k = 0
  do i=1, nb_ligneCONTACT_BODY
    if (TAB_CONTACT_BODY(i)%nature == 'CLJCx') cycle
    ! The particles ids
    l_cdt = TAB_CONTACT_BODY(i)%icdent
    l_ant = TAB_CONTACT_BODY(i)%ianent

    ! The center of each particle
    vtk_cd_center(:) = TAB_MESH(l_cdt)%center(:)
    vtk_an_center(:) = TAB_MESH(l_ant)%center(:)

    ! The branch
    Lik(:) = TAB_MESH(l_cdt)%center(:)-TAB_MESH(l_ant)%center(:)

    !Normalizing
    Lik(:) = Lik(:)/(Lik(1)**2 + Lik(2)**2)**0.5

    ! The normal and tangential vectors... are replaced following the branch
    v_l_normal(:) = Lik(:)
    v_l_t(1) = Lik(2)
    v_l_t(2) = -Lik(1)

    l_rn = TAB_CONTACT_BODY(i)%F_tot(1)*v_l_normal(1) + TAB_CONTACT_BODY(i)%F_tot(2)*v_l_normal(2)
    l_rt = abs(TAB_CONTACT_BODY(i)%F_tot(1)*v_l_t(1) + TAB_CONTACT_BODY(i)%F_tot(2)*v_l_t(2))

    k=k+1
    if(k .le. 9) then
      write(121, '(E13.6)', advance='no') l_rt
    else
      write(121, '(E13.6)') l_rt
        k=0
    end if
  end do
  ! And jump
  write(121, '(A)') ' '

  close(121)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!  CONTACT POINTS   !!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  vtk_inter = './POSTPRO/VTK/ctci_' // vtk_counter // '.vtk'

  open(unit=121, file=vtk_inter, status='replace')
  write(121,'(A)') '# vtk DataFile Version 3.0'
  write(121,'(A,I6)') 'RIGID ', compteur_clout
  write(121,'(A)') 'ASCII'
  write(121,'(A)') 'DATASET POLYDATA'

  ! Finding the number of vertices
  ! In this case each disc will have n_l_sides sides. The vertices are written counterclockwise
  n_l_sides = 6
  n_l_vertices = nb_ligneCONTACT * n_l_sides

  write(v_n_vertices, '(I8)') n_l_vertices

  vtk_n_points = 'POINTS ' // v_n_vertices // ' float'
  write(121,'(A)') vtk_n_points

  ! Writing the coordinates of the vertices
  ! 3 vertices per row
  !

  ! Defining a fix size for contact points at nodes
  rad_ctci = l_ave_rad*0.05
  k=0
  do i=1, nb_ligneCONTACT
    !write(121,*) 'Particle', i
    do j=1, n_l_sides
      cur_angle_rad = (pi/(n_l_sides/2))*(j-1)
      cur_disc_vertex(1) = TAB_CONTACT(i)%coor_ctc(1) + rad_ctci*cos(cur_angle_rad)
      cur_disc_vertex(2) = TAB_CONTACT(i)%coor_ctc(2) + rad_ctci*sin(cur_angle_rad)
      k = k + 1

      if(k .lt. 3) then
        write(121,'(E13.6,A)', advance = 'no') cur_disc_vertex(1), ' '
        write(121,'(E13.6,A)', advance = 'no') cur_disc_vertex(2), ' '
        write(121,'(E13.6,A)', advance = 'no') 0.D0 , ' '
      else
        write(121,'(E13.6,A)', advance = 'no') cur_disc_vertex(1), ' '
        write(121,'(E13.6,A)', advance = 'no') cur_disc_vertex(2), ' '
        write(121,'(E13.6,A)') 0.D0, ' '
        k = 0
      end if
    end do
  end do

  write(121, '(A)') ' '

  ! Writing the conectivity between vertices
  write(121,'(A)', advance='no') 'POLYGONS '
  write(121, '(2(I8,A))') nb_ligneCONTACT, ' ' , nb_ligneCONTACT*(n_l_sides+1)

  curr_l_faces = 0

  do i=1, nb_ligneCONTACT
    !write(121,*) 'Particle', i
    ! Write the number of vertices

    write(121, '(I2,A)', advance= 'no') n_l_sides, ' '

    ! Write its consecutive conectivity
    do j=1, n_l_sides
      ! ... writing precisely
      if (j .lt. n_l_sides) then
        write(121, '(I7,A)', advance = 'no') curr_l_faces - 1 + j, ' '
      else
        write(121, '(I7,A)') curr_l_faces - 1 + j, ' '
      end if
    end do
    curr_l_faces = curr_l_faces + n_l_sides
  end do

  ! And jump
  write(121, '(A)') ' '


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!! INDIVIDUAL CONTACT POINTS FIELDS !!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! How many fields there will be?
  n_l_fields = 4                          ! They are: Type force (Active - Float), RN, RT, M_fr
                                          ! ...

  ! A blank space
  write(121,'(A)') ''
  ! The cells begin
  write(121,'(A)', advance='no') 'CELL_DATA'

  ! Writing the number of data by field. It corresponds to the same number of faces
  write(121, '(2(I8,A))') nb_ligneCONTACT
  write(121, '(A,I4)')  'FIELD FieldData', n_l_fields

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Type force (1 Active, -1 Float)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Naming the field, the dimension of the data, and number of lines (as before FIELD), and data type
  write(121,'(A,I1,I8,A)') 'Type ', 1, nb_ligneCONTACT, ' float'
  k = 0
  do i=1, nb_ligneCONTACT

    l_rn = TAB_CONTACT(i)%rn
    l_rt = abs(TAB_CONTACT(i)%rt)

    k=k+1
    if(k .le. 9) then
      if (l_rn .lt. 1.E-9) then
        write(121, '(I3)', advance='no') -1
      else
        write(121, '(I3)', advance='no') 1
      end if
    else
      if (l_rn .lt. 1.E-9) then
         write(121, '(I3)') -1
      else
        write(121, '(I3)') 1
      end if
      k=0
    end if
  end do
  ! And jump
  write(121, '(A)') ' '

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Normal Force
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Naming the field, the dimension of the data, and number of lines (as before FIELD), and data type
  write(121,'(A,I1,I8,A)') 'RN ', 1, nb_ligneCONTACT, ' float'
  k = 0
  do i=1, nb_ligneCONTACT
    l_rn = TAB_CONTACT(i)%rn
    l_rt = abs(TAB_CONTACT(i)%rt)

    k=k+1

    if(k .le. 9) then
      write(121, '(E13.6)', advance='no') l_rn
    else
      write(121, '(E13.6)') l_rn
      k=0
    end if
  end do
  ! And jump
  write(121, '(A)') ' '

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Tangential force
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (along the perpendicular orientation of the branch)

  ! Naming the field, the dimension of the data, and number of lines (as before FIELD), and data type
  write(121,'(A,I1,I8,A)') 'RT ', 1, nb_ligneCONTACT, ' float'
  k = 0
  do i=1, nb_ligneCONTACT
    l_rn = TAB_CONTACT(i)%rn
    l_rt = abs(TAB_CONTACT(i)%rt)

    k=k+1
    if(k .le. 9) then
      write(121, '(E13.6)', advance='no') l_rt
    else
      write(121, '(E13.6)') l_rt
        k=0
    end if
  end do
  ! And jump
  write(121, '(A)') ' '

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Friction mobilization
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (along the perpendicular orientation of the branch)

  ! Naming the field, the dimension of the data, and number of lines (as before FIELD), and data type
  write(121,'(A,I1,I8,A)') 'Mfr', 1, nb_ligneCONTACT, ' float'
  k = 0
  do i=1, nb_ligneCONTACT
    l_rn = TAB_CONTACT(i)%rn
    l_rt = abs(TAB_CONTACT(i)%rt)

    k=k+1
    if(k .le. 9) then
      write(121, '(E13.6)', advance='no') l_rt/(fric_part*l_rn)
    else
      write(121, '(E13.6)') l_rt/(fric_part*l_rn)
        k=0
    end if
  end do
  ! And jump
  write(121, '(A)') ' '

  close(121)

  print*, 'Drawing in vtk       ---> Ok!'
  !stop

end subroutine draw


!==============================================================================
! Computing the deformation and stress tensor for each element of the mesh
!==============================================================================
! Careful... it is written for triangular meshing
subroutine mesh_stresses

  implicit none

  integer                         :: i, j, k
  real*8                          :: lambda_l, mu_l
  real*8                          :: det_F, I_1, I_2, I_3, tr_c2
  real*8, dimension(2)            :: node_1, node_2, node_3
  real*8, dimension(2)            :: node_1_ref, node_2_ref, node_3_ref
  real*8, dimension(2)            :: disp_1, disp_2, disp_3
  real*8, dimension(2,2)          :: F, F_t, C_local, C_local_sqr, C_inv

  ! Initializing all the tensor
  !TAB_MESH(:)%mesh_connect(:)%jacobi(:,:) = 0.0
  !TAB_MESH(:)%mesh_connect(:)%inv_jacobi(:,:) = 0.0
  !TAB_MESH(:)%mesh_connect(:)%def_grad(:,:) = 0.0
  !TAB_MESH(:)%mesh_connect(:)%cg_left(:,:) = 0.0
  !TAB_MESH(:)%mesh_connect(:)%cg_right(:,:) = 0.0

  ! Initializing the determinant of the jacobian of the transformation
  !TAB_MESH(:)%mesh_connect(:)%det_j = 0.0

  ! Initializing the area of each element in the current configuration
  !TAB_MESH(i)%mesh_connect(:)%surf_element = 0.0

  ! 1. Computing the jacobian for each finite element
  ! This is the spacial gradient of the displacement of for each element written
  ! Big advantage! We have triangles for the meshing. Then, the gradient is
  ! constant through the finite element

  ! Computing the Lame constants
  lambda_l = (poisson_ratio*E_young)/((1+poisson_ratio)*(1-2*poisson_ratio))
  mu_l     = E_young/(2*(1+poisson_ratio))

  ! So... for each particle
  do i=1, n_def_particles
    ! For each element
    !print*, 'part: ', i
    do j=1, TAB_MESH(i)%nb_mesh
      !print*, 'mesh: ', j
      ! Initializing
      TAB_MESH(i)%mesh_connect(j)%cg_right(:,:) = 0.0

      TAB_MESH(i)%mesh_connect(j)%cg_left(:,:)  = 0.0

      TAB_MESH(i)%mesh_connect(j)%pk_2nd(:,:)   = 0.0

      TAB_MESH(i)%mesh_connect(j)%def_grad(:,:) = 0.0

      TAB_MESH(i)%mesh_connect(j)%dif_disp(:,:) = 0.0

      TAB_MESH(i)%mesh_connect(j)%cauchy_s(:,:) = 0.0

      ! Nodes coordinates
      node_1(1) = TAB_MESH(i)%node(TAB_MESH(i)%mesh_connect(j)%connectivity(1),1)
      node_1(2) = TAB_MESH(i)%node(TAB_MESH(i)%mesh_connect(j)%connectivity(1),2)

      node_2(1) = TAB_MESH(i)%node(TAB_MESH(i)%mesh_connect(j)%connectivity(2),1)
      node_2(2) = TAB_MESH(i)%node(TAB_MESH(i)%mesh_connect(j)%connectivity(2),2)

      node_3(1) = TAB_MESH(i)%node(TAB_MESH(i)%mesh_connect(j)%connectivity(3),1)
      node_3(2) = TAB_MESH(i)%node(TAB_MESH(i)%mesh_connect(j)%connectivity(3),2)

      ! Nodes coordinates in ref
      node_1_ref(1) = TAB_MESH(i)%node_ref(TAB_MESH(i)%mesh_connect(j)%connectivity(1),1)
      node_1_ref(2) = TAB_MESH(i)%node_ref(TAB_MESH(i)%mesh_connect(j)%connectivity(1),2)

      node_2_ref(1) = TAB_MESH(i)%node_ref(TAB_MESH(i)%mesh_connect(j)%connectivity(2),1)
      node_2_ref(2) = TAB_MESH(i)%node_ref(TAB_MESH(i)%mesh_connect(j)%connectivity(2),2)

      node_3_ref(1) = TAB_MESH(i)%node_ref(TAB_MESH(i)%mesh_connect(j)%connectivity(3),1)
      node_3_ref(2) = TAB_MESH(i)%node_ref(TAB_MESH(i)%mesh_connect(j)%connectivity(3),2)

      ! The total displacement of each node
      disp_1(1) = node_1(1) - node_1_ref(1)
      disp_1(2) = node_1(2) - node_1_ref(2)

      disp_2(1) = node_2(1) - node_2_ref(1)
      disp_2(2) = node_2(2) - node_2_ref(2)

      disp_3(1) = node_3(1) - node_3_ref(1)
      disp_3(2) = node_3(2) - node_3_ref(2)

      ! Computing the area of the element
      !TAB_MESH(i)%mesh_connect(j)%surf_element = node_1(1)*node_2(2) - node_1(2)*node_2(1) + &
      !                                           node_2(1)*node_3(2) - node_2(2)*node_3(1) + &
      !                                           node_3(1)*node_1(2) - node_3(2)*node_1(1)
      TAB_MESH(i)%mesh_connect(j)%surf_element = node_1_ref(1)*node_2_ref(2) - node_1_ref(2)*node_2_ref(1) + &
                                                 node_2_ref(1)*node_3_ref(2) - node_2_ref(2)*node_3_ref(1) + &
                                                 node_3_ref(1)*node_1_ref(2) - node_3_ref(2)*node_1_ref(1)
      TAB_MESH(i)%mesh_connect(j)%surf_element = abs(TAB_MESH(i)%mesh_connect(j)%surf_element*0.5)

      ! The jacobian for this element is then:
      TAB_MESH(i)%mesh_connect(j)%det_j = 2*TAB_MESH(i)%mesh_connect(j)%surf_element

      ! This is the jacobian matrix of the transformation for this triangular element (J)
      ! This equations are "easily" found after working in the natural coordinates of the element
      ! For a triangular element, the gradient of displacements are constant for all the element
      ! The jacobian matrix of the tranformation is defined as follows:
      !
      ! J = |d_x/d_zeta d_y/d_zeta|
      !     |d_x/d_eta  d_y/d_eta |
      ! Being \zeta et \eta ... and the base of the natural coordinates frame
      ! The functions x and y (the displacement functions in the local frame) are
      ! built after the shape functions for this finite element
      !TAB_MESH(i)%mesh_connect(j)%jacobi(1,1) = (node_2(1)-node_1(1))/2
      !TAB_MESH(i)%mesh_connect(j)%jacobi(1,2) = (node_3(1)-node_1(1))/2
      !TAB_MESH(i)%mesh_connect(j)%jacobi(2,1) = (node_2(2)-node_1(2))/2
      !TAB_MESH(i)%mesh_connect(j)%jacobi(2,2) = (node_3(2)-node_1(2))/2
      !TAB_MESH(i)%mesh_connect(j)%jacobi(1,1) = (node_2_ref(1)-node_1_ref(1))/2
      !TAB_MESH(i)%mesh_connect(j)%jacobi(1,2) = (node_3_ref(1)-node_1_ref(1))/2
      !TAB_MESH(i)%mesh_connect(j)%jacobi(2,1) = (node_2_ref(2)-node_1_ref(2))/2
      !TAB_MESH(i)%mesh_connect(j)%jacobi(2,2) = (node_3_ref(2)-node_1_ref(2))/2

      ! Matrix left
      TAB_MESH(i)%mesh_connect(j)%dif_disp(1,1) = disp_2(1) - disp_1(1)
      TAB_MESH(i)%mesh_connect(j)%dif_disp(1,2) = disp_3(1) - disp_1(1)
      TAB_MESH(i)%mesh_connect(j)%dif_disp(2,1) = disp_2(2) - disp_1(2)
      TAB_MESH(i)%mesh_connect(j)%dif_disp(2,2) = disp_3(2) - disp_1(2)

      ! Computing the inverse of the jacobian matrix of the transformation(J^1)
      ! Having J and J^-1 allows to transform between the natural coordinates
      ! and the cartesian coordinates
      TAB_MESH(i)%mesh_connect(j)%inv_jacobi(1,1) = (node_3_ref(2)-node_1_ref(2))!/2! TAB_MESH(i)%mesh_connect(j)%jacobi(2,2)
      TAB_MESH(i)%mesh_connect(j)%inv_jacobi(1,2) = (node_1_ref(1)-node_3_ref(1))!/2!-TAB_MESH(i)%mesh_connect(j)%jacobi(1,2)
      TAB_MESH(i)%mesh_connect(j)%inv_jacobi(2,1) = (node_1_ref(2)-node_2_ref(2))!/2!-TAB_MESH(i)%mesh_connect(j)%jacobi(2,1)
      TAB_MESH(i)%mesh_connect(j)%inv_jacobi(2,2) = (node_2_ref(1)-node_1_ref(1))!/2! TAB_MESH(i)%mesh_connect(j)%jacobi(1,1)

      TAB_MESH(i)%mesh_connect(j)%inv_jacobi(:,:) = TAB_MESH(i)%mesh_connect(j)%inv_jacobi(:,:)/ &
                                                    (2*TAB_MESH(i)%mesh_connect(j)%surf_element)
      ! Finally
      !TAB_MESH(i)%mesh_connect(j)%inv_jacobi(:,:) = TAB_MESH(i)%mesh_connect(j)%inv_jacobi(:,:)/ &
      !                                              TAB_MESH(i)%mesh_connect(j)%surf_element !TAB_MESH(i)%mesh_connect(j)%det_j

      !print*, TAB_MESH(i)%mesh_connect(j)%inv_jacobi

      ! Now we can compute the deformation gradient. The jacobian matrix helps there
      ! The deformation gradient is defined as:
      !
      ! F=I + \grad U = I + 1/2|du_x/d\zeta du_x/d\eta||d_zeta/dX d_zeta/dY|=I + 0.5*|du_x/d\zeta du_x/d\eta|*(Inv_J)
      !                        |du_y/d\zeta du_x/d\eta||d\eta /dX d_eta/dY |         |du_y/d\zeta du_x/d\eta|
      ! Being u and v the displacement functions in the eulerian frame for this finite element

      !TAB_MESH(i)%mesh_connect(j)%def_grad(1,1) = (disp_2(1) - disp_1(1))*TAB_MESH(i)%mesh_connect(j)%inv_jacobi(1,1)+ &
      !                                            (disp_3(1) - disp_1(1))*TAB_MESH(i)%mesh_connect(j)%inv_jacobi(2,1)
!
      !TAB_MESH(i)%mesh_connect(j)%def_grad(1,2) = (disp_2(1) - disp_1(1))*TAB_MESH(i)%mesh_connect(j)%inv_jacobi(1,2)+ &
      !                                            (disp_3(1) - disp_1(1))*TAB_MESH(i)%mesh_connect(j)%inv_jacobi(2,2)
!
      !TAB_MESH(i)%mesh_connect(j)%def_grad(2,1) = (disp_2(2) - disp_1(2))*TAB_MESH(i)%mesh_connect(j)%inv_jacobi(1,1)+ &
      !                                            (disp_3(2) - disp_1(2))*TAB_MESH(i)%mesh_connect(j)%inv_jacobi(2,1)
!
      !TAB_MESH(i)%mesh_connect(j)%def_grad(2,2) = (disp_2(2) - disp_1(2))*TAB_MESH(i)%mesh_connect(j)%inv_jacobi(1,2)+ &
!                                                  (disp_3(2) - disp_1(2))*TAB_MESH(i)%mesh_connect(j)%inv_jacobi(2,2)

      TAB_MESH(i)%mesh_connect(j)%def_grad = matmul(TAB_MESH(i)%mesh_connect(j)%dif_disp, &
                                                    TAB_MESH(i)%mesh_connect(j)%inv_jacobi)
      ! One half
      !TAB_MESH(i)%mesh_connect(j)%def_grad = TAB_MESH(i)%mesh_connect(j)%def_grad*0.5 &
      !                                       /(2*TAB_MESH(i)%mesh_connect(j)%surf_element)

      ! Avoiding parasite behavior
      if (abs(TAB_MESH(i)%mesh_connect(j)%def_grad(1,1)) .lt. 1e-9) then
        TAB_MESH(i)%mesh_connect(j)%def_grad(1,1) = 0.0
      end if
      if (abs(TAB_MESH(i)%mesh_connect(j)%def_grad(1,2)) .lt. 1e-9) then
        TAB_MESH(i)%mesh_connect(j)%def_grad(1,2) = 0.0
      end if
      if (abs(TAB_MESH(i)%mesh_connect(j)%def_grad(2,1)) .lt. 1e-9) then
        TAB_MESH(i)%mesh_connect(j)%def_grad(2,1) = 0.0
      end if
      if (abs(TAB_MESH(i)%mesh_connect(j)%def_grad(2,2)) .lt. 1e-9) then
        TAB_MESH(i)%mesh_connect(j)%def_grad(2,2) = 0.0
      end if

      ! Adding the identity
      TAB_MESH(i)%mesh_connect(j)%def_grad(1,1) = TAB_MESH(i)%mesh_connect(j)%def_grad(1,1)+1.0
      TAB_MESH(i)%mesh_connect(j)%def_grad(2,2) = TAB_MESH(i)%mesh_connect(j)%def_grad(2,2)+1.0

      ! For practicality
      F(1,1) = TAB_MESH(i)%mesh_connect(j)%def_grad(1,1)
      F(1,2) = TAB_MESH(i)%mesh_connect(j)%def_grad(1,2)
      F(2,1) = TAB_MESH(i)%mesh_connect(j)%def_grad(2,1)
      F(2,2) = TAB_MESH(i)%mesh_connect(j)%def_grad(2,2)

      !print*, 'Mesh stress'
      !print*, i
      !print*, j
      !print*, F

      ! Computing the transpose of F
      F_t(1,1) = F(1,1)
      F_t(1,2) = F(2,1)
      F_t(2,1) = F(1,2)
      F_t(2,2) = F(2,2)

      ! Now we can compute the right and left Cauchy-Green deformation tensor
      ! Right (C)
      !TAB_MESH(i)%mesh_connect(j)%cg_right(1,1) = F_t(1,1)*F(1,1)+F_t(1,2)*F(2,1)
      !TAB_MESH(i)%mesh_connect(j)%cg_right(1,2) = F_t(1,1)*F(1,2)+F_t(1,2)*F(2,2)
      !TAB_MESH(i)%mesh_connect(j)%cg_right(2,1) = F_t(2,1)*F(1,1)+F_t(2,2)*F(2,1)
      !TAB_MESH(i)%mesh_connect(j)%cg_right(2,2) = F_t(2,1)*F(1,2)+F_t(2,2)*F(2,2)

      TAB_MESH(i)%mesh_connect(j)%cg_right = matmul(F_t,F)

      ! Left (B)
      !TAB_MESH(i)%mesh_connect(j)%cg_left(1,1) =  F(1,1)*F_t(1,1)+F(1,2)*F_t(2,1)
      !TAB_MESH(i)%mesh_connect(j)%cg_left(1,2) =  F(1,1)*F_t(1,2)+F(1,2)*F_t(2,2)
      !TAB_MESH(i)%mesh_connect(j)%cg_left(2,1) =  F(2,1)*F_t(1,1)+F(2,2)*F_t(2,1)
      !TAB_MESH(i)%mesh_connect(j)%cg_left(2,2) =  F(2,1)*F_t(1,2)+F(2,2)*F_t(2,2)

      TAB_MESH(i)%mesh_connect(j)%cg_left = matmul(F,F_t)

      ! The elastic energy density of a hyperelastic neo-Hookean material is written as follows:
      ! W = (1/2)\lambda *(ln(J))**2 - \mu ln(J) + (1/2)\mu (I_1 - 3)
      ! Being \lambda and \mu the Lame constants
      ! J the determinant of the deformation gradient = sqrt(determinant of C)
      ! and I_1 the first invariant of C => tr(C)

      ! For simplicity
      C_local(1,1) = TAB_MESH(i)%mesh_connect(j)%cg_right(1,1)
      C_local(1,2) = TAB_MESH(i)%mesh_connect(j)%cg_right(1,2)
      C_local(2,1) = TAB_MESH(i)%mesh_connect(j)%cg_right(2,1)
      C_local(2,2) = TAB_MESH(i)%mesh_connect(j)%cg_right(2,2)

      ! Computing the 1st invariant of C
      I_1 = TAB_MESH(i)%mesh_connect(j)%cg_right(1,1) + TAB_MESH(i)%mesh_connect(j)%cg_right(2,2)

      ! Computing the 2nd
      ! First we need to have C**2
      !C_local_sqr(1,1) = C_local(1,1)*C_local(1,1)+C_local(1,2)*C_local(2,1)
      !C_local_sqr(1,2) = C_local(1,1)*C_local(1,2)+C_local(1,2)*C_local(2,2)
      !C_local_sqr(2,1) = C_local(2,1)*C_local(1,1)+C_local(2,2)*C_local(2,1)
      !C_local_sqr(2,2) = C_local(2,1)*C_local(1,2)+C_local(2,2)*C_local(2,2)

      C_local_sqr = matmul(C_local,C_local)
      ! Finding its trace
      tr_c2 = C_local_sqr(1,1) + C_local_sqr(2,2)

      I_2 = 0.5*((I_1**2) - tr_c2)

      ! Computing the third invariants
      ! Determinant of C
      I_3 = C_local(1,1)*C_local(2,2) - C_local(1,2)*C_local(2,1)

      ! Now we can compute the determinant of F
      det_F = sqrt(I_3)

      ! We need to compute the inverse of C ... this is also known as Finger tensor
      C_inv(1,1) = C_local(2,2)
      C_inv(1,2) = -C_local(1,2)
      C_inv(2,1) = -C_local(2,1)
      C_inv(2,2) = C_local(1,1)

      C_inv = C_inv/(det_F**2)

      ! Now we can compute the second Piola-Kirchhoff stress tensor for this element
      TAB_MESH(i)%mesh_connect(j)%pk_2nd(1,1) = (lambda_l*log(det_F**2)-mu_l)*C_inv(1,1) + mu_l
      TAB_MESH(i)%mesh_connect(j)%pk_2nd(1,2) = (lambda_l*log(det_F**2)-mu_l)*C_inv(1,2)
      TAB_MESH(i)%mesh_connect(j)%pk_2nd(2,1) = (lambda_l*log(det_F**2)-mu_l)*C_inv(2,1)
      TAB_MESH(i)%mesh_connect(j)%pk_2nd(2,2) = (lambda_l*log(det_F**2)-mu_l)*C_inv(2,2) + mu_l

      TAB_MESH(i)%mesh_connect(j)%green_E(1,1) = 0.5*(C_local(1,1)-1)
      TAB_MESH(i)%mesh_connect(j)%green_E(1,2) = 0.5*(C_local(1,2))
      TAB_MESH(i)%mesh_connect(j)%green_E(2,1) = 0.5*(C_local(2,1))
      TAB_MESH(i)%mesh_connect(j)%green_E(2,2) = 0.5*(C_local(2,2)-1)


      ! Computing the equivalent Cauchy stress tensor
      TAB_MESH(i)%mesh_connect(j)%cauchy_s = matmul(TAB_MESH(i)%mesh_connect(j)%pk_2nd,F_t)
      TAB_MESH(i)%mesh_connect(j)%cauchy_s = matmul(F,TAB_MESH(i)%mesh_connect(j)%cauchy_s)
      TAB_MESH(i)%mesh_connect(j)%cauchy_s = TAB_MESH(i)%mesh_connect(j)%cauchy_s* &
                                             (1/(F(1,1)*F(2,2) - F(1,2)*F(2,1)))

      if (.false.) then
        ! The nodes of this finite element in the current configuration
        print*, 'connect'
        print*, TAB_MESH(i)%mesh_connect(j)%connectivity(1)
        print*, TAB_MESH(i)%mesh_connect(j)%connectivity(2)
        print*, TAB_MESH(i)%mesh_connect(j)%connectivity(3)

        print*, 'Coordi'
        print*, node_1
        print*, node_2
        print*, node_3

        print*, 'Disp'
        print*, disp_1
        print*, disp_2
        print*, disp_3

        print*, 'area'
        print*, TAB_MESH(i)%mesh_connect(j)%surf_element

        print*, 'FFFFFFFFFF'
        print*, F

        print*, 'RIGHT'
        print*, C_local

        print*, 'Cinvinvinivinivinv'
        print*, C_inv

        print*, 'Piola-Kirchhoff'
        print*, TAB_MESH(i)%mesh_connect(j)%pk_2nd(1,1), TAB_MESH(i)%mesh_connect(j)%pk_2nd(1,2)
        print*, TAB_MESH(i)%mesh_connect(j)%pk_2nd(2,1), TAB_MESH(i)%mesh_connect(j)%pk_2nd(2,2)

      end if
    end do
  end do

  ! But we also want the stress and deformation tensor at the nodes
  do i=1, n_def_particles
    !print*, 'Part: ', i
    !if (first_over_all) then
    if (allocated(TAB_MESH(i)%node_info)) deallocate(TAB_MESH(i)%node_info)
    allocate (TAB_MESH(i)%node_info(TAB_MESH(i)%nb_node))
    !end if
    do j=1, TAB_MESH(i)%nb_node
      !print*, 'Node: ', j

      ! Initializing
      TAB_MESH(i)%node_info(j)%cg_right(:,:) = 0.0

      TAB_MESH(i)%node_info(j)%cg_left(:,:)  = 0.0

      TAB_MESH(i)%node_info(j)%pk_2nd(:,:)   = 0.0

      TAB_MESH(i)%node_info(j)%def_grad(:,:) = 0.0

      TAB_MESH(i)%node_info(j)%green_E(:,:) = 0.0

      TAB_MESH(i)%node_info(j)%cauchy_s(:,:) = 0.0

      TAB_MESH(i)%node_info(j)%common_to   = 0

      ! Looking for the node in the connectivity list
      do k=1, TAB_MESH(i)%nb_mesh
        if (j == TAB_MESH(i)%mesh_connect(k)%connectivity(1) .or. &
            j == TAB_MESH(i)%mesh_connect(k)%connectivity(2) .or. &
            j == TAB_MESH(i)%mesh_connect(k)%connectivity(3)) then

          TAB_MESH(i)%node_info(j)%cg_right(:,:) = TAB_MESH(i)%node_info(j)%cg_right(:,:) + &
                                                   TAB_MESH(i)%mesh_connect(k)%cg_right(:,:)

          TAB_MESH(i)%node_info(j)%cg_left(:,:) = TAB_MESH(i)%node_info(j)%cg_left(:,:) + &
                                                  TAB_MESH(i)%mesh_connect(k)%cg_left(:,:)

          TAB_MESH(i)%node_info(j)%pk_2nd(:,:) = TAB_MESH(i)%node_info(j)%pk_2nd(:,:) + &
                                                 TAB_MESH(i)%mesh_connect(k)%pk_2nd(:,:)

          TAB_MESH(i)%node_info(j)%def_grad(:,:) = TAB_MESH(i)%node_info(j)%def_grad(:,:) + &
                                                   TAB_MESH(i)%mesh_connect(k)%def_grad(:,:)

          TAB_MESH(i)%node_info(j)%green_E(:,:) = TAB_MESH(i)%node_info(j)%green_E(:,:) + &
                                                   TAB_MESH(i)%mesh_connect(k)%green_E(:,:)

          TAB_MESH(i)%node_info(j)%cauchy_s(:,:) = TAB_MESH(i)%node_info(j)%cauchy_s(:,:) + &
                                                   TAB_MESH(i)%mesh_connect(k)%cauchy_s(:,:)

          TAB_MESH(i)%node_info(j)%common_to = TAB_MESH(i)%node_info(j)%common_to + 1
        end if
      end do

      !print*, TAB_MESH(i)%node_info(j)%common_to
      ! This is the average to the common elements to this node
      TAB_MESH(i)%node_info(j)%cg_right(:,:) = TAB_MESH(i)%node_info(j)%cg_right(:,:)/TAB_MESH(i)%node_info(j)%common_to

      TAB_MESH(i)%node_info(j)%cg_left(:,:)  = TAB_MESH(i)%node_info(j)%cg_left(:,:)/TAB_MESH(i)%node_info(j)%common_to

      TAB_MESH(i)%node_info(j)%pk_2nd(:,:)   = TAB_MESH(i)%node_info(j)%pk_2nd(:,:)/TAB_MESH(i)%node_info(j)%common_to

      TAB_MESH(i)%node_info(j)%def_grad(:,:) = TAB_MESH(i)%node_info(j)%def_grad(:,:)/TAB_MESH(i)%node_info(j)%common_to

      TAB_MESH(i)%node_info(j)%green_E(:,:) = TAB_MESH(i)%node_info(j)%green_E(:,:)/TAB_MESH(i)%node_info(j)%common_to

      TAB_MESH(i)%node_info(j)%cauchy_s(:,:) = TAB_MESH(i)%node_info(j)%cauchy_s(:,:)/TAB_MESH(i)%node_info(j)%common_to

      !print*, TAB_MESH(i)%node_info(j)%def_grad
      !stop
    end do
  end do
  print*, 'Computing mesh stresses         ---> Ok!'
end subroutine mesh_stresses

!==============================================================================
! Writing the Green-Lagrange (E) strain tensor
!==============================================================================
subroutine green_lagrange

  implicit none

  integer                         :: i, j
  real*8, dimension(2,2)          :: temp_m
  character(len=32)               :: l_file
  logical                                  :: dir_c

  if (c_mesh_stresses .ne. 1) then
    print*, 'Strains and stresses were not computed yet!!!'
    stop
  end if

  ! Cleaning or creating the folder if necessary
  if (first_over_all) then
    ! Asking if the file already exists
    inquire(file='./POSTPRO/LIST_GL', exist=dir_c)
    if(dir_c) then
      ! Cleaning
      call system('rm ./POSTPRO/LIST_GL/*')
    else
      ! Creating
      call system('mkdir ./POSTPRO/LIST_GL')
    end if
  end if

  ! The name of the file
  l_file = './POSTPRO/LIST_GL/Green_La.     '

  ! Preparing the corresponding number of the file
  if (compteur_clout<10) then
    WRITE(l_file(28:29),'(I1)')   compteur_clout
  else if ( (compteur_clout>=10) .and. (compteur_clout<100) ) then
    WRITE(l_file(28:30),'(I2)')   compteur_clout
  else if ( (compteur_clout>=100).and. (compteur_clout<1000) ) then
    WRITE(l_file(28:31),'(I3)')   compteur_clout
  else if ( (compteur_clout>=1000).and. (compteur_clout<10000) ) then
    WRITE(l_file(28:32),'(I4)')   compteur_clout
  else
    print*, "Not enough file names :: List of Green-Lagrange"
    stop
  end if

  ! Opening the file
  open(unit=123,file=l_file,status='replace')

  write(123,*) '#    E_xx    ', '    E_xy    ', '    E_yx    ', '    E_yy    '

  do i=1, n_def_particles
    temp_m(:,:) = 0.
    do j=1, TAB_MESH(i)%nb_mesh
      temp_m = temp_m + TAB_MESH(i)%mesh_connect(j)%green_E
    end do
    temp_m = temp_m/real(TAB_MESH(i)%nb_mesh)
    write(123, '(4(1X,E12.5))') temp_m(1,1), temp_m(1,2), temp_m(2,1), temp_m(2,2)
  end do

  close(123)

  print*, 'Computing strain from DOF       ---> Ok!'

end subroutine green_lagrange

!==============================================================================
! Computing the average interparticle distance
!==============================================================================
subroutine avg_c_dist

  implicit none

  integer                          :: i, cd, an, non_valid
  real*8                           :: l_avg
  real*8, dimension(3)             :: n_vec_defo

  if (first_over_all) then
    write(124,*) '# Average intercenter distance'
  end if

  l_avg = 0.0
  non_valid = 0
  do i=1, nb_ligneCONTACT_BODY
    if  (TAB_CONTACT_BODY(i)%nature == 'CLJCx') then
      non_valid = non_valid + 1
      cycle
    end if
    cd = TAB_CONTACT_BODY(i)%icdent
    an = TAB_CONTACT_BODY(i)%ianent
    ! Computing the unit branch vector
    n_vec_defo = TAB_MESH(cd)%center - TAB_MESH(an)%center
    ! Normalizing
    n_vec_defo = n_vec_defo/(n_vec_defo(1)**2 + n_vec_defo(2)**2 + n_vec_defo(3)**2)**0.5
    if(TAB_CONTACT_BODY(i)%F_tot(1)*n_vec_defo(1) + &
       TAB_CONTACT_BODY(i)%F_tot(2)*n_vec_defo(2) .lt. 1E-9) then
      non_valid = non_valid + 1
      cycle
    else
      l_avg = l_avg + ((TAB_MESH(cd)%center(1)-TAB_MESH(an)%center(1))**2 + &
                       (TAB_MESH(cd)%center(2)-TAB_MESH(an)%center(2))**2)**0.5
    end if
  end do

  print*, l_avg
  print*, nb_ligneCONTACT_BODY
  print*, non_valid

  if (non_valid == nb_ligneCONTACT_BODY) then
    print*, 'No active contacts'
  else
    l_avg = l_avg / (nb_ligneCONTACT_BODY - non_valid)
    write(124,'(1X,E12.5)') l_avg
  end if

  print*, 'Computing average interp dist   ---> Ok!'

end subroutine avg_c_dist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                         EXTERNAL SUBROUTINES                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE envelope(x, y, n, vertex, nvert, iwk)

  !  Find the vertices (in clockwise order) of a polygon enclosing
  !  the points (x(i), y(i), i=1, ..., n.

  !  On output, vertex(i), i=1, ..., nvert contains the numbers of the vertices.
  !  iwk() is an integer work array which must have dimension at least n
  !  in the calling program.

  !  There is a limit of 100 vertices imposed by the dimension of array next.

  !  Programmer: Alan Miller
  !  Latest revision - 12 September 1987
  !  Fortran 90 version - 8 August 1996

  IMPLICIT NONE
  INTEGER :: n, vertex(n), nvert, iwk(n)
  REAL*8  :: x(n), y(n)

  !       Local variables

  INTEGER :: next(100), i, i1, i2, j, jp1, jp2, i2save, i3, i2next
  REAL*8  :: xmax, xmin, ymax, ymin, dist, dmax, dmin, x1, y1, dx, dy, x2, y2, &
             dx1, dx2, dmax1, dmax2, dy1, dy2, temp, zero = 0.0

  IF (n < 2) RETURN

  !  Choose the points with smallest & largest x- values as the
  !  first two vertices of the polygon.

  IF (x(1) > x(n)) THEN
    vertex(1) = n
    vertex(2) = 1
    xmin = x(n)
    xmax = x(1)
  ELSE
    vertex(1) = 1
    vertex(2) = n
    xmin = x(1)
    xmax = x(n)
  END IF

  DO i = 2, n-1
    temp = x(i)
    IF (temp < xmin) THEN
      vertex(1) = i
      xmin = temp
    ELSE IF (temp > xmax) THEN
      vertex(2) = i
      xmax = temp
    END IF
  END DO

  !       Special case, xmax = xmin.

  IF (xmax == xmin) THEN
    IF (y(1) > y(n)) THEN
      vertex(1) = n
      vertex(2) = 1
      ymin = y(n)
      ymax = y(1)
    ELSE
      vertex(1) = 1
      vertex(2) = n
      ymin = y(1)
      ymax = y(n)
    END IF

    DO i = 2, n-1
      temp = y(i)
      IF (temp < ymin) THEN
        vertex(1) = i
        ymin = temp
      ELSE IF (temp > ymax) THEN
        vertex(2) = i
        ymax = temp
      END IF
    END DO

    nvert = 2
    IF (ymax == ymin) nvert = 1
    RETURN
  END IF

  !  Set up two initial lists of points; those points above & those below the
  !  line joining the first two vertices.    next(i) will hold the pointer to the
  !  point furthest from the line joining vertex(i) to vertex(i+1) on the left
  !  hand side.

  i1 = vertex(1)
  i2 = vertex(2)
  iwk(i1) = -1
  iwk(i2) = -1
  dx = xmax - xmin
  y1 = y(i1)
  dy = y(i2) - y1
  dmax = zero
  dmin = zero
  next(1) = -1
  next(2) = -1

  DO i = 1, n
    IF (i == vertex(1) .OR. i == vertex(2)) CYCLE
    dist = (y(i) - y1)*dx - (x(i) - xmin)*dy
    IF (dist > zero) THEN
      iwk(i1) = i
      i1 = i
      IF (dist > dmax) THEN
        next(1) = i
        dmax = dist
      END IF
    ELSE IF (dist < zero) THEN
      iwk(i2) = i
      i2 = i
      IF (dist < dmin) THEN
        next(2) = i
        dmin = dist
      END IF
    END IF
  END DO

  !  Ends of lists are indicated by pointers to -ve positions.

  iwk(i1) = -1
  iwk(i2) = -1
  nvert = 2

  j = 1

  !  Start of main process.

  !  Introduce new vertex between vertices j & j+1, if one has been found.
  !  Otherwise increase j.   Exit if no more vertices.

  40 IF (next(j) < 0) THEN
  IF (j == nvert) RETURN
  j = j + 1
  GO TO 40
  END IF

  jp1 = j + 1
  DO i = nvert, jp1, -1
    vertex(i+1) = vertex(i)
    next(i+1) = next(i)
  END DO
  jp2 = jp1 + 1
  nvert = nvert + 1
  IF (jp2 > nvert) jp2 = 1
  i1 = vertex(j)
  i2 = next(j)
  i3 = vertex(jp2)
  vertex(jp1) = i2

  !  Process the list of points associated with vertex j.   New list at vertex j
  !  consists of those points to the left of the line joining it to the new
  !  vertex (j+1).   Similarly for the list at the new vertex.
  !  Points on or to the right of these lines are dropped.

  x1 = x(i1)
  x2 = x(i2)
  y1 = y(i1)
  y2 = y(i2)
  dx1 = x2 - x1
  dx2 = x(i3) - x2
  dy1 = y2 - y1
  dy2 = y(i3) - y2
  DMAX1 = zero
  dmax2 = zero
  next(j) = -1
  next(jp1) = -1
  i2save = i2
  i2next = iwk(i2)
  i = iwk(i1)
  iwk(i1) = -1
  iwk(i2) = -1

  60 IF (i /= i2save) THEN
    dist = (y(i) - y1)*dx1 - (x(i) - x1)*dy1
    IF (dist > zero) THEN
      iwk(i1) = i
      i1 = i
      IF (dist > DMAX1) THEN
        next(j) = i
        DMAX1 = dist
      END IF
    ELSE
      dist = (y(i) - y2)*dx2 - (x(i) - x2)*dy2
      IF (dist > zero) THEN
        iwk(i2) = i
        i2 = i
        IF (dist > dmax2) THEN
          next(jp1) = i
          dmax2 = dist
        END IF
      END IF
    END IF
    i = iwk(i)
  ELSE
    i = i2next
  END IF

  !  Get next point from old list at vertex j.

  IF (i > 0) GO TO 60

  !  End lists with -ve values.

  iwk(i1) = -1
  iwk(i2) = -1

  GO TO 40
END SUBROUTINE envelope


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
!  Discussion:
!
!    CDIV computes:
!
!      (CR,CI) = (AR,AI) / (BR,BI)
!
!    using real(kind=8) arithmetic.
!
!  Modified:
!    15 March 2001
!
!  Parameters:
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
!    30 November 1998
!
!  Author:
!    John Burkardt
!
!  Parameters:
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
!    24 March 2000
!
!  Author:
!    John Burkardt
!
!  Parameters:
!    Input, integer LDA, the leading dimension of A.
!    Input, integer N, the order of A.
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

!==============================================================================
! Closing ports
!==============================================================================
subroutine close_all

  implicit none

  if (c_compacity         ==1) close(100)
  if (c_wall_pos          ==1) close(101)
  if (c_wall_f            ==1) close(102)
  if (c_coordination      ==1) close(103)
  if (c_contact_anisotropy==1) close(107)
  if (c_force_anisotropy  ==1) close(108)
  if (c_branch_anisotropy ==1) close(109)
  if (c_qoverp            ==1) close(113)
  if (c_qp_walls_nfric    ==1) close(114)
  if (c_avg_c_dist        ==1) close(124)

end subroutine close_all

end program post2D_defo
