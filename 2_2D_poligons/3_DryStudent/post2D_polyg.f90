!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! POST-PROCESSOR !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! POLYGONS!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program post2D_polyg

implicit none

type  ::  T_CORPS
   integer                                  ::  num,lieu,nb_vertex,ctc
   real(kind=8)                             ::  radius,Area,ax1,ax2, s11, s12, s21, s22
   real(kind=8)                             ::  p, q, qoverp, crit, Area_temp
   real(kind=8),dimension(3)                ::  center,center_ref
   real(kind=8),dimension(3)                ::  V
   real(kind=8),allocatable,dimension(:,:)  ::  vertex,vertex_ref
   character(len=5)                         ::  behav
end type T_CORPS

type  ::  T_CONTACT
   integer                                  ::  icdent,ianent,sect
   integer                                  ::  type, status_points
   real(kind=8),dimension(3)                ::  n
   real(kind=8),dimension(3)                ::  t
   real(kind=8),dimension(3)                ::  coor_ctc
   real(kind=8)                             ::  rn,rt,rs
   real(kind=8)                             ::  vn,vl,vs
   real(kind=8)                             ::  gap0, gap, cd_len, tang_disp
   character(len=5)                         ::  nature, status, law
end type T_CONTACT

type  ::  T_DOF
   integer                                  ::  num
   real(kind=8),dimension(3)                ::  X
   real(kind=8),dimension(3)                ::  V
end type T_DOF

type(T_CORPS),allocatable,dimension(:)      ::  TAB_POLYG,TAB_PLAN,TAB_POLYG_INI
type(T_CONTACT),allocatable,dimension(:)    ::  TAB_CONTACT,TAB_CONTACT_persiste
type(T_CONTACT),allocatable,dimension(:)    ::  TAB_CONTACT_POLYG


!DEFINITIONS Vloc-Dof-Bodies
!================================================
integer                                     ::  n_walls=0,fin=0,n_particles=0
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
integer                                    ::  calcul_coordination=0, &
                                               calcul_compacity=0,calcul_qoverp=0,calcul_anisotropy_contact=0,&
                                               calcul_anisotropy_force=0, calcul_anisotropy_branch = 0, &
                                               walls_position = 0, walls_f = 0


!Variables DavidC.
real(kind=8), allocatable, dimension(:,:)  ::  tab_temp_disp
real(kind=8), allocatable, dimension(:)    ::  tab_temp_rot
integer                                    ::  l, nb_part_invi, nb_part_visi, nb_vi_ini
integer                                    ::  N_PLPLx, nc_poly, n_ctc_cohe
integer                                    ::  N_PLJCx, nc_wall
logical                                    ::  copy_input=.false.
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
  
  if (command=='VLOC_DOF_BODIES              :') then
    read(1,*) n_walls
    read(1,*) n_particles
    cycle
  end if
  
  if (command=='COORDINATION                 :') then
    calcul_coordination=1
    open (unit=100,file='./POSTPRO/COORDINATION.DAT',status='replace')
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
    walls_position=1
    open (unit=108,file='./POSTPRO/WALLSPOS.DAT',status='replace')
    cycle
  end if
  
  if (command=='WALLS FORCES                 :') then
    walls_f=1
    open (unit=109,file='./POSTPRO/WALLSFORCE.DAT',status='replace')
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
    if (calcul_compacity            == 1)  call compacity
    if (calcul_qoverp               == 1)  call qoverp
    if (calcul_anisotropy_contact   == 1)  call anisotropy_contact
    if (calcul_anisotropy_force     == 1)  call anisotropy_force
    if (calcul_anisotropy_branch    == 1)  call anisotropy_branch

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
  
  integer                       ::  i,err,num_part,nb_vertex,j, gen 
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
  integer                       ::  icdent,ianent,cpt_type, l_stat_p
  real(kind=8)                  ::  n1,n2,n3
  real(kind=8)                  ::  rn,rt,rs,rn_temp,rt_temp,rs_temp
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
        else
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
          
        end if
        
        if (i==n_particles+n_walls) exit
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
          
        else
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
          
        end if
      end if
      
      if (i==n_particles+n_walls) exit
    end do
  end if
  
  N_PLJCx = 0
  N_PLPLx = 0
  
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
    end do
    close(4)
    
    ! Allocating the corresponding space for the number of contacts
    if (allocated(TAB_CONTACT)) deallocate(TAB_CONTACT)
    allocate(TAB_CONTACT(N_PLPLx+N_PLJCx))
    
    ! Storing the contacts information
    ! In the array the particle-particle contacts are stored first
    nb_ligneCONTACT=0
    nc_poly = 0
    nc_wall = 0
    
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
        read(4,'(8X,i5,30X,A5,9X,i5,30X,A5)') icdent, i_law, ianent,status
        read(4,'(29X, 3(5x,D14.7,2X))') rt,rn,rs
        read(4,'(29X, 3(5x,D14.7,2X))') vlt,vln,vls
        read(4,'(29X, 1(5x,D14.7,2X))') gap_c
        read(4,'(29X, 3(5x,D14.7,2X))') n1,n2,n3
        read(4,'(29X, 3(5x,D14.7,2X))') coor_ctc1,coor_ctc2,coor_ctc3

        if (status(1:1) == "W") then
          read(4,*)
          read(4,'(3(1X,D14.7))') gap0, cd_len, tang_disp
        end if

        nb_ligneCONTACT=nb_ligneCONTACT+1

        TAB_CONTACT(nc_poly)%icdent=icdent
        TAB_CONTACT(nc_poly)%ianent=ianent
        TAB_CONTACT(nc_poly)%status=status

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

        TAB_CONTACT(nc_poly)%law=i_law
        if (status(1:1) == "W") then
          TAB_CONTACT(nc_poly)%gap0=gap0
          TAB_CONTACT(nc_poly)%gap = gap_c
          TAB_CONTACT(nc_poly)%cd_len=cd_len
          TAB_CONTACT(nc_poly)%tang_disp = tang_disp
        else
          TAB_CONTACT(nc_poly)%gap0=-99.
          TAB_CONTACT(nc_poly)%gap = gap_c
          TAB_CONTACT(nc_poly)%cd_len=-99.
          TAB_CONTACT(nc_poly)%tang_disp = -99.
        end if

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
      end if
      
      if (text2 == '$icdan  PLJCx') then
        nc_wall = nc_wall + 1
        read(4,*)
        read(4,'(8X,i5,30X,A5,9X,i5,30X,A5)') icdent, i_law, ianent,status
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

        ! Default initialization
        TAB_CONTACT(N_PLPLx+nc_wall)%gap0=-99.
        TAB_CONTACT(N_PLPLx+nc_wall)%gap = gap_c
        TAB_CONTACT(N_PLPLx+nc_wall)%cd_len=-99.
        TAB_CONTACT(N_PLPLx+nc_wall)%tang_disp = -99.
      end if

      if (nc_poly+nc_wall == N_PLPLx + N_PLJCx) then
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
      TAB_CONTACT_POLYG(nb_ligneCONTACT_POLYG)%status_points   = TAB_CONTACT(i)%status_points
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
  ! print*, nb_ligneCONTACT_POLYG
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
    if (abs(Rnik) .le. 1.D-8) cycle
    
    Lik(1) = TAB_POLYG(cd)%center(1)-TAB_POLYG(an)%center(1)
    Lik(2) = TAB_POLYG(cd)%center(2)-TAB_POLYG(an)%center(2)
    Fik(1) = (Rnik*nik(1)+Rtik*tik(1))
    Fik(2) = (Rnik*nik(2)+Rtik*tik(2))
    
    Moment(1,1:2) = Fik(1)*Lik(1:2) + Moment(1,1:2)
    Moment(2,1:2) = Fik(2)*Lik(1:2) + Moment(2,1:2)
  end do
  
  Moment = Moment / (height*width)
  
  lda  = 2
  matz = 1
  
  ! Finding the eigenvalues of the global stress tensor
  call rg(lda, 2, Moment, wr, wi, matz, localframe, ierror)
  S1 = max( wr(1),wr(2) )
  S2 = min( wr(1),wr(2) )

  ! Computing the orientation of the stress tensor 
  ang_c = abs(0.5*atan2(-2*Moment(1,2),(Moment(1,1)-Moment(2,2))))

  ang_c = ang_c*180/pi
  
  if (first_over_all) then
    write(103,*) '#   time     ', '      S1     ', '      S2     ', '(S1-S2)/(S1+S2)', '  theta_s  '
  endif
  
  write(103,'(8(1X,E12.5))') time, S1,S2,(S1-S2)/(S1+S2), ang_c, Moment(1,1), Moment(1,2), Moment(2,2)
  
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
    if (abs(TAB_CONTACT_POLYG(i)%rn) .gt. 1e-9) then
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
! Fermeture des fichiers
!==============================================================================

 subroutine close_all
 implicit none
 
   if (calcul_coordination              == 1)  close (100)
   if (calcul_compacity                 == 1)  close (102)
   if (calcul_qoverp                    == 1)  close (103)
   if (calcul_anisotropy_contact        == 1)  close (104)
   if (calcul_anisotropy_force          == 1)  close (105)
   if (calcul_anisotropy_branch         == 1)  close (106)
   if (walls_position                   == 1)  close (108)
   if (walls_f                          == 1)  close (109)
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
