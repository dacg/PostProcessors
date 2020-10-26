!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!     VISU    !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!! SPHERES !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


program visu_spheres

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

end type T_CORPS

type  ::  T_CONTACT
  integer                                        ::  icdent,ianent,type
  real(kind=8),dimension(3)                      ::  n,t,s
  real(kind=8),dimension(3)                      ::  coor_ctc
  real(kind=8)                                   ::  rn,rt,rs,vln,vls,vlt,gap
  integer                                        ::  sect
  character(len=5)                               ::  nature, statut
end type T_CONTACT

type(T_CORPS),dimension(:),allocatable          ::  TAB_SPHERE,TAB_PLAN
type(T_CONTACT),dimension(:),allocatable        ::  TAB_CONTACT

!DEFINITIONS Vloc-Dof-Bodies
!================================================================
integer                                         ::  n_walls=0, n_particles=0
integer                                         ::  compteur_clout=0
logical                                         ::  fin_post=.false. , first_bodies=.true.
integer                                         ::  step,nb_ligneCONTACT
real(kind=8)                                    ::  time,Mean_total_Normal_force, pfric_coef

!DEFINITIONS BOITE
!================================================================
real(kind=8)                                    ::  H_ini,Long_ini,Larg_ini
logical                                         ::  deformation_ini = .true.


!Global variables
!================================================================
character(len=30)                               ::  command
real(kind=8)                                    ::  PI = 3.14159265359, &
                                                    width,large,height


! Variables David C.
logical                                          ::  first_over_all = .true.
logical                                          ::  v_negative_f = .false.
real*8                                           ::  f_b_scale

!================================================================
!Reading input file
!================================================================
print*,'-------------------------------------------'
print*,'!                                         !'
print*,'!       Visualization of spheres          !'
print*,'!                                         !'
print*,'-------------------------------------------'
print*,'Reading INPUT'

open(unit=1,file='POST_INPUT.DAT',status='old')
! Reading the commands
do

  read(1,'(A30)') command
  print*,command

  if (command=='INFO SAMPLE                  :') then
    read(1,*) n_walls
    read(1,*) n_particles
    read(1,*) pfric_coef
    cycle
  end if

  if (command=='NEGATIVE FORCES              :') then
    v_negative_f = .true.
    cycle
  end if

  if (command=='FORCE BARS SCALE             :') then
    read(1,*) f_b_scale
    cycle
  end if

  if (command=='END                          :') then
    close(1)
    exit
  end if
end do

!==============================================================================
! Reading OUTBOX and printing vtk files
!==============================================================================

  do
    ! Calling the Subroutines
    call read_Vloc_dof_bodies

    call box_size

    call draw

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
!Reading OUTBOX
!======================================================================
subroutine read_Vloc_dof_bodies

  implicit none

  integer                              ::  i=0,j=0,n_spsp,n_sppl,num_part,err
  integer                              ::  nb_SPSPx,nb_SPPLx,icdent,ianent
  real(kind=8)                         ::  rn,rs,rt,t1,t2,t3,n1,n2,n3,s1,s2,s3,rayon,vls,vln,vlt,gap
  character(len=6)                     ::  text
  character(len=5)                     ::  statut,behav,color
  character(len=13)                    ::  text2
  real(kind=8)                         ::  center1,center2,center3,ax1r,ax2r,ax3r,coor1,coor2,coor3
  real(kind=8)                         ::  coor_ctc1,coor_ctc2,coor_ctc3
  real(kind=8)                         ::  Volume
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
          read(2,'(29X, 3(5x,D14.7,2X))') ax1r,ax2r,ax3r
          if (first_over_all) then
            TAB_PLAN(i-n_particles)%ax1=ax1r
            TAB_PLAN(i-n_particles)%ax2=ax2r
            TAB_PLAN(i-n_particles)%ax3=ax3r
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
                                                  TAB_PLAN(i-n_particles)%vertex_ref(3,j) *rot(1,3)+&
                                                  TAB_PLAN(i-n_particles)%center(1)

            TAB_PLAN(i-n_particles)%vertex(2,j) = TAB_PLAN(i-n_particles)%vertex_ref(1,j) *rot(2,1)+&
                                                  TAB_PLAN(i-n_particles)%vertex_ref(2,j) *rot(2,2)+&
                                                  TAB_PLAN(i-n_particles)%vertex_ref(3,j) *rot(2,3)+&
                                                  TAB_PLAN(i-n_particles)%center(2)

            TAB_PLAN(i-n_particles)%vertex(3,j) = TAB_PLAN(i-n_particles)%vertex_ref(1,j) *rot(3,1)+&
                                                  TAB_PLAN(i-n_particles)%vertex_ref(2,j) *rot(3,2)+&
                                                  TAB_PLAN(i-n_particles)%vertex_ref(3,j) *rot(3,3)+&
                                                  TAB_PLAN(i-n_particles)%center(3)
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
    if (TAB_CONTACT(i)%rn == 0.D0 .and. (v_negative_f .eqv. .false.)) cycle
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
    if (TAB_CONTACT(i)%rn .le. 0. .and. (v_negative_f .eqv. .false.)) cycle

    TAB_SPHERE(icdent)%nctc = TAB_SPHERE(icdent)%nctc + 1
    TAB_SPHERE(ianent)%nctc = TAB_SPHERE(ianent)%nctc + 1
  end do

  print*, 'Nb particules :', n_particles
  print*, 'Nb Contacts :', nb_ligneCONTACT

  print*, 'Step:', step
  print*, 'Time:', time

end subroutine read_Vloc_dof_bodies


!==============================================================================
!Calcul des deformation
!==============================================================================
subroutine box_size

  implicit none

  integer               :: i
  real(kind=8)          :: Xplan_max,Xplan_min,Yplan_max,Yplan_min,Zplan_max,Zplan_min

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

end subroutine box_size


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
  n_l_vertices = n_l_vertices + (n_walls*8)

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
          ! Writing the 3 components
          write(114,'(E13.6,A)', advance = 'no') local_sphere(j,1), ' '
          write(114,'(E13.6,A)', advance = 'no') local_sphere(j,2), ' '
          write(114,'(E13.6,A)', advance = 'no') local_sphere(j,3), ' '
        else
          write(114,'(E13.6,A)', advance = 'no') local_sphere(j,1), ' '
          write(114,'(E13.6,A)', advance = 'no') local_sphere(j,2), ' '
          write(114,'(E13.6,A)') local_sphere(j,3), ' '
          k = 0
        end if
      end do
    else
      ! Writing the vertices of the walls
      do j=1, 8

        k = k +1
        if(k .lt. 3) then
          write(114,'(E13.6,A)', advance = 'no') TAB_PLAN(i-n_particles)%vertex(1,j), ' '
          write(114,'(E13.6,A)', advance = 'no') TAB_PLAN(i-n_particles)%vertex(2,j), ' '
          write(114,'(E13.6,A)', advance = 'no') TAB_PLAN(i-n_particles)%vertex(3,j), ' '
        else
          write(114,'(E13.6,A)', advance = 'no') TAB_PLAN(i-n_particles)%vertex(1,j), ' '
          write(114,'(E13.6,A)', advance = 'no') TAB_PLAN(i-n_particles)%vertex(2,j), ' '
          write(114,'(E13.6,A)') TAB_PLAN(i-n_particles)%vertex(3,j), ' '
          k = 0
        end if
      end do
    end if
  end do

  write(114, '(A)') ' '

  ! Finding the total number of faces
  n_l_faces = n_particles * sphere_n_faces

  ! Writing the conectivity between vertices
  ! Plus 6 faces each wall
  write(114,'(A)', advance='no') 'POLYGONS '
  write(114, '(2(I8,A))') (n_l_faces + n_walls*6), ' ' , (n_l_faces*4)+(n_walls*5*6)

  curr_l_faces = 0

  do i=1, n_particles
    !write(114,*) 'Particle', i
    do j=1, sphere_n_faces
      ! We have always triangles
      write(114, '(I1,A)', advance= 'no') 3, ' '

      ! First face...
      write(114, '(I7,A)', advance = 'no') sphere_connect(j,1) + curr_l_faces -1, ' '

      ! Second face
      write(114, '(I7,A)', advance = 'no') sphere_connect(j,2) + curr_l_faces -1, ' '

      ! Third face
      write(114, '(I7,A)') sphere_connect(j,3) + curr_l_faces -1, ' '
    end do
    curr_l_faces = curr_l_faces + sphere_n_vertices
  end do

  ! Writing the conectivity for the walls
  ! For all the walls
  do i=1, n_walls
    ! Each one with 6 faces

    !!!!! 1st (1-2-3-4)
    write(114, '(I1,A)', advance= 'no') 4, ' '
    write(114, '(I7,A)', advance = 'no') curr_l_faces, ' '
    write(114, '(I7,A)', advance = 'no') curr_l_faces +1, ' '
    write(114, '(I7,A)', advance = 'no') curr_l_faces +2, ' '
    write(114, '(I7)') curr_l_faces +3

    !!!!! 2nd (5-6-7-8)
    write(114, '(I1,A)', advance= 'no') 4, ' '
    write(114, '(I7,A)', advance = 'no') curr_l_faces + 4, ' '
    write(114, '(I7,A)', advance = 'no') curr_l_faces + 5, ' '
    write(114, '(I7,A)', advance = 'no') curr_l_faces + 6, ' '
    write(114, '(I7,A)') curr_l_faces +7

    !!!!! 3nd (1-5-8-4)
    write(114, '(I1,A)', advance= 'no') 4, ' '
    write(114, '(I7,A)', advance = 'no') curr_l_faces, ' '
    write(114, '(I7,A)', advance = 'no') curr_l_faces + 4, ' '
    write(114, '(I7,A)', advance = 'no') curr_l_faces + 7, ' '
    write(114, '(I7)') curr_l_faces +3

    !!!!! 4th (2-6-7-3)
    write(114, '(I1,A)', advance= 'no') 4, ' '
    write(114, '(I7,A)', advance = 'no') curr_l_faces + 1, ' '
    write(114, '(I7,A)', advance = 'no') curr_l_faces + 5, ' '
    write(114, '(I7,A)', advance = 'no') curr_l_faces + 6, ' '
    write(114, '(I7)') curr_l_faces +2

    !!!!! 5th (1-2-6-5)
    write(114, '(I1,A)', advance= 'no') 4, ' '
    write(114, '(I7,A)', advance = 'no') curr_l_faces, ' '
    write(114, '(I7,A)', advance = 'no') curr_l_faces + 1, ' '
    write(114, '(I7,A)', advance = 'no') curr_l_faces + 5, ' '
    write(114, '(I7)') curr_l_faces +4

    !!!!! 6th (4-3-7-8)
    write(114, '(I1,A)', advance= 'no') 4, ' '
    write(114, '(I7,A)', advance = 'no') curr_l_faces + 3, ' '
    write(114, '(I7,A)', advance = 'no') curr_l_faces + 2, ' '
    write(114, '(I7,A)', advance = 'no') curr_l_faces + 6, ' '
    write(114, '(I7)') curr_l_faces +7

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
  write(114,'(A,I1,I8,A)') 'Id ', 1, (n_l_faces +n_walls*6), ' float'
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
  write(114,'(A,I1,I8,A)') 'Material ', 1, (n_l_faces +n_walls*6), ' float'
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

  write(114,'(A,I1,I8,A)') 'Disp ', 3, (n_l_faces +n_walls*6), ' float'
  k = 0
  do i=1, n_particles + n_walls
    if (i .le. n_particles) then
      curr_l_vector(:) = TAB_SPHERE(i)%center(:) - TAB_SPHERE(i)%center_ref(:)
      do j=1, sphere_n_faces
        k=k+3
        if(k .le. 9) then
          write(114, '(3(E13.6,A))', advance='no') curr_l_vector(1), ' ', curr_l_vector(2), ' ', curr_l_vector(3), ' '
        else
          write(114, '(3(E13.6,A))') curr_l_vector(1), ' ', curr_l_vector(2), ' ', curr_l_vector(3)
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
          write(114, '(3(E13.6,A))', advance='no') curr_l_vector(1), ' ', curr_l_vector(2), ' ', curr_l_vector(3), ' '
        else
          write(114, '(3(E13.6,A))') curr_l_vector(1), ' ', curr_l_vector(2), ' ', curr_l_vector(3)
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

  write(114,'(A,I1,I8,A)') 'Veloc ', 3, (n_l_faces +n_walls*6), ' float'
  k = 0
  do i=1, n_particles + n_walls
    if (i .le. n_particles) then
      curr_l_vector(1) = TAB_SPHERE(i)%Vx
      curr_l_vector(2) = TAB_SPHERE(i)%Vy
      curr_l_vector(3) = TAB_SPHERE(i)%Vz

      do j=1, sphere_n_faces
        k=k+3
        if(k .le. 9) then
          write(114, '(3(E13.6,A))', advance='no') curr_l_vector(1), ' ', curr_l_vector(2), ' ', curr_l_vector(3), ' '
        else
          write(114, '(3(E13.6,A))') curr_l_vector(1), ' ', curr_l_vector(2), ' ', curr_l_vector(3)
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
          write(114, '(3(E13.6,A))', advance='no') curr_l_vector(1), ' ', curr_l_vector(2), ' ', curr_l_vector(3), ' '
        else
          write(114, '(3(E13.6,A))') curr_l_vector(1), ' ', curr_l_vector(2), ' ', curr_l_vector(3)
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

  write(114,'(A,I1,I8,A)') 'Z ', 1, (n_l_faces +n_walls*6), ' float'
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

  do i=1, nb_ligneCONTACT

    if(TAB_CONTACT(i)%nature /= 'SPSPx') cycle
    ! Removing contacts with walls
    l_cdt = TAB_CONTACT(i)%icdent
    l_ant = TAB_CONTACT(i)%ianent

    if (l_cdt .gt. n_particles .or. l_ant .gt. n_particles) cycle

    if (TAB_CONTACT(i)%rn .le. 0 .and. (v_negative_f .eqv. .false.)) cycle

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

  ! Force scale parameter --> Its a proportion of the mean radius
  l_force_scale = (ave_rad)*f_b_scale

  ! Writing
  k=0
  do i=1, nb_ligneCONTACT

    if (TAB_CONTACT(i)%nature /= 'SPSPx') cycle
    ! The particles ids
    l_cdt = TAB_CONTACT(i)%icdent
    l_ant = TAB_CONTACT(i)%ianent

    if (l_cdt .gt. n_particles .or. l_ant .gt. n_particles) cycle

    if (TAB_CONTACT(i)%rn .le. 0 .and. (v_negative_f .eqv. .false.)) cycle

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
      write(114, '(3(E13.6,A))', advance='no') vtk_cd_center(1)+l_force*v_l_t(1)+ l_force*v_l_s(1), ' ', &
                                               vtk_cd_center(2)+l_force*v_l_t(2)+ l_force*v_l_s(2), ' ', &
                                               vtk_cd_center(3)+l_force*v_l_t(3)+ l_force*v_l_s(3), ' '
      ! Candidate --- Second vertex.
      write(114, '(3(E13.6,A))', advance='no') vtk_cd_center(1)+l_force*v_l_t(1)- l_force*v_l_s(1), ' ', &
                                               vtk_cd_center(2)+l_force*v_l_t(2)- l_force*v_l_s(2), ' ', &
                                               vtk_cd_center(3)+l_force*v_l_t(3)- l_force*v_l_s(3), ' '
      ! Candidate --- Third vertex.
      write(114, '(3(E13.6,A))', advance='no') vtk_cd_center(1)-l_force*v_l_t(1)- l_force*v_l_s(1), ' ', &
                                               vtk_cd_center(2)-l_force*v_l_t(2)- l_force*v_l_s(2), ' ', &
                                               vtk_cd_center(3)-l_force*v_l_t(3)- l_force*v_l_s(3), ' '
      ! Candidate --- 4th vertex.
      write(114, '(3(E13.6,A))') vtk_cd_center(1)-l_force*v_l_t(1)+ l_force*v_l_s(1), ' ', &
                                 vtk_cd_center(2)-l_force*v_l_t(2)+ l_force*v_l_s(2), ' ', &
                                 vtk_cd_center(3)-l_force*v_l_t(3)+ l_force*v_l_s(3), ' '

      ! Antagonist --- First vertex.
      write(114, '(3(E13.6,A))', advance='no') vtk_an_center(1)+l_force*v_l_t(1)+ l_force*v_l_s(1), ' ', &
                                               vtk_an_center(2)+l_force*v_l_t(2)+ l_force*v_l_s(2), ' ', &
                                               vtk_an_center(3)+l_force*v_l_t(3)+ l_force*v_l_s(3), ' '
      ! Antagonist --- Second vertex.
      write(114, '(3(E13.6,A))', advance='no') vtk_an_center(1)+l_force*v_l_t(1)- l_force*v_l_s(1), ' ', &
                                               vtk_an_center(2)+l_force*v_l_t(2)- l_force*v_l_s(2), ' ', &
                                               vtk_an_center(3)+l_force*v_l_t(3)- l_force*v_l_s(3), ' '
      ! Antagonist --- Third vertex.
      write(114, '(3(E13.6,A))', advance='no') vtk_an_center(1)-l_force*v_l_t(1)- l_force*v_l_s(1), ' ', &
                                               vtk_an_center(2)-l_force*v_l_t(2)- l_force*v_l_s(2), ' ', &
                                               vtk_an_center(3)-l_force*v_l_t(3)- l_force*v_l_s(3), ' '
      ! Antagonist --- 4th vertex.
      write(114, '(3(E13.6,A))') vtk_an_center(1)-l_force*v_l_t(1)+ l_force*v_l_s(1), ' ', &
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
      write(114, '(3(E13.6,A))', advance='no') vtk_cd_center(1)+l_force*v_l_t(1)+ l_force*v_l_s(1), ' ', &
                                               vtk_cd_center(2)+l_force*v_l_t(2)+ l_force*v_l_s(2), ' ', &
                                               vtk_cd_center(3)+l_force*v_l_t(3)+ l_force*v_l_s(3), ' '
      ! Candidate --- Second vertex.
      write(114, '(3(E13.6,A))', advance='no') vtk_cd_center(1)+l_force*v_l_t(1)- l_force*v_l_s(1), ' ', &
                                               vtk_cd_center(2)+l_force*v_l_t(2)- l_force*v_l_s(2), ' ', &
                                               vtk_cd_center(3)+l_force*v_l_t(3)- l_force*v_l_s(3), ' '
      ! Candidate --- Third vertex.
      write(114, '(3(E13.6,A))', advance='no') vtk_cd_center(1)-l_force*v_l_t(1)- l_force*v_l_s(1), ' ', &
                                               vtk_cd_center(2)-l_force*v_l_t(2)- l_force*v_l_s(2), ' ', &
                                               vtk_cd_center(3)-l_force*v_l_t(3)- l_force*v_l_s(3), ' '
      ! Candidate --- 4th vertex.
      write(114, '(3(E13.6,A))') vtk_cd_center(1)-l_force*v_l_t(1)+ l_force*v_l_s(1), ' ', &
                                 vtk_cd_center(2)-l_force*v_l_t(2)+ l_force*v_l_s(2), ' ', &
                                 vtk_cd_center(3)-l_force*v_l_t(3)+ l_force*v_l_s(3), ' '

      ! Antagonist --- First vertex.
      write(114, '(3(E13.6,A))', advance='no') vtk_an_center(1)+l_force*v_l_t(1)+ l_force*v_l_s(1), ' ', &
                                               vtk_an_center(2)+l_force*v_l_t(2)+ l_force*v_l_s(2), ' ', &
                                               vtk_an_center(3)+l_force*v_l_t(3)+ l_force*v_l_s(3), ' '
      ! Antagonist --- Second vertex.
      write(114, '(3(E13.6,A))', advance='no') vtk_an_center(1)+l_force*v_l_t(1)- l_force*v_l_s(1), ' ', &
                                               vtk_an_center(2)+l_force*v_l_t(2)- l_force*v_l_s(2), ' ', &
                                               vtk_an_center(3)+l_force*v_l_t(3)- l_force*v_l_s(3), ' '
      ! Antagonist --- Third vertex.
      write(114, '(3(E13.6,A))', advance='no') vtk_an_center(1)-l_force*v_l_t(1)- l_force*v_l_s(1), ' ', &
                                               vtk_an_center(2)-l_force*v_l_t(2)- l_force*v_l_s(2), ' ', &
                                               vtk_an_center(3)-l_force*v_l_t(3)- l_force*v_l_s(3), ' '
      ! Antagonist --- 4th vertex.
      write(114, '(3(E13.6,A))') vtk_an_center(1)-l_force*v_l_t(1)+ l_force*v_l_s(1), ' ', &
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
      write(114, '(3(E13.6,A))', advance='no') vtk_cd_center(1)+l_force*v_l_t(1)+ l_force*v_l_s(1), ' ', &
                                               vtk_cd_center(2)+l_force*v_l_t(2)+ l_force*v_l_s(2), ' ', &
                                               vtk_cd_center(3)+l_force*v_l_t(3)+ l_force*v_l_s(3), ' '
      ! Candidate --- Second vertex.
      write(114, '(3(E13.6,A))', advance='no') vtk_cd_center(1)+l_force*v_l_t(1)- l_force*v_l_s(1), ' ', &
                                               vtk_cd_center(2)+l_force*v_l_t(2)- l_force*v_l_s(2), ' ', &
                                               vtk_cd_center(3)+l_force*v_l_t(3)- l_force*v_l_s(3), ' '
      ! Candidate --- Third vertex.
      write(114, '(3(E13.6,A))', advance='no') vtk_cd_center(1)-l_force*v_l_t(1)- l_force*v_l_s(1), ' ', &
                                               vtk_cd_center(2)-l_force*v_l_t(2)- l_force*v_l_s(2), ' ', &
                                               vtk_cd_center(3)-l_force*v_l_t(3)- l_force*v_l_s(3), ' '
      ! Candidate --- 4th vertex.
      write(114, '(3(E13.6,A))') vtk_cd_center(1)-l_force*v_l_t(1)+ l_force*v_l_s(1), ' ', &
                                 vtk_cd_center(2)-l_force*v_l_t(2)+ l_force*v_l_s(2), ' ', &
                                 vtk_cd_center(3)-l_force*v_l_t(3)+ l_force*v_l_s(3), ' '

      ! Antagonist --- First vertex.
      write(114, '(3(E13.6,A))', advance='no') vtk_an_center(1)+l_force*v_l_t(1)+ l_force*v_l_s(1), ' ', &
                                               vtk_an_center(2)+l_force*v_l_t(2)+ l_force*v_l_s(2), ' ', &
                                               vtk_an_center(3)+l_force*v_l_t(3)+ l_force*v_l_s(3), ' '
      ! Antagonist --- Second vertex.
      write(114, '(3(E13.6,A))', advance='no') vtk_an_center(1)+l_force*v_l_t(1)- l_force*v_l_s(1), ' ', &
                                               vtk_an_center(2)+l_force*v_l_t(2)- l_force*v_l_s(2), ' ', &
                                               vtk_an_center(3)+l_force*v_l_t(3)- l_force*v_l_s(3), ' '
      ! Antagonist --- Third vertex.
      write(114, '(3(E13.6,A))', advance='no') vtk_an_center(1)-l_force*v_l_t(1)- l_force*v_l_s(1), ' ', &
                                               vtk_an_center(2)-l_force*v_l_t(2)- l_force*v_l_s(2), ' ', &
                                               vtk_an_center(3)-l_force*v_l_t(3)- l_force*v_l_s(3), ' '
      ! Antagonist --- 4th vertex.
      write(114, '(3(E13.6,A))') vtk_an_center(1)-l_force*v_l_t(1)+ l_force*v_l_s(1), ' ', &
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
    write(114, '(I7,A)', advance = 'no') curr_l_faces, ' '
    write(114, '(I7,A)', advance = 'no') curr_l_faces +1, ' '
    write(114, '(I7,A)', advance = 'no') curr_l_faces +2, ' '
    write(114, '(I7)') curr_l_faces +3

    !!!!! 2nd (5-6-7-8)
    write(114, '(I1,A)', advance= 'no') 4, ' '
    write(114, '(I7,A)', advance = 'no') curr_l_faces + 4, ' '
    write(114, '(I7,A)', advance = 'no') curr_l_faces + 5, ' '
    write(114, '(I7,A)', advance = 'no') curr_l_faces + 6, ' '
    write(114, '(I7)') curr_l_faces +7

    !!!!! 3nd (1-5-8-4)
    write(114, '(I1,A)', advance= 'no') 4, ' '
    write(114, '(I7,A)', advance = 'no') curr_l_faces, ' '
    write(114, '(I7,A)', advance = 'no') curr_l_faces + 4, ' '
    write(114, '(I7,A)', advance = 'no') curr_l_faces + 7, ' '
    write(114, '(I7)') curr_l_faces +3


    !!!!! 4th (2-6-7-3)
    write(114, '(I1,A)', advance= 'no') 4, ' '
    write(114, '(I7,A)', advance = 'no') curr_l_faces + 1, ' '
    write(114, '(I7,A)', advance = 'no') curr_l_faces + 5, ' '
    write(114, '(I7,A)', advance = 'no') curr_l_faces + 6, ' '
    write(114, '(I7)') curr_l_faces +2

    !!!!! 5th (1-2-6-5)
    write(114, '(I1,A)', advance= 'no') 4, ' '
    write(114, '(I7,A)', advance = 'no') curr_l_faces, ' '
    write(114, '(I7,A)', advance = 'no') curr_l_faces + 1, ' '
    write(114, '(I7,A)', advance = 'no') curr_l_faces + 5, ' '
    write(114, '(I7)') curr_l_faces +4

    !!!!! 6th (4-3-7-8)
    write(114, '(I1,A)', advance= 'no') 4, ' '
    write(114, '(I7,A)', advance = 'no') curr_l_faces + 3, ' '
    write(114, '(I7,A)', advance = 'no') curr_l_faces + 2, ' '
    write(114, '(I7,A)', advance = 'no') curr_l_faces + 6, ' '
    write(114, '(I7)') curr_l_faces +7

    curr_l_faces = curr_l_faces + 8
  end do

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!! FORCES FIELDS !!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! How many fields there will be?
  n_l_fields = 4                          ! They are: RN, RT, Friction Mobilization Fric_M, TracComp
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

    Lik(:) = TAB_SPHERE(l_cdt)%center(:)-TAB_SPHERE(l_ant)%center(:)

    if (TAB_CONTACT(i)%rn .le. 0 .and. (v_negative_f .eqv. .false.)) cycle

    ! Being careful with the periodic case
    if (sqrt(Lik(1)**2 + Lik(2)**2 + Lik(3)**2) .gt. (TAB_SPHERE(l_cdt)%Rmax + TAB_SPHERE(l_ant)%Rmax)) then 
      do j=1, 12
        k=k+1
        if(k .le. 9) then
          write(114, '(E13.6)', advance='no') TAB_CONTACT(i)%rn
        else
          write(114, '(E13.6)') TAB_CONTACT(i)%rn
          k=0
        end if
      end do
    else
      do j=1, 6
        k=k+1
        if(k .le. 9) then
          write(114, '(E13.6)', advance='no') TAB_CONTACT(i)%rn
        else
          write(114, '(E13.6)') TAB_CONTACT(i)%rn
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

    Lik(:) = TAB_SPHERE(l_cdt)%center(:)-TAB_SPHERE(l_ant)%center(:)

    if (TAB_CONTACT(i)%rn .le. 0 .and. (v_negative_f .eqv. .false.)) cycle

    ! Being careful with the periodic case
    if (sqrt(Lik(1)**2 + Lik(2)**2 + Lik(3)**2) .gt. (TAB_SPHERE(l_cdt)%Rmax + TAB_SPHERE(l_ant)%Rmax)) then 
      do j=1, 12
        k=k+1
        if(k .le. 9) then
          write(114, '(E13.6)', advance='no') sqrt(TAB_CONTACT(i)%rt**2 + TAB_CONTACT(i)%rs**2)
        else
          write(114, '(E13.6)') sqrt(TAB_CONTACT(i)%rt**2 + TAB_CONTACT(i)%rs**2)
          k=0
        end if
      end do
    else
      do j=1, 6
        k=k+1
        if(k .le. 9) then
          write(114, '(E13.6)', advance='no') sqrt(TAB_CONTACT(i)%rt**2 + TAB_CONTACT(i)%rs**2)
        else
          write(114, '(E13.6)') sqrt(TAB_CONTACT(i)%rt**2 + TAB_CONTACT(i)%rs**2)
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

    Lik(:) = TAB_SPHERE(l_cdt)%center(:)-TAB_SPHERE(l_ant)%center(:)

    if (TAB_CONTACT(i)%rn .le. 0 .and. (v_negative_f .eqv. .false. )) cycle

    ! Attention.. periodic case
    if (sqrt(Lik(1)**2 + Lik(2)**2 + Lik(3)**2) .gt. (TAB_SPHERE(l_cdt)%Rmax + TAB_SPHERE(l_ant)%Rmax)) then
      do j=1, 12
        k=k+1
        if(k .le. 9) then
          ! Careful--> I'm writing the coefficient of friction explicitly
          write(114, '(E13.6)', advance='no') sqrt(TAB_CONTACT(i)%rt**2 + TAB_CONTACT(i)%rs**2)/(pfric_coef*TAB_CONTACT(i)%rn)
        else
          write(114, '(E13.6)') sqrt(TAB_CONTACT(i)%rt**2 + TAB_CONTACT(i)%rs**2)/(pfric_coef*TAB_CONTACT(i)%rn)
          k=0
        end if
      end do
    else
      do j=1, 6
        k=k+1
        if(k .le. 9) then
          ! Careful--> I'm writing the coefficient of friction explicitly
          write(114, '(E13.6)', advance='no') sqrt(TAB_CONTACT(i)%rt**2 + TAB_CONTACT(i)%rs**2)/(pfric_coef*TAB_CONTACT(i)%rn)
        else
          write(114, '(E13.6)') sqrt(TAB_CONTACT(i)%rt**2 + TAB_CONTACT(i)%rs**2)/(pfric_coef*TAB_CONTACT(i)%rn)
          k=0
        end if
      end do
    end if
  end do
  ! And jump
  write(114, '(A)') ' '

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Traction (-1) Compression (1)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Naming the field, the dimension of the data, and number of lines (as before FIELD), and data type
  write(114,'(A,I1,I8,A)') 'TracComp ', 1, (f_counter+f_counter_periodic*2)*6, ' float'
  k = 0
  do i=1, nb_ligneCONTACT
    if (TAB_CONTACT(i)%nature /= 'SPSPx') cycle
    l_cdt = TAB_CONTACT(i)%icdent
    l_ant = TAB_CONTACT(i)%ianent

    Lik(:) = TAB_SPHERE(l_cdt)%center(:)-TAB_SPHERE(l_ant)%center(:)

    if (TAB_CONTACT(i)%rn .le. 0 .and. (v_negative_f .eqv. .false.)) cycle

    ! Being careful with the periodic case
    if (sqrt(Lik(1)**2 + Lik(2)**2 + Lik(3)**2) .gt. (TAB_SPHERE(l_cdt)%Rmax + TAB_SPHERE(l_ant)%Rmax)) then 
      do j=1, 12
        k=k+1
        if(k .le. 9) then
          if(TAB_CONTACT(i)%rn .le. 0.) then
            write(114, '(E13.6)', advance='no') -1.0
          else
            write(114, '(E13.6)', advance='no')  1.0
          end if
        else
          if(TAB_CONTACT(i)%rn .le. 0.) then
            write(114, '(E13.6)') -1.0
          else
            write(114, '(E13.6)')  1.0
          end if
          k=0
        end if
      end do
    else
      do j=1, 6
        k=k+1
        if(k .le. 9) then
          if(TAB_CONTACT(i)%rn .le. 0.) then
            write(114, '(E13.6)', advance='no') -1.0
          else
            write(114, '(E13.6)', advance='no')  1.0
          end if
        else
          if(TAB_CONTACT(i)%rn .le. 0.) then
            write(114, '(E13.6)') -1.0
          else
            write(114, '(E13.6)')  1.0
          end if
          k=0
        end if
      end do
    end if
  end do
  ! And jump
  write(114, '(A)') ' '

  close(114)

  print*, 'Drawing in vtk  ---> Ok!'

end subroutine draw


!==================================================================================================
!Closing units
!=================================================================================================
subroutine close_all

  implicit none

  print*,'---> fermeture des fichiers'

  deallocate(TAB_SPHERE)

  print*,'--> Fin du postraitement <--'
  stop

end subroutine close_all


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

end program visu_spheres
