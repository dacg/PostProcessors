!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!  POST-PROCESSOR !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! POLYHEDRA !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Notes: - This code is build to read up to 10000 OUTBOX frames

program post3D_polyhe

implicit none

type  ::  T_CORPS
  integer                                        ::  nb_vertex,n_faces,id,nb_ctc
  integer,dimension(:,:),pointer                 ::  face
  real(kind=8)                                   ::  Rmax, volume, I1, I2, I3
  real(kind=8)                                   ::  Vx, Vy, Vz, Vrx, Vry, Vrz
  real(kind=8)                                   ::  Wx, Wy, Wz
  real(kind=8)                                   ::  ax1, ax2, ax3
  real(kind=8)                                   ::  Norme_V, Norme_Vr, Force, Sphericity, Aire
  real(kind=8)                                   ::  nctc, nctcslide, nctcstick, nctcs, nctcd, nctct
  real(kind=8),dimension(3)                      ::  center,center_ref
  real(kind=8),dimension(3,3)                    ::  momentI, Rot
  real(kind=8),dimension(:,:),pointer            ::  vertex, vertex_ref
  character(len=5)                               ::  behav,color

end type T_CORPS

type  ::  T_FACE
  integer, dimension(2)                          ::  face_id
  real(kind=8)                                   ::  face_status, f_area
  real(kind=8), dimension(3)                     ::  nface_vector
end type T_FACE

type  ::  T_CONTACT
  logical                                        ::  deja_compte
  integer                                        ::  icdent,ianent,type, cdver
  integer                                        ::  sect, loc
  real(kind=8)                                   ::  x_status
  real(kind=8),dimension(3)                      ::  n,t,s
  real(kind=8),dimension(3)                      ::  coor_ctc
  real(kind=8)                                   ::  rn,rt,rs,vls,vln,vlt, gapTT
  character(len=5)                               ::  nature,status, i_law
  character(len=62)                              ::  str_internal
end type T_CONTACT

type(T_CORPS),dimension(:),pointer               ::  TAB_POLY,TAB_PLAN
type(T_CONTACT),dimension(:),allocatable         ::  TAB_CONTACT,TAB_CONTACT_POLYR

! DEFINITIONS Vloc-Dof-Bodies
! ================================================
integer                                          ::  n_walls=0,n_particles=0
integer                                          ::  compteur_clout=0
logical                                          ::  fin_post=.false.
integer                                          ::  step,nb_ligneCONTACT,&
                                                     nb_ligneCONTACT_POLYR, &
                                                     nf_PRPRx,nf_PRPLx
real(kind=8)                                     ::  time,Mean_total_Normal_force=0


!COMMANDES D APPELLE
!================================================
character(len=30)                                ::  command
real(kind=8)                                     ::  PI = 3.14159265359
integer                                          ::  c_draw = 0, c_toppling_v

! Variables DC.
logical                                          ::  first_over_all = .true.
logical                                          ::  find_ini_vt, find_end_vt
real*8                                           ::  time_ini, time_end

!================================================
!Lecture du fichier de commande
!================================================
print*,'-------------------------------------------'
print*,'!                                         !'
print*,'!      Module de Post-traitement 3D       !'
print*,'!         from Emilien Azema Code         !'
print*,'!                 David C.                !'
print*,'-------------------------------------------'
print*,'Reading INPUT'

open(unit=1,file='POST_INPUT.DAT',status='old')
! Reading the commands
do
  read(1,'(A30)') command

  if (command=='VLOC_DOF_BODIES              :') then
    read(1,*) n_walls
    read(1,*) n_particles
    cycle
  end if

  if (command=='DRAW                         :') then
    c_draw = 1
    ! Uses port 121
    cycle
  end if

  if (command=='TOPPLING SPEED               :') then
    c_toppling_v = 1
    find_ini_vt = .false.
    find_end_vt = .false.
    open(unit=110, file='./POSTPRO/VELOC_FALL.OUT', status='replace')
    ! Uses port 110
    cycle
  end if

  if (command == 'END                          :') exit
end do

! Closing the reading
close(1)

!=================================================================
!Appel des differentes ressources lues dans le fichier de commande
!=================================================================
do

  call read_Vloc_dof_bodies

  if (c_draw                              == 1) call draw
  if (c_toppling_v                        == 1) call toppling_v

  if (fin_post) then
    call close_all
    print*,'--> End of post-processing <--'
    stop
  endif

  if(first_over_all) first_over_all = .false.

end do

!================================================
!================================================
!================================================
!================================================

 contains

!======================================================================
!======================================================================
SUBROUTINE read_Vloc_dof_bodies

  implicit none

  integer                              ::  i, j
  integer                              ::  n_prpr,n_prpl,num_part,err,cd,an,n_vertices,n_faces
  integer                              ::  n_PRPRx, n_PRPLx, icdent,ianent,nb_ctc_deja_compte,cpt,cpt_actif
  integer                              ::  icdver
  real(kind=8)                         ::  rn,rs,rt,t1,t2,t3,n1,n2,n3,s1,s2,s3,vls,vln,vlt, igap
  real(kind=8)                         ::  I1,I2,I3
  character(len=6)                     ::  text
  character(len=5)                     ::  status, law
  character(len=5)                     ::  color,behav
  character(len=13)                    ::  text2
  real(kind=8)                         ::  center1,center2,center3,ax1r,ax2r,ax3r,coor1,coor2,coor3
  real(kind=8)                         ::  coor_ctc1,coor_ctc2,coor_ctc3
  real(kind=8)                         ::  ver1,ver2,ver3,Volume_tetra,Air_triangle,demi_somme
  integer                              ::  face1,face2,face3
  real(kind=8),dimension(3)            ::  X,Y,Z,centre_gravity,ab,ac,bc
  real(kind=8),dimension(3,3)          ::  rot
  real(kind=8)                         ::  Vrx,Vry,Vrz,vx,vy,vz,dcd_ptc_contact,dan_ptc_contact,&
                                           n_ff,n_fs,n_fv,n_ss
  character(len=22)                    ::  clout_DOF
  character(len=28)                    ::  clout_Vloc
  character(len=20)                    ::  clout_Bodies
  real(kind=8)                         ::  Number_mean_Total_force=0, total_p_volume

  real(kind=8)                         ::  curr_rad
  integer                              ::  temp_i_var, total_faces
  real(kind=8)                         ::  avr_part_rad


  clout_Bodies = './OUTBOX/BODIES.OUT'
  clout_DOF  =  './OUTBOX/DOF.OUT.     '
  clout_Vloc =  './OUTBOX/Vloc_Rloc.OUT.     '

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
  else if ( (compteur_clout>=1000).and. (compteur_clout<10000) ) then
    WRITE(clout_DOF(18:22),'(I4)')   compteur_clout
    WRITE(clout_Vloc(24:28),'(I4)')  compteur_clout
  end if

  ! Reading the BODIES.OUT file
  if (first_over_all) then
    open(unit=2,file=clout_Bodies,status='old')
    print*,'--->',clout_Bodies

    temp_i_var = 0
    total_faces = 0

    allocate(TAB_POLY(n_particles))
    allocate(TAB_PLAN(n_walls))
    i=0
    do
      read(2,'(A6)') text
      if (text == '$bdyty') then
        read(2,*)
        read(2,*)
        read(2,'(22X, A5,7X,D14.7)') behav, curr_rad

        !print*, i
        ! Storing the particles information
        if (behav == 'PLEXx') then
          i=i+1
          TAB_POLY(i)%id=i

          ! Computing the average equivalent particle radius
          avr_part_rad = avr_part_rad + curr_rad

          TAB_POLY(i)%behav = behav
          read(2,'(29X, 3(5x,D14.7,2X))') I1,I2,I3
          TAB_POLY(i)%I1 = I1
          TAB_POLY(i)%I2 = I2
          TAB_POLY(i)%I3 = I3
          read(2,*)
          read(2,'(29X, 3(5x,D14.7,2X))') center1,center2,center3
          TAB_POLY(i)%center_ref(1)=center1
          TAB_POLY(i)%center_ref(2)=center2
          TAB_POLY(i)%center_ref(3)=center3
          read(2,*)
          read(2,*)
          read(2,'(22X,A5,12X,i7,13X,i7)') color,n_vertices,n_faces

          TAB_POLY(i)%color     = color
          TAB_POLY(i)%n_faces   = n_faces
          TAB_POLY(i)%nb_vertex = n_vertices
          allocate(TAB_POLY(i)%vertex_ref(3,n_vertices))
          allocate(TAB_POLY(i)%vertex(3,n_vertices))
          allocate(TAB_POLY(i)%face(3,n_faces))

          ! Counting the number of faces
          total_faces = total_faces + n_faces

          ! Reading the coordinates of the vertices
          do j=1,TAB_POLY(i)%nb_vertex
            read(2,'(29X, 3(5x,D14.7,2X))') ver1,ver2,ver3
            TAB_POLY(i)%vertex_ref(1,j) = ver1
            TAB_POLY(i)%vertex_ref(2,j) = ver2
            TAB_POLY(i)%vertex_ref(3,j) = ver3
          end do

          ! Finding the centroid of the polyhedron
          centre_gravity = 0
          do j=1,TAB_POLY(i)%nb_vertex
            centre_gravity(1)=centre_gravity(1)+TAB_POLY(i)%vertex_ref(1,j)
            centre_gravity(2)=centre_gravity(2)+TAB_POLY(i)%vertex_ref(2,j)
            centre_gravity(3)=centre_gravity(3)+TAB_POLY(i)%vertex_ref(3,j)
          enddo
          centre_gravity=centre_gravity/real(TAB_POLY(i)%nb_vertex,8)

          ! Finding the farthest vertex from the particle centroid
          TAB_POLY(i)%Rmax = 0
          do j=1,TAB_POLY(i)%nb_vertex
            X(1) = TAB_POLY(i)%vertex_ref(1,j) - centre_gravity(1)
            X(2) = TAB_POLY(i)%vertex_ref(2,j) - centre_gravity(2)
            X(3) = TAB_POLY(i)%vertex_ref(3,j) - centre_gravity(3)

            TAB_POLY(i)%Rmax = max(TAB_POLY(i)%Rmax , sqrt(X(1)**2 + X(2)**2 + X(3)**2))
          end do

          ! Finding the volume of the particle
          Volume_tetra = 0.D0
          Air_triangle = 0.D0
          do j=1,TAB_POLY(i)%n_faces
            read(2,'(29X, 3(5x,i7,9X))') face1,face2,face3
            TAB_POLY(i)%face(1,j) = face1
            TAB_POLY(i)%face(2,j) = face2
            TAB_POLY(i)%face(3,j) = face3
            X(1) = TAB_POLY(i)%vertex_ref(1,face1) - centre_gravity(1)
            X(2) = TAB_POLY(i)%vertex_ref(2,face1) - centre_gravity(2)
            X(3) = TAB_POLY(i)%vertex_ref(3,face1) - centre_gravity(3)

            Y(1) = TAB_POLY(i)%vertex_ref(1,face2) - centre_gravity(1)
            Y(2) = TAB_POLY(i)%vertex_ref(2,face2) - centre_gravity(2)
            Y(3) = TAB_POLY(i)%vertex_ref(3,face2) - centre_gravity(3)

            Z(1) = TAB_POLY(i)%vertex_ref(1,face3) - centre_gravity(1)
            Z(2) = TAB_POLY(i)%vertex_ref(2,face3) - centre_gravity(2)
            Z(3) = TAB_POLY(i)%vertex_ref(3,face3) - centre_gravity(3)

            ! Orignal
            !Volume_tetra = Volume_tetra + abs( X(1)*(Y(2)*Z(3)-Y(3)*Z(2)) &
            !                                  -X(2)*(Y(1)*Z(3)-Y(3)*Z(1)) &
            !                                  +X(3)*(Y(1)*Z(2)-Y(2)*Z(1)))
            ! New
            Volume_tetra = Volume_tetra + abs( X(1)*(Y(2)*Z(3)-Y(3)*Z(2)) &
                                              +X(2)*(Y(3)*Z(1)-Y(1)*Z(3)) &
                                              +X(3)*(Y(1)*Z(2)-Y(2)*Z(1)))

            ab(:) = X(:)-Y(:)
            ac(:) = X(:)-Z(:)
            bc(:) = Y(:)-Z(:)
            demi_somme = 0.5*( sqrt( ab(1)**2+ab(2)**2+ab(3)**2 ) +&
                               sqrt( ac(1)**2+ac(2)**2+ac(3)**2 ) +&
                               sqrt( bc(1)**2+bc(2)**2+bc(3)**2 ))

            Air_triangle = Air_triangle + sqrt(demi_somme * (demi_somme - sqrt(ab(1)**2+ab(2)**2+ab(3)**2))*&
                                                            (demi_somme - sqrt(ac(1)**2+ac(2)**2+ac(3)**2))*&
                                                            (demi_somme - sqrt(bc(1)**2+bc(2)**2+bc(3)**2)))
          end do

          TAB_POLY(i)%volume     = Volume_tetra/6.D0

          ! Computing the total volume in the box
          total_p_volume = total_p_volume + TAB_POLY(i)%volume

          TAB_POLY(i)%Aire       = Air_triangle

          TAB_POLY(i)%center(1)=0.D0
          TAB_POLY(i)%momentI=0
          TAB_POLY(i)%nctc = 0
          TAB_POLY(i)%nctcslide = 0
          TAB_POLY(i)%nctcstick = 0
          TAB_POLY(i)%nctcs = 0
          TAB_POLY(i)%nctcd = 0
          TAB_POLY(i)%nctct = 0
          TAB_POLY(i)%nb_ctc = 0
          TAB_POLY(i)%force = 0.D0

        !--Storing the walls information
        else
          i=i+1
          TAB_PLAN(1)%id=1
          read(2,*)
          read(2,*)
          read(2,*)
          read(2,*)
          read(2,*)
          read(2,'(29X, 3(5x,D14.7,2X))') center1,center2,center3
          TAB_PLAN(1)%center_ref(1)=center1
          TAB_PLAN(1)%center_ref(2)=center2
          TAB_PLAN(1)%center_ref(3)=center3
          read(2,*)
          read(2,*)
          read(2,'(29X, 3(5x,D14.7,2X))') ax1r,ax2r,ax3r
          TAB_PLAN(1)%ax1=ax1r
          TAB_PLAN(1)%ax2=ax2r
          TAB_PLAN(1)%ax3=ax3r
          read(2,*)
          TAB_PLAN(1)%center(:)=0.D0
          TAB_PLAN(1)%behav = 'WALL_'

          ! Similar to poly, i'm going to create vertices (8 vertices in each wall)
          allocate(TAB_PLAN(1)%vertex_ref(3,8))
          allocate(TAB_PLAN(1)%vertex(3,8))

          ! First sheet
          TAB_PLAN(1)%vertex_ref(1,1) =  ax1r
          TAB_PLAN(1)%vertex_ref(2,1) =  ax2r
          TAB_PLAN(1)%vertex_ref(3,1) =  ax3r

          TAB_PLAN(1)%vertex_ref(1,2) = (- ax1r)
          TAB_PLAN(1)%vertex_ref(2,2) =  ax2r
          TAB_PLAN(1)%vertex_ref(3,2) =  ax3r

          TAB_PLAN(1)%vertex_ref(1,3) = (- ax1r)
          TAB_PLAN(1)%vertex_ref(2,3) = (- ax2r)
          TAB_PLAN(1)%vertex_ref(3,3) =  ax3r

          TAB_PLAN(1)%vertex_ref(1,4) =  ax1r
          TAB_PLAN(1)%vertex_ref(2,4) = (- ax2r)
          TAB_PLAN(1)%vertex_ref(3,4) =  ax3r

          ! Second sheet
          TAB_PLAN(1)%vertex_ref(1,5) =  ax1r
          TAB_PLAN(1)%vertex_ref(2,5) =  ax2r
          TAB_PLAN(1)%vertex_ref(3,5) = (- ax3r)

          TAB_PLAN(1)%vertex_ref(1,6) = (- ax1r)
          TAB_PLAN(1)%vertex_ref(2,6) =  ax2r
          TAB_PLAN(1)%vertex_ref(3,6) = (- ax3r)

          TAB_PLAN(1)%vertex_ref(1,7) = (- ax1r)
          TAB_PLAN(1)%vertex_ref(2,7) = (- ax2r)
          TAB_PLAN(1)%vertex_ref(3,7) = (- ax3r)

          TAB_PLAN(1)%vertex_ref(1,8) =  ax1r
          TAB_PLAN(1)%vertex_ref(2,8) = (- ax2r)
          TAB_PLAN(1)%vertex_ref(3,8) = (- ax3r)

        end if
        if (i==n_particles+n_walls) exit
      end if
    end do
    ! Printing the average particle radius
    print*, 'Average radius of all particles: ', avr_part_rad/n_particles
  end if
  close(2)

  !-Reading the DOF.OUT.#
  open(unit=5,file=clout_DOF,iostat=err,status='old')
  if (err/=0) then
    fin_post=.true.
  else
    print*,'--->',clout_DOF
    fin_post=.false.
    i = 0
    read(5,*)
    read(5,*)
    read(5,*)
    read(5,'(7X,i9,19X,D14.7)') step,time

    do
      read(5,'(A6)') text
      if (text == '$bdyty') then
        i = i+1
        read(5,'(6X,i9)') num_part

        !--Storing the particles information
        if (i<n_particles+1) then
          read(5,*)
          read(5,'(29X, 3(5x,D14.7,2X))') coor1,coor2,coor3
          TAB_POLY(i)%center(1)=coor1+TAB_POLY(i)%center_ref(1)
          TAB_POLY(i)%center(2)=coor2+TAB_POLY(i)%center_ref(2)
          TAB_POLY(i)%center(3)=coor3+TAB_POLY(i)%center_ref(3)
          read(5,*)
          read(5,'(29X, 3(5x,D14.7,2X))') Vx,Vy,Vz
          read(5,'(29X, 3(5x,D14.7,2X))') Vrx,Vry,Vrz

          TAB_POLY(i)%Vx=Vx
          TAB_POLY(i)%Vy=Vy
          TAB_POLY(i)%Vz=Vz
          TAB_POLY(i)%Vrx=Vrx
          TAB_POLY(i)%Vry=Vry
          TAB_POLY(i)%Vrz=Vrz
          TAB_POLY(i)%Norme_V  = sqrt(Vx**2+Vy**2+Vz**2)
          TAB_POLY(i)%Norme_Vr = sqrt(Vrx**2+Vry**2+Vrz**2)
          read(5,'(29X, 3(5x,D14.7,2X))') rot(1,1),rot(2,1),rot(3,1)
          read(5,'(29X, 3(5x,D14.7,2X))') rot(1,2),rot(2,2),rot(3,2)
          read(5,'(29X, 3(5x,D14.7,2X))') rot(1,3),rot(2,3),rot(3,3)
          TAB_POLY(i)%Rot = rot
          TAB_POLY(i)%Wx = Vrx*rot(1,1) + Vry*rot(1,2) + Vrz*rot(1,3)
          TAB_POLY(i)%Wy = Vrx*rot(2,1) + Vry*rot(2,2) + Vrz*rot(2,3)
          TAB_POLY(i)%Wz = Vrx*rot(3,1) + Vry*rot(3,2) + Vrz*rot(3,3)

          ! On rotationne et translate
          do j=1,TAB_POLY(i)%nb_vertex
            ! histoire de tout avoir dans le repère global
            TAB_POLY(i)%vertex(1,j)=TAB_POLY(i)%vertex_ref(1,j)*rot(1,1)+&
                                    TAB_POLY(i)%vertex_ref(2,j)*rot(1,2)+&
                                    TAB_POLY(i)%vertex_ref(3,j)*rot(1,3)+&
                                    TAB_POLY(i)%center(1)

            TAB_POLY(i)%vertex(2,j)=TAB_POLY(i)%vertex_ref(1,j)*rot(2,1)+&
                                    TAB_POLY(i)%vertex_ref(2,j)*rot(2,2)+&
                                    TAB_POLY(i)%vertex_ref(3,j)*rot(2,3)+&
                                    TAB_POLY(i)%center(2)

            TAB_POLY(i)%vertex(3,j)=TAB_POLY(i)%vertex_ref(1,j)*rot(3,1)+&
                                    TAB_POLY(i)%vertex_ref(2,j)*rot(3,2)+&
                                    TAB_POLY(i)%vertex_ref(3,j)*rot(3,3)+&
                                    TAB_POLY(i)%center(3)
          end do

          TAB_POLY(i)%nctc = 0
          TAB_POLY(i)%nctcslide = 0
          TAB_POLY(i)%nctcstick = 0
          TAB_POLY(i)%nctcs = 0
          TAB_POLY(i)%nctcd = 0
          TAB_POLY(i)%nctct = 0
          TAB_POLY(i)%nb_ctc = 0
          TAB_POLY(i)%force = 0.D0

        !--Storing the walls informations
        else
          read(5,*)
          read(5,'(29X, 3(5x,D14.7,2X))') coor1,coor2,coor3
          TAB_PLAN(i-n_particles)%center(1)=coor1+TAB_PLAN(i-n_particles)%center_ref(1)
          TAB_PLAN(i-n_particles)%center(2)=coor2+TAB_PLAN(i-n_particles)%center_ref(2)
          TAB_PLAN(i-n_particles)%center(3)=coor3+TAB_PLAN(i-n_particles)%center_ref(3)
          read(5,*)
          read(5,'(29X, 3(5x,D14.7,2X))') Vx,Vy,Vz
          read(5,'(29X, 3(5x,D14.7,2X))') Vrx,Vry,Vrz
          TAB_PLAN(i-n_particles)%Vx = Vx
          TAB_PLAN(i-n_particles)%Vy = Vy
          TAB_PLAN(i-n_particles)%Vz = Vz
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
          TAB_PLAN(i-n_particles)%nctcslide = 0
          TAB_PLAN(i-n_particles)%nctcstick = 0
          TAB_PLAN(i-n_particles)%nctcs = 0
          TAB_PLAN(i-n_particles)%nctcd = 0
          TAB_PLAN(i-n_particles)%nctct = 0
          TAB_PLAN(i-n_particles)%nb_ctc = 0
          TAB_PLAN(i-n_particles)%force = 0.D0

          ! Translating and rotating the vertices
          ! There are 8 vertices per face
          do j=1, 8
            TAB_PLAN(i-n_particles)%vertex(1,j) = TAB_PLAN(i-n_particles)%vertex_ref(1,j)*rot(1,1)+&
                                                  TAB_PLAN(i-n_particles)%vertex_ref(2,j)*rot(1,2)+&
                                                  TAB_PLAN(i-n_particles)%vertex_ref(3,j)*rot(1,3)+&
                                                  TAB_PLAN(i-n_particles)%center(1)

            TAB_PLAN(i-n_particles)%vertex(2,j) = TAB_PLAN(i-n_particles)%vertex_ref(1,j)*rot(2,1)+&
                                                  TAB_PLAN(i-n_particles)%vertex_ref(2,j)*rot(2,2)+&
                                                  TAB_PLAN(i-n_particles)%vertex_ref(3,j)*rot(2,3)+&
                                                  TAB_PLAN(i-n_particles)%center(2)

            TAB_PLAN(i-n_particles)%vertex(3,j) = TAB_PLAN(i-n_particles)%vertex_ref(1,j)*rot(3,1)+&
                                                  TAB_PLAN(i-n_particles)%vertex_ref(2,j)*rot(3,2)+&
                                                  TAB_PLAN(i-n_particles)%vertex_ref(3,j)*rot(3,3)+&
                                                  TAB_PLAN(i-n_particles)%center(3)
          end do
        end if
      end if
      if (i==n_particles+n_walls) exit
    end do
  end if

  ! Variables to count the number of PRPR and PRPL contacts
  n_PRPRx=0
  n_PRPLx=0
  !-Reading the Vloc_Rloc files
  !-- Finding the number and type of contacts
  open(unit=4,file=clout_Vloc,iostat=err,status='old')
  if (err/=0) then
    ! - Void
  else

    ! Reading
    print*,'--->',clout_Vloc
    read(4,*) 
    read(4,*) 
    read(4,*) 
    read(4,*) 
    read(4,*) 
    read(4,*) 
    do
      ! Reading the type of contact
      read(4,'(A13)',iostat=err) text2
      if (err/=0) then
        close(2)
        close(4)
        close(5)
        exit
      end if

      if (text2 == '!-------------') exit                ! Me parece que hay que quitarlo
      if (text2 == '$icdan  PRPLx') then
        n_PRPLx = n_PRPLx +1
      end if
      if (text2 == '$icdan  PRPRx') then
        n_PRPRx = n_PRPRx +1
      end if
    end do
    close(4)

    ! Storing global variables
    nf_PRPLx = n_PRPLx
    nf_PRPRx = n_PRPRx

    ! Allocating the container for the contacts
    if (allocated(TAB_CONTACT)) deallocate(TAB_CONTACT)
    allocate(TAB_CONTACT(n_PRPRx+n_PRPLx+1))

    ! Initializing the container
    do i=1,(n_PRPRx+n_PRPLx)
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
      TAB_CONTACT(i)%vln=0.D0
      TAB_CONTACT(i)%vlt=0.D0
      TAB_CONTACT(i)%vls=0.D0
      TAB_CONTACT(i)%type=0
      TAB_CONTACT(i)%nature = '    '
    end do

    ! The Vloc file is re-opened to store the information
    n_prpr = 0
    n_prpl = 0
    nb_ligneCONTACT=0

    open(unit=4,file=clout_Vloc,status='old')
    read(4,*) 
    read(4,*) 
    read(4,*) 
    read(4,*) 
    read(4,*) 
    read(4,*) 

    !--Storing the contacts information
    do
      read(4,'(A13)',iostat=err) text2

      if (text2 == '$icdan  PRPRx') then
        n_prpr = n_prpr + 1
        read(4,*)
        read(4,'(7X,i6,15X,i6,2X,A5,8X,i6,16X,A5)') icdent, icdver, law, ianent, status
        read(4,'(29X, 3(5x,D14.7,2X))') rt,rn,rs
        read(4,'(29X, 3(5x,D14.7,2X))') vlt,vln,vls
        read(4,'(29X, 1(5x,D14.7,2X))') igap
        read(4,'(29X, 3(5x,D14.7,2X))') coor_ctc1,coor_ctc2,coor_ctc3
        read(4,'(29X, 3(5x,D14.7,2X))') t1,t2,t3
        read(4,'(29X, 3(5x,D14.7,2X))') n1,n2,n3
        read(4,'(29X, 3(5x,D14.7,2X))') s1,s2,s3

        nb_ligneCONTACT=nb_ligneCONTACT+1

        TAB_CONTACT(n_prpr)%icdent=icdent
        TAB_CONTACT(n_prpr)%cdver=icdver
        TAB_CONTACT(n_prpr)%i_law=law
        TAB_CONTACT(n_prpr)%ianent=ianent
        TAB_CONTACT(n_prpr)%status=status
        TAB_CONTACT(n_prpr)%gapTT=igap
        TAB_CONTACT(n_prpr)%n(1)=n1
        TAB_CONTACT(n_prpr)%n(2)=n2
        TAB_CONTACT(n_prpr)%n(3)=n3
        TAB_CONTACT(n_prpr)%t(1)=t1
        TAB_CONTACT(n_prpr)%t(2)=t2
        TAB_CONTACT(n_prpr)%t(3)=t3
        TAB_CONTACT(n_prpr)%s(1)=s1
        TAB_CONTACT(n_prpr)%s(2)=s2
        TAB_CONTACT(n_prpr)%s(3)=s3
        TAB_CONTACT(n_prpr)%coor_ctc(1)=coor_ctc1
        TAB_CONTACT(n_prpr)%coor_ctc(2)=coor_ctc2
        TAB_CONTACT(n_prpr)%coor_ctc(3)=coor_ctc3
        TAB_CONTACT(n_prpr)%rn=rn
        TAB_CONTACT(n_prpr)%rt=rt
        TAB_CONTACT(n_prpr)%rs=rs
        TAB_CONTACT(n_prpr)%vln=vln
        TAB_CONTACT(n_prpr)%vlt=vlt
        TAB_CONTACT(n_prpr)%vls=vls
        TAB_CONTACT(n_prpr)%type=1
        TAB_CONTACT(n_prpr)%nature = 'PRPRx'
        TAB_CONTACT(n_prpr)%deja_compte = .false.

    ! If the contact occurs between two polyhedra
      else if (text2 == '$icdan  PRPLx') then
        n_prpl = n_prpl + 1
        read(4,*)
        read(4,'(7X,i6,15X,i6,2X,A5,8X,i6,16X,A5)') icdent, icdver, law, ianent, status
        read(4,'(29X, 3(5x,D14.7,2X))') rs,rt,rn
        read(4,'(29X, 3(5x,D14.7,2X))') vls,vlt,vln
        read(4,'(29X, 1(5x,D14.7,2X))') igap
        read(4,'(29X, 3(5x,D14.7,2X))') coor_ctc1,coor_ctc2,coor_ctc3
        read(4,'(29X, 3(5x,D14.7,2X))') t1,t2,t3
        read(4,'(29X, 3(5x,D14.7,2X))') n1,n2,n3
        read(4,'(29X, 3(5x,D14.7,2X))') s1,s2,s3

        nb_ligneCONTACT=nb_ligneCONTACT+1

        TAB_CONTACT(n_prpl+n_PRPRx)%icdent=icdent
        TAB_CONTACT(n_prpl+n_PRPRx)%cdver=icdver
        TAB_CONTACT(n_prpl+n_PRPRx)%i_law=law
        TAB_CONTACT(n_prpl+n_PRPRx)%ianent=ianent
        TAB_CONTACT(n_prpl+n_PRPRx)%status=status
        TAB_CONTACT(n_prpl+n_PRPRx)%gapTT=igap
        TAB_CONTACT(n_prpl+n_PRPRx)%n(1)=n1
        TAB_CONTACT(n_prpl+n_PRPRx)%n(2)=n2
        TAB_CONTACT(n_prpl+n_PRPRx)%n(3)=n3
        TAB_CONTACT(n_prpl+n_PRPRx)%t(1)=t1
        TAB_CONTACT(n_prpl+n_PRPRx)%t(2)=t2
        TAB_CONTACT(n_prpl+n_PRPRx)%t(3)=t3
        TAB_CONTACT(n_prpl+n_PRPRx)%s(1)=s1
        TAB_CONTACT(n_prpl+n_PRPRx)%s(2)=s2
        TAB_CONTACT(n_prpl+n_PRPRx)%s(3)=s3
        TAB_CONTACT(n_prpl+n_PRPRx)%coor_ctc(1)=coor_ctc1
        TAB_CONTACT(n_prpl+n_PRPRx)%coor_ctc(2)=coor_ctc2
        TAB_CONTACT(n_prpl+n_PRPRx)%coor_ctc(3)=coor_ctc3
        TAB_CONTACT(n_prpl+n_PRPRx)%rn=rn
        TAB_CONTACT(n_prpl+n_PRPRx)%rt=rt
        TAB_CONTACT(n_prpl+n_PRPRx)%rs=rs
        TAB_CONTACT(n_prpl+n_PRPRx)%vln=vln
        TAB_CONTACT(n_prpl+n_PRPRx)%vlt=vlt
        TAB_CONTACT(n_prpl+n_PRPRx)%vls=vls
        TAB_CONTACT(n_prpl+n_PRPRx)%type=1
        TAB_CONTACT(n_prpl+n_PRPRx)%nature = 'PRPLx'
        TAB_CONTACT(n_prpl+n_PRPRx)%deja_compte = .false.
      end if

      if (n_prpl+n_prpr == n_PRPLx + n_PRPRx) then
        exit
      end if
    end do
  endif

  ! It is necessary to build a new summarized contact array
  nb_ctc_deja_compte = 0
  do i=1,nb_ligneCONTACT
    if ( (TAB_CONTACT(i)%deja_compte) ) cycle

    cd = TAB_CONTACT(i)%icdent
    an = TAB_CONTACT(i)%ianent

    do j=i+1,nb_ligneCONTACT
      if (TAB_CONTACT(i)%type==4) exit ! car le max se sont des contacts quadruples
      if (j-i > 5) exit     ! car en fait les contacts face/face sont ecrit les un a la suite des autres
      if ((cd == TAB_CONTACT(j)%icdent).and.(an==TAB_CONTACT(j)%ianent).or.&
          (cd == TAB_CONTACT(j)%ianent).and.(an==TAB_CONTACT(j)%icdent) ) then

        TAB_CONTACT(j)%deja_compte = .true.
        nb_ctc_deja_compte = nb_ctc_deja_compte + 1

        ! NOOOOOOOO PUTAAAAAA
        !if ((TAB_CONTACT(i)%rn>0.D0).and.(TAB_CONTACT(j)%rn>0.D0)) then
        !TAB_CONTACT(i)%rn = TAB_CONTACT(i)%rn + TAB_CONTACT(j)%rn
        !TAB_CONTACT(i)%rt = TAB_CONTACT(i)%rt + TAB_CONTACT(j)%rt
        !TAB_CONTACT(i)%rs = TAB_CONTACT(i)%rs + TAB_CONTACT(j)%rs

        !TAB_CONTACT(i)%vln = TAB_CONTACT(i)%vln + TAB_CONTACT(j)%vln
        !TAB_CONTACT(i)%vlt = TAB_CONTACT(i)%vlt + TAB_CONTACT(j)%vlt
        !TAB_CONTACT(i)%vls = TAB_CONTACT(i)%vls + TAB_CONTACT(j)%vls

        !TAB_CONTACT(i)%coor_ctc(:) = TAB_CONTACT(i)%coor_ctc(:) + TAB_CONTACT(j)%coor_ctc(:)

        TAB_CONTACT(i)%type = TAB_CONTACT(i)%type + 1
        TAB_CONTACT(j)%type = TAB_CONTACT(i)%type

        !else if ((TAB_CONTACT(i)%rn==0.D0).and.(TAB_CONTACT(j)%rn>0.D0)) then
        !  TAB_CONTACT(i)%rn = TAB_CONTACT(j)%rn
        !  TAB_CONTACT(i)%rt = TAB_CONTACT(j)%rt
        !  TAB_CONTACT(i)%rs = TAB_CONTACT(j)%rs
        !
        !  TAB_CONTACT(i)%vln = TAB_CONTACT(j)%vln
        !  TAB_CONTACT(i)%vlt = TAB_CONTACT(j)%vlt
        !  TAB_CONTACT(i)%vls = TAB_CONTACT(j)%vls
        !
        !  TAB_CONTACT(i)%coor_ctc(:) = TAB_CONTACT(j)%coor_ctc(:)
        !
        !  TAB_CONTACT(i)%type = TAB_CONTACT(j)%type
        !end if
      end if
    end do
  end do

  nb_ligneCONTACT_POLYR = nb_ligneCONTACT - nb_ctc_deja_compte
  cpt = 0
  cpt_actif = 0
  n_ff = 0
  n_fs = 0
  n_fv = 0
  n_ss = 0

  if (allocated(TAB_CONTACT_POLYR)) deallocate(TAB_CONTACT_POLYR)
  allocate(TAB_CONTACT_POLYR(nb_ligneCONTACT_POLYR))

  TAB_CONTACT_POLYR(:)%x_status = 0.0

  do i=1,nb_ligneCONTACT
    if ((TAB_CONTACT(i)%deja_compte)) cycle
    cpt = cpt + 1
    dcd_ptc_contact = 1000000000.D0
    dan_ptc_contact = 1000000000.D0

    ! Your data
    TAB_CONTACT_POLYR(cpt)%loc      = i
    TAB_CONTACT_POLYR(cpt)%icdent   = TAB_CONTACT(i)%icdent
    TAB_CONTACT_POLYR(cpt)%ianent   = TAB_CONTACT(i)%ianent
    TAB_CONTACT_POLYR(cpt)%n(:)     = TAB_CONTACT(i)%n(:)
    TAB_CONTACT_POLYR(cpt)%t(:)     = TAB_CONTACT(i)%t(:)
    TAB_CONTACT_POLYR(cpt)%s(:)     = TAB_CONTACT(i)%s(:)
    TAB_CONTACT_POLYR(cpt)%coor_ctc = TAB_CONTACT(i)%coor_ctc
    TAB_CONTACT_POLYR(cpt)%rn       = TAB_CONTACT(i)%rn
    TAB_CONTACT_POLYR(cpt)%rt       = TAB_CONTACT(i)%rt
    TAB_CONTACT_POLYR(cpt)%rs       = TAB_CONTACT(i)%rs
    TAB_CONTACT_POLYR(cpt)%vln      = TAB_CONTACT(i)%vln
    TAB_CONTACT_POLYR(cpt)%vlt      = TAB_CONTACT(i)%vlt
    TAB_CONTACT_POLYR(cpt)%vls      = TAB_CONTACT(i)%vls
    TAB_CONTACT_POLYR(cpt)%type     = TAB_CONTACT(i)%type
    TAB_CONTACT_POLYR(cpt)%nature   = TAB_CONTACT(i)%nature

    cd = TAB_CONTACT_POLYR(cpt)%icdent
    an = TAB_CONTACT_POLYR(cpt)%ianent

    ! Data of your friends
    do j=i+1, nb_ligneCONTACT
      if (j-i > 5) exit
      if ((cd == TAB_CONTACT(j)%icdent).and.(an==TAB_CONTACT(j)%ianent) .or. &
        (cd == TAB_CONTACT(j)%ianent).and.(an==TAB_CONTACT(j)%icdent) ) then

        TAB_CONTACT_POLYR(cpt)%coor_ctc = TAB_CONTACT_POLYR(cpt)%coor_ctc + TAB_CONTACT(j)%coor_ctc
        TAB_CONTACT_POLYR(cpt)%rn       = TAB_CONTACT_POLYR(cpt)%rn     + TAB_CONTACT(j)%rn
        TAB_CONTACT_POLYR(cpt)%rt       = TAB_CONTACT_POLYR(cpt)%rt     + TAB_CONTACT(j)%rt
        TAB_CONTACT_POLYR(cpt)%rs       = TAB_CONTACT_POLYR(cpt)%rs     + TAB_CONTACT(j)%rs
        TAB_CONTACT_POLYR(cpt)%vln      = TAB_CONTACT_POLYR(cpt)%vln    + TAB_CONTACT(j)%vln
        TAB_CONTACT_POLYR(cpt)%vlt      = TAB_CONTACT_POLYR(cpt)%vlt    + TAB_CONTACT(j)%vlt
        TAB_CONTACT_POLYR(cpt)%vls      = TAB_CONTACT_POLYR(cpt)%vls    + TAB_CONTACT(j)%vls
      end if
    end do

    ! Finally
    TAB_CONTACT_POLYR(cpt)%coor_ctc = TAB_CONTACT_POLYR(cpt)%coor_ctc/TAB_CONTACT(i)%type
    TAB_CONTACT_POLYR(cpt)%vln = TAB_CONTACT_POLYR(cpt)%vln/TAB_CONTACT(i)%type
    TAB_CONTACT_POLYR(cpt)%vlt = TAB_CONTACT_POLYR(cpt)%vlt/TAB_CONTACT(i)%type
    TAB_CONTACT_POLYR(cpt)%vls = TAB_CONTACT_POLYR(cpt)%vls/TAB_CONTACT(i)%type

    if (TAB_CONTACT_POLYR(cpt)%type == 2) n_fs = n_fs + 1
    if (TAB_CONTACT_POLYR(cpt)%type == 3) n_ff = n_ff + 1
    if (TAB_CONTACT_POLYR(cpt)%type == 4) n_ff = n_ff + 1

    if (TAB_CONTACT_POLYR(cpt)%type==1) then  ! Il faut distinguer les contacts face-vertex des contact side-side
      cd = TAB_CONTACT_POLYR(cpt)%icdent
      an = TAB_CONTACT_POLYR(cpt)%ianent

      if (cd>n_particles.or.an>n_particles) then
        n_fv = n_fv + 1
        TAB_CONTACT_POLYR(cpt)%type = 1
      else
        do j=1,TAB_POLY(cd)%nb_vertex
          dcd_ptc_contact = min(dcd_ptc_contact, sqrt((TAB_CONTACT_POLYR(cpt)%coor_ctc(1)-TAB_POLY(cd)%vertex(1,j))**2+&
                                                      (TAB_CONTACT_POLYR(cpt)%coor_ctc(2)-TAB_POLY(cd)%vertex(2,j))**2+&
                                                      (TAB_CONTACT_POLYR(cpt)%coor_ctc(3)-TAB_POLY(cd)%vertex(3,j))**2))
        end do
        do j=1,TAB_POLY(an)%nb_vertex
          dan_ptc_contact = min(dan_ptc_contact, sqrt((TAB_CONTACT_POLYR(cpt)%coor_ctc(1)-TAB_POLY(an)%vertex(1,j))**2+&
                                                      (TAB_CONTACT_POLYR(cpt)%coor_ctc(2)-TAB_POLY(an)%vertex(2,j))**2+&
                                                      (TAB_CONTACT_POLYR(cpt)%coor_ctc(3)-TAB_POLY(an)%vertex(3,j))**2))
        end do

        if (min( dcd_ptc_contact , dan_ptc_contact ) < 0.0001) then
          n_fv = n_fv + 1
          TAB_CONTACT_POLYR(cpt)%type = 1
        else
          n_ss = n_ss + 1
          TAB_CONTACT_POLYR(cpt)%type = 0
        end if
      end if
    end if

    if (TAB_CONTACT_POLYR(cpt)%vln**2 + TAB_CONTACT_POLYR(cpt)%vlt**2 > 0.000000001) then
      TAB_CONTACT_POLYR(cpt)%status='slide'
    end if

    if (TAB_CONTACT_POLYR(cpt)%vln**2 + TAB_CONTACT_POLYR(cpt)%vlt**2 < 0.000000001) then
      TAB_CONTACT_POLYR(cpt)%status='stick'
    end if

    if (TAB_CONTACT_POLYR(cpt)%rn>0.D0)  cpt_actif = cpt_actif +1
  end do

  Mean_total_Normal_force = 0.D0
  Number_mean_Total_force = 0.D0
  do i=1,nb_ligneCONTACT_POLYR
    if (TAB_CONTACT_POLYR(i)%rn == 0.D0 ) cycle
    Mean_total_Normal_force = Mean_total_Normal_force + TAB_CONTACT(i)%rn
    Number_mean_Total_force = Number_mean_Total_force + 1.D0
  end do
  Mean_total_Normal_force = Mean_total_Normal_force / Number_mean_Total_force

  ! Computing the number of active contacts per particle
  TAB_POLY(:)%nctc = 0.

  do i=1,nb_ligneCONTACT_POLYR
    icdent = TAB_CONTACT_POLYR(i)%icdent
    ianent = TAB_CONTACT_POLYR(i)%ianent
    ! Contacts between particles
    if (TAB_CONTACT_POLYR(i)%nature /= 'PRPRx') cycle

    ! Only active contacts
    if (TAB_CONTACT_POLYR(i)%rn .le. 0.) cycle

    TAB_POLY(icdent)%nctc = TAB_POLY(icdent)%nctc + 1
    TAB_POLY(ianent)%nctc = TAB_POLY(ianent)%nctc + 1
  end do

  print*, ''
  print*, ''
  print*,'Step:', step
  print*,'Time:', time
  print*,'Nb particles: ',n_particles
  print*,'Nb Contacts: ',nb_ligneCONTACT
  print*,'Nb real Contacts: ',nb_ligneCONTACT_POLYR
  print*,'Nb active Contacts: ',cpt_actif
  print*,' n_ff:',n_ff
  print*,' n_fs:',n_fs
  print*,' n_ss ',n_ss
  print*,' n_fv ',n_fv

END SUBROUTINE read_Vloc_dof_bodies



!================================================
! Finding the velocity of the toppling
!================================================
subroutine toppling_v

  implicit none

  integer                                       :: i
  real*8                                        :: v_fall, dist_toppling
  real*8, dimension(3)                          :: pos_c_ini, pos_c_end


  if (find_ini_vt .eqv. .false.) then
    do i=1, nb_ligneCONTACT_POLYR
      if (TAB_CONTACT_POLYR(i)%nature == 'PRPLx') cycle
      if (TAB_CONTACT_POLYR(i)%icdent == 50 .or. TAB_CONTACT_POLYR(i)%ianent == 50) then
        if (TAB_CONTACT_POLYR(i)%rn .gt. 0) then
          find_ini_vt = .true.
          time_ini = time
          pos_c_ini(:) = TAB_CONTACT_POLYR(i)%coor_ctc(:)
        end if
      end if
    end do
  end if

  if (find_end_vt .eqv. .false.) then
    do i=1, nb_ligneCONTACT_POLYR
      if (TAB_CONTACT_POLYR(i)%nature == 'PRPLx') cycle
      if (TAB_CONTACT_POLYR(i)%icdent == 149 .or. TAB_CONTACT_POLYR(i)%ianent == 149) then
        if (TAB_CONTACT_POLYR(i)%rn .gt. 0.) then
          find_end_vt = .true.
          time_end = time
          pos_c_end(:) = TAB_CONTACT_POLYR(i)%coor_ctc(:)

          dist_toppling = sqrt((pos_c_end(1) - pos_c_ini(1))**2 + (pos_c_end(3) - pos_c_ini(3))**2)
          v_fall = dist_toppling/(time_end - time_ini)

          print*, 'asdadasdadadsdasdasdasadasas'
          print*, v_fall
          print*, time_end
          print*, time_ini
          print*, dist_toppling

          write(110,*) '#The falling velocity ', compteur_clout
          write(110,'(E12.5)') v_fall
        end if
      end if
    end do
  end if

  print*, 'Searching toppling velocity ---> Ok!'

end subroutine toppling_v


!================================================
! Drawing in vtk
!================================================
subroutine draw

  implicit none

  character(:), allocatable                      ::  vtk_part, vtk_counter
  character(:), allocatable                      ::  vtk_n_points
  logical                                        ::  dir_vtk
  character(len=6)                               ::  vtk_c_temp,  v_n_vertices
  integer                                        ::  i, j, k
  integer                                        ::  n_l_vertices, n_l_faces, curr_l_faces
  integer                                        ::  n_l_fields, l_counter
  real(kind=8), dimension(3)                     ::  curr_l_vector
  integer                                        ::  l_cdt, l_ant
  ! Variables for the forces
  character(:), allocatable                      ::  vtk_forces
  character(:), allocatable                      ::  vtk_n_v_forces
  character(len=6)                               ::  v_n_v_forces
  real(kind=8)                                   ::  vtk_ave_force, l_force_scale, l_force
  real(kind=8)                                   ::  l_rn, l_rt, l_rs
  real(kind=8), dimension(3)                     ::  vtk_cd_center, vtk_an_center, v_l_normal, v_l_t, v_l_s
  integer                                        ::  f_counter

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
    n_l_vertices = n_l_vertices + TAB_POLY(i)%nb_vertex 
  end do
  
  ! Adding the vertices of the walls
  n_l_vertices = n_l_vertices + (n_walls*8)
  
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
      !print*, 'Particle: ', i
      do j=1, TAB_POLY(i)%nb_vertex
        k = k + 1
        !print*, k
        if(k .lt. 3) then
          if (TAB_POLY(i)%vertex(1,j) .lt. 0) then
            write(121,'(E13.6,A)', advance = 'no') TAB_POLY(i)%vertex(1,j), ' '
          else
            write(121,'(E13.6,A)', advance = 'no') TAB_POLY(i)%vertex(1,j), ' '
          end if

          if (TAB_POLY(i)%vertex(2,j) .lt. 0) then
            write(121,'(E13.6,A)', advance = 'no') TAB_POLY(i)%vertex(2,j), ' '
          else
            write(121,'(E13.6,A)', advance = 'no') TAB_POLY(i)%vertex(2,j), ' '
          end if

          if (TAB_POLY(i)%vertex(3,j) .lt. 0) then
            write(121,'(E13.6,A)', advance = 'no') TAB_POLY(i)%vertex(3,j), ' '
          else
            write(121,'(E13.6,A)', advance = 'no') TAB_POLY(i)%vertex(3,j), ' '
          end if
        else
          if (TAB_POLY(i)%vertex(1,j) .lt. 0) then
            write(121,'(E13.6,A)', advance = 'no') TAB_POLY(i)%vertex(1,j), ' '
          else
            write(121,'(E13.6,A)', advance = 'no') TAB_POLY(i)%vertex(1,j), ' '
          end if

          if (TAB_POLY(i)%vertex(2,j) .lt. 0) then
            write(121,'(E13.6,A)', advance = 'no') TAB_POLY(i)%vertex(2,j), ' '
          else
            write(121,'(E13.6,A)', advance = 'no') TAB_POLY(i)%vertex(2,j), ' '
          end if

          if (TAB_POLY(i)%vertex(3,j) .lt. 0) then
            write(121,'(E13.6,A)') TAB_POLY(i)%vertex(3,j), ' '
          else
            write(121,'(E13.6,A)') TAB_POLY(i)%vertex(3,j), ' '
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
            write(121,'(E13.6,A)', advance = 'no') TAB_PLAN(i-n_particles)%vertex(1,j), ' '
          else
            write(121,'(E13.6,A)', advance = 'no') TAB_PLAN(i-n_particles)%vertex(1,j), ' '
          end if

          if (TAB_PLAN(i-n_particles)%vertex(2,j) .lt. 0) then
            write(121,'(E13.6,A)', advance = 'no') TAB_PLAN(i-n_particles)%vertex(2,j), ' '
          else
            write(121,'(E13.6,A)', advance = 'no') TAB_PLAN(i-n_particles)%vertex(2,j), ' '
          end if

          if (TAB_PLAN(i-n_particles)%vertex(3,j) .lt. 0) then
            write(121,'(E13.6,A)', advance = 'no') TAB_PLAN(i-n_particles)%vertex(3,j), ' '
          else 
            write(121,'(E13.6,A)', advance = 'no') TAB_PLAN(i-n_particles)%vertex(3,j), ' '
          end if
        else
          if (TAB_PLAN(i-n_particles)%vertex(1,j) .lt. 0) then
            write(121,'(E13.6,A)', advance = 'no') TAB_PLAN(i-n_particles)%vertex(1,j), ' '
          else 
            write(121,'(E13.6,A)', advance = 'no') TAB_PLAN(i-n_particles)%vertex(1,j), ' '
          end if
          
          if (TAB_PLAN(i-n_particles)%vertex(2,j) .lt. 0) then
            write(121,'(E13.6,A)', advance = 'no') TAB_PLAN(i-n_particles)%vertex(2,j), ' '
          else 
            write(121,'(E13.6,A)', advance = 'no') TAB_PLAN(i-n_particles)%vertex(2,j), ' '
          end if
          
          if (TAB_PLAN(i-n_particles)%vertex(3,j) .lt. 0) then
            write(121,'(E13.6,A)') TAB_PLAN(i-n_particles)%vertex(3,j), ' '
          else 
            write(121,'(E13.6,A)') TAB_PLAN(i-n_particles)%vertex(3,j), ' '
          end if
          
          k = 0
          
        end if
      end do
    end if
  end do
  
  write(121, '(A)') ' '
  
  ! Finding the total number of faces
  n_l_faces = 0
  
  do i=1, n_particles
    n_l_faces = n_l_faces + TAB_POLY(i)%n_faces
  end do
  
  ! Writing the conectivity between vertices
  ! Plus 6 faces each wall (=6*6)
  write(121,'(A)', advance='no') 'POLYGONS '
  write(121, '(2(I6,A))') n_l_faces + n_walls*6, ' ' , (n_l_faces*4)+(n_walls*6*5)
  
  curr_l_faces = 0
  
  do i=1, n_particles
    !write(121,*) 'Particle', i
    do j=1, TAB_POLY(i)%n_faces
      ! We have always triangles
      write(121, '(I1,A)', advance= 'no') 3, ' '
      
      ! First face... writing precisely 
      if (TAB_POLY(i)%face(1,j) + curr_l_faces - 1.lt. 10) then
        write(121, '(I1)', advance = 'no') TAB_POLY(i)%face(1,j) + curr_l_faces -1
      else if (TAB_POLY(i)%face(1,j) + curr_l_faces -1 .lt. 100) then
        write(121, '(I2)', advance = 'no') TAB_POLY(i)%face(1,j) + curr_l_faces -1 
      else if (TAB_POLY(i)%face(1,j) + curr_l_faces -1 .lt. 1000) then
        write(121, '(I3)', advance = 'no') TAB_POLY(i)%face(1,j) + curr_l_faces -1 
      else if (TAB_POLY(i)%face(1,j) + curr_l_faces -1 .lt. 10000) then
        write(121, '(I4)', advance = 'no') TAB_POLY(i)%face(1,j) + curr_l_faces -1 
      else if (TAB_POLY(i)%face(1,j) + curr_l_faces -1 .lt. 100000) then
        write(121, '(I5)', advance = 'no') TAB_POLY(i)%face(1,j) + curr_l_faces -1
      end if
      
      write(121, '(A)', advance='no') ' '
      
      ! Second face
      if (TAB_POLY(i)%face(2,j) + curr_l_faces -1 .lt. 10) then
        write(121, '(I1)', advance = 'no') TAB_POLY(i)%face(2,j) + curr_l_faces -1
      else if (TAB_POLY(i)%face(2,j) + curr_l_faces -1 .lt. 100) then
        write(121, '(I2)', advance = 'no') TAB_POLY(i)%face(2,j) + curr_l_faces -1
      else if (TAB_POLY(i)%face(2,j) + curr_l_faces -1 .lt. 1000) then
        write(121, '(I3)', advance = 'no') TAB_POLY(i)%face(2,j) + curr_l_faces -1
      else if (TAB_POLY(i)%face(2,j) + curr_l_faces -1 .lt. 10000) then
        write(121, '(I4)', advance = 'no') TAB_POLY(i)%face(2,j) + curr_l_faces -1
      else if (TAB_POLY(i)%face(2,j) + curr_l_faces -1 .lt. 100000) then
        write(121, '(I5)', advance = 'no') TAB_POLY(i)%face(2,j) + curr_l_faces -1
      end if
      
      ! Third face
      write(121, '(A)', advance='no') ' '
      
      if (TAB_POLY(i)%face(3,j) + curr_l_faces -1 .lt. 10) then
        write(121, '(I1)') TAB_POLY(i)%face(3,j) + curr_l_faces -1
      else if (TAB_POLY(i)%face(3,j) + curr_l_faces -1 .lt. 100) then
        write(121, '(I2)') TAB_POLY(i)%face(3,j) + curr_l_faces -1
      else if (TAB_POLY(i)%face(3,j) + curr_l_faces -1 .lt. 1000) then
        write(121, '(I3)') TAB_POLY(i)%face(3,j) + curr_l_faces -1 
      else if (TAB_POLY(i)%face(3,j) + curr_l_faces -1 .lt. 10000) then
        write(121, '(I4)') TAB_POLY(i)%face(3,j) + curr_l_faces -1
      else if (TAB_POLY(i)%face(3,j) + curr_l_faces -1 .lt. 100000) then
        write(121, '(I5)') TAB_POLY(i)%face(3,j) + curr_l_faces -1
      end if
      
    end do
    curr_l_faces = curr_l_faces + TAB_POLY(i)%nb_vertex
  end do
  
  ! Writing the conectivity for the walls
  ! For all the walls 
  do i=1, n_walls
    ! Each one with 6 faces
        
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
    
    !!!!! 2nd (5-6-7-8)
    write(121, '(I1,A)', advance= 'no') 4, ' '
    
    if (curr_l_faces + 4 .lt. 10) then
      write(121, '(I1)', advance = 'no') curr_l_faces + 4
    else if (curr_l_faces + 4 .lt. 100) then
      write(121, '(I2)', advance = 'no') curr_l_faces + 4
    else if (curr_l_faces + 4 .lt. 1000) then
      write(121, '(I3)', advance = 'no') curr_l_faces + 4
    else if (curr_l_faces + 4 .lt. 10000) then
      write(121, '(I4)', advance = 'no') curr_l_faces + 4 
    else if (curr_l_faces + 4 .lt. 100000) then
      write(121, '(I5)', advance = 'no') curr_l_faces + 4 
    end if
    
    write(121, '(A)', advance='no') ' '
    
    if (curr_l_faces +5 .lt. 10) then
      write(121, '(I1)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 100) then
      write(121, '(I2)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 1000) then
      write(121, '(I3)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 10000) then
      write(121, '(I4)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 100000) then
      write(121, '(I5)', advance = 'no') curr_l_faces + 5 
    end if
    
    write(121, '(A)', advance='no') ' '
    
    if (curr_l_faces +6 .lt. 10) then
      write(121, '(I1)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 100) then
      write(121, '(I2)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 1000) then
      write(121, '(I3)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 10000) then
      write(121, '(I4)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 100000) then
      write(121, '(I5)', advance = 'no') curr_l_faces + 6
    end if
    
    write(121, '(A)', advance='no') ' '
    
    if (curr_l_faces +7 .lt. 10) then
      write(121, '(I1)') curr_l_faces +7
    else if (curr_l_faces +7 .lt. 100) then
      write(121, '(I2)') curr_l_faces +7 
    else if (curr_l_faces +7 .lt. 1000) then
      write(121, '(I3)') curr_l_faces +7 
    else if (curr_l_faces +7 .lt. 10000) then
      write(121, '(I4)') curr_l_faces +7
    else if (curr_l_faces +7 .lt. 100000) then
      write(121, '(I5)') curr_l_faces +7
    end if
    
    !!!!! 3nd (1-5-8-4)
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
    
    if (curr_l_faces +4 .lt. 10) then
      write(121, '(I1)', advance = 'no') curr_l_faces + 4
    else if (curr_l_faces +4 .lt. 100) then
      write(121, '(I2)', advance = 'no') curr_l_faces + 4
    else if (curr_l_faces +4 .lt. 1000) then
      write(121, '(I3)', advance = 'no') curr_l_faces + 4
    else if (curr_l_faces +4 .lt. 10000) then
      write(121, '(I4)', advance = 'no') curr_l_faces + 4
    else if (curr_l_faces +4 .lt. 100000) then
      write(121, '(I5)', advance = 'no') curr_l_faces + 4
    end if
    
    write(121, '(A)', advance='no') ' '
    
    if (curr_l_faces +7 .lt. 10) then
      write(121, '(I1)', advance = 'no') curr_l_faces + 7
    else if (curr_l_faces +7 .lt. 100) then
      write(121, '(I2)', advance = 'no') curr_l_faces + 7
    else if (curr_l_faces +7 .lt. 1000) then
      write(121, '(I3)', advance = 'no') curr_l_faces + 7
    else if (curr_l_faces +7 .lt. 10000) then
      write(121, '(I4)', advance = 'no') curr_l_faces + 7
    else if (curr_l_faces +7 .lt. 100000) then
      write(121, '(I5)', advance = 'no') curr_l_faces + 7
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
    
    
    !!!!! 4th (2-6-7-3)
    write(121, '(I1,A)', advance= 'no') 4, ' '
    
    if (curr_l_faces + 1 .lt. 10) then
      write(121, '(I1)', advance = 'no') curr_l_faces + 1
    else if (curr_l_faces + 1 .lt. 100) then
      write(121, '(I2)', advance = 'no') curr_l_faces + 1
    else if (curr_l_faces + 1 .lt. 1000) then
      write(121, '(I3)', advance = 'no') curr_l_faces + 1
    else if (curr_l_faces + 1 .lt. 10000) then
      write(121, '(I4)', advance = 'no') curr_l_faces + 1
    else if (curr_l_faces + 1 .lt. 100000) then
      write(121, '(I5)', advance = 'no') curr_l_faces + 1
    end if
    
    write(121, '(A)', advance='no') ' '
    
    if (curr_l_faces +5 .lt. 10) then
      write(121, '(I1)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 100) then
      write(121, '(I2)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 1000) then
      write(121, '(I3)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 10000) then
      write(121, '(I4)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 100000) then
      write(121, '(I5)', advance = 'no') curr_l_faces + 5 
    end if
    
    write(121, '(A)', advance='no') ' '
    
    if (curr_l_faces +6 .lt. 10) then
      write(121, '(I1)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 100) then
      write(121, '(I2)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 1000) then
      write(121, '(I3)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 10000) then
      write(121, '(I4)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 100000) then
      write(121, '(I5)', advance = 'no') curr_l_faces + 6
    end if
    
    write(121, '(A)', advance='no') ' '
    
    if (curr_l_faces +2 .lt. 10) then
      write(121, '(I1)') curr_l_faces +2
    else if (curr_l_faces +2 .lt. 100) then
      write(121, '(I2)') curr_l_faces +2
    else if (curr_l_faces +2 .lt. 1000) then
      write(121, '(I3)') curr_l_faces +2
    else if (curr_l_faces +2 .lt. 10000) then
      write(121, '(I4)') curr_l_faces +2
    else if (curr_l_faces +2 .lt. 100000) then
      write(121, '(I5)') curr_l_faces +2
    end if
    
    !!!!! 5th (1-2-6-5)
    write(121, '(I1,A)', advance= 'no') 4, ' '
    
    if (curr_l_faces .lt. 10) then
      write(121, '(I1)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 100) then
      write(121, '(I2)', advance = 'no') curr_l_faces
    else if (curr_l_faces  .lt. 1000) then
      write(121, '(I3)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 10000) then
      write(121, '(I4)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 100000) then
      write(121, '(I5)', advance = 'no') curr_l_faces
    end if
    
    write(121, '(A)', advance='no') ' '
    
    if (curr_l_faces +1 .lt. 10) then
      write(121, '(I1)', advance = 'no') curr_l_faces + 1
    else if (curr_l_faces +1 .lt. 100) then
      write(121, '(I2)', advance = 'no') curr_l_faces + 1
    else if (curr_l_faces +1 .lt. 1000) then
      write(121, '(I3)', advance = 'no') curr_l_faces + 1
    else if (curr_l_faces +1 .lt. 10000) then
      write(121, '(I4)', advance = 'no') curr_l_faces + 1
    else if (curr_l_faces +1 .lt. 100000) then
      write(121, '(I5)', advance = 'no') curr_l_faces + 1
    end if
    
    write(121, '(A)', advance='no') ' '
    
    if (curr_l_faces +5 .lt. 10) then
      write(121, '(I1)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 100) then
      write(121, '(I2)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 1000) then
      write(121, '(I3)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 10000) then
      write(121, '(I4)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 100000) then
      write(121, '(I5)', advance = 'no') curr_l_faces + 5
    end if
    
    write(121, '(A)', advance='no') ' '
    
    if (curr_l_faces +4 .lt. 10) then
      write(121, '(I1)') curr_l_faces +4
    else if (curr_l_faces +4 .lt. 100) then
      write(121, '(I2)') curr_l_faces +4
    else if (curr_l_faces +4 .lt. 1000) then
      write(121, '(I3)') curr_l_faces +4
    else if (curr_l_faces +4 .lt. 10000) then
      write(121, '(I4)') curr_l_faces +4
    else if (curr_l_faces +4 .lt. 100000) then
      write(121, '(I5)') curr_l_faces +4
    end if
    
    
    !!!!! 6th (4-3-7-8)
    write(121, '(I1,A)', advance= 'no') 4, ' '
    
    if (curr_l_faces + 3 .lt. 10) then
      write(121, '(I1)', advance = 'no') curr_l_faces + 3
    else if (curr_l_faces + 3 .lt. 100) then
      write(121, '(I2)', advance = 'no') curr_l_faces + 3
    else if (curr_l_faces + 3 .lt. 1000) then
      write(121, '(I3)', advance = 'no') curr_l_faces + 3
    else if (curr_l_faces + 3 .lt. 10000) then
      write(121, '(I4)', advance = 'no') curr_l_faces + 3
    else if (curr_l_faces + 3 .lt. 100000) then
      write(121, '(I5)', advance = 'no') curr_l_faces + 3
    end if
    
    write(121, '(A)', advance='no') ' '
    
    if (curr_l_faces +2 .lt. 10) then
      write(121, '(I1)', advance = 'no') curr_l_faces + 2
    else if (curr_l_faces +2 .lt. 100) then
      write(121, '(I2)', advance = 'no') curr_l_faces + 2
    else if (curr_l_faces +2 .lt. 1000) then
      write(121, '(I3)', advance = 'no') curr_l_faces + 2
    else if (curr_l_faces +2 .lt. 10000) then
      write(121, '(I4)', advance = 'no') curr_l_faces + 2
    else if (curr_l_faces +2 .lt. 100000) then
      write(121, '(I5)', advance = 'no') curr_l_faces + 2
    end if
    
    write(121, '(A)', advance='no') ' '
    
    if (curr_l_faces +6 .lt. 10) then
      write(121, '(I1)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 100) then
      write(121, '(I2)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 1000) then
      write(121, '(I3)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 10000) then
      write(121, '(I4)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 100000) then
      write(121, '(I5)', advance = 'no') curr_l_faces + 6
    end if
    
    write(121, '(A)', advance='no') ' '
    
    if (curr_l_faces +7 .lt. 10) then
      write(121, '(I1)') curr_l_faces +7
    else if (curr_l_faces +7 .lt. 100) then
      write(121, '(I2)') curr_l_faces +7 
    else if (curr_l_faces +7 .lt. 1000) then
      write(121, '(I3)') curr_l_faces +7 
    else if (curr_l_faces +7 .lt. 10000) then
      write(121, '(I4)') curr_l_faces +7
    else if (curr_l_faces +7 .lt. 100000) then
      write(121, '(I5)') curr_l_faces +7
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
  write(121,'(A)') ''
  ! The cells begin
  write(121,'(A)', advance='no') 'CELL_DATA'
  
  ! Writing the number of data by field. It corresponds to the same number of faces
  write(121, '(2(I6,A))') n_l_faces + n_walls*6
  write(121, '(A,I4)')  'FIELD FieldData', n_l_fields
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Id
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! Naming the field, the dimension of the data, and number of lines (as before FIELD), and data type
  write(121,'(A,I1,I6,A)') 'Id ', 1, (n_l_faces +n_walls*6), ' float'
  k = 0
  l_counter = 0
  do i=1, n_particles + n_walls
    l_counter = l_counter + 1
    if (i .le. n_particles) then
      do j=1, TAB_POLY(i)%n_faces
        k=k+1
        if(k .le. 9) then
          write(121, '(I6)', advance='no') l_counter
        else 
          write(121, '(I6)') l_counter
          k=0
        end if 
      end do
    else
      ! Six facer for each wall 
      do j=1, 6
        k=k+1
        if(k .le. 9) then
          write(121, '(I6)', advance='no') l_counter
        else 
          write(121, '(I6)') l_counter
          k=0
        end if 
      end do
    end if
  end do
  ! And jump
  write(121, '(A)') ' '
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Material
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(121,'(A,I1,I6,A)') 'Material ', 1, (n_l_faces + n_walls*6), ' float'
  k = 0
  l_counter = 1
  do i=1, n_particles + n_walls
    if (i .le. n_particles) then
      ! Material number 1
      do j=1, TAB_POLY(i)%n_faces
        k=k+1
        if(k .le. 9) then
          write(121, '(I6)', advance='no') l_counter
        else 
          write(121, '(I6)') l_counter
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
          write(121, '(I6)', advance='no') l_counter
        else 
          write(121, '(I6)') l_counter
          k=0
        end if 
      end do
    end if
  end do
  
  ! And jump
  write(121, '(A)') ' '
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Dispacement
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  curr_l_vector(:) = 0.D0
  
  write(121,'(A,I1,I6,A)') 'Disp ', 3, (n_l_faces + n_walls*6), ' float'
  k = 0
  do i=1, n_particles + n_walls
    if (i .le. n_particles) then
      curr_l_vector(:) = TAB_POLY(i)%center(:) - TAB_POLY(i)%center_ref
      do j=1, TAB_POLY(i)%n_faces
        k=k+3
        if(k .le. 9) then
          write(121, '(3(F12.9,A))', advance='no') curr_l_vector(1), ' ', curr_l_vector(2), ' ', curr_l_vector(3), ' '
        else 
          write(121, '(3(F12.9,A))') curr_l_vector(1), ' ', curr_l_vector(2), ' ', curr_l_vector(3)
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
          write(121, '(3(F12.9,A))', advance='no') curr_l_vector(1), ' ', curr_l_vector(2), ' ', curr_l_vector(3), ' '
        else 
          write(121, '(3(F12.9,A))') curr_l_vector(1), ' ', curr_l_vector(2), ' ', curr_l_vector(3)
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
  
  write(121,'(A,I1,I6,A)') 'Veloc ', 3, (n_l_faces + n_walls*6), ' float'
  k = 0
  do i=1, n_particles + n_walls
    if (i .le. n_particles) then
      curr_l_vector(1) = TAB_POLY(i)%Vx
      curr_l_vector(2) = TAB_POLY(i)%Vy
      curr_l_vector(3) = TAB_POLY(i)%Vz
      
      do j=1, TAB_POLY(i)%n_faces
        k=k+3
        if(k .le. 9) then
          write(121, '(3(F13.9,A))', advance='no') curr_l_vector(1), ' ', curr_l_vector(2), ' ', curr_l_vector(3), ' '
        else 
          write(121, '(3(F13.9,A))') curr_l_vector(1), ' ', curr_l_vector(2), ' ', curr_l_vector(3)
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
          write(121, '(3(F13.9,A))', advance='no') curr_l_vector(1), ' ', curr_l_vector(2), ' ', curr_l_vector(3), ' '
        else 
          write(121, '(3(F13.9,A))') curr_l_vector(1), ' ', curr_l_vector(2), ' ', curr_l_vector(3)
          k=0
        end if 
      end do
    end if
  end do
  
  ! And jump
  write(121, '(A)') ' '
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Coordination
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  l_counter = 0
  
  write(121,'(A,I1,I6,A)') 'Z ', 1, (n_l_faces + n_walls*6), ' float'
  k = 0
  do i=1, n_particles + n_walls
    l_counter = 0
    if (i .le. n_particles) then
      ! Counting the number of contacts
      do j=1, nb_ligneCONTACT_POLYR
        l_cdt = TAB_CONTACT_POLYR(j)%icdent
        l_ant = TAB_CONTACT_POLYR(j)%ianent
        
        if(TAB_POLY(l_cdt)%behav /= 'PLEXx' .or. TAB_POLY(l_ant)%behav /= 'PLEXx') cycle
        if(l_cdt == i .or. l_ant == i) then
          l_counter = l_counter + 1 
        end if
      end do
        
      do j=1, TAB_POLY(i)%n_faces
        k=k+1
        if(k .le. 9) then
          write(121, '(I4,A)', advance='no') l_counter, ' '
        else 
          write(121, '(I4,A)') l_counter
          k=0
        end if
      end do
    else
      do j=1, 6
        k=k+1
        if(k .le. 9) then
          write(121, '(I4,A)', advance='no') 0, ' '
        else 
          write(121, '(I4,A)') 0
          k=0
        end if 
      end do
    end if
  end do
  
  ! And jump
  write(121, '(A)') ' '

  
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
  ! Eight vertices for each parallelepiped
  f_counter = 0
  do i=1, nb_ligneCONTACT_POLYR
    if(TAB_CONTACT_POLYR(i)%nature /= 'PRPRx') cycle
    f_counter = f_counter + 1
  end do
  
  write(v_n_v_forces, '(I6)') f_counter * 8
  
  vtk_n_v_forces = 'POINTS' // v_n_v_forces // ' float'
  write(121,'(A)') vtk_n_v_forces
  
  ! Finding the average force
  vtk_ave_force = 0.D0
  
  do i = 1, nb_ligneCONTACT_POLYR
    if (TAB_CONTACT_POLYR(i)%nature /= 'PRPRx') cycle
    vtk_ave_force = vtk_ave_force + abs(TAB_CONTACT_POLYR(i)%rn)
  end do
  
  vtk_ave_force = vtk_ave_force/nb_ligneCONTACT_POLYR
  
  ! Force scale parameter 
  l_force_scale = 0.D0
  do i=1, n_particles
    l_force_scale = l_force_scale + TAB_POLY(i)%Rmax
  end do
  ! The scale is here!
  ! 1% of average radius
  l_force_scale = (l_force_scale/n_particles)*0.01
  
  k=0
  do i=1, nb_ligneCONTACT_POLYR
    
    if (TAB_CONTACT_POLYR(i)%nature /= 'PRPRx') cycle
    ! The particles ids
    l_cdt = TAB_CONTACT_POLYR(i)%icdent
    l_ant = TAB_CONTACT_POLYR(i)%ianent
    ! The center of each particle
    vtk_cd_center(:) = TAB_POLY(l_cdt)%center(:)
    vtk_an_center(:) = TAB_POLY(l_ant)%center(:)
    ! The normal and tangential vectors
    v_l_normal(:) = TAB_CONTACT_POLYR(i)%n(:)
    v_l_t(:) = TAB_CONTACT_POLYR(i)%t(:)
    v_l_s(:) = TAB_CONTACT_POLYR(i)%s(:)
    ! The forces
    l_rn = TAB_CONTACT_POLYR(i)%rn
    l_rt = TAB_CONTACT_POLYR(i)%rt
    l_rs = TAB_CONTACT_POLYR(i)%rs
    
    l_force = (l_rn/vtk_ave_force)*l_force_scale
    
    ! Write! ... 4 vertices for each line
    
    ! Candidate --- First vertex. 
    write(121, '(3(F15.7,A))', advance='no') vtk_cd_center(1)+l_force*v_l_t(1)+ l_force*v_l_s(1), ' ', &
                                             vtk_cd_center(2)+l_force*v_l_t(2)+ l_force*v_l_s(2), ' ', &
                                             vtk_cd_center(3)+l_force*v_l_t(3)+ l_force*v_l_s(3), ' '
    ! Candidate --- Second vertex. 
    write(121, '(3(F15.7,A))', advance='no') vtk_cd_center(1)+l_force*v_l_t(1)- l_force*v_l_s(1), ' ', &
                                             vtk_cd_center(2)+l_force*v_l_t(2)- l_force*v_l_s(2), ' ', &
                                             vtk_cd_center(3)+l_force*v_l_t(3)- l_force*v_l_s(3), ' '
    ! Candidate --- Third vertex. 
    write(121, '(3(F15.7,A))', advance='no') vtk_cd_center(1)-l_force*v_l_t(1)- l_force*v_l_s(1), ' ', &
                                             vtk_cd_center(2)-l_force*v_l_t(2)- l_force*v_l_s(2), ' ', &
                                             vtk_cd_center(3)-l_force*v_l_t(3)- l_force*v_l_s(3), ' '
    ! Candidate --- 4th vertex. 
    write(121, '(3(F15.7,A))') vtk_cd_center(1)-l_force*v_l_t(1)+ l_force*v_l_s(1), ' ', &
                               vtk_cd_center(2)-l_force*v_l_t(2)+ l_force*v_l_s(2), ' ', &
                               vtk_cd_center(3)-l_force*v_l_t(3)+ l_force*v_l_s(3), ' ' 
    
    ! Antagonist --- First vertex. 
    write(121, '(3(F15.7,A))', advance='no') vtk_an_center(1)+l_force*v_l_t(1)+ l_force*v_l_s(1), ' ', &
                                             vtk_an_center(2)+l_force*v_l_t(2)+ l_force*v_l_s(2), ' ', &
                                             vtk_an_center(3)+l_force*v_l_t(3)+ l_force*v_l_s(3), ' '
    ! Antagonist --- Second vertex.          
    write(121, '(3(F15.7,A))', advance='no') vtk_an_center(1)+l_force*v_l_t(1)- l_force*v_l_s(1), ' ', &
                                             vtk_an_center(2)+l_force*v_l_t(2)- l_force*v_l_s(2), ' ', &
                                             vtk_an_center(3)+l_force*v_l_t(3)- l_force*v_l_s(3), ' '
    ! Antagonist --- Third vertex.           
    write(121, '(3(F15.7,A))', advance='no') vtk_an_center(1)-l_force*v_l_t(1)- l_force*v_l_s(1), ' ', &
                                             vtk_an_center(2)-l_force*v_l_t(2)- l_force*v_l_s(2), ' ', &
                                             vtk_an_center(3)-l_force*v_l_t(3)- l_force*v_l_s(3), ' '
    ! Antagonist --- 4th vertex.             
    write(121, '(3(F15.7,A))') vtk_an_center(1)-l_force*v_l_t(1)+ l_force*v_l_s(1), ' ', &
                               vtk_an_center(2)-l_force*v_l_t(2)+ l_force*v_l_s(2), ' ', &
                               vtk_an_center(3)-l_force*v_l_t(3)+ l_force*v_l_s(3), ' ' 
    
  end do
  
  ! Writing the conectivity
  write(121,'(A)', advance='no') 'POLYGONS '
  ! 6 faces for each parallepiped
  write(121, '(2(I6,A))') f_counter*6, ' ' , (f_counter*6*5)
  
  curr_l_faces = 0
  
  do i=1, nb_ligneCONTACT_POLYR
    ! Each one with 6 faces
    if (TAB_CONTACT_POLYR(i)%nature /= 'PRPRx') cycle
    
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
    
    !!!!! 2nd (5-6-7-8)
    write(121, '(I1,A)', advance= 'no') 4, ' '
    
    if (curr_l_faces + 4 .lt. 10) then
      write(121, '(I1)', advance = 'no') curr_l_faces + 4
    else if (curr_l_faces + 4 .lt. 100) then
      write(121, '(I2)', advance = 'no') curr_l_faces + 4
    else if (curr_l_faces + 4 .lt. 1000) then
      write(121, '(I3)', advance = 'no') curr_l_faces + 4
    else if (curr_l_faces + 4 .lt. 10000) then
      write(121, '(I4)', advance = 'no') curr_l_faces + 4 
    else if (curr_l_faces + 4 .lt. 100000) then
      write(121, '(I5)', advance = 'no') curr_l_faces + 4 
    end if
    
    write(121, '(A)', advance='no') ' '
    
    if (curr_l_faces +5 .lt. 10) then
      write(121, '(I1)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 100) then
      write(121, '(I2)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 1000) then
      write(121, '(I3)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 10000) then
      write(121, '(I4)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 100000) then
      write(121, '(I5)', advance = 'no') curr_l_faces + 5 
    end if
    
    write(121, '(A)', advance='no') ' '
    
    if (curr_l_faces +6 .lt. 10) then
      write(121, '(I1)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 100) then
      write(121, '(I2)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 1000) then
      write(121, '(I3)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 10000) then
      write(121, '(I4)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 100000) then
      write(121, '(I5)', advance = 'no') curr_l_faces + 6
    end if
    
    write(121, '(A)', advance='no') ' '
    
    if (curr_l_faces +7 .lt. 10) then
      write(121, '(I1)') curr_l_faces +7
    else if (curr_l_faces +7 .lt. 100) then
      write(121, '(I2)') curr_l_faces +7 
    else if (curr_l_faces +7 .lt. 1000) then
      write(121, '(I3)') curr_l_faces +7 
    else if (curr_l_faces +7 .lt. 10000) then
      write(121, '(I4)') curr_l_faces +7
    else if (curr_l_faces +7 .lt. 100000) then
      write(121, '(I5)') curr_l_faces +7
    end if
    
    !!!!! 3nd (1-5-8-4)
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
    
    if (curr_l_faces +4 .lt. 10) then
      write(121, '(I1)', advance = 'no') curr_l_faces + 4
    else if (curr_l_faces +4 .lt. 100) then
      write(121, '(I2)', advance = 'no') curr_l_faces + 4
    else if (curr_l_faces +4 .lt. 1000) then
      write(121, '(I3)', advance = 'no') curr_l_faces + 4
    else if (curr_l_faces +4 .lt. 10000) then
      write(121, '(I4)', advance = 'no') curr_l_faces + 4
    else if (curr_l_faces +4 .lt. 100000) then
      write(121, '(I5)', advance = 'no') curr_l_faces + 4
    end if
    
    write(121, '(A)', advance='no') ' '
    
    if (curr_l_faces +7 .lt. 10) then
      write(121, '(I1)', advance = 'no') curr_l_faces + 7
    else if (curr_l_faces +7 .lt. 100) then
      write(121, '(I2)', advance = 'no') curr_l_faces + 7
    else if (curr_l_faces +7 .lt. 1000) then
      write(121, '(I3)', advance = 'no') curr_l_faces + 7
    else if (curr_l_faces +7 .lt. 10000) then
      write(121, '(I4)', advance = 'no') curr_l_faces + 7
    else if (curr_l_faces +7 .lt. 100000) then
      write(121, '(I5)', advance = 'no') curr_l_faces + 7
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
    
    
    !!!!! 4th (2-6-7-3)
    write(121, '(I1,A)', advance= 'no') 4, ' '
    
    if (curr_l_faces + 1 .lt. 10) then
      write(121, '(I1)', advance = 'no') curr_l_faces + 1
    else if (curr_l_faces + 1 .lt. 100) then
      write(121, '(I2)', advance = 'no') curr_l_faces + 1
    else if (curr_l_faces + 1 .lt. 1000) then
      write(121, '(I3)', advance = 'no') curr_l_faces + 1
    else if (curr_l_faces + 1 .lt. 10000) then
      write(121, '(I4)', advance = 'no') curr_l_faces + 1
    else if (curr_l_faces + 1 .lt. 100000) then
      write(121, '(I5)', advance = 'no') curr_l_faces + 1
    end if
    
    write(121, '(A)', advance='no') ' '
    
    if (curr_l_faces +5 .lt. 10) then
      write(121, '(I1)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 100) then
      write(121, '(I2)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 1000) then
      write(121, '(I3)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 10000) then
      write(121, '(I4)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 100000) then
      write(121, '(I5)', advance = 'no') curr_l_faces + 5 
    end if
    
    write(121, '(A)', advance='no') ' '
    
    if (curr_l_faces +6 .lt. 10) then
      write(121, '(I1)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 100) then
      write(121, '(I2)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 1000) then
      write(121, '(I3)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 10000) then
      write(121, '(I4)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 100000) then
      write(121, '(I5)', advance = 'no') curr_l_faces + 6
    end if
    
    write(121, '(A)', advance='no') ' '
    
    if (curr_l_faces +2 .lt. 10) then
      write(121, '(I1)') curr_l_faces +2
    else if (curr_l_faces +2 .lt. 100) then
      write(121, '(I2)') curr_l_faces +2
    else if (curr_l_faces +2 .lt. 1000) then
      write(121, '(I3)') curr_l_faces +2
    else if (curr_l_faces +2 .lt. 10000) then
      write(121, '(I4)') curr_l_faces +2
    else if (curr_l_faces +2 .lt. 100000) then
      write(121, '(I5)') curr_l_faces +2
    end if
    
    !!!!! 5th (1-2-6-5)
    write(121, '(I1,A)', advance= 'no') 4, ' '
    
    if (curr_l_faces .lt. 10) then
      write(121, '(I1)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 100) then
      write(121, '(I2)', advance = 'no') curr_l_faces
    else if (curr_l_faces  .lt. 1000) then
      write(121, '(I3)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 10000) then
      write(121, '(I4)', advance = 'no') curr_l_faces
    else if (curr_l_faces .lt. 100000) then
      write(121, '(I5)', advance = 'no') curr_l_faces
    end if
    
    write(121, '(A)', advance='no') ' '
    
    if (curr_l_faces +1 .lt. 10) then
      write(121, '(I1)', advance = 'no') curr_l_faces + 1
    else if (curr_l_faces +1 .lt. 100) then
      write(121, '(I2)', advance = 'no') curr_l_faces + 1
    else if (curr_l_faces +1 .lt. 1000) then
      write(121, '(I3)', advance = 'no') curr_l_faces + 1
    else if (curr_l_faces +1 .lt. 10000) then
      write(121, '(I4)', advance = 'no') curr_l_faces + 1
    else if (curr_l_faces +1 .lt. 100000) then
      write(121, '(I5)', advance = 'no') curr_l_faces + 1
    end if
    
    write(121, '(A)', advance='no') ' '
    
    if (curr_l_faces +5 .lt. 10) then
      write(121, '(I1)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 100) then
      write(121, '(I2)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 1000) then
      write(121, '(I3)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 10000) then
      write(121, '(I4)', advance = 'no') curr_l_faces + 5
    else if (curr_l_faces +5 .lt. 100000) then
      write(121, '(I5)', advance = 'no') curr_l_faces + 5
    end if
    
    write(121, '(A)', advance='no') ' '
    
    if (curr_l_faces +4 .lt. 10) then
      write(121, '(I1)') curr_l_faces +4
    else if (curr_l_faces +4 .lt. 100) then
      write(121, '(I2)') curr_l_faces +4
    else if (curr_l_faces +4 .lt. 1000) then
      write(121, '(I3)') curr_l_faces +4
    else if (curr_l_faces +4 .lt. 10000) then
      write(121, '(I4)') curr_l_faces +4
    else if (curr_l_faces +4 .lt. 100000) then
      write(121, '(I5)') curr_l_faces +4
    end if
    
    
    !!!!! 6th (4-3-7-8)
    write(121, '(I1,A)', advance= 'no') 4, ' '
    
    if (curr_l_faces + 3 .lt. 10) then
      write(121, '(I1)', advance = 'no') curr_l_faces + 3
    else if (curr_l_faces + 3 .lt. 100) then
      write(121, '(I2)', advance = 'no') curr_l_faces + 3
    else if (curr_l_faces + 3 .lt. 1000) then
      write(121, '(I3)', advance = 'no') curr_l_faces + 3
    else if (curr_l_faces + 3 .lt. 10000) then
      write(121, '(I4)', advance = 'no') curr_l_faces + 3
    else if (curr_l_faces + 3 .lt. 100000) then
      write(121, '(I5)', advance = 'no') curr_l_faces + 3
    end if
    
    write(121, '(A)', advance='no') ' '
    
    if (curr_l_faces +2 .lt. 10) then
      write(121, '(I1)', advance = 'no') curr_l_faces + 2
    else if (curr_l_faces +2 .lt. 100) then
      write(121, '(I2)', advance = 'no') curr_l_faces + 2
    else if (curr_l_faces +2 .lt. 1000) then
      write(121, '(I3)', advance = 'no') curr_l_faces + 2
    else if (curr_l_faces +2 .lt. 10000) then
      write(121, '(I4)', advance = 'no') curr_l_faces + 2
    else if (curr_l_faces +2 .lt. 100000) then
      write(121, '(I5)', advance = 'no') curr_l_faces + 2
    end if
    
    write(121, '(A)', advance='no') ' '
    
    if (curr_l_faces +6 .lt. 10) then
      write(121, '(I1)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 100) then
      write(121, '(I2)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 1000) then
      write(121, '(I3)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 10000) then
      write(121, '(I4)', advance = 'no') curr_l_faces + 6
    else if (curr_l_faces +6 .lt. 100000) then
      write(121, '(I5)', advance = 'no') curr_l_faces + 6
    end if
    
    write(121, '(A)', advance='no') ' '
    
    if (curr_l_faces +7 .lt. 10) then
      write(121, '(I1)') curr_l_faces +7
    else if (curr_l_faces +7 .lt. 100) then
      write(121, '(I2)') curr_l_faces +7 
    else if (curr_l_faces +7 .lt. 1000) then
      write(121, '(I3)') curr_l_faces +7 
    else if (curr_l_faces +7 .lt. 10000) then
      write(121, '(I4)') curr_l_faces +7
    else if (curr_l_faces +7 .lt. 100000) then
      write(121, '(I5)') curr_l_faces +7
    end if
    
    curr_l_faces = curr_l_faces + 8
    
  end do
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!! FORCES FIELDS !!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! How many fields there will be?
  n_l_fields = 4                          ! They are: Type force (Compression - tension), RN, RT
                                          ! ...
  
  ! A blank space
  write(121,'(A)') ''
  ! The cells begin
  write(121,'(A)', advance='no') 'CELL_DATA'
  
  ! Writing the number of data by field. It corresponds to the same number of faces
  write(121, '(2(I6,A))') f_counter*6
  write(121, '(A,I4)')  'FIELD FieldData', n_l_fields
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Type force (1 Compression, -1 Tension)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! Naming the field, the dimension of the data, and number of lines (as before FIELD), and data type
  write(121,'(A,I1,I6,A)') 'Type ', 1, f_counter*6, ' float'
  k = 0
  do i=1, nb_ligneCONTACT_POLYR
    if (TAB_CONTACT_POLYR(i)%nature /= 'PRPRx') cycle
    do j=1, 6
      k=k+1
      if(k .le. 9) then
        if (TAB_CONTACT_POLYR(i)%rn .lt. 0) then
          write(121, '(I3)', advance='no') -1
        else 
          write(121, '(I3)', advance='no') 1
        end if
      else 
        if (TAB_CONTACT_POLYR(i)%rn .lt. 0) then
          write(121, '(I3)') -1
        else 
          write(121, '(I3)') 1
        end if
        k=0
      end if 
    end do
  end do
  ! And jump
  write(121, '(A)') ' '  
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Normal Force
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! Naming the field, the dimension of the data, and number of lines (as before FIELD), and data type
  write(121,'(A,I1,I6,A)') 'RN ', 1, f_counter*6, ' float'
  k = 0
  do i=1, nb_ligneCONTACT_POLYR
    if (TAB_CONTACT_POLYR(i)%nature /= 'PRPRx') cycle
    do j=1, 6
      k=k+1
      if(k .le. 9) then
        write(121, '(F15.7)', advance='no') TAB_CONTACT_POLYR(i)%rn
      else 
        write(121, '(F15.7)') TAB_CONTACT_POLYR(i)%rn
        k=0
      end if 
    end do
  end do
  ! And jump
  write(121, '(A)') ' '
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Tangential force (Rs+Rt)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! Naming the field, the dimension of the data, and number of lines (as before FIELD), and data type
  write(121,'(A,I1,I6,A)') 'RT ', 1, f_counter*6, ' float'
  k = 0
  do i=1, nb_ligneCONTACT_POLYR
    if (TAB_CONTACT_POLYR(i)%nature /= 'PRPRx') cycle
    do j=1, 6
      k=k+1
      if(k .le. 9) then
        write(121, '(F15.7)', advance='no') sqrt(TAB_CONTACT_POLYR(i)%rt**2 + TAB_CONTACT_POLYR(i)%rs**2)
      else 
        write(121, '(F15.7)') sqrt(TAB_CONTACT_POLYR(i)%rt**2 + TAB_CONTACT_POLYR(i)%rs**2)
        k=0
      end if 
    end do
  end do
  ! And jump
  write(121, '(A)') ' '
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Contact STATUS (0: noctc, 1=stick, 2=slide, 3=Mstck)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! Naming the field, the dimension of the data, and number of lines (as before FIELD), and data type
  write(121,'(A,I1,I6,A)') 'Status ', 1, f_counter*6, ' float'
  k = 0
  do i=1, nb_ligneCONTACT_POLYR
    if (TAB_CONTACT_POLYR(i)%nature /= 'PRPRx') cycle
    do j=1, 6
      k=k+1
      if(k .le. 9) then
        write(121, '(F15.7)', advance='no') TAB_CONTACT_POLYR(i)%x_status
      else
        write(121, '(F15.7)') TAB_CONTACT_POLYR(i)%x_status
        k=0
      end if
    end do
  end do
  ! And jump
  write(121, '(A)') ' '
  
  close(121)
  
  print*, 'Drawing in vtk       ---> Ok!'

end subroutine draw


!================================================
! Closing files
!================================================
subroutine close_all

  implicit none

  ! Variables
  integer                                        :: i

  print*,'--> Closing ports <--'

  ! Deallocating the polyhedra container
  do,i=1,n_particles
    deallocate(TAB_POLY(i)%vertex_ref)
    deallocate(TAB_POLY(i)%vertex)
    deallocate(TAB_POLY(i)%face)
  end do

  if (c_toppling_v                              ==1) then
    if (find_end_vt .eqv. .false.) then
      write(110,*) '# I did not find a velocity'
      write(110,*) '-1.0'
    end if
    close(110)
  end if

  deallocate(TAB_POLY)

end subroutine close_all


end program post3D_polyhe