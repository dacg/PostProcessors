!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!  POST-PROCESSOR !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! POLYHEDRA !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Notes: - This code is build to read up to 1000 OUTBOX frames

program post3D_polyhe

implicit none

type  ::  T_CORPS
  integer                                        ::  n_vertices,n_faces,id,gmv_material,nb_ctc
  integer                                        ::  n_contactors, n_faces_cluster, n_vertices_cluster
  real(kind=8),dimension(:,:),pointer            ::  vertex,vertex_ref
  integer,dimension(:,:),pointer                 ::  face
  real(kind=8)                                   ::  Rmax,volume,I1,I2,I3,&
                                                     Sphericity,Aire
  real(kind=8),dimension(3)                      ::  center,center_ref
  real(kind=8)                                   ::  Vx,Vy,Vz,Vrx,Vry,Vrz
  real(kind=8)                                   ::  Wx,Wy,Wz
  real(kind=8)                                   ::  ax1, ax2, ax3
  real(kind=8)                                   ::  Norme_V,Norme_Vr,Force
  real(kind=8)                                   ::  av_fn, av_mob
  real(kind=8),dimension(3,3)                    ::  momentI,Rot
  character(len=5)                               ::  behav,color
  real(kind=8)                                   ::  nctc,nctcslide,nctcstick,nctcs,nctcd,&
                                                     nctct,Pressure
end type T_CORPS

type  ::  T_FACE
  integer,dimension(2)                           ::  face_id
  real(kind=8)                                   ::  face_status
  real(kind=8), dimension(3)                     ::  nface_vector
end type T_FACE

type  ::  T_CONTACT
  integer                                        ::  icdent,ianent,type, cdver, id
  real(kind=8),dimension(3)                      ::  n,t,s
  real(kind=8),dimension(3)                      ::  coor_ctc
  real(kind=8)                                   ::  rn,rt,rs,vls,vln,vlt, gapTT
  integer                                        ::  sect
  character(len=5)                               ::  nature,status, i_law
  logical                                        ::  deja_compte,cycle_gmv
  real(kind=8)                                   ::  x_status
end type T_CONTACT

type  ::  T_GEO
  integer,dimension(:),pointer                   ::  g_face
  integer                                        ::  g_n_faces
end type T_GEO

type(T_CORPS),dimension(:),pointer               ::  TAB_POLY,TAB_PLAN
type(T_CONTACT),dimension(:),allocatable         ::  TAB_CONTACT,TAB_CONTACT_POLYR

!DEFINITIONS Vloc-Dof-Bodies
!================================================
integer                                          ::  n_walls=0,n_particles=0
integer                                          ::  compteur_clout=0
logical                                          ::  fin_post=.false.
integer                                          ::  step,nb_ligneCONTACT,&
                                                     nb_ligneCONTACT_POLYR, nf_PRPRx,nf_PRPLx
real(kind=8)                                     ::  time, x_period, y_period, mean_n_force

!DEFINITIONS BOITE
!================================================
real(kind=8)                                     ::  H_ini,Long_ini,Larg_ini,&
                                                     Depsilon_p,Depsilon_q
logical                                          ::  deformation_ini = .true.

!CALCUL DU PROFIL DE VITESSE MOYEN
!================================================
type  ::  T_step_velocity
  real(kind=8),dimension(:,:),allocatable        ::  tab_vmoyen
  real(kind=8),dimension(:,:),allocatable        ::  tab_vrmoyen
end type T_step_velocity

!CALCUL DU PROFIL DE CONTRAINTE MOYEN
!================================================
type  ::  T_sigma
  real(kind=8),dimension(3,3)                    :: sigma
end type T_sigma
type  ::  T_step_contrainte
  type(T_sigma),dimension(:),allocatable         ::  tab_sigma
end type T_step_contrainte

!COMMANDES D APPELLE
!================================================
character(len=30)                                ::  command
real(kind=8)                                     ::  PI = 3.14159265359,&
                                                     width,large,height
integer                                          ::  c_avg_velocity=0,&
                                                     c_coordination=0, c_qoverp=0, c_compacity=0, &
                                                     c_walls_force=0, c_walls_position=0, &
                                                     c_anisotropy_contact=0, c_anisotropy_force=0, &
                                                     c_anisotropy_branch=0, c_draw = 0, &
                                                     c_float_particles=0, c_coordination_class=0, &
                                                     c_list_n_forces=0, c_list_branches=0, c_list_mobilization=0, &
                                                     c_fn_class=0, c_mobilization_class=0, &
                                                     c_ctc_dir=0, c_brc_dir=0, c_frc_dir=0, c_gap_evol=0

! Variables David C.
logical                                          ::  copy_input=.false.
logical                                          ::  first_over_all = .true.
integer                                          ::  qoverp_option
real(kind=8)                                     ::  total_p_volume=0

!================================================
!Lecture du fichier de commande
!================================================
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

  if (command=='COPY IMPUT                   :') then
    copy_input = .true.
    cycle
  end if

  if (command=='VLOC_DOF_BODIES              :') then
    read(1,*) n_walls
    read(1,*) n_particles
    read(1,*) x_period
    read(1,*) y_period
    cycle
  end if

  if (command=='COORDINATION                 :') then
   ! Uses port 100
   c_coordination=1
   open (unit=100,file='./POSTPRO/COORDINATION.DAT',status='replace')
   cycle
  end if

  if (command=='AVERAGE VELOCITY             :') then
    ! Uses port 101
    c_avg_velocity=1
    open (unit=101,file='VITESSE_MOYENNE.DAT',status='replace')
    cycle
  end if

  if (command=='QoverP                       :') then
    ! Uses port 102 and may use 152
    c_qoverp=1
    read(1,*) qoverp_option
    open (unit=102,file='./POSTPRO/QoverP.DAT',status='replace')
    if (qoverp_option == 1) then
      open (unit=152,file='./POSTPRO/STRESSTENSOR.DAT',status='replace')
    end if
    cycle
  end if

  if (command=='COMPACITY                    :') then
    ! Uses port 103
    c_compacity=1
    open (unit=103,file='./POSTPRO/COMPACITY.DAT',status='replace')
    cycle
  end if

  if (command=='ANISOTROPY CONTACT           :') then
    ! Uses port 104
    c_anisotropy_contact=1
    open (unit=104,file='./POSTPRO/ANISOTROPY_CONTACT.DAT',status='replace')
    cycle
  end if

  if (command=='ANISOTROPY FORCE             :') then
    ! Uses port 105
    c_anisotropy_force=1
    open (unit=105,file='./POSTPRO/ANISOTROPY_FORCE.DAT',status='replace')
    cycle
  end if

  if (command=='ANISOTROPY BRANCH            :') then
    ! Uses port 106
    c_anisotropy_branch=1
    open (unit=106,file='./POSTPRO/ANISOTROPY_BRANCH.DAT',status='replace')
    cycle
  end if

  if (command=='WALLS POSITION               :') then
    ! Uses port 108
    c_walls_position=1
    open (unit=108,file='./POSTPRO/WALLSPOS.DAT',status='replace')
    cycle
  end if

  if (command=='WALLS FORCES                 :') then
    ! Uses port 109
    c_walls_force =1
    open (unit=109,file='./POSTPRO/WALLSFORCE.DAT',status='replace')
    cycle
  end if

  !if (command=='CONTACT NUMBER PROBABILITY   :') then
  !  ! Uses port 111
  !  c_n_ctc_probability = 1
  !  open (unit=111,file='./POSTPRO/CTC_PROBABILITY.DAT',status='replace')
  !  cycle
  !end if

  if (command=='FLOATING PARTICLES           :') then
    ! Uses port 112
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

  if (command=='GAP EVOLUTION                :') then
    c_gap_evol = 1
    ! Uses port 123
    open (unit=123,file='./POSTPRO/GAP_EVOLUTION.DAT',status='replace')
    cycle
  end if

  if (command=='END                          :') exit
end do

close(1)

!=================================================================
!Appel des different ressources lues dans le fichier de commande
!=================================================================
do

  if (copy_input .and. first_over_all) then
    call system('cp DATBOX/DOF.INI OUTBOX/DOF.OUT.0')
    call system('cp DATBOX/Vloc_Rloc.INI OUTBOX/Vloc_Rloc.OUT.0')
  endif

  call read_Vloc_dof_bodies

  call Calcul_deformation

  if (c_qoverp                    == 1)  call qoverp
  if (c_coordination              == 1)  call nb_coordination
  if (c_walls_position            == 1)  call walls_position
  if (c_compacity                 == 1)  call compacity
  if (c_walls_force               == 1)  call walls_force
  if (c_draw                      == 1)  call draw
  if (c_anisotropy_contact        == 1)  call anisotropy_contact
  if (c_anisotropy_force          == 1)  call anisotropy_force
  if (c_anisotropy_branch         == 1)  call anisotropy_branch
  if (c_float_particles           == 1)  call float_particles
  if (c_coordination_class        == 1)  call coordination_class
  if (c_list_n_forces             == 1)  call list_n_forces
  if (c_list_branches             == 1)  call list_branches
  if (c_list_mobilization         == 1)  call list_mobilization
  if (c_fn_class                  == 1)  call fn_class
  if (c_mobilization_class        == 1)  call mobilization_class
  if (c_ctc_dir                   == 1)  call ctc_dir
  if (c_brc_dir                   == 1)  call branch_dir
  if (c_frc_dir                   == 1)  call forces_dir
  if (c_gap_evol                  == 1)  call gap_evol

  if (fin_post) then
    call close_all
    print*,'--> End of post-processing <--'
    exit
  endif

  if(first_over_all) first_over_all = .false.
end do

!================================================
!================================================
!================================================
!================================================

 contains

!======================================================================
!Lecture des fichiers Vloc, Dof et BODIES et on commence remplir notre type
!TAB_POLY.... c'est un peu bancale je l'accorde volontier :-)
!======================================================================
subroutine read_Vloc_dof_bodies

  implicit none

  integer                              ::  i, ii, j,n_prpr,n_prpl,num_part,err,cd,an,n_vertices,n_faces
  integer                              ::  n_PRPRx, n_PRPLx, icdent,ianent,nb_ctc_deja_compte,cpt,cpt_actif, &
                                           icdver
  real(kind=8)                         ::  rn,rs,rt,t1,t2,t3,n1,n2,n3,s1,s2,s3,vls,vln,vlt, igap
  real(kind=8)                         ::  I1,I2,I3,mean_Sphericity
  character(len=6)                     ::  text
  character(len=5)                     ::  status, law
  character(len=5)                     ::  color,behav
  character(len=13)                    ::  text2
  real(kind=8)                         ::  center1,center2,center3,ax1r,ax2r,ax3r,coor1,coor2,coor3
  real(kind=8)                         ::  coor_ctc1,coor_ctc2,coor_ctc3
  real(kind=8)                         ::  ver1,ver2,ver3,Volume_tetra,Air_triangle,demi_somme
  integer                              ::  face1,face2,face3
  real(kind=8),dimension(3)            ::  X,Y,Z,center_gravity,ab,ac,bc, center_gravity_first
  real(kind=8),dimension(3,3)          ::  rot
  real(kind=8)                         ::  Vrx,Vry,Vrz,vx,vy,vz,dcd_ptc_contact,dan_ptc_contact,&
                                           n_ff,n_fs,n_fv,n_ss
  character(len=21)                    ::  clout_DOF
  character(len=27)                    ::  clout_Vloc
  character(len=20)                    ::  clout_Bodies

  real(kind=8)                         ::  curr_rad
  integer                              ::  temp_i_var

  clout_Bodies = './OUTBOX/BODIES.OUT'
  clout_DOF  =  './OUTBOX/DOF.OUT.    '
  clout_Vloc =  './OUTBOX/Vloc_Rloc.OUT.    '

  ! Setting the index of the output files
  compteur_clout = compteur_clout+1

  if (copy_input .and. first_over_all) then
    compteur_clout = compteur_clout - 1
  endif

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
    mean_Sphericity = 0
    temp_i_var = 0

    allocate(TAB_POLY(n_particles))
    allocate(TAB_PLAN(n_walls))

    ! A first lecture in order to count the number of contactors per cluster
    i=0
    TAB_POLY(:)%n_contactors = 0
    do
      read(2,'(A6)') text
      if (text == '$bdyty') then
        i=i+1
        if (i<n_particles+1) then
          do
            read(2,'(A6)') text
            if (text == ' POLYR') then
              TAB_POLY(i)%n_contactors = TAB_POLY(i)%n_contactors + 1
            else if (text == '$$$$$$') then
              exit
            end if
          end do
        else
          exit
        end if
      end if
    end do

    ! Going all the way up to read the rest of information
    rewind(2)

    i=0
    do
      read(2,'(A6)') text
      if (text == '$bdyty') then
        i=i+1
        !print*, i
        ! Storing the particles information
        if (i<n_particles+1) then
          TAB_POLY(i)%id=i
          read(2,*)
          read(2,*)
          read(2,'(22X, A5,7X,D14.7)') behav, curr_rad

          ! Computing the average equivalent particle radius
          !avr_part_rad = avr_part_rad + curr_rad

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
          ! Reading all the contactors
          ! I assume that cluster may have different number of contactors,
          ! however every contactor have the same shape.
          read(2,'(22X,A5,12X,i7,13X,i7)') color,n_vertices,n_faces
          TAB_POLY(i)%color     = color
          TAB_POLY(i)%n_faces   = n_faces
          TAB_POLY(i)%n_vertices = n_vertices
          TAB_POLY(i)%n_faces_cluster = n_faces*TAB_POLY(i)%n_contactors
          TAB_POLY(i)%n_vertices_cluster = n_vertices*TAB_POLY(i)%n_contactors

          allocate(TAB_POLY(i)%vertex_ref(3,TAB_POLY(i)%n_vertices_cluster))
          allocate(TAB_POLY(i)%vertex(3,TAB_POLY(i)%n_vertices_cluster))
          allocate(TAB_POLY(i)%face(3,TAB_POLY(i)%n_faces_cluster))

          ! Reading all the vertices and faces of contactors
          ii=0
          do
            ii=ii+1
            !print*, ii
            if (ii .le. TAB_POLY(i)%n_contactors) then
              ! Counting the number of faces
              temp_i_var = temp_i_var + n_faces

              ! Reading the coordinates of the vertices
              do j=1,TAB_POLY(i)%n_vertices
                read(2,'(29X, 3(5x,D14.7,2X))') ver1,ver2,ver3
                TAB_POLY(i)%vertex_ref(1,j+(ii-1)*TAB_POLY(i)%n_vertices) = ver1
                TAB_POLY(i)%vertex_ref(2,j+(ii-1)*TAB_POLY(i)%n_vertices) = ver2
                TAB_POLY(i)%vertex_ref(3,j+(ii-1)*TAB_POLY(i)%n_vertices) = ver3
              end do
              do j=1,TAB_POLY(i)%n_faces
                read(2,'(29X, 3(5x,i7,9X))') face1,face2,face3
                TAB_POLY(i)%face(1,j+(ii-1)*TAB_POLY(i)%n_faces) = face1+(ii-1)*TAB_POLY(i)%n_vertices
                TAB_POLY(i)%face(2,j+(ii-1)*TAB_POLY(i)%n_faces) = face2+(ii-1)*TAB_POLY(i)%n_vertices
                TAB_POLY(i)%face(3,j+(ii-1)*TAB_POLY(i)%n_faces) = face3+(ii-1)*TAB_POLY(i)%n_vertices
              end do
              read(2,*)
            else
              exit
            end if
          end do

          !print*, TAB_POLY(i)%face

          ! Finding the centroid of the polyhedron-cluster
          center_gravity = 0
          do j=1,TAB_POLY(i)%n_vertices_cluster
           center_gravity(1)=center_gravity(1)+TAB_POLY(i)%vertex_ref(1,j)
           center_gravity(2)=center_gravity(2)+TAB_POLY(i)%vertex_ref(2,j)
           center_gravity(3)=center_gravity(3)+TAB_POLY(i)%vertex_ref(3,j)
          enddo

          center_gravity=center_gravity/real(TAB_POLY(i)%n_vertices_cluster,8)

          ! Finding the centroid of the first element in the polyhedron-cluster
          center_gravity_first = 0
          do j=1,TAB_POLY(i)%n_vertices
           center_gravity_first(1)=center_gravity_first(1)+TAB_POLY(i)%vertex_ref(1,j)
           center_gravity_first(2)=center_gravity_first(2)+TAB_POLY(i)%vertex_ref(2,j)
           center_gravity_first(3)=center_gravity_first(3)+TAB_POLY(i)%vertex_ref(3,j)
          enddo

          center_gravity_first=center_gravity_first/real(TAB_POLY(i)%n_vertices,8)

          ! Finding the farthest vertex from the centroid
          TAB_POLY(i)%Rmax = 0
          do j=1,TAB_POLY(i)%n_vertices_cluster
            X(1) = TAB_POLY(i)%vertex_ref(1,j) - center_gravity(1)
            X(2) = TAB_POLY(i)%vertex_ref(2,j) - center_gravity(2)
            X(3) = TAB_POLY(i)%vertex_ref(3,j) - center_gravity(3)

            TAB_POLY(i)%Rmax = max(TAB_POLY(i)%Rmax , sqrt(X(1)**2 + X(2)**2 + X(3)**2))
          end do

          ! Finding the volume of the particle-cluster
          Volume_tetra = 0.D0
          Air_triangle = 0.D0
          do j=1,TAB_POLY(i)%n_faces
            X(1) = TAB_POLY(i)%vertex_ref(1,face1) - center_gravity_first(1)
            X(2) = TAB_POLY(i)%vertex_ref(2,face1) - center_gravity_first(2)
            X(3) = TAB_POLY(i)%vertex_ref(3,face1) - center_gravity_first(3)

            Y(1) = TAB_POLY(i)%vertex_ref(1,face2) - center_gravity_first(1)
            Y(2) = TAB_POLY(i)%vertex_ref(2,face2) - center_gravity_first(2)
            Y(3) = TAB_POLY(i)%vertex_ref(3,face2) - center_gravity_first(3)

            Z(1) = TAB_POLY(i)%vertex_ref(1,face3) - center_gravity_first(1)
            Z(2) = TAB_POLY(i)%vertex_ref(2,face3) - center_gravity_first(2)
            Z(3) = TAB_POLY(i)%vertex_ref(3,face3) - center_gravity_first(3)

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

          ! Thinking it is a cluster with the same basic shape
          TAB_POLY(i)%volume = (Volume_tetra/6.0)*TAB_POLY(i)%n_contactors

          ! Computing the total volume in the box
          total_p_volume = total_p_volume + TAB_POLY(i)%volume

          TAB_POLY(i)%Aire       = Air_triangle
          TAB_POLY(i)%Sphericity = (PI**(0.333333333)) * ((6*TAB_POLY(i)%volume)**(0.66666666)) / &
                                        TAB_POLY(i)%Aire    !QUE PITOS ES ESTO?

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
          TAB_POLY(i)%Pressure = 0.D0

          mean_Sphericity = mean_Sphericity + TAB_POLY(i)%Sphericity

        !--Storing the walls information
        else
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
          TAB_PLAN(i-n_particles)%ax1=ax1r
          TAB_PLAN(i-n_particles)%ax2=ax2r
          TAB_PLAN(i-n_particles)%ax3=ax3r
          read(2,*)
          TAB_PLAN(i-n_particles)%center(:)=0.D0
          TAB_PLAN(i-n_particles)%behav = 'WALL_'

          ! Similar to poly, i'm going to create vertices (8 vertices in each wall)
          allocate(TAB_PLAN(i-n_particles)%vertex_ref(3,8))
          allocate(TAB_PLAN(i-n_particles)%vertex(3,8))

          ! First sheet
          TAB_PLAN(i-n_particles)%vertex_ref(1,1) =  ax1r
          TAB_PLAN(i-n_particles)%vertex_ref(2,1) =  ax2r
          TAB_PLAN(i-n_particles)%vertex_ref(3,1) =  ax3r

          TAB_PLAN(i-n_particles)%vertex_ref(1,2) = (- ax1r)
          TAB_PLAN(i-n_particles)%vertex_ref(2,2) =  ax2r
          TAB_PLAN(i-n_particles)%vertex_ref(3,2) =  ax3r

          TAB_PLAN(i-n_particles)%vertex_ref(1,3) = (- ax1r)
          TAB_PLAN(i-n_particles)%vertex_ref(2,3) = (- ax2r)
          TAB_PLAN(i-n_particles)%vertex_ref(3,3) =  ax3r

          TAB_PLAN(i-n_particles)%vertex_ref(1,4) =  ax1r
          TAB_PLAN(i-n_particles)%vertex_ref(2,4) = (- ax2r)
          TAB_PLAN(i-n_particles)%vertex_ref(3,4) =  ax3r

          ! Second sheet
          TAB_PLAN(i-n_particles)%vertex_ref(1,5) =  ax1r
          TAB_PLAN(i-n_particles)%vertex_ref(2,5) =  ax2r
          TAB_PLAN(i-n_particles)%vertex_ref(3,5) = (- ax3r)

          TAB_PLAN(i-n_particles)%vertex_ref(1,6) = (- ax1r)
          TAB_PLAN(i-n_particles)%vertex_ref(2,6) =  ax2r
          TAB_PLAN(i-n_particles)%vertex_ref(3,6) = (- ax3r)

          TAB_PLAN(i-n_particles)%vertex_ref(1,7) = (- ax1r)
          TAB_PLAN(i-n_particles)%vertex_ref(2,7) = (- ax2r)
          TAB_PLAN(i-n_particles)%vertex_ref(3,7) = (- ax3r)

          TAB_PLAN(i-n_particles)%vertex_ref(1,8) =  ax1r
          TAB_PLAN(i-n_particles)%vertex_ref(2,8) = (- ax2r)
          TAB_PLAN(i-n_particles)%vertex_ref(3,8) = (- ax3r)

        end if
        if (i==n_particles+n_walls) exit
      end if
    end do
  end if

  close(2)

  !-Reading the DOF.OUT.#
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

          ! Translating and rotating vertices
          do j=1,TAB_POLY(i)%n_vertices_cluster
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

        TAB_CONTACT(n_prpr)%id = n_prpr
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

        TAB_CONTACT(n_prpl+n_PRPRx)%id = n_prpl
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
      if ((cd == TAB_CONTACT(j)%icdent).and.(an==TAB_CONTACT(j)%ianent) .or. &
          (cd == TAB_CONTACT(j)%ianent).and.(an==TAB_CONTACT(j)%icdent)) then

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

        !  TAB_CONTACT(i)%vln = TAB_CONTACT(j)%vln
        !  TAB_CONTACT(i)%vlt = TAB_CONTACT(j)%vlt
        !  TAB_CONTACT(i)%vls = TAB_CONTACT(j)%vls

        !  TAB_CONTACT(i)%coor_ctc(:) = TAB_CONTACT(j)%coor_ctc(:)

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
    TAB_CONTACT_POLYR(cpt)%id       = TAB_CONTACT(i)%id
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

    ! Initializing the corresponding damage value as function of the contact type
    !if (TAB_CONTACT_POLYR(cpt)%type == 3) TAB_CONTACT_POLYR(cpt)%x_status = 0.33333
    !if (TAB_CONTACT_POLYR(cpt)%type == 4) TAB_CONTACT_POLYR(cpt)%x_status = 0.25

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

        ! Writing relative parameter to quantify the damage in the cohesive contact

        !if (TAB_CONTACT_POLYR(cpt)%type == 3 .and. TAB_CONTACT(j)%status == 'Mstck') then
        !  TAB_CONTACT_POLYR(cpt)%x_status = TAB_CONTACT_POLYR(cpt)%x_status + 0.33333
        !else if (TAB_CONTACT_POLYR(cpt)%type == 4 .and. TAB_CONTACT(j)%status == 'Mstck') then
        !  TAB_CONTACT_POLYR(cpt)%x_status = TAB_CONTACT_POLYR(cpt)%x_status + 0.25
        !end if
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
        do j=1,TAB_POLY(cd)%n_vertices
          dcd_ptc_contact = min(dcd_ptc_contact, sqrt((TAB_CONTACT_POLYR(cpt)%coor_ctc(1)-TAB_POLY(cd)%vertex(1,j))**2+&
                                                      (TAB_CONTACT_POLYR(cpt)%coor_ctc(2)-TAB_POLY(cd)%vertex(2,j))**2+&
                                                      (TAB_CONTACT_POLYR(cpt)%coor_ctc(3)-TAB_POLY(cd)%vertex(3,j))**2))
        end do
        do j=1,TAB_POLY(an)%n_vertices
          dan_ptc_contact = min(dan_ptc_contact, sqrt((TAB_CONTACT_POLYR(cpt)%coor_ctc(1)-TAB_POLY(an)%vertex(1,j))**2+&
                                                      (TAB_CONTACT_POLYR(cpt)%coor_ctc(2)-TAB_POLY(an)%vertex(2,j))**2+&
                                                      (TAB_CONTACT_POLYR(cpt)%coor_ctc(3)-TAB_POLY(an)%vertex(3,j))**2))
        end do

        if (min( dcd_ptc_contact , dan_ptc_contact ) < 0.0001) then
          n_fv = n_fv + 1
          TAB_CONTACT_POLYR(cpt)%type = 1
        else
          n_ss = n_ss + 1
          TAB_CONTACT_POLYR(cpt)%type = 11
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

  ! Computing the number of active contacts per particle
  ! Adding the number of contacts per particle
  TAB_POLY(:)%nctc = 0.
  TAB_POLY(:)%av_mob = 0.
  TAB_POLY(:)%av_fn = 0.

  do i=1,nb_ligneCONTACT_POLYR
    icdent = TAB_CONTACT_POLYR(i)%icdent
    ianent = TAB_CONTACT_POLYR(i)%ianent
    ! Contacts between particles
    if (TAB_CONTACT_POLYR(i)%nature /= 'PRPRx') cycle
    if ((icdent .gt. n_particles) .or. (ianent .gt. n_particles)) cycle
    ! Only active contacts
    if (TAB_CONTACT_POLYR(i)%rn .le. 0.) cycle

    TAB_POLY(icdent)%nctc = TAB_POLY(icdent)%nctc + 1
    TAB_POLY(ianent)%nctc = TAB_POLY(ianent)%nctc + 1

    if (c_fn_class == 1) then
      TAB_POLY(icdent)%av_fn = TAB_POLY(icdent)%av_fn + TAB_CONTACT_POLYR(i)%rn
      TAB_POLY(ianent)%av_fn = TAB_POLY(ianent)%av_fn + TAB_CONTACT_POLYR(i)%rn
    end if

    if (c_mobilization_class == 1) then
      TAB_POLY(icdent)%av_mob = TAB_POLY(icdent)%av_mob + &
              sqrt(TAB_CONTACT_POLYR(i)%rt**2 + TAB_CONTACT_POLYR(i)%rs**2)/(0.4*TAB_CONTACT_POLYR(i)%rn)
      TAB_POLY(ianent)%av_mob = TAB_POLY(ianent)%av_mob + &
              sqrt(TAB_CONTACT_POLYR(i)%rt**2 + TAB_CONTACT_POLYR(i)%rs**2)/(0.4*TAB_CONTACT_POLYR(i)%rn)
    end if
  end do

  ! Finding averages
  if (c_fn_class == 1) then
    do i=1, n_particles
      if (TAB_POLY(i)%nctc .gt. 0) then
        TAB_POLY(i)%av_fn = TAB_POLY(i)%av_fn/TAB_POLY(i)%nctc
      end if
    end do
  end if

  if (c_mobilization_class == 1) then
    do i=1, n_particles
      if (TAB_POLY(i)%nctc .gt. 0) then
        TAB_POLY(i)%av_mob = TAB_POLY(i)%av_mob/TAB_POLY(i)%nctc
      end if
    end do
  end if

  ! Finding the average normal force in the sample
  cpt = 0
  mean_n_force = 0.
  do i=1, nb_ligneCONTACT_POLYR
    if (TAB_CONTACT_POLYR(i)%nature /= 'PRPRx') cycle
    if (TAB_CONTACT_POLYR(i)%rn .gt. 0.) then
      mean_n_force = mean_n_force + TAB_CONTACT_POLYR(i)%rn
      cpt = cpt + 1
    end if
  end do
  mean_n_force = mean_n_force/cpt

  !print*,'Nb particles: ',n_particles
  print*,'Nb Contacts: ', nb_ligneCONTACT
  print*,'Nb real Contacts: ', nb_ligneCONTACT_POLYR
  print*,'Nb active Contacts: ', cpt_actif
  !print*,' n_ff:',n_ff
  !print*,' n_fs:',n_fs
  !print*,' n_ss ',n_ss
  !print*,' n_fv ',n_fv

  print*,'Step:', step
  print*,'Time:', time
end subroutine read_Vloc_dof_bodies


!================================================
!Computation of the box size and deformation
!================================================
subroutine calcul_deformation

  implicit none

  integer                                        :: i,j
  real*8                                         :: Xplan_max,Xplan_min,Yplan_max
  real*8                                         :: Yplan_min,Zplan_max,Zplan_min

  if (deformation_ini) then

    Xplan_max =-10000000
    Xplan_min = 10000000
    Yplan_max =-10000000
    Yplan_min = 10000000
    Zplan_max =-10000000
    Zplan_min = 10000000

    ! Finding the size of the box as function of the position of the vertices
    do i=1,n_particles
      do j=1, TAB_POLY(i)%n_vertices_cluster
        Xplan_max = max(Xplan_max, TAB_POLY(i)%vertex(1,j))
        Xplan_min = min(Xplan_min, TAB_POLY(i)%vertex(1,j))

        Yplan_max = max(Yplan_max, TAB_POLY(i)%vertex(2,j))
        Yplan_min = min(Yplan_min, TAB_POLY(i)%vertex(2,j))

        Zplan_max = max(Zplan_max, TAB_POLY(i)%vertex(3,j))
        Zplan_min = min(Zplan_min, TAB_POLY(i)%vertex(3,j))
      end do
    end do

    Long_ini = Xplan_max-Xplan_min
    Larg_ini = Yplan_max-Yplan_min
    H_ini    = Zplan_max-Zplan_min

    Depsilon_p = 0.D0
    Depsilon_q = 0.D0

    deformation_ini=.false.

  else

    Xplan_max =-10000000
    Xplan_min = 10000000
    Yplan_max =-10000000
    Yplan_min = 10000000
    Zplan_max =-10000000
    Zplan_min = 10000000

    do i=1,n_particles
      do j=1, TAB_POLY(i)%n_vertices_cluster
        Xplan_max = max(Xplan_max, TAB_POLY(i)%vertex(1,j))
        Xplan_min = min(Xplan_min, TAB_POLY(i)%vertex(1,j))

        Yplan_max = max(Yplan_max, TAB_POLY(i)%vertex(2,j))
        Yplan_min = min(Yplan_min, TAB_POLY(i)%vertex(2,j))

        Zplan_max = max(Zplan_max, TAB_POLY(i)%vertex(3,j))
        Zplan_min = min(Zplan_min, TAB_POLY(i)%vertex(3,j))
      end do
    end do
  end if

  large  = Xplan_max - Xplan_min
  width  = Yplan_max - Yplan_min
  height = Zplan_max - Zplan_min

  !print*, large
  !print*, width
  !print*, height
  !Def = (H_ini-height)/ H_ini

end subroutine calcul_deformation

!==============================================================================
! Walls positions
!==============================================================================
subroutine walls_position

  implicit none

  if (first_over_all) then
    if (n_walls == 2) then
      write(108,*) '#   time     ', '   wall 1-X  ', '   wall 1-Y  ', '   wall 1-Z  ', &
                                    '   wall 2-X  ', '   wall 2-Y  ', '   wall 2-Z  '
    else if(n_walls == 6) then
      write(108,*) '#   time     ', '   wall 1-X  ', '   wall 1-Y  ', '   wall 1-Z  ', &
                                    '   wall 2-X  ', '   wall 2-Y  ', '   wall 2-Z  ', &
                                    '   wall 3-X  ', '   wall 3-Y  ', '   wall 3-Z  ', &
                                    '   wall 4-X  ', '   wall 4-Y  ', '   wall 4-Z  ', &
                                    '   wall 5-X  ', '   wall 5-Y  ', '   wall 5-Z  ', &
                                    '   wall 6-X  ', '   wall 6-Y  ', '   wall 6-Z '
    else
      print*, "How many walls?? -- Walls position"
      stop
    end if
  endif

  if (n_walls == 2) then
    write(108,'(19(1X,E12.5))') time, TAB_PLAN(1)%center(1), TAB_PLAN(1)%center(2), TAB_PLAN(1)%center(3) &
                                    , TAB_PLAN(2)%center(1), TAB_PLAN(2)%center(2), TAB_PLAN(2)%center(3)

  else if(n_walls == 6) then
    write(108,'(19(1X,E12.5))') time, TAB_PLAN(1)%center(1), TAB_PLAN(1)%center(2), TAB_PLAN(1)%center(3) &
                                    , TAB_PLAN(2)%center(1), TAB_PLAN(2)%center(2), TAB_PLAN(2)%center(3) &
                                    , TAB_PLAN(3)%center(1), TAB_PLAN(3)%center(2), TAB_PLAN(3)%center(3) &
                                    , TAB_PLAN(4)%center(1), TAB_PLAN(4)%center(2), TAB_PLAN(4)%center(3) &
                                    , TAB_PLAN(5)%center(1), TAB_PLAN(5)%center(2), TAB_PLAN(5)%center(3) &
                                    , TAB_PLAN(6)%center(1), TAB_PLAN(6)%center(2), TAB_PLAN(6)%center(3)
  else
    print*, "How many walls?? -- Walls position"
    stop
  end if

  print*, 'Write Walls position---> Ok!'

end subroutine walls_position


!==============================================================================
! Walls forces
!==============================================================================
subroutine walls_force

  implicit none

  integer                            :: i, idw1, idw2, idw3, idw4, idw5, idw6
  real*8                             :: fx1, fx2, fx3, fx4, fx5, fx6
  real*8                             :: fy1, fy2, fy3, fy4, fy5, fy6
  real*8                             :: fz1, fz2, fz3, fz4, fz5, fz6

  ! Initializing the x forces
  fx1 = 0.
  fx2 = 0.
  fx3 = 0.
  fx4 = 0.
  fx5 = 0.
  fx6 = 0.

  ! Initializing the y forces
  fy1 = 0.
  fy2 = 0.
  fy3 = 0.
  fy4 = 0.
  fy5 = 0.
  fy6 = 0.

  ! Initializing the z forces
  fz1 = 0.
  fz2 = 0.
  fz3 = 0.
  fz4 = 0.
  fz5 = 0.
  fz6 = 0.

  ! ID for the walls
  idw1 = n_particles + 1
  idw2 = n_particles + 2
  idw3 = n_particles + 3
  idw4 = n_particles + 4
  idw5 = n_particles + 5
  idw6 = n_particles + 6

  do i=1, nb_ligneCONTACT_POLYR
    if (TAB_CONTACT_POLYR(i)%nature== 'PRPLx') then
      if(TAB_CONTACT_POLYR(i)%ianent == idw1) then
        fx1 = fx1 + TAB_CONTACT_POLYR(i)%rn*TAB_CONTACT_POLYR(i)%n(1) + &
                    TAB_CONTACT_POLYR(i)%rt*TAB_CONTACT_POLYR(i)%t(1) + &
                    TAB_CONTACT_POLYR(i)%rs*TAB_CONTACT_POLYR(i)%s(1)

        fy1 = fy1 + TAB_CONTACT_POLYR(i)%rn*TAB_CONTACT_POLYR(i)%n(2) + &
                    TAB_CONTACT_POLYR(i)%rt*TAB_CONTACT_POLYR(i)%t(2) + &
                    TAB_CONTACT_POLYR(i)%rs*TAB_CONTACT_POLYR(i)%s(2)

        fz1 = fz1 + TAB_CONTACT_POLYR(i)%rn*TAB_CONTACT_POLYR(i)%n(3) + &
                    TAB_CONTACT_POLYR(i)%rt*TAB_CONTACT_POLYR(i)%t(3) + &
                    TAB_CONTACT_POLYR(i)%rs*TAB_CONTACT_POLYR(i)%s(3)

      elseif(TAB_CONTACT_POLYR(i)%ianent == idw2) then
        fx2 = fx2 + TAB_CONTACT_POLYR(i)%rn*TAB_CONTACT_POLYR(i)%n(1) + &
                    TAB_CONTACT_POLYR(i)%rt*TAB_CONTACT_POLYR(i)%t(1) + &
                    TAB_CONTACT_POLYR(i)%rs*TAB_CONTACT_POLYR(i)%s(1)

        fy2 = fy2 + TAB_CONTACT_POLYR(i)%rn*TAB_CONTACT_POLYR(i)%n(2) + &
                    TAB_CONTACT_POLYR(i)%rt*TAB_CONTACT_POLYR(i)%t(2) + &
                    TAB_CONTACT_POLYR(i)%rs*TAB_CONTACT_POLYR(i)%s(2)

        fz2 = fz2 + TAB_CONTACT_POLYR(i)%rn*TAB_CONTACT_POLYR(i)%n(3) + &
                    TAB_CONTACT_POLYR(i)%rt*TAB_CONTACT_POLYR(i)%t(3) + &
                    TAB_CONTACT_POLYR(i)%rs*TAB_CONTACT_POLYR(i)%s(3)
      end if

      if (n_walls == 6) then
        if(TAB_CONTACT_POLYR(i)%ianent == idw3) then
          fx3 = fx3 + TAB_CONTACT_POLYR(i)%rn*TAB_CONTACT_POLYR(i)%n(1) + &
                      TAB_CONTACT_POLYR(i)%rt*TAB_CONTACT_POLYR(i)%t(1) + &
                      TAB_CONTACT_POLYR(i)%rs*TAB_CONTACT_POLYR(i)%s(1)

          fy3 = fy3 + TAB_CONTACT_POLYR(i)%rn*TAB_CONTACT_POLYR(i)%n(2) + &
                      TAB_CONTACT_POLYR(i)%rt*TAB_CONTACT_POLYR(i)%t(2) + &
                      TAB_CONTACT_POLYR(i)%rs*TAB_CONTACT_POLYR(i)%s(2)

          fz3 = fz3 + TAB_CONTACT_POLYR(i)%rn*TAB_CONTACT_POLYR(i)%n(3) + &
                      TAB_CONTACT_POLYR(i)%rt*TAB_CONTACT_POLYR(i)%t(3) + &
                      TAB_CONTACT_POLYR(i)%rs*TAB_CONTACT_POLYR(i)%s(3)

        elseif(TAB_CONTACT_POLYR(i)%ianent == idw4) then
          fx4 = fx4 + TAB_CONTACT_POLYR(i)%rn*TAB_CONTACT_POLYR(i)%n(1) + &
                      TAB_CONTACT_POLYR(i)%rt*TAB_CONTACT_POLYR(i)%t(1) + &
                      TAB_CONTACT_POLYR(i)%rs*TAB_CONTACT_POLYR(i)%s(1)

          fy4 = fy4 + TAB_CONTACT_POLYR(i)%rn*TAB_CONTACT_POLYR(i)%n(2) + &
                      TAB_CONTACT_POLYR(i)%rt*TAB_CONTACT_POLYR(i)%t(2) + &
                      TAB_CONTACT_POLYR(i)%rs*TAB_CONTACT_POLYR(i)%s(2)

          fz4 = fz4 + TAB_CONTACT_POLYR(i)%rn*TAB_CONTACT_POLYR(i)%n(3) + &
                      TAB_CONTACT_POLYR(i)%rt*TAB_CONTACT_POLYR(i)%t(3) + &
                      TAB_CONTACT_POLYR(i)%rs*TAB_CONTACT_POLYR(i)%s(3)

        elseif(TAB_CONTACT_POLYR(i)%ianent == idw5) then
          fx5 = fx5 + TAB_CONTACT_POLYR(i)%rn*TAB_CONTACT_POLYR(i)%n(1) + &
                      TAB_CONTACT_POLYR(i)%rt*TAB_CONTACT_POLYR(i)%t(1) + &
                      TAB_CONTACT_POLYR(i)%rs*TAB_CONTACT_POLYR(i)%s(1)

          fy5 = fy5 + TAB_CONTACT_POLYR(i)%rn*TAB_CONTACT_POLYR(i)%n(2) + &
                      TAB_CONTACT_POLYR(i)%rt*TAB_CONTACT_POLYR(i)%t(2) + &
                      TAB_CONTACT_POLYR(i)%rs*TAB_CONTACT_POLYR(i)%s(2)

          fz5 = fz5 + TAB_CONTACT_POLYR(i)%rn*TAB_CONTACT_POLYR(i)%n(3) + &
                      TAB_CONTACT_POLYR(i)%rt*TAB_CONTACT_POLYR(i)%t(3) + &
                      TAB_CONTACT_POLYR(i)%rs*TAB_CONTACT_POLYR(i)%s(3)

        elseif(TAB_CONTACT_POLYR(i)%ianent == idw6) then
          fx6 = fx6 + TAB_CONTACT_POLYR(i)%rn*TAB_CONTACT_POLYR(i)%n(1) + &
                      TAB_CONTACT_POLYR(i)%rt*TAB_CONTACT_POLYR(i)%t(1) + &
                      TAB_CONTACT_POLYR(i)%rs*TAB_CONTACT_POLYR(i)%s(1)

          fy6 = fy6 + TAB_CONTACT_POLYR(i)%rn*TAB_CONTACT_POLYR(i)%n(2) + &
                      TAB_CONTACT_POLYR(i)%rt*TAB_CONTACT_POLYR(i)%t(2) + &
                      TAB_CONTACT_POLYR(i)%rs*TAB_CONTACT_POLYR(i)%s(2)

          fz6 = fz6 + TAB_CONTACT_POLYR(i)%rn*TAB_CONTACT_POLYR(i)%n(3) + &
                      TAB_CONTACT_POLYR(i)%rt*TAB_CONTACT_POLYR(i)%t(3) + &
                      TAB_CONTACT_POLYR(i)%rs*TAB_CONTACT_POLYR(i)%s(3)
        end if
      end if
    end if
  end do

  if (first_over_all) then
    if (n_walls == 2) then
      write(109,*) '#   time     ', '     fx1     ', '     fy1     ', '     fz1     ', &
                                    '     fx2     ', '     fy2     ', '     fz2     '
    else if (n_walls == 6) then
      write(109,*) '#   time     ', '     fx1     ', '     fy1     ', '     fz1     ', &
                                    '     fx2     ', '     fy2     ', '     fz2     ', &
                                    '     fx3     ', '     fy3     ', '     fz3     ', &
                                    '     fx4     ', '     fy4     ', '     fz4     ', &
                                    '     fx5     ', '     fy5     ', '     fz5     ', &
                                    '     fx6     ', '     fy6     ', '     fz6    '
    else
      print*, "How many walls?? -- Walls force"
      stop
    end if
  endif

  if (n_walls == 2) then
    write(109,'(7(1X,E12.5))') time, fx1, fy1, fz1, &
                                     fx2, fy2, fz2

  else if (n_walls == 6) then
    write(109,'(19(1X,E12.5))') time, fx1, fy1, fz1, &
                                      fx2, fy2, fz2, &
                                      fx3, fy3, fz3, &
                                      fx4, fy4, fz4, &
                                      fx5, fy5, fz5, &
                                      fx5, fy6, fz6

  else
    print*, "How many walls?? -- Walls force"
    stop
  end if

  print*, 'Write Walls forces  ---> Ok!'

end subroutine walls_force

!================================================
! Computing the compacity
!================================================
subroutine compacity

  implicit none

  integer                                        ::  i
  real(kind=8)                                   ::  vol_grains, vol_box

  vol_grains = 0.D0

  do i=1,n_particles
    vol_grains = vol_grains + TAB_POLY(i)%volume
  end do

  if(first_over_all) then
    write(103,*) '#   time     ', '   height    ', '   width    ', '    large    ',' V_grain/V_box '
  endif

  vol_box = height*large*width

  write(103,'(7(1X,E12.5))') time, height, width, large, vol_grains/vol_box

  print*, 'Write Compacity     ---> Ok!'

end subroutine compacity

!================================================
! Computing the coordination
!================================================
subroutine nb_coordination

  implicit none

  integer                                  :: i, j, cd, an, float_part
  real*8                                   :: z, zc, zp, zc1, zc2, zc3

  ! Initalizing variables
  z   = 0
  zc  = 0
  zp  = 0
  float_part = 0

  zc1 = 0
  zc2 = 0
  zc3 = 0

  ! Computing the number of contacts
  do i=1, nb_ligneCONTACT_POLYR
    ! Only between particles
    if (TAB_CONTACT_POLYR(i)%nature == 'PRPLx') cycle

    ! Computing the number of contacts if they are active
    if (TAB_CONTACT_POLYR(i)%rn .gt. 0) then
      zc=zc+1
    end if
  end do

  ! Computing the number of simple and multiple contacts
  do i=1, nb_ligneCONTACT_POLYR
    ! Only contacts between particles
    if (TAB_CONTACT_POLYR(i)%nature == 'PRPLx') cycle

    ! Only active contacts
    if (TAB_CONTACT_POLYR(i)%rn .gt. 0) then

      ! One-One contacts
      ! 1 face-vertex, 11 vertex-vertex
      if (TAB_CONTACT_POLYR(i)%type == 1 .or. TAB_CONTACT_POLYR(i)%type == 11) then
        zc1 = zc1 + 1
      ! Multiple
      else if (TAB_CONTACT_POLYR(i)%type == 2) then
        zc2 = zc2 + 1
      else if (TAB_CONTACT_POLYR(i)%type .ge. 3) then
        zc3 = zc3 + 1
      else
        print*, TAB_CONTACT_POLYR(i)%type
        print*, TAB_CONTACT_POLYR(i)%icdent
        print*, TAB_CONTACT_POLYR(i)%ianent
        print*, TAB_CONTACT_POLYR(i)%id
        print*, TAB_CONTACT_POLYR(i)%gapTT
        print*, TAB_CONTACT_POLYR(i)%rn
        print*, 'the heck!?'
        stop
      end if
    end if
  end do

  ! Computing the number of floating particles
  do i=1,n_particles
    if (TAB_POLY(i)%nctc == 0) float_part = float_part + 1
  end do

  zp = n_particles - float_part

  z   = 2*zc/zp
  zc1 = 2*zc1/zp
  zc2 = 2*zc2/zp
  zc3 = 2*zc3/zp

  if (first_over_all) then
    write(100,*) '#   time     ', '      z      ', '      z1     ', '      z2     ', '      z3     '
  end if

  write(100,'(5(1X,E12.5))') time, z, zc1, zc2, zc3

  print*, 'Write Coordination  ---> Ok!'

end subroutine nb_coordination


!================================================
! Computation of q/p
!================================================
subroutine qoverp

  implicit none

  integer                                  :: i, j, cd, an, count_f
  integer                                  :: ierror,matz,lda
  real*8                                   :: S1,S2,S3
  real*8                                   :: Rtik,Rnik,Rsik,q,p, Rayon_max, ave_n_force
  real*8,dimension(3)                      :: nik,tik,sik,Lik,Fik, max_d_vert_cd, max_d_vert_an
  real*8,dimension(3)                      :: wr,wi
  real*8,dimension(3,3)                    :: Moment, M
  real*8,dimension(3,3)                    :: localframe

  real*8                                   :: pos_z_w1, pos_z_w2

  ! Initalizing
  Moment(:,:) = 0
  Lik = 0.
  Fik = 0.

  pos_z_w1 = TAB_PLAN(1)%center(3)
  pos_z_w2 = TAB_PLAN(6)%center(3)

  !print*, pos_z_w1
  !print*, pos_z_w2

  ave_n_force = 0
  count_f = 0

  do i=1,nb_ligneCONTACT_POLYR

    cd   = TAB_CONTACT_POLYR(i)%icdent
    an   = TAB_CONTACT_POLYR(i)%ianent

    ! Only for the contacts between particles. Contacts with the walls are excluded
    if (TAB_CONTACT_POLYR(i)%nature /= 'PRPRx') cycle

    ! Active contacts
    if (TAB_CONTACT_POLYR(i)%rn .le. 0.) cycle

    ! Getting rid of parasite forces
    if (TAB_CONTACT_POLYR(i)%rn .gt. 12*mean_n_force) cycle

    ! Far from the walls
    if (TAB_POLY(cd)%center(3) .lt. pos_z_w1 + 0.01*10.0) cycle
    if (TAB_POLY(cd)%center(3) .gt. pos_z_w2 - 0.01*10.0) cycle

    if (TAB_POLY(an)%center(3) .lt. pos_z_w1 + 0.01*10.0) cycle
    if (TAB_POLY(an)%center(3) .gt. pos_z_w2 - 0.01*10.0) cycle

    nik  = TAB_CONTACT_POLYR(i)%n
    tik  = TAB_CONTACT_POLYR(i)%t
    sik  = TAB_CONTACT_POLYR(i)%s
    Rtik = TAB_CONTACT_POLYR(i)%rt
    Rnik = TAB_CONTACT_POLYR(i)%rn
    Rsik = TAB_CONTACT_POLYR(i)%rs

    ! The branch vector
    Lik(1:3) = TAB_POLY(cd)%center(1:3)-TAB_POLY(an)%center(1:3)

    max_d_vert_cd(:) = 0.
    max_d_vert_an(:) = 0.

    ! Looking for the farthest vertex for each polyhedron
    ! The candidate
    do j=1, TAB_POLY(cd)%n_vertices_cluster
      max_d_vert_cd(1) = max(max_d_vert_cd(1), abs(TAB_POLY(cd)%vertex(1,j) - TAB_POLY(cd)%center(1)))
      max_d_vert_cd(2) = max(max_d_vert_cd(2), abs(TAB_POLY(cd)%vertex(2,j) - TAB_POLY(cd)%center(2)))
      max_d_vert_cd(3) = max(max_d_vert_cd(3), abs(TAB_POLY(cd)%vertex(3,j) - TAB_POLY(cd)%center(3)))
    end do

    ! The antagonist
    do j=1, TAB_POLY(an)%n_vertices_cluster
      max_d_vert_an(1) = max(max_d_vert_an(1), abs(TAB_POLY(an)%vertex(1,j) - TAB_POLY(an)%center(1)))
      max_d_vert_an(2) = max(max_d_vert_an(2), abs(TAB_POLY(an)%vertex(2,j) - TAB_POLY(an)%center(2)))
      max_d_vert_an(3) = max(max_d_vert_an(3), abs(TAB_POLY(an)%vertex(3,j) - TAB_POLY(an)%center(3)))
    end do

    ! Correction of the branch components if necessary
    ! X component
    if (Lik(1) .gt. (max_d_vert_an(1) + max_d_vert_cd(1))) then
      Lik(1) = Lik(1) - x_period*(Lik(1)/abs(Lik(1)))
    end if
    ! Y component
    if (Lik(2) .gt. (max_d_vert_an(2) + max_d_vert_cd(2))) then
      Lik(2) = Lik(2) - y_period*(Lik(2)/abs(Lik(2)))
    end if

    if (Lik(1)*nik(1) + Lik(2)*nik(2) + Lik(3)*nik(3) .gt. (TAB_POLY(cd)%n_contactors+1)*0.01) then
      print*, 'que pitos qp'
      print*, TAB_CONTACT_POLYR(i)%id
      print*, Lik
      print*, nik
      print*, max_d_vert_cd
      print*, max_d_vert_an
      !stop
      cycle
    end if

    ! The force
    Fik(1:3) = (Rnik*nik(1:3)+Rtik*tik(1:3)+Rsik*sik(1:3))

    Moment(1,1:3) = Fik(1)*Lik(1:3) + Moment(1,1:3)
    Moment(2,1:3) = Fik(2)*Lik(1:3) + Moment(2,1:3)
    Moment(3,1:3) = Fik(3)*Lik(1:3) + Moment(3,1:3)

  end do

  Moment = Moment / (height*large*width)

  lda  = 3
  matz = 1
  M = Moment

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
      do i=1,3
        if ((wr(i)<S3) .and. (wr(i)>S1)) S2=wr(i)
      end do
    end if
  end if

  q = (S3-S1)*0.5
  p = (S1+S3)*0.5

  if (first_over_all) then
    write(102,*) '#   time     ', '      S1     ', '      S2     ', '      S3     ', '     q/p     '
  end if

  write(102,'(5(1X,E12.5))') time, S1, S2, S3, q/p

  print*, 'Writing q/p          ---> Ok!'

  ! Writing the stress tensor if required
  if (qoverp_option == 1) then
    if (first_over_all) then
      write(152,*) '#   time     ', '   M( 1,1 )  ', '   M( 1,2 )  ', '   M( 1,3 )  ', &
                                    '   M( 2,1 )  ', '   M( 2,2 )  ', '   M( 2,3 )  ', &
                                    '   M( 3,1 )  ', '   M( 3,2 )  ', '   M( 3,3 )  '
    end if

    write(152,'(10(1X,E12.5))') time, Moment(1,1)/total_p_volume, Moment(1,2)/total_p_volume, Moment(1,3)/total_p_volume, &
                                      Moment(2,1)/total_p_volume, Moment(2,2)/total_p_volume, Moment(2,3)/total_p_volume, &
                                      Moment(3,1)/total_p_volume, Moment(3,2)/total_p_volume, Moment(3,3)/total_p_volume
  end if

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
  pos_z_w2 = TAB_PLAN(6)%center(3)

  ! Building the fabric tensor
  do i=1,nb_ligneCONTACT_POLYR
    ! Checking it is a contact between two bodies
    if(TAB_CONTACT_POLYR(i)%nature /= 'PRPRx') cycle

    cd   = TAB_CONTACT_POLYR(i)%icdent
    an   = TAB_CONTACT_POLYR(i)%ianent

    if ((cd .gt. n_particles) .or. (an .gt. n_particles)) cycle

    ! Counting a smaller box letting a space to the walls
    if (TAB_POLY(cd)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
    if (TAB_POLY(cd)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle

    if (TAB_POLY(an)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
    if (TAB_POLY(an)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle

    nik  = TAB_CONTACT_POLYR(i)%n
    Rnik = TAB_CONTACT_POLYR(i)%rn

    ! Checking if it is an active contact
    if (Rnik .le. 0.D0) cycle

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
  pos_z_w2 = TAB_PLAN(6)%center(3)

  ! Building the fabric forces tensor
  do i=1,nb_ligneCONTACT_POLYR

    if  (TAB_CONTACT_POLYR(i)%nature /= 'PRPRx') cycle

    cd   = TAB_CONTACT_POLYR(i)%icdent
    an   = TAB_CONTACT_POLYR(i)%ianent

    if ((cd .gt. n_particles) .or. (an .gt. n_particles)) cycle

    ! Counting a smaller box letting a space to the walls
    if (TAB_POLY(cd)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
    if (TAB_POLY(cd)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle

    if (TAB_POLY(an)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
    if (TAB_POLY(an)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle

    nik  = TAB_CONTACT_POLYR(i)%n
    tik  = TAB_CONTACT_POLYR(i)%t
    sik  = TAB_CONTACT_POLYR(i)%s
    Rnik = TAB_CONTACT_POLYR(i)%rn
    Rtik = TAB_CONTACT_POLYR(i)%rt
    Rsik = TAB_CONTACT_POLYR(i)%rs

    !print*, 'Aniso Force'
    !print*, Rtik
    !print*, Rsik

    ! Only active contacts
    if (Rnik .le. 0.D0) cycle

    ! Active contacts
    cpt = cpt + 1

    ! The force vector
    Fik(:) = Rnik*nik(:) + Rtik*tik(:) + Rsik*sik(:)

    ! Average normal force
    !av_force = av_force + Fik(1)*nik(1) + Fik(2)*nik(2) + Fik(3)*nik(3)
    ! Average total force
    av_force = av_force + sqrt(Fik(1)**2 + Fik(2)**2 + Fik(3)**2)

    ! The norm of the normal force
    norm_FN = Fik(1)*nik(1) + Fik(2)*nik(2) + Fik(3)*nik(3)

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

    !print*, 'normFT'
    !print*, norm_FT

    ! The direction of FikT
    if (norm_FT .lt. 1.0E-9) then
      dirFT = 0.
    else
      dirFT = FikT/norm_FT
    end if

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

  !print*, Fabric_N
  !print*, Fabric_F

  !print*, Fabric_N(1,1) + Fabric_N(2,2) + Fabric_N(3,3)

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
  real*8                                   :: cpt, av_length, Lnik, Ltik, Rnik, Rsik, Rtik
  real*8, dimension(3)                     :: max_d_vert_cd, max_d_vert_an
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
  Lik(:) = 0.

  pos_z_w1 = TAB_PLAN(1)%center(3)
  pos_z_w2 = TAB_PLAN(6)%center(3)

  ! Building the chi tensor with the length of branches
  do i=1,nb_ligneCONTACT_POLYR

    if  (TAB_CONTACT_POLYR(i)%nature /= 'PRPRx') cycle

    cd   = TAB_CONTACT_POLYR(i)%icdent
    an   = TAB_CONTACT_POLYR(i)%ianent
    nik  = TAB_CONTACT_POLYR(i)%n
    tik  = TAB_CONTACT_POLYR(i)%t
    sik  = TAB_CONTACT_POLYR(i)%s
    Rnik = TAB_CONTACT_POLYR(i)%rn
    Rsik = TAB_CONTACT_POLYR(i)%rs
    Rtik = TAB_CONTACT_POLYR(i)%rt

    if ((cd .gt. n_particles) .or. (an .gt. n_particles)) cycle

    ! Only active contacts
    if (Rnik .le. 0.D0) cycle

    ! Active contacts
    cpt = cpt + 1

    ! The branch
    Lik(1:3) = TAB_POLY(cd)%center(1:3)-TAB_POLY(an)%center(1:3)

    max_d_vert_cd(:) = 0.
    max_d_vert_an(:) = 0.

    ! Looking for the farthest vertex for each polyhedron
    ! The candidate
    do j=1, TAB_POLY(cd)%n_vertices_cluster
      max_d_vert_cd(1) = max(max_d_vert_cd(1), abs(TAB_POLY(cd)%vertex(1,j) - TAB_POLY(cd)%center(1)))
      max_d_vert_cd(2) = max(max_d_vert_cd(2), abs(TAB_POLY(cd)%vertex(2,j) - TAB_POLY(cd)%center(2)))
      max_d_vert_cd(3) = max(max_d_vert_cd(3), abs(TAB_POLY(cd)%vertex(3,j) - TAB_POLY(cd)%center(3)))
    end do

    ! The antagonist
    do j=1, TAB_POLY(an)%n_vertices_cluster
      max_d_vert_an(1) = max(max_d_vert_an(1), abs(TAB_POLY(an)%vertex(1,j) - TAB_POLY(an)%center(1)))
      max_d_vert_an(2) = max(max_d_vert_an(2), abs(TAB_POLY(an)%vertex(2,j) - TAB_POLY(an)%center(2)))
      max_d_vert_an(3) = max(max_d_vert_an(3), abs(TAB_POLY(an)%vertex(3,j) - TAB_POLY(an)%center(3)))
    end do

    !print*, Lik

    ! Correction of the branch components if necessary
    ! X component
    if (Lik(1) .gt. (max_d_vert_an(1) + max_d_vert_cd(1))) then
      Lik(1) = Lik(1) - x_period*(Lik(1)/abs(Lik(1)))
    end if
    ! Y component
    if (Lik(2) .gt. (max_d_vert_an(2) + max_d_vert_cd(2))) then
      Lik(2) = Lik(2) - y_period*(Lik(2)/abs(Lik(2)))
    end if

    ! The norm of the branch in the normal contact direction
    Lnik = Lik(1)*nik(1) + Lik(2)*nik(2) + Lik(3)*nik(3)

    if (Lnik .gt. (TAB_POLY(cd)%n_contactors+1.5)*0.01) then
      print*, 'que pitos'
      print*, Lnik
      print*, TAB_CONTACT_POLYR(i)%id
      print*, Lik
      print*, nik
      print*, max_d_vert_cd
      print*, max_d_vert_an
      !stop
      cycle
    end if

    ! Average branch length in the normal contact direction
    !av_length = av_length + Lnik

    ! The direction of the resultant for tangent forces and Normalizing
    ikt = tik*Rtik + sik*Rsik

    ! Computing the unitary tangential vector
    if (Rsik + Rtik .lt. 1E-9) then
      ikt = 0.
    else
      ikt = ikt/sqrt(ikt(1)**2 + ikt(2)**2 + ikt(3)**2)
    end if

    ! The branch in the tangential resultant
    Ltik = Lik(1)*ikt(1) + Lik(2)*ikt(2) + Lik(3)*ikt(3)

    ! Building the chi_n
    Fabric_N(1,1:3) = Lnik*nik(1)*nik(1:3) + Fabric_N(1,1:3)
    Fabric_N(2,1:3) = Lnik*nik(2)*nik(1:3) + Fabric_N(2,1:3)
    Fabric_N(3,1:3) = Lnik*nik(3)*nik(1:3) + Fabric_N(3,1:3)

    ! Building the tensor of tangential branches
    Fabric_T(1,1:3) = Ltik*nik(1)*ikt(1:3) + Fabric_T(1,1:3)
    Fabric_T(2,1:3) = Ltik*nik(2)*ikt(1:3) + Fabric_T(2,1:3)
    Fabric_T(3,1:3) = Ltik*nik(3)*ikt(1:3) + Fabric_T(3,1:3)
  end do

  ! Computing the average normal force
  av_length = av_length / cpt

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

  ! Initializing
  float_part = 0
  real_particles = 0

  do i=1, n_particles

    real_particles = real_particles + 1

    if (TAB_POLY(i)%nctc .lt. 1.0) then
      float_part = float_part + 1
    end if
  end do

  !print*, float_part
  !print*, real_particles

  ! Normalizing by the number of particles in the box
  float_ratio = real(float_part) / real(real_particles)

  if (first_over_all) then
    write(112,*) '#   time     ', ' float_ratio  '
  end if

  write(112,'(2(1X,E12.5))') time, float_ratio

  print*, 'Write Floating particles ---> Ok!'

end subroutine float_particles

!==============================================================================
! Computing the coordination by class
!==============================================================================
subroutine coordination_class

  implicit none

  integer                                  :: i, j, cd, an, n_bin
  real*8                                   :: r_max, r_min, interval
  real*8, dimension(:,:), allocatable      :: bins
  character(len=24)                        :: coor_c_file


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

  ! Finding the max and min equiv radius
  do i=1, n_particles
    r_min = min(r_min, TAB_POLY(i)%Rmax)
    r_max = max(r_max, TAB_POLY(i)%Rmax)
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
      if((TAB_POLY(j)%Rmax .ge. (r_min + (i-1)*interval)) .and. (TAB_POLY(j)%Rmax*0.999 .lt. (interval*i + r_min))) then 

        bins(i,2) = bins(i,2) + TAB_POLY(j)%nctc
        bins(i,3) = bins(i,3) + 1
      end if
    end do
  end do

  ! Writing
  write(113,*) '#   Class     ', '      z      '

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

  integer                                        ::  i, j, k, l_cdt, l_ant
  integer                                        ::  n_l_vertices, n_l_faces, curr_l_faces
  integer                                        ::  n_l_fields, l_counter
  real*8                                         ::  ave_rad
  real*8, dimension(3)                           ::  curr_l_vector
  logical                                        ::  dir_vtk
  character(len=8)                               ::  vtk_c_temp,  v_n_vertices
  character(:), allocatable                      ::  vtk_part, vtk_inter, vtk_counter
  character(:), allocatable                      ::  vtk_n_points

  ! Variables for the forces
  integer                                        ::  f_counter, f_counter_periodic
  real*8                                         ::  vtk_ave_force, l_force_scale, l_force
  real*8                                         ::  l_rn, l_rt, l_rs
  real*8, dimension(3)                           ::  vtk_cd_center, vtk_an_center, v_l_normal, v_l_t, v_l_s
  real*8, dimension(3)                           ::  Lik
  character(len=8)                               ::  v_n_v_forces
  character(:), allocatable                      ::  vtk_forces
  character(:), allocatable                      ::  vtk_n_v_forces


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
    ave_rad = ave_rad + TAB_POLY(i)%Rmax
  end do

  ave_rad = ave_rad/n_particles

  ! Writing the particles information
  vtk_part = './POSTPRO/VTK/rigid_' // vtk_counter // '.vtk'

  open(unit=114, file=vtk_part, status='replace')
  write(114,'(A)') '# vtk DataFile Version 3.0'
  write(114,'(A,I6)') 'RIGID ', compteur_clout
  write(114,'(A)') 'ASCII'
  write(114,'(A)') 'DATASET POLYDATA'

  ! Writing the number of vertices
  n_l_vertices = 0
  do i=1, n_particles
    n_l_vertices = n_l_vertices + TAB_POLY(i)%n_vertices_cluster
  end do

  ! Adding the vertices of the walls
  n_l_vertices = n_l_vertices + (2*8)

  write(v_n_vertices, '(I8)') n_l_vertices

  vtk_n_points = 'POINTS ' // v_n_vertices // ' float'
  write(114,'(A)') vtk_n_points

  ! Writing the coordinates of the vertices of particles and walls
  ! 3 vertices per row
  k=0
  do i=1, n_particles + n_walls
    !write(114,*) 'Particle', i
    if (i .le. n_particles) then
      do j=1, TAB_POLY(i)%n_vertices_cluster
        k = k + 1
        if(k .lt. 3) then
          if (TAB_POLY(i)%vertex(1,j) .lt. 0) then
            write(114,'(F11.8,A)', advance = 'no') TAB_POLY(i)%vertex(1,j), ' '
          else
            write(114,'(F9.7,A)', advance = 'no') TAB_POLY(i)%vertex(1,j), ' '
          end if

          if (TAB_POLY(i)%vertex(2,j) .lt. 0) then
            write(114,'(F11.8,A)', advance = 'no') TAB_POLY(i)%vertex(2,j), ' '
          else
            write(114,'(F9.7,A)', advance = 'no') TAB_POLY(i)%vertex(2,j), ' '
          end if

          if (TAB_POLY(i)%vertex(3,j) .lt. 0) then
            write(114,'(F11.8,A)', advance = 'no') TAB_POLY(i)%vertex(3,j), ' '
          else
            write(114,'(F9.7,A)', advance = 'no') TAB_POLY(i)%vertex(3,j), ' '
          end if
        else
          if (TAB_POLY(i)%vertex(1,j) .lt. 0) then
            write(114,'(F11.8,A)', advance = 'no') TAB_POLY(i)%vertex(1,j), ' '
          else
            write(114,'(F9.7,A)', advance = 'no') TAB_POLY(i)%vertex(1,j), ' '
          end if

          if (TAB_POLY(i)%vertex(2,j) .lt. 0) then
            write(114,'(F11.8,A)', advance = 'no') TAB_POLY(i)%vertex(2,j), ' '
          else
            write(114,'(F9.7,A)', advance = 'no') TAB_POLY(i)%vertex(2,j), ' '
          end if

          if (TAB_POLY(i)%vertex(3,j) .lt. 0) then
            write(114,'(F11.8,A)') TAB_POLY(i)%vertex(3,j), ' '
          else
            write(114,'(F9.7,A)') TAB_POLY(i)%vertex(3,j), ' '
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
  n_l_faces = 0

  do i=1, n_particles
    n_l_faces = n_l_faces + TAB_POLY(i)%n_faces_cluster
  end do

  ! Writing the conectivity between vertices
  ! Plus 6 faces each wall (=6*6)
  write(114,'(A)', advance='no') 'POLYGONS '
  write(114, '(2(I8,A))') n_l_faces + 12, ' ' , (n_l_faces*4)+(12*5)

  curr_l_faces = 0

  do i=1, n_particles
    !write(114,*) 'Particle', i
    do j=1, TAB_POLY(i)%n_faces_cluster
      !print*, TAB_POLY(i)%n_faces_cluster
      ! We have always triangles
      write(114, '(I1,A)', advance= 'no') 3, ' '

      ! First face... writing precisely
      if (TAB_POLY(i)%face(1,j) + curr_l_faces - 1.lt. 10) then
        write(114, '(I1)', advance = 'no') TAB_POLY(i)%face(1,j) + curr_l_faces -1
      else if (TAB_POLY(i)%face(1,j) + curr_l_faces -1 .lt. 100) then
        write(114, '(I2)', advance = 'no') TAB_POLY(i)%face(1,j) + curr_l_faces -1
      else if (TAB_POLY(i)%face(1,j) + curr_l_faces -1 .lt. 1000) then
        write(114, '(I3)', advance = 'no') TAB_POLY(i)%face(1,j) + curr_l_faces -1
      else if (TAB_POLY(i)%face(1,j) + curr_l_faces -1 .lt. 10000) then
        write(114, '(I4)', advance = 'no') TAB_POLY(i)%face(1,j) + curr_l_faces -1
      else if (TAB_POLY(i)%face(1,j) + curr_l_faces -1 .lt. 100000) then
        write(114, '(I5)', advance = 'no') TAB_POLY(i)%face(1,j) + curr_l_faces -1
      else if (TAB_POLY(i)%face(1,j) + curr_l_faces -1 .lt. 1000000) then
        write(114, '(I6)', advance = 'no') TAB_POLY(i)%face(1,j) + curr_l_faces -1
      else if (TAB_POLY(i)%face(1,j) + curr_l_faces -1 .lt. 10000000) then
        write(114, '(I7)', advance = 'no') TAB_POLY(i)%face(1,j) + curr_l_faces -1
      end if

      write(114, '(A)', advance='no') ' '

      ! Second face
      if (TAB_POLY(i)%face(2,j) + curr_l_faces -1 .lt. 10) then
        write(114, '(I1)', advance = 'no') TAB_POLY(i)%face(2,j) + curr_l_faces -1
      else if (TAB_POLY(i)%face(2,j) + curr_l_faces -1 .lt. 100) then
        write(114, '(I2)', advance = 'no') TAB_POLY(i)%face(2,j) + curr_l_faces -1
      else if (TAB_POLY(i)%face(2,j) + curr_l_faces -1 .lt. 1000) then
        write(114, '(I3)', advance = 'no') TAB_POLY(i)%face(2,j) + curr_l_faces -1
      else if (TAB_POLY(i)%face(2,j) + curr_l_faces -1 .lt. 10000) then
        write(114, '(I4)', advance = 'no') TAB_POLY(i)%face(2,j) + curr_l_faces -1
      else if (TAB_POLY(i)%face(2,j) + curr_l_faces -1 .lt. 100000) then
        write(114, '(I5)', advance = 'no') TAB_POLY(i)%face(2,j) + curr_l_faces -1
      else if (TAB_POLY(i)%face(2,j) + curr_l_faces -1 .lt. 1000000) then
        write(114, '(I6)', advance = 'no') TAB_POLY(i)%face(2,j) + curr_l_faces -1
      else if (TAB_POLY(i)%face(2,j) + curr_l_faces -1 .lt. 10000000) then
        write(114, '(I7)', advance = 'no') TAB_POLY(i)%face(2,j) + curr_l_faces -1
      end if

      ! Third face
      write(114, '(A)', advance='no') ' '

      if (TAB_POLY(i)%face(3,j) + curr_l_faces -1 .lt. 10) then
        write(114, '(I1)') TAB_POLY(i)%face(3,j) + curr_l_faces -1
      else if (TAB_POLY(i)%face(3,j) + curr_l_faces -1 .lt. 100) then
        write(114, '(I2)') TAB_POLY(i)%face(3,j) + curr_l_faces -1
      else if (TAB_POLY(i)%face(3,j) + curr_l_faces -1 .lt. 1000) then
        write(114, '(I3)') TAB_POLY(i)%face(3,j) + curr_l_faces -1
      else if (TAB_POLY(i)%face(3,j) + curr_l_faces -1 .lt. 10000) then
        write(114, '(I4)') TAB_POLY(i)%face(3,j) + curr_l_faces -1
      else if (TAB_POLY(i)%face(3,j) + curr_l_faces -1 .lt. 100000) then
        write(114, '(I5)') TAB_POLY(i)%face(3,j) + curr_l_faces -1
      else if (TAB_POLY(i)%face(3,j) + curr_l_faces -1 .lt. 1000000) then
        write(114, '(I6)') TAB_POLY(i)%face(3,j) + curr_l_faces -1
      else if (TAB_POLY(i)%face(3,j) + curr_l_faces -1 .lt. 10000000) then
        write(114, '(I7)') TAB_POLY(i)%face(3,j) + curr_l_faces -1
      end if
    end do
    curr_l_faces = curr_l_faces + TAB_POLY(i)%n_vertices_cluster
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
      do j=1, TAB_POLY(i)%n_faces_cluster
        k=k+1
        if(k .le. 9) then
          write(114, '(I6)', advance='no') l_counter
        else
          write(114, '(I6)') l_counter
          k=0
        end if
      end do
    else
      ! Six facer for each wall
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
      do j=1, TAB_POLY(i)%n_faces_cluster
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
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Dispacement
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  curr_l_vector(:) = 0.D0

  write(114,'(A,I1,I8,A)') 'Disp ', 3, (n_l_faces +12), ' float'
  k = 0
  do i=1, n_particles + n_walls
    if (i .le. n_particles) then
      curr_l_vector(:) = TAB_POLY(i)%center(:) - TAB_POLY(i)%center_ref
      do j=1, TAB_POLY(i)%n_faces_cluster
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
      curr_l_vector(1) = TAB_POLY(i)%Vx
      curr_l_vector(2) = TAB_POLY(i)%Vy
      curr_l_vector(3) = TAB_POLY(i)%Vz

      do j=1, TAB_POLY(i)%n_faces_cluster
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
      l_counter = TAB_POLY(i)%nctc

      do j=1, TAB_POLY(i)%n_faces_cluster
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
  !print*, nb_ligneCONTACT_POLYR
  do i=1, nb_ligneCONTACT_POLYR
    !print*, i
    if(TAB_CONTACT_POLYR(i)%nature /= 'PRPRx') cycle
    ! Removing contacts with walls
    l_cdt = TAB_CONTACT_POLYR(i)%icdent
    l_ant = TAB_CONTACT_POLYR(i)%ianent

    !print*, l_cdt, l_ant

    if (l_cdt .gt. n_particles .or. l_ant .gt. n_particles) cycle

    if (TAB_CONTACT_POLYR(i)%rn .le. 0) cycle

    Lik(:) = TAB_POLY(l_cdt)%center(:)-TAB_POLY(l_ant)%center(:)

    ! Forces in the periodic zone should be counted separately
    if (sqrt(Lik(1)**2 + Lik(2)**2 + Lik(3)**2) .gt. (TAB_POLY(l_cdt)%Rmax + TAB_POLY(l_ant)%Rmax)) then
      f_counter_periodic = f_counter_periodic + 1
    else
      f_counter = f_counter + 1
    end if
    vtk_ave_force = vtk_ave_force + TAB_CONTACT_POLYR(i)%rn
  end do

  vtk_ave_force = vtk_ave_force/(f_counter+f_counter_periodic)

  write(v_n_v_forces, '(I8)') (f_counter + f_counter_periodic*2) * 8

  vtk_n_v_forces = 'POINTS' // v_n_v_forces // ' float'
  write(114,'(A)') vtk_n_v_forces

  ! Force scale parameter
  ! The scale is here!
  ! 30% of average radius
  l_force_scale = (ave_rad)*0.03

  ! Writing
  k=0
  do i=1, nb_ligneCONTACT_POLYR

    if (TAB_CONTACT_POLYR(i)%nature /= 'PRPRx') cycle
    ! The particles ids
    l_cdt = TAB_CONTACT_POLYR(i)%icdent
    l_ant = TAB_CONTACT_POLYR(i)%ianent

    if (l_cdt .gt. n_particles .or. l_ant .gt. n_particles) cycle

    if (TAB_CONTACT_POLYR(i)%rn .le. 0) cycle

    Lik(:) = TAB_POLY(l_cdt)%center(:)-TAB_POLY(l_ant)%center(:)

    ! Being careful with particles in the periodic borders
    ! Periodic case
    if (sqrt(Lik(1)**2 + Lik(2)**2 + Lik(3)**2) .gt. (TAB_POLY(l_cdt)%Rmax + TAB_POLY(l_ant)%Rmax)) then
      ! The center of the candidate and computing the center of the Antagonist
      vtk_cd_center(:) = TAB_POLY(l_cdt)%center(:)

      ! The normal and tangential vectors
      v_l_normal(:) = TAB_CONTACT_POLYR(i)%n(:)
      v_l_t(:) = TAB_CONTACT_POLYR(i)%t(:)
      v_l_s(:) = TAB_CONTACT_POLYR(i)%s(:)

      ! The center of the Antagonist
      vtk_an_center(:) = TAB_POLY(l_cdt)%center(:) - (TAB_POLY(l_cdt)%Rmax+TAB_POLY(l_ant)%Rmax)*v_l_normal(:)

      ! The forces
      l_rn = TAB_CONTACT_POLYR(i)%rn
      l_rt = TAB_CONTACT_POLYR(i)%rt
      l_rs = TAB_CONTACT_POLYR(i)%rs

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
      vtk_an_center(:) = TAB_POLY(l_ant)%center(:)

      ! The normal and tangential vectors
      v_l_normal(:) = TAB_CONTACT_POLYR(i)%n(:)
      v_l_t(:) = TAB_CONTACT_POLYR(i)%t(:)
      v_l_s(:) = TAB_CONTACT_POLYR(i)%s(:)

      ! The center of the Antagonist
      vtk_cd_center(:) = TAB_POLY(l_ant)%center(:) + (TAB_POLY(l_cdt)%Rmax+TAB_POLY(l_ant)%Rmax)*v_l_normal(:)

      ! The forces
      l_rn = TAB_CONTACT_POLYR(i)%rn
      l_rt = TAB_CONTACT_POLYR(i)%rt
      l_rs = TAB_CONTACT_POLYR(i)%rs

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
  n_l_fields = 3                          ! They are: RN, RT, Mobilization
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
  do i=1, nb_ligneCONTACT_POLYR
    if (TAB_CONTACT_POLYR(i)%nature /= 'PRPRx') cycle

    ! The particles ids
    l_cdt = TAB_CONTACT_POLYR(i)%icdent
    l_ant = TAB_CONTACT_POLYR(i)%ianent

    if (l_cdt .gt. n_particles .or. l_ant .gt. n_particles) cycle

    if (TAB_CONTACT_POLYR(i)%rn .le. 0) cycle

    Lik(:) = TAB_POLY(l_cdt)%center(:)-TAB_POLY(l_ant)%center(:)

    ! Being careful with the periodic case
    if (sqrt(Lik(1)**2 + Lik(2)**2 + Lik(3)**2) .gt. (TAB_POLY(l_cdt)%Rmax + TAB_POLY(l_ant)%Rmax)) then
      do j=1, 12
        k=k+1
        if(k .le. 9) then
          write(114, '(F15.7)', advance='no') TAB_CONTACT_POLYR(i)%rn
        else
          write(114, '(F15.7)') TAB_CONTACT_POLYR(i)%rn
          k=0
        end if
      end do
    else
      do j=1, 6
        k=k+1
        if(k .le. 9) then
          write(114, '(F15.7)', advance='no') TAB_CONTACT_POLYR(i)%rn
        else
          write(114, '(F15.7)') TAB_CONTACT_POLYR(i)%rn
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
  do i=1, nb_ligneCONTACT_POLYR
    if (TAB_CONTACT_POLYR(i)%nature /= 'PRPRx') cycle

    ! The particles ids
    l_cdt = TAB_CONTACT_POLYR(i)%icdent
    l_ant = TAB_CONTACT_POLYR(i)%ianent

    if (l_cdt .gt. n_particles .or. l_ant .gt. n_particles) cycle

    if (TAB_CONTACT_POLYR(i)%rn .le. 0) cycle

    Lik(:) = TAB_POLY(l_cdt)%center(:)-TAB_POLY(l_ant)%center(:)

    ! Being careful with the periodic case
    if (sqrt(Lik(1)**2 + Lik(2)**2 + Lik(3)**2) .gt. (TAB_POLY(l_cdt)%Rmax + TAB_POLY(l_ant)%Rmax)) then
      do j=1, 12
        k=k+1
        if(k .le. 9) then
          write(114, '(F15.7)', advance='no') sqrt(TAB_CONTACT_POLYR(i)%rt**2 + TAB_CONTACT_POLYR(i)%rs**2)
        else
          write(114, '(F15.7)') sqrt(TAB_CONTACT_POLYR(i)%rt**2 + TAB_CONTACT_POLYR(i)%rs**2)
          k=0
        end if
      end do
    else
      do j=1, 6
        k=k+1
        if(k .le. 9) then
          write(114, '(F15.7)', advance='no') sqrt(TAB_CONTACT_POLYR(i)%rt**2 + TAB_CONTACT_POLYR(i)%rs**2)
        else
          write(114, '(F15.7)') sqrt(TAB_CONTACT_POLYR(i)%rt**2 + TAB_CONTACT_POLYR(i)%rs**2)
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
  do i=1, nb_ligneCONTACT_POLYR
    if (TAB_CONTACT_POLYR(i)%nature /= 'PRPRx') cycle

    l_cdt = TAB_CONTACT_POLYR(i)%icdent
    l_ant = TAB_CONTACT_POLYR(i)%ianent

    if (l_cdt .gt. n_particles .or. l_ant .gt. n_particles) cycle

    Lik(:) = TAB_POLY(l_cdt)%center(:)-TAB_POLY(l_ant)%center(:)

    if (TAB_CONTACT_POLYR(i)%rn .le. 0) cycle

    ! Attention.. periodic case
    if (sqrt(Lik(1)**2 + Lik(2)**2 + Lik(3)**2) .gt. (TAB_POLY(l_cdt)%Rmax + TAB_POLY(l_ant)%Rmax)) then
      do j=1, 12
        k=k+1
        if(k .le. 9) then
          ! Careful--> I'm writing the coeficient of friction explicitly
          write(114, '(F15.7)', advance='no') sqrt(TAB_CONTACT_POLYR(i)%rt**2 +  &
                            TAB_CONTACT_POLYR(i)%rs**2)/(0.4*TAB_CONTACT_POLYR(i)%rn)
        else
          write(114, '(F15.7)') sqrt(TAB_CONTACT_POLYR(i)%rt**2 + &
                            TAB_CONTACT_POLYR(i)%rs**2)/(0.4*TAB_CONTACT_POLYR(i)%rn)
          k=0
        end if
      end do
    else
      do j=1, 6
        k=k+1
        if(k .le. 9) then
          ! Careful--> I'm writing the coeficient of friction explicitly
          write(114, '(F15.7)', advance='no') sqrt(TAB_CONTACT_POLYR(i)%rt**2 + &
                            TAB_CONTACT_POLYR(i)%rs**2)/(0.4*TAB_CONTACT_POLYR(i)%rn)
        else
          write(114, '(F15.7)') sqrt(TAB_CONTACT_POLYR(i)%rt**2 + &
                            TAB_CONTACT_POLYR(i)%rs**2)/(0.4*TAB_CONTACT_POLYR(i)%rn)
          k=0
        end if
      end do
    end if
  end do
  ! And jump
  write(114, '(A)') ' '

  close(114)

  print*, 'Drawing in vtk       ---> Ok!'

end subroutine draw


!==============================================================================
! List of normal forces
!==============================================================================
subroutine list_n_forces

  implicit none

  integer                                  :: i, cd, an
  real*8                                   :: dist_ref
  real*8, dimension(3)                     :: Lik
  character(len=24)                        :: l_file

  real*8                                   :: pos_z_w1, pos_z_w2


  pos_z_w1 = TAB_PLAN(1)%center(3)
  pos_z_w2 = TAB_PLAN(6)%center(3)

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

  do i=1,nb_ligneCONTACT_POLYR
    ! Contacts between particles
    if (TAB_CONTACT_POLYR(i)%nature /= 'PRPRx') cycle

    ! Active contacts
    if (TAB_CONTACT_POLYR(i)%rn .le. 0.0) cycle

    ! Getting rid of parasite forces
    if (TAB_CONTACT_POLYR(i)%rn .gt. 20*mean_n_force) cycle

    cd   = TAB_CONTACT_POLYR(i)%icdent
    an   = TAB_CONTACT_POLYR(i)%ianent

    if ((cd .gt. n_particles) .or. (an .gt. n_particles)) cycle


    dist_ref = (((TAB_POLY(cd)%n_contactors - 1)/3) + 1) * (2*0.01)
    ! Counting a smaller box letting a space to the walls
    if (TAB_POLY(cd)%center(3) .lt. pos_z_w1 + dist_ref) cycle
    if (TAB_POLY(cd)%center(3) .gt. pos_z_w2 - dist_ref) cycle

    if (TAB_POLY(an)%center(3) .lt. pos_z_w1 + dist_ref) cycle
    if (TAB_POLY(an)%center(3) .gt. pos_z_w2 - dist_ref) cycle

    Lik(:) = TAB_POLY(cd)%center(:)-TAB_POLY(an)%center(:)

    ! Not joining contacts in the periodic zones
    if (sqrt(Lik(1)**2 + Lik(2)**2 + Lik(3)**2) .gt. (TAB_POLY(cd)%Rmax + TAB_POLY(an)%Rmax)) cycle

    write(115, '(E15.7)') TAB_CONTACT_POLYR(i)%rn
  end do

  close(115)

  print*, 'Write List FN     ---> Ok!'

end subroutine list_n_forces


!==============================================================================
! List of branches
!==============================================================================
subroutine list_branches

  implicit none

  integer                                  :: i, cd, an
  real*8, dimension(3)                     :: Lik
  character(len=24)                        :: l_file

  real*8                                   :: pos_z_w1, pos_z_w2


  pos_z_w1 = TAB_PLAN(1)%center(3)
  pos_z_w2 = TAB_PLAN(6)%center(3)

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

  do i=1,nb_ligneCONTACT_POLYR
    ! Contacts between particles
    if (TAB_CONTACT_POLYR(i)%nature /= 'PRPRx') cycle

    ! Active contacts
    if (TAB_CONTACT_POLYR(i)%rn .le. 0.0) cycle

    cd   = TAB_CONTACT_POLYR(i)%icdent
    an   = TAB_CONTACT_POLYR(i)%ianent

    if ((cd .gt. n_particles) .or. (an .gt. n_particles)) cycle

    ! Counting a smaller box letting a space to the walls
    if (TAB_POLY(cd)%center(3) .lt. pos_z_w1 + (TAB_POLY(cd)%n_contactors + 1.5)*0.01) cycle
    if (TAB_POLY(cd)%center(3) .gt. pos_z_w2 - (TAB_POLY(cd)%n_contactors + 1.5)*0.01) cycle

    if (TAB_POLY(an)%center(3) .lt. pos_z_w1 + (TAB_POLY(an)%n_contactors + 1.5)*0.01) cycle
    if (TAB_POLY(an)%center(3) .gt. pos_z_w2 - (TAB_POLY(an)%n_contactors + 1.5)*0.01) cycle

    Lik(:) = TAB_POLY(cd)%center(:)-TAB_POLY(an)%center(:)

    ! Not joining contacts in the periodic zones
    if (sqrt(Lik(1)**2 + Lik(2)**2 + Lik(3)**2) .gt. (TAB_POLY(cd)%Rmax + TAB_POLY(an)%Rmax)) cycle

    write(116, '(E15.7)') sqrt(Lik(1)**2 + Lik(2)**2 + Lik(3)**2)
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
  pos_z_w2 = TAB_PLAN(6)%center(3)


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

  write(117,*) '#  Mobil    '

  do i=1,nb_ligneCONTACT_POLYR
    ! Contacts between particles
    if (TAB_CONTACT_POLYR(i)%nature /= 'PRPRx') cycle

    ! Active contacts
    if (TAB_CONTACT_POLYR(i)%rn .le. 0.0) cycle

    cd   = TAB_CONTACT_POLYR(i)%icdent
    an   = TAB_CONTACT_POLYR(i)%ianent

    if ((cd .gt. n_particles) .or. (an .gt. n_particles)) cycle

    ! Counting a smaller box letting a space to the walls
    if (TAB_POLY(cd)%center(3) .lt. pos_z_w1 + (TAB_POLY(cd)%n_contactors + 1.5)*0.01) cycle
    if (TAB_POLY(cd)%center(3) .gt. pos_z_w2 - (TAB_POLY(cd)%n_contactors + 1.5)*0.01) cycle

    if (TAB_POLY(an)%center(3) .lt. pos_z_w1 + (TAB_POLY(an)%n_contactors + 1.5)*0.01) cycle
    if (TAB_POLY(an)%center(3) .gt. pos_z_w2 - (TAB_POLY(an)%n_contactors + 1.5)*0.01) cycle

    Lik(:) = TAB_POLY(cd)%center(:)-TAB_POLY(an)%center(:)

    ! Not joining contacts in the periodic zones
    if (sqrt(Lik(1)**2 + Lik(2)**2 + Lik(3)**2) .gt. (TAB_POLY(cd)%Rmax + TAB_POLY(an)%Rmax)) cycle

    ! Attention ---> coeficient of friction explicitly written
    write(117, '(E15.7)') sqrt(TAB_CONTACT_POLYR(i)%rs**2 + TAB_CONTACT_POLYR(i)%rt**2)/(0.4*TAB_CONTACT_POLYR(i)%rn)
  end do

  close(117)

  print*, 'Write List mobilization ---> Ok!'

end subroutine list_mobilization

!==============================================================================
! Computing the normal force by class
!==============================================================================
subroutine fn_class

  implicit none

  integer                                  :: i, j, n_bin
  real*8                                   :: r_max, r_min, interval
  real*8, dimension(:,:), allocatable      :: bins
  character(len=24)                        :: coor_c_file


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
    r_min = min(r_min, TAB_POLY(i)%Rmax)
    r_max = max(r_max, TAB_POLY(i)%Rmax)
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
      if((TAB_POLY(j)%Rmax .ge. (r_min + (i-1)*interval)) .and. (TAB_POLY(j)%Rmax*0.999 .lt. (interval*i + r_min))) then

        bins(i,2) = bins(i,2) + TAB_POLY(j)%av_fn
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
  
  integer                                  :: i, j, n_bin
  real*8                                   :: r_max, r_min, interval
  real*8, dimension(:,:), allocatable      :: bins
  character(len=24)                        :: coor_c_file

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
    r_min = min(r_min, TAB_POLY(i)%Rmax)
    r_max = max(r_max, TAB_POLY(i)%Rmax)
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
      if((TAB_POLY(j)%Rmax .ge. (r_min + (i-1)*interval)) .and. (TAB_POLY(j)%Rmax*0.999 .lt. (interval*i + r_min))) then 

        bins(i,2) = bins(i,2) + TAB_POLY(j)%av_mob
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
  pos_z_w2 = TAB_PLAN(6)%center(3)

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
  contotal = 0
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
    do j=1,nb_ligneCONTACT_POLYR
      cd = TAB_CONTACT_POLYR(j)%icdent
      an = TAB_CONTACT_POLYR(j)%ianent
      ! Only between discs
      if(TAB_CONTACT_POLYR(j)%nature /= 'PRPRx') cycle
      ! Only active contacts
      if(TAB_CONTACT_POLYR(j)%rn .le. 0.D0) cycle

      if (TAB_POLY(cd)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
      if (TAB_POLY(cd)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle
      if (TAB_POLY(an)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
      if (TAB_POLY(an)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle

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

  ! Computing the number of active contacts
  do i=1, ninter
    contotal = contotal + theta_interval(i,2)
  end do

  ! Normalizing
  theta_interval(:,3) = theta_interval(:,2)/(2*contotal)

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
  pos_z_w2 = TAB_PLAN(6)%center(3)

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
    do j=1,nb_ligneCONTACT_POLYR
      cd = TAB_CONTACT_POLYR(j)%icdent
      an = TAB_CONTACT_POLYR(j)%ianent
      ! Only between discs
      if(TAB_CONTACT_POLYR(j)%nature /= 'PRPRx') cycle
      ! Only active contacts
      if(TAB_CONTACT_POLYR(j)%rn .le. 0.D0) cycle
      if (TAB_POLY(cd)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
      if (TAB_POLY(cd)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle
      if (TAB_POLY(an)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
      if (TAB_POLY(an)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle

      ! The branch vector
      Lik(:) = TAB_POLY(cd)%center(:)-TAB_POLY(an)%center(:)

      ! Not joining contacts in the periodic zones
      if (sqrt(Lik(1)**2 + Lik(2)**2 + Lik(3)**2) .gt. (TAB_POLY(cd)%Rmax + TAB_POLY(an)%Rmax)) cycle

      ! The normal vector
      nik(1)  = TAB_CONTACT_POLYR(j)%n(1)
      nik(2)  = TAB_CONTACT_POLYR(j)%n(2)
      nik(3)  = TAB_CONTACT_POLYR(j)%n(3)

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
  pos_z_w2 = TAB_PLAN(6)%center(3)

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
    do j=1,nb_ligneCONTACT_POLYR
      cd = TAB_CONTACT_POLYR(j)%icdent
      an = TAB_CONTACT_POLYR(j)%ianent
      ! Only polyhedra
      if(TAB_CONTACT_POLYR(j)%nature /= 'PRPRx') cycle
      ! Only active contacts
      if(TAB_CONTACT_POLYR(j)%rn .le. 0.D0) cycle
      if (TAB_POLY(cd)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
      if (TAB_POLY(cd)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle
      if (TAB_POLY(an)%center(3) .lt. pos_z_w1 + 0.01*5.0) cycle
      if (TAB_POLY(an)%center(3) .gt. pos_z_w2 - 0.01*5.0) cycle

      Rnik = TAB_CONTACT_POLYR(j)%rn
      Rtik = TAB_CONTACT_POLYR(j)%rt
      Rsik = TAB_CONTACT_POLYR(j)%rs

      ! The normal vector
      nik = TAB_CONTACT_POLYR(j)%n

      ! The tangential vector
      tik = TAB_CONTACT_POLYR(j)%t
      sik = TAB_CONTACT_POLYR(j)%s
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
! Gap Evolution
!==============================================================================
subroutine gap_evol

  implicit none

  integer                                      :: i, l_counter
  real(kind=8)                                 :: max_gap, mean_gap, rela_gap

  ! Initializing
  l_counter = 0
  max_gap = 0.0
  mean_gap = 0.0
  rela_gap = 0.0


  ! Computing the mean gap and looking for the maximum gap in the sample between particles
  do i=1, nb_ligneCONTACT
    if (TAB_CONTACT(i)%nature /= 'PRPRx') cycle

    max_gap = min(max_gap, TAB_CONTACT(i)%gapTT)

    if (TAB_CONTACT(i)%gapTT .le. 0) then
      l_counter = l_counter + 1
      mean_gap = mean_gap + TAB_CONTACT(i)%gapTT
    end if
  end do

  do i=1, n_particles
    rela_gap = rela_gap + TAB_POLY(i)%Rmax
  end do

  mean_gap = mean_gap/l_counter
  rela_gap = rela_gap/n_particles

  ! Writing the heading
  if (first_over_all) then
    write(123,*) '   # Time   ', '    max_gap   ', '   mean_gap  ', '  mean_r_gap '
  end if

  ! Writing
  write(123,'(4(1X,E12.5))') time, abs(max_gap), abs(mean_gap), abs(mean_gap)/rela_gap

  print*, 'Gap evolution    ---> Ok!'

end subroutine gap_evol


!================================================
! Closing files
!================================================
subroutine close_all

  implicit none

  ! Variables
  integer                                        :: i

  if (c_coordination                        == 1) close (100)
  !if (c_vitesse_moyenne                     == 1) close (101)
  if (c_compacity                           == 1) close (103)
  if (c_qoverp                              == 1) close (103)
  if (c_anisotropy_contact                  == 1) close (104)
  if (c_anisotropy_force                    == 1) close (105)
  if (c_anisotropy_branch                   == 1) close (106)
  if (c_walls_position                      == 1) close (108)
  if (c_walls_force                         == 1) close (109)
  if (c_float_particles                     == 1) close (112)
  if (c_gap_evol                            == 1) close (123)
  
  ! Deallocating the polyhedra container
  do,i=1,n_particles
    deallocate(TAB_POLY(i)%vertex_ref)
    deallocate(TAB_POLY(i)%vertex)
    deallocate(TAB_POLY(i)%face)
  end do
  
  deallocate(TAB_POLY)
  
  print*,'--> Closing all <--'

  stop
end subroutine close_all


!==============================================================================
!==============================================================================
!==============================================================================
!==============================================================================
!-ADITIONAL ROUTINES USED IN THE SUBROUTINES ABOVE

!---------------------------------------------------
subroutine rg(lda, n, a, wr, wi, matz, z, ierror)
!  RG finds the eigenvalues and eigenvectors of a real(kind=8) general matrix.
!  
!  Modified:
!  01 May 2000
!  
!  Parameters:
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
  
  ! Variables
  integer                                        :: lda
  integer                                        :: n
  integer                                        :: nm
  
  real(kind=8)                                   :: a(lda,n)
  real(kind=8)                                   :: fv1(n)
  integer                                        :: ierror
  integer                                        :: is1
  integer                                        :: is2
  integer                                        :: iv1(n)
  integer                                        :: matz
  real(kind=8)                                   :: wi(n)
  real(kind=8)                                   :: wr(n)
  real(kind=8)                                   :: z(lda,n)
  
  
  ierror = 0
  
  if (n > lda) then
    ierror = 10 * n
    return
  end if
  
  ! Balance the matrix.
  call balanc(lda, n, a, is1, is2, fv1)
  
  ! Put the matrix into upper Hessenberg form.
  call elmhes( lda, n, is1, is2, a, iv1)
  
  if (matz == 0) then
    call hqr(lda, n, is1, is2, a, wr, wi, ierror)
    
    if (ierror /= 0) then
      return
    end if
    
  else
    call eltran(lda, n, is1, is2, a, iv1, z)
    call hqr2(lda, n, is1, is2, a, wr, wi, z, ierror)
    
    if (ierror /= 0) then
      return
    end if
    
    call balbak(lda, n, is1, is2, fv1, n, z)
  end if
end subroutine rg



subroutine balanc(nm, n, a, low, igh, scale)
!! BALANC balances a real(kind=8) matrix before eigenvalue calculations.
!
!  Discussion:
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
  
  ! Variables
  integer                                        :: nm
  integer                                        :: n
  
  real(kind=8)                                   :: a(nm,n)
  real(kind=8)                                   :: b2
  real(kind=8)                                   :: c
  real(kind=8)                                   :: f
  real(kind=8)                                   :: g
  
  integer                                        :: i
  integer                                        :: iexc
  integer                                        :: igh
  integer                                        :: j
  integer                                        :: k
  integer                                        :: l
  integer                                        :: low
  integer                                        :: m
  logical                                        :: noconv
  real(kind=8)                                   :: r
  real(kind=8), parameter                        :: radix = 16.0D+00
  real(kind=8) s
  real(kind=8) scale(n)
  
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

subroutine init_random_seed()

      INTEGER :: i, n, clock
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed

      CALL RANDOM_SEED(size = n)
      ALLOCATE(seed(n))

      CALL SYSTEM_CLOCK(COUNT=clock)

      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
      CALL RANDOM_SEED(PUT = seed)

      DEALLOCATE(seed)
end


end program post3D_polyhe