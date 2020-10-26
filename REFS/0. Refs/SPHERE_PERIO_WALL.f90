program grandeur3d

implicit none

type  ::  T_CORPS
real(kind=8)                               ::  Rmax,volume,Ax1,Ax2,Ax3
real(kind=8),dimension(3)                  ::  centre,centre_ref
logical                                    ::  inbox
real(kind=8)                               ::  Vx,Vy,Vz,Vrx,Vry,Vrz
real(kind=8)                               ::  Wx,Wy,Wz
real(kind=8)                               ::  Norme_V,Norme_Vr,I1,I2,I3
real(kind=8),dimension(3,3)                ::  momentI,Rot,moment_part,moment_part1,moment_part2,&
                                               moment_part3,moment_part4,moment_part5,moment_part6,moment_part7,&
                                               moment_part8,moment_part9
real(kind=8)                               ::  nctc,nctcslide,nctcstick,nb_ctc,Norm_F,ctc_float
character(len=5)                           ::  behav,color
integer                                    ::  num
real(kind=8),dimension(3,3)                ::  Fabric,Fabric_N,Fabric_T,Fabric_T2
end type T_CORPS

type  ::  T_CONTACT                                                   
integer                                   ::  icdent,ianent,type        
real(kind=8),dimension(3)                 ::  n,t,s                 
real(kind=8),dimension(3)                 ::  coor_ctc              
real(kind=8)                              ::  rn,rt,rs,vln,vls,vlt,gap
integer                                   ::  sect
character(len=5)                          ::  nature,statut  
logical                                   ::  inbox,cycle_gmv,contact_perdu
end type T_CONTACT


type(T_CORPS),dimension(:),allocatable        ::  TAB_SPHERE,TAB_PLAN
type(T_CONTACT),dimension(:),allocatable      ::  TAB_CONTACT

!DEFINITIONS Vloc-Dof-Bodies
!================================================================
integer                                     ::  nb_paroi=0,nbParticules=0
integer                                     ::  compteur_clout=-1,all_pas=0
integer                                     ::  nb_face,nb_sommet
logical                                     ::  the_first_time=.true., fin_post=.false. 
integer                                     ::  pas,nb_ligneCONTACT
real(kind=8)                                ::  temps,Mean_total_Normal_force

!DEFINITIONS BOITE
!================================================================
real(kind=8),dimension(3,8)                 ::  tab_sommet_box
real(kind=8)                                ::  Def,H_ini,Long_ini,Larg_ini,Def_H,step_pdf_start,step_pdf_stop,&
                                                step_ksi_reseau,stop_ksi_reseau,Ep_avant=0,Eq_avant=0,Ep,Eq,Depsilon_p,Depsilon_q,&
                                                pressure_in_sample=0,pressure_c_in_sample = 0.D0,&
                                                taux_def=0,sin_de_phi=0,partial_pressure_in_sample=0
logical                                     ::  deformation_ini = .true.,deformation_clout=.false.

!MOYENNE POUR KSI-RESEAU
!================================================================
real(kind=8),dimension(:,:,:),allocatable   ::  Mean_Ksi_reseau 
integer                                     ::  pas_ksi = 0

integer                                     ::  quel_pas_de_temps =0


!CALCUL DES PROFILS 
!================================================================
integer                                     ::  debut_profils=0,fin=0
logical                                     ::  the_first_time_profils

!CALCUL DU PROFIL DE VITESSE MOYEN
!================================================================
type  ::  T_step_velocity
  real(kind=8),dimension(:,:),allocatable             ::  tab_vmoyen
  real(kind=8),dimension(:,:),allocatable             ::  tab_vrmoyen
end type T_step_velocity
type(T_step_velocity),dimension(:),allocatable        ::  TAB_vitesse

!CALCUL DU PROFIL DE CONTRAINTE MOYEN
!================================================================
type  ::  T_sigma
  real(kind=8),dimension(3,3)                       :: sigma
end type
type  ::  T_step_contrainte
  type(T_sigma),dimension(:),allocatable            ::  tab_sigma
end type 
type(T_step_contrainte),dimension(:),allocatable      ::  TAB_contrainte
!CALCUL DU PROFIL DE COMPACITY MOYEN
!================================================================
type  ::  T_step_compacity
  real(kind=8),dimension(:,:),allocatable             ::  tab_compacity_z
end type T_step_compacity
type(T_step_compacity),dimension(:),allocatable       ::  TAB_compacity


!----------------Pour Persiste Contact---------------------------
logical                                     :: first_list = .true.
type(T_CONTACT),dimension(:),pointer        :: TAB_CONTACT_INITIAL
integer                                     :: nb_ligneCONTACT_first=0,debut_persistance=0


!COMMANDES D APPELLE
!================================================================
character(len=30)                           ::  command
real(kind=8)                                ::  precision, PI = 3.14159265359,Long=0,Larg=0,Haut=0,step_all_gmv

integer                                     ::  type_box=0,calcul_vitesse_moyenne=0,calcul_coordination=0,&
                                                calcul_qsurp=0,calcul_compacity=0,calcul_contact_anisotropies=0,&
                                                calcul_write_all_gmv=0,calcul_pdf=0,calcul_ksi_reseau=0,&
                                                calcul_contact_anisotropies_locale=0,calcul_construct_bodies=0,&
                                                calcul_profils=0,calcul_persiste_contacts=0


integer                                     ::  calcul_contact_orientation=0,calcul_contact_orientation3D=0,&
                                                calcul_spatial_correlation=0,calcul_radial_distribution=0,&
                                                calcul_correlation_function=0

!================================================================
!Lecture du fichier de commande
!================================================================
print*,'-------------------------------------------'
print*,'!                                         !'
print*,'!      Module de Post-traitement 3D       !'
print*,'!             Emilien Azema               !'
print*,'!                                         !'
print*,'-------------------------------------------'
print*,'Lecture des commandes'
open(unit=1,file='BOITE.DAT',status='old')
do 
  read(1,'(A30)') command
  print*,command
  if (command=='ALL BOX                      :') then
    type_box=0
    cycle
  end if

  if (command=='CALCUL DEFORMATION           :') then
    deformation_clout = .true.
    cycle
  end if

  if (command=='VLOC_DOF_BODIES              :') then
    read(1,*) nb_paroi
    read(1,*) nbParticules
    read(1,*) all_pas
    cycle
  end if
  
  if (command=='NOMBRE DE COORDINATION       :') then
    calcul_coordination=1
    open (unit=100,file='NOMBRE_COORDINATION.DAT',status='replace')     
    cycle
  end if

  if (command=='VITESSE MOYENNE              :') then
    calcul_vitesse_moyenne=1
    open (unit=101,file='VITESSE_MOYENNE.DAT',status='replace')     
    cycle
  end if

  if (command=='Q/P                          :') then
    calcul_qsurp=1
    open (unit=102,file='QsurP.DAT',status='replace')     
!    open (unit=1003,file='QsurP_moment.DAT',status='replace')
    cycle
  end if

  if (command=='COMPACITY                    :') then
    calcul_compacity=1
    open (unit=103,file='COMPACITY.DAT',status='replace')     
    cycle
  end if

  if (command=='CONTACT ANISOTROPIE          :') then
    calcul_contact_anisotropies=1
    open (unit=104,file='CONTACT_ANISOTROPIES.DAT',status='replace')     
!    open (unit=1004,file='CONTACT_ANISOTROPIES_TENSOR.DAT',status='replace')     
!    open (unit=1005,file='CONTACT_ANISOTROPIES_FORCE.DAT',status='replace')     
    cycle
  end if

  if (command=='CONTACT ANISOTROPIE LOCALE   :') then
    calcul_contact_anisotropies_locale=1
    open (unit=105,file='CONTACT_ANISOTROPIES_LOCALES.DAT',status='replace')     
    cycle
  end if

  if (command=='WRITE ALL GMV                :') then
    calcul_write_all_gmv=1
    read(1,*) step_all_gmv
    cycle
  end if

  if (command=='PDF                          :') then
    calcul_pdf=1
    read(1,*) step_pdf_start
    read(1,*) step_pdf_stop
    open (unit=106,file='PDF.DAT',status='replace')     
    open (unit=108,file='PDF_Ac.DAT',status='replace')     
    open (unit=109,file='PDF_Az2.DAT',status='replace')     
    open (unit=110,file='PDF_Az3.DAT',status='replace')     
    open (unit=111,file='PDF_Az4.DAT',status='replace')     
    open (unit=112,file='PDF_Az5.DAT',status='replace')     
    open (unit=113,file='PDF_Az6.DAT',status='replace')     
    open (unit=114,file='PDF_Az7.DAT',status='replace')     
    open (unit=115,file='PDF_Az8.DAT',status='replace')     
    open (unit=116,file='PDF_Az9.DAT',status='replace')     
    open (unit=117,file='PDF_Az10.DAT',status='replace')     
    open (unit=118,file='PDF_Az11.DAT',status='replace')     
    open (unit=119,file='PDF_Az12.DAT',status='replace')     
    cycle
  end if

  if (command=='KSI RESEAU                   :') then
    calcul_ksi_reseau=1
    read(1,*) step_ksi_reseau
    read(1,*) stop_ksi_reseau
    allocate( Mean_Ksi_reseau( int(stop_ksi_reseau-step_ksi_reseau)+1,10000,10 ) )
    open (unit=107,file='KSI_RESEAU.DAT',status='replace')     
    cycle
  end if


  if (command=='CONSTRUCT BODIES             :') then
    calcul_construct_bodies=1
    read(1,*) quel_pas_de_temps
    open (unit=128,file='NEW_BODIES.DAT',status='replace') 
    open (unit=129,file='NEW_DRV_DOF.DAT',status='replace') 
    open (unit=130,file='NEW_DOF.INI',status='replace') 
    cycle
  end if                                             


  if (command=='PROFILES                     :') then
    calcul_profils = 1
    the_first_time_profils = .true.
    read(1,*) debut_profils
    read(1,*) fin
    open (unit=131,file='PROFIL_VITESSE_MOYEN.DAT',status='replace')     
    open (unit=132,file='PROFIL_CONTRAINTE_MOYEN.DAT',status='replace')     
    open (unit=133,file='PROFIL_I_MOYEN.DAT',status='replace')     
    open (unit=134,file='PROFIL_COMPACITY_MOYEN.DAT',status='replace')     
    cycle
  end if

  if (command=='CONTACT PERSISTANCE          :') then
    calcul_persiste_contacts = 1
    read(1,*) debut_persistance
    open (unit=134,file='CONTACT_PERSISTANCE.DAT',status='replace')     
    cycle
  end if

!==============================================================================
!       Pour un seul pas de temps
!==============================================================================


  if (command=='CONTACT ORIENTATION          :') then
    calcul_contact_orientation=1
    open (unit=126,file='ORIENTATION_CONTACT_xy.DAT',status='replace') 
    open (unit=124,file='ORIENTATION_CONTACT_THETA.DAT',status='replace') 
    cycle
  end if                                             

  if (command=='CONTACT ORIENTATION 3D       :') then
    calcul_contact_orientation3D=1
    open (unit=127,file='ORIENTATION_CONTACT_3D.DAT',status='replace') 
    cycle
  end if                                             

  if (command=='SPATIAL CORRELATION          :') then
    calcul_spatial_correlation = 1
    open (unit=130,file='SPATIAL_CORRELATION.DAT',status='replace') 
    open (unit=131,file='SPATIAL_CORRELATION_PDF_RMAX.DAT',status='replace') 
    open (unit=132,file='SPATIAL_CORRELATION_RADIUS_MAX.DAT',status='replace') 
    open (unit=133,file='RADIAL_DISTRIBUTION_FLOAT.DAT',status='replace') 
    cycle
  end if   
                                            
  if (command=='RADIAL DISTRIBUTION          :') then
    calcul_radial_distribution = 1
    open (unit=133,file='RADIAL_DISTRIBUTION_FLOAT.DAT',status='replace') 
    cycle
  end if                                             
  
                                            
  if (command=='CORRELATION FUNCTION         :') then
    calcul_correlation_function = 1
    open (unit=133,file='CORRELATION_FUNCTION.DAT',status='replace') 
    cycle
  end if                                             


  if (command=='END                          :') exit
end do
close(1)

!==============================================================================
!Appel des differente ressources lues dans le fichier de commande
!==============================================================================

do 
call read_Vloc_dof_bodies

!-----Lecture de tout les pas --------------
if (all_pas>=0) then
  if (fin_post) call fermer
  call Calcul_deformation

  if (calcul_qsurp                          == 1)  then
    call qsurp    
!    call qsurp_moment
  end if
  if (calcul_coordination                   == 1)  call nb_coordination
  if (calcul_vitesse_moyenne                == 1)  call vitesse_moyenne
  if (calcul_compacity                      == 1)  call compacity
  if (calcul_contact_anisotropies           == 1)  then
    call contact_anisotropies
  end if
  if (calcul_contact_anisotropies_locale    == 1)  call contact_anisotropies_locales
  if (calcul_write_all_gmv                  == 1)  call write_all_gmv
  if (calcul_pdf                            == 1)  call pdf
  if (calcul_ksi_reseau                     == 1)  then
    pas_ksi = pas_ksi + 1
    call ksi_reseau 
  end if
  if (calcul_persiste_contacts              == 1)  call persiste_contact
  if (calcul_construct_bodies               == 1)  call construct_bodies
  if (calcul_profils                        == 1)  call profils
  if (calcul_contact_orientation            == 1)  call contact_orientation
  if (calcul_spatial_correlation            == 1)  call spatial_correlation
  if (calcul_radial_distribution            == 1)  call radial_distribution
  if (calcul_correlation_function           == 1)  call correlation_function
  

end if

!-----Lecture d'un seul pas --------------
!if (all_pas>=0) then
!  if (calcul_contact_orientation            == 1)  call contact_orientation
!  if (calcul_contact_orientation3D          == 1)  call contact_orientation_3D
!  call fermer
!end if 
end do


!==============================================================================
!==============================================================================
!==============================================================================
!==============================================================================

contains
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
    
    do i=1,nbParticules
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
    write(134,'(6(1X,D14.7))') temps,Haut,temps/(el_diameter * sqrt(1000/pressure_in_sample)),cpt/cpt_total
    print*, 0.05/Haut,cpt/cpt_total
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
real(kind=8)                            :: centre_band,xmax_box,xmin_box,ymax_box,ymin_box,zmax_box,zmin_box,&
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
  do i=1,nbParticules
    dmoyen = TAB_SPHERE(i)%Rmax + dmoyen
  end do
  dmoyen = 2*dmoyen / real(nbParticules,8)
  print*,'Mean diamater =',dmoyen
  ep = dmoyen

!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  if (the_first_time_profils) then
    nb_ligne_tab = 0
    do
      nb_ligne_tab  =   nb_ligne_tab + 1
      ep = ep + dmoyen/10
      if (ep>Haut+dmoyen) exit
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
    do i=1,nbParticules
      if (TAB_SPHERE(i)%behav/='PLEXx') cycle
      
      if ((TAB_SPHERE(i)%centre(3)<ep+dmoyen/10).and.(TAB_SPHERE(i)%centre(3)>ep)) then
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
      if ( abs(TAB_SPHERE(i)%centre(3) - ep)            <= TAB_SPHERE(i)%Rmax .and. &
           abs(TAB_SPHERE(i)%centre(3) - (ep+dmoyen/10))<= TAB_SPHERE(i)%Rmax ) then
           sphere_in_band = .true. 
           zmax_box = ep+dmoyen/10
           zmin_box = ep
      end if
      if ( abs(TAB_SPHERE(i)%centre(3) - ep)            >= TAB_SPHERE(i)%Rmax .and. &
           abs(TAB_SPHERE(i)%centre(3) - (ep+dmoyen/10))<= TAB_SPHERE(i)%Rmax ) then
           sphere_in_band = .true. 
           zmax_box = ep+dmoyen/10
           zmin_box = TAB_SPHERE(i)%centre(3) - TAB_SPHERE(i)%Rmax
      end if
           
      if ( abs(TAB_SPHERE(i)%centre(3) - ep)            <= TAB_SPHERE(i)%Rmax .and. &
           abs(TAB_SPHERE(i)%centre(3) - (ep+dmoyen/10))>= TAB_SPHERE(i)%Rmax ) then
           sphere_in_band = .true. 
           zmax_box = TAB_SPHERE(i)%centre(3) + TAB_SPHERE(i)%Rmax
           zmin_box = ep
      end if
      xmax_box = TAB_SPHERE(i)%centre(1)+TAB_SPHERE(i)%Rmax
      xmin_box = TAB_SPHERE(i)%centre(1)-TAB_SPHERE(i)%Rmax
      ymax_box = TAB_SPHERE(i)%centre(2)+TAB_SPHERE(i)%Rmax
      ymin_box = TAB_SPHERE(i)%centre(2)-TAB_SPHERE(i)%Rmax
      
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
            if ( (test_x-TAB_SPHERE(i)%centre(1))**2+&
                 (test_y-TAB_SPHERE(i)%centre(2))**2+&
                 (test_z-TAB_SPHERE(i)%centre(3))**2   <= TAB_SPHERE(i)%Rmax**2 )  nb_point_sphere = nb_point_sphere + 1.D0
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
        Lik(1:3) = TAB_SPHERE(cd)%centre(1:3)-TAB_SPHERE(an)%centre(1:3)
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
    
    vboite = Long*Larg*dmoyen/10

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
                           
    if (ep>Haut+dmoyen) exit
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
! Calcul du Ksi-reseau 
!==============================================================================
subroutine ksi_reseau
implicit none
real(kind=8),dimension(3,3)              :: Fabric,Fabric_LN,Fabric_LT,Fabric_L,Moment,M
real(kind=8),dimension(3,3)              :: Fabric_FN,Fabric_FT,Fabric_F
real(kind=8),dimension(3)                :: nik,tik,sik,Lik,Fik
real(kind=8)                             :: Rtik,Rnik,Rsik,norm_FN,FN,FT,LN,LT,Norm_F,Norm_L,cpt,epsilon
integer                                  :: i,j,cd,an,nb_ksi
real(kind=8)                             :: a,aln,alt,afn,aft,qp,ksi,PRESSION_FT,PRESSION_FN,PRESSION_P,PRESSION_L,&
                                            kmobilized,fmobilized,nb_ksi_ctc, S1,S2,S3,Rayon_max,Theta_fabric,theta_sigma,&
                                            fmax,fmin
real(kind=8),dimension(3)                :: wr,wi
real(kind=8),dimension(3,3)              :: localframe  
integer                                  :: ierror,matz,lda
real(kind=8),dimension(:,:),allocatable  :: mean_val


if ( (compteur_clout >= step_ksi_reseau).and.(compteur_clout <= stop_ksi_reseau) ) then
  print*, '---> Ksi-reseau'
  Lik       = 0.D0
  Fik       = 0.D0
  Norm_L    = 0.D0
  Norm_F    = 0.D0
  norm_FN   = 0.D0
  cpt       = 0.D0
  Fabric    = 0.D0
  Fabric_LN = 0.D0
  Fabric_LT = 0.D0
  Fabric_FN = 0.D0
  Fabric_FT = 0.D0
  Moment    = 0.D0
  fmax      = -100000000000.D0
  fmin      = 1000000000000.D0
  do i=1,nb_ligneCONTACT
    cd   = TAB_CONTACT(i)%icdent
    an   = TAB_CONTACT(i)%ianent
    if (TAB_CONTACT(i)%rn<0.000000000001*Mean_total_Normal_force) cycle
    if ((TAB_SPHERE(cd)%behav/='PLEXx').or.(TAB_SPHERE(an)%behav/='PLEXx')) cycle

    nik  = TAB_CONTACT(i)%n
    tik  = TAB_CONTACT(i)%t
    sik  = TAB_CONTACT(i)%s
    Rtik = TAB_CONTACT(i)%rt
    Rnik = TAB_CONTACT(i)%rn
    Rsik = TAB_CONTACT(i)%rs

    Lik(1:3) = TAB_SPHERE(cd)%centre(1:3)-TAB_SPHERE(an)%centre(1:3)
    if ( sqrt( Lik(1)**2 + Lik(2)**2 + Lik(3)**2) > 2*( TAB_SPHERE(cd)%Rmax + TAB_SPHERE(an)%Rmax ) ) then
      Rayon_max = sqrt( ( Lik(1)**2 + Lik(2)**2 + Lik(3)**2) ) 
      lik(:) = Lik(:) / Rayon_max
      Lik(:) = abs(TAB_SPHERE(cd)%Rmax+TAB_SPHERE(an)%Rmax)*Lik(:)
    end if

    norm_FN  = norm_FN + Rnik
    
    Fik(1:3) = (Rnik*nik(1:3)+Rtik*tik(1:3)+Rsik*sik(1:3))
    Norm_F   = Norm_F + sqrt( Fik(1)**2 + Fik(2)**2 + Fik(3)**2 )
    Norm_L   = Norm_L + sqrt( Lik(1)**2 + Lik(2)**2 + Lik(3)**2 )

    LN       =  sqrt( Lik(1)**2 + Lik(2)**2 + Lik(3)**2 ) 
    FN       =  Rnik
   
    tik(1:3) =  Fik(1:3) - FN*nik(1:3)
    FT       =  sqrt( tik(1)**2 + tik(2)**2 + tik(3)**2 )

    if ( FT  > 0.00000000001) tik(1:3) = tik(1:3) / FT
    if ( FT <= 0.00000000001) tik(1:3) = 0.D0

    FT       =  Fik(1)*tik(1) + Fik(2)*tik(2) + Fik(3)*tik(3)  

    Moment(1,1:3) = Fik(1)*Lik(1:3) + Moment(1,1:3)
    Moment(2,1:3) = Fik(2)*Lik(1:3) + Moment(2,1:3)
    Moment(3,1:3) = Fik(3)*Lik(1:3) + Moment(3,1:3)

    Fabric(1,1:3) = nik(1)*nik(1:3) + Fabric(1,1:3)
    Fabric(2,1:3) = nik(2)*nik(1:3) + Fabric(2,1:3)
    Fabric(3,1:3) = nik(3)*nik(1:3) + Fabric(3,1:3)

    Fabric_LN(1,1:3) = LN*nik(1)*nik(1:3) + Fabric_LN(1,1:3)
    Fabric_LN(2,1:3) = LN*nik(2)*nik(1:3) + Fabric_LN(2,1:3)
    Fabric_LN(3,1:3) = LN*nik(3)*nik(1:3) + Fabric_LN(3,1:3)

    Fabric_FN(1,1:3) = FN*nik(1)*nik(1:3) + Fabric_FN(1,1:3)
    Fabric_FN(2,1:3) = FN*nik(2)*nik(1:3) + Fabric_FN(2,1:3)
    Fabric_FN(3,1:3) = FN*nik(3)*nik(1:3) + Fabric_FN(3,1:3)
  
    Fabric_FT(1,1:3) = FT*nik(1)*tik(1:3) + Fabric_FT(1,1:3)
    Fabric_FT(2,1:3) = FT*nik(2)*tik(1:3) + Fabric_FT(2,1:3)
    Fabric_FT(3,1:3) = FT*nik(3)*tik(1:3) + Fabric_FT(3,1:3)

    cpt = cpt + 1.D0
    fmax = max(fmax,FN)
    fmin = min(fmin,FN)
  end do
  norm_FN     = norm_FN / cpt
  norm_F      = norm_F  / cpt
  norm_L      = norm_L  / cpt
  
  Fabric      = Fabric(:,:) / cpt
  Fabric_LN   = Fabric_LN   / (Norm_L * cpt)
  Fabric_FN   = Fabric_FN    / (Norm_F * cpt)
  Fabric_FT   = Fabric_FT    / (Norm_F * cpt)
  Fabric_F    = Fabric_FT + Fabric_FN

  PRESSION_P  = (Moment(1,1)+Moment(2,2)+Moment(3,3)) / 3
  PRESSION_L  = Fabric_LN(1,1)+Fabric_LN(2,2)+Fabric_LN(3,3)
  PRESSION_FN = Fabric_FN(1,1)+Fabric_FN(2,2)+Fabric_FN(3,3)
  PRESSION_FT = Fabric_F(1,1)+Fabric_F(2,2)+Fabric_F(3,3)

  theta_sigma = 3.141592/4.0
!-------------------------------------------------------------------!
  ksi = 0.D0
  nb_ksi = 0
  epsilon = 0.01
  do 
    ksi = ksi + epsilon
    nb_ksi = nb_ksi + 1    
    Lik       = 0.D0
    Fik       = 0.D0
    Fabric    = 0.D0
    Fabric_LN = 0.D0
    Fabric_LT = 0.D0
    Fabric_FN = 0.D0
    Fabric_FT = 0.D0
    Moment    = 0.D0
    kmobilized= 0.D0
    fmobilized= 0.D0
    nb_ksi_ctc= 0.D0

    do i=1,nb_ligneCONTACT
      cd   = TAB_CONTACT(i)%icdent
      an   = TAB_CONTACT(i)%ianent
      if (TAB_CONTACT(i)%rn<0.000000000001*Mean_total_Normal_force) cycle
      if ((TAB_SPHERE(cd)%behav/='PLEXx').or.(TAB_SPHERE(an)%behav/='PLEXx')) cycle
      if (TAB_CONTACT(i)%rn<0.000000000001*Mean_total_Normal_force) cycle
      
!      if  ( ((TAB_CONTACT(i)%rn/norm_FN) > ksi) .and.(TAB_CONTACT(i)%rn/norm_FN) < ksi+epsilon) then
      if  ( ( (TAB_CONTACT(i)%rn/norm_FN) < ksi+epsilon) ) then
        cd   = TAB_CONTACT(i)%icdent
        an   = TAB_CONTACT(i)%ianent
        nik  = TAB_CONTACT(i)%n
        tik  = TAB_CONTACT(i)%t
        sik  = TAB_CONTACT(i)%s
        Rtik = TAB_CONTACT(i)%rt
        Rnik = TAB_CONTACT(i)%rn
        Rsik = TAB_CONTACT(i)%rs

        Lik(1:3) = TAB_SPHERE(cd)%centre(1:3)-TAB_SPHERE(an)%centre(1:3)
        if ( sqrt( Lik(1)**2 + Lik(2)**2 + Lik(3)**2) > 2*( TAB_SPHERE(cd)%Rmax + TAB_SPHERE(an)%Rmax ) ) then
          Rayon_max = sqrt( ( Lik(1)**2 + Lik(2)**2 + Lik(3)**2) ) 
          lik(:) = Lik(:) / Rayon_max
          Lik(:) = abs(TAB_SPHERE(cd)%Rmax+TAB_SPHERE(an)%Rmax)*Lik(:)
        end if

        Fik(1:3) = (Rnik*nik(1:3)+Rtik*tik(1:3)+Rsik*sik(1:3))
  
        LN       =  sqrt( Lik(1)**2 + Lik(2)**2 + Lik(3)**2 ) 
        FN       =  Rnik
   
        tik(1:3) =  Fik(1:3) - FN*nik(1:3)
        FT       =  sqrt( tik(1)**2 + tik(2)**2 + tik(3)**2 )

        if ( FT  > 0.00000000001) tik(1:3) = tik(1:3) / FT
        if ( FT <= 0.00000000001) tik(1:3) = 0.D0

        FT       =  Fik(1)*tik(1) + Fik(2)*tik(2) + Fik(3)*tik(3)  

        Moment(1,1:3) = Fik(1)*Lik(1:3) + Moment(1,1:3)
        Moment(2,1:3) = Fik(2)*Lik(1:3) + Moment(2,1:3)
        Moment(3,1:3) = Fik(3)*Lik(1:3) + Moment(3,1:3)

        Fabric(1,1:3) = nik(1)*nik(1:3) + Fabric(1,1:3)
        Fabric(2,1:3) = nik(2)*nik(1:3) + Fabric(2,1:3)
        Fabric(3,1:3) = nik(3)*nik(1:3) + Fabric(3,1:3)

        Fabric_LN(1,1:3) = LN*nik(1)*nik(1:3) + Fabric_LN(1,1:3)
        Fabric_LN(2,1:3) = LN*nik(2)*nik(1:3) + Fabric_LN(2,1:3)
        Fabric_LN(3,1:3) = LN*nik(3)*nik(1:3) + Fabric_LN(3,1:3)

        Fabric_FN(1,1:3) = FN*nik(1)*nik(1:3) + Fabric_FN(1,1:3)
        Fabric_FN(2,1:3) = FN*nik(2)*nik(1:3) + Fabric_FN(2,1:3)
        Fabric_FN(3,1:3) = FN*nik(3)*nik(1:3) + Fabric_FN(3,1:3)
  
        Fabric_FT(1,1:3) = FT*nik(1)*tik(1:3) + Fabric_FT(1,1:3)
        Fabric_FT(2,1:3) = FT*nik(2)*tik(1:3) + Fabric_FT(2,1:3)
        Fabric_FT(3,1:3) = FT*nik(3)*tik(1:3) + Fabric_FT(3,1:3)

        if ( abs( sqrt(Rtik**2+Rsik**2) ) > 0.39 * Rnik ) &
                                   kmobilized = kmobilized+1.D0
                                   
        fmobilized = fmobilized + &
                      abs( sqrt(Rtik**2+Rsik**2) ) / ( 0.4 * Rnik )

        nb_ksi_ctc =  nb_ksi_ctc + 1.D0
      end if
    end do
    Fabric      = Fabric(:,:) / cpt
    Fabric_LN   = Fabric_LN   / (Norm_L * cpt)
    Fabric_FN   = Fabric_FN   / (Norm_F * cpt)
    Fabric_FT   = Fabric_FT   / (Norm_F * cpt)
    Fabric_F    = Fabric_FT + Fabric_FN

    a = 0.D0 ; aln = 0.D0 ; alt = 0.D0 ; afn = 0.D0 ; aft = 0.D0

    lda  = 3
    matz = 1

    M = Moment
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
    qp = ((S3-S1)/3) / PRESSION_P

    M = Fabric
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

    if (S3==wr(1)) Theta_fabric = asin(abs(localframe(3,1))/sqrt(localframe(1,1)**2+localframe(2,1)**2+localframe(3,1)**2) )
    if (S3==wr(2)) Theta_fabric = asin(abs(localframe(3,2))/sqrt(localframe(1,2)**2+localframe(2,2)**2+localframe(3,2)**2) )
    if (S3==wr(3)) Theta_fabric = asin(abs(localframe(3,3))/sqrt(localframe(1,3)**2+localframe(2,3)**2+localframe(3,3)**2) )
    a = 2*(S3-S1)*cos(2*abs(Theta_fabric-theta_sigma))

    aln = 0.D0

    M = Fabric_FN
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
    afn = 2*(S3-S1)/PRESSION_FN - a

    M = Fabric_F
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
    aft = 2*(S3-S1)/PRESSION_FT - a - afn
    
    if (nb_ksi_ctc==0) then
      kmobilized = 0
      fmobilized = 0
    else  
      kmobilized = kmobilized / nb_ksi_ctc
      fmobilized = fmobilized / nb_ksi_ctc
    end if

    Mean_Ksi_reseau(pas_ksi,nb_ksi,1) = ksi
    Mean_Ksi_reseau(pas_ksi,nb_ksi,2) = qp
    Mean_Ksi_reseau(pas_ksi,nb_ksi,3) = a
    Mean_Ksi_reseau(pas_ksi,nb_ksi,4) = aln
    Mean_Ksi_reseau(pas_ksi,nb_ksi,5) = afn
    Mean_Ksi_reseau(pas_ksi,nb_ksi,6) = aft
    Mean_Ksi_reseau(pas_ksi,nb_ksi,7) = (a+aln+afn+aft)/2
    Mean_Ksi_reseau(pas_ksi,nb_ksi,8) = kmobilized
    Mean_Ksi_reseau(pas_ksi,nb_ksi,9) = fmobilized
    Mean_Ksi_reseau(pas_ksi,nb_ksi,10) = nb_ksi_ctc / cpt

    if (ksi>50) exit

  end do
!  all_pas = all_pas + 1

  if ( (compteur_clout == stop_ksi_reseau ) ) then
    allocate(mean_val(10,nb_ksi))
    mean_val(:,:) = 0
  
    do i=1,pas_ksi
      do j=1,nb_ksi
        mean_val(1,j) = mean_val(1,j) + Mean_Ksi_reseau(i,j,1)
        mean_val(2,j) = mean_val(2,j) + Mean_Ksi_reseau(i,j,2)
        mean_val(3,j) = mean_val(3,j) + Mean_Ksi_reseau(i,j,3)
        mean_val(4,j) = mean_val(4,j) + Mean_Ksi_reseau(i,j,4)
        mean_val(5,j) = mean_val(5,j) + Mean_Ksi_reseau(i,j,5)
        mean_val(6,j) = mean_val(6,j) + Mean_Ksi_reseau(i,j,6)
        mean_val(7,j) = mean_val(7,j) + Mean_Ksi_reseau(i,j,7)
        mean_val(8,j) = mean_val(8,j) + Mean_Ksi_reseau(i,j,8)
        mean_val(9,j) = mean_val(9,j) + Mean_Ksi_reseau(i,j,9)
        mean_val(10,j) = mean_val(10,j) + Mean_Ksi_reseau(i,j,10)
      end do
    end do
    mean_val(:,:) = mean_val(:,:) / pas_ksi


    i = 0
    do 
      i = i + 1
      if (mean_val(1,i)>50) exit
      if (mean_val(10,i)==0) cycle
      write(107,'(16(1X,D12.5))') mean_val(1,i), mean_val(2,i), mean_val(3,i),&
                                                         mean_val(4,i), mean_val(5,i),&
                                                         mean_val(6,i), mean_val(7,i),&
                                                         mean_val(8,i), mean_val(9,i),&
                                                         mean_val(10,i)
    end do
  
    stop
  end if
end if 


end subroutine ksi_reseau


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

do i=1,nbParticules
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
  all_pas = all_pas + 1
  
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
  
    Lik(1:3) = TAB_SPHERE(cd)%centre(1:3)-TAB_SPHERE(an)%centre(1:3)
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

  do i=1,nbParticules
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
! Calcul des differentes anisotropies de contact
!==============================================================================
subroutine contact_anisotropies
implicit none
real(kind=8),dimension(3,3)              :: Fabric,Fabric_LN,Fabric_LT,Fabric_L,Matrice_temp,Moment
real(kind=8),dimension(3,3)              :: Fabric_FN,Fabric_FT,Fabric_F
real(kind=8),dimension(3)                :: nik,tik,sik,Lik,Fik
real(kind=8)                             :: Rtik,Rnik,Rsik,norm,FN,FT,LN,LT,Norm_F,Norm_L,cpt,Rayon_max
integer                                  :: i,cd,an,j
real(kind=8)                             :: ac,aln,alt,afn,aft,S1,S2,S3
real(kind=8)                             :: acnd,alnnd,altnd,afnnd,aftnd
real(kind=8)                             :: acp,alnp,altp,afnp,aftp
real(kind=8)                             :: afnndp,aftndp
real(kind=8),dimension(3)                :: wr,wi
real(kind=8),dimension(3,3)              :: localframe  
integer                                  :: ierror,matz,lda

real(kind=8),dimension(3,3)              :: InvFabric,Fabric_psiN,Fabric_psiL,InvFabric_psiN,Fabric_psiT
real(kind=8),dimension(3,3)              :: Matrix_rotation_sigma,Matrix_rotation_contact
real(kind=8),dimension(3,3)              :: Matrix_A,invMatrix_A,Matrix_rotation_force_normale,Matrix_rotation_A,&
                                            Matrix_rotation_force_tangente
real(kind=8),dimension(3,3)              :: Fabric_psiNv2,Fabric_psiNv2_temp,Fabric_psiTv2,Fabric_psiTv2_temp
real(kind=8),dimension(3,3)              :: PsiN,PsiT
real(kind=8)                             :: determinant,theta_c,theta_l,theta_n,theta_t,theta_sigma,n_theta,&
                                            theta_epsilon,theta_normale,Theta_A,Trace_FN,Trace_F,Trace_L

print*,'.... contact anisotropies ....'

Lik       = 0.D0
Fik       = 0.D0
Norm_L    = 0.D0
Norm_F    = 0.D0
cpt       = 0.D0
Fabric    = 0.D0
Fabric_LN = 0.D0
Fabric_LT = 0.D0
Fabric_FN = 0.D0
Fabric_FT = 0.D0
Moment    = 0.D0                            



do i=1,nb_ligneCONTACT
  cd   = TAB_CONTACT(i)%icdent
  an   = TAB_CONTACT(i)%ianent
  if (TAB_CONTACT(i)%rn<0.000000000001*Mean_total_Normal_force) cycle
  if ((TAB_SPHERE(cd)%behav/='PLEXx').or.(TAB_SPHERE(an)%behav/='PLEXx')) cycle

  nik  = TAB_CONTACT(i)%n
  tik  = TAB_CONTACT(i)%t
  sik  = TAB_CONTACT(i)%s
  
  Lik(1:3) = TAB_SPHERE(cd)%centre(1:3)-TAB_SPHERE(an)%centre(1:3)
  if ( sqrt( Lik(1)**2 + Lik(2)**2 + Lik(3)**2) > 2*( TAB_SPHERE(cd)%Rmax + TAB_SPHERE(an)%Rmax ) ) then
    Rayon_max = sqrt( ( Lik(1)**2 + Lik(2)**2 + Lik(3)**2) ) 
    lik(:) = Lik(:) / Rayon_max
    Lik(:) = abs(TAB_SPHERE(cd)%Rmax+TAB_SPHERE(an)%Rmax)*Lik(:)
  end if

  Rtik = TAB_CONTACT(i)%rt
  Rnik = TAB_CONTACT(i)%rn
  Rsik = TAB_CONTACT(i)%rs

  Fik(1:3) = (Rnik*nik(1:3)+Rtik*tik(1:3)+Rsik*sik(1:3))
  Norm_F   = Norm_F + Rnik !Norm_F + sqrt( Fik(1)**2 + Fik(2)**2 + Fik(3)**2 )
  Norm_L   = Norm_L + sqrt( Lik(1)**2 + Lik(2)**2 + Lik(3)**2 )

  LN       =  sqrt( Lik(1)**2 + Lik(2)**2 + Lik(3)**2 ) 
  FN       =  Rnik
   
  tik(1:3) =  Fik(1:3) - FN*nik(1:3)
  FT       =  sqrt( tik(1)**2 + tik(2)**2 + tik(3)**2 )


  if ( FT  > 0.00001*FN) tik(1:3) = tik(1:3) / FT
  if ( FT <= 0.00001*FN) tik(1:3) = 0.D0

  FT       =  Fik(1)*tik(1) + Fik(2)*tik(2) + Fik(3)*tik(3)  

  Fabric(1,1:3) = nik(1)*nik(1:3) + Fabric(1,1:3)
  Fabric(2,1:3) = nik(2)*nik(1:3) + Fabric(2,1:3)
  Fabric(3,1:3) = nik(3)*nik(1:3) + Fabric(3,1:3)

  Fabric_LN(1,1:3) = LN*nik(1)*nik(1:3) + Fabric_LN(1,1:3)
  Fabric_LN(2,1:3) = LN*nik(2)*nik(1:3) + Fabric_LN(2,1:3)
  Fabric_LN(3,1:3) = LN*nik(3)*nik(1:3) + Fabric_LN(3,1:3)

  Fabric_FN(1,1:3) = FN*nik(1)*nik(1:3) + Fabric_FN(1,1:3)
  Fabric_FN(2,1:3) = FN*nik(2)*nik(1:3) + Fabric_FN(2,1:3)
  Fabric_FN(3,1:3) = FN*nik(3)*nik(1:3) + Fabric_FN(3,1:3)
  
  Fabric_FT(1,1:3) = FT*nik(1)*tik(1:3) + Fabric_FT(1,1:3)
  Fabric_FT(2,1:3) = FT*nik(2)*tik(1:3) + Fabric_FT(2,1:3)
  Fabric_FT(3,1:3) = FT*nik(3)*tik(1:3) + Fabric_FT(3,1:3)
                
  Moment(1,1:3) = Fik(1)*Lik(1:3) + Moment(1,1:3)
  Moment(2,1:3) = Fik(2)*Lik(1:3) + Moment(2,1:3)
  Moment(3,1:3) = Fik(3)*Lik(1:3) + Moment(3,1:3)

  cpt = cpt + 1
end do

Norm_F = Norm_F / real(cpt,8)
Norm_L = Norm_L / real(cpt,8)

Fabric      = Fabric(:,:) / (cpt)
Fabric_LN   = Fabric_LN   / (Norm_L*cpt)
Fabric_FN   = Fabric_FN   / (Norm_F*cpt)
Fabric_FT   = Fabric_FT   / (Norm_F*cpt)
Fabric_F    = Fabric_FT + Fabric_FN

Trace_FN = Fabric_FN(1,1)+Fabric_FN(2,2)+Fabric_FN(3,3)
Trace_F  = Fabric_F(1,1)+Fabric_F(2,2)+Fabric_F(3,3)
Trace_L  = Fabric_LN(1,1)+Fabric_LN(2,2)+Fabric_LN(3,3)

!!!@@@@@ RECHERCHE DE THETA_SIGMA
lda  = 3
matz = 1

call rg ( lda, 3, Moment, wr, wi, matz, localframe, ierror )
if ( wr(1)==wr(2) .and. (wr(2)==wr(3)) ) then
  S1=wr(1)
  S2=S1
  S3=S1
else
  S3 = max( wr(1),max(wr(2),wr(3)))
  S1 = min( wr(1),min(wr(2),wr(3)))
  if (wr(1)==wr(2)) S2=S1
  if (wr(3)==wr(2)) S2=S3

  do i=1,3
    if ((wr(i)<S3) .and. (wr(i)>S1)) S2=wr(i)
  end do
end if

if (S3==wr(1)) Theta_sigma = asin(abs(localframe(3,1))/sqrt(localframe(1,1)**2+localframe(2,1)**2+localframe(3,1)**2) )
if (S3==wr(2)) Theta_sigma = asin(abs(localframe(3,2))/sqrt(localframe(1,2)**2+localframe(2,2)**2+localframe(3,2)**2) )
if (S3==wr(3)) Theta_sigma = asin(abs(localframe(3,3))/sqrt(localframe(1,3)**2+localframe(2,3)**2+localframe(3,3)**2) )


!!!@@@@ ANISOTROPIES CALCULEES AVEC LES TENSEURS XI
ac = 0.D0 ; aln = 0.D0 ; alt = 0.D0 ; afn = 0.D0 ; aft = 0.D0
acnd = 0.D0 ; alnnd = 0.D0 ; altnd = 0.D0 ; afnnd = 0.D0 ; aftnd = 0.D0

Matrice_temp = Fabric
call rg ( lda, 3, Matrice_temp, wr, wi, matz, localframe, ierror )
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
ac = 2*(S3-S1)/(S1+S3)

aln = 0

Matrice_temp = Fabric_FN
call rg ( lda, 3, Matrice_temp, wr, wi, matz, localframe, ierror )
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
afn = 2*( S3-S1 ) / (S1+S3)  - ac


Matrice_temp = Fabric_F
call rg ( lda, 3, Matrice_temp, wr, wi, matz, localframe, ierror )
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
aft = 2*( S3-S1 ) / (S1+S3) - ac - afn


if (type_box==0)&
write(104,'(31(1X,D12.5))') temps,Def,Def_H,&
                            ac,aln,afn,aft,(ac+aln+afn+aft)/2


end subroutine contact_anisotropies



!==============================================================================
! Calcul des differentes anisotropies de contact en derivant sur le vecteur force
!==============================================================================
subroutine contact_anisotropies_forces
implicit none
real(kind=8),dimension(3,3)              :: Fabric,Fabric_LN,Fabric_LT,Fabric_L,Matrice_temp,Moment
real(kind=8),dimension(3,3)              :: Fabric_F,Fabric_FL
real(kind=8),dimension(3)                :: nik,tik,sik,Lik,Fik
real(kind=8)                             :: Rtik,Rnik,Rsik,norm,LN,LT,Norm_F,Norm_L,cpt,Rayon_max,&
                                            Fnorm,lnorm,Trace_f,Trace_LN,Trace_L,norm_FL,Trace_FL
integer                                  :: i,cd,an,j
real(kind=8)                             :: ac,aln,alt,af,S1,S2,S3,afl
real(kind=8),dimension(3)                :: wr,wi
real(kind=8),dimension(3,3)              :: localframe  
integer                                  :: ierror,matz,lda

print*,'.... contact anisotropies ....'

Lik       = 0.D0
Fik       = 0.D0
Norm_L    = 0.D0
Norm_F    = 0.D0
norm_FL   = 0.D0
cpt       = 0.D0
Fabric    = 0.D0
Fabric_LN = 0.D0
Fabric_LT = 0.D0
Fabric_F  = 0.D0
Moment    = 0.D0                            
Fnorm     = 0.D0
Fabric_FL = 0.D0

do i=1,nb_ligneCONTACT
!  if (TAB_CONTACT(i)%rn==0.D0) cycle
  if (TAB_CONTACT(i)%rn<0.000000000001*Mean_total_Normal_force) cycle
  cd   = TAB_CONTACT(i)%icdent
  an   = TAB_CONTACT(i)%ianent
  nik  = TAB_CONTACT(i)%n
  tik  = TAB_CONTACT(i)%t
  sik  = TAB_CONTACT(i)%s
  
  Lik(1:3) = TAB_SPHERE(cd)%centre(1:3)-TAB_SPHERE(an)%centre(1:3)
  if ( sqrt( Lik(1)**2 + Lik(2)**2 + Lik(3)**2) > 2*( TAB_SPHERE(cd)%Rmax + TAB_SPHERE(an)%Rmax ) ) then
    Rayon_max = sqrt( ( Lik(1)**2 + Lik(2)**2 + Lik(3)**2) ) 
    lik(:) = Lik(:) / Rayon_max
    Lik(:) = abs(TAB_SPHERE(cd)%Rmax+TAB_SPHERE(an)%Rmax)*Lik(:)
  end if

  Rtik = TAB_CONTACT(i)%rt
  Rnik = TAB_CONTACT(i)%rn
  Rsik = TAB_CONTACT(i)%rs

  Fik(1:3) = (Rnik*nik(1:3)+Rtik*tik(1:3)+Rsik*sik(1:3))
  
  Fnorm    = sqrt( abs(Fik(1)**2 + Fik(2)**2 + Fik(3)**2) )
  lnorm    = sqrt( abs(Lik(1)**2 + Lik(2)**2 + Lik(3)**2) )
  
  Norm_F   = Norm_F  + sqrt( Fik(1)**2 + Fik(2)**2 + Fik(3)**2 )
  Norm_L   = Norm_L  + sqrt( Lik(1)**2 + Lik(2)**2 + Lik(3)**2 )
  norm_FL  = norm_FL + sqrt( abs( Fik(1)*lik(1) + Fik(2)*lik(2) + Fik(3)*lik(3) ) )
  
  nik(1:3) =  Fik(1:3) / Fnorm
  tik(1:3) =  Lik(1:3) - Fnorm*nik(1:3)
  
  LN       =  sqrt( abs(Lik(1)*nik(1) + Lik(2)*nik(2) + Lik(3)*nik(3)) )
  LT       =  sqrt( tik(1)**2 + tik(2)**2 + tik(3)**2 )
  
  if ( LT  > 0.00001*LN) tik(1:3) = tik(1:3) / LT
  if ( LT <= 0.00001*LN) tik(1:3) = 0.D0

  Fabric(1,1:3) = nik(1)*nik(1:3) + Fabric(1,1:3)
  Fabric(2,1:3) = nik(2)*nik(1:3) + Fabric(2,1:3)
  Fabric(3,1:3) = nik(3)*nik(1:3) + Fabric(3,1:3)

  Fabric_FL(1,1:3) = Fnorm*lnorm*nik(1)*nik(1:3) + Fabric_FL(1,1:3)
  Fabric_FL(2,1:3) = Fnorm*lnorm*nik(2)*nik(1:3) + Fabric_FL(2,1:3)
  Fabric_FL(3,1:3) = Fnorm*lnorm*nik(3)*nik(1:3) + Fabric_FL(3,1:3)

  Fabric_LN(1,1:3) = LN*nik(1)*nik(1:3) + Fabric_LN(1,1:3)
  Fabric_LN(2,1:3) = LN*nik(2)*nik(1:3) + Fabric_LN(2,1:3)
  Fabric_LN(3,1:3) = LN*nik(3)*nik(1:3) + Fabric_LN(3,1:3)

  Fabric_LT(1,1:3) = LT*nik(1)*tik(1:3) + Fabric_LT(1,1:3)
  Fabric_LT(2,1:3) = LT*nik(2)*tik(1:3) + Fabric_LT(2,1:3)
  Fabric_LT(3,1:3) = LT*nik(3)*tik(1:3) + Fabric_LT(3,1:3)

  Fabric_F(1,1:3) = Fnorm*nik(1)*nik(1:3) + Fabric_F(1,1:3)
  Fabric_F(2,1:3) = Fnorm*nik(2)*nik(1:3) + Fabric_F(2,1:3)
  Fabric_F(3,1:3) = Fnorm*nik(3)*nik(1:3) + Fabric_F(3,1:3)
                  
  Moment(1,1:3) = Fik(1)*Lik(1:3) + Moment(1,1:3)
  Moment(2,1:3) = Fik(2)*Lik(1:3) + Moment(2,1:3)
  Moment(3,1:3) = Fik(3)*Lik(1:3) + Moment(3,1:3)

  cpt = cpt + 1
end do

Norm_F  = Norm_F  / real(cpt,8)
Norm_L  = Norm_L  / real(cpt,8)
Norm_FL = Norm_FL / real(cpt,8)

Fabric      = Fabric(:,:) / (cpt)
Fabric_FL   = Fabric_FL(:,:) / (Norm_FL*cpt)
Fabric_LN   = Fabric_LN   / (Norm_L*cpt)
Fabric_LT   = Fabric_LT   / (Norm_L*cpt)
Fabric_L    = Fabric_LT   + Fabric_LN
Fabric_F    = Fabric_F    / (Norm_F*cpt)

Trace_F    = Fabric_F(1,1)+Fabric_F(2,2)+Fabric_F(3,3)
Trace_LN   = Fabric_LN(1,1)+Fabric_LN(2,2)+Fabric_LN(3,3)
Trace_L    = Fabric_L(1,1)+Fabric_L(2,2)+Fabric_L(3,3)
Trace_FL   = Fabric_FL(1,1)+Fabric_FL(2,2)+Fabric_FL(3,3)


!!!@@@@@ RECHERCHE DE THETA_SIGMA
lda  = 3
matz = 1

call rg ( lda, 3, Moment, wr, wi, matz, localframe, ierror )
if ( wr(1)==wr(2) .and. (wr(2)==wr(3)) ) then
  S1=wr(1)
  S2=S1
  S3=S1
else
  S3 = max( wr(1),max(wr(2),wr(3)))
  S1 = min( wr(1),min(wr(2),wr(3)))
  if (wr(1)==wr(2)) S2=S1
  if (wr(3)==wr(2)) S2=S3

  do i=1,3
    if ((wr(i)<S3) .and. (wr(i)>S1)) S2=wr(i)
  end do
end if

!if (S3==wr(1)) Theta_sigma = asin(abs(localframe(3,1))/sqrt(localframe(1,1)**2+localframe(2,1)**2+localframe(3,1)**2) )
!if (S3==wr(2)) Theta_sigma = asin(abs(localframe(3,2))/sqrt(localframe(1,2)**2+localframe(2,2)**2+localframe(3,2)**2) )
!if (S3==wr(3)) Theta_sigma = asin(abs(localframe(3,3))/sqrt(localframe(1,3)**2+localframe(2,3)**2+localframe(3,3)**2) )
                            
!!!@@@@ ANISOTROPIES CALCULEES AVEC LES TENSEURS PSI
ac = 0.D0 ; aln = 0.D0 ; alt = 0.D0 ; af  = 0.D0 ; afl = 0.D0

Matrice_temp = Fabric
call rg ( lda, 3, Matrice_temp, wr, wi, matz, localframe, ierror )
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
ac = 4*(S3-S1)
!print*, S1+S2+S3

Matrice_temp = Fabric_FL
call rg ( lda, 3, Matrice_temp, wr, wi, matz, localframe, ierror )
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
!print*, S1+S2+S3
afl  = 4*( S3-S1 ) / Trace_FL  - ac


Matrice_temp = Fabric_LN
call rg ( lda, 3, Matrice_temp, wr, wi, matz, localframe, ierror )
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
!print*, S1+S2+S3
aln  = 4*( S3-S1 ) / Trace_LN  - ac

Matrice_temp = Fabric_L
call rg ( lda, 3, Matrice_temp, wr, wi, matz, localframe, ierror )
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
!print*, S1+S2+S3
alt  = 4*( S3-S1 ) / Trace_L  - ac -aln

Matrice_temp = Fabric_F
call rg ( lda, 3, Matrice_temp, wr, wi, matz, localframe, ierror )
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
!print*, S1+S2+S3

af  = 4*( S3-S1 ) / Trace_F  - ac

!print*,ac,aln,af,afl
!print*,(ac+af+aln+alt)/4
!print*, (ac+afl)/4

if (type_box==0)&
write(1005,'(31(1X,D12.5))') temps,Def,Def_H,&
                            ac,aln,alt,af,(ac+af+aln+alt)/4,&
                            afl, (ac+afl)/4


                            
end subroutine contact_anisotropies_forces


!==============================================================================
! Calcul des differentes anisotropies de contact locales
!==============================================================================
subroutine contact_anisotropies_locales
implicit none
real(kind=8),dimension(3,3)              :: Fabric,Fabric_LN,Fabric_LT,Fabric_L,localframe
real(kind=8),dimension(3,3)              :: Fabric_FN,Fabric_FT,Fabric_F
real(kind=8),dimension(3)                :: nik,tik,sik,Lik,Fik
real(kind=8)                             :: Rtik,Rnik,Rsik,norm,FN,FT,LN,LT,Norm_F,Norm_L,cpt,&
                                            ac,acS,Sa,Sa_NS,Theta_c,Theta_fabric,&
                                            Theta_c_s,Theta_fabric_s,S1,S2,S3
integer                                  :: i,cd,an,j
real(kind=8),dimension(3,3)              :: a,aln,alt,afn,aft
real(kind=8),dimension(3)                :: wr,wi
integer                                  :: ierror,matz,lda


print*,'.... contact anisotropies locales ....'


do i=1,nbParticules
  TAB_SPHERE(i)%Fabric = 0.D0
  TAB_SPHERE(i)%Fabric_N = 0.D0
  TAB_SPHERE(i)%Fabric_T = 0.D0
  TAB_SPHERE(i)%nb_ctc = 0
  TAB_SPHERE(i)%Norm_F = 0.D0
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
  
  Lik(1:3) = TAB_SPHERE(cd)%centre(1:3)-TAB_SPHERE(an)%centre(1:3)
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

lda  = 3
matz = 1
call rg ( lda, 3, Fabric, wr, wi, matz, localframe, ierror )
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

if (S3==wr(1)) Theta_fabric = acos(localframe(3,1)/sqrt(localframe(1,1)**2+localframe(2,1)**2+localframe(3,1)**2) )
if (S3==wr(2)) Theta_fabric = acos(localframe(3,2)/sqrt(localframe(1,2)**2+localframe(2,2)**2+localframe(3,2)**2) )
if (S3==wr(3)) Theta_fabric = acos(localframe(3,3)/sqrt(localframe(1,3)**2+localframe(2,3)**2+localframe(3,3)**2) )


Sa=0
Sa_NS=0
do i=1,nbParticules
  if (TAB_SPHERE(i)%nb_ctc<1) cycle

  TAB_SPHERE(i)%Norm_F   = TAB_SPHERE(i)%Norm_F   / TAB_SPHERE(i)%nb_ctc
  TAB_SPHERE(i)%Fabric   = TAB_SPHERE(i)%Fabric   / TAB_SPHERE(i)%nb_ctc

  Fabric  = TAB_SPHERE(i)%Fabric

  lda  = 3
  matz = 1
  call rg ( lda, 3, Fabric, wr, wi, matz, localframe, ierror )
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
  if (S3==wr(1)) Theta_c = acos(localframe(3,1)/sqrt(localframe(1,1)**2+localframe(2,1)**2+localframe(3,1)**2) )
  if (S3==wr(2)) Theta_c = acos(localframe(3,2)/sqrt(localframe(1,2)**2+localframe(2,2)**2+localframe(3,2)**2) )
  if (S3==wr(3)) Theta_c = acos(localframe(3,3)/sqrt(localframe(1,3)**2+localframe(2,3)**2+localframe(3,3)**2) )

  
  ac  = 5*(S3-S1) / 2
!  acS = ( 5*(S3-S1)/2 ) *   cos(2*(Theta_c-Theta_fabric))
  acS = ( 5*(TAB_SPHERE(i)%Fabric(3,3)-TAB_SPHERE(i)%Fabric(1,1))/2 ) !*   cos(2*(Theta_c-Theta_fabric))
    
  Sa    = Sa    + TAB_SPHERE(i)%nb_ctc*acS
  Sa_NS = Sa_NS + TAB_SPHERE(i)%nb_ctc*ac
end do

Sa  =  Sa  / real(2*cpt,8)
Sa_NS = Sa_NS / real(2*cpt,8)


if (type_box==0)&
write(105,'(31(1X,D12.5))') temps,Def,Def_H,Sa,Sa_NS

print*, Sa,Sa_NS

end subroutine contact_anisotropies_locales

!==============================================================================
! Compacity
!==============================================================================
subroutine compacity
implicit none
real(kind=8)                             :: Vol,cpt
integer                                  :: i

print*,'.... Compacity ....'

Vol = 0.D0
do i=1,nbParticules
  Vol  = Vol  + TAB_SPHERE(i)%Volume
end do

write(103,'(6(1X,D12.5))') temps,Def,Def_H,Vol/(Long*Larg*Haut)
end subroutine compacity


!==============================================================================
! Calcul de q/p
!==============================================================================
subroutine qsurp
implicit none
real(kind=8),dimension(3,3)              :: Moment,M
real(kind=8),dimension(3)                :: nik,tik,sik,Lik,Fik,Vik
real(kind=8)                             :: Rtik,Rnik,Rsik
integer                                  :: i,cd,an
real(kind=8)                             :: S1,S2,S3,Rayon_max,q,p,dmean
real(kind=8),dimension(3)                :: wr,wi
real(kind=8),dimension(3,3)              :: localframe  
integer                                  :: ierror,matz,lda

print*,'.... q/p ....'

do i=1,nbParticules
  dmean = dmean + 2*TAB_SPHERE(i)%Rmax
end do
dmean = dmean /real(nbParticules,8)

Moment(:,:)   = 0.D0
Lik = 0
Fik = 0
Vik = 0
do i=1,nb_ligneCONTACT
  cd   = TAB_CONTACT(i)%icdent
  an   = TAB_CONTACT(i)%ianent
!  if (TAB_CONTACT(i)%rn<0.000000000001*Mean_total_Normal_force) cycle
  if (cd>nbParticules.or.an>nbParticules) cycle

  nik  = TAB_CONTACT(i)%n
  tik  = TAB_CONTACT(i)%t
  sik  = TAB_CONTACT(i)%s
  Rtik = TAB_CONTACT(i)%rt
  Rnik = TAB_CONTACT(i)%rn
  Rsik = TAB_CONTACT(i)%rs

  Lik(1:3) = TAB_SPHERE(cd)%centre(1:3)-TAB_SPHERE(an)%centre(1:3)
  if ( sqrt( Lik(1)**2 + Lik(2)**2 + Lik(3)**2) > 3*( TAB_SPHERE(cd)%Rmax + TAB_SPHERE(an)%Rmax ) ) then
    Rayon_max = sqrt( ( Lik(1)**2 + Lik(2)**2 + Lik(3)**2) ) 
    lik(:) = Lik(:) / Rayon_max
    Lik(:) = abs(TAB_SPHERE(cd)%Rmax+TAB_SPHERE(an)%Rmax)*Lik(:)
  end if
  
  Fik(1:3) = (Rnik*nik(1:3)+Rtik*tik(1:3)+Rsik*sik(1:3))

  Moment(1,1:3) = Fik(1)*Lik(1:3) + Moment(1,1:3)
  Moment(2,1:3) = Fik(2)*Lik(1:3) + Moment(2,1:3)
  Moment(3,1:3) = Fik(3)*Lik(1:3) + Moment(3,1:3)
end do
Moment   = Moment   / ( Long*Larg*Haut )

pressure_in_sample = (Moment(1,1) + Moment(2,2) + Moment(3,3))/3

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

q = (S3-S1)/2
p = (S1+S2+S3)/3
write(102,'(20(1X,D12.5))') temps,Def,-Moment(1,3)/ Moment(3,3),&
                                      -Moment(2,3)/ Moment(3,3),&
                                      (S3-S1)/(S3+S1),&
                                      (S3-S2)/(S3+S2)
end subroutine qsurp


!==============================================================================
! Calcul de q/p par les moments
!==============================================================================
subroutine qsurp_moment
implicit none
real(kind=8),dimension(3,3)              :: Moment,M,Moment1,Moment2,Moment3,Moment4,Moment5,Moment6,Moment7,Moment8,Moment9
real(kind=8),dimension(3)                :: nik,tik,sik,Lik,Fik,Vik
real(kind=8)                             :: Rtik,Rnik,Rsik,Masse_reduite
integer                                  :: i,cd,an
real(kind=8)                             :: S1,S2,S3,Rayon_max,q,p,qc,pc,dmean,tan_phi_contrainte,qprime,pprime,qprime2,pprime2,&
                                            sin_phi_contrainte,alpha,&
                                            tan_phi_contrainte1,tan_phi_contrainte2,tan_phi_contrainte3,tan_phi_contrainte4,&
                                            tan_phi_contrainte5,tan_phi_contrainte6,tan_phi_contrainte7,tan_phi_contrainte8,&
                                            tan_phi_contrainte9,q1p,q2p,q3p,q4p,q5p,q6p,q7p
real(kind=8),dimension(3)                :: wr,wi
real(kind=8),dimension(3,3)              :: localframe  
integer                                  :: ierror,matz,lda

print*,'.... q/p moment....'

do i=1,nbParticules
  dmean = dmean + 2*TAB_SPHERE(i)%Rmax
  TAB_SPHERE(i)%nctc = 0
  TAB_SPHERE(i)%moment_part=0
  TAB_SPHERE(i)%moment_part1=0
  TAB_SPHERE(i)%moment_part2=0
  TAB_SPHERE(i)%moment_part3=0
  TAB_SPHERE(i)%moment_part4=0
  TAB_SPHERE(i)%moment_part5=0
  TAB_SPHERE(i)%moment_part6=0
  TAB_SPHERE(i)%moment_part7=0
  TAB_SPHERE(i)%moment_part8=0
end do
dmean = dmean /real(nbParticules,8)
do i=1,nb_ligneCONTACT
  cd = TAB_CONTACT(i)%icdent
  an = TAB_CONTACT(i)%ianent
  if (TAB_CONTACT(i)%rn<0.000000000001*Mean_total_Normal_force) cycle
  if ((TAB_SPHERE(cd)%behav/='PLEXx').or.(TAB_SPHERE(an)%behav/='PLEXx')) cycle

  TAB_SPHERE(cd)%nctc = TAB_SPHERE(cd)%nctc + 1
  TAB_SPHERE(an)%nctc = TAB_SPHERE(an)%nctc + 1
end do

Moment(:,:)   = 0.D0
Moment1=0.D0
Moment2=0.D0
Moment3=0.D0
Moment4=0.D0
Moment5=0.D0
Moment6=0.D0
Moment7=0.D0
Moment8=0.D0
Moment9=0.D0
Lik = 0
Fik = 0
Vik = 0
do i=1,nb_ligneCONTACT
  cd   = TAB_CONTACT(i)%icdent
  an   = TAB_CONTACT(i)%ianent
  if (TAB_CONTACT(i)%rn<0.000000000001*Mean_total_Normal_force) cycle
  if ((TAB_SPHERE(cd)%behav/='PLEXx').or.(TAB_SPHERE(an)%behav/='PLEXx')) cycle

  nik  = TAB_CONTACT(i)%n
  tik  = TAB_CONTACT(i)%t
  sik  = TAB_CONTACT(i)%s
  Rtik = TAB_CONTACT(i)%rt
  Rnik = TAB_CONTACT(i)%rn
  Rsik = TAB_CONTACT(i)%rs


  Lik(1:3) = TAB_CONTACT(i)%coor_ctc(1:3)-TAB_SPHERE(cd)%centre(1:3)
  if ( sqrt( Lik(1)**2 + Lik(2)**2 + Lik(3)**2) > 2*( TAB_SPHERE(cd)%Rmax + TAB_SPHERE(an)%Rmax ) ) then
    Rayon_max = sqrt( ( Lik(1)**2 + Lik(2)**2 + Lik(3)**2) ) 
    lik(:) = Lik(:) / Rayon_max
    Lik(:) = abs(TAB_SPHERE(cd)%Rmax)*Lik(:)
  end if

  Fik(1:3) = (Rnik*nik(1:3)+Rtik*tik(1:3)+Rsik*sik(1:3))
  
  TAB_SPHERE(cd)%moment_part(1,1:3) = TAB_SPHERE(cd)%moment_part(1,1:3) + Fik(1)*Lik(1:3)
  TAB_SPHERE(cd)%moment_part(2,1:3) = TAB_SPHERE(cd)%moment_part(2,1:3) + Fik(2)*Lik(1:3)
  TAB_SPHERE(cd)%moment_part(3,1:3) = TAB_SPHERE(cd)%moment_part(3,1:3) + Fik(3)*Lik(1:3)

  if (TAB_SPHERE(cd)%nctc==1) then
    TAB_SPHERE(cd)%moment_part1(1,1:3) = TAB_SPHERE(cd)%moment_part1(1,1:3) + Fik(1)*Lik(1:3)
    TAB_SPHERE(cd)%moment_part1(2,1:3) = TAB_SPHERE(cd)%moment_part1(2,1:3) + Fik(2)*Lik(1:3)
    TAB_SPHERE(cd)%moment_part1(3,1:3) = TAB_SPHERE(cd)%moment_part1(3,1:3) + Fik(3)*Lik(1:3)
  end if
  if (TAB_SPHERE(cd)%nctc==2) then
    TAB_SPHERE(cd)%moment_part2(1,1:3) = TAB_SPHERE(cd)%moment_part2(1,1:3) + Fik(1)*Lik(1:3)
    TAB_SPHERE(cd)%moment_part2(2,1:3) = TAB_SPHERE(cd)%moment_part2(2,1:3) + Fik(2)*Lik(1:3)
    TAB_SPHERE(cd)%moment_part2(3,1:3) = TAB_SPHERE(cd)%moment_part2(3,1:3) + Fik(3)*Lik(1:3)
  end if
  if (TAB_SPHERE(cd)%nctc==3) then
    TAB_SPHERE(cd)%moment_part3(1,1:3) = TAB_SPHERE(cd)%moment_part3(1,1:3) + Fik(1)*Lik(1:3)
    TAB_SPHERE(cd)%moment_part3(2,1:3) = TAB_SPHERE(cd)%moment_part3(2,1:3) + Fik(2)*Lik(1:3)
    TAB_SPHERE(cd)%moment_part3(3,1:3) = TAB_SPHERE(cd)%moment_part3(3,1:3) + Fik(3)*Lik(1:3)
  end if
  if (TAB_SPHERE(cd)%nctc==4) then
    TAB_SPHERE(cd)%moment_part4(1,1:3) = TAB_SPHERE(cd)%moment_part4(1,1:3) + Fik(1)*Lik(1:3)
    TAB_SPHERE(cd)%moment_part4(2,1:3) = TAB_SPHERE(cd)%moment_part4(2,1:3) + Fik(2)*Lik(1:3)
    TAB_SPHERE(cd)%moment_part4(3,1:3) = TAB_SPHERE(cd)%moment_part4(3,1:3) + Fik(3)*Lik(1:3)
  end if
  if (TAB_SPHERE(cd)%nctc==5) then
    TAB_SPHERE(cd)%moment_part5(1,1:3) = TAB_SPHERE(cd)%moment_part5(1,1:3) + Fik(1)*Lik(1:3)
    TAB_SPHERE(cd)%moment_part5(2,1:3) = TAB_SPHERE(cd)%moment_part5(2,1:3) + Fik(2)*Lik(1:3)
    TAB_SPHERE(cd)%moment_part5(3,1:3) = TAB_SPHERE(cd)%moment_part5(3,1:3) + Fik(3)*Lik(1:3)
  end if
  if (TAB_SPHERE(cd)%nctc==6) then
    TAB_SPHERE(cd)%moment_part6(1,1:3) = TAB_SPHERE(cd)%moment_part6(1,1:3) + Fik(1)*Lik(1:3)
    TAB_SPHERE(cd)%moment_part6(2,1:3) = TAB_SPHERE(cd)%moment_part6(2,1:3) + Fik(2)*Lik(1:3)
    TAB_SPHERE(cd)%moment_part6(3,1:3) = TAB_SPHERE(cd)%moment_part6(3,1:3) + Fik(3)*Lik(1:3)
  end if
  if (TAB_SPHERE(cd)%nctc==7) then
    TAB_SPHERE(cd)%moment_part7(1,1:3) = TAB_SPHERE(cd)%moment_part7(1,1:3) + Fik(1)*Lik(1:3)
    TAB_SPHERE(cd)%moment_part7(2,1:3) = TAB_SPHERE(cd)%moment_part7(2,1:3) + Fik(2)*Lik(1:3)
    TAB_SPHERE(cd)%moment_part7(3,1:3) = TAB_SPHERE(cd)%moment_part7(3,1:3) + Fik(3)*Lik(1:3)
  end if
  if (TAB_SPHERE(cd)%nctc==8) then
    TAB_SPHERE(cd)%moment_part8(1,1:3) = TAB_SPHERE(cd)%moment_part8(1,1:3) + Fik(1)*Lik(1:3)
    TAB_SPHERE(cd)%moment_part8(2,1:3) = TAB_SPHERE(cd)%moment_part8(2,1:3) + Fik(2)*Lik(1:3)
    TAB_SPHERE(cd)%moment_part8(3,1:3) = TAB_SPHERE(cd)%moment_part8(3,1:3) + Fik(3)*Lik(1:3)
  end if
  if (TAB_SPHERE(cd)%nctc==9) then
    TAB_SPHERE(cd)%moment_part9(1,1:3) = TAB_SPHERE(cd)%moment_part9(1,1:3) + Fik(1)*Lik(1:3)
    TAB_SPHERE(cd)%moment_part9(2,1:3) = TAB_SPHERE(cd)%moment_part9(2,1:3) + Fik(2)*Lik(1:3)
    TAB_SPHERE(cd)%moment_part9(3,1:3) = TAB_SPHERE(cd)%moment_part9(3,1:3) + Fik(3)*Lik(1:3)
  end if

  Lik(1:3) = TAB_CONTACT(i)%coor_ctc(1:3)-TAB_SPHERE(an)%centre(1:3)
  if ( sqrt( Lik(1)**2 + Lik(2)**2 + Lik(3)**2) > 2*( TAB_SPHERE(cd)%Rmax + TAB_SPHERE(an)%Rmax ) ) then
    Rayon_max = sqrt( ( Lik(1)**2 + Lik(2)**2 + Lik(3)**2) ) 
    lik(:) = Lik(:) / Rayon_max
    Lik(:) = abs(TAB_SPHERE(an)%Rmax)*Lik(:)
  end if

  Fik(1:3) = -Fik(1:3)
  
  TAB_SPHERE(an)%moment_part(1,1:3) = TAB_SPHERE(an)%moment_part(1,1:3) + Fik(1)*Lik(1:3)
  TAB_SPHERE(an)%moment_part(2,1:3) = TAB_SPHERE(an)%moment_part(2,1:3) + Fik(2)*Lik(1:3)
  TAB_SPHERE(an)%moment_part(3,1:3) = TAB_SPHERE(an)%moment_part(3,1:3) + Fik(3)*Lik(1:3)

  if (TAB_SPHERE(an)%nctc==1) then
    TAB_SPHERE(an)%moment_part1(1,1:3) = TAB_SPHERE(an)%moment_part1(1,1:3) + Fik(1)*Lik(1:3)
    TAB_SPHERE(an)%moment_part1(2,1:3) = TAB_SPHERE(an)%moment_part1(2,1:3) + Fik(2)*Lik(1:3)
    TAB_SPHERE(an)%moment_part1(3,1:3) = TAB_SPHERE(an)%moment_part1(3,1:3) + Fik(3)*Lik(1:3)
  end if
  if (TAB_SPHERE(an)%nctc==2) then
    TAB_SPHERE(an)%moment_part2(1,1:3) = TAB_SPHERE(an)%moment_part2(1,1:3) + Fik(1)*Lik(1:3)
    TAB_SPHERE(an)%moment_part2(2,1:3) = TAB_SPHERE(an)%moment_part2(2,1:3) + Fik(2)*Lik(1:3)
    TAB_SPHERE(an)%moment_part2(3,1:3) = TAB_SPHERE(an)%moment_part2(3,1:3) + Fik(3)*Lik(1:3)
  end if
  if (TAB_SPHERE(an)%nctc==3) then
    TAB_SPHERE(an)%moment_part3(1,1:3) = TAB_SPHERE(an)%moment_part3(1,1:3) + Fik(1)*Lik(1:3)
    TAB_SPHERE(an)%moment_part3(2,1:3) = TAB_SPHERE(an)%moment_part3(2,1:3) + Fik(2)*Lik(1:3)
    TAB_SPHERE(an)%moment_part3(3,1:3) = TAB_SPHERE(an)%moment_part3(3,1:3) + Fik(3)*Lik(1:3)
  end if
  if (TAB_SPHERE(cd)%nctc==4) then
    TAB_SPHERE(an)%moment_part4(1,1:3) = TAB_SPHERE(an)%moment_part4(1,1:3) + Fik(1)*Lik(1:3)
    TAB_SPHERE(an)%moment_part4(2,1:3) = TAB_SPHERE(an)%moment_part4(2,1:3) + Fik(2)*Lik(1:3)
    TAB_SPHERE(an)%moment_part4(3,1:3) = TAB_SPHERE(an)%moment_part4(3,1:3) + Fik(3)*Lik(1:3)
  end if
  if (TAB_SPHERE(an)%nctc==5) then
    TAB_SPHERE(an)%moment_part5(1,1:3) = TAB_SPHERE(an)%moment_part5(1,1:3) + Fik(1)*Lik(1:3)
    TAB_SPHERE(an)%moment_part5(2,1:3) = TAB_SPHERE(an)%moment_part5(2,1:3) + Fik(2)*Lik(1:3)
    TAB_SPHERE(an)%moment_part5(3,1:3) = TAB_SPHERE(an)%moment_part5(3,1:3) + Fik(3)*Lik(1:3)
  end if
  if (TAB_SPHERE(an)%nctc==6) then
    TAB_SPHERE(an)%moment_part6(1,1:3) = TAB_SPHERE(an)%moment_part6(1,1:3) + Fik(1)*Lik(1:3)
    TAB_SPHERE(an)%moment_part6(2,1:3) = TAB_SPHERE(an)%moment_part6(2,1:3) + Fik(2)*Lik(1:3)
    TAB_SPHERE(an)%moment_part6(3,1:3) = TAB_SPHERE(an)%moment_part6(3,1:3) + Fik(3)*Lik(1:3)
  end if
  if (TAB_SPHERE(an)%nctc==7) then
    TAB_SPHERE(an)%moment_part7(1,1:3) = TAB_SPHERE(an)%moment_part7(1,1:3) + Fik(1)*Lik(1:3)
    TAB_SPHERE(an)%moment_part7(2,1:3) = TAB_SPHERE(an)%moment_part7(2,1:3) + Fik(2)*Lik(1:3)
    TAB_SPHERE(an)%moment_part7(3,1:3) = TAB_SPHERE(an)%moment_part7(3,1:3) + Fik(3)*Lik(1:3)
  end if
  if (TAB_SPHERE(an)%nctc==8) then
    TAB_SPHERE(an)%moment_part8(1,1:3) = TAB_SPHERE(an)%moment_part8(1,1:3) + Fik(1)*Lik(1:3)
    TAB_SPHERE(an)%moment_part8(2,1:3) = TAB_SPHERE(an)%moment_part8(2,1:3) + Fik(2)*Lik(1:3)
    TAB_SPHERE(an)%moment_part8(3,1:3) = TAB_SPHERE(an)%moment_part8(3,1:3) + Fik(3)*Lik(1:3)
  end if
  if (TAB_SPHERE(an)%nctc==9) then
    TAB_SPHERE(an)%moment_part9(1,1:3) = TAB_SPHERE(an)%moment_part9(1,1:3) + Fik(1)*Lik(1:3)
    TAB_SPHERE(an)%moment_part9(2,1:3) = TAB_SPHERE(an)%moment_part9(2,1:3) + Fik(2)*Lik(1:3)
    TAB_SPHERE(an)%moment_part9(3,1:3) = TAB_SPHERE(an)%moment_part9(3,1:3) + Fik(3)*Lik(1:3)
  end if

end do
lda  = 3
matz = 1
do i=1,nbParticules
  if ((TAB_SPHERE(i)%behav/='PLEXx')) cycle
  Moment(:,:) = Moment(:,:) - TAB_SPHERE(i)%moment_part(:,:)
  if (TAB_SPHERE(i)%nctc==1) Moment1(:,:) = Moment1(:,:) - TAB_SPHERE(i)%moment_part1(:,:)
  if (TAB_SPHERE(i)%nctc==2) Moment2(:,:) = Moment2(:,:) - TAB_SPHERE(i)%moment_part2(:,:)
  if (TAB_SPHERE(i)%nctc==3) Moment3(:,:) = Moment3(:,:) - TAB_SPHERE(i)%moment_part3(:,:)
  if (TAB_SPHERE(i)%nctc==4) Moment4(:,:) = Moment4(:,:) - TAB_SPHERE(i)%moment_part4(:,:)
  if (TAB_SPHERE(i)%nctc==5) Moment5(:,:) = Moment5(:,:) - TAB_SPHERE(i)%moment_part5(:,:)
  if (TAB_SPHERE(i)%nctc==6) Moment6(:,:) = Moment6(:,:) - TAB_SPHERE(i)%moment_part6(:,:)
  if (TAB_SPHERE(i)%nctc==7) Moment7(:,:) = Moment7(:,:) - TAB_SPHERE(i)%moment_part7(:,:)
  if (TAB_SPHERE(i)%nctc==8) Moment8(:,:) = Moment8(:,:) - TAB_SPHERE(i)%moment_part8(:,:)
  if (TAB_SPHERE(i)%nctc==9) Moment9(:,:) = Moment9(:,:) - TAB_SPHERE(i)%moment_part9(:,:)
end do


tan_phi_contrainte = -Moment(2,3)/ Moment(3,3)
tan_phi_contrainte1= -Moment1(2,3)/ Moment(3,3)
tan_phi_contrainte2= -Moment2(2,3)/ Moment(3,3)
tan_phi_contrainte3= -Moment3(2,3)/ Moment(3,3)
tan_phi_contrainte4= -Moment4(2,3)/ Moment(3,3)
tan_phi_contrainte5= -Moment5(2,3)/ Moment(3,3)
tan_phi_contrainte6= -Moment6(2,3)/ Moment(3,3)
tan_phi_contrainte7= -Moment7(2,3)/ Moment(3,3)
tan_phi_contrainte8= -Moment8(2,3)/ Moment(3,3)
tan_phi_contrainte9= -Moment9(2,3)/ Moment(3,3)

M = Moment
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

qprime = (S3-S1)/2
pprime = (S3+S1)/2

M = Moment1
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
q1p = (S3-S1)/2
q1p = q1p / pprime


M = Moment2
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
q2p = (S3-S1)/2
q2p = q2p / pprime

M = Moment3
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
q3p = (S3-S1)/2
q3p = q3p / pprime

M = Moment4
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
q4p = (S3-S1)/2
q4p = q4p / pprime

M = Moment5
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
q5p = (S3-S1)/2
q5p = q5p / pprime

M = Moment6
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
q6p = (S3-S1)/2
q6p = q6p / pprime

M = Moment1
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
q7p = (S3-S1)/2
q7p = q7p / pprime


if (type_box==0) write(1003,'(20(1X,D12.5))') temps,Def,qprime/pprime,&
                                              qprime/pprime,&
                                              q1p,&
                                              q2p,&
                                              q3p,&
                                              q4p,&
                                              q5p,&
                                              q6p,&
                                              q7p,&
                                              q1p+q2p+q3p+q4p+q5p+q6p+q7p
end subroutine qsurp_moment

!==============================================================================
! Vitesse moyenne
!==============================================================================
subroutine vitesse_moyenne
implicit none
real(kind=8)                             :: V,Vx,Vy,Vz,Vrx,cpt,el_diameter,zc,viscosity,&
                                            I_Calculated,Ec,V_centre_mass,&
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
V_centre_mass = 0
Nombre_Inertie = 0
mean_mass = 0

cpt = 0
do i=1,nbParticules
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
applied_pressure = 6.5790026 / (Long*Larg)

I_Calculated   = (0.047472569/Haut) * el_diameter * sqrt( 1000./ pressure_in_sample )
Nombre_Inertie = (0.047472569/Haut) * el_diameter * sqrt( 1000./ applied_pressure   )

write(101,'(25(1X,D12.5))') temps,Def,Haut,&
                            V,Vx,Vy,Vz,Vrx,&
                            Nombre_Inertie,I_Calculated

 end subroutine vitesse_moyenne

!==============================================================================
! Nombre de coordination
!==============================================================================
subroutine nb_coordination
implicit none
real(kind=8)                             :: z,zc,zp,zc2,mean_gap,min_gap,dmoyen,viscosity,mean_force_calculated,Mobilized,kslide,&
                                            P0,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,zc3,float_part,nb_part_without_paroi,&
                                            z_float,dist,&
                                            P1_all,P2_all,P3_all,P4_all,P5_all,P6_all,P7_all,P8_all,P9_all,P10_all,P11_all,P12_all
integer                                  :: i,j,cd,an
real(kind=8),dimension(:),allocatable    :: ctc_part

print*,'.... Nombre coordination ....'

z   = 0
zc  = 0
zp  = 0
zc2 = 0
mean_gap = 0.D0
min_gap  = 0.D0
dmoyen   = 0.D0
mean_force_calculated = 0.D0
Mobilized = 0.D0
kslide    = 0.D0

zc3= 0
P0 = 0
P1 = 0
P2 = 0
P2 = 0
P3 = 0
P4 = 0
P5 = 0
P6 = 0
P7 = 0
P8 = 0
P9 = 0
P10= 0
P11= 0
P12= 0

P1_all = 0
P2_all = 0
P2_all = 0
P3_all = 0
P4_all = 0
P5_all = 0
P6_all = 0
P7_all = 0
P8_all = 0
P9_all = 0
P10_all= 0
P11_all= 0
P12_all= 0

float_part = 0
nb_part_without_paroi = 0

z_float = 0

do i= 1,nb_ligneCONTACT
!  if (TAB_CONTACT(i)%rn==0.D0) cycle
  if (TAB_CONTACT(i)%rn<0.000000000001*Mean_total_Normal_force) cycle
  zc = zc+1
  mean_gap = mean_gap + abs(TAB_CONTACT(i)%gap)
  min_gap  = max( min_gap , abs(TAB_CONTACT(i)%gap) )
  mean_force_calculated = mean_force_calculated + &
                          sqrt( TAB_CONTACT(i)%rn**2 + & 
                                TAB_CONTACT(i)%rt**2 + &
                                TAB_CONTACT(i)%rs**2    )
                                
  Mobilized = Mobilized + sqrt( TAB_CONTACT(i)%rt**2 + TAB_CONTACT(i)%rs**2 )/(0.4*TAB_CONTACT(i)%rn)
  if (TAB_CONTACT(i)%statut=='slide') kslide = kslide+1.0
  
end do
mean_force_calculated = mean_force_calculated / zc

do i=1,nbParticules
  TAB_SPHERE(i)%nctc = 0
  TAB_SPHERE(i)%ctc_float = 0
  if (TAB_SPHERE(i)%behav=='PLEXx') nb_part_without_paroi = nb_part_without_paroi + 1.0
end do

do i=1,nb_ligneCONTACT
  cd = TAB_CONTACT(i)%icdent
  an = TAB_CONTACT(i)%ianent
  if (TAB_CONTACT(i)%rn<0.000000000001*Mean_total_Normal_force) cycle
  if ((TAB_SPHERE(cd)%behav/='PLEXx').or.(TAB_SPHERE(an)%behav/='PLEXx')) cycle

  TAB_SPHERE(cd)%nctc = TAB_SPHERE(cd)%nctc + 1
  TAB_SPHERE(an)%nctc = TAB_SPHERE(an)%nctc + 1
end do
    
do i=1,nbParticules
  if (TAB_SPHERE(i)%nctc > 0 ) then
    zp=zp+1
    dmoyen = dmoyen + TAB_SPHERE(i)%Rmax
  end if
end do
dmoyen = dmoyen / zp

do i=1,nbParticules
  if ((TAB_SPHERE(i)%behav/='PLEXx')) cycle
  
  if (TAB_SPHERE(i)%nctc >= 2 )  zc2  = zc2 + real(TAB_SPHERE(i)%nctc,8)
  if (TAB_SPHERE(i)%nctc == 0 )  P0   = P0  + 1.0
  if (TAB_SPHERE(i)%nctc == 1 )  P1   = P1  + 1.0
  if (TAB_SPHERE(i)%nctc == 2 )  P2   = P2  + 1.0
  if (TAB_SPHERE(i)%nctc == 3 )  P3   = P3  + 1.0
  if (TAB_SPHERE(i)%nctc == 4 )  P4   = P4  + 1.0
  if (TAB_SPHERE(i)%nctc == 5 )  P5   = P5  + 1.0
  if (TAB_SPHERE(i)%nctc == 6 )  P6   = P6  + 1.0
  if (TAB_SPHERE(i)%nctc == 7 )  P7   = P7  + 1.0
  if (TAB_SPHERE(i)%nctc == 8 )  P8   = P8  + 1.0
  if (TAB_SPHERE(i)%nctc == 9 )  P9   = P9  + 1.0
  if (TAB_SPHERE(i)%nctc == 10 ) P10  = P10 + 1.0
  if (TAB_SPHERE(i)%nctc == 11 ) P11  = P11 + 1.0
  if (TAB_SPHERE(i)%nctc == 12 ) P12  = P12 + 1.0
end do

Mobilized = Mobilized / zc
kslide    = kslide / zc

P0  =     P0  / nb_part_without_paroi

P1_all  = 1 * P1  / nb_part_without_paroi
P2_all  = 2 * P2  / nb_part_without_paroi
P3_all  = 3 * P3  / nb_part_without_paroi
P4_all  = 4 * P4  / nb_part_without_paroi
P5_all  = 5 * P5  / nb_part_without_paroi
P6_all  = 6 * P6  / nb_part_without_paroi
P7_all  = 7 * P7  / nb_part_without_paroi
P8_all  = 8 * P8  / nb_part_without_paroi
P9_all  = 9 * P9  / nb_part_without_paroi
P10_all = 10* P10 / nb_part_without_paroi
P11_all = 11* P11 / nb_part_without_paroi
P12_all = 12* P12 / nb_part_without_paroi
zc2 = P1_all+P2_all+P3_all+P4_all+P5_all+P6_all+P7_all+P8_all+P9_all+P10_all+P11_all+P12_all

P1  = 1 * P1  / zp
P2  = 2 * P2  / zp
P3  = 3 * P3  / zp
P4  = 4 * P4  / zp
P5  = 5 * P5  / zp
P6  = 6 * P6  / zp
P7  = 7 * P7  / zp
P8  = 8 * P8  / zp
P9  = 9 * P9  / zp
P10 = 10* P10 / zp
P11 = 11* P11 / zp
P12 = 12* P12 / zp
zc3 = P1+P2+P3+P4+P5+P6+P7+P8+P9+P10+P11+P12


z   = 2*zc / zp

! zp  =  nombre de particules avec au moins un contact !!!
! zc3 =  coordination sans les partiucles flottantes
! zc2 =  coordination avec les particules flottantes
! z = 2 fois le nombre de contacts divisé par le nombre particules qui ne flottent pas...


float_part = float_part /( zp+float_part )

mean_gap = (mean_gap / zc) / dmoyen
min_gap  = min_gap / dmoyen


write(100,'(40(1X,D12.5))') temps,Def,Def_H,z,zc2,zc3,kslide,Mobilized,&
                                             P0,P1,P2,P3,P4,P5,P6,P7,P8,P9,&
                                             P1_all,P2_all,P3_all,P4_all,P5_all,P6_all,P7_all,P8_all,P9_all,&
                                             mean_force_calculated,pressure_in_sample


end subroutine nb_coordination


!==================================================================================================
!Procedure pour ecrire le gmv a chaque pas de temps
!==================================================================================================
subroutine write_all_gmv
implicit none
character(len=13)                           ::  nom
integer                                     ::  i,j,cpt_force,cd,an,cpt_force_contact
real(kind=8),dimension(3)                   ::  nik,tik,sik,Lik,Fik
real(kind=8)                                ::  Rtik,Rnik,Rsik,LN,fmax,fmin,fmean,ax2,force_contact,FN,FN_contact,&
                                                fmax_contact,fmean_contact,fmin_contact,alpha,mean_diameter
print*,'....  WRITE GMV FILE....'


if ( (step_all_gmv == compteur_clout).or. (step_all_gmv < 0) ) then

  nom       =  'GMV_FILE.0000'
  if (compteur_clout<10) then
    WRITE(nom(12:13),'(I1)')   compteur_clout
  else if ( (compteur_clout>=10) .and. (compteur_clout<100) ) then
    WRITE(nom(11:13),'(I2)')   compteur_clout
  else if ( (compteur_clout>=100).and. (compteur_clout<1000) ) then
    WRITE(nom(10:13),'(I3)')   compteur_clout
  end if   

!!! -------------- CALCUL DE LA FORCE RADIAL MEAN, MAX, MIN ---------------
  fmax           = 0
  fmean          = 0
  fmin           = 100000
  cpt_force      = 0
  fmax_contact           = 0
  fmean_contact          = 0
  fmin_contact           = 100000
  cpt_force_contact      = 0
  Fik(:)         = 0
  FN_contact     = 0
  mean_diameter  = 0.D0

  do i=1,nb_ligneCONTACT
    if  (TAB_CONTACT(i)%nature == 'SPPLx') cycle 
    if  (TAB_CONTACT(i)%rn<0.000000000001*Mean_total_Normal_force) cycle
    cd   = TAB_CONTACT(i)%icdent
    an   = TAB_CONTACT(i)%ianent
    Fik(1:3) = 0.D0
    nik  = TAB_CONTACT(i)%n
    tik  = TAB_CONTACT(i)%t
    sik  = TAB_CONTACT(i)%s
    Rtik = TAB_CONTACT(i)%rt
    Rnik = TAB_CONTACT(i)%rn
    Rsik = TAB_CONTACT(i)%rs
    TAB_CONTACT(i)%cycle_gmv = .false.
    Lik(1:3) = TAB_SPHERE(cd)%centre(1:3)-TAB_SPHERE(an)%centre(1:3)
    if ( sqrt( Lik(1)**2 + Lik(2)**2 + Lik(3)**2) > 2*( TAB_SPHERE(cd)%Rmax + TAB_SPHERE(an)%Rmax ) ) then
      TAB_CONTACT(i)%cycle_gmv = .true.
!      cpt_cycle_gmv = cpt_cycle_gmv + 1
      cycle
    end if

    Fik(1:3) = Fik(1:3) + (Rnik*nik(1:3)+Rtik*tik(1:3)+Rsik*sik(1:3))

    FN       =  Rnik
    
    fmax   = max(fmax,FN)
    fmin   = min(fmin,FN)
    fmean  = fmean + FN
    cpt_force = cpt_force + 1
  end do
  fmean = fmean / real(cpt_force,8)

  do i=1,nbParticules
    TAB_SPHERE(i)%nctc = 0
    mean_diameter = mean_diameter + TAB_SPHERE(i)%Rmax
  end do
  mean_diameter = mean_diameter / real(nbParticules,8)

  do i=1,nb_ligneCONTACT
    cd = TAB_CONTACT(i)%icdent
    an = TAB_CONTACT(i)%ianent
    if (TAB_CONTACT(i)%rn<0.000000000001*Mean_total_Normal_force) cycle
    if ((TAB_SPHERE(cd)%behav/='PLEXx').or.(TAB_SPHERE(an)%behav/='PLEXx')) cycle

    TAB_SPHERE(cd)%nctc = TAB_SPHERE(cd)%nctc + 1
    TAB_SPHERE(an)%nctc = TAB_SPHERE(an)%nctc + 1
  end do

  print*,fmin
  print*,fmax
  print*,fmean

!!! -------------- FIN CALCUL DE LA FORCE RADIAL/Normale MEAN, MAX, MIN ---------------


  open(unit=30000,file=nom,status='replace') 

  write(30000,'(A14)') 'gmvinput ascii' 
  write(30000,'(A5,1X,I9)') 'nodes',nbParticules*48 + (cpt_force)*8  ! 1 sphere c'est 48 nodes

!!! ----------- ECRITURE DES NOEUDS DES SPHERES
  do i=1,nbParticules ! X Sphere 
    call write_gmv_x(30000,TAB_SPHERE(i)%centre(1),TAB_SPHERE(i)%Rmax) 
  enddo 
!!! -------------- ECRITURE DES NOEUDS DES FORCES RADIALES suivant x
  alpha = 1.5
  Ax2 = 0.D0
  cd  = 0
  an  = 0
  Fik(:) = 0.D0
  do i=1,nb_ligneCONTACT
    if  (TAB_CONTACT(i)%nature == 'SPPLx') cycle 
    if  (TAB_CONTACT(i)%cycle_gmv) cycle
    if  (TAB_CONTACT(i)%rn<0.000000000001*Mean_total_Normal_force) cycle
    cd   = TAB_CONTACT(i)%icdent
    an   = TAB_CONTACT(i)%ianent
    Fik(1:3) = 0.D0
    nik  = TAB_CONTACT(i)%n
    tik  = TAB_CONTACT(i)%t
    sik  = TAB_CONTACT(i)%s
    Rtik = TAB_CONTACT(i)%rt
    Rnik = TAB_CONTACT(i)%rn
    Rsik = TAB_CONTACT(i)%rs

    Fik(1:3) = 0.D0 
    force_contact =  Rnik
        
!    Ax2 = alpha * (TAB_SPHERE(cd)%Rmax + TAB_SPHERE(an)%Rmax) *( (force_contact - fmin)/(fmax - fmin) )
    Ax2 = mean_diameter / 4.0
    write(30000,'(8(1X,E14.7))')  TAB_SPHERE(cd)%centre(1)-Ax2,TAB_SPHERE(cd)%centre(1)+Ax2,&
                                  TAB_SPHERE(cd)%centre(1)+Ax2,TAB_SPHERE(cd)%centre(1)-Ax2,&
                                  TAB_SPHERE(an)%centre(1)-Ax2,TAB_SPHERE(an)%centre(1)+Ax2,&
                                  TAB_SPHERE(an)%centre(1)+Ax2,TAB_SPHERE(an)%centre(1)-Ax2
  end do

!!! ----------- ECRITURE DES NOEUDS DES SPHERES
  do i=1,nbParticules ! Y Sphere 
    call write_gmv_y(30000,TAB_SPHERE(i)%centre(2),TAB_SPHERE(i)%Rmax) 
  enddo 
!!! -------------- ECRITURE DES NOEUDS DES FORCES suivant y
  cd  = 0
  an  = 0
  Fik(:) = 0.D0
  do i=1,nb_ligneCONTACT
    if  (TAB_CONTACT(i)%nature == 'SPPLx') cycle 
    if  (TAB_CONTACT(i)%cycle_gmv) cycle
    if  (TAB_CONTACT(i)%rn<0.000000000001*Mean_total_Normal_force) cycle
    cd   = TAB_CONTACT(i)%icdent
    an   = TAB_CONTACT(i)%ianent
    Fik(1:3) = 0.D0
    nik  = TAB_CONTACT(i)%n
    tik  = TAB_CONTACT(i)%t
    sik  = TAB_CONTACT(i)%s
    Rtik = TAB_CONTACT(i)%rt
    Rnik = TAB_CONTACT(i)%rn
    Rsik = TAB_CONTACT(i)%rs

    Fik(1:3) = 0.D0 
    force_contact =  Rnik

!    Ax2 = alpha * (TAB_SPHERE(cd)%Rmax + TAB_SPHERE(an)%Rmax) *( (force_contact - fmin)/(fmax - fmin) )
    Ax2 = mean_diameter / 4.0
    write(30000,'(8(1X,E14.7))')  TAB_SPHERE(cd)%centre(2)-Ax2,TAB_SPHERE(cd)%centre(2)-Ax2,&
                                  TAB_SPHERE(cd)%centre(2)+Ax2,TAB_SPHERE(cd)%centre(2)+Ax2,&
                                  TAB_SPHERE(an)%centre(2)-Ax2,TAB_SPHERE(an)%centre(2)-Ax2,&
                                  TAB_SPHERE(an)%centre(2)+Ax2,TAB_SPHERE(an)%centre(2)+Ax2
  end do

!!! ----------- ECRITURE DES NOEUDS DES SPHERES
  do i=1,nbParticules ! Z Sphere 
    call write_gmv_z(30000,TAB_SPHERE(i)%centre(3),TAB_SPHERE(i)%Rmax) 
  enddo 
!!! -------------- ECRITURE DES NOEUDS DES FORCES suivant z
  cd  = 0
  an  = 0
  Fik(:) = 0.D0
  do i=1,nb_ligneCONTACT
    if  (TAB_CONTACT(i)%nature == 'SPPLx') cycle 
    if  (TAB_CONTACT(i)%cycle_gmv) cycle
    if  (TAB_CONTACT(i)%rn<0.000000000001*Mean_total_Normal_force) cycle
    cd   = TAB_CONTACT(i)%icdent
    an   = TAB_CONTACT(i)%ianent
    Fik(1:3) = 0.D0
    nik  = TAB_CONTACT(i)%n
    tik  = TAB_CONTACT(i)%t
    sik  = TAB_CONTACT(i)%s
    Rtik = TAB_CONTACT(i)%rt
    Rnik = TAB_CONTACT(i)%rn
    Rsik = TAB_CONTACT(i)%rs

    Fik(1:3) = 0.D0 
    force_contact =  Rnik

!    Ax2 = alpha * (TAB_SPHERE(cd)%Rmax + TAB_SPHERE(an)%Rmax) *( (force_contact - fmin)/(fmax - fmin) )
    Ax2 = mean_diameter / 4.0
    write(30000,'(8(1X,E14.7))')  TAB_SPHERE(cd)%centre(3),TAB_SPHERE(cd)%centre(3),&
                                  TAB_SPHERE(cd)%centre(3),TAB_SPHERE(cd)%centre(3),&
                                  TAB_SPHERE(an)%centre(3),TAB_SPHERE(an)%centre(3),&
                                  TAB_SPHERE(an)%centre(3),TAB_SPHERE(an)%centre(3)
  end do

  call write_gmv_cells(30000,nbParticules,cpt_force) 
  call write_gmv_surf(30000,nbParticules ) 
  call write_gmv_mat(30000,nbParticules*62,cpt_force,0) 

end if

end subroutine write_all_gmv

!======================================================================
!Lecture des fichiers Vloc, Dof et BODIES et on commence  remplir notre type
!TAB_SPHER.... c'est un peu bancale je l'accorde volontier :-)
!======================================================================
subroutine read_Vloc_dof_bodies
implicit none
integer                       ::  i=0,j=0,n_spsp,n_sppl,num_part,cpt_type,err
integer                       ::  nb_SPSPx,nb_SPPLx,icdent,ianent
real(kind=8)                  ::  rn,rs,rt,t1,t2,t3,n1,n2,n3,s1,s2,s3,rayon,vls,vln,vlt,gap
character(len=6)              ::  text
character(len=5)              ::  statut,behav,color
character(len=13)             ::  text2
real(kind=8)                  ::  centre1,centre2,centre3,ax1,ax2,ax3,coor1,coor2,coor3,&
                                  centrecluster1,centrecluster2,centrecluster3
real(kind=8)                  ::  coor_ctc1,coor_ctc2,coor_ctc3
real(kind=8)                  ::  som1,som2,som3,Volume
real(kind=8),dimension(3)     ::  X,Y,Z,centre_gravity
real(kind=8),dimension(3,3)   ::  rot
real(kind=8)                  ::  Vrx,Vry,Vrz,vx,vy,vz
real(kind=8)                  ::  I1,I2,I3,Number_mean_Total_force
character(len=22)             ::  clout_DOF
character(len=28)             ::  clout_Vloc


if (all_pas==0) then 
  compteur_clout = compteur_clout+1
  clout_DOF  =  '../OUTBOX/DOF.OUT.   '
  clout_Vloc =  '../OUTBOX/Vloc_Rloc.OUT.   '
else
  compteur_clout = all_pas
  clout_DOF  =  '../OUTBOX/DOF.OUT.   '
  clout_Vloc =  '../OUTBOX/Vloc_Rloc.OUT.   '   
end if

if ( deformation_clout .and. deformation_ini ) then
  clout_DOF  =  '../OUTBOX/DOF.OUT.0'
  clout_Vloc =  '../OUTBOX/Vloc_Rloc.OUT.0'
  compteur_clout = 0
else
  if (compteur_clout<10) then
    WRITE(clout_DOF(19:20),'(I1)')   compteur_clout
    WRITE(clout_Vloc(25:26),'(I1)')  compteur_clout
  else if ( (compteur_clout>=10) .and. (compteur_clout<100) ) then
    WRITE(clout_DOF(19:21),'(I2)')   compteur_clout
    WRITE(clout_Vloc(25:27),'(I2)')  compteur_clout
  else if ( (compteur_clout>=100).and. (compteur_clout<1000) ) then
    WRITE(clout_DOF(19:22),'(I3)')   compteur_clout
    WRITE(clout_Vloc(25:28),'(I3)')  compteur_clout
    WRITE(clout_Vloc(1:2),'(A2)')  '..'
    WRITE(clout_DOF(1:2),'(A2)')  '..'
  end if   
end if

if (the_first_time) then
  open(unit=2,file='../OUTBOX/BODIES.OUT',status='old')
  allocate(TAB_SPHERE(nbParticules))
  allocate(TAB_PLAN(nb_paroi))
  i=0
  do
    read(2,'(A6)') text
    if (text == '$bdyty') then
      i = i +1
      if (i<nbParticules+1) then
        TAB_SPHERE(i)%num=i
        read(2,*)
        read(2,*)
        read(2,'(22X,A5)') behav
        TAB_SPHERE(i)%behav = behav
        read(2,'(29X, 3(5x,D14.7,2X))') I1,I2,I3
        read(2,*)
        read(2,'(29X, 3(5x,D14.7,2X))') centre1,centre2,centre3
        TAB_SPHERE(i)%I1=I1
        TAB_SPHERE(i)%I2=I2
        TAB_SPHERE(i)%I3=I3
        TAB_SPHERE(i)%centre_ref(1)=centre1
        TAB_SPHERE(i)%centre_ref(2)=centre2
        TAB_SPHERE(i)%centre_ref(3)=centre3
        read(2,*)
        read(2,*)
        read(2,'(22X,A5,7X,D14.7)') color,rayon
        TAB_SPHERE(i)%color = color
        TAB_SPHERE(i)%Rmax = rayon
        TAB_SPHERE(i)%volume = 4.188790204*rayon**3
        TAB_SPHERE(i)%centre =0.D0
        TAB_SPHERE(i)%inbox=.false.
        TAB_SPHERE(i)%momentI=0
        TAB_SPHERE(i)%nctc = 0
        TAB_SPHERE(i)%nctcslide = 0
        TAB_SPHERE(i)%nctcstick = 0
      end if
      if (i>nbParticules) then
        the_first_time=.false.
        Print*, 'Lecture Bodies.out--> OK'
        exit
      end if
    end if
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
  read(5,'(7X,i9,19X,D14.7)') pas,temps
  print*,'pas de temps n:',pas, temps
  print*,'--->',clout_DOF
  read(5,*)
  read(5,*)
  do
    read(5,'(A6)') text
    if (text == '!-----') exit
    if (text == '$bdyty') then
      i = i+1
      read(5,'(6X,i9)') num_part
      if (num_part<nbParticules+1) then
        read(5,*)
        read(5,'(29X, 3(5x,D14.7,2X))') coor1,coor2,coor3
        TAB_SPHERE(i)%centre(1)=coor1+TAB_SPHERE(i)%centre_ref(1)
        TAB_SPHERE(i)%centre(2)=coor2+TAB_SPHERE(i)%centre_ref(2)
        TAB_SPHERE(i)%centre(3)=coor3+TAB_SPHERE(i)%centre_ref(3)
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
      end if
    end if
  end do
end if



nb_SPSPx=0
nb_SPPLx=0
open(unit=4,file=clout_Vloc,iostat=err,status='old')
if (err/=0) then

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
  !---------------------------------------
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
  read(4,'(29X,5x,D14.7)') gap
  read(4,'(29X, 3(5x,D14.7,2X))') s1,s2,s3
  read(4,'(29X, 3(5x,D14.7,2X))') t1,t2,t3
  read(4,'(29X, 3(5x,D14.7,2X))') n1,n2,n3
  read(4,'(29X, 3(5x,D14.7,2X))') coor_ctc1,coor_ctc2,coor_ctc3
  !--------------------------
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
  TAB_CONTACT(n_spsp)%nature = 'SPSPx'
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

print*,'Nb particules :',nbParticules,'Nb Contacts :',nb_ligneCONTACT
print*,'Mean Normal force :',Mean_total_Normal_force 

end subroutine read_Vloc_dof_bodies


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
  Lik(1:3) = TAB_SPHERE(cd)%centre(1:3)-TAB_SPHERE(an)%centre(1:3)
  norm = sqrt( Lik(1)**2 + Lik(2)**2 + Lik(3)**2 )
  Lik(1:3) = Lik(1:3) / norm
  if (Lik(1)**2+Lik(2)**2 ==0.D0 ) cycle
  cpt_phi = cpt_phi + 1
end do

do i=1,Nsect
  do j=1,nb_ligneCONTACT
    cd  = TAB_CONTACT(j)%icdent
    an  = TAB_CONTACT(j)%ianent
    Lik(1:3) = TAB_SPHERE(cd)%centre(1:3)-TAB_SPHERE(an)%centre(1:3)
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
  
!  Lik(1:3) = TAB_SPHERE(cd)%centre(1:3)-TAB_SPHERE(an)%centre(1:3)
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
! spatial_correlation
!==============================================================================
subroutine spatial_correlation
implicit none
integer                                             ::  i,cd,an,nb_part_without_paroi,particule_i,nb_ksi,j
real(kind=8)                                        ::  rmoyen,P0,compteur_z_float,ksi,epsilon,dist,compteur_part,nb_part,&
                                                        radius_max,radius_mean
real(kind=8),dimension(:,:,:),allocatable           ::  Tab_correlation,Tab_correlation_max
real(kind=8),dimension(:,:),allocatable             ::  mean_correlation
real(kind=8),dimension(:),allocatable               ::  max_correlation
logical                                             ::  test_radius_max

print*,'.... spatial_correlation ...'

rmoyen = 0
P0 = 0
nb_part_without_paroi = 0
compteur_part  = 0
nb_part = 0
radius_max = 0

do i=1,nbParticules
  TAB_SPHERE(i)%nctc = 0
  TAB_SPHERE(i)%ctc_float = 0
  if (TAB_SPHERE(i)%behav=='PLEXx') nb_part_without_paroi = nb_part_without_paroi + 1.0
  rmoyen = rmoyen + TAB_SPHERE(i)%Rmax
end do
rmoyen = rmoyen/nbParticules
print*,'R_moyen= ',rmoyen

do i=1,nb_ligneCONTACT
  cd = TAB_CONTACT(i)%icdent
  an = TAB_CONTACT(i)%ianent
  if (TAB_CONTACT(i)%rn<0.000000000001*Mean_total_Normal_force) cycle
  if ((TAB_SPHERE(cd)%behav/='PLEXx').or.(TAB_SPHERE(an)%behav/='PLEXx')) cycle

  TAB_SPHERE(cd)%nctc = TAB_SPHERE(cd)%nctc + 1
  TAB_SPHERE(an)%nctc = TAB_SPHERE(an)%nctc + 1
end do

do i=1,nbParticules
  if ((TAB_SPHERE(i)%behav/='PLEXx')) cycle
  if (TAB_SPHERE(i)%nctc == 0 )  P0   = P0  + 1.0
  nb_part = nb_part + 1
end do

allocate(Tab_correlation(int(P0),10000,2))
allocate(Tab_correlation_max(int(P0),10000,2))
allocate(max_correlation(int(P0)))
P0 = P0 / real(nb_part,8)
Tab_correlation_max(:,:,:) = 0
max_correlation(:) = 0

print*, 'P(c=0) = ', P0

!ksi = 0.0000
!nb_ksi = 0
!epsilon = rmoyen/5
!do 
!  if (ksi/rmoyen<0.1) then
!    epsilon=rmoyen/10
!  else
!    epsilon = rmoyen/5
!  end if
!  particule_i = 0
!  nb_ksi = nb_ksi + 1
!  do i=1,nbParticules
!    if ((TAB_SPHERE(i)%behav/='PLEXx')) cycle
!    if (TAB_SPHERE(i)%nctc == 0 ) then
!      particule_i = particule_i + 1
!      compteur_z_float = 0
!      compteur_part = 0
!      do j=1,nbParticules
!        if ((TAB_SPHERE(j)%behav/='PLEXx')) cycle
!        dist = sqrt( (TAB_SPHERE(i)%centre(1)-TAB_SPHERE(j)%centre(1))**2 + &
!                     (TAB_SPHERE(i)%centre(2)-TAB_SPHERE(j)%centre(2))**2 + &
!                     (TAB_SPHERE(i)%centre(3)-TAB_SPHERE(j)%centre(3))**2 )

!        if ( TAB_SPHERE(i)%Rmax+ksi > dist - TAB_SPHERE(j)%Rmax) then
!          compteur_part = compteur_part + 1.
!          if (TAB_SPHERE(j)%nctc == 0) compteur_z_float = compteur_z_float + 1.
!        end if
!      end do
!      Tab_correlation(particule_i,nb_ksi,1) = (TAB_SPHERE(i)%Rmax+ksi)/(rmoyen)
!      Tab_correlation_max(particule_i,nb_ksi,1) = (TAB_SPHERE(i)%Rmax+ksi)/(rmoyen)

!      Tab_correlation(particule_i,nb_ksi,2) = compteur_z_float/compteur_part
!    end if
!  end do

!  ksi = ksi + epsilon
!  print*,ksi/rmoyen,nb_ksi
!  if (ksi>20*rmoyen) exit
!end do

particule_i = 0
do i=1, nbParticules
  if ((TAB_SPHERE(i)%behav/='PLEXx')) cycle
  if (TAB_SPHERE(i)%nctc == 0 ) then
    ksi         = TAB_SPHERE(i)%Rmax
    epsilon     = rmoyen / 5
    particule_i = particule_i + 1
    nb_ksi      = 1
    test_radius_max = .false.
    do 
      compteur_part    = 0.
      compteur_z_float = 0.
      do j=1,nbParticules
        if ((TAB_SPHERE(j)%behav/='PLEXx')) cycle
          dist = sqrt( (TAB_SPHERE(i)%centre(1)-TAB_SPHERE(j)%centre(1))**2 + &
                       (TAB_SPHERE(i)%centre(2)-TAB_SPHERE(j)%centre(2))**2 + &
                      (TAB_SPHERE(i)%centre(3)-TAB_SPHERE(j)%centre(3))**2 )
        if ( TAB_SPHERE(i)%Rmax+ksi > dist - TAB_SPHERE(j)%Rmax) then 
          compteur_part = compteur_part + 1.
          if (TAB_SPHERE(j)%nctc == 0) compteur_z_float = compteur_z_float + 1.
        end if
      end do
      
      Tab_correlation(particule_i,nb_ksi,1)     = (TAB_SPHERE(i)%Rmax+ksi)/ TAB_SPHERE(i)%Rmax
      Tab_correlation(particule_i,nb_ksi,2)     = compteur_z_float        / compteur_part
      
!      print*,'----'
!      print*,Tab_correlation(particule_i,nb_ksi,2)
      
      if (.not. test_radius_max) then
!        if ( (Tab_correlation(particule_i,nb_ksi,2)<0.50) ) then ! Critère a 0.5
        if ( (Tab_correlation(particule_i,nb_ksi,2)<0.20) ) then ! Critère a 0.2
          max_correlation(particule_i) = Tab_correlation(particule_i,nb_ksi,1)
          test_radius_max = .true.
          exit
        end if
      end if
      ksi = ksi + epsilon
      nb_ksi = nb_ksi + 1
      if (ksi>15*rmoyen) then
        max_correlation(particule_i) = (TAB_SPHERE(i)%Rmax+ksi)/ TAB_SPHERE(i)%Rmax
        print*, i
        exit
      end if
    end do
  end if

end do



allocate(mean_correlation(nb_ksi,2))
mean_correlation = 0
do i=1,nb_ksi
  do j=1,particule_i
    mean_correlation(i,1) = mean_correlation(i,1) +  Tab_correlation(j,i,1)
    mean_correlation(i,2) = mean_correlation(i,2) +  Tab_correlation(j,i,2)
  end do
end do
mean_correlation      = mean_correlation / real(particule_i,8)

do i=1,nb_ksi
  write(130,'(16(1X,D12.5))') mean_correlation(i,1), mean_correlation(i,2)
end do

radius_mean = 0
do i=1,particule_i
  radius_max = max(radius_max,max_correlation(i))
  radius_mean = radius_mean + max_correlation(i)
  write(131,'(16(1X,D12.5))') max_correlation(i)
end do

radius_mean = radius_mean/particule_i

write(132,'(16(1X,D12.5))') radius_max,radius_mean


stop

end subroutine spatial_correlation



!==============================================================================
! radial distribution
!==============================================================================
subroutine radial_distribution
implicit none
integer                                             ::  i,cd,an,nb_part_without_paroi,particule_i,nb_ksi,j
real(kind=8)                                        ::  rmoyen,P0,compteur_z_float,ksi,epsilon,dist,compteur_part,nb_part,&
                                                        radius_max,radius_mean,nb_total_float,nb_total_no_float,compteur_z_no_float
real(kind=8),dimension(:,:,:),allocatable           ::  Tab_correlation
real(kind=8),dimension(:,:),allocatable             ::  mean_correlation


print*,'.... radial_distribution ...'

rmoyen = 0
P0 = 0
nb_part_without_paroi = 0
compteur_part  = 0
nb_part = 0
radius_max = 0

do i=1,nbParticules
  TAB_SPHERE(i)%nctc = 0
  TAB_SPHERE(i)%ctc_float = 0
  if (TAB_SPHERE(i)%behav=='PLEXx') nb_part_without_paroi = nb_part_without_paroi + 1.0
  rmoyen = rmoyen + TAB_SPHERE(i)%Rmax
end do
rmoyen = rmoyen/nbParticules
print*,'R_moyen= ',rmoyen

do i=1,nb_ligneCONTACT
  cd = TAB_CONTACT(i)%icdent
  an = TAB_CONTACT(i)%ianent
  if (TAB_CONTACT(i)%rn<0.000000000001*Mean_total_Normal_force) cycle
  if ((TAB_SPHERE(cd)%behav/='PLEXx').or.(TAB_SPHERE(an)%behav/='PLEXx')) cycle

  TAB_SPHERE(cd)%nctc = TAB_SPHERE(cd)%nctc + 1
  TAB_SPHERE(an)%nctc = TAB_SPHERE(an)%nctc + 1
end do

do i=1,nbParticules
  if ((TAB_SPHERE(i)%behav/='PLEXx')) cycle
  if (TAB_SPHERE(i)%nctc == 0 )  P0   = P0  + 1.0
  nb_part = nb_part + 1
end do

allocate(Tab_correlation(int(P0),10000,3))

nb_total_float         = P0
nb_total_no_float      = nb_part-nb_total_float

print*, 'nb_float = ',  nb_total_float , 'nb_no_float = ', nb_total_no_float

P0                     = P0 / real(nb_part,8)
Tab_correlation(:,:,:) = 0

print*, 'P(c=0) = ', P0

particule_i = 0
do i=1, nbParticules
  if ((TAB_SPHERE(i)%behav/='PLEXx')) cycle
  if (TAB_SPHERE(i)%nctc == 0 ) then
!  if (TAB_SPHERE(i)%nctc > 0 ) then
  
    if ( (TAB_SPHERE(i)%centre(2) > 0.30+0.150).and.(TAB_SPHERE(i)%centre(2) < 0.30-0.150) ) cycle
    if ( (TAB_SPHERE(i)%centre(1) > 0.25+0.125).and.(TAB_SPHERE(i)%centre(1) < 0.25-0.125) ) cycle
    if ( (TAB_SPHERE(i)%centre(3) > 0.25+0.125).and.(TAB_SPHERE(i)%centre(3) < 0.25-0.125) ) cycle
  
    ksi         = rmoyen
    epsilon     = rmoyen / 10
    particule_i = particule_i + 1
    nb_ksi      = 1
    do 
      compteur_part       = 0.
      compteur_z_float    = 0.
      compteur_z_no_float = 0.
      do j=1,nbParticules
        if (i==j) cycle 
        if ((TAB_SPHERE(j)%behav/='PLEXx')) cycle
          dist = sqrt( (TAB_SPHERE(i)%centre(1)-TAB_SPHERE(j)%centre(1))**2 + &
                       (TAB_SPHERE(i)%centre(2)-TAB_SPHERE(j)%centre(2))**2 + &
                       (TAB_SPHERE(i)%centre(3)-TAB_SPHERE(j)%centre(3))**2 )
        if ( (dist < ksi+epsilon) .and. (dist > ksi) ) then 
          compteur_part = compteur_part + 1.
!          if ( (TAB_SPHERE(j)%nctc >  0).and.(TAB_SPHERE(i)%nctc >  0) ) compteur_z_no_float = compteur_z_no_float + 1.
          if ( (TAB_SPHERE(j)%nctc == 0).and.(TAB_SPHERE(i)%nctc == 0) ) compteur_z_float    = compteur_z_float + 1.
        end if
      end do
      
      Tab_correlation(particule_i,nb_ksi,1)     = (ksi)/ (2*rmoyen)
      Tab_correlation(particule_i,nb_ksi,2)     = compteur_z_float           /  (4*3.141592*epsilon*ksi**2) 
!      Tab_correlation(particule_i,nb_ksi,3)     = compteur_z_no_float        /  (4*3.141592*epsilon*ksi**2) 
                  
      ksi = ksi + epsilon
      nb_ksi = nb_ksi + 1
            
      if (ksi>16*rmoyen) exit
    end do
  end if
end do

print*, 'particule_i = ', particule_i

allocate(mean_correlation(nb_ksi,3))
mean_correlation = 0
do i=1,nb_ksi
  do j=1,particule_i
    mean_correlation(i,1) = mean_correlation(i,1) +  Tab_correlation(j,i,1)
    mean_correlation(i,2) = mean_correlation(i,2) +  Tab_correlation(j,i,2)
    mean_correlation(i,3) = mean_correlation(i,3) +  Tab_correlation(j,i,3)
  end do
end do
mean_correlation      = mean_correlation / real(particule_i,8)

do i=1,nb_ksi
  write(133,'(16(1X,D12.5))') mean_correlation(i,1), mean_correlation(i,2),&
                              mean_correlation(i,2)  /(real(nb_total_float,8)  / (Long*Larg*Haut) )
end do

stop

end subroutine radial_distribution


!==============================================================================
! correlation function
!==============================================================================
subroutine correlation_function
implicit none
integer                                             ::  i,cd,an,nb_part_without_paroi,particule_i,nb_ksi,j
real(kind=8)                                        ::  rmoyen,P0,compteur_z_float,ksi,epsilon,dist,compteur_part,nb_part,&
                                                        mean_density_i,mean_density_j,nb_part_r
real(kind=8),dimension(:,:),allocatable             ::  mean_correlation
real(kind=8),dimension(:),allocatable               ::  density_part
real(kind=8),dimension(:,:,:),allocatable           ::  Tab_correlation

print*,'.... correlation function ...'

rmoyen = 0
P0 = 0
compteur_part  = 0
nb_part = 0

do i=1,nbParticules
  TAB_SPHERE(i)%nctc = 0
  TAB_SPHERE(i)%ctc_float = 0
  if (TAB_SPHERE(i)%behav=='PLEXx') nb_part_without_paroi = nb_part_without_paroi + 1.0
  rmoyen = rmoyen + TAB_SPHERE(i)%Rmax
end do
rmoyen = rmoyen/nbParticules
print*,'R_moyen= ',rmoyen

do i=1,nb_ligneCONTACT
  cd = TAB_CONTACT(i)%icdent
  an = TAB_CONTACT(i)%ianent
  if (TAB_CONTACT(i)%rn<0.000000000001*Mean_total_Normal_force) cycle
  if ((TAB_SPHERE(cd)%behav/='PLEXx').or.(TAB_SPHERE(an)%behav/='PLEXx')) cycle

  TAB_SPHERE(cd)%nctc = TAB_SPHERE(cd)%nctc + 1
  TAB_SPHERE(an)%nctc = TAB_SPHERE(an)%nctc + 1
end do

do i=1,nbParticules
  if ((TAB_SPHERE(i)%behav/='PLEXx')) cycle
  if (TAB_SPHERE(i)%nctc == 0 )  P0   = P0  + 1.0
  nb_part = nb_part + 1
end do

allocate(Tab_correlation(int(P0),10000,2))
allocate(density_part(nbParticules))

P0 = P0 / real(nb_part,8)
density_part(:) = 0

!---- On calcul la densité de particules flotantes autour de chaque particules eloignées d'une distance 3d -----

!do i=1,nbParticules
!  if (TAB_SPHERE(i)%behav/='PLEXx') cycle
!  compteur_part    = 0
!  compteur_z_float = 0
!  do j=1,nbParticules
!    if (i==j) cycle
!    if ((TAB_SPHERE(j)%behav/='PLEXx')) cycle
    
!    dist = sqrt( (TAB_SPHERE(i)%centre(1)-TAB_SPHERE(j)%centre(1))**2 + &
!                 (TAB_SPHERE(i)%centre(2)-TAB_SPHERE(j)%centre(2))**2 + &
!                 (TAB_SPHERE(i)%centre(3)-TAB_SPHERE(j)%centre(3))**2 )

!    if ( dist < 3*(2*rmoyen) ) then
!      compteur_part = compteur_part + 1.
!      if (TAB_SPHERE(j)%nctc == 0) compteur_z_float = compteur_z_float + 1.
!    end if
!  end do
  
!  if (compteur_part == 0) then
!    density_part(i) = 0
!  else
!    density_part(i) = compteur_z_float / compteur_part  !( (4*3.141592*(3*rmoyen)**3 )/ 3 )  
!  end if
!end do

do i=1,nbParticules
  if (TAB_SPHERE(i)%behav/='PLEXx') cycle
  if (TAB_SPHERE(i)%nctc == 0) density_part(i) = 0
  if (TAB_SPHERE(i)%nctc >  0) density_part(i) = 1
end do

mean_density_i = 0
nb_part_r      = 0
do i=1,nbParticules
  if (TAB_SPHERE(i)%behav/='PLEXx') cycle
    mean_density_i = mean_density_i + density_part(i)
    nb_part_r      = nb_part_r + 1. 
end do
mean_density_i = mean_density_i / nb_part_r
mean_density_j = mean_density_i

print*,'mean density',mean_density_i, 'P(0)',P0

particule_i = 0
do i=1, nbParticules
  if ((TAB_SPHERE(i)%behav/='PLEXx')) cycle
    ksi                                    =  TAB_SPHERE(i)%Rmax
    epsilon                                =  rmoyen / 2
    particule_i                            =  particule_i + 1
    nb_ksi                                 =  1
    do 
      Tab_correlation(particule_i,nb_ksi,2)  =  0
      nb_part_r                              =  0
      do j=1,nbParticules
        if ((TAB_SPHERE(j)%behav/='PLEXx')) cycle
        if (i==j) cycle
        
          dist = sqrt( (TAB_SPHERE(i)%centre(1)-TAB_SPHERE(j)%centre(1))**2 + &
                       (TAB_SPHERE(i)%centre(2)-TAB_SPHERE(j)%centre(2))**2 + &
                      (TAB_SPHERE(i)%centre(3)-TAB_SPHERE(j)%centre(3))**2 )
        if ( (dist < ksi+epsilon) .and. (dist > ksi) ) then 
!        if ( (dist < ksi) ) then 
          Tab_correlation(particule_i,nb_ksi,2)     = Tab_correlation(particule_i,nb_ksi,2) + &
                                                      density_part(particule_i)*density_part(j)
          nb_part_r      = nb_part_r + 1. 
        end if
      end do

      if (nb_part_r==0) then
        nb_part_r = 1
        Tab_correlation(particule_i,nb_ksi,2) = 0
      end if 

      Tab_correlation(particule_i,nb_ksi,1)     =  ( TAB_SPHERE(i)%Rmax+ksi ) 
      Tab_correlation(particule_i,nb_ksi,2)     =  ( Tab_correlation(particule_i,nb_ksi,2) / nb_part_r) / &
                                                   ( mean_density_i*mean_density_j )

      print*, ksi/rmoyen, particule_i,nbParticules
      ksi = ksi + epsilon
      nb_ksi = nb_ksi + 1
      if (ksi>20*rmoyen) exit

    end do
end do



allocate(mean_correlation(nb_ksi,2))
mean_correlation = 0
do i=1,nb_ksi
  do j=1,particule_i
    mean_correlation(i,1) = mean_correlation(i,1) +  Tab_correlation(j,i,1)
    mean_correlation(i,2) = mean_correlation(i,2) +  Tab_correlation(j,i,2)
  end do
end do
mean_correlation      = mean_correlation / real(particule_i,8)

do i=1,nb_ksi
  write(133,'(16(1X,D12.5))') mean_correlation(i,1)/(2*rmoyen), mean_correlation(i,2)
end do




stop

end subroutine correlation_function




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
  Lik(1:3) = TAB_SPHERE(cd)%centre(1:3)-TAB_SPHERE(an)%centre(1:3)
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
      Lik(1:3)  =   TAB_SPHERE(cd)%centre(1:3)-TAB_SPHERE(an)%centre(1:3)
      norm      =   sqrt( Lik(1)**2 + Lik(2)**2 + Lik(3)**2 )
      Lik(1:3)  =   Lik(1:3) / norm
      
      phi       = acos( Lik(1) / sqrt(Lik(1)**2+Lik(2)**2) )
      if (Lik(1)**2+Lik(2)**2 == 0.D0 ) cycle
      
      
      
    end do                   !  Boucle sur les contacts
    
  end do        !  Boucle sur theta

end do          !  Boucle sur phi



end subroutine contact_orientation_3D




!==================================================================================================
!Calcul des deformation
!=================================================================================================
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
  do i=1,nbParticules
    if (TAB_SPHERE(i)%behav/='PLEXx') cycle
    
    Xplan_max = max( Xplan_max, TAB_SPHERE(i)%centre(1) + TAB_SPHERE(i)%Rmax )
    Xplan_min = min( Xplan_min, TAB_SPHERE(i)%centre(1) - TAB_SPHERE(i)%Rmax)

    Yplan_max = max( Yplan_max, TAB_SPHERE(i)%centre(2) + TAB_SPHERE(i)%Rmax )
    Yplan_min = min( Yplan_min, TAB_SPHERE(i)%centre(2) - TAB_SPHERE(i)%Rmax )

    Zplan_max = max( Zplan_max, TAB_SPHERE(i)%centre(3) + TAB_SPHERE(i)%Rmax )
    Zplan_min = min( Zplan_min, TAB_SPHERE(i)%centre(3) - TAB_SPHERE(i)%Rmax )
  end do
  Long_ini = Xplan_max-Xplan_min
  Larg_ini = Yplan_max-Yplan_min
  H_ini    = Zplan_max-Zplan_min

  Depsilon_p = 0.D0
  Depsilon_q = 0.D0

  taux_def = 0.047472569 / H_ini

  deformation_ini=.false.
end if

Xplan_max =-10000000
Xplan_min = 10000000
Yplan_max =-10000000
Yplan_min = 10000000
Zplan_max =-10000000
Zplan_min = 10000000
do i=1,nbParticules
  if (TAB_SPHERE(i)%behav/='PLEXx') cycle

  Xplan_max = max( Xplan_max, TAB_SPHERE(i)%centre(1) + TAB_SPHERE(i)%Rmax )
  Xplan_min = min( Xplan_min, TAB_SPHERE(i)%centre(1) - TAB_SPHERE(i)%Rmax)

  Yplan_max = max( Yplan_max, TAB_SPHERE(i)%centre(2) + TAB_SPHERE(i)%Rmax )
  Yplan_min = min( Yplan_min, TAB_SPHERE(i)%centre(2) - TAB_SPHERE(i)%Rmax )

  Zplan_max = max( Zplan_max, TAB_SPHERE(i)%centre(3) + TAB_SPHERE(i)%Rmax )
  Zplan_min = min( Zplan_min, TAB_SPHERE(i)%centre(3) - TAB_SPHERE(i)%Rmax )
end do

Long = Xplan_max-Xplan_min
Larg = Yplan_max-Yplan_min
Haut = Zplan_max -Zplan_min

E1 = log( Long / Long_ini )
E2 = log( Larg / Larg_ini )
E3 = log( Haut / H_ini    )

Depsilon_p = ( DE1 + DE2 + DE3 )
Depsilon_q = ( DE3 - DE1 ) 

Ep = E1+E2+E3
Eq = E3-E1

Def   = 0.047472569*temps
Def_H = abs(Haut - H_ini )
print*, 'Deformation = ',Def
print*, 'Long =' , Long, 'Larg =' , Larg , 'haut =' , Haut

end subroutine calcul_deformation


!==================================================================================================
!On ferme tout les fichiers
!=================================================================================================
subroutine fermer
implicit none
integer           :: i

print*,'---> fermeture des fichiers'

if (calcul_coordination                   == 1)    close (100)
if (calcul_vitesse_moyenne                == 1)    close (101)
if (calcul_qsurp                          == 1)    close (102)
if (calcul_compacity                      == 1)    close (103)
if (calcul_contact_anisotropies           == 1)    close (104)


deallocate(TAB_SPHERE)

print*,'--> Fin du postraitement <--'
stop

end subroutine fermer

!========================================================================================================!
!========================================================================================================!
!========================================================================================================!
!========================================================================================================!
!========================================================================================================!
!========================================================================================================!

!  ROUTINE GMV CA PUT DU CUL GRAVE, NE REGARDER QUE EN CAS DE NECESSITE ABSOLUE... FAUDRA VITE 
!  TROUVER UNE AUTRE SOLUTION

!========================================================================================================!
!========================================================================================================!
!========================================================================================================!
!========================================================================================================!
!========================================================================================================!
!========================================================================================================!


!============================================================================================================ 
!                      Subroutine write_gmv_cells 
! 
! Ecrit dans le fichier gmv repere par 'a' le nombre de cells necessaires a la 
! visualisation des 'n' particules 
!============================================================================================================ 
subroutine write_gmv_cells(a,n,cpt_force_in) 
 
  integer                     :: a,n,cpt_force_in
 
  write(a,'(A5,1X,I9)') 'cells',n*62 + 6*(cpt_force_in) ! 
 
end subroutine write_gmv_cells 


!============================================================================================================ 
!                      Subroutine write_gmv_x 
! 
! Ecrit dans le fichier gmv repere par 'a' la coordonnee x des 48 points d'une 
! particule spherique, centre en x sur 'coo1' et de rayon 'ray', servant a sa 
! representation 
!============================================================================================================ 
subroutine write_gmv_x(a,coo1,ray) 
 
  integer                     :: a,i 
  real(kind=8)                :: coo1,ray 
  real(kind=8),dimension(6,8) :: x_surf 
 
  do i=1,4 
    x_surf(1,i) = coo1+ray*cos(pi/12.D0) 
  enddo 
  do i=5,8 
    x_surf(1,i) = coo1+ray*cos(pi/4.D0) 
  enddo 
  do i=1,4 
    x_surf(2,i) = coo1+ray*cos(pi/4.D0) 
  enddo 
  do i=5,8 
    x_surf(2,i) = coo1+ray*cos(5.D0*pi/12.D0) 
  enddo 
  do i=1,8 
    x_surf(3,i) = coo1+ray*cos(5.D0*pi/12.D0) 
  enddo 
  do i=1,8 
    x_surf(4,i) = coo1+ray*cos(7.D0*pi/12.D0) 
  enddo 
  do i=1,4 
    x_surf(5,i) = coo1+ray*cos(7.D0*pi/12.D0) 
  enddo 
  do i=5,8 
    x_surf(5,i) = coo1+ray*cos(3.D0*pi/4.D0) 
  enddo 
  do i=1,4 
    x_surf(6,i) = coo1+ray*cos(3.D0*pi/4.D0) 
  enddo 
  do i=5,8 
    x_surf(6,i) = coo1+ray*cos(11.D0*pi/12.D0) 
  enddo 
 
  do i=1,6 
    write(a,'(8(1X,E14.7))') x_surf(i,1),x_surf(i,2),x_surf(i,3),x_surf(i,4), & 
                             x_surf(i,5),x_surf(i,6),x_surf(i,7),x_surf(i,8) 
  enddo 
 
end subroutine write_gmv_x 

!============================================================================================================ 
!                      Subroutine write_gmv_y 
! 
! Ecrit dans le fichier gmv repere par 'a' la coordonnee y des 48 points d'une 
! particule spherique, centre en y sur 'coo2' et de rayon 'ray', servant a sa 
! representation 
!============================================================================================================ 
subroutine write_gmv_y(a,coo2,ray) 
 
  integer                     :: a,i 
  real(kind=8)                :: coo2,ray 
  real(kind=8),dimension(3)   :: r_surf 
  real(kind=8),dimension(6,8) :: y_surf 
 
  r_surf(1) = ray*sin(pi/12.D0) 
  r_surf(2) = ray*sin(pi/4.D0) 
  r_surf(3) = ray*sin(5.D0*pi/12.D0) 
 
  do i=1,4 
    y_surf(1,i) = coo2+r_surf(1)*cos((1.D0+2.D0*(i-1.D0))*pi/4.D0) 
    y_surf(6,i+4) = coo2+r_surf(1)*cos((1.D0+2.D0*(i-1.D0))*pi/4.D0) 
  enddo 
  do i=1,4 
    y_surf(1,i+4) = coo2+r_surf(2)*cos((1.D0+2.D0*(i-1.D0))*pi/8.D0) 
    y_surf(6,i) = coo2+r_surf(2)*cos((1.D0+2.D0*(i+3.D0))*pi/8.D0) 
  enddo 
  do i=1,4 
    y_surf(2,i) = coo2+r_surf(2)*cos((1.D0+2.D0*(i+3.D0))*pi/8.D0) 
    y_surf(5,i+4) = coo2+r_surf(2)*cos((1.D0+2.D0*(i-1.D0))*pi/8.D0) 
  enddo 
  do i=1,4 
    y_surf(2,i+4) = coo2+r_surf(3)*cos((1.D0+2.D0*(i-1.D0))*pi/12.D0) 
    y_surf(5,i) = coo2+r_surf(3)*cos((1.D0+2.D0*(i+7.D0))*pi/12.D0) 
  enddo 
  do i=1,8 
    y_surf(3,i) = coo2+r_surf(3)*cos((1.D0+2.D0*(i+3.D0))*pi/12.D0) 
    y_surf(4,i) = coo2+r_surf(3)*cos((1.D0+2.D0*(i-1.D0))*pi/12.D0) 
  enddo 
 
  do i=1,6 
    write(a,'(8(1X,E14.7))') y_surf(i,1),y_surf(i,2),y_surf(i,3),y_surf(i,4), & 
                             y_surf(i,5),y_surf(i,6),y_surf(i,7),y_surf(i,8) 
  enddo 
 
end subroutine write_gmv_y 

!============================================================================================================ 
!                      Subroutine write_gmv_z 
! 
! Ecrit dans le fichier gmv repere par 'a' la coordonnee z des 48 points d'une 
! particule spherique, centre en z sur 'coo3' et de rayon 'ray', servant a sa 
! representation 
!============================================================================================================ 
subroutine write_gmv_z(a,coo3,ray) 
 
  integer                     :: a,i 
  real(kind=8)                :: coo3,ray 
  real(kind=8),dimension(3)   :: r_surf 
  real(kind=8),dimension(6,8) :: z_surf 
 
  r_surf(1) = ray*sin(pi/12.D0) 
  r_surf(2) = ray*sin(pi/4.D0) 
  r_surf(3) = ray*sin(5.D0*pi/12.D0) 
 
  do i=1,4 
    z_surf(1,i) = coo3+r_surf(1)*sin((1.D0+2.D0*(i-1.D0))*pi/4.D0) 
    z_surf(6,i+4) = coo3+r_surf(1)*sin((1.D0+2.D0*(i-1.D0))*pi/4.D0) 
  enddo 
  do i=1,4 
    z_surf(1,i+4) = coo3+r_surf(2)*sin((1.D0+2.D0*(i-1.D0))*pi/8.D0) 
    z_surf(6,i) = coo3+r_surf(2)*sin((1.D0+2.D0*(i+3.D0))*pi/8.D0) 
  enddo 
  do i=1,4 
    z_surf(2,i) = coo3+r_surf(2)*sin((1.D0+2.D0*(i+3.D0))*pi/8.D0) 
    z_surf(5,i+4) = coo3+r_surf(2)*sin((1.D0+2.D0*(i-1.D0))*pi/8.D0) 
  enddo 
  do i=1,4 
    z_surf(2,i+4) = coo3+r_surf(3)*sin((1.D0+2.D0*(i-1.D0))*pi/12.D0) 
    z_surf(5,i) = coo3+r_surf(3)*sin((1.D0+2.D0*(i+7.D0))*pi/12.D0) 
  enddo 
  do i=1,8 
    z_surf(3,i) = coo3+r_surf(3)*sin((1.D0+2.D0*(i+3.D0))*pi/12.D0) 
    z_surf(4,i) = coo3+r_surf(3)*sin((1.D0+2.D0*(i-1.D0))*pi/12.D0) 
  enddo 
 
  do i=1,6 
    write(a,'(8(1X,E14.7))') z_surf(i,1),z_surf(i,2),z_surf(i,3),z_surf(i,4), & 
                             z_surf(i,5),z_surf(i,6),z_surf(i,7),z_surf(i,8) 
  enddo 
 
end subroutine write_gmv_z 

!============================================================================================================ 
!                      Subroutine write_gmv_surf 
! 
! Ecrit dans le fichier gmv repere par 'a' les definitions de toutes les 
! surfaces (tri et qua) representant les 'n' spheres 
! Une sphere = 62 faces
!============================================================================================================ 
subroutine write_gmv_surf(a,n) 
 
  integer                     :: a,n,i,j,k,cd,an
 
  do i=1,n 
    j = 48*(i-1) 
    write(a,'(A3,1X,I7)') 'tri',3 
    write(a,'(3(I11))') 1+j,5+j,6+j 
    write(a,'(A3,1X,I7)') 'tri',3 
    write(a,'(3(I11))') 45+j,37+j,38+j 
    write(a,'(A3,1X,I7)') 'tri',3 
    write(a,'(3(I11))') 6+j,5+j,14+j 
    write(a,'(A3,1X,I7)') 'tri',3 
    write(a,'(3(I11))') 38+j,37+j,26+j 
    write(a,'(A3,1X,I7)') 'tri',3 
    write(a,'(3(I11))') 5+j,13+j,14+j 
    write(a,'(A3,1X,I7)') 'tri',3 
    write(a,'(3(I11))') 6+j,14+j,15+j 
    write(a,'(A3,1X,I7)') 'tri',3 
    write(a,'(3(I11))') 37+j,25+j,26+j 
    write(a,'(A3,1X,I7)') 'tri',3 
    write(a,'(3(I11))') 38+j,26+j,27+j 
    write(a,'(A3,1X,I7)') 'tri',3 
    write(a,'(3(I11))') 2+j,7+j,8+j 
    write(a,'(A3,1X,I7)') 'tri',3 
    write(a,'(3(I11))') 46+j,39+j,40+j 
    write(a,'(A3,1X,I7)') 'tri',3 
    write(a,'(3(I11))') 8+j,7+j,17+j 
    write(a,'(A3,1X,I7)') 'tri',3 
    write(a,'(3(I11))') 40+j,39+j,29+j 
    write(a,'(A3,1X,I7)') 'tri',3 
    write(a,'(3(I11))') 7+j,16+j,17+j 
    write(a,'(A3,1X,I7)') 'tri',3 
    write(a,'(3(I11))') 8+j,17+j,18+j 
    write(a,'(A3,1X,I7)') 'tri',3 
    write(a,'(3(I11))') 39+j,28+j,29+j 
    write(a,'(A3,1X,I7)') 'tri',3 
    write(a,'(3(I11))') 40+j,29+j,30+j 
    write(a,'(A3,1X,I7)') 'tri',3 
    write(a,'(3(I11))') 3+j,9+j,10+j 
    write(a,'(A3,1X,I7)') 'tri',3 
    write(a,'(3(I11))') 47+j,41+j,42+j 
    write(a,'(A3,1X,I7)') 'tri',3 
    write(a,'(3(I11))') 10+j,9+j,20+j 
    write(a,'(A3,1X,I7)') 'tri',3 
    write(a,'(3(I11))') 42+j,41+j,32+j 
    write(a,'(A3,1X,I7)') 'tri',3 
    write(a,'(3(I11))') 9+j,19+j,20+j 
    write(a,'(A3,1X,I7)') 'tri',3 
    write(a,'(3(I11))') 10+j,20+j,21+j 
    write(a,'(A3,1X,I7)') 'tri',3 
    write(a,'(3(I11))') 41+j,31+j,32+j 
    write(a,'(A3,1X,I7)') 'tri',3 
    write(a,'(3(I11))') 42+j,32+j,33+j 
    write(a,'(A3,1X,I7)') 'tri',3 
    write(a,'(3(I11))') 4+j,11+j,12+j 
    write(a,'(A3,1X,I7)') 'tri',3 
    write(a,'(3(I11))') 48+j,43+j,44+j 
    write(a,'(A3,1X,I7)') 'tri',3 
    write(a,'(3(I11))') 12+j,11+j,23+j 
    write(a,'(A3,1X,I7)') 'tri',3 
    write(a,'(3(I11))') 44+j,43+j,35+j 
    write(a,'(A3,1X,I7)') 'tri',3 
    write(a,'(3(I11))') 11+j,22+j,23+j 
    write(a,'(A3,1X,I7)') 'tri',3 
    write(a,'(3(I11))') 12+j,23+j,24+j 
    write(a,'(A3,1X,I7)') 'tri',3 
    write(a,'(3(I11))') 43+j,34+j,35+j 
    write(a,'(A3,1X,I7)') 'tri',3 
    write(a,'(3(I11))') 44+j,35+j,36+j 
    write(a,'(A4,1X,I6)') 'quad',4 
    write(a,'(4(I11))') 1+j,2+j,3+j,4+j 
    write(a,'(A4,1X,I6)') 'quad',4 
    write(a,'(4(I11))') 45+j,46+j,47+j,48+j 
    write(a,'(A4,1X,I6)') 'quad',4 
    write(a,'(4(I11))') 1+j,6+j,7+j,2+j 
    write(a,'(A4,1X,I6)') 'quad',4 
    write(a,'(4(I11))') 45+j,38+j,39+j,46+j 
    write(a,'(A4,1X,I6)') 'quad',4 
    write(a,'(4(I11))') 6+j,15+j,16+j,7+j 
    write(a,'(A4,1X,I6)') 'quad',4 
    write(a,'(4(I11))') 38+j,27+j,28+j,39+j 
    write(a,'(A4,1X,I6)') 'quad',4 
    write(a,'(4(I11))') 3+j,10+j,11+j,4+j 
    write(a,'(A4,1X,I6)') 'quad',4 
    write(a,'(4(I11))') 47+j,42+j,43+j,48+j 
    write(a,'(A4,1X,I6)') 'quad',4 
    write(a,'(4(I11))') 10+j,21+j,22+j,11+j 
    write(a,'(A4,1X,I6)') 'quad',4 
    write(a,'(4(I11))') 42+j,33+j,34+j,43+j 
    write(a,'(A4,1X,I6)') 'quad',4 
    write(a,'(4(I11))') 1+j,4+j,12+j,5+j 
    write(a,'(A4,1X,I6)') 'quad',4 
    write(a,'(4(I11))') 45+j,48+j,44+j,37+j 
    write(a,'(A4,1X,I6)') 'quad',4 
    write(a,'(4(I11))') 12+j,24+j,13+j,5+j 
    write(a,'(A4,1X,I6)') 'quad',4 
    write(a,'(4(I11))') 37+j,44+j,36+j,25+j 
    write(a,'(A4,1X,I6)') 'quad',4 
    write(a,'(4(I11))') 2+j,8+j,9+j,3+j 
    write(a,'(A4,1X,I6)') 'quad',4 
    write(a,'(4(I11))') 46+j,40+j,41+j,47+j 
    write(a,'(A4,1X,I6)') 'quad',4 
    write(a,'(4(I11))') 8+j,18+j,19+j,9+j 
    write(a,'(A4,1X,I6)') 'quad',4 
    write(a,'(4(I11))') 40+j,30+j,31+j,41+j 
    write(a,'(A4,1X,I6)') 'quad',4 
    write(a,'(4(I11))') 13+j,24+j,36+j,25+j 
    write(a,'(A4,1X,I6)') 'quad',4 
    write(a,'(4(I11))') 14+j,13+j,25+j,26+j 
    write(a,'(A4,1X,I6)') 'quad',4 
    write(a,'(4(I11))') 15+j,14+j,26+j,27+j 
    write(a,'(A4,1X,I6)') 'quad',4 
    write(a,'(4(I11))') 16+j,15+j,27+j,28+j 
    write(a,'(A4,1X,I6)') 'quad',4 
    write(a,'(4(I11))') 17+j,16+j,28+j,29+j 
    write(a,'(A4,1X,I6)') 'quad',4 
    write(a,'(4(I11))') 18+j,17+j,29+j,30+j 
    write(a,'(A4,1X,I6)') 'quad',4 
    write(a,'(4(I11))') 19+j,18+j,30+j,31+j 
    write(a,'(A4,1X,I6)') 'quad',4 
    write(a,'(4(I11))') 20+j,19+j,31+j,32+j 
    write(a,'(A4,1X,I6)') 'quad',4 
    write(a,'(4(I11))') 21+j,20+j,32+j,33+j 
    write(a,'(A4,1X,I6)') 'quad',4 
    write(a,'(4(I11))') 22+j,21+j,33+j,34+j 
    write(a,'(A4,1X,I6)') 'quad',4 
    write(a,'(4(I11))') 23+j,22+j,34+j,35+j 
    write(a,'(A4,1X,I6)') 'quad',4 
    write(a,'(4(I11))') 24+j,23+j,35+j,36+j 
  enddo 
  
  ! FORCES AUX BRANCHES-----
  k = j + 48
  k = 48*n
  do i=1,nb_ligneCONTACT
    if  (TAB_CONTACT(i)%nature == 'SPPLx') cycle 
    if  (TAB_CONTACT(i)%cycle_gmv) cycle
    if  (TAB_CONTACT(i)%rn<0.000000000001*Mean_total_Normal_force) cycle
    cd   = TAB_CONTACT(i)%icdent
    an   = TAB_CONTACT(i)%ianent
    write(30000,'(A4,1X,I6)') 'quad',4 
    write(30000,'(4(I11))')  1+k,2+k,3+k,4+K
    write(30000,'(A4,1X,I6)') 'quad',4 
    write(30000,'(4(I11))')  2+k,6+k,7+k,3+K
    write(30000,'(A4,1X,I6)') 'quad',4 
    write(30000,'(4(I11))')  3+k,7+k,8+k,4+K
    write(30000,'(A4,1X,I6)') 'quad',4 
    write(30000,'(4(I11))')  4+k,8+k,5+k,1+K
    write(30000,'(A4,1X,I6)') 'quad',4 
    write(30000,'(4(I11))')  1+k,5+k,6+k,2+K
    write(30000,'(A4,1X,I6)') 'quad',4 
    write(30000,'(4(I11))')  8+k,7+k,6+k,5+K
    k = k + 8
  end do

end subroutine write_gmv_surf 


!============================================================================================================ 
!                      Subroutine write_gmv_mat 
! 
! Ecrit dans le fichier gmv repere par 'a' l'indice materiel 'i_mat' affecte 
! aux 'm' particules 
!============================================================================================================ 
subroutine write_gmv_mat(a,m,cpt_force_ini,cpt_force_N_ini) 
 
  integer                     :: a,m,i,j,cpt_force_ini,cpt_force_N_ini,cd,an,k
  real(kind=8),dimension(3)                   ::  nik,tik,sik,Lik,Fik
  real(kind=8)                                ::  Rtik,Rnik,Rsik,LN,force_contact
  logical                                     ::  mobilized=.false.
  
!------------------------------------------------------------------------------ 
! affectation d'un label materiau pour gmv 
!------------------------------------------------------------------------------ 
  do i=1,nbParticules
    TAB_SPHERE(i)%nctc = 0
  end do

  do i=1,nb_ligneCONTACT
    cd = TAB_CONTACT(i)%icdent
    an = TAB_CONTACT(i)%ianent
    if (TAB_CONTACT(i)%rn<0.000000000001*Mean_total_Normal_force) cycle
    if ((TAB_SPHERE(cd)%behav/='PLEXx').or.(TAB_SPHERE(an)%behav/='PLEXx')) cycle

    TAB_SPHERE(cd)%nctc = TAB_SPHERE(cd)%nctc + 1
    TAB_SPHERE(an)%nctc = TAB_SPHERE(an)%nctc + 1
  end do


  write(a,'(A8)')     'variable'
  
  call write_gmv_V(a)

  write(a,'(A7)') 'endvars'

  write(a,'(A8,1X,I4,1X,I4)') 'material',4,0 
  write(a,'(A8)') 'CTC_0   '
  write(a,'(A8)') 'CTC_UP_0'
  write(a,'(A8)') 'WALLx   '
  write(a,'(A8)') 'FORCES  ' 
  
  do i=1,nbParticules
    if (TAB_SPHERE(i)%behav=='PLEXx') then
      if (TAB_SPHERE(i)%nctc==0)  write(a,'(8(I6))') ( 1,j=1,62 )
      if (TAB_SPHERE(i)%nctc>0)   write(a,'(8(I6))') ( 2,j=1,62 )
    else
      write(a,'(8(I6))') ( 3,j=1,62 )
    end if
  end do
  
  do i=1,nb_ligneCONTACT
    cd = TAB_CONTACT(i)%icdent
    an = TAB_CONTACT(i)%ianent
    if (TAB_CONTACT(i)%rn<0.000000000001*Mean_total_Normal_force) cycle
    if ((TAB_SPHERE(cd)%behav/='PLEXx').or.(TAB_SPHERE(an)%behav/='PLEXx')) cycle

    Lik(1:3) = TAB_SPHERE(cd)%centre(1:3)-TAB_SPHERE(an)%centre(1:3)
    LN       =  sqrt(Lik(1)**2+ Lik(2)**2 + Lik(3)**2)  
    Fik(1:3) = 0.D0
    mobilized = .false.

    nik  = TAB_CONTACT(i)%n
    tik  = TAB_CONTACT(i)%t
    sik  = TAB_CONTACT(i)%s
    Rtik = TAB_CONTACT(i)%rt
    Rnik = TAB_CONTACT(i)%rn
    Rsik = TAB_CONTACT(i)%rs

    Fik(1:3) = Rnik

    write(a,'(6(I6))') (4,j=1,6)

  end do

  write(a,'(A8)') 'polygons' 
  write(a,'(A7)') 'endpoly' 
  write(a,'(A6)') 'endgmv' 
end subroutine write_gmv_mat 


subroutine write_gmv_V(a)
implicit none
integer                     :: a,i,j
write(a,'(A5,4X,I2)')  '<Vx>  ',1
do i=1,nbParticules
  write(a,'(8(E14.7,1X))') ( TAB_SPHERE(i)%Vx,j=1,48)
end do
write(a,'(A5,4X,I2)')  '<Vy>  ',1
do i=1,nbParticules
  write(a,'(8(E14.7,1X))') ( TAB_SPHERE(i)%Vy,j=1,48)
end do
write(a,'(A5,4X,I2)')  '<Vz>  ',1
do i=1,nbParticules
  write(a,'(8(E14.7,1X))') ( TAB_SPHERE(i)%Vz,j=1,48)
end do
end subroutine write_gmv_V


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
  do i=1,nbParticules
    dmean = dmean + 2*TAB_SPHERE(i)%Rmax
    H_max = max(H_max,TAB_SPHERE(i)%centre(3))
  end do
  dmean = dmean / real(nbParticules,8)
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

  do i=1,nbParticules
    TAB_SPHERE(i)%behav = 'PLEXx'
    if ( TAB_SPHERE(i)%centre(3) < 2*dmean ) TAB_SPHERE(i)%behav='wallb'
    if ( TAB_SPHERE(i)%centre(3) > H_max - 2*dmean ) TAB_SPHERE(i)%behav='wallh'
  end do

  cpt = 0
  do i=1,nbParticules
    if (TAB_SPHERE(i)%behav /= 'PLEXx')  cycle
    cpt = cpt + 1
    
    write(128,'(A6)')                                        '$bdyty'
    write(128,'(1X,A5,1X,i6)')                               'RBDY3',cpt
    write(128,'(A6)')                                        '$blmty'
    write(128,'(1X,A5,5X,i2,2X,A5,2X,A5,1X,A6,D14.7)')       'PLAIN',1,'behav','PLEXx','avrd=',0.D0
    write(128,'(29X,3(A5,D14.7,2X))')  'I1  =',TAB_SPHERE(i)%I1,'I2  =',TAB_SPHERE(i)%I2,'I3  =',TAB_SPHERE(i)%I3
    write(128,'(A6)')                                        '$nodty'
    write(128,'(1X,A5,5X,i2,16X,3(A5,D14.7,2X))')            'NO6xx',1,'coo1=',TAB_SPHERE(i)%centre(1),&
                                                                       'coo2=',TAB_SPHERE(i)%centre(2),&
                                                                       'coo3=',TAB_SPHERE(i)%centre(3)
    write(128,'(29X,3(A5,D14.7,2X))') 'coo4=',0.D0,'coo5=',0.D0,'coo6=',0.D0
    write(128,'(A6)')                                        '$tacty'
    write(128,'(1X,A5,5X,i2,2X,A5,2X,A5,2X,A5,D14.7)')'SPHER',1,'color',TAB_SPHERE(i)%color,'byrd=',TAB_SPHERE(i)%Rmax
    write(128,'(A6)')                                        '$$$$$$'
  end do

  nbPart = cpt
  nb_wallh = 0
  nb_wallb = 0
  
  do i=1,nbParticules
    if (TAB_SPHERE(i)%behav /= 'wallh')  cycle
    cpt = cpt + 1
    nb_wallh = nb_wallh + 1
    
    write(128,'(A6)')                                        '$bdyty'
    write(128,'(1X,A5,1X,i6)')                               'RBDY3',cpt
    write(128,'(A6)')                                        '$blmty'
    write(128,'(1X,A5,5X,i2,2X,A5,2X,A5,1X,A6,D14.7)')       'PLAIN',1,'behav',TAB_SPHERE(i)%behav,'avrd=',0.D0
    write(128,'(29X,3(A5,D14.7,2X))')  'I1  =',TAB_SPHERE(i)%I1,'I2  =',TAB_SPHERE(i)%I2,'I3  =',TAB_SPHERE(i)%I3
    write(128,'(A6)')                                        '$nodty'
    write(128,'(1X,A5,5X,i2,16X,3(A5,D14.7,2X))')            'NO6xx',1,'coo1=',TAB_SPHERE(i)%centre(1),&
                                                                       'coo2=',TAB_SPHERE(i)%centre(2),&
                                                                       'coo3=',TAB_SPHERE(i)%centre(3)
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
  do i=1,nbParticules
    if (TAB_SPHERE(i)%behav /= 'wallb')  cycle  
    nb_wallb = nb_wallb + 1
    write(128,'(1X,A5,5X,i2,2X,A5,2X,A5,2X,A5,D14.7)')'SPHEb',1,'color','TATAb','byrd=',TAB_SPHERE(i)%Rmax
    write(128,'(29X,3(A5,D14.7,2X))') 'coo1=',TAB_SPHERE(i)%centre(1)-0.5/2,&
                                      'coo2=',TAB_SPHERE(i)%centre(2)-0.6/2,&
                                      'coo3=',TAB_SPHERE(i)%centre(3)-0.D0
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

end program grandeur3d
