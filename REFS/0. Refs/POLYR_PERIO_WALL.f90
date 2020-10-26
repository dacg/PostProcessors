program grandeur3d

implicit none

type  ::  T_CORPS
  integer                                    ::  nb_vertex,nb_face,num,gmv_material,nb_ctc
  real(kind=8),dimension(:,:),pointer        ::  vertex,vertex_ref
  integer,dimension(:,:),pointer             ::  face
  real(kind=8)                               ::  Rmax,volume,Ax1,Ax2,Ax3,I1,I2,I3,Sphericity,Aire
  real(kind=8),dimension(3)                  ::  centre,centre_ref
  real(kind=8)                               ::  Vx,Vy,Vz,Vrx,Vry,Vrz
  real(kind=8)                               ::  Wx,Wy,Wz
  real(kind=8)                               ::  Norme_V,Norme_Vr,Force
  real(kind=8),dimension(3,3)                ::  momentI,Rot
  character(len=5)                           ::  behav,color  
  real(kind=8)                               ::  nctc,nctcslide,nctcstick,nctcs,nctcd,nctct,Pressure
end type T_CORPS

type  :: T_FACE
  integer,dimension(:,:),pointer              ::  face   
  integer                                     ::  nb_face
end type

type  ::  T_CONTACT                                                   
  integer                                   ::  icdent,ianent,type        
  real(kind=8),dimension(3)                 ::  n,t,s                 
  real(kind=8),dimension(3)                 ::  coor_ctc              
  real(kind=8)                              ::  rn,rt,rs,vls,vln,vlt
  integer                                   ::  sect
  character(len=5)                          ::  nature,statut  
  logical                                   ::  deja_compte,cycle_gmv
end type T_CONTACT


type(T_CORPS),dimension(:),pointer            ::  TAB_POLYEDRE,TAB_PLAN
type(T_FACE),dimension(:),pointer             ::  TAB_DES_FACE
type(T_CONTACT),dimension(:),allocatable      ::  TAB_CONTACT,TAB_CONTACT_POLYR

!DEFINITIONS Vloc-Dof-Bodies
!================================================================
integer                                     ::  nb_paroi=0,fin=0,nbParticules=0,numero_paroi=0
integer                                     ::  compteur_clout=0,all_pas=0
logical                                     ::  the_first_time=.true., fin_post=.false. 
integer                                     ::  pas,nb_ligneCONTACT,nb_ligneCONTACT_POLYR,step_all_gmv,step_all_vtk,&
                                                step_pdf_start,step_pdf_stop,compteur_vtk=0
real(kind=8)                                ::  temps,Mean_total_Normal_force=0
 
!DEFINITIONS BOITE
!================================================================
real(kind=8),dimension(3,8)                 ::  tab_sommet_box
real(kind=8)                                ::  Def,H_ini,Long_ini,Larg_ini,&
                                                Depsilon_p,Depsilon_q,&
                                                Haut_avant,Long_avant,Larg_avant,Ep_avant=0,Eq_avant=0,Ep,Eq,&
                                                pressure_in_sample
logical                                     ::  deformation_ini = .true.

!CALCUL DES PROFILS 
!================================================================
integer                                     ::  debut_profils=0
logical                                     ::  the_first_time_profils = .true.

!CALCUL DU PROFIL DE VITESSE MOYEN
!================================================================
type  ::  T_step_velocity
  real(kind=8),dimension(:,:),allocatable             ::  tab_vmoyen
  real(kind=8),dimension(:,:),allocatable             ::  tab_vrmoyen
end type T_step_velocity
type(T_step_velocity),dimension(:),allocatable      ::  TAB_vitesse

!CALCUL DU PROFIL DE CONTRAINTE MOYEN
!================================================================
type  ::  T_sigma
  real(kind=8),dimension(3,3)                       :: sigma
end type
type  ::  T_step_contrainte
  type(T_sigma),dimension(:),allocatable            ::  tab_sigma
end type 
type(T_step_contrainte),dimension(:),allocatable      ::  TAB_contrainte


!MOYENNE POUR KSI-RESEAU
!================================================================
real(kind=8),dimension(:,:,:),allocatable   ::  Mean_Ksi_reseau 
integer                                     ::  pas_ksi = 0,step_ksi_reseau,stop_ksi_reseau,quel_pas_de_temps


!COMMANDES D APPELLE
!================================================================
character(len=30)                           ::  command
real(kind=8)                                ::  precision, PI = 3.14159265359,Longueur,Largeur,Hauteur

integer                                     ::  type_box=0,calcul_vitesse_moyenne=0,calcul_coordination=0,&
                                                calcul_qsurp=0,calcul_compacity=0,calcul_write_all_gmv=0,&
                                                contact_anisotropes=0,calcul_write_all_vtk=0,branche_anisotropes=0,&
                                                calcul_pdf=0,calcul_ksi_reseau=0,calcul_construct_bodies=0,&
                                                calcul_profils=0
                                                     
integer                                     ::  calcul_contact_orientation=0

                                                
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

  if (command=='VLOC_DOF_BODIES              :') then
    read(1,*) nb_paroi
    read(1,*) fin
    read(1,*) nbParticules
    read(1,*) all_pas
    read(1,*) numero_paroi
    open (unit=10,file='DEFORMATION.DAT',status='replace')     
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
    cycle
  end if

  if (command=='COMPACITY                    :') then
    calcul_compacity=1
    open (unit=103,file='COMPACITY.DAT',status='replace')     
    cycle
  end if

  if (command=='CONTACT ANISOTROPIE          :') then
    contact_anisotropes=1
    open (unit=104,file='CONTACT_ANISOTROPIES.DAT',status='replace')     
    cycle
  end if

  if (command=='BRANCHE ANISOTROPIE          :') then
    branche_anisotropes=1
    open (unit=105,file='BRANCHE_ANISOTROPIES.DAT',status='replace')     
    cycle
  end if

  if (command=='PDF                          :') then
    calcul_pdf=1
    read(1,*) step_pdf_start
    read(1,*) step_pdf_stop
    open (unit=106,file='PDF.DAT',status='replace')     
    cycle
  end if
  
  if (command=='WRITE ALL GMV                :') then
    calcul_write_all_gmv=1
    read(1,*) step_all_gmv
    cycle
  end if

  if (command=='WRITE ALL VTK                :') then
    calcul_write_all_vtk=1
    read(1,*) step_all_vtk
    cycle
  end if

  if (command=='KSI RESEAU                   :') then
    calcul_ksi_reseau=1
    read(1,*) step_ksi_reseau
    read(1,*) stop_ksi_reseau
    allocate( Mean_Ksi_reseau( int(stop_ksi_reseau-step_ksi_reseau)+1,1000,15 ) )
    open (unit=107,file='KSI_RESEAU.DAT',status='replace')     
    open (unit=108,file='KSI_RESEAU_Prop_Contact.DAT',status='replace')     
    cycle
  end if
  

  if (command=='CONTACT ORIENTATION          :') then
    calcul_contact_orientation=1
    open (unit=126,file='ORIENTATION_CONTACT_xy.DAT',status='replace') 
    open (unit=124,file='ORIENTATION_CONTACT_THETA.DAT',status='replace') 
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
    read(1,*) debut_profils
    open (unit=131,file='PROFIL_VITESSE_MOYEN.DAT',status='replace')     
    open (unit=132,file='PROFIL_CONTRAINTE_MOYEN.DAT',status='replace') 
    open (unit=133,file='PROFIL_I_MOYEN.DAT',status='replace')     
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

    if (calcul_qsurp                          == 1)  call qsurp    
    if (calcul_coordination                   == 1)  call nb_coordination
    if (calcul_vitesse_moyenne                == 1)  call vitesse_moyenne
    if (calcul_compacity                      == 1)  call compacity
    if (contact_anisotropes                   == 1)  call contact_anisotropies
    if (branche_anisotropes                   == 1)  call branche_anisotropies
    if (calcul_write_all_gmv                  == 1)  call write_all_gmv
    if (calcul_write_all_vtk                  == 1)  call write_all_vtk
    if (calcul_pdf                            == 1)  call pdf
    if (calcul_ksi_reseau                     == 1)  then
      pas_ksi = pas_ksi + 1
      call ksi_reseau 
    end if
    if (calcul_construct_bodies               == 1)  call construct_bodies
    if (calcul_profils                        == 1)  call profils

    if (calcul_contact_orientation            == 1)  call contact_orientation

  end if

end do


!==============================================================================
!==============================================================================
!==============================================================================
!==============================================================================

contains

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
type(T_moyenne),dimension(:),allocatable   :: Moyenne_contrainte_moyenne
real(kind=8)                               ::  Rtik,Rnik,Rsik
real(kind=8),dimension(3)                  ::  nik,sik,tik,Lik,Fik
integer                                    :: ierror_sigma,matz_sigma,lda_sigma
integer                                    :: cd,an
real(kind=8),dimension(3)                  :: wr_sigma,wi_sigma
real(kind=8)                               :: S1_sigma,S2_sigma,S3_sigma,S_P,cpt_p
real(kind=8),dimension(3,3)                :: localframe_sigma
real(kind=8),dimension(3,3)                :: M_sigma
!! POUR LE CALCUL DU PROFIL DE I !!
real(kind=8),dimension(:,:),allocatable    :: Moyenne_I_moyen
real(kind=8)                               :: Mean_pressure

if (compteur_clout>=debut_profils) then

  print*,'.... Profils ....'

  dmoyen = 0  
  cpt_p  = 0
  do i=1,nbParticules
    if (TAB_POLYEDRE(i)%behav /= 'plexx')  cycle
    dmoyen = TAB_POLYEDRE(i)%Rmax + dmoyen
    cpt_p = cpt_p+1
  end do
  dmoyen = 2*dmoyen / real(cpt_p,8)
  ep = dmoyen
  print*,'d_moyen = ',dmoyen
!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  if (the_first_time_profils) then
    nb_ligne_tab = 0
    do
      nb_ligne_tab  =   nb_ligne_tab + 1
      ep = ep + dmoyen/10
      if (ep>hauteur+dmoyen) exit
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
      do j=1,nb_ligne_tab+100
        TAB_contrainte(i)%tab_sigma(j)%sigma = 0
      end do
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
!???

!! POUR LE CALCUL DU PROFIL DE CONTRAINTE !!
!! non, la il n y a rien a initialiser

!! POUR LE CALCUL DU PROFIL DE I !!
!! non, la il n y a rien a initialiser

    cpt_v    = 0
    cpt_c    = 0
    ep = ep+dmoyen/10
    nb_ligne_tab = nb_ligne_tab + 1

!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!! Boucle sur les particules
    do i=1,nbParticules
      if (TAB_POLYEDRE(i)%behav/='plexx') cycle
      
      if ((TAB_POLYEDRE(i)%centre(3)<ep+dmoyen/10).and.(TAB_POLYEDRE(i)%centre(3)>ep)) then
!! POUR LE CALCUL DU PROFIL DE VITESSE !!
        Vmoyen = Vmoyen + sqrt(TAB_POLYEDRE(i)%Vx**2 + TAB_POLYEDRE(i)%Vy**2 + TAB_POLYEDRE(i)%Vz**2)
        Vxmoyen = Vxmoyen + TAB_POLYEDRE(i)%Vx 
        Vymoyen = Vymoyen + TAB_POLYEDRE(i)%Vy  
        Vzmoyen = Vzmoyen + TAB_POLYEDRE(i)%Vz  

        Vrmoyen = Vrmoyen + sqrt(TAB_POLYEDRE(i)%Wx**2 + TAB_POLYEDRE(i)%Wy**2 + TAB_POLYEDRE(i)%Wz**2)
        Vrxmoyen = Vrxmoyen + (TAB_POLYEDRE(i)%Wx)
        Vrymoyen = Vrymoyen + (TAB_POLYEDRE(i)%Wy) 
        Vrzmoyen = Vrzmoyen + (TAB_POLYEDRE(i)%Wz) 
    
        cpt_v = cpt_v + 1
      end if
    end do
!! Boucle sur les contacts
    do i=1,nb_ligneCONTACT_POLYR
      if ( (TAB_CONTACT_POLYR(i)%coor_ctc(3) < ep+ dmoyen/10).and.(TAB_CONTACT_POLYR(i)%coor_ctc(3) > ep) ) then
        if  (TAB_CONTACT_POLYR(i)%nature == 'SPPLx') cycle
        cd   = TAB_CONTACT_POLYR(i)%icdent
        an   = TAB_CONTACT_POLYR(i)%ianent
        nik  = TAB_CONTACT_POLYR(i)%n
        tik  = TAB_CONTACT_POLYR(i)%t
        sik  = TAB_CONTACT_POLYR(i)%s
        Rtik = TAB_CONTACT_POLYR(i)%rt
        Rnik = TAB_CONTACT_POLYR(i)%rn
        Rsik = TAB_CONTACT_POLYR(i)%rs
        if (Rnik==0) cycle
        Lik(1:3) = TAB_POLYEDRE(cd)%centre(1:3)-TAB_POLYEDRE(an)%centre(1:3)
        Fik(1:3) = (Rnik*nik(1:3)+Rtik*tik(1:3)+Rsik*sik(1:3))
    
    
        if ( sqrt( Lik(1)**2 + Lik(2)**2 + Lik(3)**2) > 2*( TAB_POLYEDRE(cd)%Rmax + TAB_POLYEDRE(an)%Rmax ) ) then
          Rayon_max = sqrt( ( Lik(1)**2 + Lik(2)**2 + Lik(3)**2) ) 
          lik(:) = Lik(:) / Rayon_max
          Lik(:) = abs(TAB_POLYEDRE(cd)%Rmax+TAB_POLYEDRE(an)%Rmax)*Lik(:)
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
    
    vboite = longueur*largeur*dmoyen/10

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

!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    

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

    if (ep>hauteur+dmoyen) exit
  end do

!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  if (compteur_clout == fin) then
!! POUR LE CALCUL DU PROFIL DE VITESSE !!
    allocate (Moyenne_vitesse_moyenne(nb_ligne_tab+100,4))
    allocate(Moyenne_vitesse_moyenne_r(nb_ligne_tab+100,5))
    Moyenne_vitesse_moyenne(:,:) = 0.D0
    Moyenne_vitesse_moyenne_r(:,:) = 0.D0
!! POUR LE CALCUL DU PROFIL DE CONTRAINTE !!
    allocate (Moyenne_contrainte_moyenne(nb_ligne_tab+100))
    do i=1,nb_ligne_tab+100
      Moyenne_contrainte_moyenne(i)%sigma = 0.D0
    end do
!! POUR LE CALCUL DU PROFIL DE I !!
    allocate (Moyenne_I_moyen(nb_ligne_tab+100,4))

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
      Moyenne_I_moyen(i,2) = ( Moyenne_vitesse_moyenne(i+1,3)-Moyenne_vitesse_moyenne(i,3) ) / &      ! ATTENTION IL FAUT SE RAPPELER QU ON 
                             ( dmoyen*(Moyenne_vitesse_moyenne(i+1,1)-Moyenne_vitesse_moyenne(i,1)) )  
      Moyenne_I_moyen(i,3) =  Moyenne_I_moyen(i,2) * dmoyen * sqrt(2800/Mean_pressure)
      Moyenne_I_moyen(i,4) =  abs(2*Moyenne_vitesse_moyenne_r(i,2)) * dmoyen * sqrt(2800/Mean_pressure)
    end do    

!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!! POUR LE CALCUL DU PROFIL DE VITESSE !!
    do i=1,nb_ligne_tab+20
      write(131,'(11(1X,D14.7))') Moyenne_vitesse_moyenne(i,1),Moyenne_vitesse_moyenne(i,2),&
                                  Moyenne_vitesse_moyenne(i,3),Moyenne_vitesse_moyenne(i,4),&
                                  Moyenne_vitesse_moyenne_r(i,2),Moyenne_vitesse_moyenne_r(i,3),&
                                  Moyenne_vitesse_moyenne_r(i,4)
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

      write(132,'(11(1X,D14.7))') Moyenne_vitesse_moyenne(i,1), (S3_sigma-S1_sigma)/(S3_sigma+S1_sigma),S_P
      
    end do

!! POUR LE CALCUL DU PROFIL DE I !!
    do i=1,nb_ligne_tab+20
      write(133,'(11(1X,D14.7))') Moyenne_vitesse_moyenne(i,1),Moyenne_I_moyen(i,2),Moyenne_I_moyen(i,3),Moyenne_I_moyen(i,4)
    end do

!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  end if
end if  
end subroutine profils



!==============================================================================
! Orientation des contacts
!==============================================================================
subroutine contact_orientation
implicit none
integer                                             ::  i,cd,an,Nsect,j,cpt_theta,cpt_phi,cpt_zero,k,cpt
real(kind=8)                                        ::  sect,sect2,theta,phi,fmoyen,norm,SinTheta,ThetaMax,ThetaMin,ThetaMean,&
                                                        Rtik,Rnik,Rsik,FN,FT,LN,LT
real(kind=8),dimension(3)                           ::  nik,Lik,tik,Fik,sik
real(kind=8),dimension(:),allocatable               ::  tab_alpha,tab_alpha2
real(kind=8),dimension(:),allocatable               ::  tab_number_n_theta,tab_number_n_xy,tab_number_nb_in_class
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
do i=1,nb_ligneCONTACT_POLYR
  if  (TAB_CONTACT_POLYR(i)%nature == 'PRPLx') cycle 
  if  (TAB_CONTACT_POLYR(i)%rn==0.D0) cycle

  nik  = TAB_CONTACT_POLYR(i)%n
  
  Fabric(1,1:3) = nik(1)*nik(1:3) + Fabric(1,1:3)
  Fabric(2,1:3) = nik(2)*nik(1:3) + Fabric(2,1:3)
  Fabric(3,1:3) = nik(3)*nik(1:3) + Fabric(3,1:3)
  cpt=cpt+1
end do
Fabric(:,:) = Fabric(:,:) / real(cpt,8)
print*,Fabric(1,1),Fabric(2,2),Fabric(3,3),5*(Fabric(3,3) - Fabric(1,1))/2
cpt = 0

!----------------ORIENTATION DES NORMALES DANS LE PLAN XY ------------------
cpt_phi   = 0
do i=1,nb_ligneCONTACT_POLYR
  if (TAB_CONTACT_POLYR(i)%nature == 'PRPLx') cycle
  if (TAB_CONTACT_POLYR(i)%rn==0.D0)          cycle
  cd  = TAB_CONTACT_POLYR(i)%icdent
  an  = TAB_CONTACT_POLYR(i)%ianent
  Lik(1:3) = TAB_POLYEDRE(cd)%centre(1:3)-TAB_POLYEDRE(an)%centre(1:3)
  norm = sqrt( Lik(1)**2 + Lik(2)**2 + Lik(3)**2 )
!  Lik(1:3) = Lik(1:3) / norm
  Lik = TAB_CONTACT_POLYR(i)%n
!  if (Lik(1)**2+Lik(2)**2 ==0.D0 ) cycle
  cpt_phi = cpt_phi + 1
end do
print*, cpt_phi

do i=1,Nsect
  do j=1,nb_ligneCONTACT_POLYR
    if (TAB_CONTACT_POLYR(j)%nature == 'PRPLx') cycle
    if (TAB_CONTACT_POLYR(j)%rn==0.D0)          cycle

    cd  = TAB_CONTACT_POLYR(j)%icdent
    an  = TAB_CONTACT_POLYR(j)%ianent
    Lik(1:3) = TAB_POLYEDRE(cd)%centre(1:3)-TAB_POLYEDRE(an)%centre(1:3)
    norm = sqrt( Lik(1)**2 + Lik(2)**2 + Lik(3)**2 )

    nik  = TAB_CONTACT_POLYR(j)%n
!    nik  = Lik(1:3) / norm

    if (nik(1)**2+nik(2)**2 ==0.D0 ) cycle

    phi = acos( nik(1) / sqrt(nik(1)**2+nik(2)**2) )
        
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
allocate(tab_number_nb_in_class(Nsect))

tab_number_n_theta  = 0
tab_number_nb_in_class = 0
cpt_theta   = 0
do i=1,nb_ligneCONTACT_POLYR
  if (TAB_CONTACT_POLYR(i)%nature == 'PRPLx') cycle
  if (TAB_CONTACT_POLYR(i)%rn==0.D0)          cycle
  cd  = TAB_CONTACT_POLYR(i)%icdent
  an  = TAB_CONTACT_POLYR(i)%ianent
  Lik(1:3) = TAB_POLYEDRE(cd)%centre(1:3)-TAB_POLYEDRE(an)%centre(1:3)
  norm = sqrt( Lik(1)**2 + Lik(2)**2 + Lik(3)**2 )

  nik  = TAB_CONTACT_POLYR(i)%n
!  nik  = Lik(1:3) / norm
  
  
  if (sqrt(nik(1)**2+nik(2)**2) < 0.001) then
    cpt_theta = cpt_theta + 1
  else

    phi = acos( (nik(1)) / sqrt(nik(1)**2+nik(2)**2) )
    if ( (phi> Pi/2 + 2*sect2).and.(phi< 3*Pi/2 - 2*sect2)  ) cycle
    if ( (phi< Pi/2 - 2*sect2).and.(phi> 3*Pi/2 + 2*sect2)  ) cycle

    cpt_theta = cpt_theta + 1

  end if
  
end do

do i=1,Nsect
  do j=1,nb_ligneCONTACT_POLYR
    if (TAB_CONTACT_POLYR(j)%nature == 'PRPLx') cycle
    if (TAB_CONTACT_POLYR(j)%rn==0.D0)          cycle

    cd  = TAB_CONTACT_POLYR(j)%icdent
    an  = TAB_CONTACT_POLYR(j)%ianent
    Lik(1:3) = TAB_POLYEDRE(cd)%centre(1:3)-TAB_POLYEDRE(an)%centre(1:3)
    norm = sqrt( Lik(1)**2 + Lik(2)**2 + Lik(3)**2 )

    Rtik = TAB_CONTACT_POLYR(j)%rt
    Rnik = TAB_CONTACT_POLYR(j)%rn
    Rsik = TAB_CONTACT_POLYR(j)%rs
    nik  = TAB_CONTACT_POLYR(j)%n
    tik  = TAB_CONTACT_POLYR(j)%t
    sik  = TAB_CONTACT_POLYR(j)%s
    Fik(1:3) = (Rnik*nik(1:3)+Rtik*tik(1:3)+Rsik*sik(1:3))


    nik  = TAB_CONTACT_POLYR(j)%n
!    nik  = Lik(1:3) / norm

    FN       =  (Fik(1)*nik(1) + Fik(2)*nik(2) + Fik(3)*nik(3))
    tik(1:3) =  ( Fik(1:3) - FN*nik(1:3) )
    
    norm       = sqrt( tik(1)**2 + tik(2)**2 + tik(3)**2 )
!    print*, norm,tik
    if (norm>0.000001) then
      tik(1:3) = tik(1:3) / norm
      FT = Rtik + Rsik !(Fik(1)*tik(1) + Fik(2)*tik(2) + Fik(3)*tik(3))
    else
       FT =0.D0
    end if
!    if ( nik(3) < 0.D0 )  nik(:)=-nik(:)

    if (sqrt(nik(1)**2+nik(2)**2) < 0.001) then
      cycle
      theta = Pi/2-0.001
    else
      phi = acos( (nik(1)) / sqrt(nik(1)**2+nik(2)**2) )

      if ( (phi> Pi/2 + 2*sect2).and.(phi< 3*Pi/2 - 2*sect2)  ) cycle
      if ( (phi< Pi/2 - 2*sect2).and.(phi> 3*Pi/2 + 2*sect2)  ) cycle

      theta = acos( nik(2) / sqrt( nik(2)**2+nik(3)**2 ) )

    end if
    

    if (i==1) then
      if ( (theta)<=(tab_alpha(i)) ) then
        tab_number_n_theta(i) = tab_number_n_theta(i) + FT !1 
        tab_number_nb_in_class(i) = tab_number_nb_in_class(i) + 1.0
      end if
    end if

    if ( (i>1).and.(i<Nsect)) then
      if ( ((theta)>=(tab_alpha(i-1)) ) .and. ( (theta)<=(tab_alpha(i)) ) ) then
        tab_number_n_theta(i) = tab_number_n_theta(i) + FT !1 
        tab_number_nb_in_class(i) = tab_number_nb_in_class(i) + 1.0
      end if
    end if

    if (i==Nsect) then
      if ( ((theta)<=(tab_alpha(i))) .and. ((theta)>=(tab_alpha(i-1))) ) then
        tab_number_n_theta(i) = tab_number_n_theta(i) + FT !1
        tab_number_nb_in_class(i) = tab_number_nb_in_class(i) + 1.0
      end if
    end if
  end do
end do
! tab_number_n_theta(:)  = tab_number_n_theta(:) / real(cpt_theta,8)
 tab_number_n_theta(:)  = tab_number_n_theta(:) / tab_number_nb_in_class(:)


do i=1,Nsect
  write(124,'(8(1X,D14.7))') sect*(i-1) , tab_number_n_theta(i)
  write(124,'(8(1X,D14.7))') sect*(i)   , tab_number_n_theta(i)
end do
do i=1,Nsect
  write(124,'(8(1X,D14.7))') sect*(i-1)+PI , tab_number_n_theta(i)
  write(124,'(8(1X,D14.7))') sect*(i)+PI   , tab_number_n_theta(i)
end do

stop
end subroutine contact_orientation




!==============================================================================
! Calcul du Ksi-reseau 
!==============================================================================
subroutine ksi_reseau
implicit none
real(kind=8),dimension(3,3)              :: Fabric,Fabric_LN,Fabric_LT,Fabric_L,Moment
real(kind=8),dimension(3,3)              :: Fabric_FN,Fabric_FT,Fabric_F
real(kind=8),dimension(3)                :: nik,tik,sik,Lik,Fik
real(kind=8)                             :: Rtik,Rnik,Rsik,norm_FN,FN,FT,LN,LT,Norm_F,Norm_L,cpt,epsilon,cpts,cptd,cptt,&
                                            norm_FT,LN_max,Sx,Sy,Sxy,Pearson
integer                                  :: i,j,cd,an,nb_ksi
real(kind=8)                             :: a,aln,alt,afn,aft,qp,ksi,PRESSION_FT,PRESSION_FN,PRESSION_P,PRESSION_L,&
                                            ksimple,ksimple_ss,kdouble,ktriple,kmobilized,fmobilized,nb_ksi_ctc,nb_ksi_ln_ctc,&
                                            k_ss_weak,k_ss_strong,&
                                            k_fv_weak,k_fv_strong,&
                                            k_ff_weak,k_ff_strong,&
                                            k_fs_weak,k_fs_strong,&
                                            k_weak,k_strong
real(kind=8),dimension(:,:),allocatable  :: mean_val


if ( (compteur_clout >= step_ksi_reseau).and.(compteur_clout <= stop_ksi_reseau) ) then
  print*, '---> Ksi-reseau'
  Lik       = 0.D0
  Fik       = 0.D0
  Norm_L    = 0.D0
  Norm_F    = 0.D0
  norm_FN   = 0.D0
  norm_FT   = 0.D0
  cpt       = 0.D0
  cpts      = 0.D0
  cptd      = 0.D0
  cptt      = 0.D0
  Fabric    = 0.D0
  Fabric_LN = 0.D0
  Fabric_LT = 0.D0
  Fabric_FN = 0.D0
  Fabric_FT = 0.D0
  Moment    = 0.D0
  LN_max    = 0.D0
  Sx        = 0.D0
  Sy        = 0.D0
  Sxy       = 0.D0
  Pearson   = 0.D0
  nb_ksi_ln_ctc = 0.D0
  
  k_ss_weak = 0.D0
  k_ss_strong = 0.D0
  k_fv_weak  = 0.D0 
  k_fv_strong  = 0.D0
  k_ff_weak = 0.D0
  k_ff_strong  = 0.D0
  k_fs_weak  = 0.D0
  k_fs_strong  = 0.D0
  k_weak  = 0.D0
  k_strong = 0.D0

  do i=1,nb_ligneCONTACT_POLYR
    if  (TAB_CONTACT_POLYR(i)%nature == 'PRPLx') cycle 
    if  (TAB_CONTACT_POLYR(i)%rn==0.D0) cycle
    cd   = TAB_CONTACT_POLYR(i)%icdent
    an   = TAB_CONTACT_POLYR(i)%ianent
    nik  = TAB_CONTACT_POLYR(i)%n
    tik  = TAB_CONTACT_POLYR(i)%t
    sik  = TAB_CONTACT_POLYR(i)%s
    Rtik = TAB_CONTACT_POLYR(i)%rt
    Rnik = TAB_CONTACT_POLYR(i)%rn
    Rsik = TAB_CONTACT_POLYR(i)%rs

    Lik(1:3) = TAB_POLYEDRE(cd)%centre(1:3)-TAB_POLYEDRE(an)%centre(1:3)
    Fik(1:3) = (Rnik*nik(1:3)+Rtik*tik(1:3)+Rsik*sik(1:3))
    
    Norm_F   = Norm_F + sqrt( Fik(1)**2 + Fik(2)**2 + Fik(3)**2 )
    Norm_L   = Norm_L + sqrt( Lik(1)**2 + Lik(2)**2 + Lik(3)**2 )
    LN       =  sqrt( Lik(1)**2 + Lik(2)**2 + Lik(3)**2 ) 
    LN_max   = max( LN_max, 0.5*LN )
    
    nik(1:3) =  Lik(1:3)  / LN

    FN       =  abs(Fik(1)*nik(1) + Fik(2)*nik(2) + Fik(3)*nik(3))
    norm_FN  = norm_FN + FN

    tik(1:3) =  Fik(1:3) - FN*nik(1:3)
    FT       = sqrt( tik(1)**2 + tik(2)**2 + tik(3)**2 )
    norm_FT  = norm_FT + FT
    
    if ( FT > 0.00000000001) tik(1:3) = tik(1:3) / FT

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

    cpt = cpt + 1
    if (TAB_CONTACT_POLYR(i)%type == 1) cpts = cpts + 1.D0
    if (TAB_CONTACT_POLYR(i)%type == 2) cptd = cptd + 1.D0
    if (TAB_CONTACT_POLYR(i)%type >  2) cptt = cptt + 1.D0

  end do
  norm_FN     = norm_FN / cpt
  norm_FT     = norm_FT / cpt
  norm_F      = norm_F  / cpt
  norm_L      = norm_L  / cpt
  
  Fabric      = Fabric(:,:) / cpt
  Fabric_LN   = Fabric_LN   / (Norm_L * cpt)
  Fabric_FN   = Fabric_FN    / (Norm_F * cpt)
  Fabric_FT   = Fabric_FT    / (Norm_F * cpt)
  Fabric_F    = Fabric_FT + Fabric_FN

  PRESSION_P  = (Moment(1,1)+Moment(2,2)+Moment(3,3)) / 3
  PRESSION_L  = (Fabric_LN(1,1)+Fabric_LN(2,2)+Fabric_LN(3,3))
  PRESSION_FN = (Fabric_FN(1,1)+Fabric_FN(2,2)+Fabric_FN(3,3))
  PRESSION_FT = (Fabric_F(1,1)+Fabric_F(2,2)+Fabric_F(3,3))
!-------------------------------------------------------------------!
! Calcul du coefficent de correlation
  do i=1,nb_ligneCONTACT_POLYR
    if  (TAB_CONTACT_POLYR(i)%nature == 'PRPLx') cycle 
    if  (TAB_CONTACT_POLYR(i)%rn==0.D0) cycle
    cd   = TAB_CONTACT_POLYR(i)%icdent
    an   = TAB_CONTACT_POLYR(i)%ianent
    nik  = TAB_CONTACT_POLYR(i)%n
    tik  = TAB_CONTACT_POLYR(i)%t
    sik  = TAB_CONTACT_POLYR(i)%s
    Rtik = TAB_CONTACT_POLYR(i)%rt
    Rnik = TAB_CONTACT_POLYR(i)%rn
    Rsik = TAB_CONTACT_POLYR(i)%rs

    Lik(1:3) = TAB_POLYEDRE(cd)%centre(1:3)-TAB_POLYEDRE(an)%centre(1:3)
    Fik(1:3) = (Rnik*nik(1:3)+Rtik*tik(1:3)+Rsik*sik(1:3))

    Sx  = Sx  + ( sqrt( Lik(1)**2 + Lik(2)**2 + Lik(3)**2 ) - Norm_L )**2
    Sy  = Sy  + ( sqrt( Fik(1)**2 + Fik(2)**2 + Fik(3)**2 ) - Norm_F )**2
    Sxy = Sxy + ( sqrt( Lik(1)**2 + Lik(2)**2 + Lik(3)**2 ) - Norm_L )*&
                ( sqrt( Fik(1)**2 + Fik(2)**2 + Fik(3)**2 ) - Norm_F )
  end do
  Pearson = Sxy/( sqrt(Sx) * sqrt(Sy) )
  print*, 'Pearson coefficient = ',Pearson
!-------------------------------------------------------------------! Calcul de la proportion de chaque contacts dans les reseaux
  do i=1,nb_ligneCONTACT_POLYR
    if  (TAB_CONTACT_POLYR(i)%nature == 'PRPLx') cycle 
    if  (TAB_CONTACT_POLYR(i)%rn==0.D0) cycle
    cd   = TAB_CONTACT_POLYR(i)%icdent
    an   = TAB_CONTACT_POLYR(i)%ianent
    nik  = TAB_CONTACT_POLYR(i)%n
    tik  = TAB_CONTACT_POLYR(i)%t
    sik  = TAB_CONTACT_POLYR(i)%s
    Rtik = TAB_CONTACT_POLYR(i)%rt
    Rnik = TAB_CONTACT_POLYR(i)%rn
    Rsik = TAB_CONTACT_POLYR(i)%rs

    Lik(1:3) = TAB_POLYEDRE(cd)%centre(1:3)-TAB_POLYEDRE(an)%centre(1:3)
    Fik(1:3) = (Rnik*nik(1:3)+Rtik*tik(1:3)+Rsik*sik(1:3))
    
    LN       =  sqrt( Lik(1)**2 + Lik(2)**2 + Lik(3)**2 ) 

    nik(1:3) =  Lik(1:3)  / LN
    FN       =  abs(Fik(1)*nik(1) + Fik(2)*nik(2) + Fik(3)*nik(3))
    
    if  ( FN/norm_FN > 1.D0 ) then
      k_strong = k_strong + 1
      if (TAB_CONTACT_POLYR(i)%type == 0) k_ss_strong = k_ss_strong + 1.D0
      if (TAB_CONTACT_POLYR(i)%type == 1) k_fv_strong = k_fv_strong + 1.D0
      if (TAB_CONTACT_POLYR(i)%type == 2) k_fs_strong = k_fs_strong + 1.D0
      if (TAB_CONTACT_POLYR(i)%type >  2) k_ff_strong = k_ff_strong + 1.D0
    else
      k_weak   = k_weak + 1
      if (TAB_CONTACT_POLYR(i)%type == 0) k_ss_weak = k_ss_weak + 1.D0
      if (TAB_CONTACT_POLYR(i)%type == 1) k_fv_weak = k_fv_weak + 1.D0
      if (TAB_CONTACT_POLYR(i)%type == 2) k_fs_weak = k_fs_weak + 1.D0
      if (TAB_CONTACT_POLYR(i)%type >  2) k_ff_weak = k_ff_weak + 1.D0
    end if

  end do

  write(108,'(40(1X,D12.5))') temps,Def,Hauteur,k_strong/cpt,k_weak/cpt,&
                                              k_ss_strong/k_strong, k_ss_weak/k_weak,&
                                              k_fv_strong/k_strong, k_fv_weak/k_weak,&
                                              k_fs_strong/k_strong, k_fs_weak/k_weak,&
                                              k_ff_strong/k_strong, k_ff_weak/k_weak
!-------------------------------------------------------------------!
  ksi = 0.D0
  nb_ksi = 0
  epsilon = 0.05
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
    ksimple   = 0.D0
    ksimple_ss= 0.D0
    kdouble   = 0.D0
    ktriple   = 0.D0
    kmobilized= 0.D0
    fmobilized= 0.D0
    nb_ksi_ctc= 0.D0
    nb_ksi_ln_ctc = 0.D0
    do i=1,nb_ligneCONTACT_POLYR
      if  (TAB_CONTACT_POLYR(i)%nature == 'PRPLx') cycle 
      if  (TAB_CONTACT_POLYR(i)%rn==0.D0) cycle
      
      cd   = TAB_CONTACT_POLYR(i)%icdent
      an   = TAB_CONTACT_POLYR(i)%ianent
      nik  = TAB_CONTACT_POLYR(i)%n
      tik  = TAB_CONTACT_POLYR(i)%t
      sik  = TAB_CONTACT_POLYR(i)%s
      Rtik = TAB_CONTACT_POLYR(i)%rt
      Rnik = TAB_CONTACT_POLYR(i)%rn
      Rsik = TAB_CONTACT_POLYR(i)%rs

      Lik(1:3) = TAB_POLYEDRE(cd)%centre(1:3)-TAB_POLYEDRE(an)%centre(1:3)
      Fik(1:3) = (Rnik*nik(1:3)+Rtik*tik(1:3)+Rsik*sik(1:3))
    
      LN       =  sqrt( Lik(1)**2 + Lik(2)**2 + Lik(3)**2 ) 

      nik(1:3) =  Lik(1:3)  / LN
      FN       =  abs(Fik(1)*nik(1) + Fik(2)*nik(2) + Fik(3)*nik(3))

      tik(1:3) =  Fik(1:3) - FN*nik(1:3)
      FT       = sqrt( tik(1)**2 + tik(2)**2 + tik(3)**2 )
      if ( FT > 0.00000000001) tik(1:3) = tik(1:3) / FT

!      if  ( ((FN/norm_FN) > ksi).and.(FN/norm_Fn ) < ksi+ epsilon) then
      if  ( (FN/norm_Fn ) < ksi+ epsilon) then

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
        
        if (TAB_CONTACT_POLYR(i)%type == 0) ksimple_ss = ksimple_ss + 1.D0
        if (TAB_CONTACT_POLYR(i)%type == 1) ksimple = ksimple + 1.D0
        if (TAB_CONTACT_POLYR(i)%type == 2) kdouble = kdouble + 1.D0
        if (TAB_CONTACT_POLYR(i)%type >  2) ktriple = ktriple + 1.D0
        
        nb_ksi_ctc = nb_ksi_ctc + 1.D0
        
        if ( abs( sqrt(TAB_CONTACT_POLYR(i)%rt**2+TAB_CONTACT_POLYR(i)%rs**2) ) > 0.39 * TAB_CONTACT_POLYR(i)%rn ) &
                                   kmobilized = kmobilized+1.D0
                                   
        fmobilized = fmobilized + &
                      abs( sqrt(TAB_CONTACT_POLYR(i)%rt**2+TAB_CONTACT_POLYR(i)%rs**2) ) / ( 0.4 * TAB_CONTACT_POLYR(i)%rn )
      end if

      if  ( ((LN/LN_max) > ksi).and.(LN/LN_max ) < ksi+ epsilon) then

        nb_ksi_ln_ctc = nb_ksi_ln_ctc + 1.D0

      end if

    end do
    Fabric      = Fabric(:,:) / cpt
    Fabric_LN   = Fabric_LN   / (Norm_L * cpt)
    Fabric_FN   = Fabric_FN   / (Norm_F * cpt)
    Fabric_FT   = Fabric_FT   / (Norm_F * cpt)
    Fabric_F    = Fabric_FT + Fabric_FN

    nb_ksi_ctc = nb_ksi_ctc / cpt
    ksimple    = ksimple / cpt
    ksimple_ss = ksimple_ss / cpt
    kdouble    = kdouble / cpt
    ktriple    = ktriple / cpt
    kmobilized = kmobilized / cpt
    fmobilized = fmobilized / cpt

    nb_ksi_ln_ctc = nb_ksi_ln_ctc / cpt
    
    qp = ((Moment(3,3)-Moment(1,1))/3) / PRESSION_P
    a  = 5*(Fabric(3,3) - Fabric(1,1)) / 2
    aln= 5*(Fabric_LN(3,3) - Fabric_LN(1,1)) / (2*PRESSION_L) - a
    afn= 5*(Fabric_FN(3,3) - Fabric_FN(1,1)) / (2*PRESSION_FN) - a
    aft= 5*(Fabric_F(3,3) - Fabric_F(1,1)) / (2*PRESSION_FT) - a - afn

!    write(107,'(16(1X,D12.5))') ksi-epsilon, qp, a,aln,afn,aft,2*(a+aln+afn+aft)/5

    Mean_Ksi_reseau(pas_ksi,nb_ksi,1)  = ksi
    Mean_Ksi_reseau(pas_ksi,nb_ksi,2)  = qp
    Mean_Ksi_reseau(pas_ksi,nb_ksi,3)  = a
    Mean_Ksi_reseau(pas_ksi,nb_ksi,4)  = aln
    Mean_Ksi_reseau(pas_ksi,nb_ksi,5)  = afn
    Mean_Ksi_reseau(pas_ksi,nb_ksi,6)  = aft
    Mean_Ksi_reseau(pas_ksi,nb_ksi,7)  = 2*(a+aln+afn+aft)/5
    Mean_Ksi_reseau(pas_ksi,nb_ksi,8)  = ksimple
    Mean_Ksi_reseau(pas_ksi,nb_ksi,9)  = ksimple_ss
    Mean_Ksi_reseau(pas_ksi,nb_ksi,10)  = kdouble
    Mean_Ksi_reseau(pas_ksi,nb_ksi,11) = ktriple
    Mean_Ksi_reseau(pas_ksi,nb_ksi,12) = kmobilized
    Mean_Ksi_reseau(pas_ksi,nb_ksi,13) = fmobilized
    Mean_Ksi_reseau(pas_ksi,nb_ksi,14) = nb_ksi_ctc
    Mean_Ksi_reseau(pas_ksi,nb_ksi,15) = nb_ksi_ln_ctc
    
    if (ksi>15) exit

  end do
  all_pas = all_pas + 1
  if ( (compteur_clout == stop_ksi_reseau ) ) then
    allocate(mean_val(15,nb_ksi))
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
        mean_val(11,j) = mean_val(11,j) + Mean_Ksi_reseau(i,j,11)
        mean_val(12,j) = mean_val(12,j) + Mean_Ksi_reseau(i,j,12)
        mean_val(13,j) = mean_val(13,j) + Mean_Ksi_reseau(i,j,13)
        mean_val(14,j) = mean_val(14,j) + Mean_Ksi_reseau(i,j,14)
        mean_val(15,j) = mean_val(15,j) + Mean_Ksi_reseau(i,j,15)
      end do
    end do
    mean_val(:,:) = mean_val(:,:) / pas_ksi


    i = 0
    do 
      i = i + 1
      if (mean_val(1,i)>15) exit
      if (mean_val(8,i)+mean_val(9,i)+mean_val(10,i)==0) cycle
      write(107,'(20(1X,D12.5))') mean_val(1,i), mean_val(2,i), mean_val(3,i),&
                                                         mean_val(4,i), mean_val(5,i),&
                                                         mean_val(6,i), mean_val(7,i),&
                                                         mean_val(8,i), mean_val(9,i),&
                                                         mean_val(10,i), mean_val(11,i),&
                                                         mean_val(12,i),mean_val(13,i),&
                                                         mean_val(14,i),mean_val(15,i)
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
integer                                  :: i,cd,an
real(kind=8)                             :: norm,increment


if ( (compteur_clout >= step_pdf_start).and.(compteur_clout <= step_pdf_stop) ) then
  print*, '---> PDF'
  do i = 1, nb_ligneCONTACT_POLYR
    if  (TAB_CONTACT_POLYR(i)%rn     == 0.D0   ) cycle
    write(106,'(4(1X,D12.5))') TAB_CONTACT_POLYR(i)%rn,&
                                sqrt(TAB_CONTACT_POLYR(i)%rs**2+TAB_CONTACT_POLYR(i)%rt**2)
  end do
  all_pas = all_pas + 1
else if ( (compteur_clout >= step_pdf_stop ) ) then
  stop
end if 

end subroutine pdf

!==============================================================================
! Branche des differentes anisotropies
!==============================================================================
subroutine branche_anisotropies
implicit none
real(kind=8),dimension(3,3)              :: M

real(kind=8),dimension(3,3)              :: Fabric,Fabric_LN,Fabric_L
real(kind=8),dimension(3,3)              :: Fabric_FN,Fabric_FT,Fabric_F

real(kind=8),dimension(3,3)              :: Fabric_s,Fabric_LN_s,Fabric_LT_s,Fabric_L_s
real(kind=8),dimension(3,3)              :: Fabric_FN_s,Fabric_FT_s,Fabric_F_s

real(kind=8),dimension(3,3)              :: Fabric_ss,Fabric_LN_ss,Fabric_LT_ss,Fabric_L_ss
real(kind=8),dimension(3,3)              :: Fabric_FN_ss,Fabric_FT_ss,Fabric_F_ss

real(kind=8),dimension(3,3)              :: Fabric_d,Fabric_LN_d,Fabric_LT_d,Fabric_L_d
real(kind=8),dimension(3,3)              :: Fabric_FN_d,Fabric_FT_d,Fabric_F_d

real(kind=8),dimension(3,3)              :: Fabric_t,Fabric_LN_t,Fabric_LT_t,Fabric_L_t
real(kind=8),dimension(3,3)              :: Fabric_FN_t,Fabric_FT_t,Fabric_F_t

real(kind=8),dimension(3)                :: nik,tik,sik,Lik,Fik
real(kind=8)                             :: Rtik,Rnik,Rsik,norm,FN,FT,LN,Norm_F,Norm_L,cpt,Rayon_max
integer                                  :: i,cd,an
real(kind=8)                             :: a,aln,alt,afn,aft
real(kind=8)                             :: as,alns,afns,afts
real(kind=8)                             :: ass,alnss,afnss,aftss
real(kind=8)                             :: ad,alnd,afnd,aftd
real(kind=8)                             :: at,alnt,afnt,aftt

real(kind=8)                             :: S1,S2,S3
real(kind=8),dimension(3,3)              :: localframe  
integer                                  :: ierror,matz,lda
real(kind=8),dimension(3)                :: wr,wi

print*,'.... branche anisotropies ....'

Lik       = 0.D0
Fik       = 0.D0
Norm_L    = 0.D0
Norm_F    = 0.D0
cpt       = 0.D0
Fabric    = 0.D0
Fabric_LN = 0.D0
Fabric_FN = 0.D0
Fabric_FT = 0.D0

Fabric_s   = 0.D0
Fabric_LN_s = 0.D0
Fabric_FN_s = 0.D0
Fabric_FT_s = 0.D0

Fabric_ss   = 0.D0
Fabric_LN_ss = 0.D0
Fabric_FN_ss = 0.D0
Fabric_FT_ss = 0.D0

Fabric_d   = 0.D0
Fabric_LN_d = 0.D0
Fabric_FN_d = 0.D0
Fabric_FT_d = 0.D0

Fabric_t   = 0.D0
Fabric_LN_t = 0.D0
Fabric_FN_t = 0.D0
Fabric_FT_t = 0.D0

do i=1,nb_ligneCONTACT_POLYR
  if (TAB_CONTACT_POLYR(i)%rn<0.000000000001*Mean_total_Normal_force) cycle
  cd   = TAB_CONTACT_POLYR(i)%icdent
  an   = TAB_CONTACT_POLYR(i)%ianent
  nik  = TAB_CONTACT_POLYR(i)%n
  tik  = TAB_CONTACT_POLYR(i)%t
  sik  = TAB_CONTACT_POLYR(i)%s
  Rtik = TAB_CONTACT_POLYR(i)%rt
  Rnik = TAB_CONTACT_POLYR(i)%rn
  Rsik = TAB_CONTACT_POLYR(i)%rs

  if ((TAB_POLYEDRE(cd)%behav/='plexx').or.(TAB_POLYEDRE(an)%behav/='plexx')) cycle

  Lik(1:3) = TAB_POLYEDRE(cd)%centre(1:3)-TAB_POLYEDRE(an)%centre(1:3)
  if ( sqrt( Lik(1)**2 + Lik(2)**2 + Lik(3)**2) > 2*( TAB_POLYEDRE(cd)%Rmax + TAB_POLYEDRE(an)%Rmax ) ) then
    Rayon_max = sqrt( ( Lik(1)**2 + Lik(2)**2 + Lik(3)**2) ) 
    Lik(:) = Lik(:) / Rayon_max
    Lik(:) = abs(TAB_POLYEDRE(cd)%Rmax+TAB_POLYEDRE(an)%Rmax)*Lik(:)
  end if

  Fik(1:3) =  (Rnik*nik(1:3)+Rtik*tik(1:3)+Rsik*sik(1:3))
  Norm_F   =  Norm_F + sqrt( Fik(1)**2 + Fik(2)**2 + Fik(3)**2 )
  
  LN       =  sqrt(Lik(1)**2+ Lik(2)**2 + Lik(3)**2)  
  Norm_L   =  LN
  
  nik(1:3) =  Lik(1:3)  / LN
  FN       =  abs(Fik(1)*nik(1) + Fik(2)*nik(2) + Fik(3)*nik(3))

  tik(1:3) =  Fik(1:3) - FN*nik(1:3)
  FT       = sqrt( tik(1)**2 + tik(2)**2 + tik(3)**2 )

  if ( FT > 0.00000000001) tik(1:3) = tik(1:3) / FT

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

  if (TAB_CONTACT_POLYR(i)%type == 0) then
    Fabric_ss(1,1:3) = nik(1)*nik(1:3) + Fabric_ss(1,1:3)
    Fabric_ss(2,1:3) = nik(2)*nik(1:3) + Fabric_ss(2,1:3)
    Fabric_ss(3,1:3) = nik(3)*nik(1:3) + Fabric_ss(3,1:3)

    Fabric_LN_ss(1,1:3) = LN*nik(1)*nik(1:3) + Fabric_LN_ss(1,1:3)
    Fabric_LN_ss(2,1:3) = LN*nik(2)*nik(1:3) + Fabric_LN_ss(2,1:3)
    Fabric_LN_ss(3,1:3) = LN*nik(3)*nik(1:3) + Fabric_LN_ss(3,1:3)

    Fabric_FN_ss(1,1:3) = FN*nik(1)*nik(1:3) + Fabric_FN_ss(1,1:3)
    Fabric_FN_ss(2,1:3) = FN*nik(2)*nik(1:3) + Fabric_FN_ss(2,1:3)
    Fabric_FN_ss(3,1:3) = FN*nik(3)*nik(1:3) + Fabric_FN_ss(3,1:3)
  
    Fabric_FT_ss(1,1:3) = FT*nik(1)*tik(1:3) + Fabric_FT_ss(1,1:3)
    Fabric_FT_ss(2,1:3) = FT*nik(2)*tik(1:3) + Fabric_FT_ss(2,1:3)
    Fabric_FT_ss(3,1:3) = FT*nik(3)*tik(1:3) + Fabric_FT_ss(3,1:3)

  else if (TAB_CONTACT_POLYR(i)%type == 1) then
    Fabric_s(1,1:3) = nik(1)*nik(1:3) + Fabric_s(1,1:3)
    Fabric_s(2,1:3) = nik(2)*nik(1:3) + Fabric_s(2,1:3)
    Fabric_s(3,1:3) = nik(3)*nik(1:3) + Fabric_s(3,1:3)

    Fabric_LN_s(1,1:3) = LN*nik(1)*nik(1:3) + Fabric_LN_s(1,1:3)
    Fabric_LN_s(2,1:3) = LN*nik(2)*nik(1:3) + Fabric_LN_s(2,1:3)
    Fabric_LN_s(3,1:3) = LN*nik(3)*nik(1:3) + Fabric_LN_s(3,1:3)

    Fabric_FN_s(1,1:3) = FN*nik(1)*nik(1:3) + Fabric_FN_s(1,1:3)
    Fabric_FN_s(2,1:3) = FN*nik(2)*nik(1:3) + Fabric_FN_s(2,1:3)
    Fabric_FN_s(3,1:3) = FN*nik(3)*nik(1:3) + Fabric_FN_s(3,1:3)
  
    Fabric_FT_s(1,1:3) = FT*nik(1)*tik(1:3) + Fabric_FT_s(1,1:3)
    Fabric_FT_s(2,1:3) = FT*nik(2)*tik(1:3) + Fabric_FT_s(2,1:3)
    Fabric_FT_s(3,1:3) = FT*nik(3)*tik(1:3) + Fabric_FT_s(3,1:3)
  else if (TAB_CONTACT_POLYR(i)%type == 2) then
    Fabric_d(1,1:3) = nik(1)*nik(1:3) + Fabric_d(1,1:3)
    Fabric_d(2,1:3) = nik(2)*nik(1:3) + Fabric_d(2,1:3)
    Fabric_d(3,1:3) = nik(3)*nik(1:3) + Fabric_d(3,1:3)

    Fabric_LN_d(1,1:3) = LN*nik(1)*nik(1:3) + Fabric_LN_d(1,1:3)
    Fabric_LN_d(2,1:3) = LN*nik(2)*nik(1:3) + Fabric_LN_d(2,1:3)
    Fabric_LN_d(3,1:3) = LN*nik(3)*nik(1:3) + Fabric_LN_d(3,1:3)

    Fabric_FN_d(1,1:3) = FN*nik(1)*nik(1:3) + Fabric_FN_d(1,1:3)
    Fabric_FN_d(2,1:3) = FN*nik(2)*nik(1:3) + Fabric_FN_d(2,1:3)
    Fabric_FN_d(3,1:3) = FN*nik(3)*nik(1:3) + Fabric_FN_d(3,1:3)
  
    Fabric_FT_d(1,1:3) = FT*nik(1)*tik(1:3) + Fabric_FT_d(1,1:3)
    Fabric_FT_d(2,1:3) = FT*nik(2)*tik(1:3) + Fabric_FT_d(2,1:3)
    Fabric_FT_d(3,1:3) = FT*nik(3)*tik(1:3) + Fabric_FT_d(3,1:3)
  else if (TAB_CONTACT_POLYR(i)%type > 2) then
    Fabric_t(1,1:3) = nik(1)*nik(1:3) + Fabric_t(1,1:3)
    Fabric_t(2,1:3) = nik(2)*nik(1:3) + Fabric_t(2,1:3)
    Fabric_t(3,1:3) = nik(3)*nik(1:3) + Fabric_t(3,1:3)

    Fabric_LN_t(1,1:3) = LN*nik(1)*nik(1:3) + Fabric_LN_t(1,1:3)
    Fabric_LN_t(2,1:3) = LN*nik(2)*nik(1:3) + Fabric_LN_t(2,1:3)
    Fabric_LN_t(3,1:3) = LN*nik(3)*nik(1:3) + Fabric_LN_t(3,1:3)

    Fabric_FN_t(1,1:3) = FN*nik(1)*nik(1:3) + Fabric_FN_t(1,1:3)
    Fabric_FN_t(2,1:3) = FN*nik(2)*nik(1:3) + Fabric_FN_t(2,1:3)
    Fabric_FN_t(3,1:3) = FN*nik(3)*nik(1:3) + Fabric_FN_t(3,1:3)
  
    Fabric_FT_t(1,1:3) = FT*nik(1)*tik(1:3) + Fabric_FT_t(1,1:3)
    Fabric_FT_t(2,1:3) = FT*nik(2)*tik(1:3) + Fabric_FT_t(2,1:3)
    Fabric_FT_t(3,1:3) = FT*nik(3)*tik(1:3) + Fabric_FT_t(3,1:3)
  end if 

  cpt = cpt + 1
end do
Norm_L = Norm_L /cpt
Norm_F = Norm_F /cpt

Fabric      = Fabric(:,:) / cpt
Fabric_LN   = Fabric_LN   / (Norm_L * cpt)
Fabric_FN   = Fabric_FN    / (Norm_F * cpt)
Fabric_FT   = Fabric_FT    / (Norm_F * cpt)
Fabric_F    = Fabric_FT + Fabric_FN

Fabric_s      = Fabric_s(:,:) / cpt
Fabric_LN_s   = Fabric_LN_s   / (Norm_L * cpt)
Fabric_FN_s   = Fabric_FN_s    / (Norm_F * cpt)
Fabric_FT_s   = Fabric_FT_s    / (Norm_F * cpt)
Fabric_F_s    = Fabric_FT_s + Fabric_FN_s

Fabric_ss      = Fabric_ss(:,:) / cpt
Fabric_LN_ss   = Fabric_LN_ss   / (Norm_L * cpt)
Fabric_FN_ss   = Fabric_FN_ss    / (Norm_F * cpt)
Fabric_FT_ss   = Fabric_FT_ss    / (Norm_F * cpt)
Fabric_F_ss    = Fabric_FT_ss + Fabric_FN_ss

Fabric_d      = Fabric_d(:,:) / cpt
Fabric_LN_d   = Fabric_LN_d   / (Norm_L * cpt)
Fabric_FN_d   = Fabric_FN_d    / (Norm_F * cpt)
Fabric_FT_d   = Fabric_FT_d    / (Norm_F * cpt)
Fabric_F_d    = Fabric_FT_d + Fabric_FN_d

Fabric_t      = Fabric_t(:,:) / cpt
Fabric_LN_t   = Fabric_LN_t   / (Norm_L * cpt)
Fabric_FN_t   = Fabric_FN_t    / (Norm_F * cpt)
Fabric_FT_t   = Fabric_FT_t    / (Norm_F * cpt)
Fabric_F_t    = Fabric_FT_t + Fabric_FN_t


!===== Tout les contacts
a = 0.D0 ; aln = 0.D0  ; afn = 0.D0 ; aft = 0.D0

lda  = 3
matz = 1

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
a = 2*(S3-S1)/(S1+S3)

M = Fabric_LN
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
aln = 2*(S3-S1)/(S1+S3) - a

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
afn = 2*(S3-S1)/(S1+S3) - a

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
aft = 2*(S3-S1)/(S1+S3) - a - afn

!as = 0.D0 ; alns = 0.D0  ; afns = 0.D0 ; afts = 0.D0
!
!as(1,3) = 5*( Fabric_s(3,3) - Fabric_s(1,1) ) /2
!alns(1,3) = 5*( Fabric_LN_s(3,3) - Fabric_LN_s(1,1) ) / (2*( Fabric_LN(3,3) + Fabric_LN(2,2) + Fabric_LN(1,1) )) - as(1,3)
!afns(1,3) = 5*( Fabric_FN_s(3,3) - Fabric_FN_s(1,1) ) /(2*( Fabric_FN(3,3) + Fabric_FN(2,2) + Fabric_FN(1,1) )) - as(1,3)
!afts(1,3) = 5*( Fabric_F_s(3,3) - Fabric_F_s(1,1) ) /(2*( Fabric_F(3,3) + Fabric_F(2,2) + Fabric_F(1,1) )) - afns(1,3) - as(1,3)

!ass(1,3) = 5*( Fabric_ss(3,3) - Fabric_ss(1,1) ) /2
!alnss(1,3) = 5*( Fabric_LN_ss(3,3) - Fabric_LN_ss(1,1) ) / (2*( Fabric_LN(3,3) + Fabric_LN(2,2) + Fabric_LN(1,1) )) - ass(1,3)
!afnss(1,3) = 5*( Fabric_FN_ss(3,3) - Fabric_FN_ss(1,1) ) /(2*( Fabric_FN(3,3) + Fabric_FN(2,2) + Fabric_FN(1,1) )) - ass(1,3)
!aftss(1,3) = 5*( Fabric_F_ss(3,3) - Fabric_F_ss(1,1) ) /(2*( Fabric_F(3,3) + Fabric_F(2,2) + Fabric_F(1,1) )) - afnss(1,3) - ass(1,3)

!ad(1,3) = 5*( Fabric_d(3,3) - Fabric_d(1,1) ) /2
!alnd(1,3) = 5*( Fabric_LN_d(3,3) - Fabric_LN_d(1,1) ) / (2*( Fabric_LN(3,3) + Fabric_LN(2,2) + Fabric_LN(1,1) )) - ad(1,3)
!afnd(1,3) = 5*( Fabric_FN_d(3,3) - Fabric_FN_d(1,1) ) /(2*( Fabric_FN(3,3) + Fabric_FN(2,2) + Fabric_FN(1,1) )) - ad(1,3)
!aftd(1,3) = 5*( Fabric_F_d(3,3) - Fabric_F_d(1,1) ) /(2*( Fabric_F(3,3) + Fabric_F(2,2) + Fabric_F(1,1) )) - afnd(1,3) - ad(1,3)

!at(1,3) = 5*( Fabric_t(3,3) - Fabric_t(1,1) ) /2
!alnt(1,3) = 5*( Fabric_LN_t(3,3) - Fabric_LN_t(1,1) ) / (2*( Fabric_LN(3,3) + Fabric_LN(2,2) + Fabric_LN(1,1) )) - at(1,3)
!afnt(1,3) = 5*( Fabric_FN_t(3,3) - Fabric_FN_t(1,1) ) /(2*( Fabric_FN(3,3) + Fabric_FN(2,2) + Fabric_FN(1,1) )) - at(1,3)
!aftt(1,3) = 5*( Fabric_F_t(3,3) - Fabric_F_t(1,1) ) /(2*( Fabric_F(3,3) + Fabric_F(2,2) + Fabric_F(1,1) )) - afnt(1,3) - at(1,3)

write(105,'(40(1X,D12.5))') temps,Def,Hauteur,a,aln,afn,aft,(a+aln+afn+aft)/2

end subroutine branche_anisotropies


!==============================================================================
! Calcul des differentes anisotropies
!==============================================================================
subroutine contact_anisotropies
implicit none
real(kind=8),dimension(3,3)              :: Fabric,Fabric_LN,Fabric_LT,Fabric_L,M
real(kind=8),dimension(3,3)              :: Fabric_FN,Fabric_FT,Fabric_F

real(kind=8),dimension(3,3)              :: Fabric_s,Fabric_LN_s,Fabric_LT_s,Fabric_L_s
real(kind=8),dimension(3,3)              :: Fabric_FN_s,Fabric_FT_s,Fabric_F_s

real(kind=8),dimension(3,3)              :: Fabric_d,Fabric_LN_d,Fabric_LT_d,Fabric_L_d
real(kind=8),dimension(3,3)              :: Fabric_FN_d,Fabric_FT_d,Fabric_F_d

real(kind=8),dimension(3,3)              :: Fabric_t,Fabric_LN_t,Fabric_LT_t,Fabric_L_t
real(kind=8),dimension(3,3)              :: Fabric_FN_t,Fabric_FT_t,Fabric_F_t

real(kind=8),dimension(3)                :: nik,tik,sik,Lik,Fik
real(kind=8)                             :: Rtik,Rnik,Rsik,norm,FN,FT,LN,LT,Norm_F,Norm_L,cpt,Rayon_max
integer                                  :: i,cd,an
real(kind=8)                             :: a,aln,alt,afn,aft,&
                                            as,alns,alts,afns,afts,&
                                            ad,alnd,altd,afnd,aftd,&
                                            at,alnt,altt,afnt,aftt

real(kind=8)                             :: S1,S2,S3
real(kind=8),dimension(3,3)              :: localframe  
integer                                  :: ierror,matz,lda
real(kind=8),dimension(3)                :: wr,wi

real(kind=8)                             :: cross_product13_1,cross_product13_2,cross_product13
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
do i=1,nb_ligneCONTACT_POLYR
  if (TAB_CONTACT_POLYR(i)%rn<0.000000000001*Mean_total_Normal_force) cycle
  cd   = TAB_CONTACT_POLYR(i)%icdent
  an   = TAB_CONTACT_POLYR(i)%ianent
  nik  = TAB_CONTACT_POLYR(i)%n
  tik  = TAB_CONTACT_POLYR(i)%t
  sik  = TAB_CONTACT_POLYR(i)%s
  Rtik = TAB_CONTACT_POLYR(i)%rt
  Rnik = TAB_CONTACT_POLYR(i)%rn
  Rsik = TAB_CONTACT_POLYR(i)%rs

  if ((TAB_POLYEDRE(cd)%behav/='plexx').or.(TAB_POLYEDRE(an)%behav/='plexx')) cycle

  Lik(1:3) = TAB_POLYEDRE(cd)%centre(1:3)-TAB_POLYEDRE(an)%centre(1:3)
  if ( sqrt( Lik(1)**2 + Lik(2)**2 + Lik(3)**2) > 2*( TAB_POLYEDRE(cd)%Rmax + TAB_POLYEDRE(an)%Rmax ) ) then
    Rayon_max = sqrt( ( Lik(1)**2 + Lik(2)**2 + Lik(3)**2) ) 
    Lik(:) = Lik(:) / Rayon_max
    Lik(:) = abs(TAB_POLYEDRE(cd)%Rmax+TAB_POLYEDRE(an)%Rmax)*Lik(:)
  end if

  Fik(1:3) = (Rnik*nik(1:3)+Rtik*tik(1:3)+Rsik*sik(1:3))
  Norm_F   = Norm_F + sqrt( Fik(1)**2 + Fik(2)**2 + Fik(3)**2 )
  Norm_L   = Norm_L + sqrt( Lik(1)**2 + Lik(2)**2 + Lik(3)**2 )
  LN       =  Lik(1)*nik(1) + Lik(2)*nik(2) + Lik(3)*nik(3)  
  FN       =  Fik(1)*nik(1) + Fik(2)*nik(2) + Fik(3)*nik(3)  

  tik(1:3) =  Fik(1:3) - FN*nik(1:3)
  
  FT       =  sqrt( tik(1)**2 + tik(2)**2 + tik(3)**2 )

  if ( FT > 0.00000000001) tik(1:3) = tik(1:3) / FT

  LT       =  Lik(1)*tik(1) + Lik(2)*tik(2) + Lik(3)*tik(3)  

  Fabric(1,1:3) = nik(1)*nik(1:3) + Fabric(1,1:3)
  Fabric(2,1:3) = nik(2)*nik(1:3) + Fabric(2,1:3)
  Fabric(3,1:3) = nik(3)*nik(1:3) + Fabric(3,1:3)

  Fabric_LN(1,1:3) = LN*nik(1)*nik(1:3) + Fabric_LN(1,1:3)
  Fabric_LN(2,1:3) = LN*nik(2)*nik(1:3) + Fabric_LN(2,1:3)
  Fabric_LN(3,1:3) = LN*nik(3)*nik(1:3) + Fabric_LN(3,1:3)

  Fabric_LT(1,1:3) = LT*nik(1)*tik(1:3) + Fabric_LT(1,1:3)
  Fabric_LT(2,1:3) = LT*nik(2)*tik(1:3) + Fabric_LT(2,1:3)
  Fabric_LT(3,1:3) = LT*nik(3)*tik(1:3) + Fabric_LT(3,1:3)

  Fabric_FN(1,1:3) = FN*nik(1)*nik(1:3) + Fabric_FN(1,1:3)
  Fabric_FN(2,1:3) = FN*nik(2)*nik(1:3) + Fabric_FN(2,1:3)
  Fabric_FN(3,1:3) = FN*nik(3)*nik(1:3) + Fabric_FN(3,1:3)
  
  Fabric_FT(1,1:3) = FT*nik(1)*tik(1:3) + Fabric_FT(1,1:3)
  Fabric_FT(2,1:3) = FT*nik(2)*tik(1:3) + Fabric_FT(2,1:3)
  Fabric_FT(3,1:3) = FT*nik(3)*tik(1:3) + Fabric_FT(3,1:3)

  cpt = cpt + 1
end do
Norm_L = Norm_L /cpt
Norm_F = Norm_F /cpt

Fabric      = Fabric(:,:) / cpt
Fabric_LN   = Fabric_LN   / (Norm_L * cpt)
Fabric_LT   = Fabric_LT   / (Norm_L * cpt)
Fabric_L    = Fabric_LT + Fabric_LN
Fabric_FN   = Fabric_FN    / (Norm_F * cpt)
Fabric_FT   = Fabric_FT    / (Norm_F * cpt)
Fabric_F    = Fabric_FT + Fabric_FN

!===== Tout les contacts
a = 0.D0 ; aln = 0.D0  ; afn = 0.D0 ; aft = 0.D0

lda  = 3
matz = 1

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
a = 2*(S3-S1)/(S1+S3)

M = Fabric_LN
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
aln = 2*(S3-S1)/(S1+S3) - a

M = Fabric_L
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
alt = 2*(S3-S1)/(S1+S3) - a - aln

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
afn = 2*(S3-S1)/(S1+S3) - a

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
aft = 2*(S3-S1)/(S1+S3) - a - afn


if (type_box==0)&
write(104,'(31(1X,D12.5))') temps,Def,Hauteur,a,aln,alt,afn,aft,(a+aln+alt+afn+aft)/2

end subroutine contact_anisotropies



!==============================================================================
! Calcul de la compacity dans un �chantillon
!==============================================================================
subroutine compacity
implicit none
integer                                  ::  i
real(kind=8)                             ::  Volume


print*,'.... Compacity ....'

Volume = 0.D0
do i=1,nbParticules
  if (TAB_POLYEDRE(i)%behav/='plexx') cycle
  Volume = Volume + TAB_POLYEDRE(i)%volume
end do

if (type_box==0) write(103,'(7(1X,D12.5))') temps,Def,Hauteur,Volume/(longueur*largeur*hauteur)

end subroutine compacity

!==============================================================================
! Calcul de q/p
!==============================================================================
subroutine qsurp
implicit none
real(kind=8),dimension(3,3)              :: Moment,Moment_s,Moment_ss,Moment_d,Moment_t,M
real(kind=8),dimension(3)                :: nik,tik,sik,Lik,Fik
real(kind=8)                             :: Rtik,Rnik,Rsik,q,p,qs,qd,qt,qss,Rayon_max,tan_phi_contrainte
integer                                  :: i,cd,an
real(kind=8)                             :: S1,S2,S3
real(kind=8)                             :: S1s,S2s,S3s
real(kind=8)                             :: S1ss,S2ss,S3ss
real(kind=8)                             :: S1d,S2d,S3d
real(kind=8)                             :: S1t,S2t,S3t
real(kind=8),dimension(3,3)              :: localframe  
integer                                  :: ierror,matz,lda
real(kind=8),dimension(3)                :: wr,wi

print*,'.... q/p ....'

Moment(:,:)   = 0
Moment_s(:,:) = 0
Moment_ss(:,:) = 0
Moment_d(:,:) = 0
Moment_t(:,:) = 0
Lik = 0
Fik = 0
do i=1,nb_ligneCONTACT_POLYR
  if (TAB_CONTACT_POLYR(i)%rn<0.000000000001*Mean_total_Normal_force) cycle
  cd   = TAB_CONTACT_POLYR(i)%icdent
  an   = TAB_CONTACT_POLYR(i)%ianent
  nik  = TAB_CONTACT_POLYR(i)%n
  tik  = TAB_CONTACT_POLYR(i)%t
  sik  = TAB_CONTACT_POLYR(i)%s
  Rtik = TAB_CONTACT_POLYR(i)%rt
  Rnik = TAB_CONTACT_POLYR(i)%rn
  Rsik = TAB_CONTACT_POLYR(i)%rs

  if ((TAB_POLYEDRE(cd)%behav/='plexx').or.(TAB_POLYEDRE(an)%behav/='plexx')) cycle

  Lik(1:3) = TAB_POLYEDRE(cd)%centre(1:3)-TAB_POLYEDRE(an)%centre(1:3)
  if ( sqrt( Lik(1)**2 + Lik(2)**2 + Lik(3)**2) > 2*( TAB_POLYEDRE(cd)%Rmax + TAB_POLYEDRE(an)%Rmax ) ) then
    Rayon_max = sqrt( ( Lik(1)**2 + Lik(2)**2 + Lik(3)**2) ) 
    Lik(:) = Lik(:) / Rayon_max
    Lik(:) = abs(TAB_POLYEDRE(cd)%Rmax+TAB_POLYEDRE(an)%Rmax)*Lik(:)
  end if

  Fik(1:3) = (Rnik*nik(1:3)+Rtik*tik(1:3)+Rsik*sik(1:3))

  Moment(1,1:3) = Fik(1)*Lik(1:3) + Moment(1,1:3)
  Moment(2,1:3) = Fik(2)*Lik(1:3) + Moment(2,1:3)
  Moment(3,1:3) = Fik(3)*Lik(1:3) + Moment(3,1:3)

  if (TAB_CONTACT_POLYR(i)%type == 0) then
    Moment_ss(1,1:3) = Fik(1)*Lik(1:3) + Moment_ss(1,1:3)
    Moment_ss(2,1:3) = Fik(2)*Lik(1:3) + Moment_ss(2,1:3)
    Moment_ss(3,1:3) = Fik(3)*Lik(1:3) + Moment_ss(3,1:3)
  else if (TAB_CONTACT_POLYR(i)%type == 1) then
    Moment_s(1,1:3) = Fik(1)*Lik(1:3) + Moment_s(1,1:3)
    Moment_s(2,1:3) = Fik(2)*Lik(1:3) + Moment_s(2,1:3)
    Moment_s(3,1:3) = Fik(3)*Lik(1:3) + Moment_s(3,1:3)
  else if (TAB_CONTACT_POLYR(i)%type == 2) then
    Moment_d(1,1:3) = Fik(1)*Lik(1:3) + Moment_d(1,1:3)
    Moment_d(2,1:3) = Fik(2)*Lik(1:3) + Moment_d(2,1:3)
    Moment_d(3,1:3) = Fik(3)*Lik(1:3) + Moment_d(3,1:3)
  else if (TAB_CONTACT_POLYR(i)%type > 2) then
    Moment_t(1,1:3) = Fik(1)*Lik(1:3) + Moment_t(1,1:3)
    Moment_t(2,1:3) = Fik(2)*Lik(1:3) + Moment_t(2,1:3)
    Moment_t(3,1:3) = Fik(3)*Lik(1:3) + Moment_t(3,1:3)
  end if 
end do

tan_phi_contrainte = -Moment(2,3)/ Moment(3,3)

!===== Tout les contacts
lda  = 3
matz = 1

M = Moment
pressure_in_sample = abs(M(1,1)+M(2,2)+M(3,3)) / 3
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
p = (S1+S3)/2

pressure_in_sample = abs(S1+S2+S3) / 3


!===== contacts ss
M = Moment_ss
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
qss = (S3-S1)/2

!===== contacts sv
M = Moment_s
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
qs = (S3-S1)/2

!===== contacts sf
M = Moment_d
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
qd = (S3-S1)/2

!===== contacts ff
M = Moment_t
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
qt = (S3-S1)/2

if (type_box==0) write(102,'(27(1X,D12.5))') temps,Def,Hauteur,tan_phi_contrainte,&
                                             q/p ,qss/p,qs/p,qd/p,qt/p,qss/p+qs/p+qd/p+qt/p
                                        
end subroutine qsurp

!==============================================================================
! Vitesse moyenne
!==============================================================================
subroutine vitesse_moyenne
implicit none
real(kind=8)                             :: V,Vx,Vy,Vmax,cpt,el_diameter
integer                                  :: i
real(kind=8)                             :: Norm_V

print*,'.... Vitesse moyenne ....'

V   = 0.D0
Vx  = 0.D0
Vy  = 0.D0
Vmax  = 0.D0
cpt = 0.D0
el_diameter = 0
do i=1,nbParticules
    if (TAB_POLYEDRE(i)%behav /= 'plexx')  cycle
    V    = V  + TAB_POLYEDRE(i)%Norme_V
    Vx   = Vx + TAB_POLYEDRE(i)%Vx
    Vy   = Vy + TAB_POLYEDRE(i)%Vy
    Vmax = max( Vmax, TAB_POLYEDRE(i)%Norme_V ) 
    el_diameter = el_diameter + 2*TAB_POLYEDRE(i)%Rmax
    cpt  = cpt +1
end do

V  = V/cpt
Vx = Vx/cpt
Vy = Vy/cpt
Vmax = Vmax/cpt
el_diameter = el_diameter / cpt

print*, 'Taux de deformation \dot gamma = ',  Vy/Hauteur
print*, 'Pressure in sample             = ',  pressure_in_sample
print*, 'Inertie                        = ',  (0.05/Hauteur) * el_diameter * sqrt(2800/pressure_in_sample)
print*, 'dmoyen                         = ',  el_diameter 

if (type_box==0) write(101,'(10(1X,D12.5))') temps,Def,Hauteur,V,Vx,Vy,Vmax,&
                                             Vy/Hauteur,&
                                             (0.05/Hauteur) * el_diameter * sqrt(2800/pressure_in_sample)

end subroutine vitesse_moyenne

!==============================================================================
! Nombre de coordination
!==============================================================================
subroutine nb_coordination
implicit none
real(kind=8)                             :: z,N_c,Np_total,Np_connecteted,ks,kss,kd,kt,kmobilized,fmobilized,z_contact=0,&
                                            P0,P1,P2,P3,P4,P5,P6,P7,P8,P9
integer                                  :: i,j,cd,an
real(kind=8),dimension(:),allocatable    :: ctc_part

print*,'.... Nombre coordination ....'

z              = 0
Np_total       = 0
Np_connecteted = 0
N_c            = 0
ks             = 0
kss            = 0
kd             = 0
kt             = 0
kmobilized     = 0
fmobilized     = 0
z_contact      = 0
P0 = 0 ; P1 = 0 ; P2 = 0 ; P3 = 0 ; P4 = 0 ; P5 = 0 ; P6 = 0 ; P7 = 0 ; P8 = 0 ; P9 = 0

do i=1,nbParticules
  TAB_POLYEDRE(i)%nctc = 0
end do

do i= 1,nb_ligneCONTACT_POLYR
  if (TAB_CONTACT_POLYR(i)%rn==0.D0) cycle
  cd = TAB_CONTACT_POLYR(i)%icdent
  an = TAB_CONTACT_POLYR(i)%ianent
  if ( (TAB_POLYEDRE(cd)%behav /= 'plexx').and.(TAB_POLYEDRE(an)%behav /= 'plexx') ) cycle
  if (TAB_CONTACT_POLYR(i)%type == 0) kss = kss + 1.D0
  if (TAB_CONTACT_POLYR(i)%type == 1) ks  = ks  + 1.D0
  if (TAB_CONTACT_POLYR(i)%type == 2) kd  = kd  + 1.D0
  if (TAB_CONTACT_POLYR(i)%type >  2) kt  = kt  + 1.D0
  if ( abs( sqrt(TAB_CONTACT_POLYR(i)%rt**2+TAB_CONTACT_POLYR(i)%rs**2) ) > 0.39 * TAB_CONTACT_POLYR(i)%rn ) kmobilized = kmobilized+1.D0
  N_c = N_c+1.D0
  fmobilized = fmobilized + abs( sqrt(TAB_CONTACT_POLYR(i)%rt**2+TAB_CONTACT_POLYR(i)%rs**2) ) / ( 0.4 * TAB_CONTACT_POLYR(i)%rn )
end do


do i=1,nb_ligneCONTACT_POLYR
  if (TAB_CONTACT_POLYR(i)%nature=='PRPLx') cycle
  if (TAB_CONTACT_POLYR(i)%rn==0.D0) cycle
  cd = TAB_CONTACT_POLYR(i)%icdent
  an = TAB_CONTACT_POLYR(i)%ianent
  if ( (TAB_POLYEDRE(cd)%behav /= 'plexx').and.(TAB_POLYEDRE(an)%behav /= 'plexx') ) then
    cycle
  else if ( (TAB_POLYEDRE(cd)%behav == 'plexx').and.(TAB_POLYEDRE(an)%behav == 'plexx') ) then
    TAB_POLYEDRE(cd)%nctc = TAB_POLYEDRE(cd)%nctc + 1
    TAB_POLYEDRE(an)%nctc = TAB_POLYEDRE(an)%nctc + 1
  else  if ( (TAB_POLYEDRE(an)%behav == 'plexx').and.(TAB_POLYEDRE(cd)%behav /= 'plexx') ) then
    TAB_POLYEDRE(an)%nctc = TAB_POLYEDRE(an)%nctc + 1
  else  if ( (TAB_POLYEDRE(cd)%behav == 'plexx').and.(TAB_POLYEDRE(an)%behav /= 'plexx') ) then
    TAB_POLYEDRE(cd)%nctc = TAB_POLYEDRE(cd)%nctc + 1
  end if
end do
    
do i=1,nbParticules
  if (TAB_POLYEDRE(i)%behav /= 'plexx')  cycle

  Np_total=Np_total+1.D0
  
  if (TAB_POLYEDRE(i)%nctc > 0)  Np_connecteted = Np_connecteted + 1
  if (TAB_POLYEDRE(i)%nctc == 0) P0 = P0+1
  if (TAB_POLYEDRE(i)%nctc == 1) P1 = P1+1
  if (TAB_POLYEDRE(i)%nctc == 2) P2 = P2+1
  if (TAB_POLYEDRE(i)%nctc == 3) P3 = P3+1
  if (TAB_POLYEDRE(i)%nctc == 4) P4 = P4+1
  if (TAB_POLYEDRE(i)%nctc == 5) P5 = P5+1
  if (TAB_POLYEDRE(i)%nctc == 6) P6 = P6+1
  if (TAB_POLYEDRE(i)%nctc == 7) P7 = P7+1
  if (TAB_POLYEDRE(i)%nctc == 8) P8 = P8+1
  if (TAB_POLYEDRE(i)%nctc == 9) P9 = P9+1
    
end do

P0 = P0 / Np_total

P1 = P1 / Np_connecteted
P2 = P2 / Np_connecteted
P3 = P3 / Np_connecteted
P4 = P4 / Np_connecteted
P5 = P5 / Np_connecteted
P6 = P6 / Np_connecteted
P7 = P7 / Np_connecteted
P8 = P8 / Np_connecteted
P9 = P9 / Np_connecteted

ks = ks / N_c
kss = kss / N_c
kd = kd / N_c
kt = kt / N_c

z         = P1 + 2*P2 + 3*P3 + 4*P4 + 5*P5 + 6*P6 + 7*P7 + 8*P8 + 9*P9
z_contact = (kss+ks+2*kd+3*kt) * z

kmobilized = kmobilized / N_c
fmobilized = fmobilized / N_c


if (type_box==0) write(100,'(50(1X,D12.5))') temps,Def,Hauteur,z,z_contact,kss,ks,kd,kt,kmobilized,fmobilized,&
                                             P0,P1,P2,P3,P4,P5,P6,P7,P8,P9
end subroutine nb_coordination

!==================================================================================================
!Procedure pour ecrire le gmv a chaque pas de temps
!==================================================================================================
subroutine write_all_gmv
implicit none
character(len=13)                           ::  nom
integer                                     ::  i,j,l,nb_nodev,ncells,cpt_force
integer                                     ::  k=0,cd,an
character(len=5),dimension(:),allocatable   ::  list_bodie_material,list_material
integer                                     ::  nb_bodie_material,GMV_nb_materials
real(kind=8)                                ::  Ax2,fmax,fmean,fmin,force_contact
real(kind=8),dimension(3)                   ::  Lik


if ( (step_all_gmv == compteur_clout).or. (step_all_gmv == 0) ) then

  print*,'....  WRITE GMV FILE....'

  nom       =  'GMV_FILE.0000'
  nb_nodev  = 0
  ncells    = 0
  if (compteur_clout<10) then
    WRITE(nom(12:13),'(I1)')   compteur_clout
  else if ( (compteur_clout>=10) .and. (compteur_clout<100) ) then
    WRITE(nom(11:13),'(I2)')   compteur_clout
  else if ( (compteur_clout>=100).and. (compteur_clout<1000) ) then
    WRITE(nom(10:13),'(I3)')   compteur_clout
  end if   

  fmax           = 0
  fmean          = 0
  fmin           = 100000
  cpt_force      = 0
  print*, 'Mean total force = ',Mean_total_Normal_force
  do i=1,nb_ligneCONTACT_POLYR 
    TAB_CONTACT_POLYR(i)%cycle_gmv = .false.

    if  (TAB_CONTACT_POLYR(i)%nature=='PRPLx') cycle
    if  (TAB_CONTACT_POLYR(i)%rn<0.000000000001*Mean_total_Normal_force) cycle
    
    cd       =  TAB_CONTACT_POLYR(i)%icdent
    an       =  TAB_CONTACT_POLYR(i)%ianent
    Lik(1:3) =  TAB_POLYEDRE(cd)%centre(1:3)-TAB_POLYEDRE(an)%centre(1:3)
    if ((TAB_POLYEDRE(cd)%behav/='plexx').or.(TAB_POLYEDRE(an)%behav/='plexx')) cycle

    if ( sqrt( Lik(1)**2 + Lik(2)**2 + Lik(3)**2) > 2*( TAB_POLYEDRE(cd)%Rmax + TAB_POLYEDRE(an)%Rmax ) ) then
      TAB_CONTACT_POLYR(i)%cycle_gmv = .true.
      cycle
    end if
    
    force_contact = sqrt(TAB_CONTACT_POLYR(i)%rn**2 + TAB_CONTACT_POLYR(i)%rs**2 + TAB_CONTACT_POLYR(i)%rt**2 )  

    cpt_force = cpt_force + 1
    fmax   = max(fmax,force_contact)
    fmin   = min(fmin,force_contact)
    fmean  = fmean + force_contact
    Lik(1:3) = TAB_POLYEDRE(cd)%centre(1:3)-TAB_POLYEDRE(an)%centre(1:3)
  end do
  fmean = fmean / real(cpt_force,8)
print*, 'CPT FORCE = ',cpt_force
  do i=1,nbParticules
    nb_nodev = nb_nodev + TAB_POLYEDRE(i)%nb_vertex
  end do

  open(unit=30000,file=nom,status='replace') 

  write(30000,'(A14)')      'gmvinput ascii'
  !--------------------------------------- Nodev
  write(30000,'(A5,1X,I9)') 'nodev',(nb_nodev) + (cpt_force)*8

  do i=1,nbParticules
    do j=1,TAB_POLYEDRE(i)%nb_vertex
      write(30000,'(3(1X,E14.7))')  TAB_POLYEDRE(i)%vertex(1,j),&
                                    TAB_POLYEDRE(i)%vertex(2,j),&
                                    TAB_POLYEDRE(i)%vertex(3,j)
    end do
  enddo 

  Ax2 = 0.D0
  cd  = 0
  an  = 0
  do i=1,nb_ligneCONTACT_POLYR
    if  (TAB_CONTACT_POLYR(i)%nature=='PRPLx') cycle
    if  (TAB_CONTACT_POLYR(i)%cycle_gmv) cycle
    if  (TAB_CONTACT_POLYR(i)%rn<0.000000000001*Mean_total_Normal_force) cycle

    cd = TAB_CONTACT_POLYR(i)%icdent
    an = TAB_CONTACT_POLYR(i)%ianent
    if ((TAB_POLYEDRE(cd)%behav/='plexx').or.(TAB_POLYEDRE(an)%behav/='plexx')) cycle

    force_contact = sqrt(TAB_CONTACT_POLYR(i)%rn**2 + TAB_CONTACT_POLYR(i)%rs**2 + TAB_CONTACT_POLYR(i)%rt**2 )  
    Ax2 = 0.75 * (TAB_POLYEDRE(cd)%Rmax + TAB_POLYEDRE(an)%Rmax) *( (force_contact - fmin)/(fmax - fmin) )
     
    write(30000,'(3(1X,E14.7))')  TAB_POLYEDRE(cd)%centre(1)-Ax2,&
                                  TAB_POLYEDRE(cd)%centre(2)-Ax2,&
                                  TAB_POLYEDRE(cd)%centre(3)
    write(30000,'(3(1X,E14.7))')  TAB_POLYEDRE(cd)%centre(1)+Ax2,&
                                  TAB_POLYEDRE(cd)%centre(2)-Ax2,&
                                  TAB_POLYEDRE(cd)%centre(3)
    write(30000,'(3(1X,E14.7))')  TAB_POLYEDRE(cd)%centre(1)+Ax2,&
                                  TAB_POLYEDRE(cd)%centre(2)+Ax2,&
                                  TAB_POLYEDRE(cd)%centre(3)
    write(30000,'(3(1X,E14.7))')  TAB_POLYEDRE(cd)%centre(1)-Ax2,&
                                  TAB_POLYEDRE(cd)%centre(2)+Ax2,&
                                  TAB_POLYEDRE(cd)%centre(3)
    write(30000,'(3(1X,E14.7))')  TAB_POLYEDRE(an)%centre(1)-Ax2,&
                                  TAB_POLYEDRE(an)%centre(2)-Ax2,&
                                  TAB_POLYEDRE(an)%centre(3)
    write(30000,'(3(1X,E14.7))')  TAB_POLYEDRE(an)%centre(1)+Ax2,&
                                  TAB_POLYEDRE(an)%centre(2)-Ax2,&
                                  TAB_POLYEDRE(an)%centre(3)
    write(30000,'(3(1X,E14.7))')  TAB_POLYEDRE(an)%centre(1)+Ax2,&
                                  TAB_POLYEDRE(an)%centre(2)+Ax2,&
                                  TAB_POLYEDRE(an)%centre(3)
    write(30000,'(3(1X,E14.7))')  TAB_POLYEDRE(an)%centre(1)-Ax2,&
                                  TAB_POLYEDRE(an)%centre(2)+Ax2,&
                                  TAB_POLYEDRE(an)%centre(3)
  end do

  !--------------------------------------- Cells
  do i=1,nbParticules
    ncells = ncells + TAB_POLYEDRE(i)%nb_face
  end do
  write(30000,'(A5,1X,I9)') 'cells',ncells + 6*(cpt_force)

  !--------------------------------------- Surface
  k = 0
  do i=1,nbParticules
    if ( i==1 ) k = 0
    if ( i>1 )  k = k + TAB_POLYEDRE(i-1)%nb_vertex
    do j=1,TAB_POLYEDRE(i)%nb_face
      write(30000,'(A3,1X,I7)') 'tri',3 
      write(30000,'(3(I11))')    TAB_POLYEDRE(i)%face(1,j)+k,&
                                 TAB_POLYEDRE(i)%face(2,j)+k,&
                                 TAB_POLYEDRE(i)%face(3,j)+k
    end do
  end do

  k = k + TAB_POLYEDRE(nbParticules)%nb_vertex 
  do i=1,nb_ligneCONTACT_POLYR
    if  (TAB_CONTACT_POLYR(i)%nature=='PRPLx') cycle
    if  (TAB_CONTACT_POLYR(i)%cycle_gmv) cycle
    if  (TAB_CONTACT_POLYR(i)%rn<0.000000000001*Mean_total_Normal_force) cycle
    cd = TAB_CONTACT_POLYR(i)%icdent
    an = TAB_CONTACT_POLYR(i)%ianent
    if ((TAB_POLYEDRE(cd)%behav/='plexx').or.(TAB_POLYEDRE(an)%behav/='plexx')) cycle

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

  !--------------------------------------- Variables
  write(30000,'(A8)')     'variable'
  !--------------------------------------- Vitesse moyenne
  write(30000,'(A5,4X,I2)')  '<V>  ',1
  do i=1,nbParticules
    write(30000,'(8(E14.7,1X))') &
    ((sqrt( TAB_POLYEDRE(i)%Vx**2 + TAB_POLYEDRE(i)%Vy**2 + TAB_POLYEDRE(i)%Vz**2)),j=1,TAB_POLYEDRE(i)%nb_vertex) 
  end do
  write(30000,'(A5,4X,I2)')  '<Vr> ',1
  do i=1,nbParticules
    write(30000,'(8(E14.7,1X))') &
    ((sqrt( TAB_POLYEDRE(i)%Vrx**2 + TAB_POLYEDRE(i)%Vry**2 + TAB_POLYEDRE(i)%Vrz**2)),j=1,TAB_POLYEDRE(i)%nb_vertex) 
  end do
  !--------------------------------------- Pression
  write(30000,'(A5,4X,I2)')  '<P>  ',1
  do i=1,nb_ligneCONTACT_POLYR
    if (TAB_CONTACT_POLYR(i)%nature == 'PRPLx') cycle
    if (TAB_CONTACT_POLYR(i)%rn==0.D0) cycle
    cd = TAB_CONTACT_POLYR(i)%icdent
    an = TAB_CONTACT_POLYR(i)%ianent
    TAB_POLYEDRE(cd)%Force = sqrt(TAB_CONTACT_POLYR(i)%rn**2 + TAB_CONTACT_POLYR(i)%rt**2 + TAB_CONTACT_POLYR(i)%rs**2)
    TAB_POLYEDRE(an)%Force = sqrt(TAB_CONTACT_POLYR(i)%rn**2 + TAB_CONTACT_POLYR(i)%rt**2 + TAB_CONTACT_POLYR(i)%rs**2)
    TAB_POLYEDRE(cd)%nb_ctc = TAB_POLYEDRE(cd)%nb_ctc + 1
    TAB_POLYEDRE(an)%nb_ctc = TAB_POLYEDRE(an)%nb_ctc + 1
  end do
  do i=1,nbParticules
    write(30000,'(8(E14.7,1X))') &
    (TAB_POLYEDRE(i)%Force,j=1,TAB_POLYEDRE(i)%nb_vertex) 
  end do
  !--------------------------------------- Deplacement suivant x
!!  write(30000,'(A5,4X,I2)')  'DepX ',1
!!  do i=1,nbParticules
!!    write(30000,'(8(E14.7,1X))') &
!!    ( abs(TAB_POLYEDRE(i)%centre(1)-TAB_POLYEDRE(i)%centre_ref(1)),j=1,TAB_POLYEDRE(i)%nb_vertex) 
!!  end do
!!  !--------------------------------------- Deplacement suivant y
!!  write(30000,'(A5,4X,I2)')  'DepY ',1
!!  do i=1,nbParticules
!!    write(30000,'(8(E14.7,1X))') &
!!    ( abs(TAB_POLYEDRE(i)%centre(2)-TAB_POLYEDRE(i)%centre_ref(2)),j=1,TAB_POLYEDRE(i)%nb_vertex) 
!!  end do

  write(30000,'(A7)') 'endvars'

  !--------------------------------------- Temps
  write(30000,'(A8,2X,E14.7)') 'probtime',temps

  !--------------------------------------- List Material
  k                 = 0
  GMV_nb_materials  = 0
  nb_bodie_material = 0
  if ( allocated(list_bodie_material) ) deallocate(list_bodie_material)
  allocate( list_bodie_material(20) )
  call construct_list_material( nb_bodie_material , list_bodie_material )

  !  !-------> + Force simple, double, triple
  GMV_nb_materials = GMV_nb_materials + 8
  k = k + 8
    !-------> + Behav grains  (TOUJOURS EN DERNIER !!)
  GMV_nb_materials = GMV_nb_materials + nb_bodie_material

  if ( allocated(list_material) ) deallocate(list_material)
  allocate( list_material(GMV_nb_materials) )

  !--------------------------------------- Liste material
  list_material(1) = 'F_sW'
  list_material(2) = 'F_ssW'
  list_material(3) = 'F_dW'
  list_material(4) = 'F_tW'
  list_material(5) = 'F_sS'
  list_material(6) = 'F_ssS'
  list_material(7) = 'F_dS'
  list_material(8) = 'F_tS'
  do i=k+1,GMV_nb_materials
    list_material(i) = list_bodie_material(i-k)
  end do
  write(30000,'(A8,I7,2X,I1)') 'material',GMV_nb_materials,0
  do i = 1,GMV_nb_materials
    write(30000,'(A5)') list_material(i) 
  end do

  !-------> + Behav grains
  do i=1,nbParticules
    do j=1,GMV_nb_materials
      if ( TAB_POLYEDRE(i)%behav == list_material(j) ) then
        TAB_POLYEDRE(i)%gmv_material = j
        exit
      end if
    end do  
    write(30000,'(10(I6))') (TAB_POLYEDRE(i)%gmv_material,j=1,TAB_POLYEDRE(i)%nb_face)
  end do

!  !-------> + Force simple, double triple
  do i=1,nb_ligneCONTACT_POLYR  
    if  (TAB_CONTACT_POLYR(i)%nature=='PRPLx') cycle
    if  (TAB_CONTACT_POLYR(i)%cycle_gmv) cycle
    if  (TAB_CONTACT_POLYR(i)%rn<0.000000000001*Mean_total_Normal_force) cycle
    cd = TAB_CONTACT_POLYR(i)%icdent
    an = TAB_CONTACT_POLYR(i)%ianent
    if ((TAB_POLYEDRE(cd)%behav/='plexx').or.(TAB_POLYEDRE(an)%behav/='plexx')) cycle

    force_contact = sqrt(TAB_CONTACT_POLYR(i)%rn**2 + TAB_CONTACT_POLYR(i)%rs**2 + TAB_CONTACT_POLYR(i)%rt**2 )  
    if ((TAB_CONTACT_POLYR(i)%type == 0).and.(force_contact/fmean < 1.D0)) write(30000,'(6(I6))') (1,j=1,6)
    if ((TAB_CONTACT_POLYR(i)%type == 1).and.(force_contact/fmean < 1.D0)) write(30000,'(6(I6))') (2,j=1,6)
    if ((TAB_CONTACT_POLYR(i)%type == 2).and.(force_contact/fmean < 1.D0)) write(30000,'(6(I6))') (3,j=1,6)
    if ((TAB_CONTACT_POLYR(i)%type >  2).and.(force_contact/fmean < 1.D0)) write(30000,'(6(I6))') (4,j=1,6)

    if ((TAB_CONTACT_POLYR(i)%type == 0).and.(force_contact/fmean > 1.D0)) write(30000,'(6(I6))') (5,j=1,6)
    if ((TAB_CONTACT_POLYR(i)%type == 1).and.(force_contact/fmean > 1.D0)) write(30000,'(6(I6))') (6,j=1,6)
    if ((TAB_CONTACT_POLYR(i)%type == 2).and.(force_contact/fmean > 1.D0)) write(30000,'(6(I6))') (7,j=1,6)
    if ((TAB_CONTACT_POLYR(i)%type >  2).and.(force_contact/fmean > 1.D0)) write(30000,'(6(I6))') (8,j=1,6)
!     if (sqrt(TAB_CONTACT_POLYR(i)%rt**2+TAB_CONTACT_POLYR(i)%rs**2)  >  0.39999 * TAB_CONTACT_POLYR(i)%rn ) then
!       write(30000,'(6(I6))') (1,j=1,6)
!     else
!       write(30000,'(6(I6))') (2,j=1,6)
!     end if

  end do
!  !-------> + Force week and strong
!  do i=1,nb_ligneCONTACT_POLYR  
!    if (TAB_CONTACT_POLYR(i)%nature=='PRPLx') cycle
!    force_contact = sqrt(TAB_CONTACT_POLYR(i)%rn**2 + TAB_CONTACT_POLYR(i)%rs**2 + TAB_CONTACT_POLYR(i)%rt**2 )  
!    if (force_contact==0) cycle 
!    if (force_contact/fmean < 1.D0 )  write(30000,'(6(I6))') (1,j=1,6)
!    if (force_contact/fmean >= 1.D0 ) write(30000,'(6(I6))') (2,j=1,6)
!  end do

  write(30000,'(A8)') 'polygons'
  write(30000,'(A7)') 'endpoly'

  write(30000,'(A6)') 'endgmv' 
  close(30000)
end if
!stop
end subroutine write_all_gmv

subroutine construct_list_material(nb_bodie_material,tab_bodie_mat_TEMP)
implicit none
integer                                     ::  i,j,nb_bodie_material
character(len=5),dimension(:),allocatable   ::  tab_bodie_mat
character(len=5),dimension(20)              ::  tab_bodie_mat_TEMP !ON SUPPOSE QU IL N Y A PAS PLUS DE 20 BEHAV...

tab_bodie_mat_TEMP(:) = 'VIDE '
nb_bodie_material     = 0

if ( allocated(tab_bodie_mat) ) deallocate(tab_bodie_mat)
allocate( tab_bodie_mat(nbParticules) )
do i=1,nbParticules
  tab_bodie_mat(i) = TAB_POLYEDRE(i)%behav
end do

do i=1,nbParticules
  do j=1,20
    if ( tab_bodie_mat(i) == tab_bodie_mat_TEMP(j) ) then
      exit
    else
      if ( tab_bodie_mat_TEMP(j) == 'VIDE ' ) then
        tab_bodie_mat_TEMP(j) = tab_bodie_mat(i)
        exit
      else
        cycle
      end if
    end if
  end do
end do

do i=1,20
  if (tab_bodie_mat_TEMP(i) /= 'VIDE ') nb_bodie_material = nb_bodie_material + 1
end do
end subroutine construct_list_material

!==================================================================================================
!Procedure pour ecrire le vtk a chaque pas de temps
!==================================================================================================
subroutine write_all_vtk
implicit none
character(len=17)                           ::  nom
integer                                     ::  i,j,k,nb_nodev,ncells
integer                                     ::  cd,an,total_vertex,total_cells
logical,dimension(:),allocatable            ::  tab_selected,tab_selected_force
real(kind=8)                                ::  win_x_max,win_x_min

if ( (step_all_vtk == compteur_clout).or. (step_all_vtk < 0) ) then

  print*,'....  WRITE VTK FILE....'

  compteur_vtk = compteur_clout  !compteur_vtk + 1    !

  nom       =  'VTK_FILE_0000.vtk'
  nb_nodev  = 0
  ncells    = 0
  if (compteur_vtk<10) then
    WRITE(nom(12:13),'(I1)')   compteur_vtk
  else if ( (compteur_vtk>=10) .and. (compteur_vtk<100) ) then
    WRITE(nom(11:13),'(I2)')   compteur_vtk
  else if ( (compteur_vtk>=100).and. (compteur_vtk<1000) ) then
    WRITE(nom(10:13),'(I3)')   compteur_vtk
  end if   

  win_x_max =-10000000
  win_x_min = 10000000
  do i=1,nbParticules
    do j=1,TAB_POLYEDRE(i)%nb_vertex
      win_x_max = max(win_x_max,TAB_POLYEDRE(i)%centre(1))
      win_x_min = min(win_x_min,TAB_POLYEDRE(i)%centre(1))
    end do
  end do

  allocate(tab_selected(nbParticules))
  tab_selected(:) = .false.
  do i=1,nbParticules
!    if ( (TAB_POLYEDRE(i)%centre(1) < win_x_max/4).and.&
!         (TAB_POLYEDRE(i)%centre(1) > win_x_min/4) ) tab_selected(i) = .true.
    tab_selected(i) = .true.
  end do
  
  open(unit=30000,file=nom,status='replace') 

  write(30000,'(A26)') '# vtk DataFile Version 3.8'
  write(30000,'(A19)') nom
  write(30000,'(A5)') 'ASCII'
  write(30000,'(A1)') ' '
  write(30000,'(A25)') 'DATASET UNSTRUCTURED_GRID'


!! Nombre total de vertex 
  total_vertex = 0
  do i=1,nbParticules
    if ( .not. tab_selected(i) ) cycle
    total_vertex = total_vertex + TAB_POLYEDRE(i)%nb_vertex
  end do 
  write(30000,'(A7,i15,A6)') 'POINTS ',total_vertex,' float'

!! Ecriture des vertex
  do i=1,nbParticules
!    print*, i, tab_selected(i),'vertex'
    if ( .not. tab_selected(i) ) cycle
    do j=1,TAB_POLYEDRE(i)%nb_vertex
      write(30000,'(3(F14.7,2X))') TAB_POLYEDRE(i)%vertex(1,j),&
                                   TAB_POLYEDRE(i)%vertex(2,j),&
                                   TAB_POLYEDRE(i)%vertex(3,j)
    end do
  end do

!! Nombre totals de cells
  total_cells = 0
  do i=1,nbParticules
    if ( .not. tab_selected(i) ) cycle
    total_cells = total_cells + TAB_POLYEDRE(i)%nb_face
  end do 
  write(30000,'(A6,i15,1X,i15)') 'CELLS ', total_cells , total_cells * 4    ! car (Face1 S1 S2 S3 )*nbface

  k = 0
  do i=1,nbParticules
!    print*, i, tab_selected(i),'face'
    if ( .not. tab_selected(i) ) cycle
    k = k + 1
    do j=1,TAB_POLYEDRE(i)%nb_face
      write(30000,'(i1,1X,3(i10,1x))') 3, TAB_POLYEDRE(i)%face(1,j) + (k-1)*TAB_POLYEDRE(i)%nb_vertex - 1,&
                                          TAB_POLYEDRE(i)%face(2,j) + (k-1)*TAB_POLYEDRE(i)%nb_vertex - 1,&
                                          TAB_POLYEDRE(i)%face(3,j) + (k-1)*TAB_POLYEDRE(i)%nb_vertex - 1
    end do
  end do

!! Cells types
  write(30000,'(A11,i15)') 'CELL_TYPES ', total_cells 
  do i=1,total_cells
!    print*,i,total_cells
    write(30000,'(i1)') 5
  end do
  
!! Cells data velocity
  write(30000,'(A10,i15)') 'CELL_DATA ', total_cells 
  write(30000,'(3(A7,1X))') 'SCALARS '   , 'Vitesse', 'float'
  write(30000,'(A20,1X)') 'LOOKUP_TABLE default'
  do i=1,nbParticules
    if ( .not. tab_selected(i) ) cycle
    do j=1,TAB_POLYEDRE(i)%nb_face
      write(30000,'(F14.7)') sqrt(TAB_POLYEDRE(i)%Vx**2+TAB_POLYEDRE(i)%Vy**2+TAB_POLYEDRE(i)%Vz**2)
    end do
  end do

!! Cells data nb_contacts
  write(30000,'(3(A7,1X))') 'SCALARS '   , 'nb_ctc', 'float'
  write(30000,'(A20,1X)') 'LOOKUP_TABLE default'
  do i=1,nb_ligneCONTACT_POLYR
    if (TAB_CONTACT_POLYR(i)%rn == 0.D0 ) cycle
    if (TAB_CONTACT_POLYR(i)%nature == 'PRPLx') cycle
    cd = TAB_CONTACT_POLYR(i)%icdent
    an = TAB_CONTACT_POLYR(i)%ianent
    TAB_POLYEDRE(cd)%nb_ctc = TAB_POLYEDRE(cd)%nb_ctc + 1
    TAB_POLYEDRE(an)%nb_ctc = TAB_POLYEDRE(an)%nb_ctc + 1
  end do
  do i=1,nbParticules
    if ( .not. tab_selected(i) ) cycle
    do j=1,TAB_POLYEDRE(i)%nb_face
      write(30000,'(F14.7)') real(TAB_POLYEDRE(i)%nb_ctc,8)
    end do
  end do

!! Cells data Pressure
  write(30000,'(3(A7,1X))') 'SCALARS '   , 'Press', 'float'
  write(30000,'(A20,1X)') 'LOOKUP_TABLE default'
  do i=1,nb_ligneCONTACT_POLYR
    if (TAB_CONTACT_POLYR(i)%rn == 0.D0 ) cycle
    if (TAB_CONTACT_POLYR(i)%nature == 'PRPLx') cycle
    cd = TAB_CONTACT_POLYR(i)%icdent
    an = TAB_CONTACT_POLYR(i)%ianent
    TAB_POLYEDRE(cd)%Pressure = TAB_POLYEDRE(cd)%Pressure + &
    sqrt(TAB_CONTACT_POLYR(i)%rn**2 + TAB_CONTACT_POLYR(i)%rs**2 + TAB_CONTACT_POLYR(i)%rt**2 )
    TAB_POLYEDRE(an)%Pressure = TAB_POLYEDRE(an)%Pressure + &
    sqrt(TAB_CONTACT_POLYR(i)%rn**2 + TAB_CONTACT_POLYR(i)%rs**2 + TAB_CONTACT_POLYR(i)%rt**2 )
  end do
  do i=1,nbParticules
    if ( .not. tab_selected(i) ) cycle
    do j=1,TAB_POLYEDRE(i)%nb_face
      write(30000,'(F14.7)') TAB_POLYEDRE(i)%Pressure
    end do
  end do


  close(30000)

end if
end subroutine write_all_vtk









!======================================================================
!Lecture des fichiers Vloc, Dof et BODIES et on commence  remplir notre type
!TAB_POLY.... c'est un peu bancale je l'accorde volontier :-)
!======================================================================
subroutine read_Vloc_dof_bodies
implicit none
integer                       ::  i=0,j=0,n_prpr,n_plpl,num_part,err,cd,an,nb_sommet,nb_face,k=0
integer                       ::  nb_prprx,nb_prplx,icdent,ianent,nb_ctc_deja_compte,cpt,cpt_actif
real(kind=8)                  ::  centrecluster1,centrecluster2,centrecluster3
real(kind=8)                  ::  rn,rs,rt,t1,t2,t3,n1,n2,n3,s1,s2,s3,vls,vln,vlt
real(kind=8)                  ::  I1,I2,I3,mean_Sphericity
character(len=6)              ::  text
character(len=5)              ::  statut
character(len=5)              ::  color,behav
character(len=13)             ::  text2
real(kind=8)                  ::  centre1,centre2,centre3,ax1,ax2,ax3,coor1,coor2,coor3
real(kind=8)                  ::  coor_ctc1,coor_ctc2,coor_ctc3
real(kind=8)                  ::  som1,som2,som3,Volume_tetra,Air_triangle,demi_somme
integer                       ::  face1,face2,face3
real(kind=8),dimension(3)     ::  X,Y,Z,centre_gravity,ab,ac,bc
real(kind=8),dimension(3,3)   ::  rot
real(kind=8)                  ::  Vrx,Vry,Vrz,vx,vy,vz,dcd_ptc_contact,dan_ptc_contact,n_ff,n_fs,n_fv,n_ss
character(len=22)             ::  clout_DOF
character(len=28)             ::  clout_Vloc
real(kind=8)                  ::  Number_mean_Total_force=0

if (all_pas==0) then 
  compteur_clout = compteur_clout+1
  clout_DOF  =  '../OUTBOX/DOF.OUT.    '
  clout_Vloc =  '../OUTBOX/Vloc_Rloc.OUT.    '
else
  compteur_clout = all_pas
  clout_DOF  =  '../OUTBOX/DOF.OUT.    '
  clout_Vloc =  '../OUTBOX/Vloc_Rloc.OUT.    '   
end if
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

if (the_first_time) then
  open(unit=2,file='../OUTBOX/BODIES.OUT',status='old')

  print*,'---> Lecture Bodies'
  
  mean_Sphericity = 0
  allocate(TAB_POLYEDRE(nbParticules))
  allocate(TAB_PLAN(nb_paroi))
  i=0
  do
    read(2,'(A6)') text
    if (text == '$bdyty') then
      i = i +1
      if (i<nbParticules+1) then
        if (i<numero_paroi) then
          TAB_POLYEDRE(i)%num=i
          read(2,*)
          read(2,*)
          read(2,'(22X, A5)') behav
          TAB_POLYEDRE(i)%behav = behav
          read(2,'(29X, 3(5x,D14.7,2X))') I1,I2,I3
          TAB_POLYEDRE(i)%I1 = I1
          TAB_POLYEDRE(i)%I2 = I2
          TAB_POLYEDRE(i)%I3 = I3
          read(2,*)
          read(2,'(29X, 3(5x,D14.7,2X))') centre1,centre2,centre3
          TAB_POLYEDRE(i)%centre_ref(1)=centre1
          TAB_POLYEDRE(i)%centre_ref(2)=centre2
          TAB_POLYEDRE(i)%centre_ref(3)=centre3
          read(2,*)
          read(2,*)
          read(2,'(22X,A5,12X,i7,13X,i7)') color,nb_sommet,nb_face
          TAB_POLYEDRE(i)%color     = color
          TAB_POLYEDRE(i)%nb_face   = nb_face
          TAB_POLYEDRE(i)%nb_vertex = nb_sommet
          allocate(TAB_POLYEDRE(i)%vertex_ref(3,nb_sommet))
          allocate(TAB_POLYEDRE(i)%vertex(3,nb_sommet))
          allocate(TAB_POLYEDRE(i)%face(3,nb_face))
          do j=1,TAB_POLYEDRE(i)%nb_vertex
            read(2,'(29X, 3(5x,D14.7,2X))') som1,som2,som3
            TAB_POLYEDRE(i)%vertex_ref(1,j) = som1
            TAB_POLYEDRE(i)%vertex_ref(2,j) = som2
            TAB_POLYEDRE(i)%vertex_ref(3,j) = som3
          end do
      
          centre_gravity = 0
          do j=1,TAB_POLYEDRE(i)%nb_vertex
           centre_gravity(1)=centre_gravity(1)+TAB_POLYEDRE(i)%vertex_ref(1,j)
           centre_gravity(2)=centre_gravity(2)+TAB_POLYEDRE(i)%vertex_ref(2,j)
           centre_gravity(3)=centre_gravity(3)+TAB_POLYEDRE(i)%vertex_ref(3,j)
          enddo
          centre_gravity=centre_gravity/real(TAB_POLYEDRE(i)%nb_vertex,8)

          TAB_POLYEDRE(i)%Rmax = 0
          do j=1,TAB_POLYEDRE(i)%nb_vertex
            X(1) = TAB_POLYEDRE(i)%vertex_ref(1,j) - centre_gravity(1)
            X(2) = TAB_POLYEDRE(i)%vertex_ref(2,j) - centre_gravity(2)
            X(3) = TAB_POLYEDRE(i)%vertex_ref(3,j) - centre_gravity(3)

            TAB_POLYEDRE(i)%Rmax = max( TAB_POLYEDRE(i)%Rmax , sqrt( X(1)**2 +  X(2)**2 +  X(3)**2 ) )
          end do
          Volume_tetra = 0.D0
          Air_triangle = 0.D0
          do j=1,TAB_POLYEDRE(i)%nb_face
            read(2,'(29X, 3(5x,i7,9X))') face1,face2,face3
            TAB_POLYEDRE(i)%face(1,j) = face1      
            TAB_POLYEDRE(i)%face(2,j) = face2      
            TAB_POLYEDRE(i)%face(3,j) = face3      
            X(1) = TAB_POLYEDRE(i)%vertex_ref(1,face1) - centre_gravity(1)
            X(2) = TAB_POLYEDRE(i)%vertex_ref(2,face1) - centre_gravity(2)
            X(3) = TAB_POLYEDRE(i)%vertex_ref(3,face1) - centre_gravity(3)
        
            Y(1) = TAB_POLYEDRE(i)%vertex_ref(1,face2) - centre_gravity(1)
            Y(2) = TAB_POLYEDRE(i)%vertex_ref(2,face2) - centre_gravity(2)
            Y(3) = TAB_POLYEDRE(i)%vertex_ref(3,face2) - centre_gravity(3)
        
            Z(1) = TAB_POLYEDRE(i)%vertex_ref(1,face3) - centre_gravity(1)
            Z(2) = TAB_POLYEDRE(i)%vertex_ref(2,face3) - centre_gravity(2)
            Z(3) = TAB_POLYEDRE(i)%vertex_ref(3,face3) - centre_gravity(3)

            Volume_tetra = Volume_tetra + abs( X(1)*(Y(2)*Z(3)-Y(3)*Z(2)) &
                                              -X(2)*(Y(1)*Z(3)-Y(3)*Z(1)) &
                                              +X(3)*(Y(1)*Z(2)-Y(2)*Z(1))  )
          
            ab(:) = X(:)-Y(:) 
            ac(:) = X(:)-Z(:) 
            bc(:) = Y(:)-Z(:) 
            demi_somme = 0.5*( sqrt( ab(1)**2+ab(2)**2+ab(3)**2 ) +&
                               sqrt( ac(1)**2+ac(2)**2+ac(3)**2 ) +&
                               sqrt( bc(1)**2+bc(2)**2+bc(3)**2 ))
        
            Air_triangle = Air_triangle + sqrt( demi_somme * (demi_somme - sqrt( ab(1)**2+ab(2)**2+ab(3)**2 ) )*&
                                                             (demi_somme - sqrt( ac(1)**2+ac(2)**2+ac(3)**2 ) )*&
                                                             (demi_somme - sqrt( bc(1)**2+bc(2)**2+bc(3)**2 ) ) )
          end do
          TAB_POLYEDRE(i)%volume     = Volume_tetra/real(6,8)
          TAB_POLYEDRE(i)%Aire       = Air_triangle
          TAB_POLYEDRE(i)%Sphericity = (PI**(0.333333333)) * ((6*TAB_POLYEDRE(i)%volume)**(0.66666666)) / &
                                        TAB_POLYEDRE(i)%Aire
        
          TAB_POLYEDRE(i)%centre(1)=0.D0
          TAB_POLYEDRE(i)%momentI=0
          TAB_POLYEDRE(i)%nctc = 0
          TAB_POLYEDRE(i)%nctcslide = 0
          TAB_POLYEDRE(i)%nctcstick = 0
          TAB_POLYEDRE(i)%nctcs = 0
          TAB_POLYEDRE(i)%nctcd = 0
          TAB_POLYEDRE(i)%nctct = 0
          TAB_POLYEDRE(i)%nb_ctc = 0
          TAB_POLYEDRE(i)%force = 0.D0
          TAB_POLYEDRE(i)%Pressure = 0.D0
         
          mean_Sphericity = mean_Sphericity + TAB_POLYEDRE(i)%Sphericity
        else if (i==numero_paroi) then
          read(2,*)
          read(2,*)
          read(2,'(22X,A5)') behav
          read(2,*)
          read(2,*)
          read(2,'(29X, 3(5x,D14.7,2X))') centrecluster1,centrecluster2,centrecluster3
          read(2,*)
          read(2,*)
          do j=numero_paroi,nbParticules
            read(2,'(22X,A5,12X,i7,13X,i7)') color,nb_sommet,nb_face
            TAB_POLYEDRE(j)%num     = j
            TAB_POLYEDRE(j)%color     = color
            TAB_POLYEDRE(j)%nb_face   = nb_face
            TAB_POLYEDRE(j)%nb_vertex = nb_sommet
            TAB_POLYEDRE(j)%behav     = behav
            TAB_POLYEDRE(j)%centre_ref = 0.D0
            TAB_POLYEDRE(j)%centre     = 0.D0
            allocate(TAB_POLYEDRE(j)%vertex_ref(3,nb_sommet))
            allocate(TAB_POLYEDRE(j)%vertex(3,nb_sommet))
            allocate(TAB_POLYEDRE(j)%face(3,nb_face))
            do k=1,TAB_POLYEDRE(j)%nb_vertex
              read(2,'(29X, 3(5x,D14.7,2X))') som1,som2,som3
              TAB_POLYEDRE(j)%vertex_ref(1,k) = som1+centrecluster1
              TAB_POLYEDRE(j)%vertex_ref(2,k) = som2+centrecluster2
              TAB_POLYEDRE(j)%vertex_ref(3,k) = som3+centrecluster3
              TAB_POLYEDRE(j)%centre_ref(1:3)=TAB_POLYEDRE(j)%centre_ref(1:3)+TAB_POLYEDRE(j)%vertex_ref(1:3,k)

            end do
            TAB_POLYEDRE(j)%centre_ref(1:3)=TAB_POLYEDRE(j)%centre_ref(1:3)/real(TAB_POLYEDRE(j)%nb_vertex,8)

            do k=1,TAB_POLYEDRE(j)%nb_vertex
              TAB_POLYEDRE(j)%vertex_ref(1,k) = TAB_POLYEDRE(j)%vertex_ref(1,k)-TAB_POLYEDRE(j)%centre_ref(1)
              TAB_POLYEDRE(j)%vertex_ref(2,k) = TAB_POLYEDRE(j)%vertex_ref(2,k)-TAB_POLYEDRE(j)%centre_ref(2)
              TAB_POLYEDRE(j)%vertex_ref(3,k) = TAB_POLYEDRE(j)%vertex_ref(3,k)-TAB_POLYEDRE(j)%centre_ref(3)
            end do
            
            centre_gravity = TAB_POLYEDRE(j)%centre_ref
            TAB_POLYEDRE(j)%Rmax = 0
            do k=1,TAB_POLYEDRE(j)%nb_vertex
              X(1) = TAB_POLYEDRE(j)%vertex_ref(1,k) - centre_gravity(1)
              X(2) = TAB_POLYEDRE(j)%vertex_ref(2,k) - centre_gravity(2)
              X(3) = TAB_POLYEDRE(j)%vertex_ref(3,k) - centre_gravity(3)

              TAB_POLYEDRE(j)%Rmax = max( TAB_POLYEDRE(j)%Rmax , sqrt( X(1)**2 +  X(2)**2 +  X(3)**2 ) )
            end do
            Volume_tetra = 0.D0
            Air_triangle = 0.D0
            do k=1,TAB_POLYEDRE(j)%nb_face
              read(2,'(29X, 3(5x,i7,9X))') face1,face2,face3
              TAB_POLYEDRE(j)%face(1,k) = face1      
              TAB_POLYEDRE(j)%face(2,k) = face2      
              TAB_POLYEDRE(j)%face(3,k) = face3      
              X(1) = TAB_POLYEDRE(j)%vertex_ref(1,face1) - centre_gravity(1)
              X(2) = TAB_POLYEDRE(j)%vertex_ref(2,face1) - centre_gravity(2)
              X(3) = TAB_POLYEDRE(j)%vertex_ref(3,face1) - centre_gravity(3)
        
              Y(1) = TAB_POLYEDRE(j)%vertex_ref(1,face2) - centre_gravity(1)
              Y(2) = TAB_POLYEDRE(j)%vertex_ref(2,face2) - centre_gravity(2)
              Y(3) = TAB_POLYEDRE(j)%vertex_ref(3,face2) - centre_gravity(3)
        
              Z(1) = TAB_POLYEDRE(j)%vertex_ref(1,face3) - centre_gravity(1)
              Z(2) = TAB_POLYEDRE(j)%vertex_ref(2,face3) - centre_gravity(2)
              Z(3) = TAB_POLYEDRE(j)%vertex_ref(3,face3) - centre_gravity(3)

              Volume_tetra = Volume_tetra + abs( X(1)*(Y(2)*Z(3)-Y(3)*Z(2)) &
                                                -X(2)*(Y(1)*Z(3)-Y(3)*Z(1)) &
                                                +X(3)*(Y(1)*Z(2)-Y(2)*Z(1))  )
          
              ab(:) = X(:)-Y(:) 
              ac(:) = X(:)-Z(:) 
              bc(:) = Y(:)-Z(:) 
              demi_somme = 0.5*( sqrt( ab(1)**2+ab(2)**2+ab(3)**2 ) +&
                                 sqrt( ac(1)**2+ac(2)**2+ac(3)**2 ) +&
                                 sqrt( bc(1)**2+bc(2)**2+bc(3)**2 ))
        
              Air_triangle = Air_triangle + sqrt( demi_somme * (demi_somme - sqrt( ab(1)**2+ab(2)**2+ab(3)**2 ) )*&
                                                               (demi_somme - sqrt( ac(1)**2+ac(2)**2+ac(3)**2 ) )*&
                                                               (demi_somme - sqrt( bc(1)**2+bc(2)**2+bc(3)**2 ) ) )
            end do
            TAB_POLYEDRE(j)%volume     = Volume_tetra/real(6,8)
            TAB_POLYEDRE(j)%Aire       = Air_triangle
            TAB_POLYEDRE(j)%Sphericity = (PI**(0.333333333)) * ((6*TAB_POLYEDRE(j)%volume)**(0.66666666)) / &
                                          TAB_POLYEDRE(j)%Aire
        
            TAB_POLYEDRE(j)%centre(1)=0.D0
            TAB_POLYEDRE(j)%momentI=0
            TAB_POLYEDRE(j)%nctc = 0
            TAB_POLYEDRE(j)%nctcslide = 0
            TAB_POLYEDRE(j)%nctcstick = 0
            TAB_POLYEDRE(j)%nctcs = 0
            TAB_POLYEDRE(j)%nctcd = 0
            TAB_POLYEDRE(j)%nctct = 0
            TAB_POLYEDRE(j)%nb_ctc = 0
            TAB_POLYEDRE(j)%force = 0.D0
            TAB_POLYEDRE(j)%Pressure = 0.D0
          
            mean_Sphericity = mean_Sphericity + TAB_POLYEDRE(j)%Sphericity
          end do
          exit
        end if
      end if 
    end if

  end do
  the_first_time=.false.
  print*,'Sphericity moyenne = ', mean_Sphericity / real(nbParticules,8)
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
        if (num_part<numero_paroi) then
          read(5,*)
          read(5,'(29X, 3(5x,D14.7,2X))') coor1,coor2,coor3
          TAB_POLYEDRE(i)%centre(1)=coor1+TAB_POLYEDRE(i)%centre_ref(1)
          TAB_POLYEDRE(i)%centre(2)=coor2+TAB_POLYEDRE(i)%centre_ref(2)
          TAB_POLYEDRE(i)%centre(3)=coor3+TAB_POLYEDRE(i)%centre_ref(3)
          read(5,*)
          read(5,'(29X, 3(5x,D14.7,2X))') Vx,Vy,Vz
          read(5,'(29X, 3(5x,D14.7,2X))') Vrx,Vry,Vrz
          TAB_POLYEDRE(i)%Vx=Vx
          TAB_POLYEDRE(i)%Vy=Vy
          TAB_POLYEDRE(i)%Vz=Vz
          TAB_POLYEDRE(i)%Vrx=Vrx
          TAB_POLYEDRE(i)%Vry=Vry
          TAB_POLYEDRE(i)%Vrz=Vrz
          TAB_POLYEDRE(i)%Norme_V  = sqrt(Vx**2+Vy**2+Vz**2)
          TAB_POLYEDRE(i)%Norme_Vr = sqrt(Vrx**2+Vry**2+Vrz**2)
          read(5,'(29X, 3(5x,D14.7,2X))') rot(1,1),rot(2,1),rot(3,1)
          read(5,'(29X, 3(5x,D14.7,2X))') rot(1,2),rot(2,2),rot(3,2)
          read(5,'(29X, 3(5x,D14.7,2X))') rot(1,3),rot(2,3),rot(3,3)
          TAB_POLYEDRE(i)%Rot = rot
          TAB_POLYEDRE(i)%Wx = Vrx*rot(1,1) + Vry*rot(1,2) + Vrz*rot(1,3)
          TAB_POLYEDRE(i)%Wy = Vrx*rot(2,1) + Vry*rot(2,2) + Vrz*rot(2,3)
          TAB_POLYEDRE(i)%Wz = Vrx*rot(3,1) + Vry*rot(3,2) + Vrz*rot(3,3)

          do j=1,TAB_POLYEDRE(i)%nb_vertex                                                 ! On rotationne et translate   
            TAB_POLYEDRE(i)%vertex(1,j)=TAB_POLYEDRE(i)%vertex_ref(1,j)*rot(1,1)+&         ! histoire de tout avoir dans
                                        TAB_POLYEDRE(i)%vertex_ref(2,j)*rot(1,2)+&         ! le rep�re global...
                                        TAB_POLYEDRE(i)%vertex_ref(3,j)*rot(1,3)+&
                                        TAB_POLYEDRE(i)%centre(1)

            TAB_POLYEDRE(i)%vertex(2,j)=TAB_POLYEDRE(i)%vertex_ref(1,j)*rot(2,1)+&
                                        TAB_POLYEDRE(i)%vertex_ref(2,j)*rot(2,2)+&
                                        TAB_POLYEDRE(i)%vertex_ref(3,j)*rot(2,3)+&
                                        TAB_POLYEDRE(i)%centre(2)

            TAB_POLYEDRE(i)%vertex(3,j)=TAB_POLYEDRE(i)%vertex_ref(1,j)*rot(3,1)+&
                                        TAB_POLYEDRE(i)%vertex_ref(2,j)*rot(3,2)+&
                                        TAB_POLYEDRE(i)%vertex_ref(3,j)*rot(3,3)+&
                                        TAB_POLYEDRE(i)%centre(3)
          end do
          TAB_POLYEDRE(i)%nctc = 0
          TAB_POLYEDRE(i)%nctcslide = 0
          TAB_POLYEDRE(i)%nctcstick = 0
          TAB_POLYEDRE(i)%nctcs = 0
          TAB_POLYEDRE(i)%nctcd = 0
          TAB_POLYEDRE(i)%nctct = 0
          TAB_POLYEDRE(i)%nb_ctc = 0
          TAB_POLYEDRE(i)%force = 0.D0
        else if (num_part==numero_paroi) then
          read(5,*)
          read(5,'(29X, 3(5x,D14.7,2X))') coor1,coor2,coor3
          read(5,*)
          read(5,'(29X, 3(5x,D14.7,2X))') Vx,Vy,Vz
          read(5,'(29X, 3(5x,D14.7,2X))') Vrx,Vry,Vrz
          read(5,'(29X, 3(5x,D14.7,2X))') rot(1,1),rot(2,1),rot(3,1)
          read(5,'(29X, 3(5x,D14.7,2X))') rot(1,2),rot(2,2),rot(3,2)
          read(5,'(29X, 3(5x,D14.7,2X))') rot(1,3),rot(2,3),rot(3,3)
          do j=num_part,nbParticules
            TAB_POLYEDRE(j)%Vx=Vx
            TAB_POLYEDRE(j)%Vy=Vy
            TAB_POLYEDRE(j)%Vz=Vz
            TAB_POLYEDRE(j)%Vrx=Vrx
            TAB_POLYEDRE(j)%Vry=Vry
            TAB_POLYEDRE(j)%Vrz=Vrz
            TAB_POLYEDRE(j)%Norme_V  = sqrt(Vx**2+Vy**2+Vz**2)
            TAB_POLYEDRE(j)%Norme_Vr = sqrt(Vrx**2+Vry**2+Vrz**2)
            TAB_POLYEDRE(j)%Rot = rot
            TAB_POLYEDRE(j)%Wx = Vrx*rot(1,1) + Vry*rot(1,2) + Vrz*rot(1,3)
            TAB_POLYEDRE(j)%Wy = Vrx*rot(2,1) + Vry*rot(2,2) + Vrz*rot(2,3)
            TAB_POLYEDRE(j)%Wz = Vrx*rot(3,1) + Vry*rot(3,2) + Vrz*rot(3,3)
            TAB_POLYEDRE(j)%centre(1)=coor1+TAB_POLYEDRE(j)%centre_ref(1)
            TAB_POLYEDRE(j)%centre(2)=coor2+TAB_POLYEDRE(j)%centre_ref(2)
            TAB_POLYEDRE(j)%centre(3)=coor3+TAB_POLYEDRE(j)%centre_ref(3)

            do k=1,TAB_POLYEDRE(j)%nb_vertex                                                 ! On rotationne et translate   
              TAB_POLYEDRE(j)%vertex(1,k)=TAB_POLYEDRE(j)%vertex_ref(1,k)*rot(1,1)+&         ! histoire de tout avoir dans
                                          TAB_POLYEDRE(j)%vertex_ref(2,k)*rot(1,2)+&         ! le rep�re global...
                                          TAB_POLYEDRE(j)%vertex_ref(3,k)*rot(1,3)+&
                                          TAB_POLYEDRE(j)%centre(1)

              TAB_POLYEDRE(j)%vertex(2,k)=TAB_POLYEDRE(j)%vertex_ref(1,k)*rot(2,1)+&
                                          TAB_POLYEDRE(j)%vertex_ref(2,k)*rot(2,2)+&
                                          TAB_POLYEDRE(j)%vertex_ref(3,k)*rot(2,3)+&
                                          TAB_POLYEDRE(j)%centre(2)

              TAB_POLYEDRE(j)%vertex(3,k)=TAB_POLYEDRE(j)%vertex_ref(1,k)*rot(3,1)+&
                                          TAB_POLYEDRE(j)%vertex_ref(2,k)*rot(3,2)+&
                                          TAB_POLYEDRE(j)%vertex_ref(3,k)*rot(3,3)+&
                                          TAB_POLYEDRE(j)%centre(3)
            end do

            TAB_POLYEDRE(j)%nctc = 0
            TAB_POLYEDRE(j)%nctcslide = 0
            TAB_POLYEDRE(j)%nctcstick = 0
            TAB_POLYEDRE(j)%nctcs = 0
            TAB_POLYEDRE(j)%nctcd = 0
            TAB_POLYEDRE(j)%nctct = 0
            TAB_POLYEDRE(j)%nb_ctc = 0
            TAB_POLYEDRE(j)%force = 0.D0
          end do
        end if
      end if
    end if
  end do
end if
nb_PRPRx=0
nb_PRPLx=0
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
    if (text2 == '$icdan  PRPLx') then
      nb_PRPLx = nb_PRPLx +1
    end if
    if (text2 == '$icdan  PRPRx') then
      nb_PRPRx = nb_PRPRx +1
    end if
  end do
  close(4)
  if (allocated(TAB_CONTACT)) deallocate(TAB_CONTACT)
  allocate(TAB_CONTACT(nb_PRPRx+nb_PRPLx+1))
  do i=1,nb_PRPRx+nb_PRPLx
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
  n_prpr = 0
  n_plpl = 0
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
    if (text2 == '$icdan  PRPLx') then
      read(4,*)
      read(4,'(7X,i6,36X,i6,16X,A5)') icdent,ianent,statut
      read(4,'(29X, 3(5x,D14.7,2X))') rs,rt,rn
      read(4,'(29X, 3(5x,D14.7,2X))') vls,vlt,vln
      read(4,*)
      read(4,'(29X, 3(5x,D14.7,2X))') coor_ctc1,coor_ctc2,coor_ctc3
      read(4,'(29X, 3(5x,D14.7,2X))') t1,t2,t3
      read(4,'(29X, 3(5x,D14.7,2X))') n1,n2,n3
      read(4,'(29X, 3(5x,D14.7,2X))') s1,s2,s3
      !---------------------------------------
      nb_ligneCONTACT=nb_ligneCONTACT+1
      n_plpl = n_plpl + 1
      TAB_CONTACT(n_plpl+nb_PRPRx)%icdent=icdent
      TAB_CONTACT(n_plpl+nb_PRPRx)%ianent=ianent
      TAB_CONTACT(n_plpl+nb_PRPRx)%statut=statut
      TAB_CONTACT(n_plpl+nb_PRPRx)%n(1)=n1
      TAB_CONTACT(n_plpl+nb_PRPRx)%n(2)=n2
      TAB_CONTACT(n_plpl+nb_PRPRx)%n(3)=n3
      TAB_CONTACT(n_plpl+nb_PRPRx)%t(1)=t1
      TAB_CONTACT(n_plpl+nb_PRPRx)%t(2)=t2
      TAB_CONTACT(n_plpl+nb_PRPRx)%t(3)=t3
      TAB_CONTACT(n_plpl+nb_PRPRx)%s(1)=s1
      TAB_CONTACT(n_plpl+nb_PRPRx)%s(2)=s2
      TAB_CONTACT(n_plpl+nb_PRPRx)%s(3)=s3
      TAB_CONTACT(n_plpl+nb_PRPRx)%coor_ctc(1)=coor_ctc1
      TAB_CONTACT(n_plpl+nb_PRPRx)%coor_ctc(2)=coor_ctc2
      TAB_CONTACT(n_plpl+nb_PRPRx)%coor_ctc(3)=coor_ctc3
      TAB_CONTACT(n_plpl+nb_PRPRx)%rn=rn
      TAB_CONTACT(n_plpl+nb_PRPRx)%rt=rt
      TAB_CONTACT(n_plpl+nb_PRPRx)%rs=rs
      TAB_CONTACT(n_plpl+nb_PRPRx)%vln=vln
      TAB_CONTACT(n_plpl+nb_PRPRx)%vlt=vlt
      TAB_CONTACT(n_plpl+nb_PRPRx)%vls=vls
      TAB_CONTACT(n_plpl+nb_PRPRx)%type=1
      TAB_CONTACT(n_plpl+nb_PRPRx)%nature = 'PRPLx'
      TAB_CONTACT(n_plpl+nb_PRPRx)%deja_compte = .false.
    end if
    if (text2 == '$icdan  PRPRx') then
      read(4,*)
      read(4,'(7X,i6,36X,i6,16X,A5)') icdent,ianent,statut
      read(4,'(29X, 3(5x,D14.7,2X))') rt,rn,rs
      read(4,'(29X, 3(5x,D14.7,2X))') vlt,vln,vls
      read(4,*)
      read(4,'(29X, 3(5x,D14.7,2X))') coor_ctc1,coor_ctc2,coor_ctc3
      read(4,'(29X, 3(5x,D14.7,2X))') t1,t2,t3
      read(4,'(29X, 3(5x,D14.7,2X))') n1,n2,n3
      read(4,'(29X, 3(5x,D14.7,2X))') s1,s2,s3
      !--------------------------
      nb_ligneCONTACT=nb_ligneCONTACT+1
      n_prpr = n_prpr + 1 
      TAB_CONTACT(n_prpr)%icdent=icdent
      TAB_CONTACT(n_prpr)%ianent=ianent
      TAB_CONTACT(n_prpr)%statut=statut
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
    end if
  end do
!----- Ici il faut definir le type de chaque contact ----------
  print*,'     creation de la nouvelle liste de contact'
  nb_ctc_deja_compte = 0
  do i=1,nb_ligneCONTACT
    if ( (TAB_CONTACT(i)%deja_compte) ) cycle
    
    cd = TAB_CONTACT(i)%icdent
    an = TAB_CONTACT(i)%ianent
    do j=i+1,nb_ligneCONTACT
      if (TAB_CONTACT(i)%type==4) exit ! car le max se sont des contacts quadruples
      if (j-i > 5) exit                ! car en fait les contacts face/face sont ecrit les un � la suite des autres
      if ((cd == TAB_CONTACT(j)%icdent).and.(an==TAB_CONTACT(j)%ianent).or.&
          (cd == TAB_CONTACT(j)%ianent).and.(an==TAB_CONTACT(j)%icdent) ) then

        TAB_CONTACT(j)%deja_compte = .true.
        nb_ctc_deja_compte = nb_ctc_deja_compte + 1

        if (  (TAB_CONTACT(i)%rn>0.D0).and.(TAB_CONTACT(j)%rn>0.D0) ) then
          TAB_CONTACT(i)%rn = TAB_CONTACT(i)%rn + TAB_CONTACT(j)%rn
          TAB_CONTACT(i)%rt = TAB_CONTACT(i)%rt + TAB_CONTACT(j)%rt
          TAB_CONTACT(i)%rs = TAB_CONTACT(i)%rs + TAB_CONTACT(j)%rs

          TAB_CONTACT(i)%vln = TAB_CONTACT(i)%vln + TAB_CONTACT(j)%vln
          TAB_CONTACT(i)%vlt = TAB_CONTACT(i)%vlt + TAB_CONTACT(j)%vlt
          TAB_CONTACT(i)%vls = TAB_CONTACT(i)%vls + TAB_CONTACT(j)%vls
        
          TAB_CONTACT(i)%coor_ctc(:) = TAB_CONTACT(i)%coor_ctc(:) + TAB_CONTACT(j)%coor_ctc(:)
        
          TAB_CONTACT(i)%type = TAB_CONTACT(i)%type + 1
          
        else if (  (TAB_CONTACT(i)%rn==0.D0).and.(TAB_CONTACT(j)%rn>0.D0) ) then
          TAB_CONTACT(i)%rn = TAB_CONTACT(j)%rn
          TAB_CONTACT(i)%rt = TAB_CONTACT(j)%rt
          TAB_CONTACT(i)%rs = TAB_CONTACT(j)%rs

          TAB_CONTACT(i)%vln = TAB_CONTACT(j)%vln
          TAB_CONTACT(i)%vlt = TAB_CONTACT(j)%vlt
          TAB_CONTACT(i)%vls = TAB_CONTACT(j)%vls
        
          TAB_CONTACT(i)%coor_ctc(:) = TAB_CONTACT(j)%coor_ctc(:)
          
          TAB_CONTACT(i)%type = TAB_CONTACT(j)%type
        end if
        
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
  do i=1,nb_ligneCONTACT
    if ( (TAB_CONTACT(i)%deja_compte) ) cycle
    cpt = cpt + 1
    dcd_ptc_contact = 1000000000.D0
    dan_ptc_contact = 1000000000.D0

    TAB_CONTACT_POLYR(cpt)%icdent   =   TAB_CONTACT(i)%icdent
    TAB_CONTACT_POLYR(cpt)%ianent   =   TAB_CONTACT(i)%ianent
    TAB_CONTACT_POLYR(cpt)%n(:)     =   TAB_CONTACT(i)%n(:)
    TAB_CONTACT_POLYR(cpt)%t(:)     =   TAB_CONTACT(i)%t(:)
    TAB_CONTACT_POLYR(cpt)%s(:)     =   TAB_CONTACT(i)%s(:)
    TAB_CONTACT_POLYR(cpt)%coor_ctc =   TAB_CONTACT(i)%coor_ctc / TAB_CONTACT(i)%type
    TAB_CONTACT_POLYR(cpt)%rn       =   TAB_CONTACT(i)%rn
    TAB_CONTACT_POLYR(cpt)%rt       =   TAB_CONTACT(i)%rt
    TAB_CONTACT_POLYR(cpt)%rs       =   TAB_CONTACT(i)%rs
    TAB_CONTACT_POLYR(cpt)%vln      =   TAB_CONTACT(i)%vln / TAB_CONTACT(i)%type
    TAB_CONTACT_POLYR(cpt)%vlt      =   TAB_CONTACT(i)%vlt / TAB_CONTACT(i)%type
    TAB_CONTACT_POLYR(cpt)%vls      =   TAB_CONTACT(i)%vls / TAB_CONTACT(i)%type
    TAB_CONTACT_POLYR(cpt)%type     =   TAB_CONTACT(i)%type
    TAB_CONTACT_POLYR(cpt)%nature   =   TAB_CONTACT(i)%nature

    
    if (TAB_CONTACT_POLYR(cpt)%type == 2) n_fs = n_fs + 1 
    if (TAB_CONTACT_POLYR(cpt)%type == 3) n_ff = n_ff + 1 
    if (TAB_CONTACT_POLYR(cpt)%type == 4) n_ff = n_ff + 1 
    
    if (TAB_CONTACT_POLYR(cpt)%type==1) then     ! Il faut distinguer les contacts face-vertex des contacts side-side
      cd = TAB_CONTACT_POLYR(cpt)%icdent
      an = TAB_CONTACT_POLYR(cpt)%ianent
      
      if (cd>nbParticules.or.an>nbParticules) then
        n_fv = n_fv + 1
        TAB_CONTACT_POLYR(cpt)%type = 1

      else
      
        do j=1,TAB_POLYEDRE(cd)%nb_vertex
          dcd_ptc_contact = min( dcd_ptc_contact , sqrt( (TAB_CONTACT_POLYR(cpt)%coor_ctc(1)-TAB_POLYEDRE(cd)%vertex(1,j))**2+&
                                                         (TAB_CONTACT_POLYR(cpt)%coor_ctc(2)-TAB_POLYEDRE(cd)%vertex(2,j))**2+&
                                                         (TAB_CONTACT_POLYR(cpt)%coor_ctc(3)-TAB_POLYEDRE(cd)%vertex(3,j))**2 ) )
        end do 
        do j=1,TAB_POLYEDRE(an)%nb_vertex
          dan_ptc_contact = min( dan_ptc_contact , sqrt( (TAB_CONTACT_POLYR(cpt)%coor_ctc(1)-TAB_POLYEDRE(an)%vertex(1,j))**2+&
                                                         (TAB_CONTACT_POLYR(cpt)%coor_ctc(2)-TAB_POLYEDRE(an)%vertex(2,j))**2+&
                                                         (TAB_CONTACT_POLYR(cpt)%coor_ctc(3)-TAB_POLYEDRE(an)%vertex(3,j))**2 ) )
        end do 
            
        if ( min( dcd_ptc_contact , dan_ptc_contact ) < 0.0001 ) then
          n_fv = n_fv + 1
          TAB_CONTACT_POLYR(cpt)%type = 1
        else
          n_ss = n_ss + 1
          TAB_CONTACT_POLYR(cpt)%type = 0
        end if
      end if
      
    end if
    
    
    if ( TAB_CONTACT_POLYR(cpt)%vln**2 + TAB_CONTACT_POLYR(cpt)%vlt**2 > 0.000000001 )   &
         TAB_CONTACT_POLYR(cpt)%statut='slide'
    if ( TAB_CONTACT_POLYR(cpt)%vln**2 + TAB_CONTACT_POLYR(cpt)%vlt**2 < 0.000000001 )   &
         TAB_CONTACT_POLYR(cpt)%statut='stick'
    if (TAB_CONTACT_POLYR(cpt)%rn>0.D0)  cpt_actif = cpt_actif +1
  end do

  print*,'Nb particules =',nbParticules,'Nb Contact',nb_ligneCONTACT
  print*,'Nb reel Contact',nb_ligneCONTACT_POLYR
  print*,'Nb Contacts actif',cpt_actif
  print*,' n_ff ',n_ff, ' n_fs ',n_fs,' n_ss ',n_ss,' n_fv ',n_fv

  Mean_total_Normal_force = 0.D0
  Number_mean_Total_force = 0.D0
  do i=1,nb_ligneCONTACT_POLYR
    if (TAB_CONTACT_POLYR(i)%rn == 0.D0 ) cycle
    Mean_total_Normal_force = Mean_total_Normal_force + TAB_CONTACT(i)%rn
    Number_mean_Total_force = Number_mean_Total_force + 1.D0
  end do
  Mean_total_Normal_force = Mean_total_Normal_force / Number_mean_Total_force

end if

end subroutine read_Vloc_dof_bodies




!==============================================================================
! Reconstruction d'un BODIES.DAT, DOF.INI, DRV_DOF.DAT 
!==============================================================================



subroutine construct_bodies
implicit none
real(kind=8)                             :: dmean,H_max,H_ini,H_temp,H_temp2,ep,Volume_grain,Volume_wall
integer                                  :: i,cpt,nbPart,nb_wallh,nb_wallb,k,j

H_ini = 0.48
ep    = 0.D0

if  (compteur_clout == quel_pas_de_temps) then
  dmean = 0.D0
  H_max = 0.D0
  do i=1,nbParticules
    dmean = dmean + 2*TAB_POLYEDRE(i)%Rmax
    H_max = max(H_max,TAB_POLYEDRE(i)%centre(3))
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
    TAB_POLYEDRE(i)%behav = 'plexx'
    TAB_POLYEDRE(i)%color = 'TATAx'
    if ( TAB_POLYEDRE(i)%centre(3) < 2*dmean ) TAB_POLYEDRE(i)%behav='wallb'
    if ( TAB_POLYEDRE(i)%centre(3) > H_max - 2.5*dmean ) TAB_POLYEDRE(i)%behav='wallh'
  end do

  cpt    = 0.D0
  H_temp = 0.D0
  Volume_grain = 0
  H_temp2 = 10000000.D0
  
  do i=1,nbParticules
    if (TAB_POLYEDRE(i)%behav /= 'plexx')  cycle
    cpt = cpt + 1
    
    write(128,'(A6)')                                        '$bdyty'
    write(128,'(1X,A5,1X,i6)')                               'RBDY3',cpt
    write(128,'(A6)')                                        '$blmty'
    write(128,'(1X,A5,5X,i2,2X,A5,2X,A5,1X,A6,D14.7)')       'PLAIN',1,'behav','plexx','avrd=',0.D0
    write(128,'(29X,3(A5,D14.7,2X))')  'I1  =',TAB_POLYEDRE(i)%I1,'I2  =',TAB_POLYEDRE(i)%I2,'I3  =',TAB_POLYEDRE(i)%I3
    write(128,'(A6)')                                        '$nodty'
    write(128,'(1X,A5,5X,i2,16X,3(A5,D14.7,2X))')            'NO6xx',1,'coo1=',TAB_POLYEDRE(i)%centre(1),&
                                                                       'coo2=',TAB_POLYEDRE(i)%centre(2),&
                                                                       'coo3=',TAB_POLYEDRE(i)%centre(3)
    write(128,'(29X,3(A5,D14.7,2X))') 'coo4=',0.D0,'coo5=',0.D0,'coo6=',0.D0
    write(128,'(A6)')                                        '$tacty'
    write(128,'(1X,A5,5X,i2,2X,A5,2X,A5,2X,A10,i7,4X,A9,i7)') 'POLYR',1,'color','TATAx',&
                                                              'nb_vertex=',TAB_POLYEDRE(i)%nb_vertex,&
                                                              'nb_faces=',TAB_POLYEDRE(i)%nb_face
    do k=1,TAB_POLYEDRE(i)%nb_vertex
      write(128,'(29X,3(A5,D14.7,2X))') 'coo1=',TAB_POLYEDRE(i)%vertex(1,k)-TAB_POLYEDRE(i)%centre(1),&
                                        'coo2=',TAB_POLYEDRE(i)%vertex(2,k)-TAB_POLYEDRE(i)%centre(2),&
                                        'coo3=',TAB_POLYEDRE(i)%vertex(3,k)-TAB_POLYEDRE(i)%centre(3)
    end do
    do k=1,TAB_POLYEDRE(i)%nb_face
       write(128,'(29X,3(A5,i7,9X))')   'ver1=',TAB_POLYEDRE(i)%face(1,k),&
                                        'ver2=',TAB_POLYEDRE(i)%face(2,k),&
                                        'ver3=',TAB_POLYEDRE(i)%face(3,k)
    end do

    write(128,'(A6)')                                        '$$$$$$'
    H_temp = max(H_temp,TAB_POLYEDRE(i)%centre(3))
    Volume_grain = Volume_grain + 4*3.1415926*((TAB_POLYEDRE(i)%Rmax)**2)/3
  end do

  if (H_temp<H_ini) then
    ep = H_ini-H_temp
    if (ep>dmean) then
      do i=1,nbParticules
        if (TAB_POLYEDRE(i)%behav /= 'plexx')  cycle
        if (TAB_POLYEDRE(i)%centre(3)>H_temp-ep) then
          cpt = cpt + 1
          write(128,'(A6)')                                        '$bdyty'
          write(128,'(1X,A5,1X,i6)')                               'RBDY3',cpt
          write(128,'(A6)')                                        '$blmty'
          write(128,'(1X,A5,5X,i2,2X,A5,2X,A5,1X,A6,D14.7)')       'PLAIN',1,'behav','plexx','avrd=',0.D0
          write(128,'(29X,3(A5,D14.7,2X))')  'I1  =',TAB_POLYEDRE(i)%I1,'I2  =',TAB_POLYEDRE(i)%I2,'I3  =',TAB_POLYEDRE(i)%I3
          write(128,'(A6)')                                        '$nodty'
          write(128,'(1X,A5,5X,i2,16X,3(A5,D14.7,2X))')            'NO6xx',1,'coo1=',TAB_POLYEDRE(i)%centre(1),&
                                                                             'coo2=',TAB_POLYEDRE(i)%centre(2),&
                                                                             'coo3=',TAB_POLYEDRE(i)%centre(3)+ep
          write(128,'(29X,3(A5,D14.7,2X))') 'coo4=',0.D0,'coo5=',0.D0,'coo6=',0.D0
          write(128,'(A6)')                                        '$tacty'
          write(128,'(1X,A5,5X,i2,2X,A5,2X,A5,2X,A10,i7,4X,A9,i7)') 'POLYR',1,'color','TATAx',&
                                                                    'nb_vertex=',TAB_POLYEDRE(i)%nb_vertex,&
                                                                    'nb_faces=',TAB_POLYEDRE(i)%nb_face
          do k=1,TAB_POLYEDRE(i)%nb_vertex
            write(128,'(29X,3(A5,D14.7,2X))') 'coo1=',TAB_POLYEDRE(i)%vertex(1,k)-TAB_POLYEDRE(i)%centre(1),&
                                              'coo2=',TAB_POLYEDRE(i)%vertex(2,k)-TAB_POLYEDRE(i)%centre(2),&
                                              'coo3=',TAB_POLYEDRE(i)%vertex(3,k)-TAB_POLYEDRE(i)%centre(3)
          end do
          do k=1,TAB_POLYEDRE(i)%nb_face
             write(128,'(29X,3(A5,i7,9X))')   'ver1=',TAB_POLYEDRE(i)%face(1,k),&
                                              'ver2=',TAB_POLYEDRE(i)%face(2,k),&
                                              'ver3=',TAB_POLYEDRE(i)%face(3,k)
          end do
          write(128,'(A6)')                                        '$$$$$$'
          Volume_grain = Volume_grain + 4*3.1415926*((TAB_POLYEDRE(i)%Rmax)**2)/3
          H_temp2 = min(H_temp2,TAB_POLYEDRE(i)%centre(3)+ep)
        end if
      end do
!!! POUR 08 FACES ON DOIT ENCORE ET ENCORE RAJOUTER DES GRAINS...
!      do i=1,nbParticules
!!        if (TAB_POLYEDRE(i)%behav /= 'plexx')  cycle
!!        if (TAB_POLYEDRE(i)%centre(3)>H_temp .and. TAB_POLYEDRE(i)%centre(3)<H_temp2 ) then
!          cpt = cpt + 1
!          write(128,'(A6)')                                        '$bdyty'
!          write(128,'(1X,A5,1X,i6)')                               'RBDY3',cpt
!          write(128,'(A6)')                                        '$blmty'
 !         write(128,'(1X,A5,5X,i2,2X,A5,2X,A5,1X,A6,D14.7)')       'PLAIN',1,'behav','plexx','avrd=',0.D0
 !         write(128,'(29X,3(A5,D14.7,2X))')  'I1  =',TAB_POLYEDRE(i)%I1,'I2  =',TAB_POLYEDRE(i)%I2,'I3  =',TAB_POLYEDRE(i)%I3
 !         write(128,'(A6)')                                        '$nodty'
 !         write(128,'(1X,A5,5X,i2,16X,3(A5,D14.7,2X))')            'NO6xx',1,'coo1=',TAB_POLYEDRE(i)%centre(1),&
 !                                                                            'coo2=',TAB_POLYEDRE(i)%centre(2),&
 !                                                                            'coo3=',TAB_POLYEDRE(i)%centre(3)+H_temp
 !         write(128,'(29X,3(A5,D14.7,2X))') 'coo4=',0.D0,'coo5=',0.D0,'coo6=',0.D0
 !         write(128,'(A6)')                                        '$tacty'
 !         write(128,'(1X,A5,5X,i2,2X,A5,2X,A5,2X,A10,i7,4X,A9,i7)') 'POLYR',1,'color','TATAx',&
 !                                                                   'nb_vertex=',TAB_POLYEDRE(i)%nb_vertex,&
 !                                                                   'nb_faces=',TAB_POLYEDRE(i)%nb_face
 !         do k=1,TAB_POLYEDRE(i)%nb_vertex
 !           write(128,'(29X,3(A5,D14.7,2X))') 'coo1=',TAB_POLYEDRE(i)%vertex(1,k)-TAB_POLYEDRE(i)%centre(1),&
 !                                             'coo2=',TAB_POLYEDRE(i)%vertex(2,k)-TAB_POLYEDRE(i)%centre(2),&
 !                                             'coo3=',TAB_POLYEDRE(i)%vertex(3,k)-TAB_POLYEDRE(i)%centre(3)
 !         end do
 !         do k=1,TAB_POLYEDRE(i)%nb_face
 !            write(128,'(29X,3(A5,i7,9X))')   'ver1=',TAB_POLYEDRE(i)%face(1,k),&
 !                                             'ver2=',TAB_POLYEDRE(i)%face(2,k),&
 !                                             'ver3=',TAB_POLYEDRE(i)%face(3,k)
 !         end do
 !         write(128,'(A6)')                                        '$$$$$$'
 !         Volume_grain = Volume_grain + 4*3.1415926*((TAB_POLYEDRE(i)%Rmax)**2)/3
!        end if
 !     end do
    end if
  end if

  nbPart = cpt
  nb_wallh = 0
  nb_wallb = 0
  
  do i=1,nbParticules
    if (TAB_POLYEDRE(i)%behav /= 'wallh')  cycle
    cpt = cpt + 1
    nb_wallh = nb_wallh + 1
    
    write(128,'(A6)')                                        '$bdyty'
    write(128,'(1X,A5,1X,i6)')                               'RBDY3',cpt
    write(128,'(A6)')                                        '$blmty'
    write(128,'(1X,A5,5X,i2,2X,A5,2X,A5,1X,A6,D14.7)')       'PLAIN',1,'behav',TAB_POLYEDRE(i)%behav,'avrd=',0.D0
    write(128,'(29X,3(A5,D14.7,2X))')  'I1  =',TAB_POLYEDRE(i)%I1,'I2  =',TAB_POLYEDRE(i)%I2,'I3  =',TAB_POLYEDRE(i)%I3
    write(128,'(A6)')                                        '$nodty'
    write(128,'(1X,A5,5X,i2,16X,3(A5,D14.7,2X))')            'NO6xx',1,'coo1=',TAB_POLYEDRE(i)%centre(1),&
                                                                       'coo2=',TAB_POLYEDRE(i)%centre(2),&
                                                                       'coo3=',TAB_POLYEDRE(i)%centre(3)+ep
    write(128,'(29X,3(A5,D14.7,2X))') 'coo4=',0.D0,'coo5=',0.D0,'coo6=',0.D0
    write(128,'(A6)')                                        '$tacty'
    write(128,'(1X,A5,5X,i2,2X,A5,2X,A5,2X,A10,i7,4X,A9,i7)') 'POLYR',1,'color','TATAh',&
                                                              'nb_vertex=',TAB_POLYEDRE(i)%nb_vertex,&
                                                              'nb_faces=',TAB_POLYEDRE(i)%nb_face
    do k=1,TAB_POLYEDRE(i)%nb_vertex
      write(128,'(29X,3(A5,D14.7,2X))') 'coo1=',TAB_POLYEDRE(i)%vertex(1,k)-TAB_POLYEDRE(i)%centre(1),&
                                        'coo2=',TAB_POLYEDRE(i)%vertex(2,k)-TAB_POLYEDRE(i)%centre(2),&
                                        'coo3=',TAB_POLYEDRE(i)%vertex(3,k)-TAB_POLYEDRE(i)%centre(3)
    end do
    do k=1,TAB_POLYEDRE(i)%nb_face
       write(128,'(29X,3(A5,i7,9X))')   'ver1=',TAB_POLYEDRE(i)%face(1,k),&
                                        'ver2=',TAB_POLYEDRE(i)%face(2,k),&
                                        'ver3=',TAB_POLYEDRE(i)%face(3,k)
    end do

    write(128,'(A6)')                                        '$$$$$$'
  end do

  cpt = cpt + 1
  Volume_wall = 0.D0
  write(128,'(A6)')                                        '$bdyty'
  write(128,'(1X,A5,1X,i6)')                               'RBDY3',cpt
  write(128,'(A6)')                                        '$blmty'
  write(128,'(1X,A5,2X,I5,2X,A5,2X,A5,2(2X,A5,D14.7))')    'PLAIN',1,'behav','wallb','avrd=',0.D0,'gyrd=',0.D0
  write(128,'(A6)')                                        '$nodty'
  write(128,'(1X,A5,5X,i2,16X,3(A5,D14.7,2X))')            'NO6xx',1,'coo1=',0.5/2,&
                                                                     'coo2=',0.6/2,&
                                                                     'coo3=',0.D0
  write(128,'(29X,3(A5,D14.7,2X))') 'coo4=',0.D0,'coo5=',0.D0,'coo6=',0.D0
  write(128,'(A6)')                                        '$tacty'
  do i=1,nbParticules
    if (TAB_POLYEDRE(i)%behav /= 'wallb')  cycle  
    nb_wallb = nb_wallb + 1
    write(128,'(1X,A5,5X,i2,2X,A5,2X,A5,2X,A10,i7,4X,A9,i7)') 'POLYR',1,'color','TATAb',&
                                                              'nb_vertex=',TAB_POLYEDRE(i)%nb_vertex,&
                                                              'nb_faces=',TAB_POLYEDRE(i)%nb_face
    do k=1,TAB_POLYEDRE(i)%nb_vertex
      write(128,'(29X,3(A5,D14.7,2X))') 'coo1=',TAB_POLYEDRE(i)%vertex(1,k)-0.5/2,&
                                        'coo2=',TAB_POLYEDRE(i)%vertex(2,k)-0.6/2,&
                                        'coo3=',TAB_POLYEDRE(i)%vertex(3,k)
    end do
    do k=1,TAB_POLYEDRE(i)%nb_face
       write(128,'(29X,3(A5,i7,9X))')   'ver1=',TAB_POLYEDRE(i)%face(1,k),&
                                        'ver2=',TAB_POLYEDRE(i)%face(2,k),&
                                        'ver3=',TAB_POLYEDRE(i)%face(3,k)
    end do
    Volume_wall = Volume_wall + 4*3.1415926*((TAB_POLYEDRE(i)%Rmax)**2)/3
  enddo
  print*, 'Masse volumique equivalente = ', Volume_grain * 2800 / Volume_wall

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
    write(129,'(1X,A5,2X,I5,6(1X,D14.7))') 'vlocy',2,0.05,0.,0.,0.,1.,0.
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




  
!===================================================================================================
!Calcul de la deformation
!==================================================================================================
subroutine Calcul_deformation
implicit none
integer                          :: i,j
real(kind=8)                     :: Xplan_max,Xplan_min,Yplan_max,Yplan_min,Zplan_max,Zplan_min
real(kind=8)                     :: Long,Larg,Haut,E1,E2,E3,dmean,masse,Calcul_I,DE1,DE2,DE3


if (deformation_ini) then
  Xplan_max =-10000000
  Xplan_min = 10000000
  Yplan_max =-10000000
  Yplan_min = 10000000
  Zplan_max =-10000000
  Zplan_min = 10000000
  do i=1,nbParticules
    if (TAB_POLYEDRE(i)%behav/='plexx') cycle
    
    Xplan_max = max( Xplan_max, TAB_POLYEDRE(i)%centre(1) + TAB_POLYEDRE(i)%Rmax )
    Xplan_min = min( Xplan_min, TAB_POLYEDRE(i)%centre(1) - TAB_POLYEDRE(i)%Rmax)

    Yplan_max = max( Yplan_max, TAB_POLYEDRE(i)%centre(2) + TAB_POLYEDRE(i)%Rmax )
    Yplan_min = min( Yplan_min, TAB_POLYEDRE(i)%centre(2) - TAB_POLYEDRE(i)%Rmax )

    Zplan_max = max( Zplan_max, TAB_POLYEDRE(i)%centre(3) + TAB_POLYEDRE(i)%Rmax )
    Zplan_min = min( Zplan_min, TAB_POLYEDRE(i)%centre(3) - TAB_POLYEDRE(i)%Rmax )
  end do
  Long_ini = Xplan_max-Xplan_min
  Larg_ini = Yplan_max-Yplan_min
  H_ini    = Zplan_max-Zplan_min
  
  Depsilon_p = 0.D0
  Depsilon_q = 0.D0

  Long_avant = Long_ini
  Larg_avant = Larg_ini
  Haut_avant = H_ini
  
  deformation_ini=.false.
end if


Xplan_max =-10000000
Xplan_min = 10000000
Yplan_max =-10000000
Yplan_min = 10000000
Zplan_max =-10000000
Zplan_min = 10000000
do i=1,nbParticules
  if (TAB_POLYEDRE(i)%behav/='plexx') cycle
    
  Xplan_max = max( Xplan_max, TAB_POLYEDRE(i)%centre(1) + TAB_POLYEDRE(i)%Rmax )
  Xplan_min = min( Xplan_min, TAB_POLYEDRE(i)%centre(1) - TAB_POLYEDRE(i)%Rmax)

  Yplan_max = max( Yplan_max, TAB_POLYEDRE(i)%centre(2) + TAB_POLYEDRE(i)%Rmax )
  Yplan_min = min( Yplan_min, TAB_POLYEDRE(i)%centre(2) - TAB_POLYEDRE(i)%Rmax )

  Zplan_max = max( Zplan_max, TAB_POLYEDRE(i)%centre(3) + TAB_POLYEDRE(i)%Rmax )
  Zplan_min = min( Zplan_min, TAB_POLYEDRE(i)%centre(3) - TAB_POLYEDRE(i)%Rmax )
end do

Long = 0.5 !Xplan_max-Xplan_min
Larg = 0.6 !Yplan_max-Yplan_min
Haut = Zplan_max -Zplan_min

largeur  = 0.5 !Xplan_max-Xplan_min
longueur = 0.6 !Yplan_max-Yplan_min
Hauteur  = Zplan_max -Zplan_min

Def = temps * 0.05 / Hauteur


end subroutine Calcul_deformation

!==================================================================================================
!On ferme tout les fichiers
!=================================================================================================
subroutine fermer
implicit none
integer           :: i

print*,'---> fermeture des fichiers'

if (calcul_coordination               == 1)    close (100)
if (calcul_vitesse_moyenne            == 1)    close (101)
if (calcul_qsurp                      == 1)    close (102)
if (calcul_compacity                  == 1)    close (103)
if (contact_anisotropes               == 1)    close (104)
if (branche_anisotropes               == 1)    close (105)
do,i=1,nbparticules
  deallocate(TAB_POLYEDRE(i)%vertex_ref)
  deallocate(TAB_POLYEDRE(i)%vertex)
  deallocate(TAB_POLYEDRE(i)%face)
end do
deallocate(TAB_POLYEDRE)

print*,'--> Fin du postraitement <--'
stop

end subroutine fermer



!==============================================================================
!==============================================================================
!==============================================================================
!==============================================================================




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
