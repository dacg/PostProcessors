program Construct_polyr

implicit none
  
type  ::  T_POLY
  integer                                        ::  nb_vertex,nb_face
  real(kind=8),dimension(:,:),allocatable        ::  vertex
  integer,dimension(:,:),allocatable             ::  face
  logical,dimension(:),allocatable               ::  select_face
  real(kind=8),dimension(:,:),allocatable        ::  normal_face
  real(kind=8)                                   ::  Rmax
  real(kind=8),dimension(3)                      ::  centre
end type T_POLY
type(T_POLY),dimension(:),allocatable            ::  TAB_POLYEDRE_new
type(T_POLY),dimension(:),allocatable            ::  TAB_POLYEDRE_old
type(T_POLY),dimension(:),allocatable            ::  TAB_POLYEDRE_centre

integer                                          ::  nb_vertex
integer                                          ::  nb_polyedre
real(kind=8)                                     ::  HL_max,HL_min,LL_max,LL_min

integer                                          ::  K_poly
real(kind=8),dimension(:,:),allocatable          ::  vertex_List
real(kind=8)                                     ::  Total_ang_moyen = 0
real(kind=8)                                     ::  Total_var_ang_moyen = 0
real(kind=8)                                     ::  Total_LongH_moyen = 0
real(kind=8)                                     ::  Total_LargH_moyen = 0
real(kind=8)                                     ::  a_12_moyen=0,a_13_moyen=0,a_23_moyen=0
logical                                          ::  initial_position,write_gmv


!============================================================================================!
!============================================================================================!
!============================================================================================!
!                       Generation alleatoire de grains polyédriques                         !
!                                                                                            !
!  Programme generant un nombre "nb_polyedre" de polyedre comportant chacun un nombre        !
!  "nb_vertex" de sommets. Par construction les polyedres ont des faces triangulaires ce qui !
!  impose un nombre de face "nb_face = 2*nb_vertex - 4".                                     !
!  Les sommets sont tirés alleatoirement sur une ellipsoide de grand axe L=1, d'axe moyen H  !
!  et de petit axe l tels que :                                                              !
!        1) HL_min < (H/L) < HL_max    et    LL_min < (l/L) < LL_max                         !
!        2) l  <  H  <  L                                                                    !
!  Si un BODIES.DAT de spheres existe deja, il est possible d injecter les polyedres dedans  !
!  en posant "initial_position= .true." ; sinon tout les polyedres aurant pour diametre      !
!  maximun L = 1 et seront tous centrés en (0,0,0).                                          !
!  Tous les polyedres sont enregistres dans un fichier BODIES_POLYR.DAT lisible par LMGC90   !
!  En posant "write_gmv = .true." un fichier gmv sera cree pour chaque polyedre              !
!  Un fichier BODIES_PARAMETER resume pour chaque polyedre les donnees suivantes :           !
!        Colonne 1 = numero du polyedre                                                      !
!        Colonne 2 = Angularite du polyedre (K_ang = moyenne des angles entre les normales   ! 
!                                            de chaque face en contactes normalises par Pi   !
!                                            K_ang->0 avec le nombre de sommets)             !
!        Colonne 3 = la variance de l angularite                                             !
!        Colonne 4 = le rapport de l aire du polyedre sur l aire de la sphere circonscrite   !
!        Colonne 5 = L/H                                                                     !
!        Colonne 6 = l/L                                                                     !
!        A la fin du fichier, toutes ces grandeurs sont moyennees sur tout les polyedres     !
!============================================================================================!
!============================================================================================!
!============================================================================================!


!============================================================================================!
!============================================================================================!
!                               COMMANDES A REMPLIR ICI                                      !

                                nb_vertex        = 25
                                nb_polyedre      = 10
                                HL_max           = 1.0
                                LL_max           = 1.0
                                HL_min           = 0.94
                                LL_min           = 0.94
                                initial_position = .false.
                                write_gmv        = .true.

!============================================================================================!
!============================================================================================!







!============================================================================================!
!============================================================================================!
!                               DEBUT DU PROGRAMME                                           !
!============================================================================================!
!============================================================================================!



!-----------y a t il un fichier de sphere existant?-------
call read_initial_position

do K_poly=1,nb_polyedre
!------------- Construction d'une liste aléatoire de vertex----------
  allocate(vertex_List(nb_vertex,3))
  call random_List_vertex(nb_vertex,2)
!------------- Construction du polyedre (tetraedre) origine----------
  call INIT_POLY(K_poly,nb_vertex)
!------------- Construction du polyedre general----------
  call construct_polyedre(K_poly,nb_vertex)
!--------------Ecriture des donnees----------------------
  if (write_gmv) call write_new_gmv(K_poly,nb_polyedre)
  call write_bodies_aire(K_poly)
  call write_bodies(K_poly)
  deallocate(vertex_List)
end do
  call write_plan

stop

contains
!=========================================================================
! Generation du polyedre 
!=========================================================================
subroutine construct_polyedre(K_poly_ini,nb_vertex_ini)
implicit none
integer                                     :: i,j,k,K_poly_ini,nb_vertex_ini,cpt_visible,cpt,nb_new_face
real(kind=8),dimension(3)                   :: New_point
real(kind=8),dimension(3)                   :: a,b,c,g,ab,ac,bc,w
real(kind=8)                                :: GO_w,OP_w,norm,d
integer,dimension(:),allocatable            :: Liste_vertex
integer,dimension(:,:),allocatable          :: potential_edges,new_face

if (nb_vertex==4) then
  TAB_POLYEDRE_new(K_poly_ini)%nb_vertex = nb_vertex
  if (allocated(TAB_POLYEDRE_new(K_poly_ini)%vertex)) deallocate(TAB_POLYEDRE_new(K_poly_ini)%vertex)
  allocate(TAB_POLYEDRE_new(K_poly_ini)%vertex(nb_vertex,3))
  TAB_POLYEDRE_new(K_poly_ini)%vertex(1:nb_vertex,:) = TAB_POLYEDRE_old(K_poly_ini)%vertex(1:nb_vertex,:)

  TAB_POLYEDRE_new(K_poly_ini)%nb_face = TAB_POLYEDRE_old(K_poly_ini)%nb_face
  if (allocated(TAB_POLYEDRE_new(K_poly_ini)%face)) deallocate(TAB_POLYEDRE_new(K_poly_ini)%face)
  allocate(TAB_POLYEDRE_new(K_poly_ini)%face(TAB_POLYEDRE_new(K_poly_ini)%nb_face,3))
  TAB_POLYEDRE_new(K_poly_ini)%face = 0
  cpt = 0
  do j = 1,TAB_POLYEDRE_old(K_poly_ini)%nb_face
    TAB_POLYEDRE_new(K_poly_ini)%face(j,:) = TAB_POLYEDRE_old(K_poly_ini)%face(j,:)
  end do

else

  do i=5,nb_vertex_ini
    New_point(1) = vertex_List(i,1)
    New_point(2) = vertex_List(i,2)
    New_point(3) = vertex_List(i,3)  
  
!----On rempli le nouveau polyedre des nouveaux vertex---------------------
    TAB_POLYEDRE_new(K_poly_ini)%nb_vertex = i 
    if (allocated(TAB_POLYEDRE_new(K_poly_ini)%vertex)) deallocate(TAB_POLYEDRE_new(K_poly_ini)%vertex)
    allocate(TAB_POLYEDRE_new(K_poly_ini)%vertex(i,3))
    TAB_POLYEDRE_new(K_poly_ini)%vertex(1:i-1,:) = TAB_POLYEDRE_old(K_poly_ini)%vertex(1:i,:)
    TAB_POLYEDRE_new(K_poly_ini)%vertex(i,:) = New_point(:)
    TAB_POLYEDRE_new(K_poly_ini)%centre(:)=0
    do j=1,TAB_POLYEDRE_new(K_poly_ini)%nb_vertex
      TAB_POLYEDRE_new(K_poly_ini)%centre(1)=TAB_POLYEDRE_new(K_poly_ini)%centre(1)+TAB_POLYEDRE_new(K_poly_ini)%vertex(j,1)
      TAB_POLYEDRE_new(K_poly_ini)%centre(2)=TAB_POLYEDRE_new(K_poly_ini)%centre(2)+TAB_POLYEDRE_new(K_poly_ini)%vertex(j,2)
      TAB_POLYEDRE_new(K_poly_ini)%centre(3)=TAB_POLYEDRE_new(K_poly_ini)%centre(3)+TAB_POLYEDRE_new(K_poly_ini)%vertex(j,3)
    enddo
    TAB_POLYEDRE_new(K_poly_ini)%centre=TAB_POLYEDRE_new(K_poly_ini)%centre/real(TAB_POLYEDRE_new(K_poly_ini)%nb_vertex,8)

!-----Il nous faut orienter chacune des faces du old polyedre, chaque face pointant vers l'exterieur----------
    if (allocated(TAB_POLYEDRE_old(K_poly_ini)%normal_face)) deallocate(TAB_POLYEDRE_old(K_poly_ini)%normal_face)
    allocate(TAB_POLYEDRE_old(K_poly_ini)%normal_face(TAB_POLYEDRE_old(K_poly_ini)%nb_face,3))
    do j = 1,TAB_POLYEDRE_old(K_poly_ini)%nb_face
      a(1)=TAB_POLYEDRE_old(K_poly_ini)%vertex(TAB_POLYEDRE_old(K_poly_ini)%face(j,1),1)
      a(2)=TAB_POLYEDRE_old(K_poly_ini)%vertex(TAB_POLYEDRE_old(K_poly_ini)%face(j,1),2)
      a(3)=TAB_POLYEDRE_old(K_poly_ini)%vertex(TAB_POLYEDRE_old(K_poly_ini)%face(j,1),3)
      b(1)=TAB_POLYEDRE_old(K_poly_ini)%vertex(TAB_POLYEDRE_old(K_poly_ini)%face(j,2),1)
      b(2)=TAB_POLYEDRE_old(K_poly_ini)%vertex(TAB_POLYEDRE_old(K_poly_ini)%face(j,2),2)
      b(3)=TAB_POLYEDRE_old(K_poly_ini)%vertex(TAB_POLYEDRE_old(K_poly_ini)%face(j,2),3)
      c(1)=TAB_POLYEDRE_old(K_poly_ini)%vertex(TAB_POLYEDRE_old(K_poly_ini)%face(j,3),1)
      c(2)=TAB_POLYEDRE_old(K_poly_ini)%vertex(TAB_POLYEDRE_old(K_poly_ini)%face(j,3),2)
      c(3)=TAB_POLYEDRE_old(K_poly_ini)%vertex(TAB_POLYEDRE_old(K_poly_ini)%face(j,3),3)
      g(1)=TAB_POLYEDRE_old(K_poly_ini)%centre(1)
      g(2)=TAB_POLYEDRE_old(K_poly_ini)%centre(2)
      g(3)=TAB_POLYEDRE_old(K_poly_ini)%centre(3)
      
    !-------------Construction des vecteurs-------------------------------
      ab = b-a
      ac = c-a
      bc = c-b
    !-------------Calcul de la normale (produit vectoriel)-------------------------------
      w(1)=ab(2)*ac(3)-ab(3)*ac(2)
      w(2)=-ab(1)*ac(3)+ab(3)*ac(1)
      w(3)=ab(1)*ac(2)-ab(2)*ac(1)
      norm = sqrt(w(1)**2+w(2)**2+w(3)**2)
      w = w / norm
    !-------------Calcul du coefficient d-------------------------------
      d = -(w(1)*a(1)+w(2)*a(2)+w(3)*a(3))

    !-------------Calcul du produit scalaire de GO.w, sachant que O appartient au plan --------------------
      GO_w = -( d + g(1)*w(1) + g(2)*w(2) + g(3)*w(3) )

      if ( GO_w < 0 ) w = - w

      TAB_POLYEDRE_old(K_poly_ini)%normal_face(j,1) = w(1)
      TAB_POLYEDRE_old(K_poly_ini)%normal_face(j,2) = w(2)
      TAB_POLYEDRE_old(K_poly_ini)%normal_face(j,3) = w(3)
    end do
  
!-----Maintenant on cherche quelle sont les faces qui voient le nouveau point qu'on ajoute-----------------
!-----Il faut calculer le produit scalaire OP.w. s il est negatif alors oui la face est vue-----
    if (allocated(TAB_POLYEDRE_old(K_poly_ini)%select_face)) deallocate(TAB_POLYEDRE_old(K_poly_ini)%select_face)
    allocate(TAB_POLYEDRE_old(K_poly_ini)%select_face(TAB_POLYEDRE_old(K_poly_ini)%nb_face))
    cpt_visible     = 0
    do j = 1,TAB_POLYEDRE_old(K_poly_ini)%nb_face
      d    = -( TAB_POLYEDRE_old(K_poly_ini)%vertex(TAB_POLYEDRE_old(K_poly_ini)%face(j,1),1)*&
                TAB_POLYEDRE_old(K_poly_ini)%normal_face(j,1) +&
                TAB_POLYEDRE_old(K_poly_ini)%vertex(TAB_POLYEDRE_old(K_poly_ini)%face(j,1),2)*&
                TAB_POLYEDRE_old(K_poly_ini)%normal_face(j,2) +& 
                TAB_POLYEDRE_old(K_poly_ini)%vertex(TAB_POLYEDRE_old(K_poly_ini)%face(j,1),3)*&
                TAB_POLYEDRE_old(K_poly_ini)%normal_face(j,3) )
      OP_w = -( d + &
                New_point(1)*TAB_POLYEDRE_old(K_poly_ini)%normal_face(j,1) +&
                New_point(2)*TAB_POLYEDRE_old(K_poly_ini)%normal_face(j,2) +&
                New_point(3)*TAB_POLYEDRE_old(K_poly_ini)%normal_face(j,3) )
    
      if ( OP_w < 0 ) then
        TAB_POLYEDRE_old(K_poly_ini)%select_face(j) = .true.
        cpt_visible = cpt_visible + 1
      else
        TAB_POLYEDRE_old(K_poly_ini)%select_face(j) = .false.
      end if
    end do  

!---On va remplir le new polyedre des faces qui ne sont pas vues, autrement dit celles qui restent------
    TAB_POLYEDRE_new(K_poly_ini)%nb_face = TAB_POLYEDRE_old(K_poly_ini)%nb_face + 2
    if (allocated(TAB_POLYEDRE_new(K_poly_ini)%face)) deallocate(TAB_POLYEDRE_new(K_poly_ini)%face)
    allocate(TAB_POLYEDRE_new(K_poly_ini)%face(TAB_POLYEDRE_new(K_poly_ini)%nb_face,3))
    TAB_POLYEDRE_new(K_poly_ini)%face = 0
    cpt = 0
    do j = 1,TAB_POLYEDRE_old(K_poly_ini)%nb_face
      if ( .NOT. TAB_POLYEDRE_old(K_poly_ini)%select_face(j) ) then
        cpt = cpt + 1
        TAB_POLYEDRE_new(K_poly_ini)%face(cpt,:) = TAB_POLYEDRE_old(K_poly_ini)%face(j,:)
      end if
    end do

!---Il faut maintenant construire les nouvelles faces---
!---Il nous faut d'abord la liste de tout les cotes possibles---
!---Une face c est 3 sommets, c est donc 3 cotes potentiels.---

    if (allocated(potential_edges)) deallocate(potential_edges)
    allocate(potential_edges(3*cpt_visible,2))
    cpt = 0
    do j = 1,TAB_POLYEDRE_old(K_poly_ini)%nb_face
      if ( TAB_POLYEDRE_old(K_poly_ini)%select_face(j) ) then
        cpt = cpt + 1
        potential_edges(cpt,1) = TAB_POLYEDRE_old(K_poly_ini)%face(j,1)
        potential_edges(cpt,2) = TAB_POLYEDRE_old(K_poly_ini)%face(j,2)

        cpt = cpt + 1
        potential_edges(cpt,1) = TAB_POLYEDRE_old(K_poly_ini)%face(j,1)
        potential_edges(cpt,2) = TAB_POLYEDRE_old(K_poly_ini)%face(j,3)

        cpt = cpt + 1
        potential_edges(cpt,1) = TAB_POLYEDRE_old(K_poly_ini)%face(j,2)
        potential_edges(cpt,2) = TAB_POLYEDRE_old(K_poly_ini)%face(j,3)
      end if
    end do

!---Dans cette liste de cotes potentiels, il apparait des cotes qui sont en fait
!---interieur au polyedre. Il faut donc les eliminer. ---------------------------------
!---Si on a N faces visibles, on forme N+2 vrai nouveaux cotes------------------------
    if (allocated(new_face)) deallocate(new_face)
    allocate(new_face(cpt_visible+2,3))
    nb_new_face = 0
    do j = 1,3*cpt_visible
      cpt = 0
      do k = 1,3*cpt_visible
        if ( (j==k) ) cycle
        if ( ((potential_edges(j,1)==potential_edges(k,1)).and.(potential_edges(j,2)==potential_edges(k,2))).or.&
             ((potential_edges(j,1)==potential_edges(k,2)).and.(potential_edges(j,2)==potential_edges(k,1))) ) then
          cpt = cpt + 1
        end if
      end do
      if ( cpt==0 ) then
        nb_new_face = nb_new_face + 1
        new_face(nb_new_face,1) = potential_edges(j,1)
        new_face(nb_new_face,2) = potential_edges(j,2)
        new_face(nb_new_face,3) = i
      end if 
    end do

!---Il ne reste plus qu a completer le New polyedre par ces new faces ------
    cpt = 0
    do j = 1,TAB_POLYEDRE_new(K_poly_ini)%nb_face
      if ( (TAB_POLYEDRE_new(K_poly_ini)%face(j,1)==0) .and. &
           (TAB_POLYEDRE_new(K_poly_ini)%face(j,2)==0) .and. &
           (TAB_POLYEDRE_new(K_poly_ini)%face(j,3)==0) ) then
    
        cpt = cpt + 1
        TAB_POLYEDRE_new(K_poly_ini)%face(j,1) = new_face(cpt,1)
        TAB_POLYEDRE_new(K_poly_ini)%face(j,2) = new_face(cpt,2)
        TAB_POLYEDRE_new(K_poly_ini)%face(j,3) = new_face(cpt,3)
      end if
    end do

!---Une fois l iteration finie, on dit que old=new ------
    TAB_POLYEDRE_old(K_poly_ini)%nb_vertex = TAB_POLYEDRE_new(K_poly_ini)%nb_vertex
    if (allocated(TAB_POLYEDRE_old(K_poly_ini)%vertex)) deallocate(TAB_POLYEDRE_old(K_poly_ini)%vertex)
    allocate(TAB_POLYEDRE_old(K_poly_ini)%vertex(TAB_POLYEDRE_old(K_poly_ini)%nb_vertex,3))
    TAB_POLYEDRE_old(K_poly_ini)%vertex = TAB_POLYEDRE_new(K_poly_ini)%vertex
  
    TAB_POLYEDRE_old(K_poly_ini)%nb_face = TAB_POLYEDRE_new(K_poly_ini)%nb_face
    if (allocated(TAB_POLYEDRE_old(K_poly_ini)%face)) deallocate(TAB_POLYEDRE_old(K_poly_ini)%face)
    allocate(TAB_POLYEDRE_old(K_poly_ini)%face(TAB_POLYEDRE_old(K_poly_ini)%nb_face,3))
    TAB_POLYEDRE_old(K_poly_ini)%face = TAB_POLYEDRE_new(K_poly_ini)%face
  end do
end if

end subroutine construct_polyedre


!=========================================================================
! Generation du polyedre initial
!=========================================================================
subroutine INIT_POLY(K_poly_ini,nb_vertex_ini)
implicit none
integer                                     :: i,K_poly_ini,nb_vertex_ini

if (allocated(TAB_POLYEDRE_old)) deallocate(TAB_POLYEDRE_old)
allocate(TAB_POLYEDRE_old(nb_polyedre))
if (allocated(TAB_POLYEDRE_new)) deallocate(TAB_POLYEDRE_new)
allocate(TAB_POLYEDRE_new(nb_polyedre))

TAB_POLYEDRE_old(K_poly_ini)%nb_vertex = 4
if (allocated(TAB_POLYEDRE_old(K_poly_ini)%vertex)) deallocate(TAB_POLYEDRE_old(K_poly_ini)%vertex)
allocate(TAB_POLYEDRE_old(K_poly_ini)%vertex(4,3))
do i=1,TAB_POLYEDRE_old(K_poly_ini)%nb_vertex
  TAB_POLYEDRE_old(K_poly_ini)%vertex(i,1) = vertex_List(i,1)
  TAB_POLYEDRE_old(K_poly_ini)%vertex(i,2) = vertex_List(i,2)
  TAB_POLYEDRE_old(K_poly_ini)%vertex(i,3) = vertex_List(i,3)
end do

TAB_POLYEDRE_old(K_poly_ini)%nb_face   = 4
if (allocated(TAB_POLYEDRE_old(K_poly_ini)%face)) deallocate(TAB_POLYEDRE_old(K_poly_ini)%face)
allocate(TAB_POLYEDRE_old(K_poly_ini)%face(4,3))
TAB_POLYEDRE_old(K_poly_ini)%face(1,1) = 1
TAB_POLYEDRE_old(K_poly_ini)%face(1,2) = 2
TAB_POLYEDRE_old(K_poly_ini)%face(1,3) = 3
TAB_POLYEDRE_old(K_poly_ini)%face(2,1) = 1
TAB_POLYEDRE_old(K_poly_ini)%face(2,2) = 3
TAB_POLYEDRE_old(K_poly_ini)%face(2,3) = 4
TAB_POLYEDRE_old(K_poly_ini)%face(3,1) = 4
TAB_POLYEDRE_old(K_poly_ini)%face(3,2) = 3
TAB_POLYEDRE_old(K_poly_ini)%face(3,3) = 2
TAB_POLYEDRE_old(K_poly_ini)%face(4,1) = 1
TAB_POLYEDRE_old(K_poly_ini)%face(4,2) = 2
TAB_POLYEDRE_old(K_poly_ini)%face(4,3) = 4

!----- On cherche le centre du 1er old polyedre
TAB_POLYEDRE_old(K_poly_ini)%centre(:)=0
do i=1,TAB_POLYEDRE_old(K_poly_ini)%nb_vertex
  TAB_POLYEDRE_old(K_poly_ini)%centre(1)=TAB_POLYEDRE_old(K_poly_ini)%centre(1)+TAB_POLYEDRE_old(K_poly_ini)%vertex(i,1)
  TAB_POLYEDRE_old(K_poly_ini)%centre(2)=TAB_POLYEDRE_old(K_poly_ini)%centre(2)+TAB_POLYEDRE_old(K_poly_ini)%vertex(i,2)
  TAB_POLYEDRE_old(K_poly_ini)%centre(3)=TAB_POLYEDRE_old(K_poly_ini)%centre(3)+TAB_POLYEDRE_old(K_poly_ini)%vertex(i,3)
enddo
TAB_POLYEDRE_old(K_poly_ini)%centre=TAB_POLYEDRE_old(K_poly_ini)%centre/real(TAB_POLYEDRE_old(K_poly_ini)%nb_vertex,8)

end subroutine INIT_POLY


!=========================================================================
! Generation d une liste de vertex aleatoire disposes sur la sphere unite
!=========================================================================
subroutine random_List_vertex(nb_vertex_out,type_alleatoire_out)
implicit none
integer                                     :: i,j
character(len=8)                            :: date
character(len=10)                           :: time
integer                                     :: seed
real(kind=8)                                :: seed2
real(kind=8)                                :: x,y,z, p,Dac,dmax_vertex,dmin_vertex,d_vertex,k
real(kind=8)                                :: xmax,xmin,ymax,ymin,zmax,zmin,Haut,Larg,Long,norm,norm_L
real(kind=8)                                :: signe,dist,a_12,a_13,a_23,Long_temp,al_12,al_13,al_23
real(kind=8),dimension(3)                   :: ac,centre_g,r_vertex,centre_p
integer                                     :: nb_vertex_out,type_alleatoire_out,compteur,ver1_long,ver2_long
real(kind=8),dimension(3,3)                 :: fabric,fabric_L

compteur = 0

if (type_alleatoire_out==1) then
  do i=1,nb_vertex_out
    do 
      call date_and_time(date,time) 
      call random_number(seed2)
      read(time(1:6),'(I6)') seed
      x = cos(seed*seed2)

      call date_and_time(date,time) 
      call random_number(seed2)
      read(time(1:6),'(I6)') seed
      y = cos(seed*seed2)
    
      ! if (x**2+y**2 < 1.D0) exit
      if (x**2 + (y**2)/((LL_max-LL_min)**2) < 1.D0) exit
    end do 
    call date_and_time(date,time) 
    call random_number(seed2)
    read(time(1:6),'(I6)') seed
    signe = cos(seed*seed2)
!    z = sign(1.D0,signe)*sqrt(1-x**2-y**2)

    z = sign(1.D0,signe)*sqrt( (HL_max-HL_min)**2 - (x**2) *(HL_max-HL_min)**2 - (y**2) * ( (HL_max-HL_min)/(LL_max-LL_min) )**2   )
      
    vertex_List(i,1) = x
    vertex_List(i,2) = y
    vertex_List(i,3) = z
  end do
  
else if (type_alleatoire_out==2) then
1 do i=1,nb_vertex_out
10  do 
      call date_and_time(date,time) 
      call random_number(seed2)
      read(time(1:6),'(I6)') seed
      x = cos(seed*seed2)

      call date_and_time(date,time) 
      call random_number(seed2)
      read(time(1:6),'(I6)') seed
      y = cos(seed*seed2)

      if ( (x**2 + (y**2)/((LL_max)**2) < 1.D0) ) exit
!      if (x**2+y**2 < 1.D0) exit
    end do 
    call date_and_time(date,time) 
    call random_number(seed2)
    read(time(1:6),'(I6)') seed
    signe = cos(seed*seed2)
    z = sign(1.D0,signe)*sqrt( (HL_max)**2 - (x**2) *(HL_max)**2 - (y**2) * ( (HL_max)/(LL_max) )**2   )
    if (i>1) then
      do j=1,i-1
        ac(1)    = vertex_List(j,1) - x
        ac(2)    = vertex_List(j,2) - y
        ac(3)    = vertex_List(j,3) - z
        Dac = sqrt( ac(1)**2 + ac(2)**2 + ac(3)**2 )
        
        p = sqrt( 8 * 3.1415926  / ( real(2*nb_vertex-4) * sin( 3.1415926 / 3)) )
                
        if (Dac<0.45*p) then
          compteur = compteur + 1
          if (compteur==1000) then
            compteur = 0
            go to 1
          end if
          Go to 10
        end if
      end do
    end if
    vertex_List(i,1) = x
    vertex_List(i,2) = y
    vertex_List(i,3) = z
  end do

  centre_g(:) = 0
  do i=1,nb_vertex_out
    centre_g(:) = centre_g(:) + vertex_List(i,:)
  end do
  centre_g(:) = centre_g(:)/nb_vertex_out
  dmax_vertex = 0
  dmin_vertex = 1000000
  do i=1,nb_vertex_out
    d_vertex = sqrt( (centre_g(1)-vertex_List(i,1))**2 + (centre_g(2)-vertex_List(i,2))**2 + (centre_g(3)-vertex_List(i,3))**2 )
    dmin_vertex = min( dmin_vertex,d_vertex )
    dmax_vertex = max( dmax_vertex,d_vertex )
  end do
  if ( (dmin_vertex/dmax_vertex < min(HL_min,LL_min)).or.(dmin_vertex/dmax_vertex > max(HL_max,LL_max)) ) go to 1
  
  Long = 0.D0
  larg = 10000000.D0
  Haut = 0.D0
  do i=1,nb_vertex_out
    do j=1,nb_vertex_out
      if (i==j) cycle
      Long_temp = sqrt( (vertex_List(i,1)-vertex_List(j,1))**2 +&
                        (vertex_List(i,2)-vertex_List(j,2))**2 +&
                        (vertex_List(i,3)-vertex_List(j,3))**2 )
      if (Long_temp>Long) then
        Long = Long_temp
        ver1_long = i
        ver2_long = j
      end if 
    end do
  end do

  ac(:) = vertex_List(ver1_long,:) - vertex_List(ver2_long,:)
  norm  = sqrt( ac(1)**2 + ac(2)**2 + ac(3)**2 )
  ac(:) = ac(:) / norm
  centre_p(:) = (vertex_List(ver1_long,:)+vertex_List(ver2_long,:))/2
  p     = -( ac(1)*centre_p(1) +&
             ac(2)*centre_p(2) + &
             ac(3)*centre_p(3) )
             
  xmax = -1000000.D0
  xmin = 1000000  
  ymax = -1000000.D0
  ymin = 1000000             
  do i=1,nb_vertex_out
    if ( (i==ver1_Long).or.(i==ver2_Long) ) cycle
    k = -( p + ac(1)*vertex_List(i,1) + ac(2)*vertex_List(i,2) + ac(3)*vertex_List(i,3) )
    x = k*ac(1) + vertex_List(i,1)
    y = k*ac(2) + vertex_List(i,2)
    z = k*ac(3) + vertex_List(i,3)
    xmax = max( xmax, x )
    ymax = max( ymax, y )
    xmin = min( xmin, x )
    ymin = min( ymin, y )
  end do
  larg = min(abs(xmax-xmin),abs(ymax-ymin))
  Haut = max(abs(xmax-xmin),abs(ymax-ymin))

  if ( (Haut/Long > HL_max).or.(Haut/Long < HL_min) ) go to 1
  if ( (larg/Long > LL_max).or.(larg/Long < LL_min) ) go to 1
end if 
end subroutine random_List_vertex

!==================================================================================================
!Procedure pour ecrire un gmv du new polyedre
!==================================================================================================
subroutine write_new_gmv(K_poly_ini_ini,nb_part_out)
implicit none
integer                                     :: i,j,k,nb_nodev,ncells
integer                                     :: nb_part_out,K_poly_ini_ini
character(len=22)                           :: clout_DOF

print*,'....  WRITE GMV FILE....'


clout_DOF  =  'GMV_POLYEDRE_FILE.    '

if (K_poly_ini_ini<10) then
  WRITE(clout_DOF(19:20),'(I1)')   K_poly_ini_ini
else if ( (K_poly_ini_ini>=10) .and. (K_poly_ini_ini<100) ) then
  WRITE(clout_DOF(19:21),'(I2)')   K_poly_ini_ini
else if ( (K_poly_ini_ini>=100).and. (K_poly_ini_ini<1000) ) then
  WRITE(clout_DOF(19:22),'(I3)')   K_poly_ini_ini
end if   


nb_nodev  = 0
ncells    = 0

open(unit=30000,file=clout_DOF,status='replace') 
write(30000,'(A14)')      'gmvinput ascii'
!--------------------------------------- Nodev
do i=1,nb_part_out
  do j=1,TAB_POLYEDRE_new(i)%nb_vertex
    nb_nodev = nb_nodev + 1
  end do
end do
write(30000,'(A5,1X,I9)') 'nodev',(nb_nodev)

do i=1,nb_part_out
  do j=1,TAB_POLYEDRE_new(i)%nb_vertex
    write(30000,'(3(1X,E14.7))')  TAB_POLYEDRE_new(i)%vertex(j,1), TAB_POLYEDRE_new(i)%vertex(j,2),TAB_POLYEDRE_new(i)%vertex(j,3)
  end do
enddo 

!--------------------------------------- Cells
do i=1,nb_part_out
  ncells = ncells + TAB_POLYEDRE_new(i)%nb_face
end do
write(30000,'(A5,1X,I9)') 'cells',ncells

!--------------------------------------- Surface
k = 0
do i=1,nb_part_out
  do j=1,TAB_POLYEDRE_new(i)%nb_face
    write(30000,'(A3,1X,I7)') 'tri',3 
    write(30000,'(3(I11))')    TAB_POLYEDRE_new(i)%face(j,1),TAB_POLYEDRE_new(i)%face(j,2),TAB_POLYEDRE_new(i)%face(j,3)
  end do
end do


write(30000,'(A8)') 'polygons'
write(30000,'(A7)') 'endpoly'

write(30000,'(A6)') 'endgmv' 
close(30000)
end subroutine write_new_gmv

!==============================================================================
! lecture des positions initiales et du rayon 
!==============================================================================
subroutine read_initial_position
implicit none
integer                             :: i=0
character(len=6)                    ::  text

print*,'.... Lecture des positions initiales ....'

allocate(TAB_POLYEDRE_centre(nb_polyedre))
if (initial_position) then
  open(unit=130,file='BODIES.DAT',status='old')
  do
    read(130,'(A6)') text
    if (text == '$nodty') then
      i = i + 1
      read(130,'(29X, 3(5x,D14.7,2X))') TAB_POLYEDRE_centre(i)%centre(1),&
                                      TAB_POLYEDRE_centre(i)%centre(2),&
                                      TAB_POLYEDRE_centre(i)%centre(3)
      read(130,*)
      read(130,*)
      read(130,'(34X,D14.7)') TAB_POLYEDRE_centre(i)%Rmax
      if (i==nb_polyedre) exit
    end if
  end do
  close(130)
else
  do i=1,nb_polyedre
    TAB_POLYEDRE_centre(i)%centre(:) = 0
  end do
end if

end subroutine read_initial_position

!==============================================================================
! write Bodies.dat
!==============================================================================
subroutine write_plan
implicit none

  write(128,'(A6)') '$bdyty'
  write(128,'(1X,A5,1X,i6)') 'RBDY3',40001
  write(128,'(A6)') '$blmty'
  write(128,'(1X,A5,5X,i2,2X,A5,2X,A5,1X,A6,D14.7)') 'PLAIN',1,'behav','TETEz','avrd=',0.0000000D+01
  write(128,'(A6)') '$nodty'
  write(128,'(1X,A5,5X,i2,16X,3(A5,D14.7,2X))') &
  'NO6xx',1,'coo1=',0.0000000D+00,'coo2=',0.0000000D+00,'coo3=',-0.6000000D-02
  write(128,'(29X,3(A5,D14.7,2X))')'coo4=',0.D0,'coo5=',0.D0,'coo6=',0.D0
  write(128,'(A6)') '$tacty'
  write(128,'(1X,A5,5X,i2,2X,A5,2X,A5,2X,3(A5,D14.7,2X))') &
  'PLANx',1,'color','WALLx',&
  'axe1=',0.0800000D+01,'axe2=',0.0800000D+01,'axe3=',0.6000000D-02
  write(128,'(A6)') '$$$$$$'

  write(128,'(A6)') '$bdyty'
  write(128,'(1X,A5,1X,i6)') 'RBDY3',40002
  write(128,'(A6)') '$blmty'
  write(128,'(1X,A5,5X,i2,2X,A5,2X,A5,1X,A6,D14.7)') 'PLAIN',1,'behav','TETEz','avrd=',0.0000000D+01
  write(128,'(A6)') '$nodty'
  write(128,'(1X,A5,5X,i2,16X,3(A5,D14.7,2X))') &
  'NO6xx',1,'coo1=',0.0000000D+00,'coo2=',0.0000000D+00,'coo3=',0.6407966D+00
  write(128,'(29X,3(A5,D14.7,2X))')'coo4=',0.D0,'coo5=',0.D0,'coo6=',0.D0
  write(128,'(A6)') '$tacty'
  write(128,'(1X,A5,5X,i2,2X,A5,2X,A5,2X,3(A5,D14.7,2X))') &
  'PLANx',1,'color','WALLx',&
  'axe1=',0.0800000D+01,'axe2=',0.0800000D+01,'axe3=',0.6000000D-02
  write(128,'(A6)') '$$$$$$'


  write(128,'(A6)') '$bdyty'
  write(128,'(1X,A5,1X,i6)') 'RBDY3',40003
  write(128,'(A6)') '$blmty'
  write(128,'(1X,A5,5X,i2,2X,A5,2X,A5,1X,A6,D14.7)') 'PLAIN',1,'behav','TETEz','avrd=',0.0000000D+01
  write(128,'(A6)') '$nodty'
  write(128,'(1X,A5,5X,i2,16X,3(A5,D14.7,2X))') &
  'NO6xx',1,'coo1=',0.0000000D+00,'coo2=',0.3060000D+00,'coo3=',0.3173983D+00
  write(128,'(29X,3(A5,D14.7,2X))')'coo4=',0.D0,'coo5=',0.D0,'coo6=',0.D0
  write(128,'(A6)') '$tacty'
  write(128,'(1X,A5,5X,i2,2X,A5,2X,A5,2X,3(A5,D14.7,2X))') &
  'PLANx',1,'color','WALLx',&
  'axe1=',0.0800000D+01,'axe2=',0.0800000D+01,'axe3=',0.6000000D-02
  write(128,'(A6)') '$$$$$$'

  write(128,'(A6)') '$bdyty'
  write(128,'(1X,A5,1X,i6)') 'RBDY3',40004
  write(128,'(A6)') '$blmty'
  write(128,'(1X,A5,5X,i2,2X,A5,2X,A5,1X,A6,D14.7)') 'PLAIN',1,'behav','TETEz','avrd=',0.0000000D+01
  write(128,'(A6)') '$nodty'
  write(128,'(1X,A5,5X,i2,16X,3(A5,D14.7,2X))') &
  'NO6xx',1,'coo1=',0.0000000D+00,'coo2=',-0.3060000D+00,'coo3=',0.3173983D+00
  write(128,'(29X,3(A5,D14.7,2X))')'coo4=',0.D0,'coo5=',0.D0,'coo6=',0.D0
  write(128,'(A6)') '$tacty'
  write(128,'(1X,A5,5X,i2,2X,A5,2X,A5,2X,3(A5,D14.7,2X))') &
  'PLANx',1,'color','WALLx',&
  'axe1=',0.0800000D+01,'axe2=',0.0800000D+01,'axe3=',0.6000000D-02
  write(128,'(A6)') '$$$$$$'

  write(128,'(A6)') '$bdyty'
  write(128,'(1X,A5,1X,i6)') 'RBDY3',40005
  write(128,'(A6)') '$blmty'
  write(128,'(1X,A5,5X,i2,2X,A5,2X,A5,1X,A6,D14.7)') 'PLAIN',1,'behav','TETEz','avrd=',0.0000000D+01
  write(128,'(A6)') '$nodty'
  write(128,'(1X,A5,5X,i2,16X,3(A5,D14.7,2X))') &
  'NO6xx',1,'coo1=',0.3060000D+00,'coo2=',-0.0000000D+00,'coo3=',0.3173983D+00
  write(128,'(29X,3(A5,D14.7,2X))')'coo4=',0.D0,'coo5=',0.D0,'coo6=',0.D0
  write(128,'(A6)') '$tacty'
  write(128,'(1X,A5,5X,i2,2X,A5,2X,A5,2X,3(A5,D14.7,2X))') &
  'PLANx',1,'color','WALLx',&
  'axe1=',0.0800000D+01,'axe2=',0.0800000D+01,'axe3=',0.6000000D-02
  write(128,'(A6)') '$$$$$$'

  write(128,'(A6)') '$bdyty'
  write(128,'(1X,A5,1X,i6)') 'RBDY3',40006
  write(128,'(A6)') '$blmty'
  write(128,'(1X,A5,5X,i2,2X,A5,2X,A5,1X,A6,D14.7)') 'PLAIN',1,'behav','TETEz','avrd=',0.0000000D+01
  write(128,'(A6)') '$nodty'
  write(128,'(1X,A5,5X,i2,16X,3(A5,D14.7,2X))') &
  'NO6xx',1,'coo1=',-0.3060000D+00,'coo2=',-0.0000000D+00,'coo3=',0.3173983D+00
  write(128,'(29X,3(A5,D14.7,2X))')'coo4=',0.D0,'coo5=',0.D0,'coo6=',0.D0
  write(128,'(A6)') '$tacty'
  write(128,'(1X,A5,5X,i2,2X,A5,2X,A5,2X,3(A5,D14.7,2X))') &
  'PLANx',1,'color','WALLx',&
  'axe1=',0.0800000D+01,'axe2=',0.0800000D+01,'axe3=',0.6000000D-02
  write(128,'(A6)') '$$$$$$'


end subroutine write_plan


subroutine write_bodies(K_poly_ini_ini)
implicit none
integer                              :: i,k,l,j,cpt_rugu,K_poly_ini_ini
real(kind=8)                         :: K2,dmoyen
real(kind=8)                         ::  thetaX,thetaY,thetaZ
real(kind=8),dimension(3,3)          ::  localframe
real(kind=8),dimension(:),allocatable :: tab_a,tab_b,tab_g

print*,'.... write BODIES_POLYR.DAT ....', K_poly_ini_ini
if (K_poly_ini_ini==1) then
  open(unit=128,file='BODIES_POLYR.DAT',status='replace') 

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
end if

!do i=1,nb_polyedre
  write(128,'(A6)') '$bdyty'
  write(128,'(1X,A5,1X,i6)') &
  'RBDY3',K_poly_ini_ini
  write(128,'(A6)') &
  '$blmty'
  write(128,'(1X,A5,5X,i2,2X,A5,2X,A5,1X,A6,D14.7)') &
  'PLAIN',1,'behav','TDURx','avrd=',0.0000000D+01
  write(128,'(29X,3(A5,D14.7,2X))') &
  'I1  =',0.D0,'I2  =',0.D0,'I3  =',0.D0
  write(128,'(A6)') &
  '$nodty'
  write(128,'(1X,A5,5X,i2,16X,3(A5,D14.7,2X))') &
  'NO6xx',1,'coo1=',TAB_POLYEDRE_centre(K_poly_ini_ini)%centre(1),&
            'coo2=',TAB_POLYEDRE_centre(K_poly_ini_ini)%centre(2),&
            'coo3=',TAB_POLYEDRE_centre(K_poly_ini_ini)%centre(3)
  write(128,'(29X,3(A5,D14.7,2X))') &
  'coo4=',0.D0,'coo5=',0.D0,'coo6=',0.D0
  write(128,'(A6)') &
  '$tacty'
  write(128,'(1X,A5,5X,i2,2X,A5,2X,A5,2X,A10,i7,4X,A9,i7)') &
  'POLYR',1,'color','REDxx',&
  'nb_vertex=',TAB_POLYEDRE_new(K_poly_ini_ini)%nb_vertex,&
  'nb_faces=',TAB_POLYEDRE_new(K_poly_ini_ini)%nb_face
  
  if (initial_position) then
  !--------------on va calculer le centre de gravite pour pouvoir tous les translater en l'origine au cas ou
    TAB_POLYEDRE_new(K_poly_ini_ini)%centre(:) = 0.D0
    do k=1,TAB_POLYEDRE_new(K_poly_ini_ini)%nb_vertex
      TAB_POLYEDRE_new(K_poly_ini_ini)%centre(1) = &
      TAB_POLYEDRE_new(K_poly_ini_ini)%centre(1) + TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(k,1)
      
      TAB_POLYEDRE_new(K_poly_ini_ini)%centre(2) = & 
      TAB_POLYEDRE_new(K_poly_ini_ini)%centre(2) + TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(k,2)
      
      TAB_POLYEDRE_new(K_poly_ini_ini)%centre(3) = &
      TAB_POLYEDRE_new(K_poly_ini_ini)%centre(3) + TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(k,3)
    end do
    TAB_POLYEDRE_new(K_poly_ini_ini)%centre = &
    TAB_POLYEDRE_new(K_poly_ini_ini)%centre/real(TAB_POLYEDRE_new(K_poly_ini_ini)%nb_vertex,8)
   
   !--------------si ils sont pas centres on recalcule les coordonnees des sommets
    if (TAB_POLYEDRE_new(K_poly_ini_ini)%centre(1)/=0 .and. &
        TAB_POLYEDRE_new(K_poly_ini_ini)%centre(2)/=0 .and. &
        TAB_POLYEDRE_new(K_poly_ini_ini)%centre(3)/=0 ) then 
      do k=1,TAB_POLYEDRE_new(K_poly_ini_ini)%nb_vertex 
        TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(k,1) = & 
        TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(k,1)-TAB_POLYEDRE_new(K_poly_ini_ini)%centre(1)
       
        TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(k,2) = &
        TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(k,2)-TAB_POLYEDRE_new(K_poly_ini_ini)%centre(2)
       
        TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(k,3) = &
        TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(k,3)-TAB_POLYEDRE_new(K_poly_ini_ini)%centre(3)
      end do
    end if
    
    !------- On cherche le Rayon max -----
    TAB_POLYEDRE_new(K_poly_ini_ini)%Rmax = 0.D0
    do k=1,TAB_POLYEDRE_new(K_poly_ini_ini)%nb_vertex
      TAB_POLYEDRE_new(K_poly_ini_ini)%Rmax=max( TAB_POLYEDRE_new(K_poly_ini_ini)%Rmax , &
                                                 sqrt(TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(k,1)**2+&
                                                     TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(k,2)**2+&
                                                     TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(k,3)**2) )
    end do
  
    !------ On actualise la position de chaque sommet--------
    TAB_POLYEDRE_new(K_poly_ini_ini)%vertex      =   &
    (TAB_POLYEDRE_centre(K_poly_ini_ini)%Rmax/TAB_POLYEDRE_new(K_poly_ini_ini)%Rmax)*TAB_POLYEDRE_new(K_poly_ini_ini)%vertex
    
    
    do k=1,TAB_POLYEDRE_new(K_poly_ini_ini)%nb_vertex
      write(128,'(29X,3(A5,D14.7,2X))') &
      'coo1=',TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(k,1),&
      'coo2=',TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(k,2),&
      'coo3=',TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(k,3)
    end do
  else
    do k=1,TAB_POLYEDRE_new(K_poly_ini_ini)%nb_vertex
      write(128,'(29X,3(A5,D14.7,2X))') &
      'coo1=',TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(k,1),&
      'coo2=',TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(k,2),&
      'coo3=',TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(k,3)
    end do
  end if
  
  do k=1,TAB_POLYEDRE_new(K_poly_ini_ini)%nb_face
     write(128,'(29X,3(A5,i7,9X))') &
     'ver1=',TAB_POLYEDRE_new(K_poly_ini_ini)%face(k,1),&
     'ver2=',TAB_POLYEDRE_new(K_poly_ini_ini)%face(k,2),&
     'ver3=',TAB_POLYEDRE_new(K_poly_ini_ini)%face(k,3)
  end do
  write(128,'(A6)') '$$$$$$'
!end do

end subroutine write_bodies


!==============================================================================
! write BODIES AIRE 
!==============================================================================
subroutine write_bodies_aire(K_poly_ini_ini)
implicit none
integer                                     :: i,j,K_poly_ini_ini,cpt_ver,cpt_ang
integer                                     :: ver1,ver2,ver3
integer                                     :: ver1_bis,ver2_bis,ver3_bis
real(kind=8),dimension(3)                   :: a_i,b_i,c_i,a_j,b_j,c_j,g,r_vertex
real(kind=8),dimension(3)                   :: ab_i,ac_i,ab,ac,ab_j,ac_j,w_i,w_j,w,centre_p
real(kind=8)                                :: Atriangle,Apolyedre,k_aire,Asphere,d_i,d_j,GO_w_i,GO_w_j
real(kind=8)                                :: Ang,Ang_moyen,norm,k_ang,k_ang_moyen,k,p,var_k_ang
real(kind=8)                                :: xmax,xmin,ymax,ymin,zmax,zmin,Long,Haut,Larg
real(kind=8),dimension(3,3)                 :: fabric
real(kind=8)                                :: a_12,a_13,a_23,Long_temp,x,y,z,dmax_vertex,dmin_vertex,d_vertex
integer                                     :: ver1_Long,ver2_Long

if (K_poly_ini_ini==1) open(unit=145,file='BODIES_PARAMETERS.DAT',status='replace') 

!----------- Calcul de l aire -----------
Atriangle = 0
Apolyedre = 0
do j=1,TAB_POLYEDRE_new(K_poly_ini_ini)%nb_face
  ver1  =  TAB_POLYEDRE_new(K_poly_ini_ini)%face(j,1)
  ver2  =  TAB_POLYEDRE_new(K_poly_ini_ini)%face(j,2)
  ver3  =  TAB_POLYEDRE_new(K_poly_ini_ini)%face(j,3)

  ab(:) =  TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(ver1,:)-TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(ver2,:)
  ac(:) =  TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(ver1,:)-TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(ver3,:)

  w(1)=ab(2)*ac(3)-ab(3)*ac(2)
  w(2)=-ab(1)*ac(3)+ab(3)*ac(1)
  w(3)=ab(1)*ac(2)-ab(2)*ac(1)

  Atriangle = 0.5 * sqrt( w(1)**2 + w(2)**2 + w(3)**2 )
  Apolyedre = Apolyedre + Atriangle
end do
Asphere = 4 * 3.1415926 * 0.5 * 0.5
k_aire = Apolyedre / Asphere

print*,'      k_aire = ',k_aire

!----------- Calcul de l angularite ----------

g(:) = 0
do i=1,TAB_POLYEDRE_new(K_poly_ini_ini)%nb_vertex
  g(:) = g(:) + TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(i,:)
end do
g(:) = g(:) /real(TAB_POLYEDRE_new(K_poly_ini_ini)%nb_vertex,8)

Ang = 0
cpt_ang = 0
Ang_moyen = 0

do j=1,TAB_POLYEDRE_new(K_poly_ini_ini)%nb_face
  ver1  =  TAB_POLYEDRE_new(K_poly_ini_ini)%face(j,1)
  ver2  =  TAB_POLYEDRE_new(K_poly_ini_ini)%face(j,2)
  ver3  =  TAB_POLYEDRE_new(K_poly_ini_ini)%face(j,3)
  do i=1,TAB_POLYEDRE_new(K_poly_ini_ini)%nb_face
    cpt_ver = 0
    if (i==j) cycle
    ver1_bis  =  TAB_POLYEDRE_new(K_poly_ini_ini)%face(i,1)
    ver2_bis  =  TAB_POLYEDRE_new(K_poly_ini_ini)%face(i,2)
    ver3_bis  =  TAB_POLYEDRE_new(K_poly_ini_ini)%face(i,3)

    if ( (ver1==ver1_bis).or.(ver1==ver2_bis).or.(ver1==ver3_bis) ) cpt_ver = cpt_ver + 1
    if ( (ver2==ver1_bis).or.(ver2==ver2_bis).or.(ver2==ver3_bis) ) cpt_ver = cpt_ver + 1
    if ( (ver3==ver1_bis).or.(ver3==ver2_bis).or.(ver3==ver3_bis) ) cpt_ver = cpt_ver + 1
    if (cpt_ver == 2) then

      a_j(1)=TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(ver1,1)
      a_j(2)=TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(ver1,2)
      a_j(3)=TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(ver1,3)
      b_j(1)=TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(ver2,1)
      b_j(2)=TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(ver2,2)
      b_j(3)=TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(ver2,3)
      c_j(1)=TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(ver3,1)
      c_j(2)=TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(ver3,2)
      c_j(3)=TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(ver3,3)

      a_i(1)=TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(ver1_bis,1)
      a_i(2)=TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(ver1_bis,2)
      a_i(3)=TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(ver1_bis,3)
      b_i(1)=TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(ver2_bis,1)
      b_i(2)=TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(ver2_bis,2)
      b_i(3)=TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(ver2_bis,3)
      c_i(1)=TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(ver3_bis,1)
      c_i(2)=TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(ver3_bis,2)
      c_i(3)=TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(ver3_bis,3)
      
    !-------------Construction des vecteurs-------------------------------
      ab_i = b_i-a_i
      ac_i = c_i-a_i

      ab_j = b_j-a_j
      ac_j = c_j-a_j

    !-------------Calcul de la normale (produit vectoriel)-------------------------------
      w_i(1)=ab_i(2)*ac_i(3)-ab_i(3)*ac_i(2)
      w_i(2)=-ab_i(1)*ac_i(3)+ab_i(3)*ac_i(1)
      w_i(3)=ab_i(1)*ac_i(2)-ab_i(2)*ac_i(1)
      norm = sqrt(w_i(1)**2+w_i(2)**2+w_i(3)**2)
      w_i = w_i / norm

      w_j(1)=ab_j(2)*ac_j(3)-ab_j(3)*ac_j(2)
      w_j(2)=-ab_j(1)*ac_j(3)+ab_j(3)*ac_j(1)
      w_j(3)=ab_j(1)*ac_j(2)-ab_j(2)*ac_j(1)
      norm = sqrt(w_j(1)**2+w_j(2)**2+w_j(3)**2)
      w_j = w_j / norm

    !-------------Calcul du coefficient d-------------------------------
      d_i = -(w_i(1)*a_i(1)+w_i(2)*a_i(2)+w_i(3)*a_i(3))
      d_j = -(w_j(1)*a_j(1)+w_j(2)*a_j(2)+w_j(3)*a_j(3))

    !-------------Calcul du produit scalaire de GO.w, sachant que O appartient au plan --------------------
      GO_w_i = -( d_i + g(1)*w_i(1) + g(2)*w_i(2) + g(3)*w_i(3) )
      GO_w_j = -( d_j + g(1)*w_j(1) + g(2)*w_j(2) + g(3)*w_j(3) )

      if ( GO_w_i < 0 ) w_i = - w_i
      if ( GO_w_j < 0 ) w_j = - w_j
   
       Ang = abs( w_i(1)*w_j(1) + w_i(2)*w_j(2) + w_i(3)*w_j(3) ) / &
             ( sqrt(w_i(1)**2 + w_i(2)**2 + w_i(3)**2)*sqrt(w_j(1)**2 + w_j(2)**2 + w_j(3)**2) )
       
       Ang = Acos(Ang)
       Ang_moyen = Ang_moyen + Ang
       cpt_ang = cpt_ang + 1
    end if
  end do
end do
Ang = Ang_moyen /real(cpt_ang,8)
Total_ang_moyen = Total_ang_moyen + Ang / real(nb_polyedre,8)

k_ang = (Ang) / 3.1415926
print*,'       K_ang = ',K_ang

!----------- Calcul de la variance de l angularite ----------

g(:) = 0
do i=1,TAB_POLYEDRE_new(K_poly_ini_ini)%nb_vertex
  g(:) = g(:) + TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(i,:)
end do
g(:) = g(:) /real(TAB_POLYEDRE_new(K_poly_ini_ini)%nb_vertex,8)

Ang = 0
Ang_moyen = 0
var_k_ang = 0

do j=1,TAB_POLYEDRE_new(K_poly_ini_ini)%nb_face
  ver1  =  TAB_POLYEDRE_new(K_poly_ini_ini)%face(j,1)
  ver2  =  TAB_POLYEDRE_new(K_poly_ini_ini)%face(j,2)
  ver3  =  TAB_POLYEDRE_new(K_poly_ini_ini)%face(j,3)
  do i=1,TAB_POLYEDRE_new(K_poly_ini_ini)%nb_face
    cpt_ver = 0
    if (i==j) cycle
    ver1_bis  =  TAB_POLYEDRE_new(K_poly_ini_ini)%face(i,1)
    ver2_bis  =  TAB_POLYEDRE_new(K_poly_ini_ini)%face(i,2)
    ver3_bis  =  TAB_POLYEDRE_new(K_poly_ini_ini)%face(i,3)

    if ( (ver1==ver1_bis).or.(ver1==ver2_bis).or.(ver1==ver3_bis) ) cpt_ver = cpt_ver + 1
    if ( (ver2==ver1_bis).or.(ver2==ver2_bis).or.(ver2==ver3_bis) ) cpt_ver = cpt_ver + 1
    if ( (ver3==ver1_bis).or.(ver3==ver2_bis).or.(ver3==ver3_bis) ) cpt_ver = cpt_ver + 1
    if (cpt_ver == 2) then

      a_j(1)=TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(ver1,1)
      a_j(2)=TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(ver1,2)
      a_j(3)=TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(ver1,3)
      b_j(1)=TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(ver2,1)
      b_j(2)=TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(ver2,2)
      b_j(3)=TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(ver2,3)
      c_j(1)=TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(ver3,1)
      c_j(2)=TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(ver3,2)
      c_j(3)=TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(ver3,3)

      a_i(1)=TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(ver1_bis,1)
      a_i(2)=TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(ver1_bis,2)
      a_i(3)=TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(ver1_bis,3)
      b_i(1)=TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(ver2_bis,1)
      b_i(2)=TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(ver2_bis,2)
      b_i(3)=TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(ver2_bis,3)
      c_i(1)=TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(ver3_bis,1)
      c_i(2)=TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(ver3_bis,2)
      c_i(3)=TAB_POLYEDRE_new(K_poly_ini_ini)%vertex(ver3_bis,3)
      
    !-------------Construction des vecteurs-------------------------------
      ab_i = b_i-a_i
      ac_i = c_i-a_i

      ab_j = b_j-a_j
      ac_j = c_j-a_j

    !-------------Calcul de la normale (produit vectoriel)-------------------------------
      w_i(1)=ab_i(2)*ac_i(3)-ab_i(3)*ac_i(2)
      w_i(2)=-ab_i(1)*ac_i(3)+ab_i(3)*ac_i(1)
      w_i(3)=ab_i(1)*ac_i(2)-ab_i(2)*ac_i(1)
      norm = sqrt(w_i(1)**2+w_i(2)**2+w_i(3)**2)
      w_i = w_i / norm

      w_j(1)=ab_j(2)*ac_j(3)-ab_j(3)*ac_j(2)
      w_j(2)=-ab_j(1)*ac_j(3)+ab_j(3)*ac_j(1)
      w_j(3)=ab_j(1)*ac_j(2)-ab_j(2)*ac_j(1)
      norm = sqrt(w_j(1)**2+w_j(2)**2+w_j(3)**2)
      w_j = w_j / norm

    !-------------Calcul du coefficient d-------------------------------
      d_i = -(w_i(1)*a_i(1)+w_i(2)*a_i(2)+w_i(3)*a_i(3))
      d_j = -(w_j(1)*a_j(1)+w_j(2)*a_j(2)+w_j(3)*a_j(3))

    !-------------Calcul du produit scalaire de GO.w, sachant que O appartient au plan --------------------
      GO_w_i = -( d_i + g(1)*w_i(1) + g(2)*w_i(2) + g(3)*w_i(3) )
      GO_w_j = -( d_j + g(1)*w_j(1) + g(2)*w_j(2) + g(3)*w_j(3) )

      if ( GO_w_i < 0 ) w_i = - w_i
      if ( GO_w_j < 0 ) w_j = - w_j
   
       Ang = abs( w_i(1)*w_j(1) + w_i(2)*w_j(2) + w_i(3)*w_j(3) ) / &
             ( sqrt(w_i(1)**2 + w_i(2)**2 + w_i(3)**2)*sqrt(w_j(1)**2 + w_j(2)**2 + w_j(3)**2) )
       
       Ang = Acos(Ang) / 3.1415926

       var_k_ang = var_k_ang + ( ( Ang - k_ang )**2 ) / real(cpt_ang,8)

    end if
  end do
end do
Total_var_ang_moyen = Total_var_ang_moyen + var_k_ang / real(nb_polyedre,8)

print*,'  var(K_ang) = ',var_k_ang

!----- Calcul du rapport d'aspect

xmax = 0.D0
xmin = 1000000  
ymax = 0.D0
ymin = 1000000
Long = 0.D0
larg = 100000
Haut = 0.D0

do i=1,TAB_POLYEDRE_new(K_poly_ini_ini)%nb_vertex
  do j=1,TAB_POLYEDRE_new(K_poly_ini_ini)%nb_vertex
  if (i==j) cycle
  Long_temp = sqrt( (vertex_List(i,1)-vertex_List(j,1))**2 +&
                    (vertex_List(i,2)-vertex_List(j,2))**2 +&
                    (vertex_List(i,3)-vertex_List(j,3))**2 )
  if (Long_temp>Long) then
    Long = Long_temp
    ver1_long = i
    ver2_long = j
  end if 
  end do
end do

ac(:) = vertex_List(ver1_long,:) - vertex_List(ver2_long,:)
norm  = sqrt( ac(1)**2 + ac(2)**2 + ac(3)**2 )
ac(:) = ac(:) / norm
centre_p(:) = (vertex_List(ver1_long,:)+vertex_List(ver2_long,:))/2
p     = -( ac(1)*centre_p(1) +&
         ac(2)*centre_p(2) + &
         ac(3)*centre_p(3) )
         
do i=1,TAB_POLYEDRE_new(K_poly_ini_ini)%nb_vertex
  if ( (i==ver1_Long).or.(i==ver2_Long) ) cycle  
  k = -( p + ac(1)*vertex_List(i,1) + ac(2)*vertex_List(i,2) + ac(3)*vertex_List(i,3) )
  x = k*ac(1) + vertex_List(i,1)
  y = k*ac(2) + vertex_List(i,2)
  z = k*ac(3) + vertex_List(i,3)
  xmax = max( xmax, x )
  ymax = max( ymax, y )
  xmin = min( xmin, x )
  ymin = min( ymin, y )
end do
larg = min(abs(xmax-xmin),abs(ymax-ymin))
Haut = max(abs(xmax-xmin),abs(ymax-ymin))

Total_LongH_moyen = Total_LongH_moyen + Haut/Long
Total_LargH_moyen = Total_LargH_moyen + Larg/Long

print*,'         H/L = ',Haut/Long
print*,'         l/L = ',Larg/Long

!----- Calcul de l'anisotropie des particules
!fabric = 0
!do i=1,TAB_POLYEDRE_new(K_poly_ini_ini)%nb_vertex
!  r_vertex(:) = vertex_List(i,:) - g(:)
!  norm = sqrt( (r_vertex(1))**2 + (r_vertex(2))**2 + (r_vertex(3))**2 )
!  r_vertex(:) = r_vertex(:) / norm
!  Fabric(1,1:3) = Fabric(1,1:3) + r_vertex(1)*r_vertex(1:3)  / real(TAB_POLYEDRE_new(K_poly_ini_ini)%nb_vertex,8)
!  Fabric(2,1:3) = Fabric(2,1:3) + r_vertex(2)*r_vertex(1:3)  / real(TAB_POLYEDRE_new(K_poly_ini_ini)%nb_vertex,8)
!  Fabric(3,1:3) = Fabric(3,1:3) + r_vertex(3)*r_vertex(1:3)  / real(TAB_POLYEDRE_new(K_poly_ini_ini)%nb_vertex,8)
!end do
!a_13 = abs( 5*(Fabric(3,3)-Fabric(1,1))/2 )
!a_23 = abs( 5*(Fabric(3,3)-Fabric(2,2))/2 )
!a_12 = abs( 5*(Fabric(2,2)-Fabric(1,1))/2 ) 

!a_12_moyen = a_12_moyen + a_12
!a_13_moyen = a_13_moyen + a_13
!a_23_moyen = a_23_moyen + a_23

!print*,'         a12 = ',a_12
!print*,'         a13 = ',a_13
!print*,'         a23 = ',a_23

dmax_vertex = 0
dmin_vertex = 1000000
do i=1,TAB_POLYEDRE_new(K_poly_ini_ini)%nb_vertex
  d_vertex = sqrt( (g(1)-vertex_List(i,1))**2 + (g(2)-vertex_List(i,2))**2 + (g(3)-vertex_List(i,3))**2 )
  dmin_vertex = min( dmin_vertex,d_vertex )
  dmax_vertex = max( dmax_vertex,d_vertex )
end do
print*,'  dmax/dmin = ',dmax_vertex/dmin_vertex

print*,'------'
write(145,'(2X,i7,15(2X,D14.7))') K_poly_ini_ini, K_ang, var_k_ang, k_aire ,Long/Haut, Larg/Haut  ! , a_12,a_13,a_23

if (K_poly_ini_ini == nb_polyedre) then
  k_ang_moyen = (Total_ang_moyen) / 3.1415926

  print*,'K_ang_moyen      = ', k_ang_moyen
  print*,'Var(K_ang_moyen) = ',Total_var_ang_moyen
  
  Total_LongH_moyen = Total_LongH_moyen / real(nb_polyedre,8)
  Total_LargH_moyen = Total_LargH_moyen / real(nb_polyedre,8)
  print*,'(H/L)moyen       = ', Total_LongH_moyen
  print*,'(l/L)moyen       = ', Total_LargH_moyen
  
!  a_12_moyen = a_12_moyen / real(nb_polyedre,8)
!  a_13_moyen = a_13_moyen / real(nb_polyedre,8)
!  a_23_moyen = a_23_moyen / real(nb_polyedre,8)
!  print*,'(a_12)moyen      = ', a_12_moyen
!  print*,'(a_13)moyen      = ', a_13_moyen
!  print*,'(a_23)moyen      = ', a_23_moyen

  write(145,'(A10,D14.7)') 'K_ang   = ', k_ang_moyen
  write(145,'(A10,D14.7)') 'V(K_ang)= ', Total_var_ang_moyen
  write(145,'(A10,D14.7)') '(H/L)moy= ', Total_LongH_moyen
  write(145,'(A10,D14.7)') '(l/L)moy= ', Total_LargH_moyen
  write(145,'(A10,D14.7)') '(a12)moy= ', a_12_moyen
  write(145,'(A10,D14.7)') '(a13)moy= ', a_13_moyen
  write(145,'(A10,D14.7)') '(a23)moy= ', a_23_moyen

end if

end subroutine write_bodies_aire





end program Construct_polyr

