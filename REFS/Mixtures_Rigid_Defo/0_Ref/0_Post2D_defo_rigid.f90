program post2D

implicit none
 
type  ::  T_CORPS_rigid
    integer                                  ::  num,nb_vertex
    real(kind=8)                             ::  rayon,Aire
    real(kind=8),dimension(3)                ::  centre,centre_ref
    real(kind=8),dimension(3)                ::  V
    real(kind=8),allocatable,dimension(:,:)  ::  vertex,vertex_ref
    character(len=5)                         ::  behav
    character(len=6)                         ::  name
end type T_CORPS_rigid

type :: T_CORPS_deformable_mesh
  integer                                  ::  type_mesh
  integer,dimension(4)                     ::  connectivity
  character(len=5)                         ::  behav
end type T_CORPS_deformable_mesh

type  ::  T_CORPS_deformable
  integer                                                 ::  num
  real(kind=8),dimension(3)                               ::  centre
  real(kind=8)                                            ::  Aire
  real(kind=8),allocatable,dimension(:,:)                 ::  V
  character(len=6)                                        ::  name
  integer                                                 ::  nb_mesh,nb_node,nb_tacty
  real(kind=8),allocatable,dimension(:,:)                 ::  node,node_ref,tacty
  type(T_CORPS_deformable_mesh),allocatable,dimension(:)  ::  mesh_connectivities
end type T_CORPS_deformable

type  ::  T_CONTACT
    integer                                  ::  icdent,ianent,sect,iantact,ictacty
    integer                                  ::  type
    real(kind=8),dimension(3)                ::  n                    
    real(kind=8),dimension(3)                ::  t,F
    real(kind=8),dimension(3)                ::  coor_ctc                  
    real(kind=8)                             ::  rn,rt,rs             
    real(kind=8)                             ::  vln,vlt,vs,gap
    character(len=5)                         ::  nature
    character(len=5)                         ::  statut   
end type T_CONTACT

type  ::  T_DOF
    integer                                  ::  num
    real(kind=8),dimension(3)                ::  X                    
    real(kind=8),dimension(3)                ::  V                  
end type T_DOF
 
type(T_CORPS_rigid),allocatable,dimension(:)      ::  TAB_RIGID_PARTICLE
type(T_CORPS_deformable),allocatable,dimension(:) ::  TAB_MESHE_PARTICLE
type(T_CONTACT),allocatable,dimension(:)          ::  TAB_CONTACT,TAB_CONTACT_PARTICLE

!DEFINITIONS Vloc-Dof-Bodies
!================================================================
integer                                     ::  nb_paroi=0,nbParticules=0,nb_DISKx,nb_POLYG,nb_JONCx,nb_DEFOR,nb_RIGID
integer                                     ::  compteur_clout=0,all_pas=0
logical                                     ::  the_first_time=.true., fin_post=.false.
integer                                     ::  pas,nb_ligneCONTACT,nb_ligneCONTACT_PARTICLE
real(kind=8)                                ::  temps
!DEFINITIONS BOITE
!================================================================
real(kind=8)                                ::  hauteur,longueur,hauteur0,longueur0,epsilon1,epsilon2,epsilon_q
!COMMANDES D APPELLE
!================================================================
character(len=30)                           ::  command
integer                                     ::  calcul_compacity=0
integer                                     ::  calcul_ps_file=0,Trace_part_ps,Trace_contact_line_ps,Trace_contact_force_ps


!PARAMETRE
!================================================================
real(kind=8),parameter                     :: pi=3.141593D0

!================================================================
!Lecture du fichier de commande
!================================================================
print*,'-------------------------------------------'
print*,'!                                         !'
print*,'!      Module de Post-traitement 2D       !'
print*,'!             Emilien Azema               !'
print*,'!                                         !'
print*,'-------------------------------------------'
print*,'Lecture des commandes'
open(unit=1,file='BOITE.DAT',status='old')
do
  read(1,'(A30)') command
  print*,command
     
  if (command=='VLOC_DOF_BODIES              :') then
    read(1,*) nb_paroi
    read(1,*) nb_RIGID
    read(1,*) nb_DEFOR
    read(1,*) all_pas
    nbParticules = nb_DEFOR+nb_RIGID
    nb_DEFOR = 0
    nb_RIGID = 0
    cycle
  end if

!################# Over all time steps
  if (command=='COMPACITY                    :') then
    calcul_compacity=1
    open (unit=100,file='COMPACITY.DAT',status='replace')
    cycle
  end if

!################# Over One time step
  if (command=='WRITE ONE PS FILE            :') then
    read(1,*) Trace_part_ps
    read(1,*) Trace_contact_line_ps
    read(1,*) Trace_contact_force_ps
    calcul_ps_file=1
    open (unit=1008,file='PS_ONE_FILE.ps',status='replace')
    cycle
  end if

!################# End reading
  if (command=='END                          :') exit
end do
close(1)

!==============================================================================
!Appel des differente ressources lues dans le fichier de commande
!==============================================================================
do
  call read_Vloc_dof_bodies
  call compute_size_box

  if (all_pas<0) then   ! We wants to calculate each time step
    if (calcul_compacity            == 1)       call compacity
  else ! We wants to calculate a given time step
    if (calcul_ps_file              == 1)  call write_one_ps_file
    fin_post=.true.
  end if

  if (fin_post) exit
end do
!==============================================================================

 contains
!#################################
!#################################
! PART 1 : subroutines for all time steps
!#################################
!#################################

!==============================================================================
! Calcul de la compacite
!==============================================================================
subroutine compacity
implicit none
real(kind=8)                             :: Air_grain_rigid,Air_grain_defor,Air_totale,Air_boite
integer                                  :: i

Air_grain_rigid = 0.D0
Air_grain_defor = 0.D0
do i=1,nb_RIGID
  Air_grain_rigid = Air_grain_rigid + TAB_RIGID_PARTICLE(i)%Aire
end do
do i=1,nb_DEFOR
  Air_grain_defor = Air_grain_defor + TAB_MESHE_PARTICLE(i)%aire
end do

Air_totale = Air_grain_defor + Air_grain_rigid
Air_boite = longueur * hauteur

write(100,'(15(1X,D12.5))') temps,epsilon1,epsilon_q,Air_totale/Air_boite

end subroutine compacity



!#################################
!#################################
! PART 2 : subroutines for one time step
!#################################
!#################################
!==============================================================================
! Write one PS File
!==============================================================================

subroutine write_one_ps_file
implicit none
type T_xs_deformable
  real(kind=8),dimension(:),allocatable  ::  node
end type T_xs_deformable

integer                                          ::  i,cpt,j,k,nodea,nodeb,nodec,noded,an,cd
real(kind=8),dimension(:,:),allocatable          ::  x,xs,ys,x_ctc,x_def
type(T_xs_deformable),allocatable,dimension(:)   ::  xs_def,ys_def
real(kind=8),dimension(3)                        ::  centre
real(kind=8),dimension(:),allocatable            ::  xc1,yc1,xc2,yc2
real(kind=8)                                     ::  sca0,xdil,xsca,ysca,sca,xamp,yamp,xoff,yoff,&
                                                     rmin,rmax,rmean,rver,xmin,xmax,ymin,ymax
real(kind=8)                                      :: frel,Force_normal_max,Force_normal_min


print*,'-----> Write Ps File'

call psinit()
!---------------------------------------------------------------------
! Rescale... DO NOT MODIFY
!---------------------------------------------------------------------
sca0 = 2200.
xdil  = 1
xmin  = 0.
ymin  = 0.
xmax  = longueur
ymax  = hauteur

xsca=xdil*sca0/(xmax-xmin)
ysca=xdil*sca0/(ymax-ymin)
sca=min(xsca,ysca)
xamp=(xmax-xmin)*sca
yamp=(ymax-ymin)*sca
xoff=(sca0-xamp)*.5
yoff=(sca0-yamp)*.5

rmin=sca*rmin
rmax=sca*rmax
rmean=sca*rmean
rver=sca*rver

allocate(x(nb_RIGID+nb_paroi,20))
allocate(x_def(nb_DEFOR,20))
allocate(x_ctc(nb_ligneCONTACT_PARTICLE,20))
allocate(xs(nb_RIGID+nb_paroi,60))
allocate(ys(nb_RIGID+nb_paroi,60))
allocate(xs_def(nb_DEFOR))
allocate(ys_def(nb_DEFOR))
do i=1,nb_DEFOR
  allocate( xs_def(i)%node(TAB_MESHE_PARTICLE(i)%nb_node) )
  allocate( ys_def(i)%node(TAB_MESHE_PARTICLE(i)%nb_node) )
end do

do j=1,nb_RIGID+nb_paroi
  x(j,2)=(TAB_RIGID_PARTICLE(j)%centre(1)-xmin)*sca+xoff
  x(j,3)=(TAB_RIGID_PARTICLE(j)%centre(2)-ymin)*sca+yoff
  if ((TAB_RIGID_PARTICLE(j)%name == 'POLYG ').or.(TAB_RIGID_PARTICLE(j)%name == 'JONCx ')) then
    do k=1,TAB_RIGID_PARTICLE(j)%nb_vertex
      xs(j,k)=sca*(TAB_RIGID_PARTICLE(j)%vertex(1,k)-xmin)+xoff
      ys(j,k)=sca*(TAB_RIGID_PARTICLE(j)%vertex(2,k)-ymin)+yoff
    enddo
  end if
enddo

do j=1,nb_DEFOR
  x_def(j,2)=(TAB_MESHE_PARTICLE(j)%centre(1)-xmin)*sca+xoff
  x_def(j,3)=(TAB_MESHE_PARTICLE(j)%centre(2)-ymin)*sca+yoff
  do k=1,TAB_MESHE_PARTICLE(j)%nb_node
    xs_def(j)%node(k) = sca*(TAB_MESHE_PARTICLE(j)%node(k,1)-xmin)+xoff
    ys_def(j)%node(k) = sca*(TAB_MESHE_PARTICLE(j)%node(k,2)-ymin)+yoff
  end do
end do

do j=1,nb_ligneCONTACT_PARTICLE
  x_ctc(j,2)=(TAB_CONTACT_PARTICLE(j)%coor_ctc(1)-xmin)*sca+xoff
  x_ctc(j,3)=(TAB_CONTACT_PARTICLE(j)%coor_ctc(2)-ymin)*sca+yoff
end do
!---------------------------------------------------------------------
! End Rescale...
!---------------------------------------------------------------------

!---------------------------------------------------------------------
! draws particules
!---------------------------------------------------------------------
if (Trace_part_ps==1) then
  write(1008,*) '1 setlinewidth'
  do i=1,nb_RIGID+nb_paroi
    if (TAB_RIGID_PARTICLE(i)%name == 'POLYG ') then
      write(1008,*) '',(xs(i,k),ys(i,k),k=1,TAB_RIGID_PARTICLE(i)%nb_vertex),TAB_RIGID_PARTICLE(i)%nb_vertex,0.3,0.3,0.3, 'polyc'
    end if
    if (TAB_RIGID_PARTICLE(i)%name == 'DISKx ') then
      write(1008,*) x(i,2),x(i,3),TAB_RIGID_PARTICLE(i)%rayon*sca,0.3,0.3,0.3,'C2'
    end if
    if (TAB_RIGID_PARTICLE(i)%name == 'JONCx ') then
      write(1008,*) '',(xs(i,k),ys(i,k),k=1,4),4,0.8,0.8,0.8, 'polyc'
    end if
  end do
  do i=1,nb_DEFOR
    do k=1,TAB_MESHE_PARTICLE(i)%nb_mesh
      if (TAB_MESHE_PARTICLE(i)%mesh_connectivities(k)%type_mesh==4) then
        nodea = TAB_MESHE_PARTICLE(i)%mesh_connectivities(k)%connectivity(1)
        nodeb = TAB_MESHE_PARTICLE(i)%mesh_connectivities(k)%connectivity(2)
        nodec = TAB_MESHE_PARTICLE(i)%mesh_connectivities(k)%connectivity(3)
        noded = TAB_MESHE_PARTICLE(i)%mesh_connectivities(k)%connectivity(4)
        write(1008,*) '',xs_def(i)%node(nodea),ys_def(i)%node(nodea),&
                         xs_def(i)%node(nodeb),ys_def(i)%node(nodeb),&
                         xs_def(i)%node(nodec),ys_def(i)%node(nodec),&
                         xs_def(i)%node(noded),ys_def(i)%node(noded) ,4,0.3,0.6,0.8, 'polyc'
      end if
      if (TAB_MESHE_PARTICLE(i)%mesh_connectivities(k)%type_mesh==3) then
        nodea = TAB_MESHE_PARTICLE(i)%mesh_connectivities(k)%connectivity(1)
        nodeb = TAB_MESHE_PARTICLE(i)%mesh_connectivities(k)%connectivity(2)
        nodec = TAB_MESHE_PARTICLE(i)%mesh_connectivities(k)%connectivity(3)
        write(1008,*) '',xs_def(i)%node(nodea),ys_def(i)%node(nodea),&
                         xs_def(i)%node(nodeb),ys_def(i)%node(nodeb),&
                         xs_def(i)%node(nodec),ys_def(i)%node(nodec),&
                         3,0.3,0.6,0.8, 'polyc'
      end if
    end do
  end do
end if
!---------------------------------------------------------------------
! draws contact lines
!---------------------------------------------------------------------
if (Trace_contact_line_ps==1) then
  write(1008,*) 1,'setlinecap'
  do i=1,nb_ligneCONTACT_PARTICLE
    if ( TAB_CONTACT_PARTICLE(i)%rn     == 0. ) cycle
    cd = TAB_CONTACT_PARTICLE(i)%icdent
    an = TAB_CONTACT_PARTICLE(i)%ianent
    if ( (cd>nbParticules).or.(an>nbParticules) ) cycle
    frel = 0.025
    if (TAB_CONTACT_PARTICLE(i)%nature=='DKDKx') then
      write(1008,*) x(cd,2),x(cd,3),x_ctc(i,2),x_ctc(i,3),frel,0.8,0,0,' Li'
      write(1008,*) x(an,2),x(an,3),x_ctc(i,2),x_ctc(i,3),frel,0.8,0,0,' Li'
    end if
    if (TAB_CONTACT_PARTICLE(i)%nature=='DKALp') then
      write(1008,*) x(cd,2),x(cd,3),x_ctc(i,2),x_ctc(i,3),frel,0.8,0,0,' Li'
      write(1008,*) x_def(an,2),x_def(an,3),x_ctc(i,2),x_ctc(i,3),frel,0.8,0,0,' Li'
    end if
    if (TAB_CONTACT_PARTICLE(i)%nature=='CLALp') then
      write(1008,*) x_def(cd,2),x_def(cd,3),x_ctc(i,2),x_ctc(i,3),frel,0.8,0,0,' Li'
      write(1008,*) x_def(an,2),x_def(an,3),x_ctc(i,2),x_ctc(i,3),frel,0.8,0,0,' Li'
    end if
    if (TAB_CONTACT_PARTICLE(i)%nature=='DKJCx') then
      write(1008,*) x(cd,2),x(cd,3),x_ctc(i,2),x_ctc(i,3),frel,0.8,0,0,' Li'
    end if
    if (TAB_CONTACT_PARTICLE(i)%nature=='CLJCx') then
      write(1008,*) x_def(cd,2),x_def(cd,3),x_ctc(i,2),x_ctc(i,3),frel,0.8,0,0,' Li'
    end if
  end do
end if
!---------------------------------------------------------------------
! draws contact forces
!---------------------------------------------------------------------
if (Trace_contact_force_ps==1) then
  Force_normal_max = 0.
  Force_normal_min = 10000000000.
  do i=1,nb_ligneCONTACT_PARTICLE
    if ( TAB_CONTACT_PARTICLE(i)%rn     == 0. ) cycle
    Force_normal_max = max(Force_normal_max,TAB_CONTACT_PARTICLE(i)%rn)
    Force_normal_min = min(Force_normal_min,TAB_CONTACT_PARTICLE(i)%rn)
  end do

  write(1008,*) 1,'setlinecap'
  do i=1,nb_ligneCONTACT_PARTICLE
    if ( TAB_CONTACT_PARTICLE(i)%rn     == 0. ) cycle
    cd = TAB_CONTACT_PARTICLE(i)%icdent
    an = TAB_CONTACT_PARTICLE(i)%ianent
    if ( (cd>nbParticules).or.(an>nbParticules) ) cycle
    frel = 0.8*abs( (TAB_CONTACT_PARTICLE(i)%rn-Force_normal_min)/(Force_normal_max+Force_normal_min)    )
    if (TAB_CONTACT_PARTICLE(i)%nature=='DKDKx') then
      write(1008,*) x(cd,2),x(cd,3),x_ctc(i,2),x_ctc(i,3),frel,0.8,0,0,' Li'
      write(1008,*) x(an,2),x(an,3),x_ctc(i,2),x_ctc(i,3),frel,0.8,0,0,' Li'
    end if
    if (TAB_CONTACT_PARTICLE(i)%nature=='DKALp') then
      write(1008,*) x(cd,2),x(cd,3),x_ctc(i,2),x_ctc(i,3),frel,0.8,0,0,' Li'
      write(1008,*) x_def(an,2),x_def(an,3),x_ctc(i,2),x_ctc(i,3),frel,0.8,0,0,' Li'
    end if
    if (TAB_CONTACT_PARTICLE(i)%nature=='CLALp') then
      write(1008,*) x_def(cd,2),x_def(cd,3),x_ctc(i,2),x_ctc(i,3),frel,0.8,0,0,' Li'
      write(1008,*) x_def(an,2),x_def(an,3),x_ctc(i,2),x_ctc(i,3),frel,0.8,0,0,' Li'
    end if
    if (TAB_CONTACT_PARTICLE(i)%nature=='DKJCx') then
      write(1008,*) x(cd,2),x(cd,3),x_ctc(i,2),x_ctc(i,3),frel,0.8,0,0,' Li'
    end if
    if (TAB_CONTACT_PARTICLE(i)%nature=='CLJCx') then
      write(1008,*) x_def(cd,2),x_def(cd,3),x_ctc(i,2),x_ctc(i,3),frel,0.8,0,0,' Li'
    end if
  end do
end if

!---------------------------------------------------------------------
! Fermer le fichier PostScript
!---------------------------------------------------------------------
write(1008,*) 'showpage'
close(1008)
end subroutine write_one_ps_file


!==============================================================================
!==============================================================================
!==============================================================================
!==============================================================================
!==============================================================================
!==============================================================================
!==============================================================================


!#################################
! PART 3 : Reading BODIES.DAT, DOF.OUT.XXX, Vloc-Rloc.OUT.XXX
!#################################
subroutine read_Vloc_dof_bodies
implicit none
integer                       ::  i,j,k,kk,err,num_part,nb_vertex,nodea,nodeb,nodec,noded
character(len=22)             ::  clout_DOF
character(len=28)             ::  clout_Vloc
character(len=6)              ::  text,name
character(len=13)             ::  text2
character(len=5)              ::  behav,statut
real(kind=8)                  ::  centre1,centre2,centre3
real(kind=8)                  ::  vertex_coor1,vertex_coor2,vertex_coor3
real(kind=8),dimension(3)     ::  g_gravity
real(kind=8)                  ::  coor1,coor2,coor3,Norm,air_maille
real(kind=8),dimension(2)     ::  X1,X2
real(kind=8)                  ::  Vx,Vy,Vr,gap
real(kind=8)                  ::  Rayon,ax1,ax2
integer                       ::  nb_plplx,nb_pljcx,nb_dkdkx,nb_dkjcx,nb_dkplx,nb_clalp,nb_cljcx,nb_dkalp
integer                       ::  N_PLPLX,N_PLJCx,N_DKDKx,N_DKJCx,N_DKPLx,N_CLALp,N_CLJCx,N_DKALp
integer                       ::  icdent,ianent,icdtact,iantact,cpt_type,cd,an,cpt1,cpt2,cpt3,cpt4
real(kind=8)                  ::  n1,n2,n3
real(kind=8)                  ::  rn,rt,rs,rn_temp,rt_temp,rs_temp
real(kind=8)                  ::  gap_temp
real(kind=8)                  ::  vln,vlt,vls,vln_temp,vlt_temp
real(kind=8)                  ::  coor_ctc1,coor_ctc2,coor_ctc3
real(kind=8)                  ::  coor_ctc1_temp,coor_ctc2_temp,coor_ctc3_temp
integer                       ::  n_simple,n_double
logical                       ::  first_double=.true.


if (all_pas<0) then
  compteur_clout = compteur_clout+1
  clout_DOF  =  '../OUTBOX/DOF.OUT.   '
  clout_Vloc =  '../OUTBOX/Vloc_Rloc.OUT.   '
else
  compteur_clout = all_pas
  clout_DOF  =  '../OUTBOX/DOF.OUT.   '
  clout_Vloc =  '../OUTBOX/Vloc_Rloc.OUT.   '
end if
if (compteur_clout<10) then
  WRITE(clout_DOF(19:20),'(I1)')   compteur_clout
  WRITE(clout_Vloc(25:26),'(I1)')  compteur_clout
else if ( (compteur_clout>=10) .and. (compteur_clout<100) ) then
  WRITE(clout_DOF(19:21),'(I2)')   compteur_clout
  WRITE(clout_Vloc(25:27),'(I2)')  compteur_clout
else if ( (compteur_clout>=100) .and. (compteur_clout<1000) ) then
  WRITE(clout_DOF(19:22),'(I3)')   compteur_clout
  WRITE(clout_Vloc(25:28),'(I3)')  compteur_clout
end if

if (the_first_time) then
  nb_DEFOR = 0
  nb_RIGID = 0
  open(unit=2,file='../DATBOX/BODIES.DAT',status='old')
  i = 0
  do
    read(2,'(A6)') text
    if (text == '$bdyty') then
      i = i + 1
      read(2,'(A6)') text
      if (text == ' MAILx') then
        nb_DEFOR = nb_DEFOR + 1
      else if (text == ' RBDY2') then
        nb_RIGID = nb_RIGID + 1
      end if
    end if
    if (i==nbParticules+nb_paroi) exit
  end do
  close(2)
  nb_RIGID = nb_RIGID - nb_paroi
  print*, "Nb rigide =",nb_RIGID
  print*, "Nb deformable =",nb_DEFOR

! ==== LECTURE DES CORPS MAILLES
  allocate(TAB_MESHE_PARTICLE(nb_DEFOR))
  TAB_MESHE_PARTICLE(:)%nb_node  = 0
  TAB_MESHE_PARTICLE(:)%nb_mesh  = 0
  TAB_MESHE_PARTICLE(:)%nb_tacty = 0
  open(unit=2,file='../DATBOX/BODIES.DAT',status='old')
  i = 0
  do
    read(2,'(A6)') text
    if (text == '$bdyty') then
      read(2,'(A6,1X,i6)') text,num_part
      if (text == ' MAILx') then
        i = i + 1
        TAB_MESHE_PARTICLE(num_part)%num = num_part
        read(2,'(A6)') text
        if (text == '$blmty') then
          do
            read(2,'(A6)') text
            if ((text==' Q4xxx').or.(text==' T3xxx')) TAB_MESHE_PARTICLE(num_part)%nb_mesh = &
                                                      TAB_MESHE_PARTICLE(num_part)%nb_mesh + 1
            if (text=='$nodty') exit
          end do
          do
            read(2,'(A6)') text
            if (text==' NO2xx') TAB_MESHE_PARTICLE(num_part)%nb_node = &
                                TAB_MESHE_PARTICLE(num_part)%nb_node + 1
            if (text=='$tacty') exit
          end do
          do
            read(2,'(A6)') text
            if ((text==' ALpxx').or.(text=='+ALpxx')) TAB_MESHE_PARTICLE(num_part)%nb_tacty = &
                                                      TAB_MESHE_PARTICLE(num_part)%nb_tacty + 1
            if (text=='$$$$$$') exit
          end do
        end if
      end if
    end if
    if (i==nb_DEFOR) exit
  end do
  close(2)

  do i=1,nb_DEFOR
    allocate(TAB_MESHE_PARTICLE(i)%mesh_connectivities(TAB_MESHE_PARTICLE(i)%nb_mesh))
    allocate(TAB_MESHE_PARTICLE(i)%tacty(TAB_MESHE_PARTICLE(i)%nb_tacty,2))
    allocate(TAB_MESHE_PARTICLE(i)%V(TAB_MESHE_PARTICLE(i)%nb_node,2))
    allocate(TAB_MESHE_PARTICLE(i)%node(TAB_MESHE_PARTICLE(i)%nb_node,2))
    allocate(TAB_MESHE_PARTICLE(i)%node_ref(TAB_MESHE_PARTICLE(i)%nb_node,2))
  end do

  open(unit=2,file='../DATBOX/BODIES.DAT',status='old')
  i = 0
  do
    read(2,'(A6)') text
    if (text == '$bdyty') then
      read(2,'(A6,1X,i6)') text,num_part
      if (text == ' MAILx') then
        i = i + 1
        read(2,'(A6)') text
        if (text == '$blmty') then
          do j=1,TAB_MESHE_PARTICLE(num_part)%nb_mesh
            read(2,'(A6)') text
            if (text==' Q4xxx') TAB_MESHE_PARTICLE(num_part)%mesh_connectivities(j)%type_mesh = 4
            if (text==' T3xxx') TAB_MESHE_PARTICLE(num_part)%mesh_connectivities(j)%type_mesh = 3
            read(2,*)
          end do
          read(2,*)
          do j=1,TAB_MESHE_PARTICLE(num_part)%nb_node
            read(2,'(34X,D14.7,7X,D14.7)') vertex_coor1,vertex_coor2
            TAB_MESHE_PARTICLE(num_part)%node_ref(j,1) = vertex_coor1
            TAB_MESHE_PARTICLE(num_part)%node_ref(j,2) = vertex_coor2
          end do
          read(2,*)
          do j=1,TAB_MESHE_PARTICLE(num_part)%nb_tacty
            read(2,'(A6,28X,i5,7X,i5)') text,nodea,nodeb
            if ((text==' CLxxx').or.(text=='+CLxxx')) cycle
            TAB_MESHE_PARTICLE(num_part)%tacty(j,1) = nodea
            TAB_MESHE_PARTICLE(num_part)%tacty(j,2) = nodeb
          end do
        end if
      end if
    end if
    if (i == nb_DEFOR) exit
  end do
  close(2)

  open(unit=2,file='../DATBOX/BODIES.DAT',status='old')
  i = 0
  do
    read(2,'(A6)') text
    if (text == '$bdyty') then
      read(2,'(A6,1X,i6)') text,num_part
      if (text == ' MAILx') then
        i = i +1
        read(2,'(A6)') text
        if (text == '$blmty') then
          do j=1,TAB_MESHE_PARTICLE(num_part)%nb_mesh
            TAB_MESHE_PARTICLE(num_part)%mesh_connectivities(j)%connectivity(:)=0
            if (TAB_MESHE_PARTICLE(num_part)%mesh_connectivities(j)%type_mesh==4) then
              read(2,'(20X,4(i7))') nodea,nodeb,nodec,noded
              read(2,'(A36)') behav
              TAB_MESHE_PARTICLE(num_part)%mesh_connectivities(j)%connectivity(1)=nodea
              TAB_MESHE_PARTICLE(num_part)%mesh_connectivities(j)%connectivity(2)=nodeb
              TAB_MESHE_PARTICLE(num_part)%mesh_connectivities(j)%connectivity(3)=nodec
              TAB_MESHE_PARTICLE(num_part)%mesh_connectivities(j)%connectivity(4)=noded
              TAB_MESHE_PARTICLE(num_part)%mesh_connectivities(j)%behav=behav
            endif
            if (TAB_MESHE_PARTICLE(num_part)%mesh_connectivities(j)%type_mesh==3) then
              read(2,'(20X,3(i7))') nodea,nodeb,nodec
              read(2,'(A36)') behav
              TAB_MESHE_PARTICLE(num_part)%mesh_connectivities(j)%connectivity(1)=nodea
              TAB_MESHE_PARTICLE(num_part)%mesh_connectivities(j)%connectivity(2)=nodeb
              TAB_MESHE_PARTICLE(num_part)%mesh_connectivities(j)%connectivity(3)=nodec
              TAB_MESHE_PARTICLE(num_part)%mesh_connectivities(j)%behav=behav
            endif
          end do
        end if
      end if
    end if
    if (i == nb_DEFOR) exit
  end do
  close(2)

! ==== LECTURE DES CORPS RIGIDES

  allocate(TAB_RIGID_PARTICLE(nb_RIGID+nb_paroi))
  nb_DISKx = 0
  nb_POLYG = 0
  open(unit=2,file='../DATBOX/BODIES.DAT',status='old')
  i = 0
  do
    read(2,'(A6)') text
    if (text == '$bdyty') then
      read(2,'(A6)') text
      if (text == ' MAILx') cycle
      i = i + 1
      read(2,*)
      read(2,*)
      read(2,*)
      read(2,*)
      read(2,*)
      read(2,'(A6)') name
      TAB_RIGID_PARTICLE(i)%name = name   !DISKx or POLYG or JONCx?
    end if
    if (i==nb_RIGID+nb_paroi) exit
  end do
  close(2)

  open(unit=2,file='../DATBOX/BODIES.DAT',status='old')
  i=0
  do
    read(2,'(A6)') text
    if (text == '$bdyty') then
      read(2,'(A6)') text
      if (text == ' MAILx') cycle
      i = i +1
      TAB_RIGID_PARTICLE(i)%num=i
      read(2,*)
      read(2,'(22X,A5)') behav
      TAB_RIGID_PARTICLE(i)%behav = behav
      read(2,*)
      read(2,'(29X, 3(5x,D14.7,2X))') centre1,centre2,centre3
      TAB_RIGID_PARTICLE(i)%centre_ref(1)=centre1
      TAB_RIGID_PARTICLE(i)%centre_ref(2)=centre2
      TAB_RIGID_PARTICLE(i)%centre_ref(3)=centre3
      TAB_RIGID_PARTICLE(i)%centre(:)=0.D0
      read(2,*)
      if (TAB_RIGID_PARTICLE(i)%name == ' POLYG') then
        nb_POLYG = nb_POLYG + 1
        read(2,'(40X,i6,11X,D14.7)') nb_vertex
        TAB_RIGID_PARTICLE(i)%nb_vertex=nb_vertex
        allocate(TAB_RIGID_PARTICLE(i)%vertex_ref(3,nb_vertex))
        allocate(TAB_RIGID_PARTICLE(i)%vertex(3,nb_vertex))
        g_gravity(:) = 0.D0
        do j=1,nb_vertex
          read(2,'(29X, 2(5x,D14.7,2X))') vertex_coor1,vertex_coor2
          TAB_RIGID_PARTICLE(i)%vertex_ref(1,j) = vertex_coor1
          TAB_RIGID_PARTICLE(i)%vertex_ref(2,j) = vertex_coor2
          TAB_RIGID_PARTICLE(i)%vertex_ref(3,j) = 0
          g_gravity(1) = g_gravity(1) + vertex_coor1
          g_gravity(2) = g_gravity(2) + vertex_coor2
          g_gravity(3) = g_gravity(3) + 0
        end do
        g_gravity(:) = g_gravity(:) / real(nb_vertex,8)
        TAB_RIGID_PARTICLE(i)%Rayon = 0.
        do j=1,nb_vertex
          TAB_RIGID_PARTICLE(i)%Rayon = max(TAB_RIGID_PARTICLE(i)%Rayon, &
                                   sqrt( (TAB_RIGID_PARTICLE(i)%vertex_ref(1,j)-g_gravity(1))**2+&
                                         (TAB_RIGID_PARTICLE(i)%vertex_ref(2,j)-g_gravity(2))**2+&
                                         (TAB_RIGID_PARTICLE(i)%vertex_ref(3,j)-g_gravity(3))**2))
        end do
        TAB_RIGID_PARTICLE(i)%Aire = 0.D0
        do j=1,nb_vertex-2
          TAB_RIGID_PARTICLE(i)%Aire = TAB_RIGID_PARTICLE(i)%Aire + &
               0.5d0*( ( ( TAB_RIGID_PARTICLE(i)%vertex_ref(1,j+1) - TAB_RIGID_PARTICLE(i)%vertex_ref(1,1))  &
                        *( TAB_RIGID_PARTICLE(i)%vertex_ref(2,j+2) - TAB_RIGID_PARTICLE(i)%vertex_ref(2,1))) &
                        -( ( TAB_RIGID_PARTICLE(i)%vertex_ref(2,j+1) - TAB_RIGID_PARTICLE(i)%vertex_ref(2,1))  &
                        *( TAB_RIGID_PARTICLE(i)%vertex_ref(1,j+2) - TAB_RIGID_PARTICLE(i)%vertex_ref(1,1))))
        enddo
        TAB_RIGID_PARTICLE(i)%name = 'POLYG '

      else if (TAB_RIGID_PARTICLE(i)%name == ' DISKx') then
          nb_DISKx = nb_DISKx + 1
          read(2,'(34X,D14.7)') Rayon
          TAB_RIGID_PARTICLE(i)%name  = 'DISKx '
          TAB_RIGID_PARTICLE(i)%Rayon = Rayon
          TAB_RIGID_PARTICLE(i)%Aire  = pi*Rayon**2
          ! Il faut mettre toute les autres valeurs a 0
          TAB_RIGID_PARTICLE(i)%nb_vertex = 0
          allocate(TAB_RIGID_PARTICLE(i)%vertex_ref(3,3))
          allocate(TAB_RIGID_PARTICLE(i)%vertex(3,3))
          TAB_RIGID_PARTICLE(i)%vertex_ref(:,:) = 0.D0
          TAB_RIGID_PARTICLE(i)%vertex(:,:) = 0.D0

      else if (TAB_RIGID_PARTICLE(i)%name == ' JONCx') then
          nb_JONCx = nb_JONCx + 1
          TAB_RIGID_PARTICLE(i)%name = 'JONCx '
          read(2,'(29X, 2(5x,D14.7,2X))') vertex_coor1,vertex_coor2
          TAB_RIGID_PARTICLE(i)%nb_vertex=4
          allocate(TAB_RIGID_PARTICLE(i)%vertex_ref(3,4))
          allocate(TAB_RIGID_PARTICLE(i)%vertex(3,4))
          TAB_RIGID_PARTICLE(i)%vertex_ref(1,1) =  vertex_coor1
          TAB_RIGID_PARTICLE(i)%vertex_ref(2,1) =  vertex_coor2
          TAB_RIGID_PARTICLE(i)%vertex_ref(3,1) = 0.

          TAB_RIGID_PARTICLE(i)%vertex_ref(1,2) =  - vertex_coor1
          TAB_RIGID_PARTICLE(i)%vertex_ref(2,2) =  vertex_coor2
          TAB_RIGID_PARTICLE(i)%vertex_ref(3,2) = 0.

          TAB_RIGID_PARTICLE(i)%vertex_ref(1,3) =  - vertex_coor1
          TAB_RIGID_PARTICLE(i)%vertex_ref(2,3) =  - vertex_coor2
          TAB_RIGID_PARTICLE(i)%vertex_ref(3,3) = 0.

          TAB_RIGID_PARTICLE(i)%vertex_ref(1,4) =  + vertex_coor1
          TAB_RIGID_PARTICLE(i)%vertex_ref(2,4) =  - vertex_coor2
          TAB_RIGID_PARTICLE(i)%vertex_ref(3,4) = 0.
      end if
    end if
    if (i==nb_RIGID+nb_paroi) exit
  end do
  nbParticules = nb_DEFOR + nb_DISKx + nb_POLYG
  close(2)
end if

open(unit=5,file=clout_DOF,iostat=err,status='old')
if (err/=0) then
  fin_post=.true.
  print*,'--> Fin du postraitement <--'
else
  fin_post=.false.
  i = 0
  k = 0
  kk= 0
  read(5,*)
  read(5,*)
  read(5,*)
  read(5,'(7X,i8,19X,D14.7)') pas,temps
  print*,'pas de temps n:',pas, temps
  print*,'--->',clout_DOF
  do
    read(5,'(A6)') text
    if (text == '$bdyty') then
      i = i+1
      read(5,'(A6)') text
      if (text == ' MAILx') then
        kk = kk + 1
        read(5,*)
        read(5,*)
        read(5,*)
        ! CALCUL DU CENTRE DE GRAVIY
        g_gravity(:) = 0.D0
        do j=1,TAB_MESHE_PARTICLE(kk)%nb_node
          read(5,'(29X, 2(5x,D14.7,2X))') coor1,coor2
          read(5,'(29X, 2(5x,D14.7,2X))') Vx,Vy
          TAB_MESHE_PARTICLE(kk)%node(j,1) = TAB_MESHE_PARTICLE(kk)%node_ref(j,1) + coor1
          TAB_MESHE_PARTICLE(kk)%node(j,2) = TAB_MESHE_PARTICLE(kk)%node_ref(j,2) + coor2
          TAB_MESHE_PARTICLE(kk)%V(j,1)    = Vx
          TAB_MESHE_PARTICLE(kk)%V(j,2)    = Vy

          g_gravity(1) = g_gravity(1) + TAB_MESHE_PARTICLE(kk)%node(j,1)
          g_gravity(2) = g_gravity(2) + TAB_MESHE_PARTICLE(kk)%node(j,2)
          g_gravity(3) = g_gravity(3) + 0
        end do
        g_gravity(:) = g_gravity(:) / real(TAB_MESHE_PARTICLE(kk)%nb_node,8)
        TAB_MESHE_PARTICLE(kk)%centre(:) = g_gravity(:)
        ! CALCUL DE L AIRE
        TAB_MESHE_PARTICLE(i)%aire = 0.
        do j=1,TAB_MESHE_PARTICLE(kk)%nb_mesh
          air_maille = 0.
          if (TAB_MESHE_PARTICLE(i)%mesh_connectivities(j)%type_mesh==4) then
            nodea = TAB_MESHE_PARTICLE(kk)%mesh_connectivities(j)%connectivity(1)
            nodeb = TAB_MESHE_PARTICLE(kk)%mesh_connectivities(j)%connectivity(2)
            nodec = TAB_MESHE_PARTICLE(kk)%mesh_connectivities(j)%connectivity(3)
            noded = TAB_MESHE_PARTICLE(kk)%mesh_connectivities(j)%connectivity(4)

            X1(:) = TAB_MESHE_PARTICLE(kk)%node(nodea,:)-TAB_MESHE_PARTICLE(kk)%node(nodeb,:)
            X2(:) = TAB_MESHE_PARTICLE(kk)%node(nodea,:)-TAB_MESHE_PARTICLE(kk)%node(nodec,:)

            air_maille = 0.5*abs( X1(1)*X2(2)-X1(2)*X2(1) )

            X1(:) = TAB_MESHE_PARTICLE(kk)%node(nodea,:)-TAB_MESHE_PARTICLE(kk)%node(nodec,:)
            X2(:) = TAB_MESHE_PARTICLE(kk)%node(nodea,:)-TAB_MESHE_PARTICLE(kk)%node(noded,:)

            air_maille = air_maille +0.5*abs( X1(1)*X2(2)-X1(2)*X2(1) )
          end if
          if (TAB_MESHE_PARTICLE(kk)%mesh_connectivities(j)%type_mesh==3) then
            nodea = TAB_MESHE_PARTICLE(kk)%mesh_connectivities(j)%connectivity(1)
            nodeb = TAB_MESHE_PARTICLE(kk)%mesh_connectivities(j)%connectivity(2)
            nodec = TAB_MESHE_PARTICLE(kk)%mesh_connectivities(j)%connectivity(3)

            X1(:) = TAB_MESHE_PARTICLE(kk)%node(nodea,:)-TAB_MESHE_PARTICLE(kk)%node(nodeb,:)
            X2(:) = TAB_MESHE_PARTICLE(kk)%node(nodea,:)-TAB_MESHE_PARTICLE(kk)%node(nodec,:)

            air_maille = 0.5*abs( X1(1)*X2(2)-X1(2)*X2(1) )
          end if
          TAB_MESHE_PARTICLE(kk)%aire = TAB_MESHE_PARTICLE(kk)%aire + air_maille
        end do


      else if (text == ' RBDY2') then
        k = k + 1
        read(5,*)
        read(5,'(29X, 3(5x,D14.7,2X))') coor1,coor2,coor3
        read(5,'(29X, 3(5x,D14.7,2X))') Vx,Vy,Vr
        TAB_RIGID_PARTICLE(k)%V(1) = Vx
        TAB_RIGID_PARTICLE(k)%V(2) = Vy
        TAB_RIGID_PARTICLE(k)%V(3) = Vr
        TAB_RIGID_PARTICLE(k)%centre(1)=coor1+TAB_RIGID_PARTICLE(k)%centre_ref(1)
        TAB_RIGID_PARTICLE(k)%centre(2)=coor2+TAB_RIGID_PARTICLE(k)%centre_ref(2)
        TAB_RIGID_PARTICLE(k)%centre(3)=coor3+TAB_RIGID_PARTICLE(k)%centre_ref(3)
        if ( TAB_RIGID_PARTICLE(k)%name == 'POLYG '.or.&
             TAB_RIGID_PARTICLE(k)%name == 'JONCx '  ) then
          do j=1,TAB_RIGID_PARTICLE(k)%nb_vertex
            TAB_RIGID_PARTICLE(k)%vertex(1,j) = TAB_RIGID_PARTICLE(k)%vertex_ref(1,j) * cos(coor3) - &
                                                       TAB_RIGID_PARTICLE(k)%vertex_ref(2,j) * sin(coor3) + &
                                                       TAB_RIGID_PARTICLE(k)%centre(1)
            TAB_RIGID_PARTICLE(k)%vertex(2,j) = TAB_RIGID_PARTICLE(k)%vertex_ref(1,j) * sin(coor3) + &
                                                       TAB_RIGID_PARTICLE(k)%vertex_ref(2,j) * cos(coor3) + &
                                                       TAB_RIGID_PARTICLE(k)%centre(2)
          end do
        end if
      end if
      if (i==nbParticules+nb_paroi) exit
    end if
  end do
end if

N_PLJCx = 0
N_PLPLx = 0
N_DKDKx = 0
N_DKJCx = 0
N_DKPLx = 0
N_CLALp = 0
N_CLJCx = 0
N_DKALp = 0

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
    if (text2 == '$icdan  PLPLx') then
      N_PLPLx = N_PLPLx +1
    end if
    if (text2 == '$icdan  PLJCx') then
      N_PLJCx = N_PLJCx +1
    end if
    if (text2 == '$icdan  DKDKx') then
      N_DKDKx = N_DKDKx +1
    end if
    if (text2 == '$icdan  DKJCx') then
      N_DKJCx = N_DKJCx +1
    end if
    if (text2 == '$icdan  DKPLx') then
      N_DKPLx = N_DKPLx +1
    end if
    if (text2 == '$icdan  CLALp') then
      N_CLALp = N_CLALp +1
    end if
    if (text2 == '$icdan  CLJCx') then
      N_CLJCx = N_CLJCx +1
    end if
    if (text2 == '$icdan  DKALp') then
      N_DKALp = N_DKALp +1
    end if
  end do
  close(4)

  if (allocated(TAB_CONTACT)) deallocate(TAB_CONTACT)
  allocate(TAB_CONTACT( N_PLPLx+N_PLJCx+N_DKDKx+N_DKJCx+N_DKPLx+N_CLALp+N_CLJCx+N_DKALp ))
  nb_plplx=0
  nb_pljcx=0
  nb_dkdkx=0
  nb_dkjcx=0
  nb_dkplx=0
  nb_clalp=0
  nb_cljcx=0
  nb_dkalp=0
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
    if (text2 == '$icdan  CLJCx') then
      read(4,*)
      read(4,'(8X,i5,9X,i5,16X,i5,9X,i5,2X,A5)') icdent,icdtact,ianent,iantact,statut
      read(4,'(29X, 3(5x,D14.7,2X))') rt,rn,rs
      read(4,'(29X, 3(5x,D14.7,2X))') vlt,vln,vls
      read(4,*)
      read(4,'(29X, 3(5x,D14.7,2X))') n1,n2,n3
      read(4,'(29X, 3(5x,D14.7,2X))') coor_ctc1,coor_ctc2,coor_ctc3

      nb_ligneCONTACT=nb_ligneCONTACT+1
      nb_cljcx = nb_cljcx + 1
      TAB_CONTACT(nb_cljcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp+N_DKJCx+N_PLJCx)%icdent=icdent
      TAB_CONTACT(nb_cljcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp+N_DKJCx+N_PLJCx)%ianent=ianent
      TAB_CONTACT(nb_cljcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp+N_DKJCx+N_PLJCx)%n(1)=n1
      TAB_CONTACT(nb_cljcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp+N_DKJCx+N_PLJCx)%n(2)=n2
      TAB_CONTACT(nb_cljcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp+N_DKJCx+N_PLJCx)%n(3)=n3
      TAB_CONTACT(nb_cljcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp+N_DKJCx+N_PLJCx)%t(1)=n2
      TAB_CONTACT(nb_cljcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp+N_DKJCx+N_PLJCx)%t(2)=-n1
      TAB_CONTACT(nb_cljcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp+N_DKJCx+N_PLJCx)%t(3)=0.D0
      TAB_CONTACT(nb_cljcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp+N_DKJCx+N_PLJCx)%coor_ctc(1)=coor_ctc1
      TAB_CONTACT(nb_cljcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp+N_DKJCx+N_PLJCx)%coor_ctc(2)=coor_ctc2
      TAB_CONTACT(nb_cljcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp+N_DKJCx+N_PLJCx)%coor_ctc(3)=coor_ctc3
      TAB_CONTACT(nb_cljcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp+N_DKJCx+N_PLJCx)%rn=rn
      TAB_CONTACT(nb_cljcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp+N_DKJCx+N_PLJCx)%rt=rt
      TAB_CONTACT(nb_cljcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp+N_DKJCx+N_PLJCx)%rs=rs
      TAB_CONTACT(nb_cljcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp+N_DKJCx+N_PLJCx)%vlt=vlt
      TAB_CONTACT(nb_cljcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp+N_DKJCx+N_PLJCx)%nature='CLJCx'
      TAB_CONTACT(nb_cljcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp+N_DKJCx+N_PLJCx)%statut=statut
      TAB_CONTACT(nb_cljcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp+N_DKJCx+N_PLJCx)%type  = 1
    end if
    if (text2 == '$icdan  PLJCx') then
      read(4,*)
      read(4,'(8X,i5,44X,i5,30X,A5)') icdent,ianent,statut
      read(4,'(29X, 3(5x,D14.7,2X))') rt,rn,rs
      read(4,'(29X, 3(5x,D14.7,2X))') vlt,vln,vls
      read(4,*)
      read(4,'(29X, 3(5x,D14.7,2X))') n1,n2,n3
      read(4,'(29X, 3(5x,D14.7,2X))') coor_ctc1,coor_ctc2,coor_ctc3

      nb_ligneCONTACT=nb_ligneCONTACT+1
      nb_pljcx = nb_pljcx + 1
      TAB_CONTACT(nb_pljcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp+N_DKJCx)%icdent=icdent
      TAB_CONTACT(nb_pljcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp+N_DKJCx)%ianent=ianent
      TAB_CONTACT(nb_pljcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp+N_DKJCx)%n(1)=n1
      TAB_CONTACT(nb_pljcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp+N_DKJCx)%n(2)=n2
      TAB_CONTACT(nb_pljcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp+N_DKJCx)%n(3)=n3
      TAB_CONTACT(nb_pljcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp+N_DKJCx)%t(1)=n2
      TAB_CONTACT(nb_pljcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp+N_DKJCx)%t(2)=-n1
      TAB_CONTACT(nb_pljcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp+N_DKJCx)%t(3)=0.D0
      TAB_CONTACT(nb_pljcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp+N_DKJCx)%coor_ctc(1)=coor_ctc1
      TAB_CONTACT(nb_pljcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp+N_DKJCx)%coor_ctc(2)=coor_ctc2
      TAB_CONTACT(nb_pljcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp+N_DKJCx)%coor_ctc(3)=coor_ctc3
      TAB_CONTACT(nb_pljcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp+N_DKJCx)%rn=rn
      TAB_CONTACT(nb_pljcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp+N_DKJCx)%rt=rt
      TAB_CONTACT(nb_pljcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp+N_DKJCx)%rs=rs
      TAB_CONTACT(nb_pljcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp+N_DKJCx)%vlt=vlt
      TAB_CONTACT(nb_pljcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp+N_DKJCx)%nature='PLJCx'
      TAB_CONTACT(nb_pljcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp+N_DKJCx)%statut=statut
      TAB_CONTACT(nb_pljcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp+N_DKJCx)%type  = 0
    end if
    if (text2 == '$icdan  DKJCx') then
      read(4,*)
      read(4,'(8X,i5,30X,i5,16X,A5)') icdent,ianent,statut
      read(4,'(29X, 3(5x,D14.7,2X))') rt,rn,rs
      read(4,'(29X, 3(5x,D14.7,2X))') vlt,vln,vls
      read(4,*)
      read(4,'(29X, 3(5x,D14.7,2X))') n1,n2,n3
      read(4,'(29X, 3(5x,D14.7,2X))') coor_ctc1,coor_ctc2,coor_ctc3

      nb_ligneCONTACT=nb_ligneCONTACT+1
      nb_dkjcx = nb_dkjcx + 1
      TAB_CONTACT(nb_dkjcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp)%icdent=icdent
      TAB_CONTACT(nb_dkjcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp)%ianent=ianent
      TAB_CONTACT(nb_dkjcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp)%n(1)=n1
      TAB_CONTACT(nb_dkjcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp)%n(2)=n2
      TAB_CONTACT(nb_dkjcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp)%n(3)=n3
      TAB_CONTACT(nb_dkjcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp)%t(1)=n2
      TAB_CONTACT(nb_dkjcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp)%t(2)=-n1
      TAB_CONTACT(nb_dkjcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp)%t(3)=0.D0
      TAB_CONTACT(nb_dkjcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp)%coor_ctc(1)=coor_ctc1
      TAB_CONTACT(nb_dkjcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp)%coor_ctc(2)=coor_ctc2
      TAB_CONTACT(nb_dkjcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp)%coor_ctc(3)=coor_ctc3
      TAB_CONTACT(nb_dkjcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp)%rn=rn
      TAB_CONTACT(nb_dkjcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp)%rt=rt
      TAB_CONTACT(nb_dkjcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp)%rs=rs
      TAB_CONTACT(nb_dkjcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp)%vlt=vlt
      TAB_CONTACT(nb_dkjcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp)%nature='DKJCx'
      TAB_CONTACT(nb_dkjcx+N_PLPLx+N_DKDKx+N_DKPLx+N_CLALp+N_DKALp)%statut=statut
      TAB_CONTACT(nb_dkjcx+N_DKDKx+N_PLPLx+N_DKPLx+N_CLALp+N_DKALp)%type  = 1
    end if
    if (text2 == '$icdan  DKDKx') then
      read(4,*)
      read(4,'(8X,i5,30X,i5,16X,A5)') icdent,ianent,statut
      read(4,'(29X, 3(5x,D14.7,2X))') rt,rn,rs
      read(4,'(29X, 3(5x,D14.7,2X))') vlt,vln,vls
      read(4,'(29X, (5x,D14.7,2X))')  gap
      read(4,'(29X, 3(5x,D14.7,2X))') n1,n2,n3
      read(4,'(29X, 3(5x,D14.7,2X))') coor_ctc1,coor_ctc2,coor_ctc3

      nb_ligneCONTACT=nb_ligneCONTACT+1
      nb_dkdkx = nb_dkdkx + 1
      TAB_CONTACT(nb_dkdkx)%icdent=icdent
      TAB_CONTACT(nb_dkdkx)%ianent=ianent
      TAB_CONTACT(nb_dkdkx)%n(1)=n1
      TAB_CONTACT(nb_dkdkx)%n(2)=n2
      TAB_CONTACT(nb_dkdkx)%n(3)=n3
      TAB_CONTACT(nb_dkdkx)%t(1)=n2
      TAB_CONTACT(nb_dkdkx)%t(2)=-n1
      TAB_CONTACT(nb_dkdkx)%t(3)=0.D0
      TAB_CONTACT(nb_dkdkx)%coor_ctc(1)=coor_ctc1
      TAB_CONTACT(nb_dkdkx)%coor_ctc(2)=coor_ctc2
      TAB_CONTACT(nb_dkdkx)%coor_ctc(3)=coor_ctc3
      TAB_CONTACT(nb_dkdkx)%rn=rn
      TAB_CONTACT(nb_dkdkx)%rt=rt
      TAB_CONTACT(nb_dkdkx)%rs=rs
      TAB_CONTACT(nb_dkdkx)%gap=gap
      TAB_CONTACT(nb_dkdkx)%vlt=vlt
      TAB_CONTACT(nb_dkdkx)%nature='DKDKx'
      TAB_CONTACT(nb_dkdkx)%type  = 1
      TAB_CONTACT(nb_dkdkx)%statut  = statut
    end if
    if (text2 == '$icdan  PLPLx') then
      read(4,*)
      read(4,'(8X,i5,44X,i5,30X,A5)') icdent,ianent,statut
      read(4,'(29X, 3(5x,D14.7,2X))') rt,rn,rs
      read(4,'(29X, 3(5x,D14.7,2X))') vlt,vln,vls
      read(4,'(29X, (5x,D14.7,2X))')  gap
      read(4,'(29X, 3(5x,D14.7,2X))') n1,n2,n3
      read(4,'(29X, 3(5x,D14.7,2X))') coor_ctc1,coor_ctc2,coor_ctc3

      nb_ligneCONTACT=nb_ligneCONTACT+1
      nb_plplx = nb_plplx + 1
      TAB_CONTACT(nb_plplx+N_DKDKx)%icdent=icdent
      TAB_CONTACT(nb_plplx+N_DKDKx)%ianent=ianent
      TAB_CONTACT(nb_plplx+N_DKDKx)%n(1)=n1
      TAB_CONTACT(nb_plplx+N_DKDKx)%n(2)=n2
      TAB_CONTACT(nb_plplx+N_DKDKx)%n(3)=n3
      TAB_CONTACT(nb_plplx+N_DKDKx)%t(1)=n2
      TAB_CONTACT(nb_plplx+N_DKDKx)%t(2)=-n1
      TAB_CONTACT(nb_plplx+N_DKDKx)%t(3)=0.D0
      TAB_CONTACT(nb_plplx+N_DKDKx)%coor_ctc(1)=coor_ctc1
      TAB_CONTACT(nb_plplx+N_DKDKx)%coor_ctc(2)=coor_ctc2
      TAB_CONTACT(nb_plplx+N_DKDKx)%coor_ctc(3)=coor_ctc3
      TAB_CONTACT(nb_plplx+N_DKDKx)%rn=rn
      TAB_CONTACT(nb_plplx+N_DKDKx)%rt=rt
      TAB_CONTACT(nb_plplx+N_DKDKx)%rs=rs
      TAB_CONTACT(nb_plplx+N_DKDKx)%vlt=vlt
      TAB_CONTACT(nb_plplx+N_DKDKx)%gap=gap
      TAB_CONTACT(nb_plplx+N_DKDKx)%nature='PLPLx'
      TAB_CONTACT(nb_plplx+N_DKDKx)%type  = 0
      TAB_CONTACT(nb_plplx+N_DKDKx)%statut  = statut
    end if
    if (text2 == '$icdan  DKPLx') then
      read(4,*)
      read(4,'(8X,i5,30X,i5,44X,A5)') icdent,ianent,statut
      read(4,'(29X, 3(5x,D14.7,2X))') rt,rn,rs
      read(4,'(29X, 3(5x,D14.7,2X))') vlt,vln,vls
      read(4,'(29X, (5x,D14.7,2X))')  gap
      read(4,'(29X, 3(5x,D14.7,2X))') n1,n2,n3
      read(4,'(29X, 3(5x,D14.7,2X))') coor_ctc1,coor_ctc2,coor_ctc3

      nb_ligneCONTACT=nb_ligneCONTACT+1
      nb_dkplx = nb_dkplx + 1
      TAB_CONTACT(nb_dkplx+N_DKDKx+N_PLPLx)%icdent=icdent
      TAB_CONTACT(nb_dkplx+N_DKDKx+N_PLPLx)%ianent=ianent
      TAB_CONTACT(nb_dkplx+N_DKDKx+N_PLPLx)%n(1)=n1
      TAB_CONTACT(nb_dkplx+N_DKDKx+N_PLPLx)%n(2)=n2
      TAB_CONTACT(nb_dkplx+N_DKDKx+N_PLPLx)%n(3)=n3
      TAB_CONTACT(nb_dkplx+N_DKDKx+N_PLPLx)%t(1)=n2
      TAB_CONTACT(nb_dkplx+N_DKDKx+N_PLPLx)%t(2)=-n1
      TAB_CONTACT(nb_dkplx+N_DKDKx+N_PLPLx)%t(3)=0.D0
      TAB_CONTACT(nb_dkplx+N_DKDKx+N_PLPLx)%coor_ctc(1)=coor_ctc1
      TAB_CONTACT(nb_dkplx+N_DKDKx+N_PLPLx)%coor_ctc(2)=coor_ctc2
      TAB_CONTACT(nb_dkplx+N_DKDKx+N_PLPLx)%coor_ctc(3)=coor_ctc3
      TAB_CONTACT(nb_dkplx+N_DKDKx+N_PLPLx)%rn=rn
      TAB_CONTACT(nb_dkplx+N_DKDKx+N_PLPLx)%rt=rt
      TAB_CONTACT(nb_dkplx+N_DKDKx+N_PLPLx)%rs=rs
      TAB_CONTACT(nb_dkplx+N_DKDKx+N_PLPLx)%vlt=vlt
      TAB_CONTACT(nb_dkplx+N_DKDKx+N_PLPLx)%gap=gap
      TAB_CONTACT(nb_dkplx+N_DKDKx+N_PLPLx)%nature='DKPLx'
      TAB_CONTACT(nb_dkplx+N_DKDKx+N_PLPLx)%type  = 1
      TAB_CONTACT(nb_dkplx+N_DKDKx+N_PLPLx)%statut  = statut
    end if

    if (text2 == '$icdan  CLALp') then
      read(4,*)
      read(4,'(8X,i5,9X,i5,16X,i5,9X,i5,9X,A5)') icdent,icdtact,ianent,iantact,statut
      read(4,'(29X, 3(5x,D14.7,2X))') rt,rn,rs
      read(4,'(29X, 3(5x,D14.7,2X))') vlt,vln,vls
      read(4,'(29X, (5x,D14.7,2X))')  gap
      read(4,'(29X, 3(5x,D14.7,2X))') n1,n2,n3
      read(4,'(29X, 3(5x,D14.7,2X))') coor_ctc1,coor_ctc2,coor_ctc3

      nb_ligneCONTACT=nb_ligneCONTACT+1
      nb_clalp = nb_clalp + 1
      TAB_CONTACT(nb_clalp+N_DKDKx+N_PLPLx+N_DKPLx)%icdent=icdent
      TAB_CONTACT(nb_clalp+N_DKDKx+N_PLPLx+N_DKPLx)%ianent=ianent
      TAB_CONTACT(nb_clalp+N_DKDKx+N_PLPLx+N_DKPLx)%ictacty=icdtact
      TAB_CONTACT(nb_clalp+N_DKDKx+N_PLPLx+N_DKPLx)%iantact=iantact
      TAB_CONTACT(nb_clalp+N_DKDKx+N_PLPLx+N_DKPLx)%n(1)=n1
      TAB_CONTACT(nb_clalp+N_DKDKx+N_PLPLx+N_DKPLx)%n(2)=n2
      TAB_CONTACT(nb_clalp+N_DKDKx+N_PLPLx+N_DKPLx)%n(3)=n3
      TAB_CONTACT(nb_clalp+N_DKDKx+N_PLPLx+N_DKPLx)%t(1)=n2
      TAB_CONTACT(nb_clalp+N_DKDKx+N_PLPLx+N_DKPLx)%t(2)=-n1
      TAB_CONTACT(nb_clalp+N_DKDKx+N_PLPLx+N_DKPLx)%t(3)=0.D0
      TAB_CONTACT(nb_clalp+N_DKDKx+N_PLPLx+N_DKPLx)%coor_ctc(1)=coor_ctc1
      TAB_CONTACT(nb_clalp+N_DKDKx+N_PLPLx+N_DKPLx)%coor_ctc(2)=coor_ctc2
      TAB_CONTACT(nb_clalp+N_DKDKx+N_PLPLx+N_DKPLx)%coor_ctc(3)=coor_ctc3
      TAB_CONTACT(nb_clalp+N_DKDKx+N_PLPLx+N_DKPLx)%rn=rn
      TAB_CONTACT(nb_clalp+N_DKDKx+N_PLPLx+N_DKPLx)%rt=rt
      TAB_CONTACT(nb_clalp+N_DKDKx+N_PLPLx+N_DKPLx)%rs=rs
      TAB_CONTACT(nb_clalp+N_DKDKx+N_PLPLx+N_DKPLx)%vlt=vlt
      TAB_CONTACT(nb_clalp+N_DKDKx+N_PLPLx+N_DKPLx)%gap=gap
      TAB_CONTACT(nb_clalp+N_DKDKx+N_PLPLx+N_DKPLx)%nature='CLALp'
      TAB_CONTACT(nb_clalp+N_DKDKx+N_PLPLx+N_DKPLx)%type  = 1
      TAB_CONTACT(nb_clalp+N_DKDKx+N_PLPLx+N_DKPLx)%statut  = statut
    end if

    if (text2 == '$icdan  DKALp') then
      read(4,*)
      read(4,'(8X,i5,9X,i5,16X,i5,9X,i5,9X,A5)') icdent,icdtact,ianent,iantact,statut
      read(4,'(29X, 3(5x,D14.7,2X))') rt,rn,rs
      read(4,'(29X, 3(5x,D14.7,2X))') vlt,vln,vls
      read(4,'(29X, (5x,D14.7,2X))')  gap
      read(4,'(29X, 3(5x,D14.7,2X))') n1,n2,n3
      read(4,'(29X, 3(5x,D14.7,2X))') coor_ctc1,coor_ctc2,coor_ctc3

      nb_ligneCONTACT=nb_ligneCONTACT+1
      nb_dkalp = nb_dkalp + 1
      TAB_CONTACT(nb_dkalp+N_DKDKx+N_PLPLx+N_DKPLx+N_CLALp)%icdent=icdent
      TAB_CONTACT(nb_dkalp+N_DKDKx+N_PLPLx+N_DKPLx+N_CLALp)%ianent=ianent
      TAB_CONTACT(nb_dkalp+N_DKDKx+N_PLPLx+N_DKPLx+N_CLALp)%ictacty=icdtact
      TAB_CONTACT(nb_dkalp+N_DKDKx+N_PLPLx+N_DKPLx+N_CLALp)%iantact=iantact
      TAB_CONTACT(nb_dkalp+N_DKDKx+N_PLPLx+N_DKPLx+N_CLALp)%n(1)=n1
      TAB_CONTACT(nb_dkalp+N_DKDKx+N_PLPLx+N_DKPLx+N_CLALp)%n(2)=n2
      TAB_CONTACT(nb_dkalp+N_DKDKx+N_PLPLx+N_DKPLx+N_CLALp)%n(3)=n3
      TAB_CONTACT(nb_dkalp+N_DKDKx+N_PLPLx+N_DKPLx+N_CLALp)%t(1)=n2
      TAB_CONTACT(nb_dkalp+N_DKDKx+N_PLPLx+N_DKPLx+N_CLALp)%t(2)=-n1
      TAB_CONTACT(nb_dkalp+N_DKDKx+N_PLPLx+N_DKPLx+N_CLALp)%t(3)=0.D0
      TAB_CONTACT(nb_dkalp+N_DKDKx+N_PLPLx+N_DKPLx+N_CLALp)%coor_ctc(1)=coor_ctc1
      TAB_CONTACT(nb_dkalp+N_DKDKx+N_PLPLx+N_DKPLx+N_CLALp)%coor_ctc(2)=coor_ctc2
      TAB_CONTACT(nb_dkalp+N_DKDKx+N_PLPLx+N_DKPLx+N_CLALp)%coor_ctc(3)=coor_ctc3
      TAB_CONTACT(nb_dkalp+N_DKDKx+N_PLPLx+N_DKPLx+N_CLALp)%rn=rn
      TAB_CONTACT(nb_dkalp+N_DKDKx+N_PLPLx+N_DKPLx+N_CLALp)%rt=rt
      TAB_CONTACT(nb_dkalp+N_DKDKx+N_PLPLx+N_DKPLx+N_CLALp)%rs=rs
      TAB_CONTACT(nb_dkalp+N_DKDKx+N_PLPLx+N_DKPLx+N_CLALp)%vlt=vlt
      TAB_CONTACT(nb_dkalp+N_DKDKx+N_PLPLx+N_DKPLx+N_CLALp)%gap=gap
      TAB_CONTACT(nb_dkalp+N_DKDKx+N_PLPLx+N_DKPLx+N_CLALp)%nature='DKALp'
      TAB_CONTACT(nb_dkalp+N_DKDKx+N_PLPLx+N_DKPLx+N_CLALp)%type  = 1
      TAB_CONTACT(nb_dkalp+N_DKDKx+N_PLPLx+N_DKPLx+N_CLALp)%statut  = statut
    end if

  end do
end if
!----- Ici il faut definir le type de chaque contact ----------
do i=1,nb_ligneCONTACT-1
  if (TAB_CONTACT(i)%type>0) cycle

  if ( (TAB_CONTACT(i)%icdent == TAB_CONTACT(i+1)%icdent) .and.&
       (TAB_CONTACT(i)%ianent == TAB_CONTACT(i+1)%ianent)) then
       TAB_CONTACT(i)%type=2
       TAB_CONTACT(i+1)%type=2
  else
    TAB_CONTACT(i)%type=1
    if (i==nb_ligneCONTACT-1) TAB_CONTACT(i+1)%type=1
  end if
end do
!----------------------------------------------------------------
n_simple = 0
n_double = 0
first_double=.true.
nb_ligneCONTACT_PARTICLE = 0
do i=1,nb_ligneCONTACT
  if (TAB_CONTACT(i)%type==1) n_simple = n_simple + 1
  if (TAB_CONTACT(i)%type==2) n_double = n_double + 1
end do
if (allocated(TAB_CONTACT_PARTICLE)) deallocate(TAB_CONTACT_PARTICLE)
allocate(TAB_CONTACT_PARTICLE(n_simple+n_double))
do i=1,nb_ligneCONTACT
  if (TAB_CONTACT(i)%type==1) then
    nb_ligneCONTACT_PARTICLE                                = nb_ligneCONTACT_PARTICLE+1
    TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%icdent   = TAB_CONTACT(i)%icdent
    TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%ianent   = TAB_CONTACT(i)%ianent
    TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%n        = TAB_CONTACT(i)%n
    TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%t        = TAB_CONTACT(i)%t
    TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%coor_ctc = TAB_CONTACT(i)%coor_ctc
    TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%rn       = TAB_CONTACT(i)%rn
    TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%rt       = TAB_CONTACT(i)%rt
    TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%rs       = TAB_CONTACT(i)%rs
    TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%gap      = TAB_CONTACT(i)%gap
    TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%F        = &
    TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%rn * TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%n +&
    TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%rt * TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%t
    TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%vlt      = TAB_CONTACT(i)%vlt
    TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%vln      = TAB_CONTACT(i)%vln
    if (abs(TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%vlt) < 0.0000000001)&
      TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%statut   = 'stick'
    if (abs(TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%vlt) > 0.0000000001)&
      TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%statut   = 'slibw'
    TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%nature   = TAB_CONTACT(i)%nature
    TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%type     = 1
  end if
  if (TAB_CONTACT(i)%type==2) then
    if (first_double) then
      coor_ctc1_temp = TAB_CONTACT(i)%coor_ctc(1)
      coor_ctc2_temp = TAB_CONTACT(i)%coor_ctc(2)
      coor_ctc3_temp = TAB_CONTACT(i)%coor_ctc(3)
      rn_temp        = TAB_CONTACT(i)%rn
      rt_temp        = TAB_CONTACT(i)%rt
      rs_temp        = TAB_CONTACT(i)%rs
      vlt_temp       = TAB_CONTACT(i)%vlt
      vln_temp       = TAB_CONTACT(i)%vln
      gap_temp       = TAB_CONTACT(i)%gap
      first_double   = .false.
      cycle
    else if (.not. first_double) then
      nb_ligneCONTACT_PARTICLE                                    = nb_ligneCONTACT_PARTICLE+1
      TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%icdent       = TAB_CONTACT(i)%icdent
      TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%ianent       = TAB_CONTACT(i)%ianent
      TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%n            = TAB_CONTACT(i)%n
      TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%t            = TAB_CONTACT(i)%t
      TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%nature       = TAB_CONTACT(i)%nature

      if (rn_temp>0. .and. TAB_CONTACT(i)%rn >0.) then
        TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%coor_ctc(1)  = (TAB_CONTACT(i)%coor_ctc(1) + coor_ctc1_temp) /2
        TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%coor_ctc(2)  = (TAB_CONTACT(i)%coor_ctc(2) + coor_ctc2_temp) /2
        TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%coor_ctc(3)  = (TAB_CONTACT(i)%coor_ctc(3) + coor_ctc3_temp) /2
        TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%rn           = TAB_CONTACT(i)%rn + rn_temp
        TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%rt           = TAB_CONTACT(i)%rt + rt_temp
        TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%rs           = TAB_CONTACT(i)%rs + rs_temp
        TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%gap          = (TAB_CONTACT(i)%gap+gap_temp)/2
        TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%vlt          = (TAB_CONTACT(i)%vlt + vlt_temp)/2
        TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%vln          = (TAB_CONTACT(i)%vln + vln_temp)/2
        TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%type         = 2

      else if (rn_temp == 0. .and. TAB_CONTACT(i)%rn >0.) then
        TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%coor_ctc(1)  = TAB_CONTACT(i)%coor_ctc(1)
        TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%coor_ctc(2)  = TAB_CONTACT(i)%coor_ctc(2)
        TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%coor_ctc(3)  = TAB_CONTACT(i)%coor_ctc(3)
        TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%rn           = TAB_CONTACT(i)%rn
        TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%rt           = TAB_CONTACT(i)%rt
        TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%rs           = TAB_CONTACT(i)%rs
        TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%gap          = TAB_CONTACT(i)%gap
        TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%vlt          = TAB_CONTACT(i)%vlt
        TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%vln          = TAB_CONTACT(i)%vln
        TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%type         = 1

      else if (rn_temp >0.   .and. TAB_CONTACT(i)%rn ==0.) then
        TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%coor_ctc(1)  = coor_ctc1_temp
        TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%coor_ctc(2)  = coor_ctc2_temp
        TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%coor_ctc(3)  = coor_ctc3_temp
        TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%rn           = rn_temp
        TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%rt           = rt_temp
        TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%rs           = rs_temp
        TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%gap          = gap_temp
        TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%vlt          = vlt_temp
        TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%vln          = vln_temp
        TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%type         = 1
      else if (rn_temp ==0.   .and. TAB_CONTACT(i)%rn ==0.) then
        TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%coor_ctc(1)  = (TAB_CONTACT(i)%coor_ctc(1) + coor_ctc1_temp) /2
        TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%coor_ctc(2)  = (TAB_CONTACT(i)%coor_ctc(2) + coor_ctc2_temp) /2
        TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%coor_ctc(3)  = (TAB_CONTACT(i)%coor_ctc(3) + coor_ctc3_temp) /2
        TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%rn           = TAB_CONTACT(i)%rn + rn_temp
        TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%rt           = TAB_CONTACT(i)%rt + rt_temp
        TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%rs           = TAB_CONTACT(i)%rs + rs_temp
        TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%gap          = (TAB_CONTACT(i)%gap+gap_temp)/2
        TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%vlt          = (TAB_CONTACT(i)%vlt + vlt_temp)/2
        TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%vln          = (TAB_CONTACT(i)%vln + vln_temp)/2
      end if


      TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%F            = &
      TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%rn * TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%n +&
      TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%rt * TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%t
      if (abs(TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%vlt) < 0.000000000001)&
        TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%statut   = 'stick'
      if (abs(TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%vlt) > 0.000000000001)&
        TAB_CONTACT_PARTICLE(nb_ligneCONTACT_PARTICLE)%statut   = 'slibw'
      first_double                                          = .true.
    end if
  end if
end do
end subroutine read_Vloc_dof_bodies

!===========================================================================
!===========================================================================
!===========================================================================
! OTHERS ROUTINES
!===========================================================================
!===========================================================================
!===========================================================================


!==============================================================================
! Compute size box
!==============================================================================
 
subroutine compute_size_box
implicit none
integer                       ::  i,j
real(kind=8)                  ::  xmax,ymax
real(kind=8)                  ::  xmin,ymin

xmax = 0.D0
ymax = 0.D0
xmin = 1000000000.D0
ymin = 1000000000.D0
do i=1,nb_RIGID
  if (TAB_RIGID_PARTICLE(i)%name == 'DISKx ') then
    xmax = max( xmax,TAB_RIGID_PARTICLE(i)%centre(1)+TAB_RIGID_PARTICLE(i)%Rayon )
    ymax = max( ymax,TAB_RIGID_PARTICLE(i)%centre(2)+TAB_RIGID_PARTICLE(i)%Rayon )
    xmin = min( xmin,TAB_RIGID_PARTICLE(i)%centre(1)-TAB_RIGID_PARTICLE(i)%Rayon )
    ymin = min( ymin,TAB_RIGID_PARTICLE(i)%centre(2)-TAB_RIGID_PARTICLE(i)%Rayon )
  end if
  if (TAB_RIGID_PARTICLE(i)%name == 'POLYG ') then
    do j=1,TAB_RIGID_PARTICLE(i)%nb_vertex
      xmax = max( xmax,TAB_RIGID_PARTICLE(i)%vertex(1,j) )
      ymax = max( ymax,TAB_RIGID_PARTICLE(i)%vertex(2,j) )
      xmin = min( xmin,TAB_RIGID_PARTICLE(i)%vertex(1,j) )
      ymin = min( ymin,TAB_RIGID_PARTICLE(i)%vertex(2,j) )
    end do
  end if
end do

do i=1,nb_DEFOR
  do j=1,TAB_MESHE_PARTICLE(i)%nb_node
    xmax = max( xmax,TAB_MESHE_PARTICLE(i)%node(j,1) )
    ymax = max( ymax,TAB_MESHE_PARTICLE(i)%node(j,2) )
    xmin = min( xmin,TAB_MESHE_PARTICLE(i)%node(j,1) )
    ymin = min( ymin,TAB_MESHE_PARTICLE(i)%node(j,2) )
  end do
end do

if (the_first_time) then
  hauteur0  = abs(ymax-ymin)
  longueur0 = abs(xmax-xmin)
  the_first_time = .false.
end if

hauteur  = abs(ymax-ymin)
longueur = abs(xmax-xmin)

epsilon1 = log(abs( hauteur-hauteur0 )   / hauteur0 + 1)
epsilon2 = log(abs( longueur-longueur0 ) / longueur0 + 1)

epsilon_q = max(epsilon1,epsilon2)-min(epsilon1,epsilon2)
end subroutine compute_size_box

!====================================================================
! Initialize PostScript
!====================================================================
subroutine  psinit()

99    format(a23)
98    format(a30)
97    format(a7,i7,a16)
96    format(a14)
95    format(a17)
94    format(a7,i3,1x,i3)
write(1008,96) '%!PS-Adobe-2.0'
write(1008,95) '%%Creator: EA    '
write(1008,95) '%%Title: PS_ONE_FILE.ps'
write(1008,95) '%%Pages: (atend) '
write(1008,98) '%%BoundingBox:  80 230 520 720'
write(1008,95) '%%EndComments    '
write(1008,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
write(1008,*) '% Procedure cprint: texte centre'
write(1008,*) '% (texte) x y cprint  :'
write(1008,*) '/cprint{moveto dup stringwidth pop 2 div neg 0'
write(1008,*) '        rmoveto show}def'
write(1008,*) '/Times-Roman findfont 0.8 scalefont setfont'
write(1008,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
write(1008,*) '% Procedure C: trace un cercle de centre xy,'
write(1008,*) '% de rayon r et niveau de gris g  :'
write(1008,*) '% x y r g C '
write(1008,*) '/C{newpath 4 1 roll 0 360 arc gsave setgray fill'
write(1008,*) '   grestore stroke}def'
write(1008,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
write(1008,*) '% Procedure C2: trace un cercle de centre xy,'
write(1008,*) '% de rayon radius et niveau de couleur r g b :'
write(1008,*) '% x y radius r g b C2'
write(1008,*) '/C2 {6 dict begin /b exch def /g exch def /r exch def'
write(1008,*) ' /radius exch def /y exch def /x exch def gsave newpath'
write(1008,*) ' r g b setrgbcolor x y radius 0 360 arc fill grestore end} def'
write(1008,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
!      write(1008,*) '% Procedure CC: trace un cercle couleur de centre xy,'
!      write(1008,*) '% de rayon r et niveau de gris g  :'
!      write(1008,*) '% x y r r g b CC '
!      write(1008,*) '/CC{newpath 4 1 roll 0 360 arc gsave setrgbcolor fill'
!      write(1008,*) '   grestore stroke}def'
!      write(1008,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
!      write(1008,*) '% Procedure Ca: i.e. C plus un rayon'
!      write(1008,*) '% trace selon l''angle t'
!      write(1008,*) '% t x y r g Ca'
!      write(1008,*) '/Ca{4 1 roll 3 copy 7 -1 roll C'
!      write(1008,*) '   /r exch def newpath moveto dup cos r mul'
!      write(1008,*) '   exch sin r mul'
!      write(1008,*) '   rlineto stroke} def'
!      write(1008,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
!      write(1008,*) '% Procedure Lis: ligne de largeur fixe'
!      write(1008,*) '% du point 0 au point 1'
!      write(1008,*) '% x0 y0 x1 y1 Lis'
!      write(1008,*) '/tmx 20 def'
!      write(1008,*) '/Lis{gsave newpath moveto lineto'
!      write(1008,*) '    stroke grestore}def'
!      write(1008,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
write(1008,*) '% Procedure Li: ligne intercentre de largeur donnee'
write(1008,*) '% largeur t, du point 0 au point 1'
write(1008,*) '% x0 y0 x1 y1 t r g b Li'
write(1008,*) '/tmx 60 def'
write(1008,*) '/Li{gsave setrgbcolor tmx mul setlinewidth newpath moveto lineto'
write(1008,*) '    stroke grestore}def'
write(1008,*) '1 setlinecap'
write(1008,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
write(1008,*) '% Procedure polyg: grey level polygones'
write(1008,*) '/polyg{gsave setgray 1 sub /nco exch def newpath'
write(1008,*) 'moveto 1 1 nco{pop lineto}for closepath fill'
write(1008,*) 'grestore}def'
write(1008,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
write(1008,*) '% Procedure polyc: coloured polygones'
write(1008,*) '/polyc{gsave setrgbcolor 1 sub /nco exch def newpath'
write(1008,*) 'moveto 1 1 nco{pop lineto}for closepath fill'
write(1008,*) 'grestore}def'
write(1008,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
write(1008,*) '% Procedure poly: polygone'
write(1008,*) '/poly{1 sub /nco exch def newpath'
write(1008,*) ' moveto 1 1 nco{pop lineto}for closepath stroke}def'
write(1008,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
!      write(1008,*) '% Procedure Ve: Vecteur'
!      write(1008,*) '% module argument x0 y0 Ve'
!      write(1008,*) '% x0 y0 depart du vecteur'
!      write(1008,*) '/fref{2 setlinewidth'
!      write(1008,*) '      newpath 0 0 moveto 86 0 lineto stroke'
!      write(1008,*) '      newpath 86 0 moveto 80 -6 lineto 100 0 lineto'
!      write(1008,*) '      80 6 lineto closepath fill}def'
!      write(1008,*) '/Ve{gsave translate rotate vesca mul '
!      write(1008,*) '    dup scale fref grestore}def'
!      write(1008,*) '/vesca 5 def'
!      write(1008,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
!      write(1008,*) '% Procedure Cg: Contact glissant'
!      write(1008,*) '% ang x0 y0 Cg'
!      write(1008,*) '% x0 y0 point de contact'
!      write(1008,*) '% ang: angle de la normale au contact'
!      write(1008,*) '%  Le symbole est une bande blanche de largeur'
!      write(1008,*) '%  2*cgwid et d''espacement 2*cgspa'
!      write(1008,*) '/Cg{gsave translate rotate 0.5 dup scale'
!      write(1008,*) ' newpath cgspa neg dup cgwid neg moveto cgwid lineto'
!      write(1008,*) ' cgspa dup cgwid lineto cgwid neg lineto closepath'
!      write(1008,*) ' 1 setgray fill 0 setgray 0 setlinewidth newpath'
!      write(1008,*) ' cgspa neg dup cgwid neg moveto cgwid lineto'
!      write(1008,*) ' cgspa dup cgwid moveto cgwid neg lineto stroke'
!      write(1008,*) ' grestore} def'
!      write(1008,*) '/cgwid{50}def'
!      write(1008,*) '/cgspa{10}def'
!      write(1008,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
!      write(1008,*) '% Procedure Cng: Contact non-glissant'
!      write(1008,*) '% ang x0 y0 Cng'
!      write(1008,*) '% x0 y0 point de contact'
!      write(1008,*) '% ang: angle de la normale au contact'
!      write(1008,*) '%  Le symbole est un disque noir de rayon cngsiz'
!      write(1008,*) '/Cng{gsave translate rotate newpath 0 0 cngsiz'
!      write(1008,*) ' 0 360 arc fill grestore} def'
!      write(1008,*) '/cngsiz{0.2}def'
!      write(1008,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
write(1008,*) '() 300 250 cprint'
write(1008,*) 90,' ',280,' ','translate 0.2 dup scale ',0.,' rotate'
write(1008,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'

return
endsubroutine psinit

!===========================================================================
!That's all folks!
!===========================================================================
end program post2D













