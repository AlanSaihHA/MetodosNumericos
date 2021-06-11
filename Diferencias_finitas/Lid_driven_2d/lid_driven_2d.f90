program LD2DNS
implicit none

integer :: i,j,nx,ny,nx1,nx2,ny1,ny2,t,tmax,nplot
real :: dx,dy,lx,ly,dt,Re,Ga,Delta,gma
real, dimension (:), allocatable :: x,y
real, dimension (:,:), allocatable :: u1,u2,v1,v2,p1,p2,F,G

real, dimension (:,:), allocatable :: du2dx,duvdy,d2udx2,d2udy2
real, dimension (:,:), allocatable :: dv2dy,dvudx,d2vdx2,d2vdy2

real, dimension (:,:), allocatable :: flx,fly

real, dimension (:), allocatable :: x1,y1

!real, dimension(:,:), allocatable :: psi, zeta

lx=1.0 ; ly=1.0   !Dimensionless size of the domain
nx=64
ny=64
dx=lx/float(nx) ; dy=ly/float(ny)

nx1=1 ; nx2=nx ; ny1=1 ; ny2=ny     !Definición de las fronteras

Re=70.0
Ga=0.7
gma=0.1

allocate(x(0:nx+1),y(0:ny+1))
allocate(u1(0:nx+1,0:ny+1),u2(0:nx+1,0:ny+1),v1(0:nx+1,0:ny+1),v2(0:nx+1,0:ny+1), &
	F(0:nx+1,0:ny+1),G(0:nx+1,0:ny+1),p1(0:nx+1,0:ny+1),p2(0:nx+1,0:ny+1), &
	flx(0:nx+1,0:ny+1),fly(0:nx+1,0:ny+1))
allocate(du2dx(0:nx+1,0:ny+1),duvdy(0:nx+1,0:ny+1), &
	d2udx2(0:nx+1,0:ny+1),d2udy2(0:nx+1,0:ny+1))
allocate(dv2dy(0:nx+1,0:ny+1),dvudx(0:nx+1,0:ny+1), &
	d2vdx2(0:nx+1,0:ny+1),d2vdy2(0:nx+1,0:ny+1))


endprogram LD2DNS
