! Program written to simulate a 2D lid driven case

include 'pressure.f90' ! pressure subroutine
include 'boundary_conditions.f90' ! boundary conditions subroutine
include 'animation.f90' ! animation subroutine
include 'output.f90' ! file writing subroutine

program LD2DNS
	implicit none

	integer :: i,j,nx,ny,nx1,nx2,ny1,ny2,t,tmax,nplot
	real :: dx,dy,lx,ly,dt,Re,Ga,Delta,gma
	real, dimension(:), allocatable :: x,y
	real, dimension(:,:), allocatable :: u1,u2,v1,v2,p1,p2,F,G

	real, dimension(:,:), allocatable :: du2dx,duvdy,d2udx2,d2udy2
	real, dimension(:,:), allocatable :: dv2dy,dvudx,d2vdx2,d2vdy2

	real, dimension(:,:), allocatable :: flx,fly

	real, dimension(:), allocatable :: x1,y1

	real, dimension(:,:), allocatable :: psi,zeta

	lx=1.0 ; ly=1.0 ! Dimensionless size of the domain
	nx=64
	ny=64
	dx=lx/float(nx) ; dy=ly/float(ny)

	nx1=1 ; nx2=nx ; ny1=1 ; ny2=ny !Definición de las fronteras

	Re=70.0
	Ga=0.7
	gma=0.1

	allocate(x(0:nx+1),y(0:ny+1))
	allocate(u1(0:nx+1,0:ny+1),u2(0:nx+1,0:ny+1),v1(0:nx+1,0:ny+1),v2(0:nx+1,0:ny+1),&
					 F(0:nx+1,0:ny+1),G(0:nx+1,0:ny+1),p1(0:nx+1,0:ny+1),p2(0:nx+1,0:ny+1),&
					 flx(0:nx+1,0:ny+1),fly(0:nx+1,0:ny+1),psi(0:nx+1,0:ny+1),zeta(0:nx+1,0:ny+1))
	
	allocate(du2dx(0:nx+1,0:ny+1),duvdy(0:nx+1,0:ny+1),&
					 d2udx2(0:nx+1,0:ny+1),d2udy2(0:nx+1,0:ny+1))
	allocate(dv2dy(0:nx+1,0:ny+1),dvudx(0:nx+1,0:ny+1),&
					 d2vdx2(0:nx+1,0:ny+1),d2vdy2(0:nx+1,0:ny+1))


! Construction of the numerical grid
	x(nx1-1)=0.0
	x(nx1)=0.5*dx
	do i=(nx1+1),nx2
		x(i)=0.5*dx+(i-nx1)*dx
	enddo
	x(nx2+1)=lx

	y(ny1-1)=0.0
	y(ny1)=0.5*dy
	do j=(ny1+1),ny2
		y(j)=0.5*dy+(j-ny1)*dy
	enddo
	y(ny2+1)=Ly

! Temporal step size definition
	dt=gma*(Re/2.)*( ( (1./dx)**2 + (1./dy)**2 )**(-1) )
	print *, 'step size dt=',dt
	tmax=10000 ! Number of steps during the simulation
	nplot=500 ! animation parameter

! Initial conditions
	u1=0.0 ; v1=0.0 ; p1=0.0
	u2=u1 ; v2=v1 ; p2=p1

! Main time loop begins
DO t=1,Tmax !==========
!======================

	call boundary_conditions(u1,u2,v1,v2,nx,ny)

! F function
	forall(i=nx1:nx2-1,j=ny1:ny2)

	du2dx(i,j)=(0.25/dx)*( (u1(i,j)+u1(i+1,j))**2 &
		- (u1(i-1,j)+u1(i,j))**2) &
		+(0.25*Ga/dx)* &
		( ABS(u1(i,j)+u1(i+1,j))*(u1(i,j)-u1(i+1,j))&
		 - ABS(u1(i-1,j)+u1(i,j))*(u1(i-1,j)-u1(i,j)) )

	duvdy(i,j)=(0.25/dy)*( (v1(i,j)+v1(i+1,j))*(u1(i,j)+u1(i,j+1)) &
		- (v1(i,j-1)+v1(i+1,j-1))*(u1(i,j-1)+u1(i,j)) ) &
		+(0.25*Ga/dy)* &
		( ABS(v1(i,j)+v1(i+1,j))*(u1(i,j)-u1(i,j+1)) &
		- ABS(v1(i,j-1)+v1(i+1,j-1))*(u1(i,j-1)-u1(i,j)) )

	d2udx2(i,j)=(1.0/(dx*dx))*(u1(i+1,j)-2.0*u1(i,j)+u1(i-1,j))

	d2udy2(i,j)=(1.0/(dy*dy))*(u1(i,j+1)-2.0*u1(i,j)+u1(i,j-1))

	flx(i,j)=0.0

	F(i,j)=u1(i,j)+dt*( (1.0/Re)*(d2udx2(i,j) + d2udy2(i,j))&
		 - du2dx(i,j) - duvdy(i,j) + flx(i,j) )
	end forall

! G function
	forall(i=nx1:nx2,j=ny1:ny2-1)

	dv2dy(i,j)=(0.25/dy)*((v1(i,j)+v1(i,j+1))**2 &
		- (v1(i,j-1)+v1(i,j))**2) &
		+ (0.25*Ga/dy)* &
		( ABS(v1(i,j)+v1(i,j+1))*(v1(i,j)-v1(i,j+1)) &
		- ABS(v1(i,j-1)+v1(i,j))*(v1(i,j-1)-v1(i,j)) )

	dvudx(i,j)=(0.25/dx)*((u1(i,j)+u1(i,j+1))*(v1(i,j)+v1(i+1,j))&
		 - (u1(i-1,j)+u1(i-1,j+1))*(v1(i-1,j)+v1(i,j)))&
		 + (0.25*Ga/dx)*( ABS(u1(i,j)+u1(i,j+1))*(v1(i,j)-v1(i+1,j)) &
		 - ABS(u1(i-1,j)+u1(i-1,j+1))*(v1(i-1,j)-v1(i,j)) )

	d2vdx2(i,j)=(1.0/(dx*dx))*(v1(i+1,j)-2.0*v1(i,j)+v1(i-1,j))

	d2vdy2(i,j)=(1.0/(dy*dy))*(v1(i,j+1)-2.0*v1(i,j)+v1(i,j-1))

	fly(i,j)=0.0

	G(i,j)=v1(i,j)+dt*( (1.0/Re)*(d2vdx2(i,j)+d2vdy2(i,j))&
		 -dvudx(i,j)-dv2dy(i,j)+fly(i,j) )
	end forall

! Calculating F and G at the boundaries
	forall(j=ny1:ny2)
		F(0,j)=u1(0,j) ; F(nx2,j)=u1(nx2,j)
	end forall

	forall(i=nx1:nx2)
		G(i,0)=v1(i,0) ; G(i,ny2)=v1(i,ny2)
	end forall

! Solving the Poisson equation for the pressure
	call pressure(p1,p2,F,G,nx,ny,dx,dy,dt)

! Calculating future velocities
	forall(i=nx1:nx2-1,j=ny1:ny2)
		u2(i,j)=F(i,j)-(dt/dx)*(p2(i+1,j)-p2(i,j))
	end forall

	forall(i=nx1:nx2,j=ny1:ny2-1)
		v2(i,j)=G(i,j)-(dt/dy)*(p2(i,j+1)-p2(i,j))
	end forall

! Computing flow properties
	! streamlines
		do i=0,nx
			psi(i,0) = 0.0
			do j=1,ny-1
				psi(i,j)=psi(i,j-1) + u2(i,j)*dy
			end do
		end do

	! vorticity
	 	do i=1,nx
			do j=1,ny
				zeta(i,j) = (1.0/dy)*(u2(i,j+1)-u2(i,j)) - (1.0/dx)*(v2(i+1,j)-v2(i,j))
			end do
		end do

! Printing values on terminal
	if(mod(t,nplot)==0) THEN
		Delta=MAXVAL(ABS(u2-u1))
		write(*,*) 'T=',t,'U1=', u1(nx/2,ny/2),'V1=', v1(nx/2,ny/2), 'E=', Delta
		call animation(x,y,u2,v2,p2,nx,ny,nplot,t,lx,ly,psi,zeta)
	end if

! Assigning old velocities to the new velocities
	u1=u2
	v1=v2

! Main time loop ends
END DO !=============
!====================

! Computing velocities inside the centers of the cell
	forall(i=nx1:nx2,j=ny1:ny2)
		u1(i,j)=0.5*(u2(i,j)+u2(i-1,j))
    v1(I,J)=0.5*(v2(i,j)+v2(i,j-1))
	end forall

	call output(x,y,nx,ny,nx1,nx2,ny1,ny2,u1,u2,v1,v2,p1,p2,&
	                  psi,zeta,lx,ly)

! Main program ends =====
end program LD2DNS !=====
!========================
