! Program written to simulate a 2D lid driven case

!include 'pressure.f90' ! pressure subroutine
!include 'boundary_conditions.f90' ! boundary conditions subroutine
!include 'animation.f90' ! animation subroutine
!include 'output.f90' ! file writing subroutine

program LD2DNS
	implicit none

	integer :: i,j,nx,ny,nx1,nx2,ny1,ny2,t,tmax,nplot
	real :: dx,dy,lx,ly,dt,Re,Ga,Delta,gma,xc,yc
	real, dimension(:), allocatable :: x,y
	real, dimension(:,:), allocatable :: u1,u2,v1,v2,p1,p2,F,G

	real, dimension(:,:), allocatable :: du2dx,duvdy,d2udx2,d2udy2
	real, dimension(:,:), allocatable :: dv2dy,dvudx,d2vdx2,d2vdy2

	real, dimension(:,:), allocatable :: flx,fly

	real, dimension(:), allocatable :: x1,y1

	real, dimension(:,:), allocatable :: psi,zeta, B0z, Bz

	lx=5.0 ; ly=5.0 ! Dimensionless size of the domain
	nx=64
	ny=64
	dx=lx/float(nx) ; dy=ly/float(ny)

	nx1=1 ; nx2=nx ; ny1=1 ; ny2=ny !Definición de las fronteras

	Re=1.0
	Ga=0.7
	gma=0.1

	allocate(x(0:nx+1),y(0:ny+1),x1(0:nx+1),y1(0:ny+1))
	allocate(u1(0:nx+1,0:ny+1),u2(0:nx+1,0:ny+1),v1(0:nx+1,0:ny+1),v2(0:nx+1,0:ny+1),&
					 F(0:nx+1,0:ny+1),G(0:nx+1,0:ny+1),p1(0:nx+1,0:ny+1),p2(0:nx+1,0:ny+1),&
					 flx(0:nx+1,0:ny+1),fly(0:nx+1,0:ny+1),psi(0:nx+1,0:ny+1),zeta(0:nx+1,0:ny+1))
	
	allocate(du2dx(0:nx+1,0:ny+1),duvdy(0:nx+1,0:ny+1),&
					 d2udx2(0:nx+1,0:ny+1),d2udy2(0:nx+1,0:ny+1))
	allocate(dv2dy(0:nx+1,0:ny+1),dvudx(0:nx+1,0:ny+1),&
					 d2vdx2(0:nx+1,0:ny+1),d2vdy2(0:nx+1,0:ny+1), &
					 B0z(0:nx+1,0:ny+1), Bz(0:nx+1,0:ny+1))


! Construction of the numerical grid
	xc = lx/2.0
	x(nx1-1)=0.0
	x(nx1)=0.5*dx
	x1(nx1-1) = -xc
	x1(nx1) = -xc + 0.5*dx
	do i=(nx1+1),nx2
		x(i)=0.5*dx+(i-nx1)*dx
		x1(i) = x(i) -xc
	enddo
	x(nx2+1)=lx
	x1(nx2+1) = lx/2

	yc = ly/2.0
	y(ny1-1)=0.0
	y(ny1)=0.5*dy
	y1(ny1-1) = -yc
	y1(ny1) = -yc +0.5*dy
	do j=(ny1+1),ny2
		y(j)=0.5*dy+(j-ny1)*dy
		y1(j) = y(j) - yc
	enddo
	y(ny2+1)=ly
	y1(ny2+1) = ly/2

! Temporal step size definition
	dt=gma*(Re/2.)*( ( (1./dx)**2 + (1./dy)**2 )**(-1) )
	print *, 'step size dt=',dt
	tmax=50000 ! Number of steps during the simulation
	nplot=500 ! animation parameter

!======================
!======================
! Calling magnetic fiel subroutine
call MagField(x,y,x1,y1,nx,ny,nx1,nx2,ny1,ny2,Bz,B0z,lx,ly)


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

	F(i,j)=u1(i,j)+dt*((d2udx2(i,j) + d2udy2(i,j))&
		 - du2dx(i,j) - duvdy(i,j) + Re*flx(i,j) )
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

	fly(i,j)=-0.5*(Bz(i,j+1)+Bz(i,j))

	G(i,j)=v1(i,j)+dt*( (d2vdx2(i,j)+d2vdy2(i,j))&
		 -dvudx(i,j)-dv2dy(i,j)+Re*fly(i,j) )
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
! ------------------------------------------------------------
! Pressure subroutine ----------------------------------------
subroutine PRESSURE(p1,p2,F,G,nx,ny,dx,dy,dt)
implicit none

	integer, value :: nx,ny
	real, value :: dx,dy,dt
	real, dimension(0:nx+1,0:ny+1),intent(inout) :: p1,p2,F,G

	integer :: pt,pj,i,j
	real :: pdt

	pdt=0.000001
	pt=100

	DO pj=1, pt

	forall(i=1:nx,j=1:ny)
		p2(i,j)=p1(i,j)+pdt*( (p1(i+1,j)-2.0*p1(i,j)+p1(i-1,j))/(dx*dx) &
			+ (p1(i,j+1)-2.0*p1(i,j)+p1(i,j-1))/(dy*dy)&
			- (1.0/dt)*( (F(i,j)-F(i-1,j))/(dx) + (G(i,j)-G(i,j-1))/(dy) ) )
	end forall

	do j=0,ny+1
		p1(0,j)=p1(1,j)
		p2(0,j)=p2(1,j)
		p1(nx+1,j)=p1(nx,j)
		p2(nx+1,j)=p2(nx,j)
	enddo

	do i=0,nx+1
		p1(i,0)=p1(i,1)
		p2(i,0)=p2(i,1)
		p1(i,ny+1)=p1(i,ny)
		p2(i,ny+1)=p2(i,ny)
	enddo

	p1=p2
	END DO

return
end subroutine

! Boundary conditions subroutine
! This subroutine contains four boundary conditions for each boundary:
! No slip, free slip, outflow, and inflow.

subroutine boundary_conditions(u1,u2,v1,v2,nx,ny)
implicit none

	integer, value :: nx,ny
	real, dimension(0:nx+1,0:ny+1), intent(inout) :: u1,u2,v1,v2
	integer :: i,j

! x=0,y - left
	do j=0,ny
		u1(0,j)=0.0
		u2(0,j)=0.0
		v1(0,j)=-v1(1,j)
		v2(0,j)=-v2(1,j)
	end do

  ! x=nx,y - right
	do j=0,ny+1
		u1(nx,j)=0.0
		u2(nx,j)=0.0
		v1(nx+1,j)=-v1(nx,j)
		v2(nx+1,j)=-v2(nx,j)
	enddo

	!x,y=ny - top
	do i=0,nx+1
		u1(i,ny+1)=-u1(i,ny)
		u2(i,ny+1)=-u2(i,ny)
		v1(i,ny)=0.0
		v2(i,ny)=0.0
	enddo

	!x,y=0 - bottom
	do i=0,nx+1
		u1(i,0)=-u1(i,1) 		!-u1(i,1) Condiciones para pared inferior fija (NO-SLIP)
		u2(i,0)=-u2(i,1) 		!-u2(i,1)
		v1(i,0)=0.0!-v1(i,0) 	!0.0
		v2(i,0)=0.0!-v2(i,0) 	!0.0
	enddo

return
end subroutine

subroutine animation(x,y,u2,v2,p2,nx,ny,nplot,t,lx,ly,psi,zeta)
implicit none

	integer :: i,j
	integer,intent(inout) :: nx,ny,nplot,t
	real,intent(inout) :: lx,ly
	real,intent(inout) :: x(0:nx+1),y(0:ny+1)
	real,dimension(0:nx+1,0:ny+1),intent(inout) :: u2,v2,p2
  real,dimension(0:nx+1,0:ny+1),intent(inout) :: psi,zeta

	integer :: ncount
	character :: ext*4,destiny*512,nctext*4,outfile*512,fname*16
	character :: anim_file*512,anim_file_name*16

	ext='.dat'
	destiny='./anim/'
	fname='frame_'
	ncount=t/nplot
	write(nctext,'(I4.4)') ncount

	outfile=trim(destiny)//trim(fname)//trim(nctext)//ext

	open(20,file=outfile)
	200 format(E12.4,1x,E12.4,1x,E12.4,1x,E12.4,1x,E12.4,1x,E12.4,1x,E12.4,1x,E12.4)

	do i=1,nx
	do j=1,ny
		write(20,200) x(i),y(j),u2(i,j),v2(i,j),p2(i,j),&
									sqrt(u2(i,j)*u2(i,j) + v2(i,j)*v2(i,j)),&
                  psi(i,j),zeta(i,j)
	end do
	end do

	anim_file_name='anim.gp'
	anim_file=trim(destiny)//trim(anim_file_name)
	open(30,file=anim_file)
	write(30,*) "p '"//trim(fname)//trim(nctext)//ext//"' every e u 1:2:(f*$3):(f*$4) w vec filled lt -1 "
	write(30,*) "pause 0.05"

return
end subroutine animation

! ------------------------------------------------------------
! Magnetic Flied subroutine ------------------------
subroutine MagField(x,y,x1,y1,nx,ny,nx1,nx2,ny1,ny2,Bz,B0z,lx,ly)

implicit none
integer :: i,j,n,m,k
integer, intent(inout) :: nx,ny,nx1,nx2,ny1,ny2
real, intent(inout) :: lx,ly
real,dimension(:), intent(inout) :: x(0:nx+1), y(0:ny+1),x1(0:nx+1), y1(0:ny+1)
real, dimension(0:nx+1,0:ny+1), intent(inout) :: Bz,B0z
real :: di,dd,a,b,c,d,E1,z,g
real, dimension(:), allocatable :: xn,ym,zk

allocate(xn(2),ym(2),zk(2))

! Parameters of the magnet
  di = 0.015
  a = 0.5*di
  b = 0.5*di
  c = 0.5*di

  xn(1) = -a; xn(2) = a
  ym(1) = -b; ym(2) = b
  zk(1) = -c; zk(2) = c

  d = 6.e-3

  E1 = 2.9e-2

  Bz = 0.0; B0z = 0.0

  do n = 1,2
  do m = 1,2
  do k = 1,2

  do i=nx1-1,nx2+1
  do j=ny1-1,ny2+1

  z = c+d
  g = 1.0/sqrt((di*x1(i)-xn(n))**2+(di*y1(j)-ym(m))**2+(z-zk(k))**2)
  B0z(i,j) = ((-1.0)**(n+m+k))*atan(((di*x1(i)-xn(n))*(di*y1(j)-ym(m)))/(z-zk(k))*g)

  end do
  end do

  Bz = Bz+B0z

  end do
  end do
  end do

  Bz = E1*Bz
  Bz=Bz/maxval(Bz)
  
  forall(i=(nx1-1):(nx2+1),j=(ny1-11):(ny2+1))
  	Bz(i,j) = Bz(i,j)
  end forall
  
! Magnetic flied profile
  do i=nx1-1,nx2+1
  do j=ny1-1,ny2+1
    open(1,file='B0z_profile_one.dat')
    write(1,*) x(i)-lx/2, Bz(i,ny/2)
  end do
  end do
  close(1)
  
! Magnetic field map
  2 format(1x,E12.5,2x,E12.5,2x,E12.5)
  open(2,file='B0z_field.dat')
  do j=ny1,ny2
    write(2,2)
  do i=nx1,nx2
    write(2,2) (x(i)-lx/2.0), (y(j)-ly/2.0), Bz(i,j)
   !	 write(1,2) x1(i)*dd, y1(j)*dd, Bz(i,j)
  end do
  end do
  close(2)


return 
end subroutine MagField

! ------------------------------------------------------------
! Output subroutine ------------------------------------------
subroutine output(x,y,nx,ny,nx1,nx2,ny1,ny2,u1,u2,v1,v2,p1,p2,&
                  psi,zeta,lx,ly)
implicit none

  integer :: i,j
  integer,intent(inout) :: nx,ny,nx1,nx2,ny1,ny2
  real,intent(inout) :: lx,ly
  real,dimension(:,:),intent(inout) :: x(0:nx+1),y(0:ny+1)
  real,dimension(0:nx+1,0:ny+1),intent(inout) :: u1,u2,v1,v2,p1,p2
  real,dimension(0:nx+1,0:ny+1),intent(inout) :: psi,zeta

  real :: mag_u ! velocity magnitude

! 1. Vector velocity field ============

  open(10,file='vec_field.dat')
    do i=nx1-1,nx2+1
    write(10,100)
    do j=ny1-1,ny2+1
      mag_u = sqrt(u1(i,j)*u1(i,j) + v1(i,j)*v1(i,j))
      write(10,100) x(i)-(lx/2.0),y(j)-(ly/2.0),u1(i,j),v1(i,j),p2(i,j),psi(i,j),mag_u,zeta(i,j)
    end do
    end do
  close(10)

  100 format(2x,E12.4,1x,E12.4,1x,E12.4,1x,E12.4,1x,E12.4,1x,E12.4,1x,E12.4,1x,E12.4)

return
end subroutine output
