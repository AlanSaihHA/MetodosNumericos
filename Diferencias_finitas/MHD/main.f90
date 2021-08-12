! Program written to simulate a 2D MHD dipole

!include 'animation.f90' ! animation subroutine
!include 'mag_fields.f90' ! Magnetic field generation subroutine

program LD2DNS
	implicit none

	integer :: i,j,nx,ny,nx1,nx2,ny1,ny2,t,tmax,nplot
	real :: dx,dy,lx,ly,dt,Re,Ga,Delta,gma

  integer :: mag_field,magnet_no

	real, dimension(:), allocatable :: x,y
	real, dimension(:,:), allocatable :: u1,u2,v1,v2,p1,p2,F,G

	real, dimension(:,:), allocatable :: du2dx,duvdy,d2udx2,d2udy2
	real, dimension(:,:), allocatable :: dv2dy,dvudx,d2vdx2,d2vdy2

	real, dimension(:,:), allocatable :: flx,fly

	real, dimension(:), allocatable :: x1,y1
	real, dimension(:,:), allocatable :: psi,zeta,sigma,lambda
  real, dimension(:,:), allocatable :: B0z

	lx=10.0 ; ly=10.0 ! Dimensionless size of the domain
	nx=100
	ny=100
	dx=lx/float(nx) ; dy=ly/float(ny)

	nx1=1 ; nx2=nx ; ny1=1 ; ny2=ny !Definición de las fronteras

	Re=3.0
	Ga=0.7
	gma=0.1

	allocate(x(0:nx+1),y(0:ny+1))
	allocate(u1(0:nx+1,0:ny+1),u2(0:nx+1,0:ny+1),v1(0:nx+1,0:ny+1),v2(0:nx+1,0:ny+1))
	allocate(F(0:nx+1,0:ny+1),G(0:nx+1,0:ny+1))
	allocate(p1(0:nx+1,0:ny+1),p2(0:nx+1,0:ny+1))

	allocate(du2dx(0:nx+1,0:ny+1),duvdy(0:nx+1,0:ny+1),d2udx2(0:nx+1,0:ny+1),d2udy2(0:nx+1,0:ny+1))
	allocate(dv2dy(0:nx+1,0:ny+1),dvudx(0:nx+1,0:ny+1),d2vdx2(0:nx+1,0:ny+1),d2vdy2(0:nx+1,0:ny+1))

	allocate(flx(0:nx+1,0:ny+1),fly(0:nx+1,0:ny+1))

	allocate(psi(0:nx+1,0:ny+1),zeta(0:nx+1,0:ny+1),sigma(0:nx+1,0:ny+1),lambda(0:nx+1,0:ny+1))
  allocate(B0z(0:nx+1,0:ny+1))

! ------------------------------------
! Parameters of the magnetic field  --
! ------------------------------------

  ! Field by: McCaig = 1, Furlani = 2
  mag_field = 1

  ! Number of magnets
  magnet_no = 1

! ------------------------------------
! ------------------------------------

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
	tmax=50000 ! Number of steps during the simulation
	nplot=50 ! animation parameter

! Calling magnetic field subroutine

  call magnetic_field(nx,ny,nx1,nx2,ny1,ny2,dx,dy,lx,ly,magnet_no,&
                            mag_field,B0z)
  B0z=B0z/maxval(B0z)

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

	fly(i,j)=(B0z(i,j+1)+B0z(i,j))

	G(i,j)=v1(i,j)+dt*( (1.0/Re)*(d2vdx2(i,j)+d2vdy2(i,j))&
		 -dvudx(i,j)-dv2dy(i,j) + Re*fly(i,j) )
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
!		do i=0,nx
!			psi(i,0) = 0.0
!			do j=1,ny-1
!				psi(i,j)=psi(i,j-1) + u2(i,j)*dy
!			end do
!		end do

! Okubo-Weiss (O-K) parameter
! 	computes the vorticity, strain field and the O-K
!    do i=1,nx
!      do j=1,ny
!        zeta(i,j) = (1.0/dy)*(u1(i,j+1)-u1(i,j)) - (1.0/dx)*(v1(i+1,j)-v1(i,j))
!				sigma(i,j) = (1.0/dy)*(u1(i,j+1)-u1(i,j)) + (1.0/dx)*(v1(i+1,j)-v1(i,j))
!				lambda(i,j) = (1.0/4.0)*(zeta(i,j)*zeta(i,j) - sigma(i,j)*sigma(i,j))
!      end do
!    end do

! Printing values on terminal
	if(mod(t,nplot)==0) THEN
		Delta=MAXVAL(ABS(u2-u1))
		write(*,*) 'T=',t,'U1=', u1(nx/2,ny/2),'V1=', v1(nx/2,ny/2), 'E=', Delta
		call animation(x,y,u2,v2,p2,nx,ny,nplot,t,lx,ly,psi,zeta,sigma,lambda)
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
	            psi,zeta,sigma,lambda,lx,ly,magnet_no)

! Main program ends =====
end program LD2DNS !=====
!========================

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
		u1(i,ny+1)=-u1(i,ny)!1.0
		u2(i,ny+1)=-u1(i,ny)!1.0
		v1(i,ny)=0.0!-v1(i,ny+1)
		v2(i,ny)=0.0!-v2(i,ny+1)
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

! Magnetic field subroutines:
! mag_field = {1,2}
! a) McCaig
! b) Furlani

subroutine magnetic_field(nx,ny,nx1,nx2,ny1,ny2,dx,dy,lx,ly,magnet_no,&
                          mag_field,B0z)
implicit none

!--------------------
! Global variables --
!--------------------

  integer :: i,j
  integer,intent(inout) :: nx,ny,nx1,nx2,ny1,ny2
  integer,intent(inout) :: magnet_no,mag_field
  real,intent(inout) :: lx,ly,dx,dy
  character :: magnets*16, outfile*64, file_name*64, ext*4

  ! Variables used to organize an array of any number of magnets
  ! li and lj correspond to the displacement distance in x and y
  ! ci and cj correspond to the implemented displacement in terms
  !   of the numerical grid
  ! remember that the influence of the magnetic field varies
  ! with the size of the magnets, in this case, with the size
  ! of a and b. Be careful to let the field fit inside the domain,
  ! or there will be trouble. This can be done by restricting the
  ! influence with "ifs" or by resizing either the domain or a and b.

  real :: li1,li2,li3,&
          lj1,lj2,lj3
  integer :: ci1,ci2,ci3,&
             cj1,cj2,cj3

!-------------------------------------
! McCaig's magnetic field variables --
!-------------------------------------

  real :: E1, xc, yc, z0, z01, z02, Bmax, a, b
  real :: z0s, z01s, z02s
  real, dimension(:), allocatable :: x, x1, y, y1
  real, dimension(:,:), allocatable :: B0z1, B0z2
  real, intent(inout) :: B0z(0:nx+1,0:ny+1)

!---------------------------------------
! Furlanis's magnetic field variables --
!---------------------------------------

! %%%%%%%%%%%%%%%

!---------------------
! Vector allocation --
!---------------------

  allocate(x(0:nx+1),y(0:ny+1),x1(0:nx+1),y1(0:ny+1))
  allocate(B0z1(0:nx+1,0:ny+1),B0z2(0:nx+1,0:ny+1))

! ---------------------------------------------
! ---------------------------------------------

!-------------------
! Grid generation --
!-------------------

  xc = lx/2.0
  x(nx1-1) = 0.0
  x(nx1) = 0.5*dx
  x1(nx1-1) = -xc
  x1(nx1) = -xc + 0.5*dx
  do i = (nx1+1), nx2
    x(i) = 0.5*dx + (i-nx1)*dx
    x1(i) = x(i) - xc
  end do
  x(nx2+1) = lx
  x1(nx2+1) = lx/2

  yc = ly/2.0
  y(ny1-1) = 0.0
  y(ny1) = 0.5*dy
  y1(ny1-1) = -yc
  y1(ny1) = -yc + 0.5*dy
  do j = (ny1+1), ny2
    y(j) = 0.5*dy + (j-ny1)*dy
    y1(j) = y(j) - yc
  end do
  y(ny2+1) = ly
  y1(ny2+1) = ly/2

! ---------------------------------------------------
! McCaig --------------------------------------------
if(mag_field==1) then

  a = 0.5
  b = 0.5 ! a and b represent the length of the sides of the magnet
  z0   = 0.00005 ! Height of the top magnetic surface
  z01  = 0.5   ! Height where the first magnetic field is taken
  z02  = 1.0   ! Height where the second magnetic field is taken
  Bmax = 4.5   ! Maximum magnetic field value

  z0s  = z0**2
  z01s = z01**2
  z02s = z02**2

! --------------------------
! Normalization parameter --
! --------------------------

  E1 = Bmax/(atan(((0.0+a)*(0.0+b))/(z0*sqrt((0.0+a)**2+(0.0+b)**2+z0s)))&
            +atan(((0.0-a)*(0.0-b))/(z0*sqrt((0.0-a)**2+(0.0-b)**2+z0s)))&
            -atan(((0.0+a)*(0.0-b))/(z0*sqrt((0.0+a)**2+(0.0-b)**2+z0s)))&
            -atan(((0.0-a)*(0.0+b))/(z0*sqrt((0.0-a)**2+(0.0+b)**2+z0s)))&
            )
  write(*,*) 'E1 =', E1

  forall(i=(nx1-1):(nx2+1),j=(ny1-1):(ny2+1))
    B0z1(i,j)=E1*&
              (atan(((x1(i)+a)*(y1(j)+b))/(z01*sqrt((x1(i)+a)**2+(y1(j)+b)**2+z01s)))&
              +atan(((x1(i)-a)*(y1(j)-b))/(z01*sqrt((x1(i)-a)**2+(y1(j)-b)**2+z01s)))&
              -atan(((x1(i)+a)*(y1(j)-b))/(z01*sqrt((x1(i)+a)**2+(y1(j)-b)**2+z01s)))&
              -atan(((x1(i)-a)*(y1(j)+b))/(z01*sqrt((x1(i)-a)**2+(y1(j)+b)**2+z01s)))&
              )
    B0z2(i,j)=E1*&
              (atan(((x1(i)+a)*(y1(j)+b))/(z02*sqrt((x1(i)+a)**2+(y1(j)+b)**2+z02s)))&
              +atan(((x1(i)-a)*(y1(j)-b))/(z02*sqrt((x1(i)-a)**2+(y1(j)-b)**2+z02s)))&
              -atan(((x1(i)+a)*(y1(j)-b))/(z02*sqrt((x1(i)+a)**2+(y1(j)-b)**2+z02s)))&
              -atan(((x1(i)-a)*(y1(j)+b))/(z02*sqrt((x1(i)-a)**2+(y1(j)+b)**2+z02s)))&
              )
  end forall

! Magnetic field calculation ----------------

! Magnetic field generated by 1 magnet

  if(magnet_no.eq.1) then
  write(magnets,'(I0.0)') magnet_no

  forall(i=(nx1-1):(nx2+1),j=(ny1-1):(ny2+1))
    B0z(i,j) = B0z1(i,j) - B0z2(i,j)
  end forall

! Magnetic field generated by 2 magnets
  else if(magnet_no.eq.2) then
  write(magnets,'(I0.0)') magnet_no

    li1 = 0.1
    ci1 = nint((nx*li1)/lx)

  forall(i=(nx1-1):(nx2+1),j=(ny1-1):(ny2+1))
    B0z(i,j) = (B0z1(i+ci1,j) - B0z2(i+ci1,j))&
              -(B0z1(i-ci1,j) - B0z2(i-ci1,j))
  end forall

! Magnetic field generated by 4 magnets

  else if(magnet_no.eq.4) then
  write(magnets,'(I0.0)') magnet_no

    li1 = 1.0
    ci1 = nint((nx*li1)/lx)

    lj1 = 1.0
    cj1 = nint((ny*lj1)/ly)

    do i=(nx1-1),(nx2+1)
    do j=(ny1-1),(ny2+1)
      if((i+ci1).lt.(nx1+1).or.(i+ci1).gt.(nx2-1)) then
        B0z(i,j)=0.0
      else if((j+cj1).lt.(ny1+1).or.(j+cj1).gt.(ny2-1)) then
        B0z(i,j)=0.0
      else if((i-ci1).lt.(nx1+1).or.(i-ci1).gt.(nx2-1)) then
        B0z(i,j)=0.0
      else if((j-cj1).lt.(ny1+1).or.(j-cj1).gt.(ny2-1)) then
        B0z(i,j)=0.0
      else
        B0z(i,j) = (B0z1(i-ci1,j-cj1) - B0z2(i-ci1,j-cj1))&
                  -(B0z1(i+ci1,j-cj1) - B0z2(i+ci1,j-cj1))&
                  +(B0z1(i+ci1,j+cj1) - B0z2(i+ci1,j+cj1))&
                  -(B0z1(i-ci1,j+cj1) - B0z2(i-ci1,j+cj1))
      end if
    end do
    end do

  end if
! ---------------------------------------------------
! ---------------------------------------------------

! ---------------------------------------------------
! Furlani -------------------------------------------
else if(mag_field==2) then
  B0z(i,j)=0.0
end if

! Magnetic field file writing -----------------------

  file_name='B0z_field_'
  ext = '.dat'
  outfile=trim(file_name)//trim(magnets)//ext
  open(1,file=outfile)
  2 format(1x,E12.5,2x,E12.5,2x,E12.5)
  do i=nx1,nx2
    write(1,2)
  do j=ny1,ny2
    write(1,2) x1(i), y1(j), B0z(i,j)
  end do
  end do
  close(1)


end subroutine magnetic_field

!--------------------------------------------------------
!--------------------------------------------------------

subroutine animation(x,y,u2,v2,p2,nx,ny,nplot,t,lx,ly,psi,zeta,sigma,lambda)
implicit none

	integer :: i,j
	integer,intent(inout) :: nx,ny,nplot,t
	real,intent(inout) :: lx,ly
	real,intent(inout) :: x(0:nx+1),y(0:ny+1)
	real,dimension(0:nx+1,0:ny+1),intent(inout) :: u2,v2,p2
  real,dimension(0:nx+1,0:ny+1),intent(inout) :: psi,zeta,sigma,lambda

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
	200 format(E12.4,1x,E12.4,1x,E12.4,1x,E12.4,1x,E12.4,1x,E12.4,1x,E12.4,1x,E12.4,1x,E12.4,1x,E12.4)

	do i=1,nx
	do j=1,ny
		write(20,200) x(i),y(j),u2(i,j),v2(i,j)
	end do
	end do

	anim_file_name='anim.gp'
	anim_file=trim(destiny)//trim(anim_file_name)
	open(30,file=anim_file)
	write(30,*) "p '"//trim(fname)//trim(nctext)//ext//"' every e u 1:2:(f*$3):(f*$4) w vec filled lt -1"
	write(30,*) "pause 0.05"

return
end subroutine animation

! ------------------------------------------------------------
! Output subroutine ------------------------------------------
subroutine output(x,y,nx,ny,nx1,nx2,ny1,ny2,u1,u2,v1,v2,p1,p2,&
                  psi,zeta,sigma,lambda,lx,ly,magnet_no)
implicit none

  integer :: i,j
  integer,intent(inout) :: nx,ny,nx1,nx2,ny1,ny2
	integer,intent(inout) :: magnet_no
  real,intent(inout) :: lx,ly
  real,dimension(:,:),intent(inout) :: x(0:nx+1),y(0:ny+1)
  real,dimension(0:nx+1,0:ny+1),intent(inout) :: u1,u2,v1,v2,p1,p2
  real,dimension(0:nx+1,0:ny+1),intent(inout) :: psi,zeta,sigma,lambda

  real :: mag_u ! velocity magnitude

	character :: magnets*16,outfile*64,filename*64,ext*4

! 1. Vector velocity field ============

	write(magnets,'(I0.0)') magnet_no
	ext='.dat'
	filename='vec_field_'
	outfile=trim(filename)//trim(magnets)//ext

  open(10,file=outfile)
    do i=nx1-1,nx2+1
    do j=ny1-1,ny2+1
      mag_u = sqrt(u1(i,j)*u1(i,j) + v1(i,j)*v1(i,j))
      write(10,100) x(i)-(lx/2.0),y(j)-(ly/2.0),u1(i,j),v1(i,j),p2(i,j),&
										psi(i,j),mag_u,zeta(i,j),sigma(i,j),lambda(i,j)
    end do
    end do
  close(10)

	write(*,*) 'NX*NY=',nx,'U1=', u1(nx/2,ny/2),'V1=', v1(nx/2,ny/2)

  100 format(2x,E12.4,1x,E12.4,1x,E12.4,1x,E12.4,1x,E12.4,&
						 1x,E12.4,1x,E12.4,1x,E12.4,1x,E12.4,1x,E12.4)

return
end subroutine output
