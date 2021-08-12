program furlani_magnetic_field

implicit none
integer :: i,j
integer :: nx, ny, nx1, nx2, ny1, ny2
integer :: n,m,k
real :: lx,ly,dx,dy,xc,yc,di
real :: a,b,c,d,g,z,E1
real :: dlx1,dlx2,dly1,dly2,dlx,dly,dd
real, dimension(:), allocatable :: x,x1,y,y1
real, dimension(:,:), allocatable :: Bz,B0z
real, dimension(:), allocatable :: xn,ym,zk

dd = 0.00072

dlx = 0.02
dly = 0.02
lx = real(nint(dlx/dd))
ly = real(nint(dly/dd))
nx = (nint(lx))*5
ny = (nint(ly))*5
dx = lx/float(nx)
dy = ly/float(ny)

nx1 = 1 ; nx2 = nx
ny1 = 1 ; ny2 = ny

allocate(x(0:nx+1),y(0:ny+1),x1(0:nx+1),y1(0:ny+1))
allocate(B0z(0:nx+1,0:ny+1),Bz(0:nx+1,0:ny+1))
allocate(xn(2),ym(2),zk(2))

!Construction of the grid
xc = lx/2.0
x(nx1-1) = 0.0
x(nx1) = 0.5*dx
x1(nx1-1) = -xc
x1(nx1) = -xc +0.5*dx

do i = (nx1+1), nx2
   x(i) = 0.5*dx + (i-nx1)*dx
   x1(i) = x(i)-xc
enddo
x(nx2+1) = lx
x1(nx2+1) = lx/2

yc = ly/2.0
y(ny1-1) = 0.0
y(ny1) = 0.5*dy
y1(ny1-1) = -yc
y1(ny1) = -yc +0.5*dy

do j = (ny1+1), ny2
   y(j) = 0.5*dy + (j-ny1)*dy
   y1(j) = y(j)-yc
enddo
y(ny2+1) = ly
y1(ny2+1) = ly/2

!Parameters of the magnet
di = 1.0e-3
a = 0.5*di
b = 0.5*di
c = 0.5*di

xn(1) = -a; xn(2) = a
ym(1) = -b; ym(2) = b
zk(1) = -c; zk(2) = c

d = 6.e-3
E1 = 2.9e-2
Bz = 0.0 ; B0z = 0.0

do n = 1,2
do m = 1,2
do k = 1,2

do i = nx1-1, nx2+1
do j = ny1-1, ny2+1

   z = c+d
   g = 1.0/sqrt((di*x1(i)-xn(n))**2+(di*y1(j)-ym(m))**2+(z-zk(k))**2)
   B0z(i,j) = ((-1.0)**(n+m+k))*atan(((di*x1(i)-xn(n))*(di*y1(j)-ym(m)))/(z-zk(k))*g)

enddo
enddo
Bz = Bz+B0z
enddo
enddo
enddo

Bz = E1*Bz
Bz = Bz/maxval(Bz)

do i=nx1-1, nx2+1
do j=ny1-1, ny2+1
   open(20, file='B0z_profile_one.dat')
   write(20,*) x(i)-lx/2, Bz(i,ny/2)
enddo
enddo

2 format(1x,E12.5,2x,E12.5,2x,E12.5)
open(1, file='B0z_field.dat')
do j=ny1,ny2
  write(1,2)
  do i=nx1,nx2
     write(1,2) (x(i)-lx/2.0), (y(j)-ly/2.0), Bz(i,j)
  enddo
enddo
close(1)

end program furlani_magnetic_field 
















