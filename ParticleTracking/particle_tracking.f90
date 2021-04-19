! Lagrangian particle tacking program. This code reads four columns from the files
! and performs a particle tracking subroutine based on the integration of the
! position of the components over time.

program ptracking
implicit none

  integer :: i,j,k,nx,ny,n,m,nfiles,tmax
  real :: lx,ly,x0,y0,u,v,press,magu,psi,zeta,dx,dy,dt,lm
  real :: thrd,rrd,Up,Vp,x10,x20,y10,y20,pi
  real :: gma,Re,nplot

  real, dimension(:,:), allocatable :: u1,v1
  real, dimension(:), allocatable :: x1,y1,x2,y2

  character :: ext*4,fname*16,ifile*1024,nctext*4,flink*1024

  Re = 3.0; gma = 0.1
  nplot = 50.0
  nx = 100
  ny = 100
  lx = 10.0
  ly = 10.0
  dx = lx/float(nx); dy = ly/float(ny)
  dt = nplot*(gma*(Re/2.)*( ( (1./dx)**2 + (1./dy)**2 )**(-1) ))
  tmax = 50000
  nfiles = tmax/int(nplot)
  write(*,*) 'time step dt= ',dt
  write(*,*) 'total number of files = ',nfiles

  n = 10000
  lm = 1.0
  pi = 3.1415926535

  ext='.dat'
  fname='frame_'
  ifile='./VecFieldAnim/'

  allocate(u1(1:nx,1:ny),v1(1:nx,1:ny))
  allocate(x1(1:n),y1(1:n),x2(1:n),y2(1:n))

! Creating particle circle at the center
  do m=1,n
    call RANDOM_NUMBER(thrd)
    call RANDOM_NUMBER(rrd)
    x1(m)=(rrd*lm)*cos(2.0*pi*thrd)+5.0
    y1(m)=(rrd*lm)*sin(2.0*pi*thrd)+5.0
  end do

!-----------------------------------------------
! Reading files --------------------------------
DO k=1,nfiles

  write(nctext,'(I4.4)') k
  flink=trim(ifile)//trim(fname)//nctext//ext
  open(1,file=flink)
  print*, 'Loaded: ', trim(fname)//nctext//ext

  do i=1,nx
  do j=1,ny
    read(1,*) x0,y0,u,v
      u1(i,j)=u
      v1(i,j)=v
  end do
  end do

  Do m=1,n

! ----------------------------------------------
! dx/dt = u ------------------------------------

  i = int(x1(m)/dx)+1 ; j = int((y1(m)+dy/2.0)/dy)+1

  x10=(i-1)*dx ; y10=((j-1)-0.5)*dy
  x20=i*dx     ; y20=(j-0.5)*dy

  Up=(1.0/(dx*dy))*((x20-x1(m))*(y20-y1(m))*u1(i-1,j-1)&
                  + (x1(m)-x10)*(y20-y1(m))*u1(i,j-1)&
                  + (x20-x1(m))*(y1(m)-y10)*u1(i-1,j)&
                  + (x1(m)-x10)*(y1(m)-y10)*u1(i,j) )

  x2(m) = x1(m) + Up*dt
  if((x2(m).lt.0).or.(x2(m).gt.lx)) x2(m)=-1.0

! dy/dt = v ------------------------------------

  i = int((x1(m)+dx/2.0)/dx)+1; j = int(y1(m)/dy)+1

  x10=((i-1)-0.5)*dx ; y10=(j-1)*dy
  x20=(i-0.5)*dx     ; y20=j*dy

  Vp=(1.0/(dx*dy))*((x20-x1(m))*(y20-y1(m))*v1(i-1,j-1)&
                  + (x1(m)-x10)*(y20-y1(m))*v1(i,j-1)&
                  + (x20-x1(m))*(y1(m)-y10)*v1(i-1,j)&
                  + (x1(m)-x10)*(y1(m)-y10)*v1(i,j) )

  y2(m) = y1(m) + Vp*dt
  if((y2(m).lt.0).or.(y2(m).gt.ly)) y2(m)=0.0

  i=int(x2(m)/dx)+1; j=int(y2(m)/dy)+1

  End do

  call animation(n,x2,y2,k)

  x1=x2
  y1=y2

END DO

end program ptracking

! Animation subroutine ----------------------------------

subroutine animation(n,x2,y2,k)
implicit none

  integer,intent(inout) :: n,k
  real,dimension(1:n),intent(inout) :: x2,y2
  integer :: m

  integer :: ncount
  character :: ext*4,fname*16,destiny*512,nctext*4,outfile*512
  character :: anim_file*512,anim_file_name*16

  ext='.dat'
  destiny='./PTAnim/'
  fname='ptframe_'
  ncount=k
  write(nctext,'(I4.4)') ncount

  outfile=trim(destiny)//trim(fname)//trim(nctext)//ext

  1 format(2x,F8.5,1x,F8.5) ! gnuplot
!  2 format(I0,',',F8.5,',',F8.5) ! paraview
  open(10,file=outfile)
!  write(10,*) 'PART,XPOS,YPOS'
  do m=1,n
    write(10,1) x2(m),y2(m)
  end do

  anim_file_name='anim.gp'
	anim_file=trim(destiny)//trim(anim_file_name)
	open(30,file=anim_file)
	write(30,*) "p '"//trim(fname)//trim(nctext)//ext//"' w p pt 7 ps 0.5 lt -1 t 't="//trim(nctext)//"'"
	write(30,*) "pause 0.05"

return
end subroutine animation
