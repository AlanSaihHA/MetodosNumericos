! program written to replicate an experimental magnetic field distribution,
! which is calculated based on the anaylitcal expression from
! Permanent Magnets in Theory and Practice, McCaig 1977, pp. 180-189.
! Jose Olvera
! 22/July/2020

program mccaig_magnetic_field

  implicit  none

  integer :: i, j
  integer :: nx, ny, nx1, nx2, ny1, ny2
  integer :: magnet_no
  real :: lx, ly, dx, dy
  real :: E1, xc, yc, z0, z01, z02, Bmax, a, b
  real :: z0s, z01s, z02s
  real, dimension(:), allocatable :: x, x1, y, y1
  real, dimension(:,:), allocatable :: B0z, B0z1, B0z2

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

! Parameters of the domain and the grid
  lx = 5.0; ly = 5.0
  nx = 128; ny = 128
  dx = lx/float(nx); dy = ly/float(ny)

  nx1 = 1; nx2 = nx
  ny1 = 1; ny2 = ny

  allocate(x(0:nx+1),y(0:ny+1),x1(0:nx+1),y1(0:ny+1))
  allocate(B0z(0:nx+1,0:ny+1),B0z1(0:nx+1,0:ny+1),B0z2(0:nx+1,0:ny+1))

! Construction of the grid

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

! Parameters of the magnetic field

  magnet_no = 1
  a = 0.5
  b = 0.5 ! a and b represent the length of the sides of the magnet
  z0   = 0.00005 ! Height of the top magnetic surface
  z01  = 0.5   ! Height where the first magnetic field is taken
  z02  = 1.0   ! Height where the second magnetic field is taken
  Bmax = 4.5   ! Maximum magnetic field value

  z0s  = z0**2
  z01s = z01**2
  z02s = z02**2

! Normalization parameter

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

! Magnetic field generated by 1 magnet
  if(magnet_no.eq.1) then
  write(magnets,'(I0.0)') magnet_no
  forall(i=(nx1-1):(nx2+1),j=(ny1-1):(ny2+1))
    B0z(i,j) = B0z1(i,j) - B0z2(i,j)
  end forall

! Writing file with the magnetic field profile computed
! along all x at y=ny/2.

  do i=nx1-1,nx2+1
  do j=ny1-1,ny2+1
    open(20,file='B0z_profile_one.dat')
    write(20,*) x(i)-lx/2, B0z(i,ny/2)
  end do
  end do

! Magnetic field generated by 2 magnets
  else if(magnet_no.eq.2) then
  write(magnets,'(I0.0)') magnet_no

    li1 = 0.5
    ci1 = nint((nx*li1)/lx)

  forall(i=(nx1-1):(nx2+1),j=(ny1-1):(ny2+1))
    B0z(i,j) = (B0z1(i+ci1,j) - B0z2(i+ci1,j))&
              -(B0z1(i-ci1,j) - B0z2(i-ci1,j))
  end forall

! Magnetic field generated by 4 magnets

  else if(magnet_no.eq.4) then
  write(magnets,'(I0.0)') magnet_no

    li1 = 0.5
    ci1 = nint((nx*li1)/lx)

    lj1 = 0.5
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

!  IF ((J-CJ.lt.NY1+1).or.(J-CJ.gt.NY2-1)) THEN
!   B0z(I,J)=0.
!  ELSE IF ((J+CJ.lt.NY1+1).or.(J+CJ.gt.NY2-1)) THEN
!   B0z(I,J)=0.
!  ELSE IF ((I+CI.lt.NX1+1).or.(I+CI.gt.NX2-1)) THEN
!   B0z(I,J)=0.
!  ELSE IF ((I+CI.lt.NX1+1).or.(I+CI.gt.NX2-1)) THEN

  file_name='B0z_field_'
  ext = '.dat'
  outfile=trim(file_name)//trim(magnets)//ext
  open(1,file=outfile)
  2 format(1x,E12.5,2x,E12.5,2x,E12.5)
  do i=nx1,nx2
    write(1,2)
  do j=ny1,ny2
    write(1,2) x(i)-lx/2.0, y(j)-ly/2.0, B0z(i,j)
  end do
  end do
  close(1)

  !call gnuplot(a,b)

end program mccaig_magnetic_field

! Subroutines

! GNUPLOT file subroutine

subroutine gnuplot(a, b)

  open(10, file="multiplot.gp")
  write(10,*) "reset; set size square; unset key; set ticscale 0"
  write(10,*) "set multiplot layout 1,2 columnsfirst"
  write(10,*) "a=",b,"; b=",b
  write(10,*) "set object rect at 0,0 size ",a,",",b,"front fs empty"
  write(10,*) "p 'B0z_field.dat' w image"
  write(10,*) "set hidden3d; set xyplane at -0.1"
  write(10,*) "sp 'B0z_field.dat' w l lt -1"
  write(10,*) "unset multiplot"
  close(10)

return
end subroutine
