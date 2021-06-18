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
