! Boundary conditions subroutine
! This subroutine contains four boundary conditions for each boundary:
! No slip, free slip, outflow, and inflow.

subroutine boundary_conditions(u1,u2,v1,v2,nx,ny)
implicit none

	integer, value :: nx,ny
	real, dimension(0:nx+1,0:ny+1), intent(inout) :: u1,u2,v1,v2
	integer :: i,j

! x=0,y - left
	do j=0,ny+1
		u1(0,j)=1.0
		u2(0,j)=1.0
		v1(0,j)=-v1(1,j)
		v2(0,j)=-v2(1,j)
	end do

  ! x=nx,y - right
	do j=0,ny+1
		u1(nx,j)=u1(nx-1,j)
		u2(nx,j)=u2(nx-1,j)
		v1(nx+1,j)=v1(nx,j)
		v2(nx+1,j)=v2(nx,j)
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
