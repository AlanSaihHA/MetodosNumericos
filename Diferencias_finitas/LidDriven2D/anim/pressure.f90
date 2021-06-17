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
