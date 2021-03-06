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
