program ubound
integer i,ei
real x(10)

do i=1,10
x(i)=i
enddo


!do i=1,10
ei=UBOUND(x,1)-1
write(*,*) x, ei
!enddo


end program 
