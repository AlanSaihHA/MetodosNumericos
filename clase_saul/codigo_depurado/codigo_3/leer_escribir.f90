program leer_escribir
implicit none
real a(1000,2)
integer i,b, REASON, nps,v
real xcc, ycc
real, allocatable :: pt(:,:), pts(:,:)
xcc=4.0
ycc=2.0


open(20,file='NACA060006_int.txt')
i=1
DO 
   READ(20,*,IOSTAT=REASON)  a(i,1), a(i,2)
   IF (REASON>0) THEN
      print*, '...algo anda mal...'
        exit
   ELSE IF (REASON<0) THEN
      print*, '...fin del archivo...'
         exit
   ELSE
      print*, 'linea: ',i,'leída'
      i=i+1
   END IF
END DO
close(20)

nps=i-1

do i=1, nps
   write(*,*) i, a(i,1), a(i,2)
enddo

write(*,*) 'nps = ',nps
!print*, a(38,1), a(38,2)

allocate(pt(nps,2))
do i=1, nps      !Se colocan las coordenadas en el lugar del objeto
   pt(i,1) = xcc + a(i,1)
   pt(i,2) = ycc + a(i,2)
enddo

allocate(pts(nps/2,4))
v=0
do i=1, nps/2
   write(*,*) i, pt(i,1),pt(i,2),pt(nps/2+i,1) ,pt(nps-v,2)
   pts(i,1)= pt(i,1)
   pts(i,2)=pt(i,2)
   pts(i,3)=pt(nps/2+i,1)
   pts(i,4)=pt(nps-v,2)
   v=i
enddo

write(*,*) 'espacio'
do i=1, nps/2
write (*,*) pts(i,1),pts(i,2),pts(i,3),pts(i,4)
enddo

!pts(i,1),pts(i,2),pts(i,3),pts(i,4)
!derecha,arriba,izquieda,abajo

r=r0
do i=1, nxp2
	do j=1,nyp2
	        x=float(i-1)/hxi
       	y=float(j-1)/hyi      	
       	do m=1,nps/2
		if ((x .gt. pts(i,3)) .and. (x .lt. pts(i,1)) )  then
			if ((y .gt. pts(i,4)) .and. (x .lt. pts(i,2)) ) then
				r(i,j)=r(i,j)+prop1
			endif
		endif
		enddo
	enddo
enddo

end program leer_escribir

!------------------------------------------------------------------
      subroutine setdens(r0,nfronts,xcc,ycc,r,nxp2,nyp2,prop1) 

	  use grid; use bounds
	  use fparam
	  implicit none 
	  
	  integer nfronts,nxp2,nyp2,i,j,l
	  real*8  x,y,xx,yy,xfl,pi,r0,xcc,ycc,prop1
	  real*8:: r(nxp2,nyp2)

      pi=4.*ATAN(1.)

      do j=1,nyp2
      do i=1,nxp2
  	    r(i,j)=r0
      enddo
      enddo
      do j=1,nyp2
      do i=1,nxp2 
       x=float(i-1)/hxi
       y=float(j-1)/hyi
             xx=abs(x-xcc)
             yy=abs(y-ycc)

            xfl=sqrt(xx**2+yy**2)-0.5d0

      if(xfl.le. 0d0)r(i,j)=r(i,j)+prop1
      enddo
      enddo

