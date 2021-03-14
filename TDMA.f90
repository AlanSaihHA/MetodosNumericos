Subroutine TDMA(nx,Tp,Tf,aE,aW,aE,aP,sP,T)
INTEGER nx,i
REAL Tp,Tf
REAL, allocatable:: A(0:nx+1),C(0:nx+1),T(0:nx+1),aW(nx),aE(nx),aP(nx),sP(nx)
allocate(A(0:nx+1),C(0:nx+1),T(0:nx+1)
A(0)=0; C(0)=Tp
A(nx+1)=0; C(nx+1)=0
do i=1, nx
	A(i)=aE(i)/(aP(i)-aW(i)*A(i-1))
	C(i)=(aW(i)*C(i-1)+sP(i))/(aP(i)-aW(i)*A(i-1))
enddo

do i=nx,1,-1
	write(*,*) i
	T(i)=A(i)*T(i+1)+C(i)
enddo

End Subroutine
