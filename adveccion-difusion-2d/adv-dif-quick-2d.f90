Program Laplace 
  real k,dx,Se,Sw,Sn,Ss,q,dv,ue,uw,un,us
  real x0,xl,y0,yl,tolerance,residual
  integer nx,i,j,ny,max_iter,m
  real pi
  real, allocatable::T(:,:),aP(:,:),aE(:,:),aW(:,:),aS(:,:),aN(:,:),SP(:,:),xc(:),x(:),yc(:),y(:),uc(:,:),vc(:,:)
  real, allocatable:: aWW(:,:),aSS(:,:),aEE(:,:),aNN(:,:)
  real, allocatable :: H(:,:)
!Inicialización de variables
x0 = 0.0; xl = 1.0
y0 = 0.0; yl = 1.0

nx = 50; ny = 50

Pi = acos(-1.0)
Pe = 200.0
gamma = 1.0/Pe

dx = (xl-x0)/float(nx)
dy = (yl-y0)/float(ny)

Se = dy; Sw = dy
Sn = dx; Ss = dx
dv = dx*dy
  !definiendo las constantes
  q=0.
  k=1.0 !es la función gamma, que se lee en  gamma*S/delta
  max_iter=1000
  tolerance= 1e-4
  
  pi =acos(-1.0)

  !definiendo el tamaño del vector
  allocate(T(0:nx+1,0:ny+1),aP(nx,ny),aE(nx,ny),aW(nx,ny),aS(nx,ny),aN(nx,ny),SP(nx,ny))
  allocate(aWW(nx,ny), aSS(nx,ny),aEE(nx,ny),aNN(nx,ny))
  allocate(x(0:nx),xc(0:nx+1),y(0:ny),yc(0:ny+1),uc(0:nx+1,0:ny+1),vc(0:nx+1,0:ny+1))
  allocate(H(0:nx+1,0:ny+1))
  !llamando la rubrutina que genera la malla
  call Mesh1D(xc,x,x0,xl,nx)
  call Mesh1D(yc,y,y0,yl,ny)
  
  !Escritura de resultados
open(1,file="velc-quick.txt", status="replace")

 !calculando los vectores velocidades
  do i=0,nx+1
     do j=0,ny+1
        uc(i,j)=-sin(pi*xc(i))*cos(pi*yc(j))
        vc(i,j)=cos(pi*xc(i))*sin(pi*yc(j))
     enddo
  enddo

do i = 0,nx+1  
do j=0,ny+1
	write(1,*)xc(i),yc(j),uc(i,j),vc(i,j)
end do
end do
close(1)

  !precolocando los vectores
  aP=0.;aE=0.;aW=0.; aS=0.; aN=0.; aWW=0.; aSS=0.; aEE=0.; aNN=0.; SP=0.; T=0.
  
!Condiciones de frontera
	!norte
	T(:,ny+1) = 1.0
	!sur
	T(:,0) = 0.0

  !definiendo las áreas
  Se=dy; Sw=dy; Sn=dx; Ss=dx
  
  !determinando los coeficientes 
  do i=1,ny
     do j=1,nx
		uw = 0.5*(uc(i,j) + uc(i-1,j))
		ue = 0.5*(uc(i,j) + uc(i+1,j))
		un = 0.5*(vc(i,j) + vc(i,j+1))
		us = 0.5*(vc(i,j) + vc(i,j-1))
		
		if (uw .gt. 0.0) then
			 aW(i,j) = gamma*Sw/dx + (1.0/8.0)*ue + (6.0/8.0)*uw
			 aWW(i,j) = -(1.0/8.0)*uw
		else 
			 aW(i,j) = gamma*Sw/dx + (3.0/8.0)*uw
		end if
		
		if (ue .gt. 0.0) then
			 aE(i,j) = gamma*Se/dx - (3.0/8.0)*ue
		else 
			 aE(i,j) = gamma*Se/dx - (6.0/8.0)*ue -(1.0/8.0)*uw
			 aEE(i,j) = (1.0/8.0)*ue
		end if
		
		if (us .gt. 0.0) then
			 aS(i,j) = gamma*Ss/dx + (1.0/8.0)*un + (6.0/8.0)*us
			 aSS(i,j) = - (1.0/8.0)*us
		else 
			 aS(i,j) = gamma*Ss/dx + (3.0/8.0)*us
		end if
		
		if (un .gt. 0.0) then
			 aN(i,j) = gamma*Sn/dx - (3.0/8.0)*un
		else 
			 aN(i,j) = gamma*Sn/dx - (6.0/8.0)*un -(1.0/8.0)*us
			 aNN(i,j) = (1.0/8.0)*un
		end if

		aP(i,j) = aE(i,j) + aW(i,j) + aN(i,j) + aS(i,j) &
			+ aWW(i,j) + aSS(i,j) + aEE(i,j) + aNN(i,j)
		Sp(i,j) = 0.0  !Termino fuente es cero no especificado en notas
     enddo    
  enddo
  
  
!Correción de condiciones de frontera
	!Este
	aP(nx,1:ny) = aP(nx,1:ny) - aE(nx,1:ny)
	aE(nx,1:ny) = 0.0
	aEE(nx,1:ny)=0.0
	!Oeste
	aP(1,1:ny) = aP(1,1:ny) - aW(1,1:ny)
	aW(1,1:ny) = 0.0
	aWW(1,1:ny) = 0.0
	!Norte
	aP(1:nx,ny) = aP(1:nx,ny) + aN(1:nx,ny)
	sP(1:nx,ny) = sP(1:nx,ny) + 2.0*aN(1:nx,ny)*T(1:nx,ny+1)
	aN(1:nx,ny) = 0.0
	aNN(1,1:ny) = 0.0
	!Sur
	aP(1:nx,1) = aP(1:nx,1) + aS(1:nx,1)
	sP(1:nx,1) = sP(1:nx,1) + 2.0*aS(1:nx,1)*T(1:nx,0)
	aS(1:nx,1) = 0.0
	aSS(1,1:ny) = 0.0

  !Para optimizar los códigos, una propuesta de saul fue asistir a un curso de ciencias computacionales

!Esta parte corresponde al método por Jacobi
count_iter=0
residual=1.0
H=0.0
do while ((count_iter <= max_iter).and.(residual > tolerance))  
	do j = 1,ny
		do i = 1,nx
			H(i,j) = (aE(i,j)*T(i+1,j) + aW(i,j)*T(i-1,j) &
			+ aN(i,j)*T(i,j+1) + aS(i,j)*T(i,j-1) &
			+ aWW(i,j)*T(i-2,j) + aEE(i,j)*T(i+2,j) &
 		        + aSS(i,j)*T(i,j-2) + aNN(i,j)*T(i,j+2) + sP(i,j))/aP(i,j)
		enddo
	enddo
	T=H
       residual = calcResidual(T,nx,ny,aP,aE,aW,aN,aS,sP,nx,ny)
	count_iter=count_iter+1
	write(*,*) count_iter
enddo 	
 
 
  !Gauss TDMA2D para un arreglo [A]{T}={S}  
  !call Gauss_TDMA2D(T,nx,ny,aP,aE,aW,aN,aS,sP,nx,ny,max_iter,tolerance,residual)

  call  WriteScalarField2D('Temp2D-quick',0,T,xc,yc,nx,ny)
  

!archivo paraview variable escalar
open(3,file='data-quick-esc.dat')                !Aquì se nombre el archivo en donde se escribiran los datos
write(3,*) 'TITLE= "TESTPLOT" '
write(3,*) 'VARIABLES="x","y","T"'              !Se dan coordenadas x,y, ademas la componente fx y fy. Ademas la magnitu de esas componentes
write(3,*) 'ZONE T="1", I=',nx+1,', J=',ny+1

do j=0, ny              !Se evalua para cada coordenada (x,y). Primerotodas las equis con cada ye. Se pasa a la siguiente equis y se vuelve a repetir
do i=0, nx 
	write(3,*) xc(i), yc(j), T(i,j)          !Evaluacion de cada dato
enddo 
enddo
close(3)
  
!archivo paraview variable velocidad
open(4,file='data-quick-vec.dat')                !Aquì se nombre el archivo en donde se escribiran los datos
write(4,*) 'TITLE= "TESTPLOT" '
write(4,*) 'VARIABLES="x","y","fx","fy","magf" '              !Se dan coordenadas x,y, ademas la componente fx y fy. Ademas la magnitu de esas componentes
write(4,*) 'ZONE T="1", I=',nx+1,', J=',ny+1

do j=0, ny              !Se evalua para cada coordenada (x,y). Primerotodas las equis con cada ye. Se pasa a la siguiente equis y se vuelve a repetir
do i=0, nx 
	write(4,*) xc(i), yc(j), uc(i,j), vc(i,j), sqrt(uc(i,j)*uc(i,j)+vc(i,j)*vc(i,j))       !Evaluacion de cada dato
enddo 
enddo
close(4)
  
  
endProgram Laplace

!!****************************************************************************************************************
  Subroutine Mesh1D(xc,x,x0,xl,nx)
    integer i,j,nx
    real*4 x0,xl,dx
    real*4 x(0:nx),xc(0:nx+1)
    dx=(1.0)/dfloat(nx)
    do i=0,nx
       x(i)=dfloat(i)*dx
       x(i)=x0+(xl-x0)*x(i)
    end do
    xc(0)=x(0); xc(nx+1)=x(nx)
    do i=1,nx
       xc(i)=(x(i)+x(i-1))*0.5
    end do
  End Subroutine Mesh1D
!!*******************************************************************************************************************

subroutine TDMA(x,a,b,c,d,n)

implicit none
 integer n,k
 real a(n),b(n),c(n),d(n),x(n),m

 do k=2,N
  m=a(k)/b(k-1)
  b(k)=b(k)-m*c(k-1)
  d(k)=d(k)-m*d(k-1)
 end do

 x(n)=d(n)/b(n)

 do k=n-1,1,-1
  x(k)=(d(k)-c(k)*x(k+1))/b(k)
 end do

end subroutine

subroutine lineX_2D(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny)

implicit none
integer bi,ei,bj,ej,i,j,nx,ny
real phi(0:ei+1,0:ej+1),aP(nx,ny),aE(nx,ny),aW(nx,ny),aN(nx,ny),aS(nx,ny),sP(nx,ny)
real a(ei),b(ei),c(ei),d(ei)
bi=1; bj=1
do j=bj,ej
do i=bi,ei
	a(i)=-aW(i,j)
	b(i)=aP(i,j)
	c(i)=-aE(i,j)
	d(i)=sp(i,j) + aN(i,j) * phi(i,j+1) + aS(i,j) * phi(i,j-1)
end do
 call TDMA(phi(bi:ei,j), a, b, c ,d ,ei)
end do

end subroutine

subroutine lineY_2D(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny)

implicit none
integer bi,ei,bj,ej,i,j,nx,ny
real phi(0:ei+1,0:ej+1),aP(nx,ny),aE(nx,ny),aW(nx,ny),aN(nx,ny),aS(nx,ny),sP(nx,ny)
real a(ej),b(ej),c(ej),d(ej)
bi=1; bj=1

do i=bi,ei
do j=bj,ej
	a(j)=-aS(i,j)
	b(j)=aP(i,j)
	c(j)=-aN(i,j)
	d(j)=sp(i,j) + aE(i,j) * phi(i+1,j) + aW(i,j) * phi(i-1,j) 
end do
 call TDMA(phi(i,bj:ej), a, b, c ,d ,ej)
end do

end subroutine

!****************************
!****************************

real*4 function calcResidual(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny)
implicit none
integer bi,ei,bj,ej,i,j,nx,ny
real*4 phi(0:ei+1,0:ej+1),aP(nx,ny),aE(nx,ny),aW(nx,ny),aN(nx,ny),aS(nx,ny),sP(nx,ny)
real*4 acum(ei,ej),residual,NINV
bi=1; bj=1
acum=0
NINV = 1.0 / dfloat(ei*ej)
do i=bi,ei
do j=bj,ej
acum(i,j) = aE(i,j) * phi(i+1,j) + aW(i,j) * phi(i-1,j) + aN(i,j) * phi(i,j+1) + &
aS(i,j) * phi(i,j-1) + sp(i,j) - aP(i,j) * phi(i,j)
end do
end do
residual = sqrt( NINV * sum(acum * acum) )
calcResidual=residual
end function


subroutine Gauss_TDMA2D(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny,max_iter,tolerance,residual)

implicit none

integer bi,ei,bj,ej,i,j,nx,ny,count_iter,max_iter
real*4 phi(0:ei+1,0:ej+1),aP(nx,ny),aE(nx,ny),aW(nx,ny),aN(nx,ny),aS(nx,ny),sP(nx,ny)
real*4 residual,tolerance
	
	interface
		real*4 function calcResidual(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny)
		implicit none
		integer bi,ei,bj,ej,i,j,nx,ny
		real*4 phi(0:ei+1,0:ej+1),aP(nx,ny),aE(nx,ny),aW(nx,ny),aN(nx,ny),aS(nx,ny),sP(nx,ny)
		end function
	end interface

count_iter=0;  residual=1.0

do while((count_iter <= max_iter).and.(residual > tolerance)) 
call lineX_2D(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny)
call lineY_2D(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny)
residual = calcResidual(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny)
count_iter=count_iter+1
end do

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!********************************************************************

Subroutine WriteScalarField2D(Name,kx,T,xc,yc,nx,ny)
integer i,j,nx,ny,kx
real*4 T(0:nx+1,0:ny+1),xc(0:nx+1),yc(0:ny+1)
character*(*)Name
character*50 txt,Filename
write(txt,'(i6)')kx
txt=ADJUSTL(txt)
Filename=name//txt(1:len_trim(txt))//".txt"

open(10,file=Filename(1:len_trim(Filename)))
	do j=0,ny+1
	do i=0,nx+1
	write(10,*)xc(i),yc(j),T(i,j)
	end do
		write(10,*)''
	end do
close(10)
End Subroutine
