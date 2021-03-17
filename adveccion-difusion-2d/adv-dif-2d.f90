Program Laplace 
  real k,dx,Se,Sw,Sn,Ss,q,dv,ue,uw,un,us
  real x0,xl,y0,yl,tolerance,residual
  integer nx,i,j,ny,max_iter
  real pi
  real, allocatable::T(:,:),aP(:,:),aE(:,:),aW(:,:),aS(:,:),aN(:,:),SP(:,:),xc(:),x(:),yc(:),y(:),uc(:,:),vc(:,:)

!Inicialización de variables
x0 = 0.0; xl = 1.0
y0 = 0.0; yl = 1.0

nx = 50; ny = 50

Pi = acos(-1.0)
Pe = 100.0
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
  allocate(x(0:nx),xc(0:nx+1),y(0:ny),yc(0:ny+1),uc(0:nx+1,0:ny+1),vc(0:nx+1,0:ny+1))

  !llamando la rubrutina que genera la malla
  call Mesh1D(xc,x,x0,xl,nx)
  call Mesh1D(yc,y,y0,yl,ny)
  
  !Escritura de resultados
open(1,file="velc.txt", status="replace")

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
  aP=0.;aE=0.;aW=0.; aS=0.; aN=0.; SP=0.; T=0.
  
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
		
		aE(i,j) = gamma*Se/dx - 0.5*ue*Se
		aW(i,j) = gamma*Sw/dx + 0.5*uw*Sw
		aN(i,j) = gamma*Sn/dy - 0.5*un*Sn
		aS(i,j) = gamma*Ss/dy + 0.5*us*Ss
		
		aP(i,j) = aE(i,j) + aW(i,j) + aN(i,j) + aS(i,j)
		Sp(i,j) = 0.0  !Termino fuente es cero no especificado en notas
     enddo    
  enddo
  

  !corrección de las condiciones de frontera
  
!Correción de condiciones de frontera
	!Este
	aP(nx,1:ny) = aP(nx,1:ny) - aE(nx,1:ny)
	aE(nx,1:ny) = 0.0
	!Oeste
	aP(1,1:ny) = aP(1,1:ny) - aW(1,1:ny)
	aW(1,1:ny) = 0.0
	
	!Norte
	aP(1:nx,ny) = aP(1:nx,ny) + aN(1:nx,ny)
	sP(1:nx,ny) = sP(1:nx,ny) + 2.0*aN(1:nx,ny)*T(1:nx,ny+1)
	aN(1:nx,ny) = 0.0
	!Sur
	aP(1:nx,1) = aP(1:nx,1) + aS(1:nx,1)
	sP(1:nx,1) = sP(1:nx,1) + 2.0*aS(1:nx,1)*T(1:nx,0)
	aS(1:nx,1) = 0.0

  !Para optimizar los códigos, una propuesta de saul fue asistir a un curso de ciencias computacionales

  !Gauss TDMA2D para un arreglo [A]{T}={S} 
  call Gauss_TDMA2D(T,nx,ny,aP,aE,aW,aN,aS,sP,nx,ny,max_iter,tolerance,residual)

  call  WriteScalarField2D('Temp2D',0,T,xc,yc,nx,ny)
  
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
