Program Laplace
  implicit none
  real k,Pe,dx,dy,Se,Sw,Sn,Ss,q,dv,ue,uw,un,us
  real x0,xl,y0,yl,tolerance,residual, dt, time
  real pi
  integer nx,i,j,ny,max_iter,it,itmax
  real, allocatable::T(:,:),Ta(:,:),aP(:,:),aE(:,:),aW(:,:),aS(:,:),aN(:,:),SP(:,:),xc(:),x(:),yc(:),y(:),u(:,:),v(:,:)
  character*50 itchar, name
  !definiendo el tamaño de la malla
  !Ocupar mallas siempre de potencia de 2
  nx=40
  ny=40
  x0=0.0
  xl=1.0
  y0=0.0
  yl=1.0
  dx=(xl-x0)/float(nx)
  dy=(yl-y0)/float(ny)
  dv=dx*dy

  !definiendo las constantes
  q=0.
  pi =acos(-1.0)
  time=10.0
  dt=0.005
  itmax=int(time/dt)+1
  Pe=50.0
  k=1.0/Pe !es la función gamma, que se lee en  gamma*S/delta
  max_iter=1000
  tolerance= 1e-4

  !definiendo el tamaño del vector
  allocate(T(0:nx+1,0:ny+1),Ta(0:nx+1,0:ny+1),aP(nx,ny),aE(nx,ny),aW(nx,ny),aS(nx,ny),aN(nx,ny))
  allocate(SP(nx,ny),u(0:nx+1,0:ny+1),v(0:nx+1,0:ny+1),xc(0:nx+1),x(0:nx),yc(0:ny+1),y(0:ny))

  !llamando la rubrutina que genera la malla
  call Mesh1D(xc,x,x0,xl,nx)
  call Mesh1D(yc,y,y0,yl,ny)

  !precolocando los vectores
  aP=0.;aE=0.;aW=0.; aS=0.; aN=0.; SP=0.;  T=0.;u=0.;v=0.

  !calculando los vectores velocidades
  do i=0,nx+1
     do j=0,ny+1
        u(i,j)=-sin(pi*xc(i))*cos(pi*yc(j))
        v(i,j)=cos(pi*xc(i))*sin(pi*yc(j))
     enddo
  enddo
    
  !determinado las condiciones de frontera, el manejo de los arreglos es parecido a Matlab 
  T(:,0)=0.
  T(:,ny+1)=1.0
  

  !definiendo las áreas
  Se=dy; Sw=dy; Sn=dx; Ss=dx
  
  !creación del archivo de animación
  open(2, file='anim.gnp',status='replace')
  write(2,*)'set xrange[0:1.0];set yrange[0:1.0]; set view map; set size square'
  
  !determinando los coeficientes y agregando el paso temporal
  do it=1,itmax
     do i=1,nx
        do j=1,ny
        ue=0.5*(u(i,j)+u(i+1,j))
        uw=0.5*(u(i,j)+u(i-1,j))
        
        un=0.5*(v(i,j)+v(i,j+1))
        us=0.5*(v(i,j)+v(i,j-1))
        
        aE(i,j)=k*Se/dx -0.5*(ue*Se)
        aW(i,j)=k*Sw/dx +0.5*(uw*Sw)
        aN(i,j)=k*Sn/dy -0.5*(un*Sn)
        aS(i,j)=k*Ss/dy +0.5*(us*Ss)      
        
        SP(i,j)=(T(i,j))*(dv/dt) !En este caso dv=dx*dy
        aP(i,j)=aE(i,j) + aW(i,j) +  aN(i,j) +  aS(i,j) + dv/dt
     enddo    
  enddo
  !corrección de las condiciones de frontera
  !Cara Este
  !SP(nx,:)=0. !condición de frontera Newmann Sp=Sp=aw*phi'*dx  pero la derivada de phi (phi'=0) es igual con cero
  aP(nx,:)=aP(nx,:)-aE(nx,:)
  aE(nx,:)=0.0

  !Cara Oeste
  !SP(1,:)=0.
  aP(1,:)=aP(1,:)-aW(1,:)
  aW(1,:)=0.0

  !Cara Norte
  SP(:,ny)=SP(:,ny)+ 2.0*aN(:,ny)*T(1:nx,ny+1)
  aP(:,ny)=aP(:,ny)+aN(:,ny)
  aN(:,ny)=0.0
  
  !Cara Sur
  SP(:,1)=SP(:,1)+ 2.0*aS(:,1)*T(1:nx,0)
  aP(:,1)=aP(:,1)+aS(:,1)
  aS(:,1)=0.0
  !Para optimizar los códigos, una propuesta de saul fue asistir a un curso de ciencias computacionales

  !Gauss TDMA2D para un arreglo [A]{T}={S} 
  call Gauss_TDMA2D(T,nx,ny,aP,aE,aW,aN,aS,sP,nx,ny,max_iter,tolerance,residual)
 
 ! Condición de frontera Newmann dPhi/dx=phiE-phiP = 0 => PhiE=PhiP
  T(0,:)=T(1,:)
  T(nx+1,:)=T(nx,:)
 
 if (mod(it,100) == 0)then 
     call  WriteScalarField2D('Temp2D',it,T,xc,yc,nx,ny)
     write(itchar, '(i6)')it !Convertir a un entero a un caracter
     itchar=adjustl(itchar)
     name='Temp2D'//itchar(1:len_trim(itchar))//'.txt'
     write(2,*)'sp "'//name(1:len_trim(name))//'" w pm3d'
     write(2,*)"pause 0.5"
  endif

enddo
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
