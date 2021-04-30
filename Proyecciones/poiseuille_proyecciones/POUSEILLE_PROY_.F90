program Lid_driven_proyeccion

implicit none

integer nx,ny,i,j,k,max_iter,itmax,it,ei,ej,c,max_iter2,kk
real x0,xl,y0,yl,dx,dy,dz,dlv,dt,time
real ue,uw,us,un,dbkx,dbky,Div
real Se,Sw,Sn,Ss,Re,residual,tolerance,gama
character*100 itchar

real, allocatable:: u(:,:),v(:,:),p(:,:),pp(:,:),aP(:,:),aE(:,:),aW(:,:),aN(:,:),aS(:,:),sP(:,:),x(:),y(:),yc(:),xc(:)
real, allocatable:: du(:,:),dv(:,:),uc(:,:),vc(:,:),u1(:,:),v1(:,:),um(:,:),vm(:,:)

nx=128
ny=64

max_iter=5000

allocate(u(0:nx,0:ny+1),v(0:nx+1,0:ny),p(0:nx+1,0:ny+1),pp(0:nx+1,0:ny+1),aP(nx,ny),aE(nx,ny),aW(nx,ny),aN(nx,ny),aS(nx,ny),	&
		sP(nx,ny),x(0:nx),y(0:ny),yc(0:ny+1),xc(0:nx+1),uc(0:nx+1,0:ny+1),vc(0:nx+1,0:ny+1),									&
		u1(0:nx,0:ny+1),v1(0:nx+1,0:ny),um(0:nx,0:ny+1),vm(0:nx+1,0:ny))

!Input
xl=10.0
x0=0.0
yl=1.0
y0=0.0

tolerance=1.0E-5
Re=100.0
gama=1.0/Re

!Construcción de la Malla

dx=(xl-x0)/float(nx)
dy=(yl-y0)/float(ny)
dz=1.0
dlv=dx*dy*dz

!Condicion de estabilidad


!dt=0.1*(0.5*Re)*(1./(1./dx**2.+1/dy**2.))S
time=10.0
dt=0.005
itmax=int(time/dt)+1   

Sw=dy*dz
Se=dy*dz
Ss=dx*dz
Sn=dx*dz


call Mesh1D(xc,x,x0,xl,nx)
call Mesh1D(yc,y,y0,yl,ny)

!Conciciones iniciales y de frontera
u=1.0;	v=0.0;	p=0.0;	!um=0.0;	vm=0.0

u(:,0)=0.0
u(:,ny+1)=0.0

!Checar
um=u
vm=v

!archivo de animación
OPEN(1,FILE='anim.gnp',STATUS='REPLACE')
WRITE(1,*)'set size square; set xrange[0:10]; set yrange[0:1]; unset key'

aP=0.0;	aE=0.0;	aW=0.0;	aN=0.0;	aS=0.0;	sP=0.0

do it=1, itmax

!************************************************************************************************************************
!Ecuación de u

aP=0.0;	aW=0.0;	aE=0.0;	aS=0.0;	aN=0.0;	sP=0.0

ei=ubound(u,1)-1
ej=ubound(u,2)-1

do i=1, ei
	do j=1, ej

		ue=0.5*(u(i,j)+u(i+1,j))
		uw=0.5*(u(i-1,j)+u(i,j))
		un=0.5*(v(i,j)+v(i+1,j))
		us=0.5*(v(i,j-1)+v(i+1,j-1))


		aE(i,j)=gama*Se/dx	-ue*Se/2.
		aW(i,j)=gama*Sw/dx	+uw*Sw/2.
		aN(i,j)=gama*Sn/dy	-un*Sn/2.
		aS(i,j)=gama*Ss/dy	+us*Ss/2.

		aP(i,j)=aE(i,j)+aW(i,j)+aN(i,j)+aS(i,j)
		sP(i,j)=0.0

	enddo
enddo

!Condiciones de frontera

dbkx=0.0
dbky=1.0

!Este
aP(ei,1:ej)=aP(ei,1:ej)-aE(ei,1:ej)
!sP(ei,1:ej)=sP(ei,1:ej)+(1.0+dbkx)*aE(ei,1:ej)*u(ei+1,1:ej)
aE(ei,1:ej)=0.0

!Oeste
aP(1,1:ej)=aP(1,1:ej)+dbkx*aW(1,1:ej)
sP(1,1:ej)=sP(1,1:ej)+(1.0+dbkx)*aW(1,1:ej)*u(0,1:ej)
aW(1,1:ej)=0.0

!Norte
aP(1:ei,ej)=aP(1:ei,ej)+dbky*aN(1:ei,ej)
sP(1:ei,ej)=sP(1:ei,ej)+(1.0+dbky)*aN(1:ei,ej)*u(1:ei,ej+1)
aN(1:ei,ej)=0.0

!Sur
aP(1:ei,1)=aP(1:ei,1)+dbky*aS(1:ei,1)
sP(1:ei,1)=sP(1:ei,1)+(1.0+dbky)*aS(1:ei,1)*u(1:ei,0)
aS(1:ei,1)=0.0

do i=1,ei
	do j=1,ej

	um(i,j)=u(i,j)+(dt/dlv)*(-aP(i,j)*u(i,j)+aE(i,j)*u(i+1,j)+aW(i,j)*u(i-1,j)+aN(i,j)*u(i,j+1)+aS(i,j)*u(i,j-1)+sP(i,j))

	enddo
enddo
!u=um
write(*,*)it,maxval(um)

!*****************************************************************************************************************************
!Ecuación para v

aP=0.0;	aW=0.0;	aE=0.0;	aS=0.0;	aN=0.0;	sP=0.0
ei=ubound(v,1)-1
ej=ubound(v,2)-1

!Cálculo de coeficientes

do i=1, ei
	do j=1, ej

		ue=0.5*(u(i,j)+u(i,j+1))
		uw=0.5*(u(i-1,j)+u(i-1,j+1))
		un=0.5*(v(i,j)+v(i,j+1))
		us=0.5*(v(i,j)+v(i,j-1))


		aE(i,j)=gama*Se/dx-ue*Se/2.
		aW(i,j)=gama*Sw/dx+uw*Sw/2.
		aN(i,j)=gama*Sn/dy-un*Sn/2.
		aS(i,j)=gama*Ss/dy+us*Ss/2.

		aP(i,j)=aE(i,j)+aW(i,j)+aN(i,j)+aS(i,j)
		sP(i,j)=0.0

	enddo
enddo

!Condiciones de frontera

dbkx=1.0
dbky=0.0

!Este
aP(ei,1:ej)=aP(ei,1:ej)-aE(ei,1:ej)
!sP(ei,1:ej)=sP(ei,1:ej)+(1.0+dbkx)*aE(ei,1:ej)*v(ei+1,1:ej)
aE(ei,1:ej)=0.0

!Oeste
aP(1,1:ej)=aP(1,1:ej)+dbkx*aW(1,1:ej)
sP(1,1:ej)=sP(1,1:ej)+(1.0+dbkx)*aW(1,1:ej)*v(0,1:ej)
aW(1,1:ej)=0.0

!Norte
aP(1:ei,ej)=aP(1:ei,ej)+dbky*aN(1:ei,ej)
sP(1:ei,ej)=sP(1:ei,ej)+(1.0+dbky)*aN(1:ei,ej)*v(1:ei,ej+1)
aN(1:ei,ej)=0.0

!Sur
aP(1:ei,1)=aP(1:ei,1)+dbky*aS(1:ei,1)
sP(1:ei,1)=sP(1:ei,1)+(1.0+dbky)*aS(1:ei,1)*v(1:ei,0)
aS(1:ei,1)=0.0


do i=1,ei
	do j=1,ej

	vm(i,j)=v(i,j)+(dt/dlv)*(-aP(i,j)*v(i,j)+aE(i,j)*v(i+1,j)+aW(i,j)*v(i-1,j)+aN(i,j)*v(i,j+1)+aS(i,j)*v(i,j-1)+sP(i,j))

	enddo
enddo

!********************************************************************************************************************************
!Ecuación para p


aP=0.0;	aW=0.0;	aE=0.0;	aS=0.0;	aN=0.0;	sP=0.0
ei=ubound(p,1)-1
ej=ubound(p,2)-1

um(nx,:)=um(nx-1,:)
vm(nx+1,:)=vm(nx,:)

!Cálculo de coeficientes

do i=1, ei
	do j=1, ej

		ue=um(i,j)
		uw=um(i-1,j)
		un=vm(i,j)
		us=vm(i,j-1)
		

		aE(i,j)=Se/dx
		aW(i,j)=Sw/dx
		aN(i,j)=Sn/dy
		aS(i,j)=Ss/dy

		aP(i,j)=aE(i,j)+aW(i,j)+aN(i,j)+aS(i,j)
		sP(i,j)=-(ue*Se-uw*Sw+un*Sn-us*Ss)/dt

	enddo
enddo

!Condiciones de frontera
p(ei+1,1:ej)=p(ei,1:ej)
p(0,1:ej)=p(1,1:ej)
p(1:ei,ej+1)=p(1:ei,ej)
p(1:ei,0)=p(1:ei,1)

dbkx=1.0
dbky=1.0


!Este
aP(ei,1:ej)=aP(ei,1:ej)+dbkx*aE(ei,1:ej)
sP(ei,1:ej)=sP(ei,1:ej)+(1.0+dbkx)*aE(ei,1:ej)*p(ei+1,1:ej)
aE(ei,1:ej)=0.0

!Oeste
aP(1,1:ej)=aP(1,1:ej)+dbkx*aW(1,1:ej)
sP(1,1:ej)=sP(1,1:ej)+(1.0+dbkx)*aW(1,1:ej)*p(0,1:ej)
aW(1,1:ej)=0.0

!Norte
aP(1:ei,ej)=aP(1:ei,ej)+dbky*aN(1:ei,ej)
sP(1:ei,ej)=sP(1:ei,ej)+(1.0+dbky)*aN(1:ei,ej)*p(1:ei,ej+1)
aN(1:ei,ej)=0.0

!Sur
aP(1:ei,1)=aP(1:ei,1)+dbky*aS(1:ei,1)
sP(1:ei,1)=sP(1:ei,1)+(1.0+dbky)*aS(1:ei,1)*p(1:ei,0)
aS(1:ei,1)=0.0




call Gauss_TDMA2D(p,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny,max_iter,tolerance,residual)

! p(nx+1,:)=p(nx,:)
! p(0,:)=p(1,:)
! p(:,ny+1)=p(:,ny)
! p(:,0)=p(:,1)

!*****************************************************************************************************************************
!Corrección de las velocidades

do i=1, nx-1
	do j=1, ny
	u(i,j)=um(i,j)-dt*(p(i+1,j)-p(i,j))/dx
	enddo
enddo
u(nx,:)=u(nx-1,:)

do i=1, nx
	do j=1, ny-1
	v(i,j)=vm(i,j)-dt*(p(i,j+1)-p(i,j))/dy
	enddo
enddo
v(nx+1,:)=v(nx,:)
write(*,*)it,residual
!Escritura de campo de velocidades 
	IF (MOD(it,100) .EQ. 0) THEN
  
		CALL interpolateToNodesUs(uc,u,nx,ny)
		CALL interpolateToNodesVs(vc,v,nx,ny)
  
		CALL WriteVectorField('vel',it,uc,vc,xc,yc,nx,ny)
  
		!para el archivo de animación
		WRITE(itchar,'(i6)')it
		itchar=ADJUSTL(itchar)
		itchar=itchar(1:LEN_TRIM(itchar))
		WRITE(1,*)"p 'vel"//itchar(1:LEN_TRIM(itchar))//".txt' u 1:2:(0.25*$3):(0.25*$4) w vec"
		WRITE(1,*)'pause 0.05'
	end if
enddo	!Fin del ciclo de tiempo

! call interpolateToNodesUs(uc,u,nx,ny)
! call interpolateToNodesVs(vc,v,nx,ny)
! 
! call WriteVectorField('uv',10,uc,vc,xc,yc,nx,ny)

end program


!*********************************************************************************************************************
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

End Subroutine

!*************************************************************************************************************
subroutine interpolateToNodesUs(uc,us,nx,ny)

implicit none
integer nx,ny
real:: us(0:nx,0:ny+1),uc(0:nx+1,0:ny+1)
integer bi,ei,bj,ej,i,j
bi=1; bj=1;ej=ny; ei=nx

! Internal points
do i=bi,ei
do j=bj-1,ej+1
	uc(I,J) = ( us(I-1,J) + us(I,J) ) * 0.5
end do 
end do

uc(bi-1,:)=us(0,:)
uc(ei+1,:)=us(nx,:)

end subroutine


subroutine interpolateToNodesVs(vc,vs,nx,ny)
implicit none
integer nx,ny
real:: vs(0:nx+1,0:ny),vc(0:nx+1,0:ny+1)
integer bi,ei,bj,ej,i,j
bi=1; bj=1;ej=ny; ei=nx

do i=bi-1,ei+1
do j=bj,ej
	vc(I,J) = ( vs(I,J) + vs(I,J-1) ) * 0.5
end do 
end do

vc(:,bj-1)=vs(:,0)
vc(:,ej+1)=vs(:,ny)
end subroutine

!************************************************************************************************************************

Subroutine WriteVectorField(Name,kx,uc,vc,xc,yc,nx,ny)
integer i,j,nx,ny,kx
real*4 uc(0:nx+1,0:ny+1),vc(0:nx+1,0:ny+1),xc(0:nx+1),yc(0:ny+1)
character*(*)Name
character*50 txt,Filename
write(txt,'(i6)')kx
txt=ADJUSTL(txt)
Filename=name//txt(1:len_trim(txt))//".txt"

open(11,file=Filename(1:len_trim(Filename)))
	do i=0,nx+1
	do j=0,ny+1
	write(11,*)xc(i),yc(j),uc(i,j),vc(i,j)
	end do
	write(11,*)''
	end do
close(11)
End Subroutine

!****************************************************************************************************************
!Solver

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

!***********************************************************************************************************************
