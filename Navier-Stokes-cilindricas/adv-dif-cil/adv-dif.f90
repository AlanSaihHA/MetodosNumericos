Program adv_dif_cil
implicit none
integer n0,nr,i,j,max_iter,it,itmax
real*4 dr,d0,r0,rl,th0,thl,pi,Se,Sw,Sn,Ss,dv,rn,rs,rp,residual,tolerance
real*4 ue,uw,vn,vs,Pe,gamma,time,dt
real*4, allocatable:: aE(:,:),aW(:,:),aN(:,:),aS(:,:),aP(:,:),sP(:,:)
real*4, allocatable:: r(:),rc(:),th(:),thc(:),T(:,:),u0(:,:),ur(:,:)
character*50 itchar,name

pi = acos(-1.0)
n0 = 100; nr = 50

max_iter = 1000
tolerance = 1e-5

r0 = 0.1; rl = 1.0
th0 = 0.0; thl = 2.0*pi
d0 = (thl-th0)/float(n0)
dr = (rl-r0)/float(nr)
Se = dr; Sw = dr

Pe = 10.0
gamma = 1/Pe

time = 3.0
dt = 0.0001
itmax = int(time/dt) + 1

allocate(aE(n0,nr),aW(n0,nr),aN(n0,nr),aS(n0,nr),aP(n0,nr),sP(n0,nr))
allocate(T(0:n0+1,0:nr+1),r(0:nr),rc(0:nr+1),th(0:n0),thc(0:n0+1))
allocate(u0(0:n0+1,0:nr+1),ur(0:n0+1,0:nr+1))

!Creación de malla
call Mesh_1D(n0,th0,thl,th,thc)
call Mesh_1D(nr,r0,rl,r,rc)

!Campo de velocidades
do i=0,n0+1
	do j=0,nr+1
		u0(i,j) = 0.1/rc(j) + rc(j)
		ur(i,j) = 0.0
	end do
end do

!Condición de frontera
T = 0.0				!Cilindro interno
T(:,nr+1) = 0.5*(cos(thc(:) + 0.5*pi) + 1.0)		!Cilindro externo
! cos(thc(:) + 0.5*pi) =- sin(thc(:))
!archivo de animación
OPEN(2,FILE='anim.gnp',STATUS='REPLACE')
WRITE(2,*)'set size square; set xrange[-1:1]; set yrange[-1:1]; set view map'

!#######################################################################
!Ciclo temporal
!#######################################################################
do it=1, itmax

	!Cálculo de coeficientes
	do i=1,n0
		do j=1,nr
			rp = rc(j)
			rn = r(j)
			rs = r(j-1)
			
			!áreas y volumen
			Sn = rn*d0
			Ss = rs*d0
			dv = rp*dr*d0
			
			ue = 0.5*(u0(i,j) + u0(i+1,j))
			uw = 0.5*(u0(i,j) + u0(i-1,j))
			vn = 0.5*(ur(i,j) + ur(i,j+1))
			vs = 0.5*(ur(i,j) + ur(i,j-1))
			
			aE(i,j) = gamma*Se/(rp*d0) - 0.5*ue*Se
			aW(i,j) = gamma*Sw/(rp*d0) + 0.5*uw*Sw
			aN(i,j) = gamma*Sn/dr - 0.5*vn*Sn
			aS(i,j) = gamma*Ss/dr + 0.5*vs*Ss
			aP(i,j) = aE(i,j)+aW(i,j)+aN(i,j)+aS(i,j) + dv/dt
			sP(i,j) = T(i,j)*dv/dt
		end do
	end do

	!Empalme por periodicidad
	T(n0+1,:) = T(1,:)
	T(0,:) = T(n0,:)
	
	!Corrección de condiciones de frontera
	!ESTE
	aP(n0,:) = aP(n0,:) + aE(n0,:)
	Sp(n0,1:nr) = Sp(n0,1:nr) + 2.0*aE(n0,1:nr)*T(n0+1,1:nr)
	aE(n0,:) = 0.0
	!OESTE
	aP(1,:) = aP(1,:) + aW(1,:)
	Sp(1,1:nr) = Sp(1,1:nr) + 2.0*aW(1,1:nr)*T(0,1:nr)
	aW(1,:) = 0.0
	!NORTE
	aP(:,nr) = aP(:,nr) + aN(:,nr)
	Sp(1:n0,nr) = Sp(1:n0,nr) + 2.0*aN(1:n0,nr)*T(1:n0,nr+1)
	aN(:,nr) = 0.0
	!SUR
	aP(1:n0,1) = aP(1:n0,1) + aS(1:n0,1)
	Sp(1:n0,1) = Sp(1:n0,1) + 2.0*aS(1:n0,1)*T(1:n0,0)
	aS(:,1) = 0.0

	call Gauss_TDMA2D(T,n0,nr,aP,aE,aW,aN,aS,Sp,n0,nr,max_iter,tolerance,residual)
	
	!Reasignamos valores de T, para deshacer desempalme
	T(0,:) = 0.5*(T(1,:) + T(n0,:))
	T(n0+1,:) = T(0,:)
	
	write(*,*)it
	
	!Escritura de resultados
	if(mod(it,100).eq.0) then
		call WriteScalarField2D('T',it,T,thc,rc,n0,nr)
		!Escritura en el archivo de gnuplot
		write(itchar,"(i6)")it
		itchar = adjustl(itchar)
		name = "T"//itchar(1:len_trim(itchar))//".txt"
		write(2,*)"sp '"//name(1:len_trim(name))//"' w pm3d"
		write(2,*)"pause 0.5"
	end if
end do
!#######################################################################
!Termina ciclo temporal
!#######################################################################
close(2)
end Program


!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
SUBROUTINE MESH_1D(nx,x0,xl,x,xc)
IMPLICIT NONE
INTEGER nx,i
REAL*4 dx,x0,xl
REAL*4 x(0:nx),xc(0:nx+1)

dx=1.0/FLOAT(nx)

DO i=0,nx
  x(i)=x0+FLOAT(i)*(xl-x0)*dx
END DO

xc(0)=x(0); xc(nx+1)=x(nx)

DO i=1,nx
  xc(i)=0.5*(x(i)+x(i-1))
END DO
END SUBROUTINE


!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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


!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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


!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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


!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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


!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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


!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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
			write(10,*)yc(j)*cos(xc(i)),yc(j)*sin(xc(i)),T(i,j)
		end do
		write(10,*)''
	end do
close(10)
End Subroutine
