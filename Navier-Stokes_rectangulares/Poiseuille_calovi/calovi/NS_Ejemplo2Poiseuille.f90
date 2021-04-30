Program NavierStokesLIDCAVITY
!!!LID DRIVEN CAVITY
  implicit none
  CHARACTER*100 itchar
  real gama,Re,dx,dy,Se,Sw,Sn,Ss,q,dv,ue,uw,un,us,ot
  real x0,xl,y0,yl,tolerance,residual,dt, time,Div,dbkx,dbky
  real pi,a
  integer nx,i,j,c,ny,max_iter,it,itmax,ei,ej,bi,bj
  real, allocatable::Pp(:,:),P(:,:),de(:,:),dn(:,:),aP(:,:),aE(:,:),aW(:,:),aS(:,:),aN(:,:),SP(:,:)
  real, allocatable::xc(:),x(:),yc(:),y(:),u1(:,:),v1(:,:),u(:,:),v(:,:),uc(:,:),vc(:,:)
  !Pp,u0,v0 valores conocidos de las iteraciones anteriores 

  !definiendo el tamaño de la malla
  !Ocupar mallas siempre de potencia de 2
  nx=200
  ny=60
  x0=0.0
  xl=10.0
  y0=0.0
  yl=1.0
  dx=(xl-x0)/float(nx)
  dy=(yl-y0)/float(ny)
  dv=dx*dy

  !definiendo las constantes
  !  q=0.
  pi =acos(-1.0)
  time=10.0
  dt=0.005
  itmax=int(time/dt)+1   
  Re=100.0
  gama=1.0/Re !es la función gamma, que se lee en  gamma*S/delta
  max_iter=1000
  tolerance= 1e-5

  !definiendo el tamaño del vector
  allocate(Pp(0:nx+1,0:ny+1),P(0:nx+1,0:ny+1),aP(nx,ny),aE(nx,ny),aW(nx,ny),aS(nx,ny),aN(nx,ny))
  allocate(SP(nx,ny),u(0:nx,0:ny+1),v(0:nx+1,0:ny),u1(0:nx,0:ny+1),v1(0:nx+1,0:ny),xc(0:nx+1),x(0:nx),yc(0:ny+1),y(0:ny))
  allocate(de(0:nx+1,0:ny+1),dn(0:nx+1,0:ny+1),uc(0:nx+1,0:ny+1),vc(0:nx+1,0:ny+1))

  !llamando la rubrutina que genera la malla
  call Mesh1D(xc,x,x0,xl,nx)
  call Mesh1D(yc,y,y0,yl,ny)

  !precolocando los vectores
  u=1.0;v=0;P=0.0;Pp=0.0;
  dn=0.0;de=0.0
  u(:,ny+1)=0.0;u(:,0)=0.0;
  !definiendo las áreas
  Se=dy; Sw=dy; Sn=dx; Ss=dx
  u1=u;v1=v
  !archivo de animación
  OPEN(1,FILE='anim.gnp',STATUS='REPLACE')
  WRITE(1,*)'set size square; set xrange[0:10]; set yrange[0:1]; unset key'

  do it=1,itmax !ciclo temporal
     Div=1.0
     c=1
     do while((Div .GE. tolerance) .AND. (c .LT. max_iter)) !ciclo simplec
        aP=0.0;aE=0.0;aW=0.0; aS=0.0; aN=0.0; SP=0.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !
        !*******************************Empieza ecuación de u1*******************************************
        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ei=ubound(u,1)-1
        ej=ubound(u,2)-1
        bi=lbound(u,1)+1
        bj=lbound(u,2)+1
write(*,*)ei,ej,bi,bj
        do j=bj,ej
           do i=bi,ei
              ue=0.5*(u(i,j)+u(i+1,j))
              uw=0.5*(u(i,j)+u(i-1,j))
              un=0.5*(v(i,j)+v(i+1,j))
              us=0.5*(v(i,j-1)+v(i+1,j-1))
              !utilizando esquema central sólo por este momento
              aE(i,j)=gama*Se/dx-0.5*(ue*Se)
              aW(i,j)=gama*Sw/dx +0.5*(uw*Sw)              
              aN(i,j)=gama*Sn/dy -0.5*(un*Sn)              
              aS(i,j)=gama*Ss/dy +0.5*(us*Ss)              
              SP(i,j)=(u(i,j))*(dv/dt) -(P(i+1,j)-P(i,j))*dv/dx  !En este caso dv=dx*dy              
              aP(i,j)=aE(i,j) + aW(i,j) +  aN(i,j) +  aS(i,j) + dv/dt
           enddo
        enddo
        !correción a la condición de frontera
        dbkx=0.0
        dbky=1.0
        !Cara Este
        aP(ei,1:ej)=aP(ei,1:ej)-aE(ei,1:ej) !condición Neumann
       ! SP(ei,1:ej)=SP(ei,1:ej)+(1.0+dbkx)*aE(ei,1:ej)*u(ei+1,1:ej)
        aE(ei,1:ej)=0.0

        !Cara Oeste  
        aP(bi,1:ej)=aP(1,1:ej)+dbkx*aW(1,1:ej)
        SP(bi,1:ej)=SP(1,1:ej)+(1.0+dbkx)*aW(1,1:ej)*u(0,1:ej)
        aW(bi,1:ej)=0.0

        !Cara Norte
        SP(1:ei,ej)=SP(1:ei,ej)+ (1+dbky)*aN(1:ei,ej)*u(1:ei,ej+1)
        aP(1:ei,ej)=aP(1:ei,ej)+dbky*aN(1:ei,ej)
        aN(1:ei,ej)=0.0

        !Cara Sur
        SP(1:ei,1)=SP(1:ei,1)+ (1+dbky)*aS(1:ei,1)*u(bi:ei,0)
        aP(1:ei,1)=aP(1:ei,1)+dbky*aS(1:ei,1)
        aS(1:ei,1)=0.0
        
        do j=bj,ej           
           do i=bi,ei
              de(i,j)=Se/(ap(i,j)-(aE(i,j)+aW(i,j)+aN(i,j)+aS(i,j)))
           enddo
        enddo

        !Gauss TDMA2D para un arreglo [A]{T}={S}
        call Gauss_TDMA2D(u1,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny,max_iter,tolerance,residual)
!###########################################################################################################
  aP=0.0; aE=0.0; aW=0.0; aN=0.0; aS=0.0; Sp=0.0         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !
        !*******************************Empieza ecuación de v1*******************************************
        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ei=ubound(v,1)-1
        ej=ubound(v,2)-1
        bi=lbound(v,1)+1
        bj=lbound(v,2)+1
        do j=bj,ej
           do i=bi,ei
              ue=0.5*(u(i,j)+u(i,j+1))
              uw=0.5*(u(i-1,j)+u(i-1,j+1))
              un=0.5*(v(i,j)+v(i,j+1))
              us=0.5*(v(i,j)+v(i,j-1))
              !utilizando esquema central sólo por este momento
              aE(i,j)=gama*Se/dx-0.5*(ue*Se)
              aW(i,j)=gama*Sw/dx +0.5*(uw*Sw)              
              aN(i,j)=gama*Sn/dy -0.5*(un*Sn)              
              aS(i,j)=gama*Ss/dy +0.5*(us*Ss)              
              SP(i,j)=(v(i,j))*(dv/dt) -(P(i,j+1)-P(i,j))*dv/dy              !En este caso dv=dx*dy              
              aP(i,j)=aE(i,j) + aW(i,j) +  aN(i,j) +  aS(i,j) + dv/dt
           enddo
        enddo
        !correción a la condición de frontera
        dbkx=1.0
        dbky=0.0
        !Cara Este
        aP(ei,1:ej)=aP(ei,1:ej)-aE(ei,1:ej)
       ! SP(ei,1:ej)=SP(ei,1:ej)+(1.0+dbkx)*aE(ei,1:ej)*v(ei+1,1:ej)
        aE(ei,1:ej)=0.0

        !Cara Oeste  
        aP(bi,1:ej)=aP(1,1:ej)+dbkx*aW(1,1:ej)
        SP(bi,1:ej)=SP(1,1:ej)+(1.0+dbkx)*aW(1,1:ej)*v(bj-1,bj:ej)
        aW(bi,1:ej)=0.0

        !Cara Norte
        SP(1:ei,ej)=SP(1:ei,ej)+ (1+dbky)*aN(1:ei,ej)*v(bi:ei,ej+1)
        aP(1:ei,ej)=aP(1:ei,ej)+dbky*aN(1:ei,ej)
        aN(1:ei,ej)=0.0

        !Cara Sur
        SP(1:ei,bj)=SP(1:ei,1)+ (1+dbky)*aS(1:ei,1)*v(bi:ei,0)
        aP(1:ei,bj)=aP(1:ei,bj)+dbky*aS(1:ei,1)
        aS(1:ei,bj)=0.0
        !Gauss TDMA2D para un arreglo [A]{T}={S}
        call Gauss_TDMA2D(v1,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny,max_iter,tolerance,residual)  
        do j=bj,ej           
           do i=bi,ei
              dn(i,j)=Sn/(ap(i,j)-(aE(i,j)+aW(i,j)+aN(i,j)+aS(i,j)))
           enddo
        enddo
	!CORRECIÓN DE LA PRESIÓN
	aP=0.0; aE=0.0; aW=0.0; aN=0.0; aS=0.0; Sp=0.0; Pp=0.0 ! Pp=p' P=P* + P'
        ei=ubound(Pp,1)-1
        ej=ubound(Pp,2)-1
        bi=lbound(Pp,1)+1
        bj=lbound(Pp,2)+1

	!Cálculo de coeficientes
        do j=bj,ej           
           do i=bi,ei
              !Flujos en las caras
              ue=u1(i,j)
              uw=u1(i-1,j)
              un=v1(i,j)
              us=v1(i,j-1)
              aE(i,j)=de(i,j)*Se
              aW(i,j)=de(i-1,j)*Sw
              aN(i,j)=dn(i,j)*Sn
              aS(i,j)=dn(i,j-1)*Ss
              ap(i,j)=aE(i,j) + aW(i,j) +  aN(i,j) +  aS(i,j)
              SP(i,j)=-(ue*Se-uw*Sw+un*Sn-us*Ss)
           enddo
        enddo
        !dado que estamos resolviendo para la corrección de la presión no es necesario correjir en las fronteras
        call Gauss_TDMA2D(Pp,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny,max_iter,tolerance,residual)
        !Correción de la presión
        P=P+Pp !P = P* + P'
        !esta parte se entiende como u=u*+d()*(P()-P()) y v=v*+d()*(P()-P())
        do i=1,nx-1
           do j=1,ny
              u1(i,j)=u1(i,j) + de(i,j)*(Pp(i,j)-Pp(i+1,j))
           enddo
        enddo
        !corrección de velocidades en la dirección u para condición de frontera Newmann        
        u1(nx,:)=u1(nx-1,:)
        do i=1,nx
           do j=1,ny-1
              v1(i,j)=v1(i,j) + dn(i,j)*(Pp(i,j)-Pp(i,j+1))
           enddo
        enddo
        !corrección de velocidades en la dirección u para condición de frontera Newmann        
        v1(nx+1,:)=v1(nx,:)
        
        Div=maxval(abs(Sp(1:ei,1:ej)))
        c=c+1
        !finaliza el ciclo SIMPLEC
     end do 
     u=u1; v=v1
     write(*,*)it,c,Div
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
  END IF
  enddo !ciclo temporal
! SOLUCIÓN ANÁLITICA

  !Escritura de resultados
open(1,file="analitica.txt", status="replace")

! do i = 0,nx+1  
do j=0,ny+1
! write(*,*)'Esto es la presión dp/dx',(Pp(nx-6,j)-Pp(nx-7,j))/(2*dx)
	write(1,*) xc(nx+1),yc(j),yc(j)*((yl-y0)-yc(j)), 0.0
end do
! end do
close(1)
!  call  interpolateToNodesUs(uc,u,nx,ny)
!  call interpolateToNodesVs(vc,v,nx,ny)
                       !(name,un entero cualsea,uc,vc,xc,yc,nx,ny) 
!  call WriteVectorField('vel',0,uc,vc,xc,yc,nx,ny)
  end Program NavierStokesLIDCAVITY
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
End Subroutine WriteScalarField2D
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

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

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

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
