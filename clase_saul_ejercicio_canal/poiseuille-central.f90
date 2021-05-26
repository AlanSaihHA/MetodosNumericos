PROGRAM LID_DRIVEN
IMPLICIT NONE
!Declaración de variables
CHARACTER*100 itchar
INTEGER i,j,nx,ny,maxiter_S,itmax,ei,ej,it,maxiter,c
REAL*4 dx,dy,x0,xl,y0,yl,dv,Se,Sw,Sn,Ss,tolerance,residual,time,dt,ue,uw,vn,vs,Div,re,gamma,dbkx,dbky,time_count
REAL*4, ALLOCATABLE :: aE(:,:),aW(:,:),aN(:,:),aS(:,:),aP(:,:),Sp(:,:),u(:,:),u1(:,:),v(:,:),v1(:,:)
REAL*4, ALLOCATABLE :: x(:),xc(:),y(:),yc(:),P(:,:),Pp(:,:),d_h(:,:),d_v(:,:),uc(:,:),vc(:,:)

!Inicialización de variables
nx=80; ny=20
x0=0.0; xl=10.0
y0=0.0; yl=4.0
Re=100.0
time=10.
dt=0.005
maxiter_S=1000
maxiter=1000
tolerance=1e-3

ALLOCATE(aE(nx,ny),aW(nx,ny),aN(nx,ny),aS(nx,ny),aP(nx,ny),Sp(nx,ny),P(0:nx+1,0:ny+1),Pp(0:nx+1,0:ny+1))
ALLOCATE(u(0:nx,0:ny+1),u1(0:nx,0:ny+1),v(0:nx+1,0:ny),v1(0:nx+1,0:ny),x(0:nx),xc(0:nx+1),y(0:ny),yc(0:ny+1))
ALLOCATE(uc(0:nx+1,0:ny+1),vc(0:nx+1,0:ny+1),d_h(0:nx+1,0:ny+1),d_v(0:nx+1,0:ny+1))

gamma=1.0/Re
dx=(xl-x0)/FLOAT(nx)
dy=(yl-y0)/FLOAT(ny)
itmax=INT(time/dt)+1
dv=dx*dy
Se=dy; Sw=dy
Sn=dx; Ss=dx

CALL MESH_1D(nx,x0,xl,x,xc)
CALL MESH_1D(ny,y0,yl,y,yc)

!Condiciones de frontera e iniciales
!calovi
u=1.0; v=0.0; P=0.0; Pp=0.0

u(:,ny+1)=0.0;u(:,0)=0.0;
d_h=0.0; d_v=0.0

! u(0,:)=1.0		!Entrada de flujo


u1=u; v1=v		!u1 y v1 son al tiempo actual

! write(*,*) ny,UBOUND(u,2)-1
! archivo de animación
OPEN(1,FILE='./datos/anim.gnp',STATUS='REPLACE')
WRITE(1,*)'set size square; set xrange[0:10]; set yrange[0:1]; unset key'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!INICIA CICLO TEMPORAL
! 
DO it=1,itmax

  Div=1.0
  c=1
  time_count=it*dt-dt
  !INICIA CICLO SIMPLEC!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  DO WHILE ((c .LT. maxiter_S) .AND. (Div .GT. tolerance))
  aE=0.0; aW=0.0; aN=0.0; aS=0.0; aP=0.0; Sp=0.0
  
  !ECUACIÓN COMPONENTE U #######################################################################################
  ei=UBOUND(u,1)-1
  ej=UBOUND(u,2)-1
  
  !Cálculo de coeficientes
  DO i=1,ei
    DO j=1,ej
      !Flujos en las caras
      ue=(u(i,j)+u(i+1,j))*0.5
      uw=(u(i,j)+u(i-1,j))*0.5
      vn=(v(i,j)+v(i+1,j))*0.5
      vs=(v(i,j-1)+v(i+1,j-1))*0.5
      
      aE(i,j)=gamma*Se/dx-0.5*ue*Se
      aW(i,j)=gamma*Sw/dx+0.5*uw*Sw
      aN(i,j)=gamma*Sn/dy-0.5*vn*Sn
      aS(i,j)=gamma*Ss/dy+0.5*vs*Ss
      aP(i,j)=aE(i,j)+aW(i,j)+aN(i,j)+aS(i,j)+dv/dt
      Sp(i,j)=u(i,j)*dv/dt-(P(i+1,j)-P(i,j))*dv/dx
    END DO
  END DO
  
  !Corrección por condiciones de frontera
  dbkx=0.0; dbky=1.0
  
  !ESTE
  aP(ei,1:ej)=aP(ei,1:ej)-aE(ei,1:ej) !condición Neumann
  !SP(ei,1:ej)=SP(ei,1:ej)+(1.0+dbkx)*aE(ei,1:ej)*u(ei+1,1:ej)
  aE(ei,1:ej)=0.0
  
  !OESTE
  aP(1,1:ej)=aP(1,1:ej)+dbkx*aW(1,1:ej)
  Sp(1,1:ej)=Sp(1,1:ej)+(1+dbkx)*aw(1,1:ej)*u(0,1:ej)
  aw(1,1:ej)=0.0
  
  !NORTE
  aP(1:ei,ej)=aP(1:ei,ej)+dbky*aN(1:ei,ej)
  Sp(1:ei,ej)=Sp(1:ei,ej)+(1+dbky)*aN(1:ei,ej)*u(1:ei,ej+1)
  aN(1:ei,ej)=0.0
  
  !SUR
  aP(1:ei,1)=aP(1:ei,1)+dbky*aS(1:ei,1)
  Sp(1:ei,1)=Sp(1:ei,1)+(1+dbky)*aS(1:ei,1)*u(1:ei,0)
  aS(1:ei,1)=0.0
  
  !Cálculo de coeficientes para corrección de la presión
  d_h(1:ei,1:ej)=Se/(aP(1:ei,1:ej)-(aE(1:ei,1:ej)+aW(1:ei,1:ej)+aN(1:ei,1:ej)+aS(1:ei,1:ej)))
  
  !Solución del sistema de ecuaciones
  CALL Gauss_TDMA2D(u1,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny,maxiter,tolerance,residual)
!###########################################################################################################
  aP=0.0; aE=0.0; aW=0.0; aN=0.0; aS=0.0; Sp=0.0 ;
  !ECUACIÓN COMPONENTE V

  ei=UBOUND(v,1)-1
  ej=UBOUND(v,2)-1
  
  !Cálculo de coeficientes
  DO i=1,ei
    DO j=1,ej
      !Flujos en las caras
      ue=(u(i,j)+u(i,j+1))*0.5
      uw=(u(i-1,j)+u(i-1,j+1))*0.5
      vn=(v(i,j)+v(i,j+1))*0.5
      vs=(v(i,j)+v(i,j-1))*0.5
      
      aE(i,j)=gamma*Se/dx-0.5*ue*Se
      aW(i,j)=gamma*Sw/dx+0.5*uw*Sw
      aN(i,j)=gamma*Sn/dy-0.5*vn*Sn
      aS(i,j)=gamma*Ss/dy+0.5*vs*Ss
      aP(i,j)=aE(i,j)+aW(i,j)+aN(i,j)+aS(i,j)+dv/dt
      Sp(i,j)=v(i,j)*dv/dt-(P(i,j+1)-P(i,j))*dv/dy
    END DO
  END DO
  
  !Corrección por condiciones de frontera
  dbkx=1.0; dbky=0.0
  
  !ESTE
  aP(ei,1:ej)=aP(ei,1:ej)-aE(ei,1:ej)
  !SP(ei,1:ej)=SP(ei,1:ej)+(1.0+dbkx)*aE(ei,1:ej)*v(ei+1,1:ej)
  aE(ei,1:ej)=0.0
  
  !OESTE
  aP(1,1:ej)=aP(1,1:ej)+dbkx*aW(1,1:ej)
  Sp(1,1:ej)=Sp(1,1:ej)+(1+dbkx)*aW(1,1:ej)*v(0,1:ej)
  aw(1,1:ej)=0.0
  
  !NORTE
  aP(1:ei,ej)=aP(1:ei,ej)+dbky*aN(1:ei,ej)
  Sp(1:ei,ej)=Sp(1:ei,ej)+(1+dbky)*aN(1:ei,ej)*v(1:ei,ej+1)
  aN(1:ei,ej)=0.0
  
  !SUR
  aP(1:ei,1)=aP(1:ei,1)+dbky*aS(1:ei,1)
  Sp(1:ei,1)=Sp(1:ei,1)+(1+dbky)*aS(1:ei,1)*v(1:ei,0)
  aS(1:ei,1)=0.0
  
  !Cálculo de coeficientes para corrección de la presión
  d_v(1:ei,1:ej)=Sn/(aP(1:ei,1:ej)-(aE(1:ei,1:ej)+aW(1:ei,1:ej)+aN(1:ei,1:ej)+aS(1:ei,1:ej)))
  
  !Solución del sistema de ecuaciones
  CALL Gauss_TDMA2D(v1,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny,maxiter,tolerance,residual)
  
!##########################################################################################################

  !CORRECIÓN DE LA PRESIÓN
aP=0.0; aE=0.0; aW=0.0; aN=0.0; aS=0.0; Sp=0.0; Pp=0.0
  ei=UBOUND(Pp,1)-1
  ej=UBOUND(Pp,2)-1
  
  !Cálculo de coeficientes
  DO i=1,ei
    DO j=1,ej
      !Flujos en las caras
      ue=u1(i,j)
      uw=u1(i-1,j)
      vn=v1(i,j)
      vs=v1(i,j-1)
      
      aE(i,j)=d_h(i,j)*Se
      aW(i,j)=d_h(i-1,j)*Sw
      aN(i,j)=d_v(i,j)*Sn
      aS(i,j)=d_v(i,j-1)*Ss
      aP(i,j)=aE(i,j)+aW(i,j)+aN(i,j)+aS(i,j)
      Sp(i,j)=-ue*Se+uw*Sw-vn*Sn+vs*Ss
    END DO
  END DO
  
  !Solución del sistema de ecuaciones
  CALL Gauss_TDMA2D(Pp,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny,maxiter,tolerance,residual)
  
  !Correción de la presión
  P=P+Pp
  
  !CORRECIÓN DE VELOCIDADES
    !COMPONENTE U
    DO i=1, nx-1
      DO j=1,ny
	u1(i,j)=u1(i,j)+d_h(i,j)*(Pp(i,j)-Pp(i+1,j))
      END DO
    END DO

    ! calovi
    !corrección de velocidades en la dirección u para condición de frontera Newmann        
     u1(nx,:)=u1(nx-1,:)
     
    !COMPONENTE V
    DO i=1,nx
      DO j=1,ny-1
	v1(i,j)=v1(i,j)+d_v(i,j)*(Pp(i,j)-Pp(i,j+1))
      END DO
    END DO

    ! calovi
    !corrección de velocidades en la dirección u para condición de frontera Newmann        
    v1(nx+1,:)=v1(nx,:)
    
    
    !Para la divergencia y el criterio de paro
    Div=MAXVAL(ABS(Sp(1:ei,1:ej))) !Divergencia de u=0, ecuación de continuidad
    c=c+1 !Aumenta el contador de iteraciones SIMPLEC
    
  END DO !TERMINA CICLO SIMPLEC &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  
  u=u1; v=v1
  
  WRITE(*,*)it,c,Div, time_count
  
  !Escritura de campo de velocidades 
  IF (MOD(it,100) .EQ. 0) THEN
    CALL interpolateToNodesUs(uc,u,nx,ny)
    CALL interpolateToNodesVs(vc,v,nx,ny)
  
   CALL WriteVectorField('./datos/vel',it,uc,vc,xc,yc,nx,ny)
  
  !para el archivo de animación
    WRITE(itchar,'(i6)')it
    itchar=ADJUSTL(itchar)
    itchar=itchar(1:LEN_TRIM(itchar))
    WRITE(1,*)"p 'vel"//itchar(1:LEN_TRIM(itchar))//".txt' u 1:2:(0.25*$3):(0.25*$4) w vec"
    WRITE(1,*)'pause 0.5'
    
   !archivo propio para el tiempo indicado
    OPEN(5,FILE="./datos/"//itchar(1:LEN_TRIM(itchar))//".gp",STATUS='REPLACE')
    WRITE(5,*)'set size square; set xrange[0:10]; set yrange[0:1]; unset key'
    WRITE(5,*)"p 'vel"//itchar(1:LEN_TRIM(itchar))//".txt' u 1:2:(0.2*$3):(0.2*$4) w vec"
    
   !archivo para paraview
   !open(10,file='./datos_paraview/data'//itchar(1:LEN_TRIM(itchar))//'.dat')       
   !write(10,*) 'TITLE= "TESTPLOT" '
   !write(10,*) 'VARIABLES="x","y","fx","fy","magf" '              !
   !write(10,*) 'ZONE T="1", I=',nx+1,', J=',ny+1

   !do i=0,nx+1
   !do j=0,ny+1
   !write(10,*)xc(i),yc(j),uc(i,j),vc(i,j),SQRT(uc(i,j)**2+vc(i,j)**2)
   !end do
   !end do
  END IF
  
END DO !TERMINA CICLO TEMPORAL

!      CALL interpolateToNodesUs(uc,u,nx,ny)
!      CALL interpolateToNodesVs(vc,v,nx,ny)
!    
!      CALL WriteVectorField('vel',0,uc,vc,xc,yc,nx,ny)

CLOSE(1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !Solucion analitica estacionaria
   OPEN(15,FILE="Sol_analitica.txt",STATUS='REPLACE')
   do i=0,0
   do j=0,ny+1
   if (yc(j) .ge. (yl/2.)) then
   write(15,*)xc(i),yc(j),1.*(1.-((yc(j)-(yl/2.))**2/yl**2)),0.
   else
   write(15,*)xc(i),yc(j),1.*(1.-(((yl/2.)-yc(j))**2/yl**2)),0.
   end if
   end do
   end do
  
END PROGRAM

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
SUBROUTINE MESH_1D(nx,x0,xl,x,xc)
IMPLICIT NONE
INTEGER i, nx
REAL*4 x0,xl,dx
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

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

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

!***************************************************************

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

!****************************************************************

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

!************************************************************************************************************

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
	write(11,*)xc(i),yc(j),uc(i,j),vc(i,j),SQRT(uc(i,j)**2+vc(i,j)**2)
	end do
	write(11,*)''
	end do
close(11)
End Subroutine
