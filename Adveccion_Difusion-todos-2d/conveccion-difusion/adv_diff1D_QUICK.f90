Program adv_diff1D
implicit none
integer nx, i, j
real*4 dx, dv, T0, T1, Se, Sw, Pe, u, gama, L,rho, x0,xl
real*4, allocatable::T(:),Ta(:),aP(:),aW(:),aWW(:),aE(:),Sp(:),xc(:),x(:)

!Inicialización

nx=25; x0=0.0; xl=1;
L = xl-x0; dx = L /float(nx)

u = 0.2;T0 = 1.0; T1 = 0.0;gama = 0.1;rho=1.0

!Alojamiento de memoria dinámica
allocate(T(0:nx+1),Ta(0:nx+1), aP(nx), aE(nx), aW(nx), aWW(nx), Sp(nx),xc(0:nx+1),x(0:nx))

!Áreas
Se = 1.0; Sw = 1.0

!precolocando los vectores
  aP=0.;aE=0.;aW=0.;aWW=0.; SP=0.; T=0.


!Valores inciales en memoria
T(0) = T0; T(nx+1) = T1

call Mesh1D(xc,x,x0,xl,nx)

!Solución numérica
do i=2,nx-1
!Esquema Central**
    if (i ==2) then 
    aE(i) = -(3./8.)*u*rho*Se + gama*Se/dx
    aW(i) = (1./8.)*u*rho*Se + (7./8.)*u*rho*Sw + gama*Sw/dx 
    aP(i) = aE(i) + aW(i) - 10./8 *u*rho*Sw+u*rho*Se
    Sp(i) = -(2./8.)*u*rho*Sw*T(0)
    else    
    aE(i) = -(3./8.)*u*rho*Se + gama*Se/dx
    aW(i) = (1./8.)*u*rho*Se + gama*Sw/dx + (6./8.)*u*rho*Sw
    aWW(i) = -(1./8.)*u*rho*Sw    
    aP(i) = aE(i) + aW(i) + aWW(i)
    Sp(i) = 0.0
    
    endif

end do

! !Condiciones de frontera
! !OEste
aE(1)=-(3./8.)*u*rho*Se +  gama*Se/dx + (1./3.)*gama*Se/dx
aP(1) = aE(1) + (10./8.)*u*rho*Se +(8./3.)*gama*Se/dx
Sp(1) = ((2./8.)*u*rho*Se + u*rho*Sw + (8./3.)*gama*Se/dx)*T(0)
aW(1) = 0.0
aWW(1) = 0.0

! !Este
aW(nx) = (6./8.)*u*rho*Sw + gama*Sw/dx + gama*Sw/dx/3.
aWW (nx)= -(1./8.)*u*rho*Sw
aP(nx) = aW(nx) + aWW(nx) - u*rho*Sw + (8./3.)*gama*Sw/dx
sP(nx) =- u*rho*Se*T(nx+1)
aE(nx) = 0.0

! write(*,*)aWW

! write(*,*) aWW,aW,aE,sP,aP
write(*,*)aWW,aW,aE,aP,sP
! !Método de Jacobi
do j=1,1000
	do i=1,nx
		T(i) = (aE(i)*T(i+1) + aW(i)*T(i-1) + aWW(i)*T(i-2) + Sp(i)) / aP(i)
	end do
end do
! 
!Solución analítica 
do i=0,nx+1
    Ta(i) = (T1-T0)*( (exp(rho*u*xc(i)/gama)-1)/(exp(rho*u*L/gama)-1) ) + T0
enddo
!Escritura de resultados
open(1,file="temp.txt", status="replace")
do i=0,nx+1
	write(1,*)xc(i),T(i),Ta(i)
end do
close(1)
end program
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
!!****************************************************************************************************************

