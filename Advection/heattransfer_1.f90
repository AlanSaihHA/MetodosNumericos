PROGRAM ADVECTION

REAL k, dx,dv,L, Se,Sw, x0,xl
REAL Th, Tc
INTEGER nx, i,m
REAL,allocatable::T(:),Ta(:),aP(:),aE(:),aW(:),sP(:),xc(:),x(:)

!Definir el tamaño de malla
nx=100; x0=0.0; xl=0.5;
L = xl-x0; dx = L /float(nx)

!Condiciones de frontera
Th = 500; Tc = 100

!Definiendo constantes
k = 1000

! Definiendo el tamaño de los arreglos
allocate(T(0:nx+1),Ta(0:nx+1),aP(nx),aE(nx),aW(nx),sP(nx),xc(0:nx+1),x(0:nx))

! Llamando la Subroutine Mesh1D

call Mesh1D(xc,x,x0,xl,nx)

!Precolar los arreglos
aP=0; aE=0; aW=0; sP=0; T=0; Ta=0

!Asignando las Condiciones de frontera
T(0) = Tc; T(nx+1) = Th

!Definiendo las áreas
Se = 0.010; Sw = 0.010

!Determinando los coeficientes

do i = 1,nx
    aE(i) = k*Se/dx
    aW(i) = k*Sw/dx
    aP(i) = aE(i) + aW(i)
    sP(i)= 0.0
enddo

!correción de las condiciones a la frontera

!Cara Este 
aP(nx) = aP(nx) + aE(nx)
sP(nx) = sP(nx) + 2.0*aE(nx)*T(nx+1)
aE(nx) = 0.0

!Cara Oeste 
aP(1) = aP(1) + aW(1)
sP(1) = sP(1) + 2.0*aW(1)*T(0)
aW(1) = 0.0

! do i = 1,nx
!     write(*,*) aW(i),aE(i),sP(i), ap(i)
! enddo   

do m = 1,100000
    do i = 1,nx
    T (i) = (aE(i)*T(i+1) + aW(i)*T(i-1)+sP(i))/aP(i)
    enddo
enddo

! do i = 1,nx
!     write(*,*) T(i)
! enddo

!Agregando la solución analítica

Ta(0) = Tc; Ta(nx+1) = Th
do i=0,nx+1
    Ta(i) = ((Th-Tc)/L)*xc(i) + Tc
enddo

open(1,file='temp.txt',status='replace')
do i=0,nx+1
    write(1,*) xc(i),T(i),Ta(i)
enddo
ENDPROGRAM ADVECTION

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
