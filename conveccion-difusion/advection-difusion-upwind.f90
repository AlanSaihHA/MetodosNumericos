PROGRAM ADV_DIFF

REAL k,dx, dv, L,Se,Sw,x0,xl   !k=constante,dx=espacio en x, dv=diferencial de volumen,L=longitud de barra,Se=Superficie este,Sw=Superficie oeste,x0=punto inicial, xl=punto final
REAL Th, Tc, u, gama, rho,T0,T1   !Th=temperatura de extremo caliente, Tc=temperatura de extremo fría
INTEGER nx,i,m   !nx=numero de volumenes, i=variable de iteracion, m=variable de iteración en matriz
REAL, allocatable:: T(:),Ta(:),aP(:),aE(:),aW(:),sP(:),xc(:),x(:)     !T(:)=Temperatura de puntos P por metodo numerico, Ta(:)=Temperatura de punto P por analitica, aP(:)=coeficiente de P, aE(:)=coeficiente este, aW(:)=coeficeiente oeste, sP(:)=término fuente, xc(:)=centros de volumenes de control, x(:)=caras de volumenes de control


!Definir el tamaño de malla
nx=25; x0=0.0; xl=1;
L=xl-x0; dx = L/float(nx)

!Definiedo constantees
u=2.5; T0=1.0; T1=0.0; gama=0.1; rho=1.0

!Definiendo el tamaño de los arreglos
allocate(T(0:nx+1),Ta(0:nx+1),aP(nx),aE(nx),aW(nx),sP(nx),xc(0:nx+1),x(0:nx))

!Llamando a la subroutine Mesh1D
call Mesh1D(xc,x,x0,xl,nx)

!Precolocar los arreglos
aP=0; aE=0; aW=0; sP=0; T=0; Ta=0

!Asignando las condiciones de frontera
T(0)=T0; T(nx+1)=T1

!Definiendo las áreas
Se=1.0; Sw=1.0

!Determinando los coeficientes (Esquema Upwind)
do i = 1, nx
	aW(i) = gama*Sw/dx + max(u*Sw*rho,0.0)
	aE(i) = gama*Se/dx - min(u*Se*rho,0.0)
	aP(i) = aE(i) + aW(i)
	sP(i) = 0.0
enddo

!Corección de las condiciones a la frontera

!Cara este
aP(nx) = aP(nx) + aE(nx)
sP(nx) = sP(nx) + 2.0*aE(nx)*T(nx+1)   !El valor fuente más dos veces el coeficiente de la cara por la temperatura de esa frontera
aE(nx) = 0.0

!Cara oeste
aP(1) = aP(1) + aW(1)
sP(1) = sP(1) + 2.0*aW(1)*T(0)
aW(1) = 0.0

do i = 1,nx
	write(*,*) aW(i),aE(i),sP(i),aP(i)
enddo

do m = 1,1000
	do i = 1,nx
	t (i) = (aE(i)*T(i+1) + aW(i)*T(i-1)+sP(i))/aP(i)
	enddo
enddo

do i = 1,nx
	write(*,*) T(i)
enddo

!Solución analítica 
do i=0,nx+1
    Ta(i) = (T1-T0)*( (exp(rho*u*xc(i)/gama)-1)/(exp(rho*u*L/gama)-1) ) + T0
enddo


open(1,file='temp.txt', status='replace')
do i=0,nx+1
	write(1,*) xc(i), T(i), Ta(i)
enddo
close(1)

ENDPROGRAM ADV_DIFF



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
