PROGRAM problema1

REAL k1,k2, dx, dv, L,Se,Sw,Sf,x0,xl   !k=constante,dx=espacio en x, dv=diferencial de volumen,L=longitud de barra,Se=Superficie este,Sw=Superficie oeste,x0=punto inicial, xl=punto final, Sf= Area asignada a la fuente
REAL Tp, Tf, q  !Tp=temperatura extremo del principio, Tf=temperatura extremo final
INTEGER nx,i,m   !nx=numero de volumenes, i=variable de iteracion, m=variable de iteración en matriz
REAL, allocatable:: T(:),Ta(:),aP(:),aE(:),aW(:),sP(:),xc(:),x(:)    !T(:)=Temperatura de puntos P por metodo numerico,Ta(:)=Temperatura de punto P por analitica, aP(:)=coeficiente de P, aE(:)=coeficiente este, aW(:)=coeficeiente oeste, sP(:)=término fuente, xc(:)=centros de volumenes de control, x(:)=caras de volumenes de control



!Definir el tamaño de malla
nx=11; x0=0.0; xl=0.02;
L=xl-x0; dx = L/float(nx)

!Condiciones de frontera
Tp= -1; Tf=0.5

!Definiendo generación de calor
q=100000  !W/m³

!Definiedo la constante k
k1=1
k2=10

!Definiendo el tamaño de los arreglos
allocate(T(0:nx+1),Ta(0:nx+1),aP(nx),aE(nx),aW(nx),sP(nx),xc(0:nx+1),x(0:nx))

!Llamando a la subroutine Mesh1D
call Mesh1D(xc,x,x0,xl,nx)

!Precolocar los arreglos
aP=0; aE=0; aW=0; sP=0; T=0; Ta=0

!Asignando las condiciones de frontera
T(0)=Tp; T(nx+1)=Tf

!Definiendo las áreas
Se=1.0; Sw=1.0; Sf=1.0

!Determinando los coeficientes
do i = 1, nx/2-1
	aE(i) = k1*Se/dx
	aW(i) = k1*Sw/dx
	aP(i) = aE(i) + aW(i)
	sP(i) = q*Sf*dx
enddo

do i = nx/2+1, nx
	aE(i) = k2*Se/dx
	aW(i) = k2*Sw/dx
	aP(i) = aE(i) + aW(i)
	sP(i) = q*Sf*dx
enddo

!Corrigiendo union
aW(nx/2)=k1*Sw/dx
aE((nx/2)=k2*Se/dx
Sp(nx/2))=aW(nx/2)+aE((nx/2)



!Corección de las condiciones a la frontera

!Cara este
aP(nx) = aP(nx) + aE(nx)
sP(nx) = sP(nx) + 2.0*aE(nx)*T(nx+1)   !El valor fuente más dos veces el coeficiente de la cara por la temperatura de esa frontera
aE(nx) = 0.0

!Cara oeste
aP(1) = aP(1) + aW(1)
sP(1) = sP(1) + 2.0*aW(1)*T(0)
aW(1) = 0.0

!Muestra de coeficientes
do i = 1,nx
	write(*,*) aW(i),aE(i),sP(i),aP(i)
enddo

!Solucion de matriz por metodo de Jacobi
do m = 1,1000000
	do i = 1,nx
	t (i) = (aE(i)*T(i+1) + aW(i)*T(i-1)+sP(i))/aP(i)
	enddo
enddo

!Mostrar resultados numericos
do i = 0,nx+1
	write(*,*) T(i)
enddo

!Agregando la solución analítica
Ta(0) = Tp; Ta(nx+1) = Tf
do i=0,nx/2
	Ta(i) = (k1*Tp+k2*Tf)/(k1+k2) & 
	+((Tf-Tp)/(L/2))*(k2/(k1+k2))*(xc(i)-(L/2)) &
	-(q/2)*(xc(i)-L/2+L/2)*((xc(i)-L/2)/k1-((2*(L/2))/(k1+k2)))
enddo

do i=nx/2+1,nx+1
	Ta(i) = (k1*Tp+k2*Tf)/(k1+k2) &
	+((Tf-Tp)/(L/2))*(k1/(k1+k2))*(xc(i)-(L/2)) &
	+(q/2)*(L/2-(xc(i)-L/2))*((xc(i)-L/2)/k2+((2*(L/2))/(k1+k2)))
enddo


open(1,file='problema1-jacobi.txt', status='replace')
do i=0,nx+1
	write(1,*) xc(i), T(i), Ta(i)
enddo

ENDPROGRAM problema1


!-----------------------------------------------------------
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

