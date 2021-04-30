PROGRAM AdvectionGeneration

REAL dx, dv, L,Se,Sw,Sf,x0,xl   !k=constante,dx=espacio en x, dv=diferencial de volumen,L=longitud de barra,Se=Superficie este,Sw=Superficie oeste,x0=punto inicial, xl=punto final, Sf= Area asignada a la fuente
REAL Tp, Tinf, k  !Tp=temperatura extremo del principio, Tf=temperatura extremo final
REAL tolerance, residual   !para gauss-seidel
INTEGER nx,i,m   !nx=numero de volumenes, i=variable de iteracion, m=variable de iteración en matriz
INTEGER max_iter, n
REAL, allocatable:: T(:),Ta(:),aP(:),aE(:),aW(:),sP(:),xc(:),x(:)     !T(:)=Temperatura de puntos P por metodo numerico,Ta(:)=Temperatura de punto P por analitica, aP(:)=coeficiente de P, aE(:)=coeficiente este, aW(:)=coeficeiente oeste, sP(:)=término fuente, xc(:)=centros de volumenes de control, x(:)=caras de volumenes de control
REAL, allocatable:: phi(:,:), phi2(:,:)
REAL, allocatable:: A(:,:), phi_v(:), sp_v(:) !Matriz A


!Definir el tamaño de malla
nx=20; x0=0.0; xl=1.0
L=xl-x0; dx = L/float(nx)

!Tolerancias
max_iter=1000
tolerance=1e-10
n=nx


!Condiciones de Temperatura
Tp= 100; Tinf=20

!Definiendo constantes
k=5 

!Definiendo el tamaño de los arreglos
allocate(T(0:nx+1),Ta(0:nx+1),aP(nx),aE(nx),aW(nx),sP(nx),xc(0:nx+1),x(0:nx))
allocate(A(nx,nx), phi_v(0:nx+1), sp_v(nx))

!Llamando a la subroutine Mesh1D
call Mesh1D(xc,x,x0,xl,nx)

!Precolocar los arreglos
aP=0; aE=0; aW=0; sP=0; T=0; Ta=0

!Asignando las condiciones de frontera
T(0)=Tp; T(nx+1)=0

!Definiendo las áreas
Se=1.0; Sw=1.0; Sf=1.0

!Determinando los coeficientes
do i = 1, nx
	aE(i) = 1/dx
	aW(i) = 1/dx
	aP(i) = aE(i) + aW(i)+k*k*dx
	sP(i) = k*k*Tinf*dx
enddo

!Corección de las condiciones a la frontera

!Cara este
aP(nx) = aP(nx) - aE(nx)
!sP(nx) = sP(nx)   !El valor fuente más dos veces el coeficiente de la cara por la temperatura de esa frontera
aE(nx) = 0.0

!Cara oeste
aP(1) = aP(1) + aW(1)
sP(1) = sP(1) + 2.0*aW(1)*T(0)
aW(1) = 0.0

!Muestra de coeficientes
do i = 1,nx
	write(*,*) aW(i),aE(i),sP(i),aP(i)
enddo

!Solucion por Guss-Seidel
!Creación de la matriz
	A=0.0
	phi_v=T
	
	do j=1, nx
		do i=1, nx
		if(i==j) then
		A(i,i)=ap(i)
		endif
		sp_v(i)=sp(i)
		end do
	end do

		!ae
	do i=1, nx-1
			A(i+1,i)=-aE(i)
		enddo
	
	!aw
		do i=2, nx
			A(i-1,i)=-aW(i)
		enddo
! 	write(*,*) A
! 	write(*,*) sp_v
call GS(max_iter, n, tolerance, phi_v, sp_v, A)
T = phi_v

!Correciòn de valor en la frontera
T(nx+1)=T(nx)

do i = 1,nx
	write(*,*) T(i)
enddo

!Agregando la solución analítica
!Ta(0) = Tp; Ta(nx+1) = Tf
do i=0,nx+1
	Ta(i) = (Tp-Tinf)*(cosh(k*(L-xc(i))))/(cosh(k*L))+Tinf
enddo


open(1,file='adv-gen-bond-seidel.txt', status='replace')
do i=0,nx+1
	write(1,*) xc(i), T(i), Ta(i)
enddo

ENDPROGRAM AdvectionGeneration


!--------------------------------------------
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

!----------------------------------------------
Subroutine GS(max_iter, n, tolerance, phi_v, sp_v, A)
	implicit none
	integer max_iter,n, i,j,k
	real tolerance, error, x_aux
	real  phi_v(0:n+1), sp_v(n), A(n,n), phi_vold(n)
	
	k=1
	error=1.0
	
	do i=1, n
		phi_v(i)=sp_v(i)/A(i,i)
	enddo

	do while ((k .LT. max_iter) .AND. (error .GT. tolerance))
		error=0.0
		do i=1, n
			x_aux=0.0
			do j=1, n
				if (i .NE. j)then 
					x_aux=A(i,j)*phi_v(j)+x_aux
				endif
			enddo
			phi_v(i) = (sp_v(i)-x_aux)/A(i,i)
			error = error+ abs(phi_v(i)-phi_vold(i))
			phi_vold(i) = phi_v(i)
		enddo
		error=error/float(n)
! 		write(*,*) error,k	
		k=k+1
	enddo
	
End subroutine




