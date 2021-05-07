PROGRAM AdvectionGeneration

REAL dx, dv, L,Se,Sw,Sf,x0,xl   !k=constante,dx=espacio en x, dv=diferencial de volumen,L=longitud de barra,Se=Superficie este,Sw=Superficie oeste,x0=punto inicial, xl=punto final, Sf= Area asignada a la fuente
REAL Tp, Tinf, k,h ,q !Tp=temperatura extremo del principio, Tf=temperatura extremo final
REAL tolerance, residual   !para gauss-seidel
INTEGER nx,i,m   !nx=numero de volumenes, i=variable de iteracion, m=variable de iteración en matriz
INTEGER max_iter, n
REAL, allocatable:: T(:),Ta(:),aP(:),aE(:),aW(:),sP(:),xc(:),x(:)     !T(:)=Temperatura de puntos P por metodo numerico,Ta(:)=Temperatura de punto P por analitica, aP(:)=coeficiente de P, aE(:)=coeficiente este, aW(:)=coeficeiente oeste, sP(:)=término fuente, xc(:)=centros de volumenes de control, x(:)=caras de volumenes de control
REAL, allocatable:: phi(:,:), phi2(:,:)
REAL, allocatable:: A(:,:), phi_v(:), sp_v(:) !Matriz A


!Definir el tamaño de malla
nx=30; x0=0.0; xl=0.2
L=xl-x0; dx = L/float(nx)

!Tolerancias
max_iter=1000
tolerance=1e-10
n=nx

!Condiciones de Temperatura
Tp= 300; Tinf=100

!Definiendo constantes
k=2; h=10; q=1.5e5

!Definiendo el tamaño de los arreglos
allocate(T(0:nx+1),Ta(0:nx+1),aP(nx),aE(nx),aW(nx),sP(nx),xc(0:nx+1),x(0:nx))
allocate(A(nx,nx), phi_v(0:nx+1), sp_v(nx))

!Llamando a la subroutine Mesh1D
call Mesh1D(xc,x,x0,xl,nx)

!Precolocar los arreglos
aP=0; aE=0; aW=0; sP=0; T=0; Ta=0

!Asignando las condiciones de frontera, T(nx+1) se corrige posteriormente
T(0)=Tp; T(nx+1)=0    

!Definiendo las áreas
Se=1.0; Sw=1.0; Sf=1.0

!Determinando los coeficientes
do i = 1, nx
	aE(i) = k*Se/dx
	aW(i) = k*Sw/dx
	aP(i) = aE(i) + aW(i)
	sP(i) = q*Sf*dx
enddo

!Corección de las condiciones a la frontera

!Cara este
aP(nx) = aP(nx) + aE(nx)*((h*dx)/k-1)
sP(nx) = sP(nx) + aE(nx)*((h*Tinf*dx)/k)
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

!Correciòn de la frontera
T(nx+1)=T(nx)


!Agregando la solución analítica
!Ta(0) = Tp; Ta(nx+1) = Tf
do i=0,nx+1
	Ta(i) = Tp-(q/(2*k))*(xc(i)*xc(i))-((Tp-Tinf-q*L*(L/(2*k)+1/h))/((k/h)+L))*xc(i)
enddo

do i = 0,nx+1
	write(*,*) T(i), Ta(i)
enddo

open(1,file='problema2-seidel-3.txt', status='replace')
do i=0,nx+1
	write(1,*) xc(i), T(i), Ta(i)
enddo

ENDPROGRAM AdvectionGeneration


!-----------------------------------------------
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