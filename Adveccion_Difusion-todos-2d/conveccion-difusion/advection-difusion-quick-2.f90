PROGRAM ADV_DIFF

REAL k,dx, dv, L,Se,Sw,x0,xl   !k=constante,dx=espacio en x, dv=diferencial de volumen,L=longitud de barra,Se=Superficie este,Sw=Superficie oeste,x0=punto inicial, xl=punto final
REAL Th, Tc, u,ue,uw, gama, rho,T0,T1,alfae,alfaw   !Th=temperatura de extremo caliente, Tc=temperatura de extremo fría
INTEGER nx,i,m   !nx=numero de volumenes, i=variable de iteracion, m=variable de iteración en matriz
REAL, allocatable:: T(:),Ta(:),aP(:),aE(:),aW(:),aWW(:),aEE(:) ,sP(:),xc(:),x(:)     !T(:)=Temperatura de puntos P por metodo numerico, Ta(:)=Temperatura de punto P por analitica, aP(:)=coeficiente de P, aE(:)=coeficiente este, aW(:)=coeficeiente oeste, sP(:)=término fuente, xc(:)=centros de volumenes de control, x(:)=caras de volumenes de control


!Definir el tamaño de malla
nx=25; x0=0.0; xl=1;
L=xl-x0; dx = L/float(nx)

!Definiedo constantees
u=0.2; uw=0.2; ue=0.2; T0=1.0; T1=0.0; gama=0.1; rho=1.0

!Definiendo el tamaño de los arreglos
allocate(T(0:nx+1),Ta(0:nx+1),aP(nx),aE(nx),aW(nx),aWW(nx),aEE(nx) ,sP(nx),xc(0:nx+1),x(0:nx))

!Llamando a la subroutine Mesh1D
call Mesh1D(xc,x,x0,xl,nx)

!Precolocar los arreglos
aP=0; aE=0; aW=0; aWW=0; aEE=0; sP=0; T=0; Ta=0

!Asignando las condiciones de frontera
T(0)=T0; T(nx+1)=T1

!Definiendo las áreas
Se=1.0; Sw=1.0

!Determinando los coeficientes (Esquema QUICK)
do i = 1, nx
		if (uw .gt. 0.0 .and. ue .gt. 0.0) then
			 alfaw=1.0 ; alfae=1.0
		else 
			 alfaw=0.0 ; alfae=0.0
		end if
		
if (i==2 .and. i==nx-1) then
		!Corrigiendo el nodo 2 o 4
if (uw .gt. 0.0 .and. ue .gt. 0.0) then
aP(i) = (6.0/8.0)*ue - (3.0/8.0)*uw + gama*Se/dx + gama*Sw/dx
aW(i) = (1.0/8.0)*ue + (7.0/8.0)*uw + gama*Sw/dx
aE(i) = - (3.0/8.0)*ue + gama*Se/dx
sP(i) = - (2.0/8.0)*uw*T(0)
aWW(i) = 0.0
aEE(i) = 0.0
else
aP(i) = - (6.0/8.0)*uw + (3.0/8.0)*ue + gama*Se/dx + gama*Sw/dx
aE(i) = - (1.0/8.0)*uw - (7.0/8.0)*ue + gama*Sw/dx
aW(i) = (3.0/8.0)*uw + gama*Se/dx
sP(i) = (2.0/8.0)*ue*T(nx+1)
aWW(i) = 0.0
aEE(i) = 0.0
endif
else
		aW(i) = gama*Sw/dx + (6.0/8.0)*alfaw*uw &
		+ (1.0/8.0)*alfae*ue + (3.0/8.0)*(1.0-alfaw)*uw
		
		aE(i) = gama*Se/dx - (3.0/8.0)*alfae*ue &
		- (6.0/8.0)*(1.0-alfae)*ue - (1.0/8.0)*(1.0-alfaw)*uw
		
		aWW(i) = - (1.0/8.0)*alfaw*uw
		
		aEE(i) = (1.0/8.0)*(1.0-alfae)*ue

		
	aP(i) = aE(i) + aW(i) + aEE(i) + aWW(i)
	sP(i) = 0.0
endif	
enddo

!Corección de las condiciones a la frontera

!Cara oeste
if (uw .gt. 0.0 .and. ue .gt. 0.0) then		
aP(1) = gama*Sw/dx + 3.0*gama*Sw/dx + (7.0/8.0)*ue        
aE(1) = gama*Se/dx + (1.0/3.0)*gama*Sw/dx - (3.0/8.0)*ue  
sP(1) = ((8.0/3.0)*gama*Sw/dx + (2.0/8.0)*ue + uw)*T(0)    
aW(1) = 0.0
aWW(1) = 0.0
aEE(1) = 0.0
else 
aP(1) = gama*Se/dx - 3.0*gama*Sw/dx + (3.0/8.0)*ue        
aE(1) = - (6.0/8.0)*ue + gama*Se/dx - (1.0/3.0)*gama*Sw/dx  
sP(1) = (-(8.0/3.0)*gama*Sw/dx + uw)*T(0)    
aW(1) = 0.0
aWW(1) = 0.0
aEE(1) = (1.0/8.0)*ue
end if

!Cara este
if (uw .gt. 0.0 .and. ue .gt. 0.0) then
aP(nx) = 3.0*gama*Se/dx + gama*Sw/dx - (3.0/8.0)*uw
aW(nx) = (6.0/8.0)*uw + gama*Sw/dx + (1.0/3.0)*gama*Se/dx 
aWW(nx) = - (1.0/8.0)*uw                                     
sP(nx) = ((8.0/3.0)*gama*Se/dx-ue)*T(nx+1)
aE(nx) = 0.0
aEE(nx) = 0.0
else
aP(nx) = gama*Sw/dx - 3.0*gama*Sw/dx - (7.0/8.0)*uw        
aW(nx) = gama*Sw/dx - (1.0/3.0)*gama*Se/dx + (3.0/8.0)*uw  
sP(nx) = (-(8.0/3.0)*gama*Se/dx - (2.0/8.0)*uw - ue)*T(nx+1)    
aE(nx) = 0.0
aWW(nx) = 0.0
aEE(nx) = 0.0
endif



do i = 1,nx
	write(*,*) aW(i),aE(i),sP(i),aP(i)
enddo

do m = 1,1000
	do i = 1,nx
	t (i) = (aE(i)*T(i+1) + aW(i)*T(i-1) &
		+ aWW(i)*T(i-2) + aEE(i)*T(i+2) +sP(i))/aP(i)
	enddo
enddo

do i = 1,nx
	write(*,*) T(i)
enddo

!Solución analítica 
do i=0,nx+1
    Ta(i) = (T1-T0)*( (exp(rho*u*xc(i)/gama)-1)/(exp(rho*u*L/gama)-1) ) + T0
enddo


open(1,file='temp-quick.txt', status='replace')
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
