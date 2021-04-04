Program Laplace
  implicit none
  real k,Pe,dx,dy,Se,Sw,Sn,Ss,q,dv,ue,uw,alfaw,alfae,gama,rho,u
  real x0,xl,y0,yl,tolerance,residual, dt, time,T0,T1,f
  real pi
  integer nx,i,j,ny,max_iter,it,itmax
  real,  allocatable::T(:),Ta(:),aP(:),aE(:),aW(:),SP(:),xc(:),x(:)
  real, allocatable:: aWW(:),aEE(:)
  character*50 itchar, name
  !definiendo el tamaño de la malla
  !Ocupar mallas siempre de potencia de 2
  nx=50
  x0=0.0
  xl=1.5
  dx=(xl-x0)/float(nx)


  !definiendo las constantes
  rho=1.
  time=10.0
  dt=0.01
  itmax=int(time/dt)+1
  gama=0.03 !es la función gamma, que se lee en  gamma*S/delta
  max_iter=1000
  tolerance= 1e-4
  
  !Definiendo velocidades
  uw=2.
  ue=2.
  u=2.
  !Definiendo valores de frontera
  T0=0
  T1=0

  !definiendo el tamaño del vector
  allocate(T(0:nx+1),Ta(0:nx+1), aP(nx), aE(nx), aW(nx), aWW(nx),aEE(nx), Sp(nx),xc(0:nx+1),x(0:nx))

  !llamando la rubrutina que genera la malla
  call Mesh1D(xc,x,x0,xl,nx)

  !precolocando los vectores
  aP=0.;aE=0.;aW=0.; aWW=0.; aEE=0.; SP=0.; T=0.; dv=dx
    
  !determinado las condiciones de frontera, el manejo de los arreglos es parecido a Matlab 
  T(0)=T0
  T(nx+1)=T1
  

  !definiendo las áreas
  Se=1.0; Sw=1.0; 

  
  !determinando los coeficientes y agregando el paso temporal
do it=1,itmax
	do i = 1, nx
		
		    if (xc(i)<=0.6) then
    			f = (-200.*xc(i)+100.)*dv
    		    elseif(xc(i)>0.6 .and. xc(i)<=0.82) then
     			f = (100.*xc(i)-80.)*dv
     		    endif
				
		

		if (uw .gt. 0.0 .and. ue .gt. 0.0) then
			 alfaw=1.0 ; alfae=1.0
		else 
			 alfaw=0.0 ; alfae=0.0
		end if
		
		!segundo volumen o penultimo volumen
		if (uw .gt. 0.0 .and. ue .gt. 0.0 .and. i==2) then
		aP(i) = (6.0/8.0)*ue - (3.0/8.0)*uw + gama*Se/dx + gama*Sw/dx + dv/dt
		aW(i) = (1.0/8.0)*ue + (7.0/8.0)*uw + gama*Sw/dx
		aE(i) = - (3.0/8.0)*ue + gama*Se/dx
		sP(i) = - (2.0/8.0)*uw*T(0) + f
		aWW(i) = 0.0
		aEE(i) = 0.0
		else if (uw .lt. 0.0 .and. ue .lt. 0.0 .and. i==(nx-1)) then
		aP(i) = - (6.0/8.0)*uw + (3.0/8.0)*ue + gama*Se/dx + gama*Sw/dx
		aE(i) = - (1.0/8.0)*uw - (7.0/8.0)*ue + gama*Sw/dx
		aW(i) = (3.0/8.0)*uw + gama*Se/dx
		sP(i) = (2.0/8.0)*ue*T(nx+1) + f
		aWW(i) = 0.0
		aEE(i) = 0.0
		endif
		
		aW(i) = gama*Sw/dx + (6.0/8.0)*alfaw*uw &
		+ (1.0/8.0)*alfae*ue + (3.0/8.0)*(1.0-alfaw)*uw
		
		aE(i) = gama*Se/dx - (3.0/8.0)*alfae*ue &
		- (6.0/8.0)*(1.0-alfae)*ue - (1.0/8.0)*(1.0-alfaw)*uw
		
		aWW(i) = - (1.0/8.0)*alfaw*uw
		
		aEE(i) = (1.0/8.0)*(1.0-alfae)*ue
		
	aP(i) = aE(i) + aW(i) + aEE(i) + aWW(i) + dv/dt
	sP(i) = T(i)*dv/dt +f	
enddo 

!Correcion de valores 
!Cara oeste
if (uw .gt. 0.0 .and. ue .gt. 0.0) then		
aP(1) = gama*Sw/dx + 3.0*gama*Sw/dx + (7.0/8.0)*ue + dv/dt     !   
aE(1) = gama*Se/dx + (1.0/3.0)*gama*Sw/dx - (3.0/8.0)*ue     !
sP(1) = ((8.0/3.0)*gama*Sw/dx + (2.0/8.0)*ue + uw)*T(0) +T(0)*dv/dt + (-200.*xc(1)+100.)*dv
aW(1) = 0.0      !
aWW(1) = 0.0      !
aEE(1) = 0.0       !
else 
aP(1) = gama*Se/dx - 3.0*gama*Sw/dx + (3.0/8.0)*ue + dv/dt        
aE(1) = - (6.0/8.0)*ue + gama*Se/dx - (1.0/3.0)*gama*Sw/dx  
sP(1) = (-(8.0/3.0)*gama*Sw/dx + uw)*T(0)    
aW(1) = 0.0
aWW(1) = 0.0
aEE(1) = (1.0/8.0)*ue
end if

!Cara este
if (uw .gt. 0.0 .and. ue .gt. 0.0) then
aP(nx) = -(3./8.)*uw + gama*Sw/dx + dv/dt    !
aW(nx) = (6./8.)*uw + gama*Sw/dx     !
aWW(nx) = - (1.0/8.0)*uw              !                             
sP(nx) = T(nx+1)*dv/dt             !uw*T(nx+1) + sP(nx)
aE(nx) = 0.0                 !
aEE(nx) = 0.0                !
else
aP(nx) = gama*Sw/dx - 3.0*gama*Sw/dx - (7.0/8.0)*uw        
aW(nx) = gama*Sw/dx - (1.0/3.0)*gama*Se/dx + (3.0/8.0)*uw  
sP(nx) = (-(8.0/3.0)*gama*Se/dx - (2.0/8.0)*uw - ue)*T(nx+1)    
aE(nx) = 0.0
aWW(nx) = 0.0
aEE(nx) = 0.0
endif


! !Condiciones de frontera
! !OEste
!aE(1)=-(3./8.)*u*rho*Se +  gama*Sw/dx + (1./3.)*gama*Se/dx
!aP(1) = aE(1) + (10./8.)*u*rho*Se +(8./3.)*gama*Se/dx + dv/dt
!Sp(1) = ((2./8.)*u*rho*Se + u*rho*Sw + (8./3.)*gama*Se/dx)*T(0) + T(0)*dv/dt  + (-200.*xc(1)+100.)*dv
!aW(1) = 0.0
!aWW(1) = 0.0

! !Este
!aW(nx) = (6./8.)*u*rho*Sw + gama*Sw/dx !+ gama*Sw/dx/3.
!aWW (nx)= -(1./8.)*u*rho*Sw
!aP(nx) = aW(nx) + aWW(nx) + dv/dt - u*rho*Sw !+ (8./3.)*gama*Sw/dx
!sP(nx) =T(nx+1)*dv/dt  !+ (100.*xc(nx)-80.)*dv  !- u*rho*Se*T(nx+1) +(8./3.)*T(nx+1)*gama*Sw/dx
!aE(nx) = 0.0




 do i=1,nx
		T(i) = (aE(i)*T(i+1) + aW(i)*T(i-1) + aWW(i)*T(i-2) + Sp(i)) / aP(i)
 end do

 ! Condición de frontera Newmann dPhi/dx=phiE-phiP = 0 => PhiE=PhiP
  T(nx+1)=T(nx)
  
 if (mod(it,100) == 0)then 
    call  WriteScalarField2D('Temp2D',it,T,xc,nx)
!      write(itchar, '(i6)')it !Convertir a un entero a un caracter
!      itchar=adjustl(itchar)
!      name='Temp2D'//itchar(1:len_trim(itchar))//'.txt'
!      write(2,*)'sp "'//name(1:len_trim(name))//'" w pm3d'
!      write(2,*)"pause 0.5"
 endif

enddo
endProgram Laplace

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!********************************************************************
Subroutine WriteScalarField2D(Name,kx,T,xc,nx)
integer i,j,nx,kx
real*4 T(0:nx+1),xc(0:nx+1)
character*(*)Name
character*50 txt,Filename
write(txt,'(i6)')kx
txt=ADJUSTL(txt)
Filename=name//txt(1:len_trim(txt))//".txt"

open(10,file=Filename(1:len_trim(Filename)))
	do i=0,nx+1
	write(10,*)xc(i),T(i)
    end do
close(10)
End Subroutine
