Program adv_diff1D
implicit none
integer nx, i, j,it,itmax
real*4 dx, dv, T0, T1, Se, Sw, Pe, u, gama, L,rho, x0,xl
real*4, allocatable::T(:),Ta(:),aP(:),aW(:),aWW(:),aE(:),Sp(:),xc(:),x(:)
real*4 dt, time,f

!Inicialización

nx=46; x0=0.0; xl=1.5;
L = xl-x0; dx = L /float(nx)
! write(*,*)dx

u = 1.7;T0 = 0.0; T1 = 0.0;gama = 0.03;rho=1.0

time=10.0
dt=0.01
itmax=int(time/dt)+1
!Alojamiento de memoria dinámica
allocate(T(0:nx+1),Ta(0:nx+1), aP(nx), aE(nx), aW(nx), aWW(nx), Sp(nx),xc(0:nx+1),x(0:nx))

!Áreas
Se = 1.0; Sw = 1.0; dv=dx

!precolocando los vectores
aP=0.;aE=0.;aW=0.;aWW=0.; SP=0.; T=0.


!Valores inciales en memoria
T(0) = T0; !T(nx+1) = T1

call Mesh1D(xc,x,x0,xl,nx)

do it=1,itmax
!Solución numérica
do i=2,nx-1
!Esquema Central**
    if (xc(i)<=0.6) then
    f = (-200.*xc(i)+100.)*dv
    elseif(xc(i)>0.6 .and. xc(i)<=0.82) then
     f = (100.*xc(i)-80.)*dv
!      write(*,*)xc(i)
    endif
    if (i ==2) then 
    aE(i) = -(3./8.)*u*rho*Se + gama*Se/dx
    aW(i) = (1./8.)*u*rho*Se + (7./8.)*u*rho*Sw + gama*Sw/dx 
    aP(i) = aE(i) + aW(i) - 10./8 *u*rho*Sw + u*rho*Se + dv/dt
    Sp(i) = -(2./8.)*u*rho*Sw*T(0) + f
    else    
    aE(i) = -(3./8.)*u*rho*Se + gama*Se/dx
    aW(i) = (1./8.)*u*rho*Se + gama*Sw/dx + (6./8.)*u*rho*Sw
    aWW(i) = -(1./8.)*u*rho*Sw    
    aP(i) = aE(i) + aW(i) + aWW(i) + dv/dt
    Sp(i) = T(i)*dv/dt +f
    
    endif

end do

! !Condiciones de frontera
! !OEste
aE(1)=-(3./8.)*u*rho*Se +  gama*Sw/dx + (1./3.)*gama*Se/dx
aP(1) = aE(1) + (10./8.)*u*rho*Se +(8./3.)*gama*Se/dx + dv/dt
Sp(1) = ((2./8.)*u*rho*Se + u*rho*Sw + (8./3.)*gama*Se/dx)*T(0) + T(0)*dv/dt  + (-200.*xc(1)+100.)*dv
aW(1) = 0.0
aWW(1) = 0.0

! !Este
aW(nx) = (6./8.)*u*rho*Sw + gama*Sw/dx !+ gama*Sw/dx/3.
aWW (nx)= -(1./8.)*u*rho*Sw
aP(nx) = aW(nx) + aWW(nx) + dv/dt - u*rho*Sw !+ (8./3.)*gama*Sw/dx
sP(nx) =T(nx+1)*dv/dt  !+ (100.*xc(nx)-80.)*dv  !- u*rho*Se*T(nx+1) +(8./3.)*T(nx+1)*gama*Sw/dx
aE(nx) = 0.0

do j=1,1000
	do i=1,nx
		T(i) = (aE(i)*T(i+1) + aW(i)*T(i-1) + aWW(i)*T(i-2) + Sp(i)) / aP(i)
	end do
end do


T(nx+1) = T(nx)

! !Solución analítica 
! do i=0,nx+1
!     Ta(i) = (T1-T0)*( (exp(rho*u*xc(i)/gama)-1)/(exp(rho*u*L/gama)-1) ) + T0
! enddo

!Escritura de resultados
 if (mod(it,100) == 0)then 
    call  WriteScalarField2D('Temp2D',it,T,xc,nx)
!      write(itchar, '(i6)')it !Convertir a un entero a un caracter
!      itchar=adjustl(itchar)
!      name='Temp2D'//itchar(1:len_trim(itchar))//'.txt'
!      write(2,*)'sp "'//name(1:len_trim(name))//'" w pm3d'
!      write(2,*)"pause 0.5"
 endif
! open(1,file="temp.txt", status="replace")
! do i=0,nx+1
! 	write(1,*)xc(i),T(i),Ta(i)
! end do

enddo
! close(1)
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
