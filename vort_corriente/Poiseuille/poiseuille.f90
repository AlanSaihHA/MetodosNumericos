Program Lidcavity
  implicit none
  real::dx,dy,dv,x0,xl,y0,yl,tolerance1,tolerance,residual,Div
  real::Se,Sw,Sn,Ss,ue,uw,un,us
  integer:: nx,ny,i,j,max_iter,itmax,it,flag
  real::  pi,time,dt,Re,eps
  real, allocatable::vor(:,:),psi(:,:),u(:,:),v(:,:),u1(:,:),v1(:,:),aE(:,:),aW(:,:),aN(:,:),aS(:,:),SP(:,:),aP(:,:)
  real, allocatable::xc(:),x(:),yc(:),y(:)  
  nx=150
  ny=30
  x0=0.0
  xl=10.0
  y0=0.0
  yl=1.0
  time=25.0
  dt=0.005
  Re=100.0
  eps=1.0/Re
  itmax=int(time/dt)+1
  max_iter=1000
  tolerance= 1.0e-5
  tolerance1=1.0e-4 !para la condición de paro en estado estacionario du/dt=0
  dx=(xl-x0)/float(nx)
  dy=(yl-y0)/float(ny)
  dv=dx*dy
  allocate(aE(nx,ny),aW(nx,ny),aN(nx,ny),aS(nx,ny),SP(nx,ny),aP(nx,ny))
  allocate(vor(0:nx+1,0:ny+1),psi(0:nx+1,0:ny+1),u(0:nx+1,0:ny+1),v(0:nx+1,0:ny+1),u1(0:nx+1,0:ny+1),v1(0:nx+1,0:ny+1))  
  allocate(xc(0:nx+1),x(0:nx),yc(0:ny+1),y(0:ny))    
  
  vor=0.0; psi=0.0; u=0.0;v=0.0; 
  aE=0.0; aW=0.0; aS=0.0; aN=0.0; SP=0.0; aP=0.0;
!  ue=0.0;uw=0.0;un=0.0;us=0.0;

  call Mesh1D(xc,x,x0,xl,nx)
  call Mesh1D(yc,y,y0,yl,ny)
  !definiendo las áreas  
  Se=dy; Sw=dy; Sn=dx; Ss=dx;  
  u = 1.0
  u(:,ny+1)=0.0;
  u(:,0)=0.0;
  u1=0.0;v1=0.0
!  do it=1,itmax !ciclo temporal
  flag=0
  it=1
  do while ((it.LE.itmax).AND.(flag.EQ.0))
     !!!!!!!!!!!!!!!!!!!!!!!!!!!! Vorticidad !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !cara norte
     do i=0,nx+1
         vor(i,ny+1)=-8.0*(psi(i,ny)-psi(i,ny+1))/(dy*dy) !-4.0*u(i,ny+1)/dy  
        !vor(i,ny+1)=-8.0*psi(i,ny)/(dy*dy)
     enddo
     !cara sur
     do i=0,nx+1
        vor(i,0)=-8.0*(-psi(i,0)+psi(i,1))/(dy*dy)  !-4.0*u(i,ny+1)/dy iria pero la vel es cero
        !vor(i,0)=-8.0*psi(i,1)/(dy*dy)
     enddo
     !cara este
     do j=0,ny+1
        vor(nx+1,j)=0.0                        !-8.0*(psi(nx,j)-psi(nx+1,j))/(dx*dx)   tenemos dvor/dy
     enddo
     !cara oeste
     !do j=0,ny+1
     !    vor(0,j)= -8.0*(-psi(0,j)+psi(1,j))/(dx*dx)                             !-8.0*(-psi(0,j)+psi(1,j))/(dx*dx)
        !vor(0,j)=-8.0*psi(1,j)/(dx*dx)
     !enddo

     do i=1,nx
        do j=1,ny
           ue=0.5*(u(i,j)+u(i+1,j))
           uw=0.5*(u(i,j)+u(i-1,j))
           un=0.5*(v(i,j)+v(i,j+1))
           us=0.5*(v(i,j)+v(i,j-1))
           aE(i,j)=eps*(Se/dx)-0.50*(ue*Se)
           aW(i,j)=eps*(Sw/dx)+0.50*(uw*Sw)
           aN(i,j)=eps*(Sn/dy)-0.50*(un*Sn)
           aS(i,j)=eps*(Ss/dy)+0.50*(us*Ss)      
           SP(i,j)=vor(i,j)*dv/dt !En este caso dv=dx*dy
           aP(i,j)=aE(i,j) + aW(i,j) +  aN(i,j) +  aS(i,j) + dv/dt
           u(i,j)=0.5*(ue+uw)
           v(i,j)=0.5*(un+us)
           Div=ue*Se-uw*Sw+un*Sn-us*Ss
        enddo
     enddo
     !corrección a las fronteras
     !Cara Este 
     !SP(nx,:)=SP(nx,:) - aE(nx,:)*vor(nx+1,1:ny)      !Condicion Neumann
     aP(nx,:)=aP(nx,:) - aE(nx,:) 
     aE(nx,:)=0.0 

     !Cara Oeste
     SP(1,:)=sp(1,:)+2.0*aW(1,:)*vor(0,1:ny)   !Condicion Dirichlet
     aP(1,:)=aP(1,:)+aW(1,:)
     aW(1,:)=0.0

     !Cara Norte  
     SP(:,ny)=SP(:,ny)+2.0*aN(:,ny)*vor(1:nx,ny+1)    !Condicion Dirichlet
     aP(:,ny)=aP(:,ny)+aN(:,ny)
     aN(:,ny)=0.0

     !Cara Sur
     SP(:,1)=SP(:,1)+2.0*aS(:,1)*vor(1:nx,0)        !Condicion Dirichlet
     aP(:,1)=aP(:,1)+aS(:,1)
     aS(:,1)=0.0

     call Gauss_TDMA2D(vor,nx,ny,aP,aE,aW,aN,aS,sP,nx,ny,max_iter,tolerance,residual) 
     call Stream(psi,nx,ny,x0,xl,y0,yl,yc,vor,tolerance,max_iter,residual)
     !esto es valido solo para los nodos internos 2,nx-1 y 2,ny-1.
     !falta correjir para las fronteras 
     do i=2,nx-1
        do j=2,ny-1
           u(i,j) =(psi(i,j+1)-psi(i,j-1))/(dy*2.0)           
           v(i,j) =-(psi(i+1,j)-psi(i-1,j))/(dx*2.0)
        enddo
     enddo
     !cálculo de las componentes de la velocidad para las fronteras
     !Sur
     do i=1,nx+1
        u(i,1)=(u(i,2)-u(i,0))/(1.5*dy)*(yc(1)-yc(0))+u(i,0)
        v(i,1)=(v(i,2)-v(i,0))/(1.5*dy)*(yc(1)-yc(0))+v(i,0)
     enddo
     !Oeste     
     do j=1,ny+1
        u(1,j)=(u(2,j)-u(0,j))/(1.5*dx)*(xc(1)-xc(0))+u(0,j)
        v(1,j)=(v(2,j)-v(0,j))/(1.5*dx)*(xc(1)-xc(0))+v(0,j)
     enddo
     !Norte     
     do i=1,nx+1
        u(i,ny)=(u(i,ny+1)-u(i,ny-1))/(1.5*dy)*(yc(ny)-yc(ny-1))+u(i,ny-1)
        v(i,ny)=(v(i,ny+1)-v(i,ny-1))/(1.5*dy)*(yc(ny)-yc(ny-1))+v(i,ny-1)
     enddo
     !Este     
     do j=1,ny+1
        u(nx,j)=(u(nx+1,j)-u(nx-1,j))/(1.5*dx)*(xc(nx)-xc(nx-1))+u(nx-1,j)
        v(nx,j)=(v(nx+1,j)-v(nx-1,j))/(1.5*dx)*(xc(nx)-xc(nx-1))+v(nx-1,j)
     enddo
          
     if (maxval(abs((u1-u)/dt)).LT.tolerance1) then
        flag=1
        itmax=it
        call  WriteScalarField2D('psi',it,psi,xc,yc,nx,ny)
        call  WriteVectorField('vel',it,u,v,xc,yc,nx,ny)
     else
        flag=0
        it=it+1
     endif
     write(*,*)it,Div,maxval(abs((u1-u)/dt))
     u(nx+1,:)=u(nx,:);  v(nx+1,:)=v(nx,:);       !corregir la frontera oeste
     u1=u;v1=v        
!     write(*,*)it,c,Div,maxval(abs((u1-u)/dt))
     if (mod(it,100) == 0)then
        call  WriteScalarField2D('psi',it,psi,xc,yc,nx,ny)
        call  WriteVectorField('vel',it,u,v,xc,yc,nx,ny)
     end if
  enddo
!enddo
end program Lidcavity
!!****************************************************************************************************************
subroutine  Stream(psi,nx,ny,x0,xl,y0,yl,yc,vor,tolerance,max_iter,residual)
  implicit none
  real dx,dy,Se,Sw,Sn,Ss,dv
  real x0,xl,y0,yl,tolerance,residual
  integer nx,i,j,ny,max_iter
  Real:: psi(0:nx+1,0:ny+1),vor(0:nx+1,0:ny+1),aP(nx,ny),aE(nx,ny),aW(nx,ny),aS(nx,ny),aN(nx,ny)
  Real:: SP(nx,ny),yc(0:ny+1)
  !determinado las condiciones de frontera, el manejo de los arreglos es parecido a Matlab
  dx=(xl-x0)/float(nx)
  dy=(yl-y0)/float(ny)
  dv=dx*dy
  !*************************condiciones de frontera**********************************
  !cara oeste
  psi(0,:)= yc(:)
  !cara este	
  !psi(nx+1,:)=0.0   !Condicion Neumann
  !cara sur
  psi(:,0)=0.0
  !cara norte
  psi(:,ny+1)=1.0
  !definiendo las áreas
  Se=dy; Sw=dy; Sn=dx; Ss=dx
  !**********Comienza a resolverse la ecuación de Poisson para vorticidad***********
  do i=1,nx
     do j=1,ny
        aE(i,j)=Se/dx
        aW(i,j)=Sw/dx
        aN(i,j)=Sn/dy
        aS(i,j)=Ss/dy
        SP(i,j)=vor(i,j)*dv !En este caso dv=dx*dy
        aP(i,j)=aE(i,j)+aW(i,j)+ aN(i,j)+aS(i,j)
     enddo
  enddo
  !!*********************************************************************************
  !
  !*******************Corrección de las condiciones de frontera**********************
  !
  !!*********************************************************************************
  !Cara Este
  !SP(nx,:)=SP(nx,:) - aE(nx,:)*psi(nx+1,1:ny)   !Condicion Neumann
  aP(nx,:)=aP(nx,:)- aE(nx,:)
  aE(nx,:)=0.0

  !Cara Oeste
  SP(1,:)=SP(1,:)+2.0*aW(1,:)*psi(0,1:ny)
  aP(1,:)=aP(1,:)+aW(1,:)
  aW(1,:)=0.0

  !Cara Norte
  SP(:,ny)=SP(:,ny)+ 2.0*aN(:,ny)*psi(1:nx,ny+1)
  aP(:,ny)=aP(:,ny)+aN(:,ny)
  aN(:,ny)=0.0

  !Cara Sur
  SP(:,1)=SP(:,1)+ 2.0*aS(:,1)*psi(1:nx,0)
  aP(:,1)=aP(:,1)+aS(:,1)
  aS(:,1)=0.0
  !!*********************************************************************************  
  !
  !*********************REsolviendo el sistema de ecuaciones*************************
  !
  !!********************************************************************************* 
  call Gauss_TDMA2D(psi,nx,ny,aP,aE,aW,aN,aS,sP,nx,ny,max_iter,tolerance,residual)  
end subroutine Stream
!!****************************************************************************************************************
Subroutine Mesh1D(xc,x,x0,xl,nx)
  implicit none
  integer i,j,nx
  real x0,xl,dx
  real x(0:nx),xc(0:nx+1)
  dx=(1.0)/float(nx)
  do i=0,nx
     x(i)=float(i)*dx
     x(i)=x0+(xl-x0)*x(i)
  end do
  xc(0)=x(0); xc(nx+1)=x(nx);
  do i=1,nx
     xc(i)=(x(i)+x(i-1))*0.5
  end do
End Subroutine Mesh1D
!!*******************************************************************************************************************

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

end subroutine TDMA


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

end subroutine lineX_2D


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

end subroutine lineY_2D


!****************************
!****************************

real function calcResidual(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny)
  implicit none
  integer bi,ei,bj,ej,i,j,nx,ny
  real phi(0:ei+1,0:ej+1),aP(nx,ny),aE(nx,ny),aW(nx,ny),aN(nx,ny),aS(nx,ny),sP(nx,ny)
  real acum(ei,ej),residual,NINV
  bi=1; bj=1
  acum=0.
  NINV = 1.0 / float(ei*ej)
  do i=bi,ei
     do j=bj,ej
        acum(i,j) = aE(i,j) * phi(i+1,j) + aW(i,j) * phi(i-1,j) + aN(i,j) * phi(i,j+1) + &
             aS(i,j) * phi(i,j-1) + sp(i,j) - aP(i,j) * phi(i,j)
     end do
  end do
  residual = sqrt( NINV * sum(acum * acum) )
  calcResidual=residual
end function calcResidual



subroutine Gauss_TDMA2D(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny,max_iter,tolerance,residual)

  implicit none

  integer bi,ei,bj,ej,i,j,nx,ny,count_iter,max_iter
  real phi(0:ei+1,0:ej+1),aP(nx,ny),aE(nx,ny),aW(nx,ny),aN(nx,ny),aS(nx,ny),sP(nx,ny)
  real residual,tolerance

  interface
     real function calcResidual(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny)
       implicit none
       integer bi,ei,bj,ej,i,j,nx,ny
       real phi(0:ei+1,0:ej+1),aP(nx,ny),aE(nx,ny),aW(nx,ny),aN(nx,ny),aS(nx,ny),sP(nx,ny)
     end function calcResidual
  end interface

  count_iter=0;  residual=1.0

  do while((count_iter <= max_iter).and.(residual > tolerance)) 
     call lineX_2D(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny)
     call lineY_2D(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny)
     residual = calcResidual(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny)
     count_iter=count_iter+1
  end do
  !write(*,*)count_iter
end subroutine Gauss_TDMA2D


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!********************************************************************

Subroutine WriteScalarField2D(Name,kx,T,xc,yc,nx,ny)
  implicit none  
  integer i,j,nx,ny,kx
  real T(0:nx+1,0:ny+1),xc(0:nx+1),yc(0:ny+1)
  character*(*)Name
  character*50 txt,Filename
  write(txt,'(i6)')kx
  txt=ADJUSTL(txt)
  Filename=name//txt(1:len_trim(txt))//".txt"

  open(10,file=Filename(1:len_trim(Filename)))
  do j=0,ny+1
     do i=0,nx+1
	write(10,*)xc(i),yc(j),T(i,j)
     end do
     write(10,*)''
  end do
  close(10)
End Subroutine WriteScalarField2D
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

Subroutine WriteVectorField(Name,kx,uc,vc,xc,yc,nx,ny)
  implicit none
  integer i,j,nx,ny,kx
  real uc(0:nx+1,0:ny+1),vc(0:nx+1,0:ny+1),xc(0:nx+1),yc(0:ny+1)
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
End Subroutine WriteVectorField

