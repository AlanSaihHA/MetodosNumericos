Program Lidcavity
  implicit none
  real::dx,dr,dv,x0,xl,r0,rl,tolerance,tolerance1,residual,Div
  real::Se,Sw,Sn,Ss,ue,uw,un,us,rp,rn,rs
  integer:: nx,nr,i,j,max_iter,itmax,it,c,flag
  real::  pi,time,dt,Re,eps
  real, allocatable::vor(:,:),psi(:,:),u(:,:),v(:,:),u1(:,:),v1(:,:),aE(:,:),aW(:,:),aN(:,:),aS(:,:),SP(:,:),aP(:,:),res(:,:)
  real, allocatable::xc(:),x(:),rc(:),r(:),W(:,:)  
  nx=200
  nr=60
  x0=0.0
  xl=10.0
  r0=0.01
  rl=1.0
pi=acos(-1.0)
  time=10.0
  dt=0.005
  Re=100.0
  eps=1.0/Re
  itmax=int(time/dt)+1
  max_iter=1000
  tolerance= 1.0e-5
  tolerance1=1.0e-4 !para la condición de paro en estado estacionario du/dt=0
  dx=(xl-x0)/float(nx)
  dr=(rl-r0)/float(nr)
  dv=dx*dr
  allocate(aE(nx,nr),aW(nx,nr),aN(nx,nr),aS(nx,nr),SP(nx,nr),aP(nx,nr))
  allocate(vor(0:nx+1,0:nr+1),psi(0:nx+1,0:nr+1),u(0:nx+1,0:nr+1),&
  v(0:nx+1,0:nr+1),W(0:nx+1,0:nr+1),u1(0:nx+1,0:nr+1),v1(0:nx+1,0:nr+1))  
  allocate(xc(0:nx+1),x(0:nx),rc(0:nr+1),r(0:nr),res(0:nx+1,0:nr+1))    
  
  vor=0.0; psi=0.0; u=1.0;v=0.0; res=0.0;W=0.0;
  aE=0.0; aW=0.0; aS=0.0; aN=0.0; SP=0.0; aP=0.0;
!  ue=0.0;uw=0.0;un=0.0;us=0.0;

  call Mesh1D(xc,x,x0,xl,nx)
  call Mesh1D(rc,r,r0,rl,nr)
  !definiendo las áreas  
  Se=dr; Sw=dr; Sn=dx; Ss=dx; dv=dx*dr 
  u(:,nr+1)=0.0;u(:,0)=u(:,1);
  u1=u;v1=v
  flag=0
  it=1
  do while((it.LE.itmax).AND.(flag.EQ.0))!ciclo temporal     
!!!!****************TRANSPORTE DE VORTICIDAD****************************************
     aP=0.0;aE=0.0;aW=0.0;aN=0.0; aS=0.0; SP=0.0; !vor=0.0
     !cara norte
     do i=0,nx+1
         vor(i,nr+1)=-2.0*(psi(i,nr)-psi(i,nr+1))/(rc(nr+1)*dr*dr) 
        !vor(i,nr+1)=-8.0*psi(i,nr)/(dy*dy)
     enddo
     !cara sur
     do i=0,nx+1
        vor(i,0)=0.0!-2.0*(-psi(i,0)+psi(i,1))/(rc(0)*dr*dr)
        !vor(i,0)=-8.0*psi(i,1)/(dy*dy)
     enddo
     !cara este
     !do j=0,nr+1
      !  vor(nx+1,j)=-8.0*psi(nx,j)/(dx*dx)
     !enddo
     !cara oeste
     do j=0,nr+1
         vor(0,j)=-2.0*(psi(1,j)-psi(0,j))/(rc(j)*dx*dx) + 1.0/rc(j)
        !vor(0,j)=-8.0*psi(1,j)/(dx*dx)
     enddo
     do i=1,nx
        do j=1,nr
           ue=0.5*(u(i,j)+u(i+1,j))
           uw=0.5*(u(i,j)+u(i-1,j))
           un=0.5*(v(i,j)+v(i,j+1))
           us=0.5*(v(i,j)+v(i,j-1))
           aE(i,j)=eps*(Se/dx)-0.50*(ue*Se)
           aW(i,j)=eps*(Sw/dx)+0.50*(uw*Sw)
           aN(i,j)=eps*(Sn/dr)-0.50*(un*Sn)
           aS(i,j)=eps*(Ss/dr)+0.50*(us*Ss)      
           SP(i,j)=vor(i,j)*dv/dt + (2.0*dv*W(i,j)/rc(j))*(W(i+1,j)-W(i-1,j))/(2.0*dx) - (eps*vor(i,j)*dv)/(rc(j)*rc(j)) 
           !SP(i,j)=vor(i,j)*dv/dt - (eps*vor(i,j)*dv)/(rc(j)*rc(j)) + vor(i,j)*v(i,j)*dv/rc(j) - (dv*W(i,j)/rc(j))*(W(i+1,j)-W(i-1,j))/(2.0*dx)
           aP(i,j)=aE(i,j) + aW(i,j) +  aN(i,j) +  aS(i,j) + dv/dt
           !Div=ue*Se-uw*Sw+un*Sn-us*Ss
           u(i,j)=0.5*(ue+uw)
           v(i,j)=0.5*(un+us)		
        enddo
     enddo
     !corrección a las fronteras
     !Cara Este 
     !SP(nx,:)=SP(nx,:)+2.0*aE(nx,:)*vor(nx+1,1:nr)  
     aP(nx,:)=aP(nx,:)-aE(nx,:) 
     aE(nx,:)=0.0 

     !Cara Oeste
     SP(1,:)=sp(1,:)+2.0*aW(1,:)*vor(0,1:nr)   
     aP(1,:)=aP(1,:)+aW(1,:)
     aW(1,:)=0.0

     !Cara Norte  
     SP(:,nr)=SP(:,nr)+2.0*aN(:,nr)*vor(1:nx,nr+1)
     aP(:,nr)=aP(:,nr)+aN(:,nr)
     aN(:,nr)=0.0

     !Cara Sur
     SP(:,1)=SP(:,1)+2.0*aS(:,1)*vor(1:nx,0)
     aP(:,1)=aP(:,1)+aS(:,1)
     aS(:,1)=0.0

     call Gauss_TDMA2D(vor,nx,nr,aP,aE,aW,aN,aS,sP,nx,nr,max_iter,tolerance,residual) 
     vor(nx+1,:)=vor(nx,:) 
     !******************************CÁLCULO DE LA VELOCIDAD SWIRL ************************
     aP=0.0;aE=0.0;aW=0.0;aN=0.0; aS=0.0; SP=0.0;!W=0.0
     do i=1,nx
        do j=1,nr
           ue=0.5*(u(i,j)+u(i+1,j))
           uw=0.5*(u(i,j)+u(i-1,j))
           un=0.5*(v(i,j)+v(i,j+1))
           us=0.5*(v(i,j)+v(i,j-1))
           aE(i,j)=eps*(Se/dx)-0.50*(ue*Se)
           aW(i,j)=eps*(Sw/dx)+0.50*(uw*Sw)
           aN(i,j)=eps*(Sn/dr)-0.50*(un*Sn)
           aS(i,j)=eps*(Ss/dr)+0.50*(us*Ss)      
           SP(i,j)=W(i,j)*dv/dt - 2.0*v(i,j)*W(i,j)*dv/rc(j) - eps*W(i,j)*dv/(rc(j)*rc(j))  
           aP(i,j)=aE(i,j) + aW(i,j) +  aN(i,j) +  aS(i,j) + dv/dt 
        enddo
     enddo
     !CORRECCIÓN A LAS FRONTERAS
     !Cara Este 
     !SP(nx,:)=SP(nx,:)+2.0*aE(nx,:)*W(nx+1,1:nr)  
     aP(nx,:)=aP(nx,:)-aE(nx,:) 
     aE(nx,:)=0.0 

     !Cara Oeste
     SP(1,:)=sp(1,:)+2.0*aW(1,:)*W(0,1:nr)   
     aP(1,:)=aP(1,:)+aW(1,:)
     aW(1,:)=0.0

     !Cara Norte  
     SP(:,nr)=SP(:,nr)+2.0*aN(:,nr)*W(1:nx,nr+1)
     aP(:,nr)=aP(:,nr)+aN(:,nr)
     aN(:,nr)=0.0

     !Cara Sur
     SP(:,1)=SP(:,1)+2.0*aS(:,1)*W(1:nx,0)
     aP(:,1)=aP(:,1)+aS(:,1)
     aS(:,1)=0.0
     call Gauss_TDMA2D(W,nx,nr,aP,aE,aW,aN,aS,sP,nx,nr,max_iter,tolerance,residual)      
	W(nx+1,:)=W(nx,:)
!!!!*********************************************************************************

!!****************************************** FUNCIÓN DE CORRIENTE **************************************************
     call Stream(psi,nx,nr,x0,xl,r0,rl,rc,r,vor,tolerance,max_iter,residual)           
  !!!!!!!!!!!!!!!!!!!!CÁLCULO DEL CAMPO DE VELOCIDADES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11111
     !esto es valido solo para los nodos internos 2,nx-1 y 2,nr-1.
     !falta corregir para las fronteras 
     do i=2,nx-1
        do j=2,nr-1
           u(i,j) =(1.0/rc(j))*(psi(i,j+1)-psi(i,j-1))/(dr*2.0)           
           v(i,j) =-(1.0/rc(j))*(psi(i+1,j)-psi(i-1,j))/(dx*2.0)
        enddo
     enddo
     !cálculo de las componenetes de la velocidad para las fronteras
     !Sur
     do i=1,nx+1
        u(i,1)=(u(i,2)-u(i,0))/(1.5*dr)*(rc(1)-rc(0))+u(i,0)
        v(i,1)=(v(i,2)-v(i,0))/(1.5*dr)*(rc(1)-rc(0))+v(i,0)
     enddo
     !Oeste     
     do j=1,nr+1
        u(1,j)=(u(2,j)-u(0,j))/(1.5*dx)*(xc(1)-xc(0))+u(0,j)
        v(1,j)=(v(2,j)-v(0,j))/(1.5*dx)*(xc(1)-xc(0))+v(0,j)
     enddo
     !Norte     
     do i=1,nx+1
        u(i,nr)=(u(i,nr+1)-u(i,nr-1))/(1.5*dr)*(rc(nr)-rc(nr-1))+u(i,nr-1)
        v(i,nr)=(v(i,nr+1)-v(i,nr-1))/(1.5*dr)*(rc(nr)-rc(nr-1))+v(i,nr-1)
     enddo
     !Este     
     do j=1,nr+1
        u(nx,j)=(u(nx+1,j)-u(nx-1,j))/(1.5*dx)*(xc(nx)-xc(nx-1))+u(nx-1,j)
        v(nx,j)=(v(nx+1,j)-v(nx-1,j))/(1.5*dx)*(xc(nx)-xc(nx-1))+v(nx-1,j)
     enddo
     u(nx+1,:)=u(nx,:);  v(nx+1,:)=v(nx,:);u(:,0)=u(:,1)
     if (maxval(abs((u1-u)/dt)).LT.tolerance1) then        
        flag=1        
        itmax=it
        call  WriteScalarField2D('psi',it,psi,xc,rc,nx,nr)
        call  WriteVectorField('vel',it,u,v,xc,rc,nx,nr)
     else
        flag=0
        it=it+1
     endif
     write(*,*)it,maxval(abs((u1-u)/dt))
     !u(nx+1,:)=u(nx,:);  v(nx+1,:)=v(nx,:);
     u1=u;v1=v        
     if (mod(it,100) == 0)then
        call  WriteScalarField2D('psi',it,psi,xc,rc,nx,nr)
        call  WriteVectorField('vel',it,u,v,xc,rc,nx,nr)
     end if
  enddo
end program Lidcavity
!!****************************************************************************************************************
subroutine  Stream(psi,nx,nr,x0,xl,r0,rl,rc,r,vor,tolerance,max_iter,residual)
  implicit none
  real dx,dr,Se,Sw,Sn,Ss,dv,rp,rn,rs
  real x0,xl,r0,rl,tolerance,residual
  integer nx,i,j,nr,max_iter
  Real:: psi(0:nx+1,0:nr+1),vor(0:nx+1,0:nr+1),aP(nx,nr),aE(nx,nr),aW(nx,nr),aS(nx,nr),aN(nx,nr)
  Real:: SP(nx,nr),rc(0:nr+1),r(0:nr)
  !determinado las condiciones de frontera, el manejo de los arreglos es parecido a Matlab
  aP=0.0;aE=0.0;aW=0.0;aN=0.0; aS=0.0; SP=0.0; psi=0.0
  dx=(xl-x0)/float(nx)
  dr=(rl-r0)/float(nr)
  dv=dx*dr
  !*************************condiciones de frontera**********************************
  !cara oeste
  psi(0,:)=0.5*rc(:)*rc(:)    !condición dirichlet 
  !cara este	
  !psi(nx+1,:)=0.0 ! Condición neumman
  !cara sur
  psi(1:nx,0)=0.5*rc(0)*rc(0) !condición dirichlet
  !cara norte
  psi(:,nr+1)=0.5*rc(nr+1)*rc(nr+1)  !condición dirichlet
  !definiendo las áreas
    Se=dr; Sw=dr; Sn=dx; Ss=dx;  
  !**********Comienza a resolverse la ecuación de Poisson para vorticidad***********
  do i=1,nx
     do j=1,nr
	rp=rc(j);rn=r(j);rs=r(j-1)   
        aE(i,j)=Se/(rp*dx)
        aW(i,j)=Sw/(rp*dx)
        aN(i,j)=Sn/(rn*dr)
        aS(i,j)=Ss/(rs*dr)
        SP(i,j)=vor(i,j)*dv - psi(i,j)*dv/rc(j)**3
        aP(i,j)=aE(i,j)+aW(i,j)+ aN(i,j)+aS(i,j)
     enddo
  enddo
  
  !*******************Corrección de las condiciones de frontera**********************
  !Cara Este
  !SP(nx,:)=SP(nx,:)+2.0*aE(nx,:)*psi(nx+1,1:ny)
  aP(nx,:)=aP(nx,:)-aE(nx,:)
  aE(nx,:)=0.0

  !Cara Oeste
  SP(1,:)=SP(1,:)+2.0*aW(1,:)*psi(0,1:nr)
  aP(1,:)=aP(1,:)+aW(1,:)
  aW(1,:)=0.0

  !Cara Norte
  SP(:,nr)=SP(:,nr)+ 2.0*aN(:,nr)*psi(1:nx,nr+1)
  aP(:,nr)=aP(:,nr)+aN(:,nr)
  aN(:,nr)=0.0

  !Cara Sur
  SP(:,1)=SP(:,1)+ 2.0*aS(:,1)*psi(1:nx,0)
  aP(:,1)=aP(:,1)+aS(:,1)
  aS(:,1)=0.0
  !*********************REsolviendo el sistema de ecuaciones*************************
  call Gauss_TDMA2D(psi,nx,nr,aP,aE,aW,aN,aS,sP,nx,nr,max_iter,tolerance,residual)  
  ! Condición de frontera Newmann dPhi/dx=phiE-phiP = 0 => PhiE=PhiP
  psi(nx+1,:)=psi(nx,:)
end subroutine Stream
!!*****************************************************************************************
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

