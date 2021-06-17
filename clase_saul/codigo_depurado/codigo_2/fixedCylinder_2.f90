
!  A TWO-DIMENSIONAL CODE FOR FLOW OVER A FIXED CYLINDER

                            
      program main
	use grid; use bounds
	use front; use press
	use flprop; use fparam
	use source
	implicit none

	integer i,j,k,kk,ii,jj,iter,nprh
	integer nfronts,nstep,lp,ll
	integer iem,nps
	integer maxstep,ibd,ibdrx,jbdrx,npvel,nppr,nback
	integer lfnm,iord
	integer	p1,p2
	integer is,npsx,npsy,nptrack,nxp2,nyp2,maxpt,maxel

	real*8 ptime,time,endtime,r0,a1,a2,a3,a4,a5,a6,a7
	real*8 am,close,am1
	real*8 bm,sm,dtvisc,pi,pd2,cnst,up,vp,drx,dry,xt,yt,Fhx,Fhy
	real*8 vmax,dt,dtadv,xc,yc,vv,velm
	real*8 Di,Re1,Re2,ucm,ucc,vcc,xcc,ycc,centx1,centy1,cent11,cent22,xxt,yyt,prop1,prop2
	real*8 ucco,vcco,Div,gama,Re
!---------------grid code-----------------------------------------------!
      	real*8, allocatable::u(:,:), v(:,:)
	real*8, allocatable::u1(:,:), v1(:,:),uc(:,:), vc(:,:)
      	real*8, allocatable::ut(:,:),vt(:,:)
      	real*8, allocatable::r(:,:), ro(:,:),p(:,:)
      	real*8, allocatable::tmp1(:,:), tmp2(:,:)
      	real*8, allocatable::tmp3(:,:), tmp4(:,:)
      	real*8, allocatable::color(:,:)
	real*8, allocatable::fx(:,:), fy(:,:)			!,u2(:,:), v2(:,:)
!---------------grid/front code-----------------------------------------!
      	real*8, allocatable::temp1(:,:), temp2(:,:)
      	real*8, allocatable::rtmp(:,:)
!---------------front tracking code-------------------------------------!
     
      	real*8, allocatable::pt(:,:),t1(:,:),fp(:,:),t2(:,:), pto(:,:),norm(:,:)
      	real*8, allocatable::cv(:,:),cnt(:,:), elprop(:,:)
      	integer, allocatable::icp(:,:), ine(:,:),iptmp(:),ip(:),ptcon(:),elcon(:),bptcon(:),belcon(:),istore(:,:)

!----------------------------------------------------------------------------
           character*20 fname, xname, statfile, datafile, filmfile,command

	  
	    logical restart, movie
 
      integer:: noutf=1,noutv=1,noutp=1,nbkg=1,nbkf=1,noutr=1,nouto=1,nouts=1,noutu=1,noutpt=1,noutphi=1
      nprh=1
	  
	nxp2=514
        nyp2=258
	maxel=1000
	maxpt=1000

	allocate(u(nxp2,nyp2), v(nxp2,nyp2))
	allocate(u1(nxp2,nyp2), v1(nxp2,nyp2),uc(nxp2,nyp2), vc(nxp2,nyp2))
      	allocate(ut(nxp2,nyp2),vt(nxp2,nyp2))

      	allocate(r(nxp2,nyp2), ro(nxp2,nyp2),p(nxp2,nyp2))
      	allocate(tmp1(nxp2,nyp2), tmp2(nxp2,nyp2))
      	allocate(tmp3(nxp2,nyp2), tmp4(nxp2,nyp2))
      	allocate(color(nxp2,nyp2))
	allocate(fx(nxp2,nyp2), fy(nxp2,nyp2))			!,u2(nxp2,nyp2), v2(nxp2,nyp2)
!---------------grid/front code-----------------------------------------!
      	allocate(temp1(nxp2+2,nyp2+2), temp2(nxp2+2,nyp2+2))
      	allocate(rtmp(nxp2,nyp2))
!---------------front tracking code-------------------------------------!
     
      	allocate(pt(maxpt,2),t1(maxel,2),fp(20,6),t2(maxel,2), pto(maxpt,2),norm(maxel,2))
      	allocate(cv(maxel,2),cnt(maxel,2), elprop(maxel,3),istore(maxel,2))
      	allocate(icp(maxel,2), ine(maxel,2),iptmp(maxpt),ip(20),ptcon(maxpt),elcon(maxel),bptcon(maxpt),belcon(maxel))

	  Open(7003,File='Force.txt',Status='replace')
	  Open(1020,File='normals.txt',Status='replace') 
!  THESE PARAMETERS SHOULD BE INPUTED 
 
!=======================================================================================
!				Initial data
!======================================================================================	  
	  maxstep=1000000
	  endtime=1000.0d0            !Se utiliza 0.0d0 para doble precision, ejemplo 5.33d6=5360000.00000000
	  xl=16.0d0
	  yl=8.0d0
	  r0=1.0d0
	  r1=2.6d0       !densidad del solido
	  r2=1.0d0       !densidad del liquido o fluido
	  Ifi1=0.0d0      !Campo indicador en el solido
	  Ifi2=1.0d0      !Campo indicador en el fluido
	  gx=0.0d0
	  gy=0.0d0
	  dt=0.0005d0
	  Re=100.0d0
	  gama=1.0d0/Re
	  xcc=4.0d0      !centro posicion del objeto solido
	  ycc=4.0d0       !centro posicion del objeto solido
	  

!=======================================================================================
!				Pressure Solution data
!======================================================================================	  
	  itmax=1000
	  xerr=0.000001d0
	  beta=1.2
!=======================================================================================
!				Boundary conditions
!========================================================================================
	  ibdry=0    !0 ciclicas
	  jbdry=0    ! 1 de frontera	  

	ucco=0.0d0   !velocidad inicial del cuerpo
	vcco=0.0d0   !velocidad inicial del cuerpo

	npvel=1
	nppr=1
	nback=6000000	  
          write(*,*)' TYPE name of output files'    !Escritura de directorio donde se guardaran los datos
	     read (5,'(a10)')fname    !nombre de carpeta donde se guarda
         lfnm=len_trim(fname)+1     !long. de cadena de nombre
         
       command='mkdir '
       command(7:7+lfnm)=fname  !Se escribe en espacio 7 nombre de carpeta 
       call system(command)     !se ejecuta en consola
       xname=fname
       fname(lfnm:lfnm+1)='/'   !resulta en "NombreCarpeta"/
       lfnm=len_trim(fname)+1 !ltr(fname,20)+1
                
      fxl=xl    !Se iguala longitud en x
      fyl=yl    ! Se iguala longiutd en y
      xmv=30.0d0
      mv=30

      call gsetup(ut,vt,ro,p,tmp1,tmp2,tmp3,nxp2,nyp2)        !inicia las variables (las coloca en cero)

! SET PARAMETERS FOR FRONT RESTRUCTURING
      amax=0.8/hxi              !maximo tamanyo de un elemento   hxi=1/dx
      amin=0.2/hxi              !minimo tamanyo de un elemento   hxi=1/dx
      am=0.5*(amax+amin)        !distancia "optima"

      time=0.0d0
      
      call finit(nfronts,fp,ip,pt,icp,ine,ptcon,elcon,bptcon,belcon,elprop,maxpt,maxel,am,prop1,prop2,xcc,ycc)    !Crear el círculo
      call fcurv(cv,t1,t2,cnt,pt,icp,ine,ptcon,elcon,maxel,maxpt,norm)               !Calcula vectores normales al contorno

	k=ffe
	do ll=1,ne	!escribe los vectores normales al contorno
		write(1020,*)cnt(k,1),cnt(k,2),norm(k,1),norm(k,2)
		k=elcon(k)
	enddo					  			  
!========================================================================================================
        call setdens(r0,nfronts,xcc,ycc,r,nxp2,nyp2,prop1)
        write(*,*)maxval(r),minval(r)

	call dens2(r,rtmp,temp1,temp2,nxp2,nyp2,pt,ptcon,icp,elcon,elprop,maxpt,maxel)   !calculo de funcion marcador con gradiente
	call gifield(color,r,nxp2,nyp2)
	
!Initial conditions	
	u=1.0d0 !todas vel axiales en 0
	v=0.0d0
	do i=1,nxp2
	do j=1,nyp2
	u(i,j)=u(i,j)*0.5d0*(color(i,j)+color(i+1,j))   !vel en solido =0 y en fluido =1
	v(i,j)=v(i,j)*0.5d0*(color(i,j)+color(i,j+1))
	enddo
	enddo	
	
        u1=u
	v1=v


!    { START TIME INTEGRATION} 
 
      do 100 nstep=0,maxstep


!================================================================
	
	time=float(nstep)*dt
	call Forceh(pt,icp,elcon,maxel,maxpt,norm,u,v,p,Fhx,Fhy,cnt,gama,nxp2,nyp2,r,xcc,ycc)
	write(7003,*)time,Fhx,Fhy
! PRINTOUT IF WANTED: 

      if (mod(nstep,1000) .eq. 0.0d0) then
      	write(*,3000)time
3000      format(' Printout at Time: ',e12.5)

	call printout(pt,ptcon,iptmp,istore,icp,elcon,maxpt,maxel,color,p,u,v,time,fname,noutp,lfnm,nxp2,nyp2)		  
      end if
!============================================================

!--------------printout done--------------------

      If(time .gt. endtime)go to 101
! solving N-S equations
      call simplec(p,u(1:nx1,1:ny1+1),v(1:nx1+1,1:ny1),xerr,dt,r,ro,fx,fy,u1(1:nx1,1:ny1+1),v1(1:nx1+1,1:ny1),  &
      			Div,xcc,ycc,color,time,gama)
      u=u1
      v=v1
      write(*,*)nstep,time

	
      time=time+dt

100   continue 
!                             { END Time Integration}
101   continue
	  
      close(1020,status='keep')	  
      close(7003,status='keep')

end program   
!-----------------------------------------------------------------------!
!-----------------------------------------------------------------------! 
! subroutine to set grid density from initial front.
      subroutine setdens(r0,nfronts,xcc,ycc,r,nxp2,nyp2,prop1) 

	  use grid; use bounds
	  use fparam
	  implicit none 
	  
	  integer nfronts,nxp2,nyp2,i,j,l
	  real*8  x,y,xx,yy,xfl,pi,r0,xcc,ycc,prop1
	  real*8:: r(nxp2,nyp2)

      pi=4.*ATAN(1.)

      do j=1,nyp2    !r=densidad, r0=1.0, prop1=r1-r2=1.6
      do i=1,nxp2
  	    r(i,j)=r0
      enddo
      enddo
      do j=1,nyp2
      do i=1,nxp2 
       x=float(i-1)/hxi   !obtengo distancia del inicio al punto
       y=float(j-1)/hyi
             xx=abs(x-xcc)    !diferencia de distancia del punto al centro del cilindro
             yy=abs(y-ycc)

            xfl=sqrt(xx**2+yy**2)-0.5d0    !reviso si estoy dentro del cilindro

         !pt(i+nptot,1)=xcc + rad*cos(th)   !coordenada en x
         !pt(i+nptot,2)=ycc + rad*sin(th)   !coordenada en y
                     
                     
            	!xfl = sqrt((x-pt(i-1,1))**2+(y-pt(m,2))**2)
            	!xfl = xfl+xfl
            
      if(xfl.le. 0d0) r(i,j)=r(i,j)+prop1   !si estoy dentro la densidad es 2.6
      
      
      enddo
      enddo

      return
      end

!-----------------------------------------------------------------------!
      subroutine gsetup(ut,vt,ro,p,tmp1,tmp2,tmp3,nxp2,nyp2)
	  
	  use grid
	  use bounds

	  implicit none

      integer i,j,nxp2,nyp2

      real*8:: ut(nxp2,nyp2),vt(nxp2,nyp2)
      real*8:: tmp1(nxp2,nyp2), tmp2(nxp2,nyp2)
      real*8:: tmp3(nxp2,nyp2),ro(nxp2,nyp2)
      real*8:: p(nxp2,nyp2)
   
      nxp1=nxp2-1    !todas las n menos 1
      nyp1=nyp2-1 
      nx=nxp2-2      !todas las n sin dos, sin celdas de fronteras
      ny=nyp2-2
      nxlast=nxp1   
      nylast=nyp1
	  
      nx1=nxp1
      ny1=nyp1

      hxi=float(nxp1-1)/xl    !xl=largo=16
      hyi=float(nyp1-1)/yl    !hxi es nx/yl, es decir, el inverso de dy 1/dy
                    
      do i=1,nxp2    !todas las varaiables son puestas en cero
      do j=1,nyp2
       ut(i,j)=0.0d0
       vt(i,j)=0.0d0
       ro(i,j)=0.0d0 
       p(i,j)=0.0d0
       tmp1(i,j)=0.0d0
       tmp2(i,j)=0.0d0
       tmp3(i,j)=0.0d0
      enddo
      enddo

      return
      end
!-----------------------------------------------------------------------!        
subroutine gifield(color,r,nxp2,nyp2)
	use flprop

	implicit none
	integer nxp2,nyp2,i,j
	real*8:: r(nxp2,nyp2),color(nxp2,nyp2) 
 
!Ifi1=0.0d0      !Campo indicador en el solido
!Ifi2=1.0d0      !Campo indicador en el fluido
!r1=2.6d0       !densidad del solido
!r2=1.0d0       !densidad del liquido o fluido 
	do i=1,nxp2
	do j=1,nyp2
	color(i,j)= (Ifi1*(r(i,j)-r2)-Ifi2*(r(i,j)-r1) )/(r1-r2)  !0 en soldio, 1 en fluido

	if (color(i,j) .gt. Ifi2) then
		color(i,j)=Ifi2
	end if 

	if (color(i,j) .lt. Ifi1) then
		color(i,j)=Ifi1
	end if 
	
	enddo
	enddo

end
!-----------------------------------------------------------------------!
 subroutine dens2(r,rtmp,tmp1,tmp2,nxp2,nyp2,pt,ptcon,icp,elcon,elprop,maxpt,maxel)  !campo indicador

      use grid !/hxi,hyi,nxp1,nyp1,nx,ny,xl,yl
      use bounds !/ibdry,jbdry,nxlast,nylast
      use front !/np,ffp,lfp,fep, ne,ffe,lfe,fee,amax,amin 
      use fparam !/xmv,mv,fxl,fyl
      use flprop !/r1,r2,m1,m2,gx,gy,rro  

	  implicit none

	  

      integer nxp2,nyp2,i,j,ii,jj,ir,jr,i1,j1,maxpt,maxel,k,ll

      real*8 tmp1(nxp2+2,nyp2+2), tmp2(nxp2+2,nyp2+2)
      real*8 rtmp(nxp2,nyp2)
      real*8 r(nxp2,nyp2)
      integer iii(4000),jjj(4000)
      real*8 pt(maxpt,2) 
      integer ptcon(maxpt)
      integer icp(maxel,2),elcon(maxel)
      real*8 elprop(maxel,3)

      real*8 xx,yy,p0x,p0y,p1x,p1y,pi,pd2,cnstx,cnsty,sc,xc,yc,x1,y1,rx,ry,xt,yt
      real*8 drx,dry,err,cf,errmax,cc,rold

      integer itmax,numnd,node,bet,iter,jo,io

      xx=(xmv+0.5d0)*xl    
      yy=(xmv+0.5d0)*yl

! generate gradient vector
 
      sc=0.5d0
      pi=4.0d0*ATAN(1.0d0)
      pd2=0.5d0*pi
      cnstx=(.25d0*hxi)**2 
      cnsty=(.25d0*hyi)**2

      do i=1,nxp2+2
      do j=1,nyp2+2  
       tmp1(i,j)=0.0d0
       tmp2(i,j)=0.0d0
      enddo
      enddo

      k=ffe
      do ll=1,ne 
       p0x=pt(icp(k,1),1)
       p0y=pt(icp(k,1),2)
       p1x= p0x+dmod((pt(icp(k,2),1)-p0x+xx),xl)-0.5d0*xl
       p1y= p0y+dmod((pt(icp(k,2),2)-p0y+yy),yl)-0.5d0*yl

       xc=0.5*(p1x+p0x)
       yc=0.5*(p1y+p0y)
       x1=p1x-p0x
       y1=p1y-p0y

! calculate components of n.dA
      rx=-y1 *elprop(k,1)
      ry= x1 *elprop(k,1)
                      
! generate a grid value   
       xt=dmod((xc+0.5d0/hxi)+xmv*xl,xl)
       yt=dmod((yc+0.5d0/hyi)+xmv*yl,yl)

       ir=1+int(nx*xt/xl)
       jr=1+int(ny*yt/yl)
                    
       do i1=1,4
       do j1=1,4
        ii=ir-2+i1   
        jj=jr-2+j1

         drx = 1. + cos((xt*hxi - float(ii-1))*pd2)
         dry = 1. + cos((yt*hyi - float(jj-1))*pd2)

        tmp1(ii+1,jj+1)=tmp1(ii+1,jj+1)+rx*drx*dry*cnstx
        tmp2(ii+1,jj+1)=tmp2(ii+1,jj+1)+ry*drx*dry*cnsty
       enddo
       enddo

      k=elcon(k)
      enddo   

! correct boundaries:
      call gbd1(tmp1,tmp2,nxp2,nyp2)

       do i=1,nxp2
       do j=1,nyp2
        rtmp(i,j)=0.0d0
       enddo
       enddo

       do i=2,nxp1
       do j=2,nyp1
        rtmp(i,j)= 0.5d0*(hxi*(tmp1(i+1+1,j+1)-tmp1(i-1+1,j+1))+		&
             hyi*(tmp2(i+1,j+1+1)-tmp2(i+1,j-1+1)))
       enddo
       enddo

! set connection
      node=0
       do 100 i=2,nxp1
       do 100 j=2,nyp1
         cc=0.0d0
         do 105 io=1,3
         do 105 jo=1,3
105      cc=cc+rtmp(i+io-2,j+jo-2)
           if(abs(cc) .gt. 1.0e-04)then
            node=node+1
            iii(node)=i
            jjj(node)=j
           end if
100      continue 
       numnd=node


      itmax=5000
      bet=1.2
      errmax=1.0e-07 
!------------------------------------------
      cf=0.5/(hxi*hxi+hyi*hyi)
      do iter=1,itmax
       err=0.0d0 
         do node=1,numnd
           i=iii(node)
           j=jjj(node) 
           rold=r(i,j)
           r(i,j)=(1.0d0-bet)*r(i,j)+bet*cf*		&
             (hxi*hxi*(r(i+1,j)+r(i-1,j))+		&
            hyi*hyi*(r(i,j+1)+r(i,j-1))-rtmp(i,j))
           err=err+abs(rold-r(i,j))
         enddo


! set boundary terms...
        do j=1,nyp2
          r(nxp2,j)=r(nxp1,j)
          r(1,j)=r(2,j)
        enddo
        do i=1,nxp2
          r(i,nyp2)=r(i,nyp1)
          r(i,1)=r(i,2)
        enddo
         

! check convergence......
      if(err.lt.(errmax*float(numnd)))go to 300
     enddo
300   continue
        write(*,*)' DENSITY: ',err,' after ',iter,' iterations'
                

      return
      end
!-----------------------------------------------------------------------!
      subroutine gbd1(tmp1,tmp2,nxp2,nyp2)
	  use grid

	  implicit none

	  integer i,j,nxp2,nyp2
	  real*8:: tmp1(nxp2+2,nyp2+2),tmp2(nxp2+2,nyp2+2)

       do i=1,nxp2+2
!              { solid wall in y}
         tmp1(i,3)=tmp1(i,3)+tmp1(i,2)
         tmp1(i,4)=tmp1(i,4)+tmp1(i,1)
         tmp1(i,nyp1)=tmp1(i,nyp1)+tmp1(i,nyp2+2)
         tmp1(i,nyp2)=tmp1(i,nyp2)+tmp1(i,nyp2+1)

         tmp2(i,3)=tmp2(i,3)-tmp2(i,2)
         tmp2(i,4)=tmp2(i,4)-tmp2(i,1)
         tmp2(i,nyp1)=tmp2(i,nyp1)-tmp2(i,nyp2+2)
         tmp2(i,nyp2)=tmp2(i,nyp2)-tmp2(i,nyp2+1)
      enddo

       do i=1,nxp2+2
         tmp1(i,2)=tmp1(i,3)
         tmp1(i,1)=tmp1(i,4)
         tmp1(i,nyp2+2)=tmp1(i,nyp1)
         tmp1(i,nyp2+1)=tmp1(i,nyp2)
         tmp2(i,2)=-tmp2(i,3)
         tmp2(i,1)=-tmp2(i,4)
         tmp2(i,nyp2+2)=-tmp2(i,nyp1)
         tmp2(i,nyp2+1)=-tmp2(i,nyp2)
       enddo

       do j=1,nyp2+2
!              { solid wall in x}
         tmp1(3,j)=tmp1(3,j)-tmp1(2,j)
         tmp1(4,j)=tmp1(4,j)-tmp1(1,j)
         tmp1(nxp1,j)=tmp1(nxp1,j)-tmp1(nxp2+2,j)
         tmp1(nxp2,j)=tmp1(nxp2,j)-tmp1(nxp2+1,j)

         tmp2(3,j)=tmp2(3,j)+tmp2(2,j)
         tmp2(4,j)=tmp2(4,j)+tmp2(1,j)
         tmp2(nxp1,j)=tmp2(nxp1,j)+tmp2(nxp2+2,j)
         tmp2(nxp2,j)=tmp2(nxp2,j)+tmp2(nxp2+1,j)
       enddo

       do j=1,nyp2+2
         tmp1(2,j)=-tmp1(3,j)
         tmp1(1,j)=-tmp1(4,j)
         tmp1(nxp2+2,j)=-tmp1(nxp1,j)
         tmp1(nxp2+1,j)=-tmp1(nxp2,j)

         tmp2(2,j)=tmp2(3,j)
         tmp2(1,j)=tmp2(4,j)
         tmp2(nxp2+2,j)=tmp2(nxp1,j)
         tmp2(nxp2+1,j)=tmp2(nxp2,j)
       enddo


      return
      end

!---------------------------------------------------------------------------------------------------------------------
subroutine finit(nfronts,fp,ip,pt,icp,ine,ptcon,elcon,bptcon,belcon,elprop,maxpt,maxel,am,prop1,prop2,xcc,ycc)
	use front
	use flprop; use fparam
	use source

	  implicit none
! a subroutine to initialize nsphr spheres by dividing the
! sphere into four ``strips'' and setting each strip separately. 

      integer is,nfronts,nptot,netot,maxpt,maxel,itp,nps
	  integer i,j,k
	  real*8 pi,radin,prop1,prop2,dth,rad,th,am,xcc,ycc

	  integer:: icp(maxel,2), ine(maxel,2),bptcon(maxpt),belcon(maxel),ip(20)
      integer:: ptcon(maxpt),elcon(maxel)	
      real*8::	pt(maxpt,2),elprop(maxel,3),fp(20,6)       

      pi=4.*ATAN(1.)


      nfronts=1
      
      num_srce=0
      nptot=0
      netot=0
      do is=1,nfronts
      radin=0.5d0
      prop1=r1-r2   !prop1=2.6-1.0=1.6
      prop2=0.0d0
          fp(is,1)=radin
          fp(is,2)=xcc
          fp(is,3)=ycc
          fp(is,4)=prop1
          fp(is,5)=0.0d0
          fp(is,6)=0.0d0
  
!                         { itp=1 for a cylinder}  !Comienzo de archivo a leer
        nps=int(2.0d0*pi*radin/am)    !Se divide perimetro entre longitud ideal am
        dth=2.0d0*pi/float(nps)       !diametro entre puntos a tener, arcos
        rad=radin 
        do i=1,nps                                 
         th=dth*(float(i)-0.5d0) 
         pt(i+nptot,1)=xcc + rad*cos(th)   !coordenada en x
         pt(i+nptot,2)=ycc + rad*sin(th)   !coordenada en y
        enddo                                 !Final archivo de leer
 
       do i=1,nps                          ! Se comienzan los elementos del circulo
           icp(i+netot,1)=i    +nptot
           icp(i+netot,2)=i+1  +nptot
           ine(i+netot,1)=i-1  +netot
           ine(i+netot,2)=i+1  +netot 
             elprop(i+netot,1)=prop1 
             elprop(i+netot,2)=prop2
       enddo  
           icp(nps+netot,2)=1  +nptot
           ine(1+netot,1)=nps  +netot
           ine(nps+netot,2)=1  +netot

      ne=nps                   +netot 
      netot=ne
      np=nps                   +nptot 
      nptot=np
     
     enddo

      ffp=1
      lfp=np
      fep=np+1

      ffe=1
      lfe=ne
      fee=ne+1

! set connectivity

      do i = 1,maxpt-1       !SE comienza conexion entre puntos y elementos, maxpt=1000
       ptcon(i)=i+1        
       bptcon(i+1)=i
     enddo
       ptcon(lfp)=1
       bptcon(1)=lfp

      do i = 1,maxel-1        !maxel=1000
       elcon(i)=i+1
       belcon(i+1)=i
      enddo
       elcon(lfe)=1
       belcon(1)=lfe

      return
      end
!-----------------------------------------------------------------------!
      subroutine fstatis(pt,icp,ine,ptcon,elcon,bptcon,belcon,		&
                             cv,elprop,maxpt,maxel,area,xc,yc)

      use front
	  use fparam

	  implicit none


	  integer p1,p2,k,kk,maxel,maxpt
	  real*8 area,xm,ym,dx,dy,xa,ya,xc,yc
      integer:: icp(maxel,2), ine(maxel,2),ptcon(maxpt),elcon(maxel)
      real*8:: cv(maxel,2), pt(maxpt,2), elprop(maxel,3)
      integer bptcon(maxpt),belcon(maxel)
      
      area=0.0d0 
!      cx=0.0d0
!      cy=0.0d0
      xm=0.0d0
      ym=0.0d0
      k=ffe 
      do 10 kk=1,ne
      k=elcon(k) 
      p1=icp(k,1)
      p2=icp(k,2)   
      area=area+0.25*((pt(p2,2)-pt(p1,2))*(pt(p2,1)+pt(p1,1))-		&
                     (pt(p2,1)-pt(p1,1))*(pt(p2,2)+pt(p1,2)) )    

       dy=pt(p2,2)-pt(p1,2)
       dx=pt(p2,1)-pt(p1,1)
       xa=0.5*(pt(p2,1)+pt(p1,1))
       ya=0.5*(pt(p2,2)+pt(p1,2))
      xm=xm+0.5*dy*xa**2
      ym=ym-0.5*dx*ya**2

!      cx=cx+cv(k,1)
!      cy=cy+cv(k,2)
10    continue 
      xc=xm/area
      yc=ym/area
!      write(*,*)'Area: ',area,' xc,yc: ',xc,yc
      return
      end
!-----------------------------------------------------------------------!

Subroutine simplec(p,u,v,xerr,dt,r,ro,fx,fy,u1,v1,Div,xcc,ycc,color,time,gama)

use grid
use bounds
use flprop
implicit none

integer i,j,k,max_iter,itmax,it,ei,ej,c,max_iter2,kk,nminx,nmaxx,nminy,nmaxy
real*8 dx,dy,dz,dlv,dt
real*8 ue,uw,us,un,dbkx,dbky,Div,xerr,pi
real*8 Se,Sw,Sn,Ss,Re,residual,tolerance,tolerance2,gama
real*8 Fw,Fe,Fs,Fn,De,Ds,Dn,Dw,dF,xcc,ycc,fac,itau,Ifm,ms,Im,ucs,vcs,time,omegas

real*8  u(nx1,ny1+1),v(nx1+1,ny1),p(nx1+1,ny1+1),Bf(nx1+1,ny1+1),color(nx1+1,ny1+1)
real*8  u1(nx1,ny1+1),v1(nx1+1,ny1),r(nx1+1,ny1+1),ro(nx1+1,ny1+1),fx(nx1+1,ny1+1),fy(nx1+1,ny1+1)
real*8, allocatable:: aP(:,:),aE(:,:),pp(:,:)
real*8, allocatable:: aW(:,:),aN(:,:),aS(:,:)
real*8, allocatable:: sP(:,:),du(:,:),dv(:,:),uso(:,:),vso(:,:)
real*8, allocatable:: x(:),xc(:)
real*8, allocatable:: y(:),yc(:)

allocate(aP(2:nx1,2:ny1),aE(2:nx1,2:ny1),pp(nx1+1,ny1+1))
allocate(aW(2:nx1,2:ny1),aN(2:nx1,2:ny1),aS(2:nx1,2:ny1))
allocate(sP(2:nx1,2:ny1),du(nx1+1,ny1+1),dv(nx1+1,ny1+1))
allocate(uso(nx1+1,ny1+1),vso(nx1+1,ny1+1))
allocate(x(nx1),xc(nx1+1),y(ny1),yc(ny1+1))

max_iter=500
max_iter2=100
tolerance=xerr*100.0d0
tolerance2=xerr

!Construcci\F3n de la Malla

dx=xl/dfloat(nx1-1)
dy=yl/dfloat(ny1-1)
dz=1.0d0
dlv=dx*dy*dz


Sw=dy*dz
Se=dy*dz
Ss=dx*dz
Sn=dx*dz

pi=acos(-1.0d0)


du=0.0d0
dv=0.0d0
call Mesh1D(xc,x,0.d0,xl,nx1)
call Mesh1D(yc,y,0.d0,yl,ny1)


nminx=int((xcc-1.0d0)/dx)
if (nminx < 1)nminx=1
nmaxx=int((xcc+1.0d0)/dx)
if (nmaxx > nx1+1)nmaxx=nx1+1
nminy=int((ycc-1.0d0)/dy)
if (nminy < 1)nminy=1
nmaxy=int((ycc+1.0d0)/dy)
if (nmaxy > ny1+1)nmaxy=ny1+1

ms=0.0d0
Im=0.0d0
do i=nminx,nmaxx
do j=nminy,nmaxy
Im=Im+dx*dy*((xc(i)-xcc)**2.0d0+(yc(j)-ycc)**2.0d0)*(1.0d0-color(i,j))*r(i,j)
ms=ms+dx*dy*(1.0d0-color(i,j))*r(i,j)
enddo
enddo
!ms=0.5d0**2*pi*r1
!Im=0.5d0*(pi*0.5d0**2)*0.5d0**2*r1
write(*,*)ms,Im!maxval((1.0d0-color)*(xc-xcc)),maxval((1.0d0-color)*(yc-ycc)),maxval(abs(color))
write(*,*)xcc,ycc,maxval(r),minval(r)
c=1
Div=1.0d0

itau=-2.0
Ifm=1.0
!Ciclo de m\E9todo simplec
do while ((c < max_iter2) .and. (Div > tolerance))

!*******************************************************************************************************************************
!Ecuaci\F3n para u

aP=0.0d0;	aW=0.0d0;	aE=0.0d0;	aS=0.0d0;	aN=0.0d0;	sP=0.0d0
ei=ubound(u,1)-1
ej=ubound(u,2)-1


!C\E1lculo de coeficientes

do i=2, ei
	do j=2, ej

		ue=0.5d0*(u(i,j)+u(i+1,j))*Se*Ifm
		uw=0.5d0*(u(i-1,j)+u(i,j))*Sw*Ifm
		un=0.5d0*(v(i,j)+v(i+1,j))*Sn*Ifm
		us=0.5d0*(v(i,j-1)+v(i+1,j-1))*Ss*Ifm

		De=gama*Se/dx
		Dw=gama*Sw/dx
		Dn=gama*Sn/dy
		Ds=gama*Ss/dx

		aE(i,j)=De-0.5d0*ue
		aW(i,j)=Dw+0.5d0*uw
		aN(i,j)=Dn-0.5d0*un
		aS(i,j)=Ds+0.5d0*us
		
		aP(i,j)=aE(i,j)+aW(i,j)+aN(i,j)+aS(i,j)+dlv/dt
		sP(i,j)=u(i,j)*dlv/dt+((p(i,j)-p(i+1,j))/dx)*dlv							
						
						
	enddo
enddo


!Condiciones de frontera

dbkx=0.0d0
dbky=1.0d0

!Este
aP(ei,2:ej)=aP(ei,2:ej)+dbkx*aE(ei,2:ej)
aE(ei,2:ej)=0.0d0

!Oeste
aP(2,2:ej)=aP(2,2:ej)+dbkx*aW(2,2:ej)
aW(2,2:ej)=0.0d0

!Norte
aP(2:ei,ej)=aP(2:ei,ej)-dbky*aN(2:ei,ej)
aN(2:ei,ej)=0.0d0

!Sur
aP(2:ei,2)=aP(2:ei,2)-dbky*aS(2:ei,2)
aS(2:ei,2)=0.0d0

do i=2,ei
	do j=2, ej
	du(i,j)=Se/(aP(i,j)-aE(i,j)-aW(i,j)-aN(i,j)-aS(i,j))
	enddo
enddo

call Gauss_TDMA2D(u1,ei,ej,aP,aE,aW,aN,aS,sP,nx1,ny1,max_iter,tolerance2,residual)


u1(:,ny1+1)=u1(:,ny1)
u1(:,1)=u1(:,2)
          

!*******************************************************************************************************************************
!Ecuaci\F3n para v


aP=0.0d0;	aW=0.0d0;	aE=0.0d0;	aS=0.0d0;	aN=0.0d0;	sP=0.0d0
ei=ubound(v,1)-1
ej=ubound(v,2)-1


!C\E1lculo de coeficientes

do i=2, ei
	do j=2, ej

		ue=0.5d0*(u(i,j)+u(i,j+1))*Se*Ifm
		uw=0.5d0*(u(i-1,j)+u(i-1,j+1))*Sw*Ifm
		un=0.5d0*(v(i,j)+v(i,j+1))*Sn*Ifm
		us=0.5d0*(v(i,j)+v(i,j-1))*Ss*Ifm

		
		De=gama*Se/dx
		Dw=gama*Sw/dx
		Dn=gama*Sn/dy
		Ds=gama*Ss/dx

		aE(i,j)=De-0.5d0*ue
		aW(i,j)=Dw+0.5d0*uw
		aN(i,j)=Dn-0.5d0*un
		aS(i,j)=Ds+0.5d0*us
		
		aP(i,j)=aE(i,j)+aW(i,j)+aN(i,j)+aS(i,j)+dlv/dt
		sP(i,j)=v(i,j)*dlv/dt+((p(i,j)-p(i,j+1))/dy)*dlv					

	enddo
enddo

!Condiciones de frontera

dbkx=1.0d0
dbky=0.0d0

!Este
aP(ei,2:ej)=aP(ei,2:ej)+dbkx*aE(ei,2:ej)
aE(ei,2:ej)=0.0d0

!Oeste
aP(2,2:ej)=aP(2,2:ej)+dbkx*aW(2,2:ej)
aW(2,2:ej)=0.0d0

!Norte
aP(2:ei,ej)=aP(2:ei,ej)+dbky*aN(2:ei,ej)
aN(2:ei,ej)=0.0d0

!Sur
aP(2:ei,2)=aP(2:ei,2)+dbky*aS(2:ei,2)
aS(2:ei,2)=0.0d0



do i=2,ei
	do j=2, ej
	dv(i,j)=Sn/(aP(i,j)-aE(i,j)-aW(i,j)-aN(i,j)-aS(i,j))
	enddo
enddo

call Gauss_TDMA2D(v1,ei,ej,aP,aE,aW,aN,aS,sP,nx1,ny1,max_iter,tolerance2,residual)


!*************************************************************************************************************
!Ecuaci\F3n para pp

aP=0.0d0;	aW=0.0d0;	aE=0.0d0;	aS=0.0d0;	aN=0.0d0;	sP=0.0d0; pp=0.0d0
ei=ubound(p,1)-1
ej=ubound(p,2)-1

!C\E1lculo de coeficientes

do i=2, ei
	do j=2, ej

		ue=u1(i,j)
		uw=u1(i-1,j)
		un=v1(i,j)
		us=v1(i,j-1)
			
		aE(i,j)=du(i,j)*Se
		aW(i,j)=du(i-1,j)*Sw
		aN(i,j)=dv(i,j)*Sn
		aS(i,j)=dv(i,j-1)*Ss

		aP(i,j)=aE(i,j)+aW(i,j)+aN(i,j)+aS(i,j)
		sP(i,j)=-ue*Se+uw*Sw-un*Sn+us*Ss

	enddo
enddo

call Gauss_TDMA2D(pp,ei,ej,aP,aE,aW,aN,aS,sP,nx1,ny1,max_iter,tolerance2,residual)


!*********************************************************************************************************************
!Correcci\F3n de la presi\F3n

p=p+pp


!Correcci\F3n de las velocidades
do i=2,nx1-1
do j=2,ny1
u1(i,j)=u1(i,j)+du(i,j)*(pp(i,j)-pp(i+1,j))
enddo
enddo

u1(:,ny1+1)=u1(:,ny1)
u1(:,1)=u1(:,2)

do i=2,nx1
do j=2,ny1-1
v1(i,j)=v1(i,j)+dv(i,j)*(pp(i,j)-pp(i,j+1))
enddo
enddo


!C\E1lculo de la divergencia
Div=maxval(abs(sp(2:ei,2:ej)))
c=c+1


enddo   !Fin del ciclo de simplec
uso=0.0d0
vso=0.0d0


do i=nminx,nmaxx
	do j=nminy,nmaxy
		u1(i,j)=0.5d0*(color(i,j)+color(i+1,j))*u1(i,j)+(1.0d0-0.5d0*(color(i,j)+color(i+1,j)))*uso(i,j)
		v1(i,j)=0.5d0*(color(i,j)+color(i,j+1))*v1(i,j)+(1.0d0-0.5d0*(color(i,j)+color(i,j+1)))*vso(i,j)
	enddo
enddo

write(33,*)time,ucs,vcs,omegas

deallocate(uso,vso)
deallocate(aP,aE,pp)
deallocate(aW,aN,aS)
deallocate(sP,du,dv)
deallocate(x,xc,y,yc)
write(*,*)'DIVERGENCE: ',c,Div

end subroutine 


!************************************************


real*8 function calcResidual(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx1,ny1)
implicit none
integer bi,ei,bj,ej,i,j,nx1,ny1
real*8 phi(ei+1,ej+1),aP(2:nx1,2:ny1),aE(2:nx1,2:ny1),aW(2:nx1,2:ny1),aN(2:nx1,2:ny1),aS(2:nx1,2:ny1),sP(2:nx1,2:ny1)
real*8,allocatable :: acum(:,:)
real*8 residual,NINV
allocate(acum(2:ei,2:ej))
bi=2; bj=2
acum=0
NINV = 1.0d0 / dfloat(ei*ej)
do i=bi,ei
do j=bj,ej
acum(i,j) = aE(i,j) * phi(i+1,j) + aW(i,j) * phi(i-1,j) + aN(i,j) * phi(i,j+1) + &
aS(i,j) * phi(i,j-1) + sp(i,j) - aP(i,j) * phi(i,j)
end do
end do
residual = sqrt( NINV * sum(acum * acum) )
calcResidual=residual
end function


!************************************************************************************************************************

!Gauss-tdma

!***********************************************************************************************************************


subroutine TDMA(x,a,b,c,d,n)

implicit none
 integer n,k
 real*8 a(2:n),b(2:n),c(2:n),d(2:n),x(2:n),m

 do k=3,N
  m=a(k)/b(k-1)
  b(k)=b(k)-m*c(k-1)
  d(k)=d(k)-m*d(k-1)
 end do

 x(n)=d(n)/b(n)

 do k=n-1,2,-1
  x(k)=(d(k)-c(k)*x(k+1))/b(k)
 end do

end subroutine

subroutine lineX_2D(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx1,ny1)

implicit none
integer bi,ei,bj,ej,i,j,nx1,ny1
real*8 phi(ei+1,ej+1),aP(2:nx1,2:ny1),aE(2:nx1,2:ny1),aW(2:nx1,2:ny1),aN(2:nx1,2:ny1),aS(2:nx1,2:ny1),sP(2:nx1,2:ny1)
real*8, allocatable:: a(:),b(:),c(:),d(:)
bi=2; bj=2
allocate(a(2:ei),b(2:ei),c(2:ei),d(2:ei))

do j=bj,ej
do i=bi,ei
	a(i)=-aW(i,j)
	b(i)=aP(i,j)
	c(i)=-aE(i,j)
	d(i)=sp(i,j) + aN(i,j) * phi(i,j+1) + aS(i,j) * phi(i,j-1)
end do
 call TDMA(phi(bi:ei,j), a, b, c ,d ,ei)
end do

end subroutine

subroutine lineY_2D(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx1,ny1)

implicit none
integer bi,ei,bj,ej,i,j,nx1,ny1
real*8 phi(ei+1,ej+1),aP(2:nx1,2:ny1),aE(2:nx1,2:ny1),aW(2:nx1,2:ny1),aN(2:nx1,2:ny1),aS(2:nx1,2:ny1),sP(2:nx1,2:ny1)
real*8, allocatable:: a(:),b(:),c(:),d(:)
bi=2; bj=2

allocate(a(2:ej),b(2:ej),c(2:ej),d(2:ej))


do i=bi,ei
do j=bj,ej
	a(j)=-aS(i,j)
	b(j)=aP(i,j)
	c(j)=-aN(i,j)
	d(j)=sp(i,j) + aE(i,j) * phi(i+1,j) + aW(i,j) * phi(i-1,j) 
end do
 call TDMA(phi(i,bj:ej), a, b, c ,d ,ej)
end do

end subroutine


subroutine Gauss_TDMA2D(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx1,ny1,max_iter,tolerance,residual)

implicit none

integer bi,ei,bj,ej,i,j,nx1,ny1,count_iter,max_iter
real*8 phi(ei+1,ej+1),aP(2:nx1,2:ny1),aE(2:nx1,2:ny1),aW(2:nx1,2:ny1),aN(2:nx1,2:ny1),aS(2:nx1,2:ny1),sP(2:nx1,2:ny1)
real*8 residual,tolerance
	
	interface
		real*8 function calcResidual(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx1,ny1)
		implicit none
		integer bi,ei,bj,ej,i,j,nx1,ny1
		real*8 phi(ei+1,ej+1),aP(2:nx1,2:ny1),aE(2:nx1,2:ny1),aW(2:nx1,2:ny1),aN(2:nx1,2:ny1),aS(2:nx1,2:ny1),sP(2:nx1,2:ny1)
		end function
	end interface

count_iter=0;  residual=1.0d0

do while((count_iter <= max_iter).and.(residual > tolerance)) 
call lineX_2D(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx1,ny1)
call lineY_2D(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx1,ny1)
residual = calcResidual(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx1,ny1)
count_iter=count_iter+1
end do

end subroutine
!***********************************************************************************************************************

Subroutine Forceh(pt,icp,elcon,maxel,maxpt,norm,u,v,p,Fhx,Fhy,cnt,gama,nxp2,nyp2,r,xcc,ycc)

use grid; use bounds
use front; 
use flprop; use fparam

implicit none


real*8, allocatable:: dudx(:,:),dvdx(:,:),dudy(:,:),dvdy(:,:)
real*8 cnt(maxel,2), norm(maxel,2)
integer i,j,iih,jjh,irh,jrh,i1,j1,k,kk,maxel,maxpt,nxp2,nyp2
real*8 p0x,p0y,pi,pd2,cnst,xc,yc,xth,yth
real*8 drxh,dryh,dx,dy,ds,gama
real*8 sigma11,sigma12,sigma22,v1,v2,u1,u2,mu,Fhx,Fhy,xcc,ycc
real*8 u(nxp2,nyp2),v(nxp2,nyp2),p(nxp2,nyp2),r(nxp2,nyp2)
real*8 pt(maxpt,2)
integer icp(maxel,2), elcon(maxel)
allocate(dudx(nxp2,nyp2),dvdx(nxp2,nyp2),dudy(nxp2,nyp2),dvdy(nxp2,nyp2))

pi=acos(-1.0d0)

pd2=0.5*pi
cnst=(0.25)*(0.25)

do j=2,nyp1
	do i=2,nxp1
		dudx(i,j)=0.0d0
		dvdx(i,j)=0.0d0
		dudy(i,j)=0.0d0
		dvdy(i,j)=0.0d0
	enddo
enddo

do j=2,nyp1
	do i=2,nxp1
		dudx(i,j)=gama*hxi*(u(i,j)-u(i-1,j))
		dvdy(i,j)=gama*hyi*(v(i,j)-v(i,j-1))
		u1=0.5*(0.5*(u(i-1,j+1)+u(i-1,j))+0.5*(u(i,j+1)+u(i,j)))
		u2=0.5*(0.5*(u(i-1,j-1)+u(i-1,j))+0.5*(u(i,j-1)+u(i,j)))
		v1=0.5*(0.5*(v(i,j)+v(i+1,j))+0.5*(v(i,j-1)+v(i+1,j-1)))
		v2=0.5*(0.5*(v(i-1,j)+v(i,j))+0.5*(v(i-1,j-1)+v(i,j-1)))		
		dvdx(i,j)=gama*hxi*(v1-v2)
		dudy(i,j)=gama*hyi*(u1-u2)
	enddo
enddo

Fhx=0.0d0
Fhy=0.0d0

k=ffe
do kk=1,ne
	xc=cnt(k,1)
	yc=cnt(k,2)
	p0x=pt(icp(k,1),1)
	p0y=pt(icp(k,1),2)
	dx=pt(icp(k,2),1)-p0x
	dy=pt(icp(k,2),2)-p0y
	ds=sqrt(dx*dx+dy*dy)

	if (xc <= xcc)xth=xc-0.5/hxi
	if (xc > xcc)xth=xc+0.5/hxi
	if (yc <= ycc)yth=yc-0.5/hyi
	if (yc > ycc)yth=yc+0.5/hyi
	irh=1+int(nx*xth/xl)
	jrh=1+int(ny*yth/yl)

	sigma11=						&
					(2.0d0*dudx(irh,jrh)-p(irh,jrh))
	sigma22=							&
					(2.0d0*dvdy(irh,jrh)-p(irh,jrh))
	sigma12=							&
					(dudy(irh,jrh)+dvdx(irh,jrh))							

	
	Fhx=Fhx+(sigma11*norm(k,1)+sigma12*norm(k,2))*ds
	Fhy=Fhy+(sigma12*norm(k,1)+sigma22*norm(k,2))*ds
	
	k=elcon(k)

enddo

end subroutine

!************************************************************
Subroutine Mesh1D(xc,x,x0,xl,nx)
implicit none
integer i,j,nx
real*8 x0,xl,dx
real*8 x(1:nx),xc(1:nx+1)
dx=(1.0d0)/dfloat(nx-1)
do i=1,nx
	x(i)=dfloat(i-1)*dx
	x(i)=x0+(xl-x0)*x(i)
end do


xc(1)=x(1); xc(nx+1)=x(nx)
do i=2,nx
xc(i)=(x(i)+x(i-1))*0.5
end do

End Subroutine


subroutine printout(pt,ptcon,iptmp,istore,icp,elcon,				&
    maxpt,maxel,r,p,u,v,time,nnm,nout,lm,nxp2,nyp2)

	use grid; use bounds
	use front; 
	use flprop; use fparam
    real*8 time	
    real*8  pt(maxpt,2)
    integer ptcon(maxpt),elcon(maxel),iptmp(maxpt),icp(maxel,2),istore(maxel,2)
    real*8 r(nxp2,nyp2),p(nxp2,nyp2),u(nxp2,nyp2),v(nxp2,nyp2)
    character*20 :: outfile, nnm
    character*1 :: num(10)
    
	num=(/'0','1','2','3','4','5','6','7','8','9'/)       
    outfile=nnm     
    m1=mod(nout,10)+1
    m0=mod((nout/10),10)+1
    outfile(lm:lm+3)='frd'
    outfile(lm+3:lm+4)=num(1)     
    outfile(lm+4:lm+5)=num(m0) 
    outfile(lm+5:lm+6)=num(m1)
    outfile(lm+6:lm+10)='.vtk'
    open(8,file=outfile,status='new') 
    
! Write the Front File
    write(8,10)
    write(8,11)time
    write(8,12)
    write(8,13)
    write(8,14)np
    10 format('# vtk DataFile Version 3.0d0')
    11 format('grid, time ',F16.8)
    12 format('ASCII')
    13 format('DATASET POLYDATA')
    14 format('POINTS ',I17,' float' )
    
    zzero=0.0d0
    k=ffp
    do kk=1,np
        xx=pt(k,1)
        pt(k,1)=xx-xl*floor(xx/xl)
        yy=pt(k,2)
        pt(k,2)=yy-yl*floor(yy/yl)          
        write(8,255) pt(k,1),pt(k,2),zzero
        255 format(' ',3e13.5)
        k=ptcon(k)
    end do
    
    k=ffp
    do kk=1,np
        iptmp(k)=kk-1
        k=ptcon(k)
    end do
    
    numLinks=0
    k=ffe
    do 288 kk=1,ne
        i=icp(k,1)
        ii=icp(k,2)
        s=sqrt((pt(i,1)-pt(ii,1))**2+(pt(i,2)-pt(ii,2))**2)
        if(s < 0.25*xl)then
            numLinks=numLinks+1
            istore(numLinks,1)=iptmp(i)
            istore(numLinks,2)=iptmp(ii)
        endif
        k=elcon(k)
        continue
    288 END DO
    
    write(8,16) numLinks,3*numLinks
    16 format('LINES ',2I17)
    
    do kk=1,numLinks        
        write(8,266) istore(kk,1),istore(kk,2)
        266 format(' 2',2I6)
    end do
    
    close(8,status = 'keep')
!---------------------------------------------
    outfile=nnm 
    outfile(lm:lm+3)='vel'
    outfile(lm+3:lm+4)=num(1)     
    outfile(lm+4:lm+5)=num(m0) 
    outfile(lm+5:lm+6)=num(m1)
    outfile(lm+6:lm+10)='.vtk'
    open(8,file=outfile,status='new')       
    write(8,10)
    write(8,11)time
    write(8,12)
    
    write(8,23)
    23 format('DATASET RECTILINEAR_GRID')
    write(8,24)(nx+1),(ny+1),1
    24 format('DIMENSIONS ',I5,I5,I5)
    write(8,111)(nx+1)
    111 format('X_COORDINATES',I5,' float')
    write(8,110)((float(i-1)/hxi), i=1,nx+1)
    write(8,112)(ny+1)
    112 format('Y_COORDINATES',I5,' float')
    write(8,110)((float(j-1)/hyi), j=1,ny+1)
    write(8,113)1
    113 format('Z_COORDINATES',I5,' float')
    write(8,110)0.0d0
    110 format(e14.5)
    
    write(8,26)(nx+1)*(ny+1)
    26 format('POINT_DATA ',I17)
    write(8,21)
    21 format('VECTORS velocity float')
    do j=1,ny+1       
        write(8,310)( (0.5*(u(i,j)+u(i,j+1))),			&
        (0.5*(v(i,j)+v(i+1,j))),i=1,nx+1)
    end do
    310 format(e14.5,e14.5,' 0')
    
    write(8,66)nx*ny
    66 format('CELL_DATA ',I17)
    
    write(8,67)'density'
    67 format('SCALARS',A20,' float 1')
    write(8,99)
    99 format('LOOKUP_TABLE default')
    do j=2,ny+1       
        write(8,210)(r(i,j),i=2,nx+1)
    end do
    210 format(e14.5)
    
    write(8,67)'pressure'
    write(8,99)
    do j=2,ny+1       
        write(8,210)(p(i,j),i=2,nx+1)
    end do
    
    close(8,status='keep')
    
    nout=nout+1
    return
    end


subroutine fcurv(cv,t1,t2,cnt,pt,icp,ine,ptcon,elcon,maxel,maxpt,norm)

	  use front
	  use fparam
	  
	  implicit none

	  integer k,kk,maxpt,maxel
	  real*8 s1x,s1y,s2y,s2x,p0x,p0y,p1x,p1y,pm1x,pm1y,xx,yy,pp1x,pp1y,s1,s2,tgx,tgy

	  real*8::pt(maxpt,2),t1(maxel,2),t2(maxel,2),cv(maxel,2),cnt(maxel,2),norm(maxel,2)
      integer:: icp(maxel,2), ine(maxel,2) 
      integer ptcon(maxpt),elcon(maxel)   
  

      k=ffe    !ffe=1
      do 10 kk=1,ne 
		xx=(xmv+0.5)*fxl    !fxl=xl
		yy=(xmv+0.5)*fyl
		p0x=pt(icp(k,1),1)   !icp(k,1) devuelve un i
		p0y=pt(icp(k,1),2)
		p1x= p0x+dmod((pt(icp(k       ,2),1)-p0x+xx),fxl)-0.5*fxl
		p1y= p0y+dmod((pt(icp(k       ,2),2)-p0y+yy),fyl)-0.5*fyl
		pm1x=p0x+dmod((pt(icp(ine(k,1),1),1)-p0x+xx),fxl)-0.5*fxl
		pm1y=p0y+dmod((pt(icp(ine(k,1),1),2)-p0y+yy),fyl)-0.5*fyl
		pp1x=p0x+dmod((pt(icp(ine(k,2),2),1)-p0x+xx),fxl)-0.5*fxl
		pp1y=p0y+dmod((pt(icp(ine(k,2),2),2)-p0y+yy),fyl)-0.5*fyl
		cnt(k,1)=(-pm1x+9.0d0*p0x+9.0d0*p1x-pp1x)*0.0d0625
		cnt(k,2)=(-pm1y+9.0d0*p0y+9.0d0*p1y-pp1y)*0.0d0625

		s1x=(-2.0d0*pm1x-3.0d0*p0x+6.0d0*p1x-pp1x)/6.0d0
		s1y=(-2.0d0*pm1y-3.0d0*p0y+6.0d0*p1y-pp1y)/6.0d0
		s2x=(pm1x-6.0d0*p0x+3.0d0*p1x+2.0d0*pp1x)/6.0d0
		s2y=(pm1y-6.0d0*p0y+3.0d0*p1y+2.0d0*pp1y)/6.0d0
		s1=sqrt(s1x**2+s1y**2)
		s2=sqrt(s2x**2+s2y**2)
		t1(k,1)=s1x/s1
		t1(k,2)=s1y/s1  
		t2(k,1)=s2x/s2
		t2(k,2)=s2y/s2 
	   
		tgx=0.5*(s1x/s1+s2x/s2)
		tgy=0.5*(s1y/s1+s2y/s2)
		norm(k,1)=tgy
		norm(k,2)=-tgx  
      k=elcon(k)
10    continue
            
      k=ffe    
      do 20 kk=1,ne
       cv(k,1)=0.5*(t2(k,1)+t1(ine(k,2),1)-t1(k,1)-t2(ine(k,1),1) )
       cv(k,2)=0.5*(t2(k,2)+t1(ine(k,2),2)-t1(k,2)-t2(ine(k,1),2) ) 

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
!       ds=sqrt((pt(icp(k,2),1)-pt(icp(k,1),1))**2 +
!     1                (pt(icp(k,2),2)-pt(icp(k,1),2))**2) 
!       cu=sqrt( cv(k,1)**2 +cv(k,2)**2)
!       curv=cu/ds
!       xxh=0.5*(pt(icp(k,1),1)+pt(icp(k,2),1)) 
!       yyh=0.5*(pt(icp(k,1),2)+pt(icp(k,2),2))
!       write(*,100)k,cu,ds,curv,cnt(k,1),xxh,cnt(k,2),yyh 
!c100    format(' !: ',i5,7e13.5)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      k=elcon(k)
20    continue 

       return
       end                                 

