Program mesh


integer nx,nxm,n1,n2
real, allocatable :: xc(:), x(:)
real x0,xl
real xm1,xm2

  n1=5         !Número de celdas en el primer tramo
  nxm=10        !Número de celdas en el tercer tramo
  n2=5          !Número de celdas en el segundo tramo

  nt=n1+nxm+n2     !Celdas totales
  x0=0.0      !punto inicial
  xl=10.0     !punto final
  xm1= 5.     !punto de inicia de seccion estrecha
  xm2= 7.     !punto final de seccion estrecha

  allocate (xc(0:nt+1),x(0:nt)) 
  
  call Mesh1D_unequal(xc,x,x0,xl,xm1,xm2,n1,n2,nxm,nt)

      
    do i=0,nt
      write(*,*) i,x(i), xc(i)
    enddo

    do i=0,nt+1
      write(*,*) i,xc(i)
    enddo


end program mesh


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

  Subroutine Mesh1D_unequal(xc,x,x0,xl,xm1,xm2,n1,n2,nxm,nt)
    integer i,j,n1,n2,nxm,nt
    real*4 x0,xl,dx,dxm,xm1,xm2,q
    real*4 x(0:nt),xc(0:nt+1)
    
    dx1=(1.0)/dfloat(n1)
    dxm=(1.0)/dfloat(nxm)
    dx2=(1.0)/dfloat(n2)
    
    x=0.
    do i=0,n1
       q=dfloat(i)*dx1
       x(i)=x0+(xm1-x0)*q
    end do
    
    do i=1,nxm
       q=dfloat(i)*dxm
       x(i+n1)=xm1+(xm2-xm1)*q
    end do  
    
    do i=1,n2
       q=dfloat(i)*dx2
       x(i+n1+nxm)=xm2+(xl-xm2)*q
    enddo
    
    xc(0)=x(0); xc(nt+1)=x(nt)
    
    do i=1,nt
       xc(i)=(x(i)+x(i-1))*0.5
    end do
  End Subroutine Mesh1D_unequal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
