
!!!!*******************************************************************************************************************************
Program Main
  implicit none
 ! integer:: n,i,col,maxcols
 ! integer::error
  ! real,allocatable::a(:,:),alpha(:),Lift(:),Drag(:), z(:)
  real,allocatable::U(:),Pitch(:),omgmax(:),R_bp(:),chord(:),twist(:),dr(:)
  real::inSpeed,delSpeed,outspeed,OMax,Rwt,rt,r_t,Nb
  real::Pi
  integer::n,i,NumCases,NE
  CHARACTER(len=50)::Name
  CHARACTER(len=100)::Dir
  Dir='/home/calovi/Documentos/Alan_BEM_Method/BEM_Method/ ' !Dirección del archivo .dat
  !Es importante que el archivo .dat se elimine la última linea que dice EOT, Ya que puede traer problemas con la extracción de la información
  Name='perfil.dat'                    !Nombre del archivo .dat con las características aerodinámicas del perfil
  Name=ADJUSTL(Name(1:LEN_TRIM(Name)))  
  Dir=ADJUSTL(Dir(1:LEN_TRIM(Dir)))
  !DEfieniendo paramétros
  Pi=acos(-1.0)

  !características de la turbina eólica
  inSpeed=5.0 !velocidad de entrada m/s
  delSpeed=0.5 !incremento de velocidad  m/s
  outSpeed=20.0!velocidad de salida m/s
  OMax=71.03 !Velocidad de rotación máxima
  n=int((outSpeed-inSpeed)/delSpeed)+1 !Para el número de casos
 
  NE=10  !Número de elementos
  allocate(U(n),Pitch(n),omgmax(n),R_bp(NE+1),chord(NE+1),twist(NE+1),dr(NE+1))  
  U(1)=inSpeed
  Pitch=0.
  omgmax(1:n)=OMax
  do i=1,n+1
     U(i+1)=U(i)+delSpeed
 enddo
! NumCases=Int((outSpeed - inSpeed)/delSpeed + 1)
 !Geometría del rotor de la turbina eólica
 Nb=3. !número de aspas
 Rwt=5.023 !Radio del rotor  [m]
 rt=0.84 !radio del hub [m]
 r_t=rt/Rwt !radio adimensional, miu=r/R
 R_bp(1)=r_t
 do i=2,NE+1
    R_bp(i)= R_bp(i)+((2./Pi)*acos(1-(i-1.)/float(NE))-0.1)*(1.-r_t)/(0.9)
 enddo
 dr(1)=r_t
 do i=1,NE+1
    dr(i+1)=R_bp(i+1)-R_bp(i)
 enddo
 !utilizando curve fitting de Gnuplot, se puede parametrizar linealmente (f(x)=a*x+b) la distribución de Chord (c,d) y twist (a,b) 
 twist(:)= -0.8596424209*(Rwt*R_bp(:))**3 + 10.2546334225*(Rwt*R_bp(:))**2 - 44.1851148314*Rwt*R_bp(:)+70.9152327719 !(-9.168+.7734)*Rwt*R_bp(:)+ (39.5373+2.423)
 !La cuerda tambien se puede representar como una función lineal del tipo f(x)=c*x+d, para propósitos de simplicidad se mantiene constante
 chord(:)=0.457
 call Bem_method(Name,Dir,Nb,R_bp,Rwt,NE,twist,Pitch,chord,dr,U,omgmax,n)
EndProgram Main
!!!!*******************************************************************************************************************************
subroutine Bem_method(Name,Dir,Nb,R_bp,Rwt,NE,twist,Pitch,chord,dr,U,omgmax,nz)
  CHARACTER(len=50)::Name
  CHARACTER(len=100)::Dir
  real::Rwt,rt,r_t,L,D,Caxi,Ct,difa,difap,pitch1,twist1,lamda
  real::a,ap,a1,ap1,f,CN,Nb
  real::Pi
  integer::nz,j,k,m,NE,axial1,tangen1

  real:: U(nz),Pitch(nz),omgmax(nz),Q(nz),CQ(nz),R_bp(NE+1),chord(NE+1),twist(NE+1),dr(NE+1),axial(NE+1,nz)
  real::tangen(NE+1,nz), CP(nz), dQ(NE+1,nz)
  axial=0.
  tangen=0.
  dQ=0.
  Cp_local=0.
  Power_coeff=0.
  !DEfieniendo paramétros
  Pi=acos(-1.0)
  alpha = 10.
  do j=1,nz
     lamda=(Rwt*omgmax(j)*2.*Pi)/(U(j)*60.)
     pitch1=pitch(j)
     
     do k=1,NE+1
        twist1=twist(k)
        !constantes
        a=0.4 !0.4        
        ap=1.0/4.0 !0.3
        difa=1.0
        difap=1.0
        m=1
        do while ((difa .GE. 1.0e-6) .and. (difap .GE. 1.0e-6) .and. (m.LT.2500))
           !ángulo de flujo
          ! phi=(180.0/Pi)*atan((1.0-a)/(lamda*R_bp(k)*(1.0+ap)))
           phi=(180.0/Pi)*atan((1.0/(lamda*R_bp(k)))*((1.0-a)/(1.0-ap)))
           !Pérdida en punta
           !f=(2.0/Pi)*acos(exp(-((Nb*(1.0-R_bp(k)))/(2.0*R_bp(k)*sin(phi*Pi/180.0)))))
           f=(2.0/Pi)*acos(exp(-(Nb/2.0)*((1.0-R_bp(k))/R_bp(k))*sqrt(1.0 + (lamda*R_bp(k))/(1.0-a)**2)))           
           !ángulo de ataque
           alpha=phi-twist1-pitch1
           !solidez       
           sigma=(Nb*chord(k))/(2.0*Pi*R_bp(k)*Rwt)           
           !Cl y Cd interpolación       
           call readtable(Name,Dir,alpha,L,D)   
!            write(*,*)L,D,name
           if (a .GT. 0.4) then
              ! if (Ct.LT.0.96) then
              CN=(8.0/9.0) + (4.0*f - (40.0/9.0))*a + ((50.0/9.0) - 4.0*f)*a**2.0
              a1=(18.0*f - 20.0 - 3.0*sqrt(CN*(50.0 - 36.0*f) + 12.0*f*(3.0*f - 4.0)))/(36.0*f - 50.0)
              ap1=-(1.0/2.0) + (1.0/2.0)*sqrt(1.0 + (4.0/(lamda*R_bp(k))**2)*a1*(1.0 - a1))                           
           else
              a1= 1.0/(1.0+((4.0*f*(sin(phi*pi/180.0))**2.0)/(sigma*L*cos(phi*Pi/180))))
              ap1=-(1.0/2.0)+(1.0/2.0)*sqrt(1.0+(4.0/(lamda*R_bp(k))**2)*a1*(1.0-a1))                         
           endif
           difa=(a1-a)
           difap=(ap1-ap)           
           a=a1          
           ap=ap1           
           m=m+1
        enddo
        dQ(k,j)=((lamda*R_bp(k))**3.0)*f*ap*(1.0-a)*(1.0-(D/L)*(1.0/tan(phi*Pi/180)))*lamda*dr(k)
!       dQ(k,j)=sigma*L*((1.0-a)**2)*(1.0/sin(phi*Pi/180))*(1.0-(D/L)*(1.0/tan(phi*Pi/180.0)))*((lamda*R_bp(k))**2.0)*lamda*dr(k)

     enddo
  enddo
  do j=1,nz
     Q(j)=sum(dQ(:,j))
  end do
 !CQ(:)=Q(:)/(0.5*1.225*(U(:)**2)*Pi*Rwt**3)
 !CP(:)=omgmax(:)*Q(:)/(0.5*1.225*(U(:)**3.)*Pi*Rwt**2.)
  CP(:)=(8.0/((Rwt*omgmax(:)*2.0*Pi)/(U(:)*60.0))**2.0)*Q(:)  
 !CP(:)=(2.0/((Rwt*omgmax(:)*2.0*Pi)/(U(:)*60.0))**2.0)*Q(:)  
 !write(*,'(1F7.4)')Q,(Rwt*omgmax(:)*2.0*Pi)/(U(:)*60.0)
!  write(*,'(<NE>F10.3)')CP
 !write(*,'(<NE>F10.5)')
  open(1,file='dQ.txt',status='replace')
write(1,'(A)')"Los resultados fueron obtenidos de acuerdo a Manwell JF et al Wind energy explained Second Edition "
  write(1,'(5X,a6,3X,a3,6X,a3,7X,a5,a11,4X,a10)')"WSpeed","TSR","Omg","Pitch","Cp","Aero_Power"
  do j=1,nz
     write(1,'(6F10.3)') U(j),(Rwt*omgmax(j)*2.0*Pi)/(U(j)*60.0),omgmax(j),pitch(j),CP(j),CP(j)*0.5*1.225*Pi*Rwt**2*U(j)**3/1000
  end do
  close(1)
endsubroutine Bem_method

!!!!*******************************************************************************************************************************
!!!!*******************************************************************************************************************************
function axial1(a,sigma,phi,L,D,Pi,f)
    integer::axial
  real::a,sigma,phi,L,D,Pi,f,fa
  axial1=sigma*(L*cos(phi*Pi/180.0)+D*sin(phi*Pi/180))*(1.0-a)-8.0*f*a*(sin(phi*Pi/180.0))**2.0
end function axial1
!!!!*******************************************************************************************************************************
function tangen1(ap,sigma,phi,L,D,Pi,f)
  real::ap,sigma,phi,L,D,Pi,f,fap
    integer::tangen
 tangen1  =sigma*(L*sin(phi*Pi/180.0)-D*cos(phi*Pi/180.0))*(1.0+ap)-8.0*f*ap*sin(phi*Pi/180.0)*cos(phi*Pi/180.0)
end function tangen1
!!!!*******************************************************************************************************************************
subroutine readtable(Name,Dir,x,L,D)
  implicit none
!  use bspline2
  integer:: n,i,col,maxcols
  integer::error
  real::x,L,D
  real,allocatable::a(:,:),alpha(:),Lift(:),Drag(:)
  CHARACTER(len=50)::Name
  CHARACTER(len=100)::Dir
  !Name='AspaDos0160.dat'
  Name=ADJUSTL(Name(1:LEN_TRIM(Name)))
  !Dir='/home/carlos/Documentos/BEM_Method/'
  Dir=ADJUSTL(Dir(1:LEN_TRIM(Dir)))

  call Ndimension(Name,Dir,n)
!  write(n)
  allocate(a(n,4),alpha(n-1),Lift(n-1),Drag(n-1))
  open(unit=10,file=Dir(1:LEN_TRIM(Dir))//Name(1:LEN_TRIM(Name)),iostat=error)
   do i=1,12
      read(10,*)
  end do
  ! Then start reading data
  a=0.
  maxcols=4
  do i=1,n-1
     read(10,*)a(i,:)
  enddo
  close(10)  
  do i=1,n-1
     alpha(i)=a(i,1)
     Lift(i)=a(i,2)
     Drag(i)=a(i,3)
 enddo

 call interpolacion(alpha,lift,n-2,x,L)
 call interpolacion(alpha,Drag,n-2,x,D) 
endsubroutine readtable
!!!!****************************************************************************************
subroutine Ndimension(Name,Dir,n)
  implicit none
  integer:: n,i,col,max_cols
  integer::error
  CHARACTER(len=50)::Name
  CHARACTER(len=100)::Dir
  open(unit=10,file=Dir(1:LEN_TRIM(Dir))//Name(1:LEN_TRIM(Name)),iostat=error)
!   write(*,*)Dir(1:LEN_TRIM(Dir))//Name(1:LEN_TRIM(Name))
  do i=1,12
      read(10,*)
  end do
  ! Then start reading data
  n=1
  max_cols=4
    do
     read(10,*, iostat = error) 
     !write(*,*)error
     IF (error > 0) then
      !  write(*,*)'check input. Something was wrong',error
        EXIT
     else if (error < 0) then
        !write(*,*)'the total is', n
        exit
     else
        n = n + 1
     endif             
  enddo
  close(10)
endsubroutine Ndimension
!***************************************************************************************************************************************
SUBROUTINE interpolacion(t,y,n,x,z)
      integer n
      real, dimension (0:n):: t,y
      real, dimension (1:10)::z
      real, dimension (0:n+1) ::a
      real, dimension (0:n+1):: h 

      interface
      subroutine bspline2_coef(n,t,y,a,h)              
      integer, intent(in) :: n                         
      real, dimension (0:n), intent(in):: t, y         
      real, dimension (0:n+1), intent(out):: a         
      real, dimension (0:n+1), intent(out):: h          
      end subroutine bspline2_coef
      function bspline2_eval(n,t,a,h,x)                 
      real, dimension (0:n+1), intent(in) :: a          
      real, dimension (0:n+1), intent(in) :: h          
      real, dimension (0:n), intent(in) :: t            
      integer, intent(in) :: n                          
      real, intent(in) :: x                             
      end function bspline2_eval
      end interface
        
      call bspline2_coef(n,t,y,a,h) 
      z(1) = bspline2_eval(n,t,a,h,x)      
     ! print*, z(1)
    end subroutine interpolacion
!!*******************************************************************************************************************
    subroutine bspline2_coef(n,t,y,a,h)              
      integer, intent(in) :: n                         
      real, dimension (0:n), intent(in):: t, y         
      real, dimension (0:n+1), intent(out):: a         
      real, dimension (0:n+1), intent(out):: h          
      integer :: i
      real :: delta, gamma, p,q
      do i = 1,n
      h(i) = t(i) - t(i-1)
      end do
      h(0) = h(1) 
      h(n+1) = h(n) 
      delta = -1.0
      gamma = 2.0*y(0)
      p = delta*gamma
      q = 2.0
      do  i = 1,n      
         r  = h(i+1)/h(i)
         delta = -r*delta
         gamma = -r*gamma + (r + 1.0)*y(i)
         p = p + gamma*delta 
         q = q + delta*delta 
      end do
      a(0) = -p/q 
      do i = 1,n+1
         a(i) = ((h(i-1)+h(i))*y(i-1)-h(i)*a(i-1))/h(i-1)
      end do
end subroutine bspline2_coef
!!!*************************************************************************************************************************  
function bspline2_eval(n,t,a,h,x)                 
      real, dimension (0:n+1), intent(in) :: a          
      real, dimension (0:n+1), intent(in) :: h          
      real, dimension (0:n), intent(in) :: t            
      integer, intent(in) :: n                          
      real, intent(in) :: x                                               
      integer :: i
      do i = n-1,0,-1      
         if ( x - t(i) >= 0.0) exit
      end do
      i = i + 1
     ! print *,"i,x",i,x
       d  =(a(i+1)*(x - t(i-1)) + a(i)*(t(i) - x + h(i+1)))/(h(i) + h(i+1)) 
 e = (a(i)*(x - t(i-1) + h(i-1)) + a(i-1)*(t(i-1) - x + h(i)))/(h(i-1) + h(i))
      bspline2_eval = ((d*(x - t(i-1)) + e*(t(i) - x)))/h(i)      
    end function bspline2_eval
!!!*************************************************************************************************************************  

                                                                       



