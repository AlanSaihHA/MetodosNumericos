Program Gauss_seidel
  implicit none
  integer max_iter,n, i,j,k
  real tolerance, error, x_aux
  real  phi_v(3), sp_v(3), A(3,3), phi_vold(3)  
  a = 0.0
  a(1,1)=300;a(2,2)=200;a(3,3)=200;a(4,4)=200;a(5,5)=300;
  a(1,2)=-100;a(2,3)=-100;a(3,4)=-100;a(4,5)=-100
  a(2,1)=-100;a(3,2)=-100;a(4,3)=-100;a(5,4)=-100
  
  sp_v(1)=20000;sp_v(2)=0;sp_v(3)=0;sp_v(4)=0;sp_v(5)=100000;
  !  x=0.0
  max_iter=1000
  tolerance=1e-5
  n=3
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
     write(*,*) error,k	
     k=k+1
  enddo
!  write(*,*)phi_vold
end Program Gauss_seidel

