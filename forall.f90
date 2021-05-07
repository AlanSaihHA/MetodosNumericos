Program forall

Integer N,i,j
REAL, allocatable :: A(:,:), B(:,:)

N=5

allocate (A(N,N), B(N,N))

do i=1,5
do j=1,5
B(i,j)=j
enddo
enddo

do i=1,n
write(*,*) B(i,1), B(i,2), B(i,3), B(i,4), B(i,5)
enddo

write(*,*) ''


do i=1,n
write(*,*) A(i,1), A(i,2), A(i,3), A(i,4), A(i,5)
enddo

write(*,*) ''

FORALL (I=1:N) A(I,I) = 1


do i=1,n
write(*,*) A(i,1), A(i,2), A(i,3), A(i,4), A(i,5)
enddo

write(*,*) ''

forall (I=1:N, J=1:N, A(I,J) .NE. 0.0 ) B(I,J) = 1.0/A(I,J)

do i=1,n
write(*,*) B(i,1), B(i,2), B(i,3), B(i,4), B(i,5)
enddo


end program forall
