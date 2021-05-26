program leer_escribir
implicit none
real a(1000,2)
integer i,b, REASON


open(20,file='NACA060006_int.txt')

i=1
DO 
   READ(20,*,IOSTAT=REASON)  a(i,1), a(i,2)
   IF (REASON>0) THEN
      print*, '...algo anda mal...'
        exit
   ELSE IF (REASON<0) THEN
      print*, '...fin del archivo...'
         exit
   ELSE
      print*, 'linea: ',i,'leída'
      i=i+1
   END IF
END DO

close(20)

b=i-1

do i=1, b
   print*, a(i,1), a(i,2)
enddo

print*,b
!print*, a(38,1), a(38,2)
end program leer_escribir

