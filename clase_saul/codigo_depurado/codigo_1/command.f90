program comando
implicit none
character*20 fname,command, xname
integer lfnm

write(*,*)' TYPE name of output files'    !Escritura de directorio donde se guardaran los datos
read (5,'(a10)')fname
lfnm=len_trim(fname)+1					
         
command='mkdir '
write(*,*) command
command(7:7+lfnm)=fname
write(*,*) command
call system(command)

xname=fname
fname(lfnm:lfnm+1)='/'
write(*,*) fname

endprogram comando
