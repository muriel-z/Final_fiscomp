program prueba_tita
use precision
use mzranmod
implicit none
real(np),parameter      ::pi=acos(-1._np)
real(np)                ::tita,gx,gy
integer                 ::i

tita=0._np
gx=0._np
gy=0._np

open(21,file='tita.dat')
do i=1,3
    tita = 2._np*pi*rmzran()
    gx = cos(tita)
    gy = sin(tita)
    write(21,*)tita,gx,gy
    write(*,*)tita,gx,gy
enddo
close(21)

end program prueba_tita
