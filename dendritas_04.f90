program dendritas_04
use precision
use mzranmod
implicit none

!Definicion de variables
real(np)                                ::ti,tf,t,dt !Variables temporales
real(np)                                ::p0,pm !particulas
real(np)                                ::D,q,Long !coeficientes
real(np)                                ::rli0,rlim,dx,dy,Lx,Ly,x,y,datt,distx,disty,dist
real(np),dimension(:),allocatable       ::ex,ey,e0x,e0y !vectores espaciales de los iones
real(np)                                ::exs,eys !escalares
real(np),dimension(:),allocatable       ::li_xd,li_yd !Arreglo del litio depositado sobre el anodo
real(np),dimension(:),allocatable       ::li0x,li0y !li0
real(np),dimension(:),allocatable       ::li_aux_x, li_aux_y
real(np)                                ::E0_x,E0_y,u,rx,ry
real(np),parameter                      ::pi = -acos(-1._np)
real(np)                                ::gx,gy,tita !vector unitario aleatorio.
real(np)                                ::alfa,mod_li0
real(np)                                ::Rli_mx,Rli_my,Rli_0x,Rli_0y,modd
integer                                 ::n0,n0max,nm !numeros de particulas 
integer                                 ::nt !numeros de pasos temporales 
integer                                 ::i,j,k,m,l
real(np)                                ::lambda,rho

!Inicializacion de variables
dt = 1.0e-6_np
n0 = 77 !comienzo con este numero de li0 depositados sobre la sup
n0max = 600 !numero maximo de particulas li0
rli0 = 1.67e-10_np  !radio atomico del li0
rlim = 1.2e-10_np   !radio atomico del li+
nm = 50 !Este numero debe mantenerse a lo largo de la ev temporal
D = 1.4e-14_np !coef de dif del Li+ en el electrolito
tita = 0._np
gx = 0._np
gy = 0._np
x=0._np
y=0._np
Long = 16.7e-9_np
datt = 1.3_np*rli0/Long !que es rli0/(pi/4)
q = sqrt(2._np*D*dt)/Long !desplazamiento medio debido a la difusion

alfa = 0._np
mod_li0 = 0._np

Rli_mx = 0._np
Rli_my = 0._np
Rli_0x = 0._np
Rli_0y = 0._np
modd = 0._np

lambda = 0._np
rho = 0._np

u = 5.6e-13_np
E0_x = 0._np
E0_y = -1.7e7_np

rx = u*E0_x*dt/Long
ry = u*E0_y*dt/Long !mu+*E*dt/Long es el desplazamiento debido al campo electrico

allocate(ex(1:nm),ey(1:nm),e0x(1:nm),e0y(1:nm)) !son las variales espaciales que guardan informacion de la posciion de todas las particulas a un dado tiempo t !OJO....si se hace esto el numero nm se debe mantener constante en cada paso temporal. Entonces se debe mantener constante el numero de iones
allocate(li_xd(1:n0),li_yd(1:n0))
allocate(li_aux_x(1:600),li_aux_y(1:600))

!creacion de los Li0 depositados uniformemente sobre el anodo--------------------------------------------------------------
open(21,file='init_li0_04.dat',status='replace')
do i = 1,n0
    !vector que guarda la posicion del litio ordenado a lo largo de x normalizado por Long
    li_xd(i) = 1.3_np*rli0*real(i,np)/Long
    li_yd(i) = 0._np/Long 
    write(21,*)li_xd(i),li_yd(i)
enddo

!Creacion del electrolito con iones*******************************
dx = rli0/(2._np*Long)
dy = rli0/(2._np*Long)
!Creacion de la posicion inicial de los 50 iones-------------------
open(22,file='init_ion_04.dat',status='replace')
do i = 1,nm
!Posiciones iniciales normalizadas por la longitud Long
    e0x(i) = dnint(rmzran()*200._np)*dx
    e0y(i) = dnint(rmzran()*200._np)*dy
    write(22,*)e0x(i),e0y(i)
enddo


li_aux_x(:) = 0._np
li_aux_y(:) = 0._np

m = n0 !numero inicial de li0
nt =5800
open(23,file='evol_li0_04.dat',status='replace')
open(24,file='evol_lim_04.dat',status='replace')
open(25,file='observables.dat',status='replace')
!Definicion de las ecuaciones de movimiento browniano
!Evoluion temporal
do j = 1,nt
    do i = 1,nm
        !Definicion del vector unitario aleatorio g
        tita = 2._np*pi*rmzran()
        gx = cos(tita)
        gy = sin(tita)
        !Evolucion de la posicion de los iones
        ex(i) = e0x(i) + q*gx + rx
        ey(i) = e0y(i) + q*gy + ry
        !Meto las PBC en una caja de tamaño 1x1
        ex(i) = ex(i) - dnint(ex(i)-0.5_np)
        ey(i) = ey(i) - dnint(ey(i)+0.5_np) + 1._np
        !Definicion de la condicion Li+-->Li0
        do k = 1,m
            if (m==600) stop
            if ( m.le.n0) then
                distx = ex(i)-li_xd(k) 
                disty = ey(i)-li_yd(k)
                dist = sqrt( distx*distx + disty*disty )
            elseif ( m.gt.n0) then
                distx = ex(i)-li0x(k) 
                disty = ey(i)-li0y(k)
                dist = sqrt( distx*distx + disty*disty )
            endif
            !alfa = atan(disty/distx)
            if (dist<datt) then 
                !Calculo el modulo del nuevo li0 para que este a una distancia datt del li0(k)
                if ( m.le.n0) then
                    modd = sqrt((ex(i)-li_xd(k))*(ex(i)-li_xd(k))+(ey(i)-li_yd(k))*(ey(i)-li_yd(k)))
                    exs = (ex(i)-li_xd(k))*datt/modd + li_xd(k)
                    eys = (ey(i)-li_yd(k))*datt/modd + li_yd(k)
                    !mod_li0 = li_xd(k)*cos(alfa) + li_yd(k)*sin(alfa) &
                             ! + 0.5_np*sqrt(4._np*li_xd(k)*li_xd(k)*(cos(alfa)*cos(alfa)+1._np) &
                             ! + 4._np*li_yd(k)*li_yd(k)*(sin(alfa)*sin(alfa)+1._np) &
                             ! + 2._np*li_xd(k)*li_yd(k)*cos(alfa)*sin(alfa) - 4._np*datt)
                elseif ( m.gt.n0) then
                    modd = sqrt((ex(i)-li0x(k))*(ex(i)-li0x(k))+(ey(i)-li0y(k))*(ey(i)-li0y(k)))
                    exs = (ex(i)-li0x(k))*datt/modd + li0x(k)
                    eys = (ey(i)-li0y(k))*datt/modd + li0y(k)
                    !mod_li0 = li0x(k)*cos(alfa) + li0y(k)*sin(alfa) &
                             ! + 0.5_np*sqrt(4._np*li0x(k)*li0x(k)*(cos(alfa)*cos(alfa)+1._np) &
                             ! + 4._np*li0y(k)*li0y(k)*(sin(alfa)*sin(alfa)+1._np) &
                             ! + 2._np*li0x(k)*li0y(k)*cos(alfa)*sin(alfa) - 4._np*datt)
                endif
                !write(*,*)modd
                !write(*,*) exs,eys
                m = m+1 !Aumenta la dimension del arreglo del li0
                !exs = mod_li0*cos(alfa)
                !eys = mod_li0*sin(alfa)
                call save_li0(m,li_xd,li_yd,exs,eys,li0x,li0y,li_aux_x,li_aux_y) !Guardo la nueva posicion del li0
                !Tengo que reponer un ion en el espacio de forma aleatoria
                ex(i) = anint(rmzran()*200._np)*dx
                ey(i) = anint(rmzran()*200._np)*dy
                !Meto las PBC
                ex(i) = ex(i) - dnint(ex(i)-0.5_np)
                ey(i) = ey(i) - dnint(ey(i)+0.5_np) + 1._np
            else
                !Necesito poner un else porque de todas formas tengo que allocatear el li0x y li0y mas alla de si se le pega un ion o no !se supone que la linea siguiente llama a la subrutina pero conserva el numero m del litio depositado entonces no deberia cambiar
                call save_li0(m,li_xd,li_yd,exs,eys,li0x,li0y,li_aux_x,li_aux_y)
            endif
        enddo
    enddo
    !Renombre de posiciones para t+dt
    t = real(j,np)*dt
    lambda = maxval(li0y)
    rho = m/dx
    e0x = ex
    e0y = ey
    write(25,*)t,lambda,rho
    write(*,*) real(j,np)/real(nt,np)*100._np
enddo
write(*,*)m
write(*,*)nm

do l = 1,nm
    write(24,*)ex(l),ey(l)
enddo

do l = 1,m
    write(23,*)li0x(l),li0y(l)
enddo
    
deallocate(ex,ey)
close(21)
close(22)
close(23)
close(24)

!*******************************************************************
contains

Subroutine save_li0(mm,li_xxd,li_yyd,exx,eyy,li0xx,li0yy,li_auxx,li_auxy)
integer,intent(in)                              ::mm
real(np),intent(in)                             ::exx,eyy
real(np),dimension(1:77),intent(in)             ::li_xxd, li_yyd
real(np),dimension(:),allocatable,intent(out)   ::li0xx, li0yy
real(np),dimension(1:600),intent(inout)         ::li_auxx, li_auxy
integer                                         ::i,nn
allocate(li0xx(1:mm),li0yy(1:mm))

!En los primeros 100 lugares guardo el litio depositado sobre el anodo
do i = 1,77
    li_auxx(i) = li_xxd(i)
    li_auxy(i) = li_yyd(i)
enddo
!En los siguientes lugares tengo que guardar las nuevas posiciones de los li+ que pasan a ser li0 

!me creo un vector auxiliar li_auxx,li_auxy. El vector auxiliar tiene en total 600 entradas, las cuales son nulas y las va llenando a medida que m>100 y se va guardando justamente las posiciones de cada li+ que paso a li0. De esta forma tengo una especie de lista donde guardo estos valores. y en li0xx y li0yy me guardo las coordenadas no nulas :D

!Se llena la entrada m
li_auxx(mm) = exx 
li_auxy(mm) = eyy 

do i = 1,mm
    li0xx(i) = li_auxx(i)
    li0yy(i) = li_auxy(i)
enddo
!Doy las PBC
do i = 1,mm
    li0xx(i) = li0xx(i) - dnint(li0xx(i)-0.5_np)
    li0yy(i) = li0yy(i) - dnint(li0yy(i)+0.5_np) + 1._np
enddo
    
end Subroutine save_li0


end program dendritas_04
