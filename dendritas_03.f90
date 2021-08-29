program dendritas_03
use precision
use mzranmod
implicit none

!Definicion de variables
real(np)                                ::ti,tf,t,dt !Variables temporales
real(np)                                ::D,q,Long !coeficientes
real(np)                                ::E0_x,E0_y,u,rx,ry
real(np)                                ::rli0,rlim,dx,dy,Lx,Ly,x,y,datt,distx,disty,dist
integer,dimension(:),allocatable       ::ex,ey,e0x,e0y !vectores espaciales de los iones
integer                                ::exs,eys !escalares
integer,dimension(:),allocatable       ::li_xd,li_yd !Arreglo del litio depositado sobre el anodo
integer,dimension(:),allocatable       ::li0x,li0y !li0
integer,dimension(:),allocatable       ::li_aux_x, li_aux_y
integer                                 ::gx,gy !vector unitario aleatorio.
integer                                 ::n0,n0max,nm !numeros de particulas 
integer                                 ::nn !numero de voxels
integer                                 ::nt !numeros de pasos temporales 
integer                                 ::i,j,k,m,l
integer                                 ::nq,nrx,nry
integer                                 ::ndistx,ndisty,ndist,ndatt

!Inicializacion de variables
dt = 1.0e-6_np
n0 = 77 !comienzo con este numero de li0 depositados sobre la sup
n0max = 600 !numero maximo de particulas li0
rli0 = 1.67e-10_np  !radio atomico del li0
rlim = 1.2e-10_np   !radio atomico del li+
nm = 50 !Este numero debe mantenerse a lo largo de la ev temporal
D = 1.4e-10_np !coef de dif del Li+ en el electrolito
gx = 0
gy = 0
x=0._np
y=0._np
Long = 16.7e-9_np
datt = 1.3_np*rli0/Long !que es rli0/(pi/4)
q = sqrt(2._np*D*dt)/Long !desplazamiento medio debido a la difusion

u = 5.6e-9_np
E0_x = 0._np
E0_y = -1.7e5_np
rx = u*E0_x*dt/Long
ry = u*E0_y*dt/Long !mu+*E*dt/Long es el desplazamiento debido al campo electrico

nq = 154 !q*154/Long
nrx = 0 !rx*154/Long
nry = -9 !ry*154/Long

nn = 154

ndistx = 0
ndisty = 0
ndist = 0
ndatt = 2

dx = 0.65_np*rli0/Long !1.3*rli0/2
dy = 0.65_np*rli0/Long

allocate(ex(1:nm),ey(1:nm),e0x(1:nm),e0y(1:nm)) !son las variales espaciales que guardan informacion de la posciion de todas las particulas a un dado tiempo t !OJO....si se hace esto el numero nm se debe mantener constante en cada paso temporal. Entonces se debe mantener constante el numero de iones
allocate(li_xd(1:n0),li_yd(1:n0))
allocate(li_aux_x(1:600),li_aux_y(1:600))

!creacion de los Li0 depositados uniformemente sobre el anodo--------------------------------------------------------------
open(21,file='init_li0_03.dat',status='replace')
do i = 1,n0
    !vector que guarda la posicion del litio ordenado a lo largo de x normalizado por Long
    li_xd(i) = 2*real(i,np)
    li_yd(i) = 0
    write(21,*)li_xd(i),li_yd(i)
enddo

!Creacion del electrolito con iones*******************************

!Creacion de la posicion inicial de los 50 iones-------------------
open(22,file='init_ion_03.dat',status='replace')
do i = 1,nm
!Posiciones iniciales normalizadas por la longitud Long
    e0x(i) = anint(rmzran()*154._np)
    e0y(i) = anint(rmzran()*154._np)
    write(22,*)e0x(i),e0y(i)
enddo


li_aux_x(:) = 0
li_aux_y(:) = 0
ex(:) = 0
ey(:) = 0

m = n0 !numero inicial de li0
nt = 190
open(23,file='evol_li0_03.dat',status='replace')
open(24,file='evol_lim_03.dat',status='replace')
!Definicion de las ecuaciones de movimiento browniano
!Evoluion temporal
do j = 1,nt
    do i = 1,nm
        !Definicion del vector unitario aleatorio g
        gx = anint(rmzran()*3._np)-1
        gy = anint(rmzran()*3._np)-1
        !Evolucion de la posicion de los iones
        ex(i) = e0x(i) + nq*gx + nrx
        ey(i) = e0y(i) + nq*gy + nry
        !Meto las PBC en una caja de tamaño nnxnn
        ex(i) = ex(i) - nn*anint(real(ex(i),np)/nn-0.5_np)
        ey(i) = ey(i) - nn*anint(real(ey(i),np)/nn+0.5_np) + nn
        !Definicion de la condicion Li+-->Li0
        do k = 1,m
            if ( m.le.n0) then
                ndistx = ex(i)-li_xd(k) 
                ndisty = ey(i)-li_yd(k)
                ndist = anint(sqrt( real(ndistx*ndistx + ndisty*ndisty ,np)))
            elseif ( m.gt.n0) then
                ndistx = ex(i)-li0x(k) 
                ndisty = ey(i)-li0y(k)
                ndist = anint(sqrt( real(ndistx*ndistx + ndisty*ndisty ,np)))
            endif
            if (ndist.le.ndatt) then 
                m = m+1 !Aumenta la dimension del arreglo del li0
                exs = ex(i)
                eys = ey(i)
                call save_li0(m,li_xd,li_yd,exs,eys,li0x,li0y,li_aux_x,li_aux_y) !Guardo la nueva posicion del li0
                !Tengo que reponer un ion en el espacio de forma aleatoria
                ex(i) = anint(rmzran()*154)
                ey(i) = anint(rmzran()*154)
                !Meto las PBC
                ex(i) = ex(i) - nn*anint(real(ex(i),np)/nn-0.5_np)
                ey(i) = ey(i) - nn*anint(real(ey(i),np)/nn+0.5_np) + nn
            else
                !Necesito poner un else porque de todas formas tengo que allocatear el li0x y li0y mas alla de si se le pega un ion o no !se supone que la linea siguiente llama a la subrutina pero conserva el numero m del litio depositado entonces no deberia cambiar
                call save_li0(m,li_xd,li_yd,exs,eys,li0x,li0y,li_aux_x,li_aux_y)
            endif
        enddo
    enddo
    !Renombre de posiciones para t+dt
    t = real(j,np)*dt
    e0x = ex
    e0y = ey
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
integer,intent(in)                             ::mm
integer,intent(in)                             ::exx,eyy
integer,dimension(1:77),intent(in)             ::li_xxd, li_yyd
integer,dimension(:),allocatable,intent(out)   ::li0xx, li0yy
integer,dimension(1:600),intent(inout)         ::li_auxx, li_auxy
integer                                        ::i,nn
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
nn=154
do i = 1,mm
    li0xx(i) = li0xx(i) - 154*anint(real(li0xx(i),np)/154._np-0.5_np)
    li0yy(i) = li0yy(i) - 154*anint(real(li0yy(i),np)/154._np+0.5_np)+154
enddo
    
end Subroutine save_li0


end program dendritas_03
