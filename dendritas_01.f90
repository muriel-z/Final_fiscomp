program dendritas_01
use precision
use mzranmod
implicit none

!Definicion de variables
real(np)                                ::ti,tf,t,dt !Variables temporales
real(np)                                ::p0,pm !particulas
real(np)                                ::D,q,r !coeficientes
real(np)                                ::rli0,rlim,dx,dy,Lx,Ly,x,y,datt,distx,disty,dist
real(np),dimension(:),allocatable       ::ex,ey,e0x,e0y !vectores espaciales de los iones
real(np)                                ::exs,eys !escalares
real(np),dimension(:),allocatable       ::li_xd,li_yd !Arreglo del litio depositado sobre el anodo
real(np),dimension(:),allocatable       ::li0x,li0y !Arreglo del litio depositado sobre el anodo
real(np)                                ::gx,gy !vector unitario aleatorio.
integer                                 ::n0,n0max,nm !numeros de particulas 
integer                                 ::nt !numeros de pasos temporales 
integer                                 ::i,j,k,m

!Inicializacion de variables
dt = 1.0e-6_np
n0 = 100 !comienzo con este numero de li0
n0max = 600 !numero maximo de particulas li0
rli0 = 1.67e-10_np  !radio atomico del li0
rlim = 1.2e-10_np   !radio atomico del li+
nm = 50 !Este numero debe mantenerse a lo largo de la ev temporal
D = 1.4e-10_np !coef de dif del Li+ en el electrolito
q = sqrt(2._np*D*dt) !desplazamiento medio debido a la difusion
r = 0._np !mu+*E*dt es el desplazamiento debido al campo electrico
gx = 0._np
gy = 0._np
datt = 1.3_np*rli0 !que es rli0/(pi/4)
x=0._np
y=0._np

allocate(ex(1:nm),ey(1:nm),e0x(1:nm),e0y(1:nm)) !son las variales espaciales que guardan informacion de la posciion de todas las particulas a un dado tiempo t !OJO....si se hace esto el numero nm se debe mantener constante en cada paso temporal. Entonces se debe mantener constante el numero de iones
allocate(li_xd(1:100),li_yd(1:100))

!creacion de los Li0 depositados uniformemente sobre el anodo--------------------------------------------------------------
open(21,file='init_li0.dat',status='replace')
do i = 1,100
    !vector que guarda la posicion del litio ordenado a lo largo de x
    li_xd(i) = rli0*real(i,np)
    li_yd(i) = 0._np 
    write(21,*)li_xd(i),li_yd(i)
enddo

!Creacion del electrolito con iones*******************************
dx = rli0/2._np
dy = rli0/2._np
Lx = 16.7e-9_np
Ly = 16.7e-9_np
!Creacion de la posicion inicial de los 50 iones-------------------
open(22,file='init_ion.dat',status='replace')
do i = 1,nm
    e0x(i) = anint(rmzran()*200._np)*dx
    e0y(i) = anint(rmzran()*200._np)*dy
    write(22,*)e0x(i),e0y(i)
enddo

m = n0 !numero inicial de li0
nt = 10
open(23,file='evol_li0.dat',status='replace')
open(24,file='evol_lim.dat',status='replace')
!Definicion de las ecuaciones de movimiento browniano
!Evoluion temporal
do j = 1,nt
    do i = 1,nm
        !Definicion del vector unitario aleatorio g
        gx = anint(rmzran()*3._np)-1
        gy = anint(rmzran()*3._np)-1
    
        ex(i) = e0x(i) + q*gx + r
        ey(i) = e0y(i) + q*gy + r
        !Definicion de la condicion Li+-->Li0
        !En cada t+dt tengo que actualizar una lista con las posiciones de los li0 y en base a eso tambien pedir (1) la actualizacion de particulas ion, es decir que se mantenga cte su densidad cada vez que pierden uno y (2) que si el ion se acerca a una cierta distancia datt se vuelva li0
        do k = 1,m
            distx = ex(i)-li_xd(k) 
            disty = ey(i)-li_yd(k)
            dist = sqrt( distx*distx + disty*disty )
            if (dist<datt) then 
                !Aparece un nuevo punto en li_xd y li_yd que va a ser igual que las coordenadas que el ion viejo
                m = m+1
                exs = ex(i)
                eys = ey(i)
                call save_li0(m,li_xd,li_yd,exs,eys,li0x,li0y) !Guardo la nueva posicion del li0
                !Tengo que reponer un ion en el espacio en la parte superior por eso le doy en los 50 primeros lugares
                ex(i) = anint(rmzran()*50._np)*dx
                ey(i) = anint(rmzran()*50._np)*dy
            endif
            write(23,*)li0x(k),li0y(k)
        enddo
        write(24,*)ex(i),ey(i)
    enddo
    !Renombre de posiciones para t+dt
    t = real(j,np)*dt
    e0x = ex
    e0y = ey
enddo
    
deallocate(ex,ey)
close(21)
close(22)
close(23)
close(24)

contains

Subroutine save_li0(mm,li_dxx,li_dyy,exx,eyy,li0xx,li0yy)
integer,intent(in)                              ::mm
real(np),intent(in)                             ::exx,eyy
real(np),dimension(:),allocatable,intent(in)    ::li_dxx, li_dyy
real(np),dimension(:),allocatable,intent(out)   ::li0xx, li0yy
real(np),dimension(:),allocatable               ::li_auxx, li_auxy
integer                                         ::i,nn
allocate(li0xx(1:mm),li0yy(1:mm))
allocate(li_auxx(101:600),li_auxy(101:600))

li_auxx = 0._np
li_auxy = 0._np

!En los primeros 100 lugares guardo el litio depositado sobre el anodo
do i = 1,100
    li_auxx(i) = li_dxx(i)
    li_auxy(i) = li_dyy(i)
enddo
!En los siguientes lugares tengo que guardar las nuevas posiciones de los li+ que pasan a ser li0 
!--> No se muy bien como hacer esto... necesito ir guardando las coordenadas x e y de la entrada 101,102,...,m. Cuando m es el mas grande deacuerdo a los pasos temporales. ya que m va subiendo con el tiempo... 
!me puedo crear un vector auxiliar li_auxx,li_auxy. El vector auxiliar tiene en total 600 entradas, las cuales son nulas y las va llenando a medida que m>100 y se va guardando justamente las posiciones de cada li+ que paso a li0. De esta forma tengo una especie de lista donde guardo estos valores. y en li0xx y li0yy me guardo las coordenadas no nulas :D

!A ver, de alguna forma aca dentro tengo que poder sistematizar guardar todas las entradas de los vectores espaciales x e y cada vez que li+-->li0

li_auxx(m) = exx !Esto llena la entrada m
li_auxy(m) = eyy 

do i = 1,mm
    li0xx(i) = li_auxx(i)
    li0yy(i) = li_auxy(i)
enddo
    
end Subroutine save_li0


end program dendritas_01
