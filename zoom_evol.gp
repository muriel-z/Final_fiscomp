#Grafico de las posiciones evolucionadas
set xlabel ''
set ylabel ''

set xzeroaxis

set title 'Zoom en la formacion de dendritas con E=0, dt=1ms y t=500ms'

set key box top right

set grid xtic
set grid ytic

set xrange[0:0.3]
set yrange[0:0.3]

set size square

plot 'evol_li0.dat' t 'li0' w p pt 7 ps 1,'evol_lim.dat' t 'ion' w p pt 7 ps 1

set terminal pdf color
set output 'grafico_zoom_evol.pdf'
replot 
exit
