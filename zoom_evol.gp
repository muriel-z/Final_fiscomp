#Grafico de las posiciones evolucionadas
set xlabel ''
set ylabel ''

set xzeroaxis

set title 'Zoom en la formacion de dendritas dt=1us y t=240us'

set key box top right

set grid xtic
set grid ytic

set xrange[0:50]
set yrange[0:50]

set size square

plot 'evol_li0_03.dat' t 'li0' w p pt 7 ps 1.2,'evol_lim_03.dat' t 'ion' w p pt 7 ps 1

set terminal pdf color
set output 'grafico_zoom_evol.pdf'
replot 
exit
