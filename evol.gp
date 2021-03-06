#Grafico de las posiciones evolucionadas
set xlabel ''
set ylabel ''

set xzeroaxis

set title ''

set key box top right

set grid xtic
set grid ytic

#set xrange[0:154]
#set yrange[0:154]

set size square

plot 'evol_li0_04.dat' t 'li0' w p pt 7 ps 0.4,'evol_lim_04.dat' t 'ion' w p pt 7 ps 0.5

set terminal pdf color
set output 'grafico_evol.pdf'
replot 
exit
