#Grafico de las posiciones evolucionadas
set xlabel ''
set ylabel ''

set xzeroaxis

set title ''

set key box top right

set grid xtic
set grid ytic

set size square

plot 'evol_li0.dat' t 'li0' w p pt 7 ps 0.5,'evol_lim.dat' t 'ion' w p pt 7 ps 0.5

set terminal pdf color
set output 'grafico_evol.pdf'
replot 
exit
