#Grafico del observable rho
set xlabel ''
set ylabel ''

set xzeroaxis

set title ''

set key box top right

set grid xtic
set grid ytic

plot 'observables.dat' u 1:3 t 'rho(t)' w p pt 7 ps 0.4

set terminal pdf color
set output 'grafico_rho.pdf'
replot 
exit
