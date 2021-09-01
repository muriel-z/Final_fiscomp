#Grafico del observable lambda
set xlabel ''
set ylabel ''

set xzeroaxis

set title ''

set key box top right

set grid xtic
set grid ytic

plot 'observables.dat' u 1:2 t 'lambda(t)' w p pt 7 ps 0.4

set terminal pdf color
set output 'grafico_lambda.pdf'
replot 
exit
