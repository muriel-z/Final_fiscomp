#Grafico de las posiciones iniciales
set xlabel ''
set ylabel ''

set xzeroaxis

set title ''

#set key top right

set grid xtic
set grid ytic

set size square

plot 'init_li0.dat' t 'li0' w p pt 7 ps 0.5,'init_ion.dat' t 'ion' w p pt 7 ps 0.5

set terminal pdf color
set output 'grafico_init.pdf'
replot 
exit
