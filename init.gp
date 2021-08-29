#Grafico de las posiciones iniciales
set xlabel ''
set ylabel ''

set xzeroaxis

set title ''

#set key top right

set grid xtic
set grid ytic

set size square

#set xrange[0:1]
#set yrange[0:1]


plot 'init_li0_03.dat' t 'li0' w p pt 7 ps 0.5,'init_ion_03.dat' t 'ion' w p pt 7 ps 0.5

set terminal pdf color
set output 'grafico_init.pdf'
replot 
exit
