set size square; set xrange[0:1];set yrange[0:2.5]; set ylabel 'Temperatura [°C]'; set xlabel 'xc'; set key left
plot 'temp-quick.txt' u 1:2 w p t 'Sol. Numerica', 'temp.txt' u 1:3 w l t 'Sol. Analitica'
