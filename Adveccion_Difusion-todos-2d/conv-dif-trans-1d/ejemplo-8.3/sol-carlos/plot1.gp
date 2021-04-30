set size square; set xrange[0:1.5];set yrange[0:14]; set ylabel 'phi'; set xlabel 'xc'; set key left
set grid ytics mytics
plot 'Temp2D1000.txt' u 1:2 w p
