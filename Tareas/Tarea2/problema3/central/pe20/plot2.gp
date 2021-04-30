 set xrange[0:10.0];set yrange[0:1.0]; set view map;
 set contour base; unset key; unset surface
 set ylabel 'yc'
set xlabel 'xc'
set cntrparam level incremental 0,0.1,10
 sp 'Temp2D3700.txt' u 1:2:3 w l lt -1 lw 1.5


