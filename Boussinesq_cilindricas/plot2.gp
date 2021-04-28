 set xrange[-1:1.0];set yrange[-1:1.0]; set view map;
 set contour base; unset key; unset surface
set cntrparam level incremental 0,0.1,10
 sp 'T3000.txt' u 1:2:3 w l lt -1 lw 1.5


