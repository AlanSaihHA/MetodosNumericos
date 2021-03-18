set xrange[0:1];set yrange[0:1];set view map;set size square
set contour base
unset surface
unset key
set cntrparam level discrete  0.0,0.0000000000001,0.000001,0.00000001,0.00000000000000001,0.0001,0.02,0.03,0.04,0.001,0.1,0.05,0.2,0.21,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1
 sp 'Temp2D-upwind0.txt' u 1:2:3 w l lt -1 lw 1.5

