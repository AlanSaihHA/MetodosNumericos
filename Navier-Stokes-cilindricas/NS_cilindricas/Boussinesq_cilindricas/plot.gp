set xrange[-1:1];set yrange[-1:1];set size square
set view map;
set contour base
set palette rgbformulae 33,13,10
unset surface
unset key
set clabel
set cntrparam bspline order 10 levels incremental 0,0.05,1 

sp 'T1000.txt' u 1:2:3 w l palette lw 1.5

#sp 'T1000.txt' w pm3d

