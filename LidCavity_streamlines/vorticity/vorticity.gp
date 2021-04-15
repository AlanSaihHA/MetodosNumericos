set palette @MATLAB
unset surface; unset key; set contour base; set cntrparam level 160;unset colorbox;set view map
set tics font ", 10"
set xlabel "X=1.0"
set ylabel "$ \\omega $ ,Re=100.0"
set size square
sp[0.0:1.0][0.0:1.0] 'omega2000.txt' w pm3d
