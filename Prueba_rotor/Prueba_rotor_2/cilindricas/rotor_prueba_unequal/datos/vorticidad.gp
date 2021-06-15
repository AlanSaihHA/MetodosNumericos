#set palette @MATLAB
unset surface; unset key; set contour base; set cntrparam level 50;unset colorbox;set view map
#set tics font ", 10"
#set xlabel "X=1.0"
set cntrlabel onecolor
set ylabel "$ \\psi $ ,Re=1000.0"
set xlabel 'x'
set xrange [0:7]
set yrange [0:2]
sp 'omega1000.txt' w pm3d lc 2
