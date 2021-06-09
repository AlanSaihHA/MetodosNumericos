set palette @MATLAB
unset key
set tics font ", 10"
set xlabel "X=1.0"
set ylabel "Velocidad,Re=100.0"
set size square
p[0.0:10.0][0.0:1.0] 'vel2000.txt' u 1:2:(0.15*$3):(0.15*$4)  w vec
