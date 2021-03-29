set palette @MATLAB
unset key
set tics font ", 10"
set xlabel "Velocidad,m/s"
set ylabel "Potencia de salida, kW"
set size square
p[0.0:25.0][0.0:13.0] 'dQ.txt' u 1:2 w l,'dQ.txt' u 3:4 w l
