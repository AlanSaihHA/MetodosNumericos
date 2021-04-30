 set size square; set xrange[0:1]; set yrange[0:1]; unset key
 p 'Sol_analitica.txt' u 1:2:(0.5*$3):(0.5*$4) w vec
 set title 'Sol. Analitica'
