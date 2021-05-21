 set size square; set xrange[0:10]; set yrange[0:4]; unset key
 p 'vel1600.txt' u 1:2:(0.2*$3):(0.2*$4) w vec
