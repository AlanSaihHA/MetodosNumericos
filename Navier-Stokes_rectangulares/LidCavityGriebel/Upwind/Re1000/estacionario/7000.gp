 set size square; set xrange[0:1]; set yrange[0:1]; unset key
 p 'vel7000.txt' u 1:2:(0.25*$3):(0.25*$4) w vec
