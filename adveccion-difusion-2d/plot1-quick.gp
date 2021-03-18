set xrange[0:1];set yrange[0:1];set view map;set size square
unset key
p 'velc-quick.txt' u 1:2:(0.1*$3):(0.1*$4) w vec

