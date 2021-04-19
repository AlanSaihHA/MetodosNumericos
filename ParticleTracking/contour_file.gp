reset
set size ratio 1

set contour base
set cntrparam bspline order 10 levels incr -5.0,0.08,5.0
unset surface
set table "contours.dat"
sp 'VecField.dat' u 1:2:6
unset table


