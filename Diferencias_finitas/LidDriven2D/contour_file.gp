reset
set size ratio 1
set contour base
set cntrparam bspline order 10 levels incr -0.5,0.01,0.5
unset surface
set table 'streamlines_data.dat'
sp 'vec_field.dat' u 1:2:6 
unset table
