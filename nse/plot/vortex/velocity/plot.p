#set size square
#unset title
#set pm3d map
#set format x""
#set format y""
#set cbrange[-0.1:1]
#set palette rgbformulae 22,13,10
#splot './u'


set view map
set samples 25, 25
set isosamples 26, 26
unset surface
set contour base
set cntrparam bspline
set cntrparam levels auto 10
set style data lines
set title "3D gnuplot demo - 2D contour projection of last plot" 
set xlabel "X axis" 
set xrange [ 0.00000 : 15.0000 ] noreverse nowriteback
set ylabel "Y axis" 
set yrange [ 0.00000 : 15.0000 ] noreverse nowriteback
set zlabel "Z axis" 
set zlabel  offset character 1, 0, 0 font "" textcolor lt -1 norotate
set zrange [ -1.20000 : 1.20000 ] noreverse nowriteback
splot "psi" using 1
