set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels 
   set xtic auto                          # set xtics automatically 
   set ytic auto                          # set ytics automatically
   set xlabel "y direction"
   set ylabel "Temperature Gradient"
#   set key 0.01,100
#   set label "Slope    = 1.9906" at 0.02, 0.005
#   set log xy
   set grid
   set term post "Helvetica" 18
   set output 'gradientx.eps'
   plot "grad" using 1:2 title 'dT/dy(y)' with lines
#       plot "error_ie" using 1:2 title 'Implicit Euler' with linespoints pointtype 12
