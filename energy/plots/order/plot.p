set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels 
   set xtic auto                          # set xtics automatically 
   set ytic auto                          # set ytics automatically
#   set title "Error Norms with varying Mesh Size"
   set xlabel "Cell Size, dx"
   set ylabel "Error Norm"
#   set key 0.01,100
   set label "Slope    = 1.9906" at 0.02, 0.005
   set log xy
   set grid
   set term post "Helvetica" 18
   set output 'eeOrder.eps'
      plot "error_ee" using 1:2 title 'Explicit Euler' with linespoints pointtype 4
#       plot "error_ie" using 1:2 title 'Implicit Euler' with linespoints pointtype 12
