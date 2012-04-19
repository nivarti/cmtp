set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels 
   set xtic auto                          # set xtics automatically 
   set ytic auto                          # set ytics automatically
#   set title "Fully developed Velocity Profile"
   set xlabel "y"
   set ylabel "Velocity"
   set grid
#   set key 0.01,100	
   set term post "Helvetica" 18
   set output 'xvel.eps'
#   set hidden3d
      plot    "xvel" title 'u(y)' with lines
