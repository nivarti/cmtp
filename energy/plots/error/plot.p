set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels 
   set xtic auto                          # set xtics automatically 
   set ytic auto                          # set ytics automatically
#   set title "Error Profile for Solution"
   set xlabel "x"
   set ylabel "y"
#   set key 0.01,100	
   set term post "Helvetica" 18
   set output 'error.eps'
   set hidden3d
   splot    "T_error_ee" title 'Explicit Euler Error(x,y)' with lines, \
      "T_error_ie" title 'Implicit Euler Error(x,y)' with lines
