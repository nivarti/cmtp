set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels 
   set xtic auto                          # set xtics automatically 
   set ytic auto                          # set ytics automatically
   set title "Convergence Behaviour of Gauss-Seidel Scheme"
   set xlabel "Iterations"
   set ylabel "Maximum Change in Solution"
#   set key 0.01,100
#   set label "Yield Point" at 0.003,260	
#      set arrow from 0.0028,250 to 0.003,280
#      set xr [0.0:0.022]
#      set yr [0:325]
   set log y
   set term post "Helvetica" 18
   set output 'sor.eps'
      plot    "SOR" using 1:2 title 'No Relaxation' with lines , \
	 "SOR" using 1:3 title 'Relaxation Coefficient: 1.3' with lines,\
	    "SOR" using 1:4 title 'Relaxation Coefficient: 1.5' with lines
