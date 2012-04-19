set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels 
   set xtic auto                          # set xtics automatically 
   set ytic auto                          # set ytics automatically
#   set title "Error Norms for Flux and Source"
   set xlabel "Cell Size, dx"
   set ylabel "Error Norm"
#   set key 0.01,100
   set label "Slope    = 1.9955" at 0.012,0.04
      set label "Slope = 1.9941" at 0.05,0.00005
#      set arrow from 0.0028,250 to 0.003,280
#      set xr [0.0:0.022]
#      set yr [0:325]
   set log xy
   set grid
   set term post "Helvetica" 18
   set output 'order.eps'
      plot    "Ef" using 1:2 title 'Flux Calculations' with linespoints pointtype 4, \
	 "Es" using 1:2 title 'Source Calculations' with linespoints pointtype 12
#	    sin(2*pi*x) title 'Exact Solution' with lines
