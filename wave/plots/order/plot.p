set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels 
   set xtic auto                          # set xtics automatically 
   set ytic auto                          # set ytics automatically
   set title "Error Norms for Different Mesh Sizes"
   set xlabel "Cell Size, dx"
   set ylabel "Error Norm"
#   set key 0.01,100
   set label "Slope    = 1.99939233" at 0.015,0.01
      set label "Slope = 1.47052215" at 0.003,0.05
#      set arrow from 0.0028,250 to 0.003,280
#      set xr [0.0:0.022]
#      set yr [0:325]
   set log xy
   set grid
   set term post "Helvetica" 18
   set output 'Order.eps'
      plot    "N" using 1:2 title 'Regular Trend' with points pointtype 4, \
	 "N2" using 1:2 title 'Odd Trend' with points pointtype 12
#	    sin(2*pi*x) title 'Exact Solution' with lines
