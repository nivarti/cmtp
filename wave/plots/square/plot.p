set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels 
   set xtic auto                          # set xtics automatically 
   set ytic auto                          # set ytics automatically
   set title "Discontinuous Solution and Flux Limiting"
   set xlabel "Position, x"
   set ylabel "Temperature"
#   set key 0.01,100
#   set label "Slope = 1.9858979" at 0.1,0.001	
#      set arrow from 0.0028,250 to 0.003,280
#      set xr [0.0:0.022]
#      set yr [0:325]
#   set log xy
   set grid
   set term post "Helvetica" 18
   set output 'Square.eps'
   plot    "Ted" using 1:2 title 'Exact Solution' with lines, \
            "Tuld" using 1:2 title 'Unlimited' with linesp pointtype 4, \
	       "Tld" using 1:2 title 'Superbee' with linespoints pointtype 12
