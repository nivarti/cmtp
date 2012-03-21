set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels 
   set xtic auto                          # set xtics automatically 
   set ytic auto                          # set ytics automatically
   set title "Temperature (Solution) Profile for Various CFL Numbers"
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
   set output 'Stability.eps'
   plot    "Ts" using 1:2 title 'Stable (CFL = 0.3)' with lines, \
            "Tb" using 1:2 title 'Borderline (CFL = 0.525)' with linesp pointtype 6, \
	       "Tb2" using 1:2 title 'Slightly Unstable (CFL = 0.528)' with linesp pointtype 12, \
		  "Tu" using 1:2 title 'Unstable (CFL = 0.54)' with linesp pointtype 4
