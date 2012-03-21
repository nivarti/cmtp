set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels 
   set xtic auto                          # set xtics automatically 
   set ytic auto                          # set ytics automatically
   set title "Temperature (Solution) Profile for Different Mesh Sizes"
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
   set output 'Temp.eps'
      plot    "T1" using 1:2 title '20 Cells' with linesp pointtype 4, \
	 "T2" using 1:2 title '40 Cells' with linesp pointtype 12, \
	    sin(2*pi*x) title 'Exact Solution' with lines
