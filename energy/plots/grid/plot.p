set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels 
   set xtic auto                          # set xtics automatically 
   set ytic auto                          # set ytics automatically
   set xlabel "x direction"
   set ylabel "Temperature Gradient"
#   set key 0.01,100
#   set label "Slope    = 1.9906" at 0.02, 0.005
#   set log xy
   set grid
   set term post "Helvetica" 18
   set output 'grid.eps'
   plot "dTdy0" using 1:2 title 'Mesh: 25x10' with lines, \
       "dTdy1" using 1:2 title 'Mesh: 50x20' with lines, \
	  "dTdy2" using 1:2 title 'Mesh: 100x40' with lines, \
	    "dTdy3" using 1:2 title 'Mesh: 200x80' with lines
