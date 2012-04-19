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
   set output 'size.eps'
   plot "dTdy0" using 1:2 title 'L = 5.0' with points pointtype 4, \
      "dTdy1" using 1:2 title 'L = 10.0' with points pointtype 12, \
	 "dTdy2" using 1:2 title 'L = 20.0' with points pointtype 10, \
	    "dTdy3" using 1:2 title 'L = 40.0' with points pointtype 6
