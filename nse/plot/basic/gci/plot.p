set   autoscale                        # scale axes automatically
#unset log                              # remove any log-scaling
unset label                            # remove any previous labels 
   set xtic auto                          # set xtics automatically 
   set ytic auto                          # set ytics automatically
#   set title "Error Norms with Iteration"
   set xlabel "x Velocity"
   set ylabel "y direction"
#   set key 0.01,100
#   set label "Slope = 1.9906" at 0.02, 0.005
#   set log xy
#   set log x
   set key right bottom
set grid
   set term post "Helvetica" 18
   set output 'gci.eps'
#   set terminal epslatex
#  set output 'gci.tex'
   plot "mesh0" using 2:1 title '10 x 10' with points pointtype 10,\
      "mesh1" using 2:1 title '20 x 20' with points pointtype 8,\
	 "mesh2" using 2:1 title '40 x 40' with points pointtype 6,\
	    "mesh3" using 2:1 title '80 x 80' with points pointtype 4
