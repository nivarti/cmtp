set   autoscale                        # scale axes automatically
#unset log                              # remove any log-scaling
unset label                            # remove any previous labels 
   set xtic auto                          # set xtics automatically 
   set ytic auto                          # set ytics automatically
#   set title "Error Norms with Iteration"
   set xlabel "Change in Solution"
   set ylabel "y direction"
#   set key 0.01,100
#   set label "Slope = 1.9906" at 0.02, 0.005
#   set log xy
#   set log x
   set key left bottom
set grid
   set term post "Helvetica" 18
   set output 'gci2.eps'
#   set terminal epslatex
#  set output 'gci.tex'
   plot "gci" using 2:1 title '20 x 20 vs. 10 x 10' with points pointtype 10,\
      "gci" using 3:1 title '40 x 40 vs. 20 x 20' with points pointtype 8,\
	 "gci" using 4:1 title '80 x 80 vs. 40 x 40' with points pointtype 6,\
	    "gci" using 5:1 title '160 x 16 vs. 80 x 80' with points pointtype 4
