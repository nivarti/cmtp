set   autoscale                        # scale axes automatically
#unset log                              # remove any log-scaling
unset label                            # remove any previous labels 
   set xtic auto                          # set xtics automatically 
   set ytic auto                          # set ytics automatically
#   set title "Error Norms with Iteration"
   set xlabel "Iterations"
   set ylabel "Change in Solution"
#   set key 0.01,100
#   set label "Slope = 1.9906" at 0.02, 0.005
#   set log xy
   set log x
   set grid
   set term post "Helvetica" 18
   set output 'conv4.eps'
#   set terminal epslatex
#  set output 'conv.tex'
   plot "conv4" using 1:2 title 'P' with points pointtype 4,\
      "conv4" using 1:3 title 'u' with points pointtype 6,\
	 "conv4" using 1:4 title 'v' with points pointtype 8

