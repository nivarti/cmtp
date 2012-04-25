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
   set output 'convU.eps'
#   set terminal epslatex
#  set output 'conv.tex'
   plot "conv" using 1:3 title 'SOR = 1.0, dt = 0.05' with linespoints pointtype 4,\
      "conv2" using 1:3 title 'SOR = 1.5, dt = 0.05' with linespoints pointtype 6,\
	 "conv3" using 1:3 title 'SOR = 1.0, dt = 0.25' with linespoints pointtype 8,\
	    "conv4" using 1:3 title 'SOR = 1.5, t = 0.25' with linespoints pointtype 10

