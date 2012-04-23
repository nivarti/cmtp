set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels 
   set xtic auto                          # set xtics automatically 
   set ytic auto                          # set ytics automatically
#   set title "Error Profile for Solution"
   set xlabel "x"
   set ylabel "y"
#   set key 0.01,100	
   set term post "Helvetica" 18
   set output 'v.eps'
   set hidden3d
	splot    "v" title 'v' with lines
