set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels 
   set xtic auto                          # set xtics automatically 
   set ytic auto                          # set ytics automatically
   set title "Error profile of a point Gauss Seidel solution of Laplace equation"
   set xlabel "x"
   set ylabel "y"
#   set key 0.01,100
#   set label "Yield Point" at 0.003,260	
#      set arrow from 0.0028,250 to 0.003,280
#      set xr [0.0:0.022]
#      set yr [0:325]
#   set log y
   set term post "Helvetica" 18
   set output 'errors.eps'
   set hidden3d
      splot    "E0" title 'Error (Computed - Exact)' with lines
