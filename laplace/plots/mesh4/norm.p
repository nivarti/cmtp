set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels 
   set xtic auto                          # set xtics automatically 
   set ytic auto                          # set ytics automatically
   set title "Error Norm variation with Mesh Size"
   set xlabel "Cell Size"
   set ylabel "Error Norm"
#   set key 0.01,100
#   set label "Yield Point" at 0.003,260	
#      set arrow from 0.0028,250 to 0.003,280
#      set xr [0.0:0.022]
#      set yr [0:325]
   set log xy
   set grid
   set term post "Helvetica" 14
   set output 'norm.eps'
      plot    "phi" using 1:2 title 'L2 Norm' with linespoints
