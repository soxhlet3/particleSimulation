reset
set grid
set terminal png enhanced
set logscale x
set format x "10^{%L}"
set xrange [1:100]
set xlabel "D/d [-]"

ws(x) = 2<x ? sqrt(1/(1+2.104/x)) : 1/0



set output "ws-dD-sim.png"
set ylabel "w_s/w_s_,_0 [-]"
set yrange [0.3:1.3]
set ytics 0.3,0.1,1.3
set key right bottom
set style line 1 lc rgb 'black' pt 5   # square
set style line 2 lc rgb 'red' pt 6   # square
 plot ws(x) with lines lw 2 title 'Empirical correlation by Ladenburg' axes x1y1, \
 	'-' w p ls 1 lw 4 title 'Resolved' axes x1y1, '-' w p ls 1 lw 4 notitle axes x1y1, '-' w p ls 1 lw 4 notitle axes x1y1, '-' w p ls 1 lw 4 notitle axes x1y1, '-' w p ls 2 lw 4 title 'Sub-grid' axes x1y1, '-' w p ls 2 lw 4 notitle axes x1y1, '-' w p ls 2 lw 4 notitle axes x1y1, '-' w p ls 2 lw 4 notitle axes x1y1
 2 0.3287836246
 e
 5 0.83073917
 e
 10 1.024849
 e
 30 1.1467
 e
 
 2 0.99978
 e
 5 0.99978
 e
 10 0.99978
 e
 30 0.99978
 e





     
  

