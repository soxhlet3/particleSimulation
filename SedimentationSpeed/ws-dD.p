reset
set grid
set terminal png enhanced
set logscale x
set format x "10^{%L}"
set xrange [1:1000]
set xlabel "D/d [-]"

ws(x) = 2<x ? sqrt(1/(1+2.104/x)) : 1/0



set output "ws-dD.png"
set ylabel "w_s/w_s_,_0 [-]"
set yrange [0.6:1]
set ytics 0.6,0.05,1
plot ws(x) with lines lw 2 notitle axes x1y1 \
     
	  



     
  

