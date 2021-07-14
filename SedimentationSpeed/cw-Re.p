reset
set grid
set terminal png enhanced
set logscale xy
set format xy "10^{%L}"
set xrange [0.01:200000]
set xlabel "Re [-]"

cw(x) = x<0.25 ? (24/x) : 0.25<=x<1000 ? (24/x + 4/sqrt(x) + 0.4) : 1000<=x<200000 ? (24/x + 5.66/sqrt(x) + 0.33) : 1/0



set output "cw.png"
set ylabel "c_w(Re) [-]"
set yrange [0.01:10000]
plot cw(x) with lines lw 2 notitle axes x1y1 \
     
	  



     
  

