reset
set grid
set terminal png enhanced
set xrange [10:1000]
set xlabel "T [K]"

k = 1.3804e-23
p = 1.0e5
dM = 3.62e-10
d = 7.0e-6

ws(x) = 1+(k*x/(pi*sqrt(2)*dM*dM*p)/d)*(2.514+0.8*exp(-0.55/(k*x/(pi*sqrt(2)*dM*dM*p)/d)))

set output "ws-T.png"
set ylabel "w_s/w_s_,_0 [-]"
set yrange [0.8:1.2]

plot ws(x) with lines lw 2 notitle axes x1y1 \
