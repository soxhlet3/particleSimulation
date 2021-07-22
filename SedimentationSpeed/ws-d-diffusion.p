reset
set grid
set terminal png enhanced
set logscale x
set format x "10^{%L}"
set xrange [0.001:1000]
set xlabel "d [μm]"

k = 1.3804e-23
p = 1.0e5
dM = 3.62e-10
T = 293


ws(x) = 1+(k*T/(pi*sqrt(2)*dM*dM*p)/(x*1.0e-6))*(2.514+0.8*exp(-0.55/(k*T/(pi*sqrt(2)*dM*dM*p)/(x*1.0e-6))))

set output "ws-d-diffusion.png"
set logscale y
set format y "10^{%L}"
set ylabel "w_s/w_s_,_0 [-]"
set yrange [0.1:1000]

plot ws(x) with lines lw 2 notitle axes x1y1 \
