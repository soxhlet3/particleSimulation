reset
set grid
set terminal png enhanced
set xrange [0:0.1]
set xtics  0,0.02,0.1
set xlabel "Time [s]"

set output "v_a_t.png"
set ylabel "Speed in z [m/s]"
set yrange [-0.1:0]
set ytics  -0.1,0.01,0
set y2label "Acceleration in z [m/s²]"
set y2range [-10:5]
set y2tics  -10,1.5,5
set key right top
plot "gnuplot.dat" using 1:3 with lines lw 3 title "v_z" , \
     "gnuplot.dat" using 1:5 with lines lw 3 title "a_z" axes x1y2, \
	  0 with lines lw 2 axes x1y2 title "0 m/s²"


set xrange [0:0.01]
set xtics  0,0.002,0.01
set xlabel "z [m]"

set output "v_a_z.png"
set ylabel "Speed in z [m/s]"
set yrange [-0.1:0]
set ytics  -0.1,0.01,0
set y2label "Acceleration in Z [m/s²]"
set y2range [-10:5]
set y2tics  -10,1.5,5
set key right top
plot "gnuplot.dat" using 2:3 with lines lw 3 title "v_z" , \
     "gnuplot.dat" using 2:5 with lines lw 3 title "a_z" axes x1y2, \
     0 with lines lw 2 axes x1y2 title "0 m/s²"
  

