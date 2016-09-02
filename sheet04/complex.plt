set xlabel "Re(z)"
set ylabel "Im(z)"
set zrange [0.9:1.1]
set view 0,0,1
set size square
splot 'data.txt' with dots title "Convergence to (1,0) depending on starting values"
