set encoding iso_8859_1
set terminal  postscript enhanced color font ",30"
set output 'wcc.eps'
unset key 
set border lw 3 
set xtics offset 0, 0.0
set xtics format "%4.1f" nomirror out 
set xlabel "k" 
set xlabel offset 0, 0.7 
set ytics 0.5 
unset xtics 
set ytics format "%4.1f" nomirror out
set ylabel "WCC"
set ylabel offset 0, 0.0 
set xrange [0: 0.5]
set yrange [0:1]
plot for [i=4:    11] 'wcc.dat' u 1:i w p  pt 7  ps 1.1 lc 'red'
