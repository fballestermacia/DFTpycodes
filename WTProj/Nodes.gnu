set encoding iso_8859_1
#set terminal  postscript enhanced color
#set output 'gap.eps'
set terminal  png      truecolor enhanced size 1920, 1680 font ",36"
set output 'Nodes.png'
unset ztics
unset key
set ticslevel 0
#set view equal xyz
set border 4095 front lt black linewidth 2 dashtype solid
set xlabel "k_1"
set ylabel "k_2"
set zlabel "k_3"
set xtics nomirror scale 0.5
set ytics nomirror scale 0.5
set ztics nomirror scale 0.5
set xtics offset  0.0,-1.0 , 0
set ytics offset -2.5,   0 , 0
set size ratio -1
#set view 60, 140, 1, 1
set xrange [  -0.50000:   0.50000]
set yrange [  -0.50000:   0.50000]
set zrange [  -0.50000:   0.50000]
splot 'Nodes.dat' u 6:7:8 w p pt 7 ps 2
