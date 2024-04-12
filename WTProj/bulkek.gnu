set terminal pdf enhanced color font ",30"
set palette defined ( 0  "green", 5 "yellow", 10 "red" )
set output 'bulkek.pdf' 
set style data linespoints
unset key
set pointsize 0.8
#set xtics font ",24"
#set ytics font ",24"
set ylabel font ",24"
set ylabel offset 1.5,0
set xrange [0:    5.90622]
emin=   -0.500000
emax=   14.712297
set yrange [0: emax]
set ylabel "Frequency (THz)"
set xtics ("E  "    0.00000,"A  "    0.62648,"G  "    1.46715,"B  "    1.90978,"D  "    2.53626,"C  "    3.09287,"Z  "    3.65175,"G  "    4.27822,"Y  "    4.83711,"C  "    5.46358,"E  "    5.90622)
set arrow from    0.62648,  0.0 to    0.62648, emax nohead
set arrow from    1.46715,  0.0 to    1.46715, emax nohead
set arrow from    1.90978,  0.0 to    1.90978, emax nohead
set arrow from    2.53626,  0.0 to    2.53626, emax nohead
set arrow from    3.09287,  0.0 to    3.09287, emax nohead
set arrow from    3.65175,  0.0 to    3.65175, emax nohead
set arrow from    4.27822,  0.0 to    4.27822, emax nohead
set arrow from    4.83711,  0.0 to    4.83711, emax nohead
set arrow from    5.46358,  0.0 to    5.46358, emax nohead
# please comment the following lines to plot the fatband 
plot 'bulkek.dat' u 1:2  w lp lw 2 pt 7  ps 0.2 lc rgb 'black', 0 w l lw 2
 
# uncomment the following lines to plot the fatband 
#plot 'bulkek.dat' u 1:2:3  w lp lw 2 pt 7  ps 0.2 lc palette, 0 w l lw 2
# uncomment the following lines to plot the spin if necessary
#plot 'bulkek.dat' u 1:2 w lp lw 2 pt 7  ps 0.2, \
     'bulkek.dat' u 1:2:($3/6):($4/6) w vec
