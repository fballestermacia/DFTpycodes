set encoding iso_8859_1
#set terminal  postscript enhanced color
#set output 'surfdos_r.eps'
#set terminal pngcairo truecolor enhanced  font ", 60" size 1920, 1680
set terminal png truecolor enhanced font ", 60" size 1920, 1680
set output 'surfdos_r.png'
set palette defined (-10 "#194eff", 0 "white", 10 "red" )
#set palette rgbformulae 33,13,10
set style data linespoints
unset ztics
unset key
set pointsize 0.8
set pm3d
set border lw 3
set size 0.8, 1
set origin 0.1, 0
#set size ratio -1
#set view equal xyz
set view map
#set cbtics font ",48"
#set xtics font ",48"
#set ytics font ",48"
#set ylabel font ",48"
set ylabel "Energy (eV)"
#set xtics offset 0, -1
#set ylabel offset -6, 0 
set xrange [0:            3.08590]
set yrange [         100.00000:         109.80000]
set xtics ("G"  0.00000,"X"  0.51195,"M"  1.13843,"Y"  1.65038,"G"  2.27685,"M"  3.08590)
set arrow from  0.51195, 100.00000 to  0.51195, 109.80000 nohead front lw 3
set arrow from  1.13843, 100.00000 to  1.13843, 109.80000 nohead front lw 3
set arrow from  1.65038, 100.00000 to  1.65038, 109.80000 nohead front lw 3
set arrow from  2.27685, 100.00000 to  2.27685, 109.80000 nohead front lw 3
set pm3d interpolate 2,2
splot 'dos.dat_r' u 1:2:3 w pm3d
set terminal png truecolor enhanced font ", 30" size 1920, 1680
set output 'spindos_r.png'
set multiplot layout 3, 1
set title 'sx'
splot 'spindos.dat_r' u 1:2:3 w pm3d 
set title 'sy'
splot 'spindos.dat_r' u 1:2:4 w pm3d
set title 'sz'
splot 'spindos.dat_r' u 1:2:5 w pm3d
