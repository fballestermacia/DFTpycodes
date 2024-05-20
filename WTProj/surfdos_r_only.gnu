set encoding iso_8859_1
#set terminal  postscript enhanced color
#set output 'surfdos_r.eps'
#set terminal pngcairo truecolor enhanced  font ", 60" size 1920, 1680
set terminal png truecolor enhanced font ", 60" size 1920, 1680
set output 'surfdos_r_only.png'
set palette defined (0  "white", 6 "red", 20 "black" )
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
unset cbtics
#set xtics offset 0, -1
#set ylabel offset -6, 0 
set xrange [0:            1.63299]
set yrange [         100.00000:         109.80000]
set xtics ("G"  0.00000,"X"  0.27091,"M"  0.60243,"Y"  0.87334,"G"  1.20486,"M"  1.63299)
set arrow from  0.27091, 100.00000 to  0.27091, 109.80000 nohead front lw 3
set arrow from  0.60243, 100.00000 to  0.60243, 109.80000 nohead front lw 3
set arrow from  0.87334, 100.00000 to  0.87334, 109.80000 nohead front lw 3
set arrow from  1.20486, 100.00000 to  1.20486, 109.80000 nohead front lw 3
set pm3d interpolate 2,2
splot 'dos.dat_r' u 1:2:(exp($4)) w pm3d
