set encoding iso_8859_1
#set terminal  postscript enhanced color
#set output 'slabek.eps'
#set terminal  pngcairo truecolor enhanced  font ",60" size 1920, 1680
set terminal  png truecolor enhanced  font ",60" size 1920, 1680
set output 'slabek.png'
set palette defined ( 0  "green", 5 "yellow", 10 "red" )
set style data linespoints
unset ztics
unset key
set pointsize 0.8
set border lw 3 
set view 0,0
#set xtics font ",36"
#set ytics font ",36"
#set ylabel font ",36"
#set xtics offset 0, -1
set ylabel offset -1, 0 
set xrange [0:    3.08590]
set ylabel "Energy (eV)"
set yrange [  -5.64012: 202.45723]
set xtics ("  G"    0.00000,"  X"    0.27091,"  M"    0.60243,"  Y"    0.87334,"  G"    1.20486,"  M"    1.63299)
set arrow from    0.27091,  -5.64012 to    0.27091, 202.45723 nohead
set arrow from    0.60243,  -5.64012 to    0.60243, 202.45723 nohead
set arrow from    0.87334,  -5.64012 to    0.87334, 202.45723 nohead
set arrow from    1.20486,  -5.64012 to    1.20486, 202.45723 nohead
#rgb(r,g,b) = int(r)*65536 + int(g)*256 + int(b)
#plot 'slabek.dat' u 1:2:(rgb(255,$3, 3)) w lp lw 2 pt 7  ps 1 lc rgb variable
# (a) 
# plot the top and bottom surface's weight together
#plot 'slabek.dat' u 1:2:($3+$4) w lp lw 2 pt 7  ps 1 lc palette
# (b) 
# plot top and bottom surface's weight with red and blue respectively
set palette defined ( -1  "blue", 0 "grey", 1 "red" )
plot 'slabek.dat' u 1:2:($4-$3) w lp lw 2 pt 7  ps 1 lc palette
