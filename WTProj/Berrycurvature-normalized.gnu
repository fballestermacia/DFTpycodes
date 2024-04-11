set encoding iso_8859_1
#set terminal  pngcairo  truecolor enhanced size 1920, 1680 font ",40"
set terminal  png       truecolor enhanced size 1920, 1680 font ",40"
set output 'Berrycurvature-normalized.png'
unset ztics
unset key
set border lw 3
set xlabel 'k (1/{\305})'
set ylabel 'k (1/{\305})'
unset colorbox
#unset xtics
#unset xlabel
set xrange [] noextend
set yrange [] noextend
set ytics 0.5 nomirror scale 0.5
set pm3d interpolate 2,2
set title '({/Symbol W}_x, {/Symbol W}_y)'
plot 'Berrycurvature-normalized.dat' u 1:2:($4/10):($5/10) w vec
