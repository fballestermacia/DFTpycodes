set encoding iso_8859_1
set terminal  pngcairo  truecolor enhanced size 3680, 1920 font ",40"
#set terminal  png       truecolor enhanced size 3680, 1920 font ",40"
set output 'Orbitalmagnetization.png'
if (!exists("MP_LEFT"))   MP_LEFT = .12
if (!exists("MP_RIGHT"))  MP_RIGHT = .92
if (!exists("MP_BOTTOM")) MP_BOTTOM = .12
if (!exists("MP_TOP"))    MP_TOP = .88
if (!exists("MP_GAP"))    MP_GAP = 0.08
set multiplot layout 1,3 rowsfirst \
              margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_GAP
 
set palette rgbformulae 33,13,10
unset ztics
unset key
set pm3d
#set zbrange [ -10: 10] 
#set cbrange [  -100 : 100 ] 
set view map
set size ratio -1
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
set title 'Orbital magnetization m_x ({/Symbol m}_B)'
splot 'Orbitalmagnetization.dat' u 1:2:4 w pm3d
unset ylabel
unset ytics
set title 'Orbital magnetization m_y ({/Symbol m}_B)'
splot 'Orbitalmagnetization.dat' u 1:2:5 w pm3d
set title 'Orbital magnetization m_{z} ({/Symbol m}_B)'
set colorbox
splot 'Orbitalmagnetization.dat' u 1:2:6 w pm3d
