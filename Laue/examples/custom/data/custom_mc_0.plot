set pm3 at s
set size square

set format z ""
set format cb ""

unset ztics
unset cbtics
unset cbmtics

set xlabel "y - cm"
set ylabel "x - cm"
set cblabel "" 100,0.

set ticslevel 0
unset surface

set yrange [*:*]# reverse
set border 0

set view 180,90

sp [][][:] 'custom_mc_0.dat' u ($1*0.2-(5.0)):($2*0.2-(5.0)):3 matr w l
pause -1

set ter png enh font 'Times' 20
set ou "custom_mc_0.png"
rep
