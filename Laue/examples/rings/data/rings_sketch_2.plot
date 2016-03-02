set data sty dots
set bmargin 4
set size square
set xlabel 'cm' font 'Times, 24'
set ylabel 'cm' font 'Times, 24'
set title '' font 'Times, 24'
set ter pos enh colo solid 'Times' 20
set polar
set ou "rings_xtalinfo_2.eps"
pl  "rings_xtalinfo_2.dat" u 4:3 t "" lw 3
!gv rings_xtalinfo_2.eps &