set data sty dots
set bmargin 4
set size square
set xlabel 'cm' font 'Times, 24'
set ylabel 'cm' font 'Times, 24'
set title '' font 'Times, 24'
set ter pos enh colo solid 'Times' 20
set polar
set ou "lens_xtalinfo_0.eps"
pl  "lens_xtalinfo_0.dat" u 4:3 t "" lw 3
!gv lens_xtalinfo_0.eps &