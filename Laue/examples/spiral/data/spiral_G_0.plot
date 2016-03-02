set sty data histeps
set bmargin 4
set xlabel 'Energy - keV' font 'Times, 24'
set ylabel 'Focusing factor' font 'Times, 24'
set title '' font 'Times, 24'
set ter pos enh colo solid 'Times' 20
set ou "spiral_G_0.eps"
pl "spiral_G_0.dat" u 1:2 t "" lw 3
!gv spiral_G_0.eps &
