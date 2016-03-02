set sty data histeps
set bmargin 4
set xlabel 'Energy - keV' font 'Times, 24'
set ylabel 'Focusing factor' font 'Times, 24'
set title '' font 'Times, 24'
set ter pos enh colo solid 'Times' 20
set ou "rings_G_5.eps"
pl "rings_G_5.dat" u 1:2 t "" lw 3
!gv rings_G_5.eps &
