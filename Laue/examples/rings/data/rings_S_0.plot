set sty data histeps
set bmargin 4
set xlabel 'Energy - keV' font 'Times, 24'
set ylabel 'Equivalent Surface - cm^2' font 'Times, 24'
set title '' font 'Times, 24'
set ter pos enh colo solid 'Times' 20
set ou "rings_S_0.eps"
pl "rings_S_0.dat" u 1:2 t "" lw 3
!gv rings_S_0.eps &
