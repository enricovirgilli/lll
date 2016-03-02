set data sty lin
set bmargin 4
set log
set xlabel 'Energy - keV' font 'Times, 24'
set ylabel 'counts/(sec cm^3 keV)' font 'Times, 24'
set title 'SPI background' font 'Times, 24'
set ter pos enh colo solid 'Times' 20
set ou "spibg.eps"
pl [:2400] \
"spibg.dat" t "" lw 3
!gv spibg.eps &