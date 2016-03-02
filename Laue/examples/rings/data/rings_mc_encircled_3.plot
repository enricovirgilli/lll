set sty data histeps
set xlabel 'cm' font 'Times, 24'
set ylabel 'Encircled photon fraction' font 'Times, 24'
set title '' font 'Times, 24'
set ter pos enh colo solid 'Times' 20
set ou "rings_mc_encircled_3.eps"
pl [:] "rings_mc_encircled_3.dat" u 1:($2/1530.42041016) t "" lw 3
!gv rings_mc_encircled_3.eps &