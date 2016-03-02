set sty data histeps
set xlabel 'cm' font 'Times, 24'
set ylabel 'Normalized counts' font 'Times, 24'
set title '' font 'Times, 24'
set ter pos enh colo solid 'Times' 20
set ou "rings_mc_y_3.eps"
pl [:] "rings_mc_y_3.dat" u 1:($2/1530.42053223) t "" lw 3
!gv rings_mc_y_3.eps &