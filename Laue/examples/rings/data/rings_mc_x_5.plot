set sty data histeps
set xlabel 'cm' font 'Times, 24'
set ylabel 'Normalized counts' font 'Times, 24'
set title '' font 'Times, 24'
set ter pos enh colo solid 'Times' 20
set ou "rings_mc_x_5.eps"
pl [:] "rings_mc_x_5.dat" u 1:($2/153.193862915) t "" lw 3
!gv rings_mc_x_5.eps &