set sty data histeps
set xlabel 'cm' font 'Times, 24'
set ylabel 'Normalized counts' font 'Times, 24'
set title '' font 'Times, 24'
set ter pos enh colo solid 'Times' 20
set ou "rings_mc_xy_3.eps"
pl [:] "rings_mc_x_3.dat" u 1:($2/1530.42053223) t "x" lw 3,"rings_mc_y_3.dat" u 1:($2/1530.42053223) t "y" lw 3 lt 3
!gv rings_mc_xy_3.eps &