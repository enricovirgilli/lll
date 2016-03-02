set sty data histeps
set xlabel 'cm' font 'Times, 24'
set ylabel 'Normalized counts' font 'Times, 24'
set title '' font 'Times, 24'
set ter pos enh colo solid 'Times' 20
set ou "rings_mc_xy_0.eps"
pl [:] "rings_mc_x_0.dat" u 1:($2/21.1658496857) t "x" lw 3,"rings_mc_y_0.dat" u 1:($2/21.1658496857) t "y" lw 3 lt 3
!gv rings_mc_xy_0.eps &