set sty data histeps
set xlabel 'cm' font 'Times, 24'
set ylabel 'Normalized counts' font 'Times, 24'
set title '' font 'Times, 24'
set ter pos enh colo solid 'Times' 20
set ou "rings_mc_xy_5.eps"
pl [:] "rings_mc_x_5.dat" u 1:($2/153.193862915) t "x" lw 3,"rings_mc_y_5.dat" u 1:($2/153.193862915) t "y" lw 3 lt 3
!gv rings_mc_xy_5.eps &