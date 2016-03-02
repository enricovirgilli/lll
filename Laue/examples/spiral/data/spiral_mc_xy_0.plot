set sty data histeps
set xlabel 'cm' font 'Times, 24'
set ylabel 'Normalized counts' font 'Times, 24'
set title '' font 'Times, 24'
set ter pos enh colo solid 'Times' 20
set ou "spiral_mc_xy_0.eps"
pl [:] "spiral_mc_x_0.dat" u 1:($2/22.0227279663) t "x" lw 3,"spiral_mc_y_0.dat" u 1:($2/22.0227279663) t "y" lw 3 lt 3
!gv spiral_mc_xy_0.eps &