set sty data histeps
set xlabel 'cm' font 'Times, 24'
set ylabel 'Normalized counts' font 'Times, 24'
set title '' font 'Times, 24'
set ter pos enh colo solid 'Times' 20
set ou "lens_mc_r_0.eps"
pl [:] "lens_mc_r_0.dat" u 1:($2/105.347206116) t "" lw 3
!gv lens_mc_r_0.eps &