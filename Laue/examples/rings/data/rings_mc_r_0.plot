set sty data histeps
set xlabel 'cm' font 'Times, 24'
set ylabel 'Normalized counts' font 'Times, 24'
set title '' font 'Times, 24'
set ter pos enh colo solid 'Times' 20
set ou "rings_mc_r_0.eps"
pl [:] "rings_mc_r_0.dat" u 1:($2/21.165851593) t "" lw 3
!gv rings_mc_r_0.eps &