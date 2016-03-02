set data sty lin
set bmargin 4
# set log
set xlabel '' font 'Times, 24'
set ylabel '' font 'Times, 24'
set title '' font 'Times, 24'
set ter pos enh colo solid 'Times' 20
set ou "table.eps"
pl  \
"table.dat" t "" lw 3
!gv table.eps &