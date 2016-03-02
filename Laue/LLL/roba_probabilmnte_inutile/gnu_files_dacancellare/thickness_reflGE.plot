set sty data histeps
set bmargin 4
set xlabel 'Thickness (cm)' font 'Times, 24'
set ylabel 'Reflectivity' font 'Times, 24'
set title '' font 'Times, 24'
set ter pos enh colo solid 'Times' 20
set ou "EtbestSIGE.eps"
pl[0:4][]\
"thickness_reflGe220_100keV.dat" u 1:2 t "Bend Ge(220) 100 keV"  with points,\
"thickness_reflGe220_200keV.dat" u 1:2 t "Bend Ge(220) 200 keV"  with points,\
"thickness_reflGe220_300keV.dat" u 1:2 t "Bend Ge(220) 300 keV"  with points,\
"thickness_reflGe111_100keV.dat" u 1:2 t "Bend Ge(111) 100 keV"  with points,\
"thickness_reflGe111_200keV.dat" u 1:2 t "Bend Ge(111) 200 keV"  with points,\
"thickness_reflGe111_300keV.dat" u 1:2 t "Bend Ge(111) 300 keV"  with points
!evince EtbestSIGE.eps &
