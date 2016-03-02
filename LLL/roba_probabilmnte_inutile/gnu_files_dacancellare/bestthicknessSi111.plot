set sty data histeps
set bmargin 4
set xlabel 'thickness (cm)' font 'Times, 24'
set ylabel 'Reflectivity' font 'Times, 24'
set title '' font 'Times, 24'
set ter pos enh colo solid 'Times' 20
set ou "tbestSi111.eps"
pl [0:6]"bendSi111_100keV.dat" u 1:2 t "Bend Si(111) -  100 keV"  with points,\
"bendSi111_200keV.dat" u 1:2 t "Bend Si(111) -  200 keV"  with points,\
"bendSi111_300keV.dat" u 1:2 t "Bend Si(111) -  300 keV"  with points
!evince tbestSi111.eps &
