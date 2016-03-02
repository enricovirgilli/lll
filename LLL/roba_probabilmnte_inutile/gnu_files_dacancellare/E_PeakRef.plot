set sty data histeps
set bmargin 4
set xlabel 'Energy (keV)' font 'Times, 24'
set ylabel 'Reflectivity' font 'Times, 24'
set title '' font 'Times, 24'
set ter pos enh colo solid 'Times' 20
set ou "tbestSi111.eps"
pl[][0.62:0.83]\
"E_PeakReflSi111.dat" u 1:2 t "Bend Si(111)"  with points,\
"E_PeakReflSi220.dat" u 1:2 t "Bend Si(220)"  with points
!evince tbestSi111.eps &
