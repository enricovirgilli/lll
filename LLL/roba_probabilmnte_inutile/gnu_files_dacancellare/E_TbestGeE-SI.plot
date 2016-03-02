set sty data histeps
set bmargin 4
set xlabel 'Energy (keV)' font 'Times, 24'
set ylabel 'T_{best} cm' font 'Times, 24'
set title '' font 'Times, 24'
set ter pos enh colo solid 'Times' 20
set ou "EtbestSIGE.eps"
pl[100:350][0.1:2.5]\
"EtbestSi111.dat" u 1:2 t "Bend Si(111)"  with points,\
"EtbestSi220.dat" u 1:2 t "Bend Si(220)"  with points,\
"EtbestGe111.dat" u 1:2 t "Bend Ge(111)"  with points,\
"EtbestGe220.dat" u 1:2 t "Bend Ge(220)"  with points
!evince EtbestSIGE.eps &
