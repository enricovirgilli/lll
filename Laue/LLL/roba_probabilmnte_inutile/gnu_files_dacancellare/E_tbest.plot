set sty data histeps
set bmargin 4
set xlabel 'Energy (keV)' font 'Times, 24'
set ylabel 'T_{best} (cm)' font 'Times, 24'
set title '' font 'Times, 24'
set ter pos enh colo solid 'Times' 20
set ou "Etbest.eps"
pl [90:500][0.65:0.9]"E_tbestSi220.dat" u 1:2 t "Bend Si(220)"  with points,\
"E_tbestSi111.dat" u 1:2 t "Bend Si(111)"  with points
!evince Etbest.eps &
