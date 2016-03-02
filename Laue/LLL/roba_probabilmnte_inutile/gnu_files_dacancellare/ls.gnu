#!/usr/local/bin/gnuplot -persist
#
#    
#    	G N U P L O T
#    	Version 4.2 patchlevel 5 
#    	last modified Mar 2009
#    	System: Linux 2.6.22.19-server-2mdv
#    
#    	Copyright (C) 1986 - 1993, 1998, 2004, 2007 - 2009
#    	Thomas Williams, Colin Kelley and many others
#    
#    	Type `help` to access the on-line reference manual.
#    	The gnuplot FAQ is available from http://www.gnuplot.info/faq/
#    
#    	Send bug reports and suggestions to <http://sourceforge.net/projects/gnuplot>
#    
set terminal pos col enhanced font "Times New Roman, 17"
set output "ls105s.ps"
set size 1.1, 1.1
set origin 14.0, 2.0

set key inside right top vertical Right noreverse enhanced autotitles nobox
set key noinvert samplen 4 spacing 1 width 0 height 0 
set title "Line Sensitivity"
set xlabel "Energy (keV)" 
set xlabel  offset character 0, 0, 0 font "Times New Roman, 20" textcolor lt -1 norotate
set x2label "" 
set x2label  offset character 0, 0, 0 font "Times New Roman, 20" textcolor lt -1 norotate
set xrange [ 95.6482 : 303.756 ] noreverse nowriteback
set x2range [ 84.7169 : 269.040 ] noreverse nowriteback
set ylabel "Photons cm^{-2} s^{-1}" 
set ylabel  offset character 0, 0, 0 font "Times New Roman, 20" textcolor lt -1 rotate by 90
set y2label "" 
set y2label  offset character 0, 0, 0 font "" textcolor lt -1 rotate by 90
set yrange [ -0.000107168 : 0.000191849 ] noreverse nowriteback
set y2range [ -8.89319e-05 : 0.000195527 ] noreverse nowriteback

set lmargin  13
set bmargin  5
set rmargin  13
set tmargin  5
set locale "C"
set pm3d explicit at s
set pm3d scansautomatic
set pm3d interpolate 1,1 flush begin noftriangles nohidden3d corners2color mean
set palette positive nops_allcF maxcolors 0 gamma 1.5 color model RGB 
set palette rgbformulae 7, 5, 15
set colorbox default
set colorbox vertical origin screen 0.9, 0.2, 0 size screen 0.05, 0.6, 0 front bdefault
set grid
set loadpath 
set fontpath 
set fit noerrorvariables
plot [][1e-6:4e-5]"linesens_binned_105s_err.dat" using 1:3:2 w xerrorbars  t "7 bins" pt 7 lt 1 lw 3, "ls_lens.dat" using 1:2  t "" w l lt -1 lw 2.5
#    EOF
