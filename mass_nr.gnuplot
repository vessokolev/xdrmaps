#!/usr/bin/env gnuplot

reset

set encoding iso_8859_1
set terminal postscript eps enhanced colour lw 1.425 font "Times-Bold,20"
set size square

set style increment userstyles
set style line 2 linewidth 4 
set border linewidth 5

unset key

# border
set style line 11 lc rgb '#000000' lt 1
set border 3 front ls 11
set tics nomirror out scale 1.5

set pm3d map
set palette color

set palette defined(\
0.00 "#000000", \
0.05 "#0000ff", \
0.10 "#00aa00", \
0.15 "#ffff00", \
0.20 "#ff0000", \
0.25 "#ff00ff", \
0.30 "#ffffff")


set xrange [0:10.1335] noreverse nowriteback
set yrange [0:10.1335] noreverse nowriteback
set cbrange [0.0:0.16]

set xlabel '{/Times-Bold-Italic x} [nm]' offset 0,0.2
set ylabel '{/Times-Bold-Italic y} [nm]' offset 1.34
set cblabel 'Atomic Mass [g]' offset 0.0

set mxtics 5
set mytics 5

#set xtics format " "
#set ytics format " "
set tics scale 1.5,0.5

#set xtics("0" 260.0,"20" 438.757375911,"40" 617.514751822,"60" 796.272127733,"80" 975.029503645,"100" 1153.78687956,277.875737591 1,295.751475182 1,313.627212773 1,331.502950364 1,349.378687956 1,367.254425547 1,385.130163138 1,403.005900729 1,420.88163832 1,456.633113502 1,474.508851093 1,492.384588684 1,510.260326276 1,528.136063867 1,546.011801458 1,563.887539049 1,581.76327664 1,599.639014231 1,635.390489413 1,653.266227005 1,671.141964596 1,689.017702187 1,706.893439778 1,724.769177369 1,742.64491496 1,760.520652551 1,778.396390142 1,814.147865325 1,832.023602916 1,849.899340507 1,867.775078098 1,885.650815689 1,903.52655328 1,921.402290871 1,939.278028462 1,957.153766053 1,992.905241236 1,1010.78097883 1,1028.65671642 1,1046.53245401 1,1064.4081916 1,1082.28392919 1,1100.15966678 1,1118.03540437 1,1135.91114196 1,1171.66261715 1,1189.53835474 1,1207.41409233 1,1225.28982992 1,1243.16556751 1,1261.0413051 1,1278.91704269 1)

#set ytics("0" 0.0,"20" 178.624590652,"40" 357.249181304,"60" 535.873771956,"80" 714.498362608,"100" 893.12295326,17.8624590652 1,35.7249181304 1,53.5873771956 1,71.4498362608 1,89.312295326 1,107.174754391 1,125.037213456 1,142.899672522 1,160.762131587 1,196.487049717 1,214.349508782 1,232.211967848 1,250.074426913 1,267.936885978 1,285.799345043 1,303.661804108 1,321.524263174 1,339.386722239 1,375.111640369 1,392.974099434 1,410.8365585 1,428.699017565 1,446.56147663 1,464.423935695 1,482.28639476 1,500.148853826 1,518.011312891 1,553.736231021 1,571.598690086 1,589.461149152 1,607.323608217 1,625.186067282 1,643.048526347 1,660.910985412 1,678.773444478 1,696.635903543 1,732.360821673 1,750.223280738 1,768.085739804 1,785.948198869 1,803.810657934 1,821.673116999 1,839.535576064 1,857.39803513 1,875.260494195 1)

#set cbtics ("0.0" 0.0,"0.4" 0.4e-5, "0.8" 0.8e-5, "1.2" 1.2e-5, "1.6" 1.6e-5, "2.0" 2.0e-5)

set title 'Atomic Mass' offset 0,-0.5

set lmargin 1.0
set rmargin 1.0

set bmargin 3.5
set tmargin 1.75

set output "mass_nr.eps" # at z-slice 1300

plot 'mass_nr.dat' u ($1*0.0200267):($2*0.0200267):($3/3.0) matrix with image
