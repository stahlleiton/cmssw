set key top right Left reverse samplen 2 width +1 box noauto 
set bar small

# Margin
set lmargin at screen 0.15
set rmargin at screen 0.95
set bmargin at screen 0.175
set tmargin at screen 0.95

# Terminal
set term post eps enh color dashed dl 2 "Helvetica" 25 size 5in,4.5in

set pointsize 1.5
set tics scale 2,1

set mxtics 5
set mytics 5

set macros
set fit logfile "/dev/null" errorvariables

# Styles, colors
green = "rgb \"dark-green\""

ltz = "rgb \"dark-green\""


set style line  1 lt 1 lw 3 lc 1      pt 6
set style line  2 lt 2 lw 3 lc @green pt 8 
set style line  3 lt 3 lw 3 lc 3      pt 4
set style line  4 lt 4 lw 3 lc 4      pt 12

PXB  = 51; set style line PXB  lt 1 lw 3 lc  0 pt 13
PXF  = 52; set style line PXF  lt 2 lw 3 lc  0 pt 12

TIB  = 53; set style line TIB  lt 1 lw 3 lc  1 pt 7
TID  = 54; set style line TID  lt 2 lw 3 lc @ltz pt 6
TOB  = 55; set style line TOB  lt 3 lw 3 lc  3 pt 9
TEC3 = 56; set style line TEC3 lt 4 lw 3 lc  4 pt 5
TEC5 = 57; set style line TEC5 lt 5 lw 3 lc  5 pt 8

