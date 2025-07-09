file(post) = sprintf("../%s/%s", energy, post)

load "myStyleEps.gnu"

set mxtics 5
unset key

set view map

set bmargin 4

set tics out

set xtics out offset 0,0.5
set ytics out offset 1.0,0

#
color1(gray) = sqrt(gray)     + exp(-gray/0.01)
color2(gray) = gray**3        + exp(-gray/0.01)
color3(gray) = sin(2*pi*gray) + exp(-gray/0.01)

set palette model RGB functions color1(gray), color2(gray), color3(gray)

set xlabel "Main deposit [ADC]" offset 0,1
set ylabel "Coupled deposit [ADC]" offset 1,0

set xrange [0:254]
set yrange [0:60]

set parametric
set urange [0:255]

r(a) = a/(1-2*a)

call "_stripProps.gnu" "TIB"
call "_stripProps.gnu" "TID"
call "_stripProps.gnu" "TOB"
call "_stripProps.gnu" "TEC3"
call "_stripProps.gnu" "TEC5"
