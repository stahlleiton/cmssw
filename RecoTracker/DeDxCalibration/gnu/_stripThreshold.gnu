load "myStyleEps.gnu"

set key top left

set style data his

###############################################
set output "../eps/stripProps/threshold.eps"

set xlabel "Cluster deposit [ADC]"
set ylabel "Fraction"

set xrange [-0.5:50.5]
set ytics 0,0.001

plot_arrow(thr,lt) = \
 sprintf("set arrow from %f,graph 0.6 to %f,graph 0.4 lt %d",\
 thr+0.01,thr+0.01,lt)

thr = system(sprintf("grep TIB ../data/stripProps.par | awk '{print $2}'"))+0
eval plot_arrow(thr, 1) # TIB

thr = system(sprintf("grep TID ../data/stripProps.par | awk '{print $2}'"))+0
eval plot_arrow(thr, 2) # TID

thr = system(sprintf("grep TOB ../data/stripProps.par | awk '{print $2}'"))+0
eval plot_arrow(thr, 3) # TOB

thr = system(sprintf("grep TEC3 ../data/stripProps.par | awk '{print $2}'"))+0
eval plot_arrow(thr, 4) # TEC3

thr = system(sprintf("grep TEC5 ../data/stripProps.par | awk '{print $2}'"))+0
eval plot_arrow(thr, 5) # TEC5

gunzip(det) = sprintf("<zcat ../out/coupling_%s.dat.gz | \
   awk '{if(NF==8) { a[$5]+=$3; s+=$3 } } END {for(x=0; x<254; x++) \
   print x,(a[x]+0)/s}'" , det)

plot \
 gunzip("TIB")  u 1:2 t "TIB"  ls TIB, \
 gunzip("TID")  u 1:2 t "TID"  ls TID, \
 gunzip("TOB")  u 1:2 t "TOB"  ls TOB, \
 gunzip("TEC3") u 1:2 t "TEC3" ls TEC3, \
 gunzip("TEC5") u 1:2 t "TEC5" ls TEC5

unset arrow

########################################
set output "../eps/stripProps/total.eps"
set xrange [-0.5:255.5]
set ytics 0,0.002
replot
