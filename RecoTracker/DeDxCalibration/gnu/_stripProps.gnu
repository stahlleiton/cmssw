set datafile missing "?"

set output sprintf("../eps/stripProps/%s.eps",ARG1)

set label 1 ARG1 at graph 0.9,0.075 right front # tc rgb "#808080"

thr   = system(sprintf("grep %s ../data/stripProps.par | awk '{print $2}'",ARG1)) + 0
alpha = system(sprintf("grep %s ../data/stripProps.par | awk '{print $3}'",ARG1)) + 0
sigma = system(sprintf("grep %s ../data/stripProps.par | awk '{print $4}'",ARG1)) + 0

sum = 1e+9
if(ARG1 eq "TIB" ) sum = 110
if(ARG1 eq "TID" ) sum =  90
if(ARG1 eq "TEC3") sum =  90
if(ARG1 eq "TEC5") sum = 150

aa = 0
if(ARG1 eq "TOB" || ARG1 eq "TEC5") aa = 30
if(ARG1 eq "TIB"                  ) aa = 25
if(ARG1 eq "TID" || ARG1 eq "TEC3") aa = 20

set style fill solid 1.0 noborder
splot sprintf("<zcat ../out/coupling_%s.dat.gz",ARG1) u 1:2:3 w image, \
       u, r(alpha)*u, 1, \
       u, thr+u*0,    1, \
       u,(sum-u < 30 ? sum-u : 1/0), 1 lt 0, \
       u,aa,          1 lt 3
