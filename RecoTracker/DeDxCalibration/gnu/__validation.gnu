bg = ARG2 + 0 # $1
l  = ARG3 + 0 # $2

set xrange [0:5]

Delta = epsilon(bg) * l

ok = system(sprintf("zcat ../out/depositCurves/sig_%s_%s_unknown_%s.dat.gz | awk '{if($5==%s) print $2}' | wc -l", ARG1,ARG4,ARG2,ARG3))

if(ok == 1) \
 sn = system(sprintf("zcat ../out/depositCurves/sig_%s_%s_unknown_%s.dat.gz | awk '{if($5==%s) print $2}'", ARG1,ARG4,ARG2,ARG3)); \
 else sn = 0.

over = system(sprintf("zcat ../out/depositCurves/his_%s_%s_unknown_%s.dat.gz | awk '{if($1==%s && $3>0) s++} END {print s+0}'", ARG1,ARG4,ARG2,ARG3))

mov = 0

amp = 0.001

gunzip(gg) = sprintf("<zcat ../out/depositCurves/his_%s_%s_unknown_%s.dat.gz | awk '{if($1 == %s) print}'", ARG1,ARG4,ARG2,ARG3)

set datafile missing "0"

set fit quiet
if(over > mov) \
 fit [0.05:] prob(x) gunzip(0) u 2:3 via amp

unset key

if(over > mov && first == 1 && take == 1) \
 set label 2 sprintf("%s   %s   {/Symbol bg} = %.2f",ARG1,ARG4,bg) at graph 0.75,0.2 center; first = 0

if(over > mov && take == 1) \
       set key at graph 0.95,1.025-0.075*ls Right nobox width -5 ; \
       else set nokey

set xrange [0:1]

if(over > mov && ok == 1 && ("$4" ne "x") && take == 1) \
 plot [0:0.6] "" u 2:3:4 w e t \
  sprintf("l = %.0f {/Symbol m}m, {/Symbol s}_n = %.1f keV",l*1e+4,sn*1e+3) \
  ls ls, \
  prob(x) ls ls

if(take) if(ls == 1) unset label 2

if(take) ls = ls + 1
