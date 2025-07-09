set label 1 "" at graph 0.6, graph 0.95 right front

unset colorbox

gunzip(charge) = \
 sprintf("<zcat ../out/histos_%s_%s_2.dat.gz", ARG1,charge)

#
set output sprintf("../eps/histo_%s_%s_pos.eps",ARG1,ARG2)
set label 1 "positives" front

splot gunzip("pos") \
  u (($1)):2:(ARG2 eq "lin" ? $3 : ($3 > 0 ? log($3) : 0)) w image, \
  (u),log(epsilon(exp(u)/mel)),1 t @elp w l lt 5 lc @green, \
  (u),log(epsilon(exp(u)/mpi)),1 t @pip w l lt 1, \
  (u),log(epsilon(exp(u)/mka)),1 t @kap w l lt 4, \
  (u),log(epsilon(exp(u)/mpr)),1 t @prp w l lt 3

#
set output sprintf("../eps/histo_%s_%s_neg.eps",ARG1,ARG2)
set label 1 "negatives" front

splot gunzip("neg") \
  u (($1)):2:(ARG2 eq "lin" ? $3 : ($3 > 0 ? log($3) : 0)) w image, \
  (u),log(epsilon(exp(u)/mel)),1 t @elm w l lt 5 lc @green, \
  (u),log(epsilon(exp(u)/mpi)),1 t @pim w l lt 1, \
  (u),log(epsilon(exp(u)/mka)),1 t @kam w l lt 4, \
  (u),log(epsilon(exp(u)/mpr)),1 t @prm w l lt 3
