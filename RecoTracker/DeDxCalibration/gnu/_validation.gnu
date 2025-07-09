ls = 1

ymax = 0.15
if(ARG2 == 0.51) { ymax = 0.05 }
if(ARG2 == 0.70) { ymax = 0.07 }
if(ARG2 == 0.86) { ymax = 0.08 }
if(ARG2 == 0.93) { ymax = 0.09 }
if(ARG2 == 1.19) { ymax = 0.10 }
if(ARG2 == 1.39) { ymax = 0.10 }
if(ARG2 == 2.08) { ymax = 0.12 }
if(ARG2 == 3.49) { ymax = 0.14 }

set yrange [1e-5:ymax]

set output sprintf("../eps/validation/%s_%s_%s.eps",ARG1,ARG4,ARG3)

set multiplot

first = 1

pixel = (("$0" eq "PXB") || ("$0" eq "PXF"))

take = 1
 call "__validation.gnu" ARG1 ARG2  300e-4 ARG4 ARG3
take = 0

take = 1
if(("$0" eq "PXF") || ("$0" eq "TID") || ("$0" eq "TEC3")) \
 call "__validation.gnu" "$0" $1  350e-4 "$3" "$2"
take = 0

unset label

take = 1
 call "__validation.gnu" ARG1 ARG2  400e-4 ARG4 ARG3
take = 0

take = 1
 call "__validation.gnu" ARG1 ARG2  500e-4 ARG4 ARG3
take = 0

take = 1
 call "__validation.gnu" ARG1 ARG2  600e-4 ARG4 ARG3
take = 0

if(!pixel) take = 1
 call "__validation.gnu" ARG1 ARG2  750e-4 ARG4 ARG3
take = 0

if(!pixel) take = 1
 call "__validation.gnu" ARG1 ARG2  900e-4 ARG4 ARG3
take = 0

unset label
set nomultiplot
