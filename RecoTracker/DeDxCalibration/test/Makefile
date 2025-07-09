# from 1 to npro

list = 1 2 3 4 5 6
nfiles = $(shell wc -l < list.txt)
args = -nfiles $(nfiles) -nchunks 6

#
all : .extractDeposit .extractHistos

#
.extractDeposit : ../out/data_2.root
	extractDeposit
	gzip -f ../out/depositCurves/*.dat
	touch $@

#
.extractHistos : .extractHistos_-1 \
                 .extractHistos_0 \
                 .extractHistos_1 \
                 .extractHistos_2
	gzip -f ../out/coupling_*.dat
	gzip -f ../out/histos_*.dat
	touch $@

.extractHistos_% : ../out/data_%.root
	extractHistos -pass $*
	touch $@

#############
# Pass 2

# Process for surveyMap and depositMap and fitProbs
../out/data_2_%.root : ../data/gain_1.dat
	processData $(args) \
          -chunk $(word 3,$(subst _, ,$(patsubst %.root,%,$@))) \
          -pass 2 -applyGain -surveyMap -depositMap -fitProbs

#############
# Pass 1

# Re-calibrate
../data/gain_1.dat : ../out/data_1.root
	calibrateGain $(args) -pass 1

# Re-process for calibGain and surveyMap
../out/data_1_%.root : ../data/stripProps.par ../data/gain_0.dat
	processData $(args) \
          -chunk $(word 3,$(subst _, ,$(patsubst %.root,%,$@))) \
          -pass 1 -applyGain -calibGain -surveyMap

#############
# Pass 0

# Calibrate
../data/gain_0.dat : ../out/data_0.root
	calibrateGain $(args) -pass 0

# Process for calibGain and surveyMap
../out/data_0_%.root : ../data/stripProps.par
	processData $(args) \
          -chunk $(word 3,$(subst _, ,$(patsubst %.root,%,$@))) \
          -pass 0 -calibGain -surveyMap

#############
# Pass -1

# Strip props
../data/stripProps.par : ../data/stripThresholds.par ../out/data_-1.root
	evalStripProps

#############
# Process for coupledHits
../out/data_-1_%.root : list.txt
	processData $(args) \
          -chunk $(word 3,$(subst _, ,$(patsubst %.root,%,$@))) \
          -pass -1 -coupledHits

#############
# Root
../out/data_%.root : \
       $(foreach i,$(list),../out/data_%_$(i).root)
	hadd -f $@ $^

clean :
	rm -f ../out/*.dat.gz
	rm -f ../out/*.root
	rm -f ../out/*.bin
	rm -f ../out/depositCurves/*.dat.gz
	rm -f .extractDeposit .extractHistos*
	rm -f list.txt
