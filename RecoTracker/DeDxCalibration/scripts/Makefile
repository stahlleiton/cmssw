all :

# source /cvmfs/cms.cern.ch/crab3/crab.csh

sub :
	rm -f trees/hadronTree_*.root
	rm -rf crab_reco
	crab submit -c conf.py

stat :
	crab status -d crab_reco --verboseErrors

kill :

get :
	crab getoutput -d crab_reco --checksum=no --outputpath /eos/user/s/sikler/DeDxCalib
	make list

list :
	ls -d -1 /afs/cern.ch/user/s/sikler/DeDxCalibration/scripts/trees/* | \
	  awk '{print "file://"$$1}' > ../test/list.txt
	wc -l ../test/list.txt

clean :
	rm -f trees/*.root
	rm -rf crab_reco
