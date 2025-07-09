all :

src ::
	scram b -j4

try :
	cd scripts ; cmsRun reco.py

cal :
	cd test; nice make -j3

doc :
	cd gnu ; nice make -j3
	cd eps ; nice make -j3
	cd tex ; make pdf

conf :
	cd scripts ; \
	cmsDriver.py try.py \
	 -s RAW2DIGI,RECO \
	 --conditions=132X_dataRun3_Prompt_v4 \
	 --eventcontent=RECO \
	 --filein='/store/hidata/HIRun2023A/HIForward4/RAW/v1/000/374/634/00000/fd7ffc7d-846f-410c-8996-73c6afcf38fd.root' \
	 --fileout=a.root \
	 -n 10 --data --no_exec --repacked --era=Run3_2023
	cd scripts ; cmsRun try_py_RAW2DIGI_RECO.py 

clean :
	#cd scripts; make clean
	cd data; make clean
	cd test; make clean
	cd gnu ; make clean
	cd eps ; make clean
	cd tex ; make clean
