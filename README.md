# VertexCompositeAnalysis

Example of setting up and running gamma+gamma to dimuon tree

cmsrel CMSSW_12_5_3

cd CMSSW_12_5_3/src

cmsenv

git clone -b ParticleFitter_12_5_X_UPC https://github.com/stahlleiton/VertexCompositeAnalysis

cd VertexCompositeAnalysis

scram b -j8

cd VertexCompositeProducer/test
