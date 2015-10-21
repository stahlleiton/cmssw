import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')

process.source = cms.Source("EmptySource")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )    

PDFName = "GaussPlusPoly"
#PDFName = "cbPlusPoly"

process.TagProbeFitTreeAnalyzer = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
    # IO parameters:
    # InputFileNames = cms.vstring("tnp_Prod_PbPb_RD_STA_Cent40100.root"),
    InputFileNames = cms.vstring("tnp_Prod_PbPb_RD_STA_AllMB_24092015.root"),
    InputDirectoryName = cms.string("fitTkFromSta"),
    InputTreeName = cms.string("fitter_tree"),
    #numbrer of CPUs to use for fitting
    # OutputFileName = cms.string("tnp_Ana_PbPb_RD_STA_Cent40100.root"),
    OutputFileName = cms.string("tnp_Ana_PbPb_RD_STA_AllMB.root"),
    NumCPU = cms.uint32(1),
    # specifies wether to save the RooWorkspace containing the data for each bin and
    # the pdf object with the initial and final state snapshots
    SaveWorkspace = cms.bool(True),
    binsForMassPlots = cms.uint32(40),
    #binnedFit = cms.bool(True),
    #binsForFit = cms.uint32(45),
    #WeightVariable = cms.string("weight"),
    
    # defines all the real variables of the probes available in the input tree and intended for use in the efficiencies
    Variables = cms.PSet(
        mass = cms.vstring("Tag-Probe Mass", "2.8", "3.4", "GeV/c^{2}"),
        pt = cms.vstring("Probe p_{T}", "0", "1000", "GeV/c"),
        p = cms.vstring("Probe p", "0", "1000", "GeV/c"),
        eta = cms.vstring("Probe #eta", "-2.4", "2.4", ""),
        abseta = cms.vstring("Probe |#eta|", "0", "2.4", ""),
        tag_pt = cms.vstring("Tag p_{T}", "0", "1000", "GeV/c"),
        #weight = cms.vstring("weight","0.0","10000.0",""),
    ),
    # defines all the discrete variables of the probes available in the input tree and intended for use in the efficiency calculations
    Categories = cms.PSet(
        #mcTrue = cms.vstring("MC true", "dummy[true=1,false=0]"),
        passingAtVtx = cms.vstring("passingAtVtx", "dummy[true=1,false=0]"),
    ),

    # defines all the PDFs that will be available for the efficiency calculations; uses RooFit's "factory" syntax;
    # each pdf needs to define "signal", "backgroundPass", "backgroundFail" pdfs, "efficiency[0.9,0,1]" and "signalFractionInPassing[0.9]" are used for initial values  
    PDFs = cms.PSet(
      cbPlusExpo = cms.vstring(
        "CBShape::signal(mass, mean[3.1,3.0,3.2], sigma[0.5], alpha[2.0, 0.2, 4.0], n[4, 0.5, 100.])",
       # "CBShape::signal(mass, mean[3.1,3.0,3.2], sigma[0.02,0.02,0.1], alpha[1.0, 0.2, 3.0], n[4, 0.5, 100.])",
        "Exponential::backgroundPass(mass, lp[0,-5,5])",
        "Exponential::backgroundFail(mass, lf[0,-5,5])",
        "efficiency[0.9,0,1]",
        "signalFractionInPassing[0.9]"
      ),
      cbGausPlusExpo = cms.vstring(
       "Gaussian::signal1(mass, mean[3.15,3.05,3.25], sigma[0.02, 0.01, 0.1])",
       "CBShape::signal2(mass, mean1[3.1,3.0,3.2], sigma2[0.02, 0.02, 0.3], alpha[2.0, 1.0, 10.0], n[4, 0.5, 100.])",
        "SUM::signalPass(fracP[0.8,0,1]*signal1,signal2)",
        "SUM::signalFail(fracF[0.8,0,1]*signal1,signal2)",
        "Exponential::backgroundPass(mass, lp[0,-5,5])",
        "Exponential::backgroundFail(mass, lf[0,-5,5])",
        "efficiency[0.9,0,1]",
        "signalFractionInPassing[0.9]"
      ),
      GaussPlusPoly = cms.vstring(
          "Gaussian::signal(mass, mean[3.09,3.0,3.2], sigma[0.1,0.0,1.])",
          "Gaussian::signal2(mass, mean1[3.1,3.05,3.2], sigma2[0.01,0.05,.50])",
          "SUM::signal(vFrac[0.8,0,1]*signal1, signal2)",
          "Chebychev::backgroundPass(mass, {cP[0,-.5,.5]})",
          "Chebychev::backgroundFail(mass, {cF[0,-.5,.5]})",
          "efficiency[0.9,0,1]",
          "signalFractionInPassing[0.9]"
      ),
       cbPlusPoly = cms.vstring(
        "CBShape::signal(mass, mean[3.08,3.0,3.2], sigma[0.01, 0.0, 0.3], alpha[2.0, 1.0, 10.0], n[4, 0.5, 100.])",
        #"CBShape::signal(mass, mean[3.1,3.0,3.2], sigma[0.02, 0.01, 0.1], alpha[1.0, 0.5, 4.0], n[1, 0.5, 100.])",
        "Chebychev::backgroundPass(mass, {cPass[0,-.5,.5], cPass2[0,-.5,.5]})",
        "Chebychev::backgroundFail(mass, {cFail[0,-.5,.5], cFail2[0,-.5,.5]})",
        #"Chebychev::backgroundPass(mass, {cPass[0,-0.5,0.5], cPass2[0,-0.5,0.5]})",
        #"Chebychev::backgroundFail(mass, {cFail[0,-0.5,0.5], cFail2[0,-0.5,0.5]})",
        "efficiency[0.9,0,1]",
        "signalFractionInPassing[0.9]"
      ),
       cbGausPlusPoly = cms.vstring(
        "CBShape::signal1(mass, mean[3.08,3.0,3.15], sigma[0.015, 0.01, 0.1], alpha[1.5, 1.0, 10.0], n[2, 0.5, 100.])",
        "Gaussian::signal2(mass, mean2[3.1,3.05,3.2], sigma2[0.10,0.05,1.00])",
        "SUM::signalPass(fracP[0.8,0,1]*signal1, signal2)",
        "SUM::signalFail(fracF[0.8,0,1]*signal1,signal2)",
        "Chebychev::backgroundPass(mass, {cPass[0,-.5,.5], cPass2[0,-.5,.5],cPass3[0,-0.5,0.5]})",
        "Chebychev::backgroundFail(mass, {cFail[0,-.5,.5], cFail2[0,-.5,.5],cFail3[0,-0.5,0.5]})",
        #"Chebychev::backgroundPass(mass, {cPass[0,-0.5,0.5], cPass2[0,-0.5,0.5]})",
        #"Chebychev::backgroundFail(mass, {cFail[0,-0.5,0.5], cFail2[0,-0.5,0.5]})",
        "efficiency[0.9,0,1]",
        "signalFractionInPassing[0.9]"
      ),
    ),
    # defines a set of efficiency calculations, what PDF to use for fitting and how to bin the data;
    # there will be a separate output directory for each calculation that includes a simultaneous fit, side band subtraction and counting. 
    Efficiencies = cms.PSet(
        #the name of the parameter set becomes the name of the directory
        PassingSTA_pt1 = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("passingAtVtx","true"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet(
              # pt = cms.vdouble(2.1,4, 4.5, 5, 5.5, 8.5, 30),
              pt = cms.vdouble(2.1,3, 4, 5, 7, 30),
              abseta = cms.vdouble(0, 1.6), 
            ),
            BinToPDFmap = cms.vstring(PDFName)
        ),
        #  PassingSTA_pt2 = cms.PSet(
        #     EfficiencyCategoryAndState = cms.vstring("passingAtVtx","true"),
        #     UnbinnedVariables = cms.vstring("mass"),
        #     BinnedVariables = cms.PSet(
        #        pt = cms.vdouble(2.1,3.,4,5,6., 30),
        #        abseta = cms.vdouble(0.9, 1.6),    
        # ),
        #     BinToPDFmap = cms.vstring(PDFName)
        # ),
         PassingSTA_pt2 = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("passingAtVtx","true"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet(
               # pt = cms.vdouble(1.5,2.5,3.5, 4.5, 5.5, 30),
               pt = cms.vdouble(1.5,2.5,3.5, 5, 30),
               abseta = cms.vdouble(1.6, 2.4),
            ),
            BinToPDFmap = cms.vstring(PDFName)
        ),
       # PassingSTA_pt4 = cms.PSet(
       #      EfficiencyCategoryAndState = cms.vstring("passingAtVtx","true"),
       #      UnbinnedVariables = cms.vstring("mass"),
       #      BinnedVariables = cms.PSet(
       #        pt = cms.vdouble(1.5,2.5,3.5,5, 30),
       #        abseta = cms.vdouble(2.1, 2.4),
       #      ),
       #      BinToPDFmap = cms.vstring(PDFName)
       #  ),

        PassingSTA_1bin = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("passingAtVtx","true"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet(
                pt = cms.vdouble(1.5,30.),
                eta = cms.vdouble(-2.4,2.4),
            ),
            BinToPDFmap = cms.vstring(PDFName)
        ),
        PassingSTA_eta = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("passingAtVtx","true"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet(
               eta = cms.vdouble(-2.4,-1.6,0,1.6,2.4),
               pt = cms.vdouble(1.5,30.),
            ),
            BinToPDFmap = cms.vstring(PDFName)
        ),
        PassingSTA_abseta = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("passingAtVtx","true"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet(
                abseta = cms.vdouble(0.,1.6,2.4),
                pt = cms.vdouble(1.5,30),
            ),
            BinToPDFmap = cms.vstring(PDFName)
        ),

    )
)

process.fitness = cms.Path(
    process.TagProbeFitTreeAnalyzer
)

