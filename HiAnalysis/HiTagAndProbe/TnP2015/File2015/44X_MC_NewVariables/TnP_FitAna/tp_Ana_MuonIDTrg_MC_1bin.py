import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("TagProbe")

# setup 'analysis'  options
options = VarParsing.VarParsing ('analysis')
options.outputFile = 'output.root'
options.inputFiles = 'input.root'
options.parseArguments()

process.load('FWCore.MessageService.MessageLogger_cfi')

process.source = cms.Source("EmptySource")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )    

PDFName = "cbGausPlusExpo"
#PDFName = "cbPlusPoly2nd"

process.TagProbeFitTreeAnalyzer = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
    # IO parameters:
    InputFileNames = cms.vstring(options.inputFiles),
    InputDirectoryName = cms.string("MuonIDTrg"),
    InputTreeName = cms.string("fitter_tree"),
    OutputFileName = cms.string(options.outputFile),
    #numbrer of CPUs to use for fitting
    NumCPU = cms.uint32(25),
    # specifies wether to save the RooWorkspace containing the data for each bin and
    # the pdf object with the initial and final state snapshots
    #binnedFit = cms.bool(True),
    #binsForFit = cms.uint32(30),
    #binsForMassPlots = cms.uint32(50),
    binsForMassPlots = cms.uint32(50),
    #binnedFit = cms.bool(True),
    #binsForFit = cms.uint32(50),
    WeightVariable = cms.string("weight"),
    SaveWorkspace = cms.bool(True),
    
    # defines all the real variables of the probes available in the input tree and intended for use in the efficiencies
    Variables = cms.PSet(
                         mass             = cms.vstring("Tag-Probe Mass", "2.6", "3.5", "GeV/c^{2}"),
                         pt               = cms.vstring("Probe p_{T}", "0.0", "1000", "GeV/c"),
                         p                = cms.vstring("Probe p", "0", "1000", "GeV/c"),
                         eta              = cms.vstring("Probe #eta", "-2.4", "2.4", ""),
                         abseta           = cms.vstring("Probe |#eta|", "0", "2.5", ""),
                         glbChi2          = cms.vstring("Probe glbChi2", "0", "80", ""),
                         glbPtError       = cms.vstring("Probe glbPtError", "0", "0.5", ""),
                         glbSigmaPtOverPt = cms.vstring("Probe #sigma(p_{T})/p_{T}", "0", ".03", ""),
                         l1dr             = cms.vstring("Probe l1dr", "0", "15", ""),
                         l1pt             = cms.vstring("Probe l1pt", "0", "140", ""),
                         l1q              = cms.vstring("Probe l1dr", "0", "8", ""),
                         tkChi2           = cms.vstring("Probe tkChi2", "0", "2", ""),
                         tkExpHitIn       = cms.vstring("Probe tkExpHitIn", "0", "2", ""),
                         tkExpHitOut      = cms.vstring("Probe tkExpHitOut", "0", "6", ""),
                         tkHitFract       = cms.vstring("Probe tkHitFract", "0", "1.5", ""),
                         tkPixelLay       = cms.vstring("Probe tkPixelLay", "0", "5", ""),
                         tkPtError        = cms.vstring("Probe tkPtError", "0", ".3", ""),
                         tkSigmaPtoverPt  = cms.vstring("Probe sigmaPt/pt", "0", ".03", ""),
                         tkValidHits      = cms.vstring("Probe tkValidHits", "0", "25", ""),
                         tkValidPixelHits = cms.vstring("Probe tkValPixHit", "0", "6", ""),
                         tmaVar           = cms.vstring("Probe tmaVar", "0", "2", ""),
                         tmlsatvar        = cms.vstring("Probe tmlsatvar", "0", "2", ""),
                         tackerMuVar      = cms.vstring("Probe trackerMuVar", "0", "2", ""),
                         
                         tag_pt     = cms.vstring("Tag p_{T}", "0.0", "1000", "GeV/c"),
                         tag_eta    = cms.vstring("Tag #eta", "-2.4", "2.4", ""),
                         tag_abseta = cms.vstring("Tag |#eta|", "0", "2.5", ""),
                         pair_pt    = cms.vstring("Pair p_{T}", "0.0", "1000", "GeV/c"),
                         pair_absy  = cms.vstring("Pair |Y|", "-2.4", "2.4", ""),
                         pair_y     = cms.vstring("Pair Y", "0", "2.5", ""),
                         weight     = cms.vstring("weight","0","1000",""),

    ),
    # defines all the discrete variables of the probes available in the input tree and intended for use in the efficiency calculations
    Categories = cms.PSet(
                          GlobalCuts= cms.vstring("GlobalCuts", "dummy[true=1,false=0]"),
                          GlobalMu  = cms.vstring("GlobalMu", "dummy[true=1,false=0]"),
                          QualityCuts = cms.vstring("QualityCuts", "dummy[true=1,false=0]"),
                          HLTL1v0   = cms.vstring("HLTL1v0", "dummy[true=1,false=0]"),
                          HLTL1v1   = cms.vstring("HLTL1v1", "dummy[true=1,false=0]"),
                          HLTL1v2   = cms.vstring("HLTL1v2", "dummy[true=1,false=0]"),
                          TMA       = cms.vstring("TMA", "dummy[true=1,false=0]"),
                          TMLSAT    = cms.vstring("TMLSAT", "dummy[true=1,false=0]"),
                          Tk        = cms.vstring("Tk", "dummy[true=1,false=0]"),
                          TrackCuts = cms.vstring("TrackCuts", "dummy[true=1,false=0]"),
                          TrackerMu = cms.vstring("TrackerMu", "dummy[true=1,false=0]"),
    ),

    # defines all the PDFs that will be available for the efficiency calculations; uses RooFit's "factory" syntax;
    # each pdf needs to define "signal", "backgroundPass", "backgroundFail" pdfs, "efficiency[0.9,0,1]" and "signalFractionInPassing[0.9]" are used for initial values  
    PDFs = cms.PSet(
      cbPlusExpo = cms.vstring(
        "CBShape::signal(mass, mean[3.1,3.0,3.2], sigma[0.5], alpha[2.0, 0.2, 4.0], n[4, 0.5, 100.])",
        "Exponential::backgroundPass(mass, lp[0,-5,5])",
        "Exponential::backgroundFail(mass, lf[0,-5,5])",
        "efficiency[0.9,0,1]",
        "signalFractionInPassing[0.9]"
      ),
      cbPlusPoly = cms.vstring(
        "CBShape::signal(mass, mean[3.1,3.0,3.2], sigma[0.02, 0.01, 0.1], alpha[2.0, 1.0, 5.0], n[1, 0.5, 100.])",
        #"CBShape::signal(mass, mean[3.1,3.0,3.2], sigma[0.02, 0.01, 0.1], alpha[1.0, 0.5, 4.0], n[1, 0.5, 100.])",
        "Chebychev::backgroundPass(mass, {cPass[0,-5.0,5.0], cPass2[0,-5.0,5.0]})",
        "Chebychev::backgroundFail(mass, {cFail[0,-5.0,5.0], cFail2[0,-5.0,5.0]})",
        #"Chebychev::backgroundPass(mass, {cPass[0,-0.5,0.5], cPass2[0,-0.5,0.5]})",
        #"Chebychev::backgroundFail(mass, {cFail[0,-0.5,0.5], cFail2[0,-0.5,0.5]})",
        "efficiency[0.9,0,1]",
        "signalFractionInPassing[0.9]"
      ),
      cbPlusPoly2 = cms.vstring(
        "Gaussian::signal1(mass, mean[3.1,3.0,3.2], sigma[0.02, 0.01, 0.1])",
        "CBShape::signal2(mass, mean, sigma2[0.02, 0.02, 0.3], alpha[2.0, 1.0, 10.0], n[4, 0.5, 100.])",
        "SUM::signalPass(fracP[0.8,0,1]*signal1,signal2)",
        "SUM::signalFail(fracF[0.8,0,1]*signal1,signal2)",
        #"CBShape::signal1(mass, mean[3.1,3.0,3.2], sigma[0.02, 0.01, 0.1], alpha[2.0, 1.0, 5.0], n[1, 0.5, 100.])",
        #"Gaussian::signal2(mass, mean, sigma2[0.02, 0.01, 0.1])",
        #"SUM::signal(vFrac[0.8,0,1]*signal1,signal2)",
        #"CBShape::signal(mass, mean[3.1,3.0,3.2], sigma[0.02, 0.01, 0.1], alpha[1.0, 0.5, 4.0], n[1, 0.5, 100.])",
        "Chebychev::backgroundPass(mass, {cPass[0,-5.0,5.0], cPass2[0,-5.0,5.0]})",
        "Chebychev::backgroundFail(mass, {cFail[0,-5.0,5.0], cFail2[0,-5.0,5.0]})",
        #"Chebychev::backgroundPass(mass, {cPass[0,-0.5,0.5], cPass2[0,-0.5,0.5]})",
        #"Chebychev::backgroundFail(mass, {cFail[0,-0.5,0.5], cFail2[0,-0.5,0.5]})",
        "efficiency[0.9,0,1]",
        "signalFractionInPassing[0.9]"
      ),
      cbGausPlusExpo = cms.vstring(
        "Gaussian::signal1(mass, mean[3.1,3.0,3.2], sigma[0.02, 0.01, 0.1])",
        "CBShape::signal2(mass, mean, sigma2[0.02, 0.02, 0.3], alpha[2.0, 1.0, 10.0], n[4, 0.5, 100.])",
        "SUM::signalPass(fracP[0.8,0,1]*signal1,signal2)",
        "SUM::signalFail(fracF[0.8,0,1]*signal1,signal2)",
        "Exponential::backgroundPass(mass, lp[0,-5.0,5.0])",
        "Exponential::backgroundFail(mass, lf[0,-5.0,5.0])",  # same slope, they're both muons
        "efficiency[0.9,0,1]",
        "signalFractionInPassing[0.9]"
      ),
      cbPlusPoly2nd = cms.vstring(
        "CBShape::signal(mass, mean[3.1,3.0,3.2], sigma[0.02, 0.01, 0.1], alpha[2.0, 1.0, 10.0], n[4, 0.5, 100.])",
        ##"Gaussian::signal(mass, mean[3.1,3.0,3.2], sigma[0.02, 0.01, 0.1])",
        "Chebychev::backgroundPass(mass, {cPass[0,-0.5,0.5], cPass2[0,-0.5,0.5]})",
        "Chebychev::backgroundFail(mass, {cFail[0,-0.5,0.5], cFail2[0,-0.5,0.5]})",
        #"Exponential::backgroundPass(mass, lp[0,-5,5])",
        #"Exponential::backgroundFail(mass, lf[0,-5,5])",  # same slope, they're both muons
        "efficiency[0.9,0,1]",
        "signalFractionInPassing[0.9]"
      ),
      gaussPlusExpo = cms.vstring(
        "Gaussian::signal(mass, mean[3.1,3.0,3.2], sigma[0.02, 0.01, 0.1])",
        "Chebychev::backgroundPass(mass, {cPass[0,-0.5,0.5], cPass2[0,-0.5,0.5]})",
        "Chebychev::backgroundFail(mass, {cFail[0,-0.5,0.5], cFail2[0,-0.5,0.5]})",
        "efficiency[0.9,0,1]",
        "signalFractionInPassing[0.9]"
      ),
    ),
    # defines a set of efficiency calculations, what PDF to use for fitting and how to bin the data;
    # there will be a separate output directory for each calculation that includes a simultaneous fit, side band subtraction and counting. 
    Efficiencies = cms.PSet(
          #the name of the parameter set becomes the name of the directory
          PassingGlb_1bin = cms.PSet(
             EfficiencyCategoryAndState = cms.vstring("QualityCuts","true","HLTL1v0","true","HLTL1v1","true","HLTL1v2","true"),
             UnbinnedVariables = cms.vstring("mass","weight"),
             BinnedVariables = cms.PSet(
                pt = cms.vdouble(1.5,30),
                eta = cms.vdouble(-2.4, 2.4),
                ),
             BinToPDFmap = cms.vstring(PDFName)
             ),
          )
    )

process.fitness = cms.Path(
    process.TagProbeFitTreeAnalyzer
)

