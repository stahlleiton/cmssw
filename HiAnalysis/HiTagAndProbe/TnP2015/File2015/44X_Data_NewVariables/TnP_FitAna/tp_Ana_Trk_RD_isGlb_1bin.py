import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')

process.source = cms.Source("EmptySource")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )    

PDFName = "twoGaussPlusPoly6v1"

process.TagProbeFitTreeAnalyzer = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
    # IO parameters:
    InputFileNames = cms.vstring(" /afs/cern.ch/user/c/chflores/work/public/TnP_2015/TP_Prod_Samples/Data/tnp_Prod_Data_PbPb_AllMB_25082015.root "),
    InputDirectoryName = cms.string("MuonTrk"),
    InputTreeName = cms.string("fitter_tree"),
    #numbrer of CPUs to use for fitting
    OutputFileName = cms.string("tnp_Ana_RD_PbPb_MuonTrk_AllMB_20150901_isGlb_1bin.root"),
    NumCPU = cms.uint32(1),
    # specifies wether to save the RooWorkspace containing the data for each bin and
    # the pdf object with the initial and final state snapshots
    SaveWorkspace = cms.bool(True),
    binsForMassPlots = cms.uint32(45),
    binnedFit = cms.bool(True),
    binsForFit = cms.uint32(45),
    #WeightVariable = cms.string("weight"),
    
    # defines all the real variables of the probes available in the input tree; can be used to select a subset of the probes
    Variables = cms.PSet(
        mass             = cms.vstring("Tag-Probe Mass", "2.0", "5.0", "GeV/c^{2}"),
        pt               = cms.vstring("Probe p_{T}", "0.0", "1000", "GeV/c"),
        p                = cms.vstring("Probe p", "0", "1000", "GeV/c"),
        eta              = cms.vstring("Probe #eta", "-2.4", "2.4", ""),
        abseta           = cms.vstring("Probe |#eta|", "0", "2.5", ""),
        staQoverP        = cms.vstring("Probe Q/p", "-200", "200", ""),
        staQoverPerror   = cms.vstring("Probe Q/p error", "0", "650", ""),
        staValidStations = cms.vstring("Probe ValidStations", "0", "5", ""),

        tag_pt           = cms.vstring("Tag p_{T}", "0.0", "1000", "GeV/c"),
        tag_eta          = cms.vstring("Tag #eta", "-2.4", "2.4", ""),
        tag_abseta       = cms.vstring("Tag |#eta|", "0", "2.5", ""),
        pair_pt          = cms.vstring("Pair p_{T}", "0.0", "1000", "GeV/c"),
        pair_absy        = cms.vstring("Pair |Y|", "-2.4", "2.4", ""),
        pair_y           = cms.vstring("Pair Y", "0", "2.5", ""),
    ),
    # defines all the Flags on which one can test the probe against (if true, is 'pass', if false is 'failed')
    Categories = cms.PSet(
        GlobalCuts      = cms.vstring("GlobalCuts", "dummy[true=1,false=0]"),
        GlobalCuts2      = cms.vstring("GlobalCuts2", "dummy[true=1,false=0]"),
        GlobalMu        = cms.vstring("GlobalMu", "dummy[true=1,false=0]"),
        PassingSta      = cms.vstring("PassingSta", "dummy[true=1,false=0]"),
        QualityMu       = cms.vstring("QualityMu", "dummy[true=1,false=0]"),
        StaTkSameCharge = cms.vstring("StaTkSameCharge", "dummy[true=1,false=0]"),
        TMA             = cms.vstring("TMA", "dummy[true=1,false=0]"),
        TMLSAT          = cms.vstring("TMLSAT", "dummy[true=1,false=0]"),
        Tk              = cms.vstring("Tk", "dummy[true=1,false=0]"),
        TrackCuts       = cms.vstring("TrackCuts", "dummy[true=1,false=0]"),
        TrackerMu       = cms.vstring("TrackerMu", "dummy[true=1,false=0]"),
        OuterValidHits  = cms.vstring("OuterValidHits", "dummy[true=1,false=0]"),
    ),

    # defines all the PDFs that will be available for the efficiency calculations; uses RooFit's "factory" syntax;
    # each pdf needs to define "signal", "backgroundPass", "backgroundFail" pdfs, "efficiency[0.9,0,1]" and "signalFractionInPassing[0.9]" are used for initial values  
    PDFs = cms.PSet(
      cbPlusExpo = cms.vstring(
        #"CBShape::signal(mass, mean[3.1,3.0,3.2], sigma[0.5], alpha[2.0, 0.2, 4.0], n[4, 0.5, 100.])",
        "CBShape::signal(mass, mean[3.1,3.0,3.2], sigma[0.02,0.02,0.1], alpha[1.0, 0.2, 3.0], n[4, 0.5, 100.])",
        "Exponential::backgroundPass(mass, lp[0,-5,5])",
        "Exponential::backgroundFail(mass, lf[0,-5,5])",
        "efficiency[0.9,0,1]",
        "signalFractionInPassing[0.9]"
      ),
      cbGausPlusExpo = cms.vstring(
        "Gaussian::signal1(mass, mean[3.1,3.0,3.2], sigma[0.02, 0.01, 0.1])",
        "CBShape::signal2(mass, mean, sigma2[0.02, 0.02, 0.3], alpha[2.0, 1.0, 10.0], n[4, 0.5, 100.])",
        "SUM::signalPass(fracP[0.8,0,1]*signal1,signal2)",
        "SUM::signalFail(fracF[0.8,0,1]*signal1,signal2)",
        "Exponential::backgroundPass(mass, lp[0,-5,5])",
        "Exponential::backgroundFail(mass, lf[0,-5,5])",
        "efficiency[0.9,0,1]",
        "signalFractionInPassing[0.9]"
      ),
      twoGaussPlusPoly6v1 = cms.vstring(
          "Gaussian::signal1(mass, mean[3.1,3.0,3.2], sigma[0.10,0.05,1.0])",
          "Gaussian::signal2(mass, mean1[3.7,3.5,3.9], sigma2[0.15,0.05,1.00])",
          "SUM::signal(vFrac[0.8,0,1]*signal1, signal2)",
          "Chebychev::backgroundPass(mass, {cP[0,-0.4,0.4], cP2[0.0,-0.04,0.04], cP3[-0.031,-0.3,0.3]})", ### good
          "Chebychev::backgroundFail(mass, {cF[-0.33,-1.0,1.0], cF2[0.05,-1.0,1.0], cF3[0.02,-1.0,1.0]})", ### good
          "efficiency[0.9,0,1]",
          "signalFractionInPassing[0.9]"
      ),
    ),
   # defines a set of efficiency calculations, what PDF to use for fitting and how to bin the data;
    # there will be a separate output directory for each calculation that includes a simultaneous fit, side band subtraction and counting. 
    Efficiencies = cms.PSet(
            #the name of the parameter set becomes the name of the directory
            # for multiple passing flags in EfficiencyCategorAndState = cms.vstring("flag1","state","flag2","state",...),
            isGlb_1bin = cms.PSet(
                EfficiencyCategoryAndState = cms.vstring("GlobalMu","true","GlobalCuts2","true"),
                UnbinnedVariables = cms.vstring("mass"),
                BinnedVariables = cms.PSet(
                    pt  = cms.vdouble(0,20),
                    eta = cms.vdouble(-2.4,2.4),
                ),
                BinToPDFmap = cms.vstring(PDFName)
            ),
    )
)

process.fitness = cms.Path(
        process.TagProbeFitTreeAnalyzer
)

