import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('python')

options.register('inputFilename', 'DataC_1.txt', #HTauTauAnalysis_1_1_Sl2.root',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Input file name"
)

options.register('outFilename', 'Output.root',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Output file name"
)

options.parseArguments()

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

#process.TFileService=cms.Service("TFileService",fileName=cms.string("Dec31_OUT_"+options.outFilename))

process.TFileService=cms.Service("TFileService",fileName=cms.string("1JAN_DATA_ZMU_OUT_"+options.outFilename))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.source = cms.Source("EmptySource")

process.tauhist = cms.EDAnalyzer('MuJetAnalysis',
                                 MCHistos = cms.string("data/MC25ns.root"),
                                 DataHistos = cms.string("data/Data25ns.root"),
                                 isolationTau = cms.string("LooseDB"),
                                 BTagInputFile = cms.string("data/CSVv250ns.csv"), 
                                 InputFile = cms.string(options.inputFilename),
                                 DYSample = cms.string(options.inputFilename),
                                 isData = cms.bool(True),    
                                 isSignalZ = cms.bool(False), # for Gen cuts
                                 isBkgZ = cms.bool(False),   # for gen Cuts *False for other DY bkgs
                                 L1Path = cms.string("data/Spring16_25nsV6_DATA_Uncertainty_AK4PFchs.txt"),  
                                 L2Path = cms.string("data/Spring16_25nsV6_DATA_Uncertainty_AK4PFchs.txt"),  
                                 L3Path = cms.string("data/Spring16_25nsV6_DATA_Uncertainty_AK4PFchs.txt"),
                                 L1DATAPath = cms.string("data/Spring16_25nsV6_DATA_Uncertainty_AK4PFchs.txt"),
                                 L2DATAPath = cms.string("data/Spring16_25nsV6_DATA_Uncertainty_AK4PFchs.txt"),
                                 L3DATAPath = cms.string("data/Spring16_25nsV6_DATA_Uncertainty_AK4PFchs.txt"),
                                 L2L3ResidualPath = cms.string("data/Spring16_25nsV6_DATA_Uncertainty_AK4PFchs.txt"), 
                                 JecFileUncMC =  cms.string("data/Spring16_25nsV6_MC_Uncertainty_AK4PFchs.txt"),
                                 JecFileUncData =  cms.string("data/Spring16_25nsV6_DATA_Uncertainty_AK4PFchs.txt"),                     
                                 MuonIDIsoSF = cms.string("data/Muon_IdIso0p15_eff.root"),
                                 MuonTriggerSF  = cms.string("data/Muon_IsoMu22_eff.root"), 
		)


process.p = cms.Path(process.tauhist)
