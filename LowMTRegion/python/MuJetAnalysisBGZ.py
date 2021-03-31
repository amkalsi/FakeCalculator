import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('python')

options.register('inputFilename', 'DYJets_inc_small', #HTauTauAnalysis_1_1_Sl2.root',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Input file name"
)


options.register('outFilename', 'Output_default.root', #HTauTauAnalysis_1_1_Sl2.root',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Output file name"
)

options.register('sample', 'DYinc', #HTauTauAnalysis_1_1_Sl2.root',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Name of sample you are running on"
)

options.parseArguments()

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

#process.TFileService=cms.Service("TFileService",fileName=cms.string("Dec31_OUT_"+options.outFilename))

process.TFileService=cms.Service("TFileService",fileName=cms.string(options.outFilename))



process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.source = cms.Source("EmptySource")

process.tauhist = cms.EDAnalyzer('LowMTRegion',
                                 isolationTau = cms.string("LooseDB"),
                                 InputFile = cms.string(options.inputFilename),
				 DYSample = cms.string(options.sample),
#                                 MuonIDIsoSF = cms.string("/user/amkalsi/CMSSW_8_0_17/src/Plots/LowMTRegion/python/data/TriggerEff.root"), #2017Feb-A-Popov_TriggerSF_Run2016All_v1.root"),
                                 MuonIDIsoSF = cms.string("/user/amkalsi/CMSSW_8_0_17/src/Plots/LowMTRegion/python/data/egammaEffi.txt_EGM2D_2016.root"),
				 MuonTriggerSF  = cms.string("/user/amkalsi/CMSSW_8_0_17/src/Plots/LowMTRegion/python//data/Muon_IsoMu22_eff.root"),
				 isOS = cms.bool(False),
                                 isSS = cms.bool(False),
				 FakeRateSSfile = cms.string("/user/amkalsi/CMSSW_8_0_17/src/Plots/LowMTRegion/python/data/fakerate_SSMtLow.root"),
				 FakeRateDYJetsfile = cms.string("/user/amkalsi/CMSSW_8_0_17/src/Plots/LowMTRegion/python/data/fakerate.root"),
#                                 LooseIsoWPFile  = cms.string("/user/amkalsi/CMSSW_8_0_17/src/Plots/ETauAnalysis_2017/python/data/TauIDSF/TauID_SF_pt_MVAoldDM2017v2_2016.root")

                                 LooseIsoWPFile  = cms.string("/user/amkalsi/CMSSW_8_0_17/src/Plots/LowMTRegion/python/data/TauID_SF_pt_DeepTau2017v2p1VSjet_2016Legacy.root"),
                                 LowPtTES = cms.string("/user/amkalsi/CMSSW_8_0_17/src/Plots/LowMTRegion/python/data/TauES_dm_DeepTau2017v2p1VSjet_2016Legacy.root"),
                                 HighPtTES  = cms.string("/user/amkalsi/CMSSW_8_0_17/src/Plots/LowMTRegion/python/data/TauES_dm_DeepTau2017v2p1VSjet_2016Legacy_ptgt100.root"),
                                 FESDM = cms.string("/user/amkalsi/CMSSW_8_0_17/src/Plots/LowMTRegion/python/data/TauFES_eta-dm_DeepTau2017v2p1VSe_2016Legacy.root"),
		)


process.p = cms.Path(process.tauhist)
