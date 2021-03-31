# This file skims the data and saves the output to ./tmp
# Do not combine files across runs, otherwise you may get inconsistent TTree structures!
# Doing things file by file is the safest way to avoid this problem, and comes at almost
# no extra cost.
# You can copy and paste json sources directly from https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/

import os
import ROOT


#path = '/pnfs/iihe/cms/store/user/rgoldouz/RunIIFall15MiniAODv2/GJets_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'
#path = '/pnfs/iihe/cms/store/user/rgoldouz/RunIIFall15MiniAODv2/ZToEE_NNPDF30_13TeV-powheg_M_50_120'
#path = '/pnfs/iihe/cms/store/user/rgoldouz/RunIIFall15MiniAODv2/ZToEE_NNPDF30_13TeV-powheg_M_120_200'
#path = '/pnfs/iihe/cms/store/user/rgoldouz/RunIIFall15MiniAODv2/ZToEE_NNPDF30_13TeV-powheg_M_200_400'
#path = '/pnfs/iihe/cms/store/user/rgoldouz/RunIIFall15MiniAODv2/ZToEE_NNPDF30_13TeV-powheg_M_400_800'
#path = '/pnfs/iihe/cms/store/user/rgoldouz/RunIIFall15MiniAODv2/ZToEE_NNPDF30_13TeV-powheg_M_800_1400'
#path = '/pnfs/iihe/cms/store/user/rgoldouz/RunIIFall15MiniAODv2/ZToEE_NNPDF30_13TeV-powheg_M_1400_2300'
#path = '/pnfs/iihe/cms/store/user/rgoldouz/RunIIFall15MiniAODv2/ZToEE_NNPDF30_13TeV-powheg_M_2300_3500'
#path = '/pnfs/iihe/cms/store/user/rgoldouz/RunIIFall15MiniAODv2/ZToEE_NNPDF30_13TeV-powheg_M_3500_4500'
#path = '/pnfs/iihe/cms/store/user/rgoldouz/RunIIFall15MiniAODv2/'
#path = '/pnfs/iihe/cms/store/user/rgoldouz/RunIIFall15MiniAODv2/'
#path = '/pnfs/iihe/cms/store/user/rgoldouz/RunIIFall15MiniAODv2/'
#path = '/pnfs/iihe/cms/store/user/rgoldouz/RunIIFall15MiniAODv2/'
#path = '/pnfs/iihe/cms/store/user/rgoldouz/RunIIFall15MiniAODv2/'

#path = '/pnfs/iihe/cms/store/user/rgoldouz/RunIIFall15MiniAODv2/GJets_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'
#path = '/pnfs/iihe/cms/store/user/rgoldouz/RunIIFall15MiniAODv2/ZToEE_NNPDF30_13TeV-powheg_M_50_120'
#path = '/pnfs/iihe/cms/store/user/rgoldouz/RunIIFall15MiniAODv2/ZToEE_NNPDF30_13TeV-powheg_M_120_200'
#path = '/pnfs/iihe/cms/store/user/rgoldouz/RunIIFall15MiniAODv2/ZToEE_NNPDF30_13TeV-powheg_M_200_400'
#path = '/pnfs/iihe/cms/store/user/rgoldouz/RunIIFall15MiniAODv2/ZToEE_NNPDF30_13TeV-powheg_M_400_800'
#path = '/pnfs/iihe/cms/store/user/rgoldouz/RunIIFall15MiniAODv2/ZToEE_NNPDF30_13TeV-powheg_M_800_1400'
#path = '/pnfs/iihe/cms/store/user/rgoldouz/RunIIFall15MiniAODv2/ZToEE_NNPDF30_13TeV-powheg_M_1400_2300'
#path = '/pnfs/iihe/cms/store/user/rgoldouz/RunIISpring16DR80/ZToEE_NNPDF30_13TeV-powheg_M_2300_3500'
#path = '/pnfs/iihe/cms/store/user/rgoldouz/RunIIFall15MiniAODv2/ZToEE_NNPDF30_13TeV-powheg_M_3500_4500'
#path = '/pnfs/iihe/cms/store/user/rgoldouz/RunIIFall15MiniAODv2/'
#path = '/pnfs/iihe/cms/store/user/rgoldouz/RunIIFall15MiniAODv2/'
#path = '/pnfs/iihe/cms/store/user/rgoldouz/RunIIFall15MiniAODv2/'
#path = '/pnfs/iihe/cms/store/user/rgoldouz/RunIIFall15MiniAODv2/'
#path = '/pnfs/iihe/cms/store/user/rgoldouz/RunIIFall15MiniAODv2/'
#path = '/pnfs/iihe/cms/store/user/wenxing/SingleElectron/crab_SingleElectron_Run2016B-PromptReco-v2_AOD_golden_0623/160623_122348/0001/'
#path =  'root://131.225.204.161:1094/' 

#path.append( '/pnfs/iihe/cms/store/user/rgoldouz/DYToEE_NNPDF30_13TeV-powheg-pythia8/crab_DYToEE_NNPDF30_13TeV-powheg-pythia8/160712_190134/0000')
#path.append( '/pnfs/iihe/cms/store/user/rgoldouz/DYToEE_NNPDF30_13TeV-powheg-pythia8/crab_DYToEE_NNPDF30_13TeV-powheg-pythia8/160712_190134/0001')
#path.append( '/pnfs/iihe/cms/store/user/rgoldouz/DYToEE_NNPDF30_13TeV-powheg-pythia8/crab_DYToEE_NNPDF30_13TeV-powheg-pythia8/160712_190134/0002')

#path.append( '/pnfs/iihe/cms/store/user/rgoldouz/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1/160712_190218/0000')
#path.append( '/pnfs/iihe/cms/store/user/rgoldouz/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1/160712_190218/0001')
#path.append( '/pnfs/iihe/cms/store/user/rgoldouz/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1/160712_190218/0002')

#path.append( '/pnfs/iihe/cms/store/user/rgoldouz/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/160712_190253/0000')
#path.append( '/pnfs/iihe/cms/store/user/rgoldouz/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/160712_190253/0001')
#path.append( '/pnfs/iihe/cms/store/user/rgoldouz/ZToEE_NNPDF30_13TeV-powheg_M_50_120/crab_ZToEE_NNPDF30_13TeV-powheg_M_50_120/160712_190042/0000')
#path.append( '/pnfs/iihe/cms/store/user/rgoldouz/ZToEE_NNPDF30_13TeV-powheg_M_120_200/crab_ZToEE_NNPDF30_13TeV-powheg_M_120_200/160712_185449/0000')
#path.append( '/pnfs/iihe/cms/store/user/rgoldouz/ZToEE_NNPDF30_13TeV-powheg_M_200_400/crab_ZToEE_NNPDF30_13TeV-powheg_M_200_400/160712_185522/0000')
#path.append( '/pnfs/iihe/cms/store/user/rgoldouz/ZToEE_NNPDF30_13TeV-powheg_M_400_800/crab_ZToEE_NNPDF30_13TeV-powheg_M_400_800/160712_185552/0000')
#path.append( '/pnfs/iihe/cms/store/user/rgoldouz/ZToEE_NNPDF30_13TeV-powheg_M_800_1400/crab_ZToEE_NNPDF30_13TeV-powheg_M_800_1400/160712_185625/0000')
#path.append( '/pnfs/iihe/cms/store/user/rgoldouz/ZToEE_NNPDF30_13TeV-powheg_M_1400_2300/crab_ZToEE_NNPDF30_13TeV-powheg_M_1400_2300/160712_185654/0000')
#path.append( '/pnfs/iihe/cms/store/user/rgoldouz/ZToEE_NNPDF30_13TeV-powheg_M_2300_3500/crab_ZToEE_NNPDF30_13TeV-powheg_M_2300_3500/160712_185726/0000')
#path.append( '/pnfs/iihe/cms/store/user/rgoldouz/ZToEE_NNPDF30_13TeV-powheg_M_3500_4500/crab_ZToEE_NNPDF30_13TeV-powheg_M_3500_4500/160712_185759/0000')
#path.append( '/pnfs/iihe/cms/store/user/rgoldouz/ZToEE_NNPDF30_13TeV-powheg_M_4500_6000/crab_ZToEE_NNPDF30_13TeV-powheg_M_4500_6000/160712_185838/0000')
#path.append( '/pnfs/iihe/cms/store/user/rgoldouz/ZToEE_NNPDF30_13TeV-powheg_M_6000_Inf/crab_ZToEE_NNPDF30_13TeV-powheg_M_6000_Inf/160712_185908/0000')

#path.append( '/pnfs/iihe/cms/store/user/rgoldouz/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/160731_122856/0000')
#path.append( '/pnfs/iihe/cms/store/user/rgoldouz/DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/160731_122954/0000')
#path.append('/pnfs/iihe/cms/store/user/rgoldouz/DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/160731_123019/0000')
#path.append('/pnfs/iihe/cms/store/user/rgoldouz/DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/160731_123058/0000')
#path.append('/pnfs/iihe/cms/store/user/rgoldouz/DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/160731_123137/0000')
#path.append('/pnfs/iihe/cms/store/user/rgoldouz/DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/160731_123209/0000')
#path.append('/pnfs/iihe/cms/store/user/rgoldouz/DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/160731_123234/0000')
#path.append('/pnfs/iihe/cms/store/user/rgoldouz/DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/160731_123258/0000')
path='root://131.225.204.161:1094//store/group/lpctau/amkalsi/MC_Run17/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-/171117_192554/0000/'

nEventsraw = 0
neventsweight = 0
nEventsStored = 0
nEventsiihe = 0
for a in path:
    print a
    filenames = os.listdir(a)
    for fname in filenames:
        filename = a + '/' + fname
        print fname 
        f = ROOT.TFile.Open(filename)
#        if not f:
#            print 'rm -rf '+fname
        tree_in = f.Get('IIHEAnalysis')
#        tree_meta = f.Get('meta')
        nEventsiihe += tree_in.GetEntries()
#        tree_meta.GetEntry(0)    
#        print tree_meta.nEventsRaw
#        nEventsraw += tree_meta.nEventsRaw
#        nEventsStored += tree_meta.nEventsStored
#        neventsweight += tree_meta.mc_nEventsWeighted
        f.Close()
#print 'nEventsraw %d   '%(nEventsraw)
#print 'neventsweight %d   '%(neventsweight)
#print 'nEventsStored %d   '%(nEventsStored)
print 'nEventsiihe %d   '%(nEventsiihe)
