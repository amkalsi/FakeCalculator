//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Jul  1 23:06:25 2018 by ROOT version 6.06/01
// from TTree IIHEAnalysis/IIHEAnalysis
// found on file: /pnfs/iihe/cms/store/user/amkalsi/Moriond2018_New_2016_June13_2018/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYJetsToLL_M-50_mcatnlo_ext2/180613_214737/0002/out_file_2008.root
//////////////////////////////////////////////////////////

#ifndef IIHEAnalysis_h
#define IIHEAnalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
using namespace std;

class IIHEAnalysis {
	public :
		TTree          *_tree;   //!pointer to the analyzed TTree or TChain
		Int_t           fCurrent; //!current Tree number in a TChain

		// Fixed size dimensions of array or collections stored in the TTree if any.
   ULong64_t       ev_event;
   ULong64_t       ev_run;
   ULong64_t       ev_luminosityBlock;
   UInt_t          ev_time;
   UInt_t          ev_time_unixTime;
   UInt_t          ev_time_microsecondOffset;
   Float_t         ev_fixedGridRhoAll;
   Float_t         ev_fixedGridRhoFastjetAll;
   Float_t         ev_fixedGridRhoFastjetAllCalo;
   Float_t         ev_fixedGridRhoFastjetCentralCalo;
   Float_t         ev_fixedGridRhoFastjetCentralChargedPileUp;
   Float_t         ev_fixedGridRhoFastjetCentralNeutral;
   vector<float>   *LHE_Pt;
   vector<float>   *LHE_Eta;
   vector<float>   *LHE_Phi;
   vector<float>   *LHE_E;
   vector<int>     *LHE_pdgid;
   vector<int>     *LHE_status;
   Float_t         LHE_weight_nominal;
   vector<float>   *LHE_weight_sys;
   UInt_t          mc_n;
   Float_t         mc_weight;
   Float_t         mc_w_sign;
   Int_t           mc_id_first;
   Int_t           mc_id_second;
   Float_t         mc_x_first;
   Float_t         mc_x_second;
   Float_t         mc_xPDF_first;
   Float_t         mc_xPDF_second;
   Float_t         mc_scalePDF;
   vector<int>     *mc_index;
   vector<int>     *mc_pdgId;
   vector<int>     *mc_charge;
   vector<int>     *mc_status;
   vector<int>     *mc_status_flags;
   vector<int>     *mc_status_tau_flags;
   vector<int>     *mc_tau_charge;
   vector<int>     *mc_tau_pdgId;
   vector<int>     *mc_tau_decay;
   vector<int>     *mc_tau_had_status;
   vector<int>     *mc_tau_had_charge;
   vector<int>     *mc_tau_had_pdgId;
   vector<float>   *mc_mass;
   vector<float>   *mc_px;
   vector<float>   *mc_py;
   vector<float>   *mc_pz;
   vector<float>   *mc_pt;
   vector<float>   *mc_eta;
   vector<float>   *mc_phi;
   vector<float>   *mc_energy;
   vector<float>   *mc_tau_pt;
   vector<float>   *mc_tau_eta;
   vector<float>   *mc_tau_phi;
   vector<float>   *mc_tau_energy;
   vector<float>   *mc_tau_had_pt;
   vector<float>   *mc_tau_had_eta;
   vector<float>   *mc_tau_had_phi;
   vector<float>   *mc_tau_had_energy;
   vector<unsigned int> *mc_numberOfDaughters;
   vector<unsigned int> *mc_numberOfMothers;
   vector<vector<int> > *mc_mother_index;
   vector<vector<int> > *mc_mother_pdgId;
   vector<vector<float> > *mc_mother_px;
   vector<vector<float> > *mc_mother_py;
   vector<vector<float> > *mc_mother_pz;
   vector<vector<float> > *mc_mother_pt;
   vector<vector<float> > *mc_mother_eta;
   vector<vector<float> > *mc_mother_phi;
   vector<vector<float> > *mc_mother_energy;
   vector<vector<float> > *mc_mother_mass;
   Int_t           mc_trueNumInteractions;
   Int_t           mc_PU_NumInteractions;
   vector<float>   *genjet_pt;
   vector<float>   *genjet_eta;
   vector<float>   *genjet_phi;
   vector<float>   *genjet_energy;
   vector<float>   *gen_weight_sys;
   UInt_t          pv_n;
   vector<float>   *pv_x;
   vector<float>   *pv_y;
   vector<float>   *pv_z;
   vector<float>   *pv_ndof;
   vector<float>   *pv_normalizedChi2;
   vector<int>     *pv_isValid;
   vector<int>     *pv_isFake;
   UInt_t          gsf_n;
   vector<int>     *gsf_classification;
   vector<float>   *gsf_ecalTrkEnergyPostCorr;
   vector<float>   *gsf_ecalEnergyPostCorr;
   vector<float>   *gsf_energy;
   vector<float>   *gsf_p;
   vector<float>   *gsf_pt;
   vector<float>   *gsf_et;
   vector<float>   *gsf_scE1x5;
   vector<float>   *gsf_scE5x5;
   vector<float>   *gsf_scE2x5Max;
   vector<float>   *gsf_full5x5_e5x5;
   vector<float>   *gsf_full5x5_e1x5;
   vector<float>   *gsf_full5x5_e2x5Max;
   vector<float>   *gsf_full5x5_sigmaIetaIeta;
   vector<float>   *gsf_full5x5_hcalOverEcal;
   vector<float>   *gsf_eta;
   vector<float>   *gsf_phi;
   vector<float>   *gsf_theta;
   vector<float>   *gsf_px;
   vector<float>   *gsf_py;
   vector<float>   *gsf_pz;
   vector<float>   *gsf_caloEnergy;
   vector<float>   *gsf_deltaEtaSuperClusterTrackAtVtx;
   vector<float>   *gsf_deltaPhiSuperClusterTrackAtVtx;
   vector<float>   *gsf_hadronicOverEm;
   vector<float>   *gsf_hcalDepth1OverEcal;
   vector<float>   *gsf_hcalDepth2OverEcal;
   vector<float>   *gsf_dr03TkSumPt;
   vector<float>   *gsf_heepTrkPtIso;
   vector<float>   *gsf_dr03EcalRecHitSumEt;
   vector<float>   *gsf_dr03HcalDepth1TowerSumEt;
   vector<float>   *gsf_dr03HcalDepth2TowerSumEt;
   vector<int>     *gsf_charge;
   vector<float>   *gsf_sigmaIetaIeta;
   vector<int>     *gsf_ecaldrivenSeed;
   vector<int>     *gsf_trackerdrivenSeed;
   vector<int>     *gsf_isEB;
   vector<int>     *gsf_isEE;
   vector<int>     *gsf_passConversionVeto;
   vector<int>     *gsf_VID_cutBasedElectronID_Fall17_94X_V1_loose;
   vector<int>     *gsf_VID_cutBasedElectronID_Fall17_94X_V1_medium;
   vector<int>     *gsf_VID_cutBasedElectronID_Fall17_94X_V1_tight;
   vector<int>     *gsf_VID_cutBasedElectronID_Fall17_94X_V1_veto;
   vector<int>     *gsf_VID_cutBasedElectronID_Summer16_80X_V1_loose;
   vector<int>     *gsf_VID_cutBasedElectronID_Summer16_80X_V1_medium;
   vector<int>     *gsf_VID_cutBasedElectronID_Summer16_80X_V1_tight;
   vector<int>     *gsf_VID_cutBasedElectronID_Summer16_80X_V1_veto;
   vector<int>     *gsf_VID_heepElectronID_HEEPV70;
   vector<int>     *gsf_VID_mvaEleID_Fall17_iso_V1_wp80;
   vector<int>     *gsf_VID_mvaEleID_Fall17_iso_V1_wp90;
   vector<int>     *gsf_VID_mvaEleID_Fall17_iso_V1_wpLoose;
   vector<int>     *gsf_VID_mvaEleID_Fall17_noIso_V1_wp80;
   vector<int>     *gsf_VID_mvaEleID_Fall17_noIso_V1_wp90;
   vector<int>     *gsf_VID_mvaEleID_Fall17_noIso_V1_wpLoose;
   vector<int>     *gsf_VID_mvaEleID_Spring16_GeneralPurpose_V1_wp80;
   vector<int>     *gsf_VID_mvaEleID_Spring16_GeneralPurpose_V1_wp90;
   vector<int>     *gsf_VID_mvaEleID_Spring16_HZZ_V1_wpLoose;
   vector<float>   *gsf_deltaEtaSeedClusterTrackAtCalo;
   vector<float>   *gsf_deltaPhiSeedClusterTrackAtCalo;
   vector<float>   *gsf_ecalEnergy;
   vector<float>   *gsf_eSuperClusterOverP;
   vector<float>   *gsf_dxy;
   vector<float>   *gsf_dxy_beamSpot;
   vector<float>   *gsf_dxy_firstPVtx;
   vector<float>   *gsf_dxyError;
   vector<float>   *gsf_dz;
   vector<float>   *gsf_dz_beamSpot;
   vector<float>   *gsf_dz_firstPVtx;
   vector<float>   *gsf_dzError;
   vector<float>   *gsf_vz;
   vector<int>     *gsf_numberOfValidHits;
   vector<int>     *gsf_nLostInnerHits;
   vector<int>     *gsf_nLostOuterHits;
   vector<int>     *gsf_convFlags;
   vector<float>   *gsf_convDist;
   vector<float>   *gsf_convDcot;
   vector<float>   *gsf_convRadius;
   vector<float>   *gsf_fBrem;
   vector<float>   *gsf_e1x5;
   vector<float>   *gsf_e2x5Max;
   vector<float>   *gsf_e5x5;
   vector<float>   *gsf_r9;
   vector<float>   *gsf_deltaEtaSeedClusterTrackAtVtx;
   vector<float>   *gsf_effArea;
   vector<float>   *gsf_sumChargedHadronPt;
   vector<float>   *gsf_sumNeutralHadronEt;
   vector<float>   *gsf_sumPhotonEt;
   vector<float>   *gsf_sumPUPt;
   vector<float>   *gsf_ooEmooP;
   vector<vector<int> > *gsf_hitsinfo;
   vector<float>   *gsf_pixelMatch_dPhi1;
   vector<float>   *gsf_pixelMatch_dPhi2;
   vector<float>   *gsf_pixelMatch_dRz1;
   vector<float>   *gsf_pixelMatch_dRz2;
   vector<int>     *gsf_pixelMatch_subDetector1;
   vector<int>     *gsf_pixelMatch_subDetector2;
   vector<float>   *gsf_mc_bestDR;
   vector<int>     *gsf_mc_index;
   vector<float>   *gsf_mc_ERatio;
   vector<float>   *gsf_sc_energy;
   vector<float>   *gsf_sc_seed_eta;
   vector<float>   *gsf_sc_eta;
   vector<float>   *gsf_sc_etacorr;
   vector<float>   *gsf_sc_theta;
   vector<float>   *gsf_sc_thetacorr;
   vector<float>   *gsf_sc_et;
   vector<float>   *gsf_sc_phi;
   vector<float>   *gsf_sc_px;
   vector<float>   *gsf_sc_py;
   vector<float>   *gsf_sc_pz;
   vector<float>   *gsf_sc_x;
   vector<float>   *gsf_sc_y;
   vector<float>   *gsf_sc_z;
   vector<float>   *gsf_sc_phiWidth;
   vector<float>   *gsf_sc_etaWidth;
   vector<int>     *gsf_sc_seed_rawId;
   vector<int>     *gsf_sc_seed_ieta;
   vector<int>     *gsf_sc_seed_iphi;
   vector<int>     *gsf_sc_seed_kHasSwitchToGain6;
   vector<int>     *gsf_sc_seed_kHasSwitchToGain1;
   vector<float>   *gsf_swissCross;
   vector<float>   *gsf_sc_rawEnergy;
   vector<float>   *gsf_sc_preshowerEnergy;
   vector<float>   *gsf_sc_lazyTools_e2x5Right;
   vector<float>   *gsf_sc_lazyTools_e2x5Left;
   vector<float>   *gsf_sc_lazyTools_e2x5Top;
   vector<float>   *gsf_sc_lazyTools_e2x5Bottom;
   vector<float>   *gsf_sc_lazyTools_eMax;
   vector<float>   *gsf_sc_lazyTools_e2nd;
   vector<float>   *gsf_sc_lazyTools_eRight;
   vector<float>   *gsf_sc_lazyTools_eLeft;
   vector<float>   *gsf_sc_lazyTools_eTop;
   vector<float>   *gsf_sc_lazyTools_eBottom;
   vector<float>   *gsf_sc_lazyTools_e2x2;
   vector<float>   *gsf_sc_lazyTools_e3x3;
   vector<float>   *gsf_sc_lazyTools_e4x4;
   vector<float>   *gsf_sc_lazyTools_e5x5;
   vector<float>   *gsf_sc_lazyTools_e1x3;
   vector<float>   *gsf_sc_lazyTools_e3x1;
   vector<float>   *gsf_sc_lazyTools_e1x5;
   vector<float>   *gsf_sc_lazyTools_e5x1;
   vector<float>   *gsf_sc_lazyTools_eshitsixix;
   vector<float>   *gsf_sc_lazyTools_eshitsiyiy;
   vector<float>   *gsf_sc_lazyTools_eseffsixix;
   vector<float>   *gsf_sc_lazyTools_eseffsiyiy;
   vector<float>   *gsf_sc_lazyTools_eseffsirir;
   vector<float>   *gsf_sc_lazyTools_BasicClusterSeedTime;
   Int_t           EHits_isSaturated;
   vector<int>     *EBHits_rawId;
   vector<int>     *EBHits_iRechit;
   vector<int>     *EBHits_ieta;
   vector<int>     *EBHits_iphi;
   vector<int>     *EBHits_RecoFlag;
   vector<int>     *EBHits_kSaturated;
   vector<int>     *EBHits_kLeadingEdgeRecovered;
   vector<int>     *EBHits_kNeighboursRecovered;
   vector<int>     *EBHits_kWeird;
   vector<float>   *EBHits_energy;
   vector<int>     *EEHits_rawId;
   vector<int>     *EEHits_iRechit;
   vector<int>     *EEHits_ieta;
   vector<int>     *EEHits_iphi;
   vector<int>     *EEHits_RecoFlag;
   vector<int>     *EEHits_kSaturated;
   vector<int>     *EEHits_kLeadingEdgeRecovered;
   vector<int>     *EEHits_kNeighboursRecovered;
   vector<int>     *EEHits_kWeird;
   vector<float>   *EEHits_energy;
   UInt_t          mu_n;
   vector<float>   *mu_gt_qoverp;
   vector<int>     *mu_gt_charge;
   vector<float>   *mu_gt_pt;
   vector<float>   *mu_gt_eta;
   vector<float>   *mu_gt_phi;
   vector<float>   *mu_gt_p;
   vector<float>   *mu_gt_px;
   vector<float>   *mu_gt_py;
   vector<float>   *mu_gt_pz;
   vector<float>   *mu_gt_theta;
   vector<float>   *mu_gt_lambda;
   vector<float>   *mu_gt_d0;
   vector<float>   *mu_gt_dz;
   vector<float>   *mu_gt_dz_beamspot;
   vector<float>   *mu_gt_dz_firstPVtx;
   vector<float>   *mu_gt_dxy;
   vector<float>   *mu_gt_dxy_beamspot;
   vector<float>   *mu_gt_dxy_firstPVtx;
   vector<float>   *mu_gt_dsz;
   vector<float>   *mu_gt_vx;
   vector<float>   *mu_gt_vy;
   vector<float>   *mu_gt_vz;
   vector<float>   *mu_gt_qoverpError;
   vector<float>   *mu_gt_ptError;
   vector<float>   *mu_gt_thetaError;
   vector<float>   *mu_gt_lambdaError;
   vector<float>   *mu_gt_phiError;
   vector<float>   *mu_gt_dxyError;
   vector<float>   *mu_gt_d0Error;
   vector<float>   *mu_gt_dszError;
   vector<float>   *mu_gt_dzError;
   vector<float>   *mu_gt_etaError;
   vector<float>   *mu_gt_chi2;
   vector<float>   *mu_gt_ndof;
   vector<float>   *mu_gt_normalizedChi2;
   vector<float>   *mu_ot_qoverp;
   vector<int>     *mu_ot_charge;
   vector<float>   *mu_ot_pt;
   vector<float>   *mu_ot_eta;
   vector<float>   *mu_ot_phi;
   vector<float>   *mu_ot_p;
   vector<float>   *mu_ot_px;
   vector<float>   *mu_ot_py;
   vector<float>   *mu_ot_pz;
   vector<float>   *mu_ot_theta;
   vector<float>   *mu_ot_lambda;
   vector<float>   *mu_ot_d0;
   vector<float>   *mu_ot_dz;
   vector<float>   *mu_ot_dz_beamspot;
   vector<float>   *mu_ot_dz_firstPVtx;
   vector<float>   *mu_ot_dxy;
   vector<float>   *mu_ot_dxy_beamspot;
   vector<float>   *mu_ot_dxy_firstPVtx;
   vector<float>   *mu_ot_dsz;
   vector<float>   *mu_ot_vx;
   vector<float>   *mu_ot_vy;
   vector<float>   *mu_ot_vz;
   vector<float>   *mu_ot_qoverpError;
   vector<float>   *mu_ot_ptError;
   vector<float>   *mu_ot_thetaError;
   vector<float>   *mu_ot_lambdaError;
   vector<float>   *mu_ot_phiError;
   vector<float>   *mu_ot_dxyError;
   vector<float>   *mu_ot_d0Error;
   vector<float>   *mu_ot_dszError;
   vector<float>   *mu_ot_dzError;
   vector<float>   *mu_ot_etaError;
   vector<float>   *mu_ot_chi2;
   vector<float>   *mu_ot_ndof;
   vector<float>   *mu_ot_normalizedChi2;
   vector<float>   *mu_it_qoverp;
   vector<int>     *mu_it_charge;
   vector<float>   *mu_it_pt;
   vector<float>   *mu_it_eta;
   vector<float>   *mu_it_phi;
   vector<float>   *mu_it_p;
   vector<float>   *mu_it_px;
   vector<float>   *mu_it_py;
   vector<float>   *mu_it_pz;
   vector<float>   *mu_it_theta;
   vector<float>   *mu_it_lambda;
   vector<float>   *mu_it_d0;
   vector<float>   *mu_it_dz;
   vector<float>   *mu_it_dz_beamspot;
   vector<float>   *mu_it_dz_firstPVtx;
   vector<float>   *mu_it_dxy;
   vector<float>   *mu_it_dxy_beamspot;
   vector<float>   *mu_it_dxy_firstPVtx;
   vector<float>   *mu_it_dsz;
   vector<float>   *mu_it_vx;
   vector<float>   *mu_it_vy;
   vector<float>   *mu_it_vz;
   vector<float>   *mu_it_qoverpError;
   vector<float>   *mu_it_ptError;
   vector<float>   *mu_it_thetaError;
   vector<float>   *mu_it_lambdaError;
   vector<float>   *mu_it_phiError;
   vector<float>   *mu_it_dxyError;
   vector<float>   *mu_it_d0Error;
   vector<float>   *mu_it_dszError;
   vector<float>   *mu_it_dzError;
   vector<float>   *mu_it_etaError;
   vector<float>   *mu_it_chi2;
   vector<float>   *mu_it_ndof;
   vector<float>   *mu_it_normalizedChi2;
   vector<float>   *mu_ibt_qoverp;
   vector<int>     *mu_ibt_charge;
   vector<float>   *mu_ibt_pt;
   vector<float>   *mu_ibt_eta;
   vector<float>   *mu_ibt_phi;
   vector<float>   *mu_ibt_p;
   vector<float>   *mu_ibt_px;
   vector<float>   *mu_ibt_py;
   vector<float>   *mu_ibt_pz;
   vector<float>   *mu_ibt_theta;
   vector<float>   *mu_ibt_lambda;
   vector<float>   *mu_ibt_d0;
   vector<float>   *mu_ibt_dz;
   vector<float>   *mu_ibt_dz_beamspot;
   vector<float>   *mu_ibt_dz_firstPVtx;
   vector<float>   *mu_ibt_dxy;
   vector<float>   *mu_ibt_dxy_beamspot;
   vector<float>   *mu_ibt_dxy_firstPVtx;
   vector<float>   *mu_ibt_dsz;
   vector<float>   *mu_ibt_vx;
   vector<float>   *mu_ibt_vy;
   vector<float>   *mu_ibt_vz;
   vector<float>   *mu_ibt_qoverpError;
   vector<float>   *mu_ibt_ptError;
   vector<float>   *mu_ibt_thetaError;
   vector<float>   *mu_ibt_lambdaError;
   vector<float>   *mu_ibt_phiError;
   vector<float>   *mu_ibt_dxyError;
   vector<float>   *mu_ibt_d0Error;
   vector<float>   *mu_ibt_dszError;
   vector<float>   *mu_ibt_dzError;
   vector<float>   *mu_ibt_etaError;
   vector<float>   *mu_ibt_chi2;
   vector<float>   *mu_ibt_ndof;
   vector<float>   *mu_ibt_normalizedChi2;
   vector<int>     *mu_isGlobalMuon;
   vector<int>     *mu_isStandAloneMuon;
   vector<int>     *mu_isTrackerMuon;
   vector<int>     *mu_isPFMuon;
   vector<int>     *mu_isPFIsolationValid;
   vector<int>     *mu_isGoodMuonTMLastStationLoose;
   vector<int>     *mu_isGoodMuonTMLastStationTight;
   vector<int>     *mu_isGoodMuonTM2DCompatibilityLoose;
   vector<int>     *mu_isGoodMuonTM2DCompatibilityTight;
   vector<int>     *mu_isGoodMuonTMOneStationLoose;
   vector<int>     *mu_isGoodMuonTMOneStationTight;
   vector<int>     *mu_isGoodMuonTMLastStationOptimizedLowPtLoose;
   vector<int>     *mu_isGoodMuonTMLastStationOptimizedLowPtTight;
   vector<int>     *mu_isTightMuon;
   vector<int>     *mu_isMediumMuon;
   vector<int>     *mu_isLooseMuon;
   vector<int>     *mu_isSoftMuon;
   vector<int>     *mu_isHighPtMuonOld;
   vector<int>     *mu_isHighPtMuon;
   vector<int>     *mu_isTrackerHighPtMuon;
   vector<int>     *mu_CutBasedIdLoose;
   vector<int>     *mu_CutBasedIdMedium;
   vector<int>     *mu_CutBasedIdMediumPrompt;
   vector<int>     *mu_CutBasedIdTight;
   vector<int>     *mu_CutBasedIdGlobalHighPt;
   vector<int>     *mu_CutBasedIdTrkHighPt;
   vector<int>     *mu_PFIsoVeryLoose;
   vector<int>     *mu_PFIsoLoose;
   vector<int>     *mu_PFIsoMedium;
   vector<int>     *mu_PFIsoTight;
   vector<int>     *mu_PFIsoVeryTight;
   vector<int>     *mu_PFIsoVeryVeryTight;
   vector<int>     *mu_TkIsoLoose;
   vector<int>     *mu_TkIsoTight;
   vector<int>     *mu_SoftCutBasedId;
   vector<int>     *mu_SoftMvaId;
   vector<int>     *mu_MvaLoose;
   vector<int>     *mu_MvaMedium;
   vector<int>     *mu_MvaTight;
   vector<int>     *mu_MiniIsoLoose;
   vector<int>     *mu_MiniIsoMedium;
   vector<int>     *mu_MiniIsoTight;
   vector<int>     *mu_MiniIsoVeryTight;
   vector<int>     *mu_TriggerIdLoose;
   vector<int>     *mu_InTimeMuon;
   vector<int>     *mu_MultiIsoLoose;
   vector<int>     *mu_MultiIsoMedium;
   vector<int>     *mu_numberOfMatchedStations;
   vector<int>     *mu_numberOfValidPixelHits;
   vector<int>     *mu_trackerLayersWithMeasurement;
   vector<int>     *mu_numberOfValidMuonHits;
   vector<int>     *mu_pixelLayersWithMeasurement;
   vector<float>   *mu_innerTrack_validFraction;
   vector<float>   *mu_combinedQuality_trkKink;
   vector<float>   *mu_combinedQuality_chi2LocalPosition;
   vector<float>   *mu_segmentCompatibility;
   vector<float>   *mu_dB;
   vector<float>   *mu_pt_default;
   vector<float>   *mu_isolationR03_sumPt;
   vector<float>   *mu_isolationR03_trackerVetoPt;
   vector<float>   *mu_isolationR03_emEt;
   vector<float>   *mu_isolationR03_emVetoEt;
   vector<float>   *mu_isolationR03_hadEt;
   vector<float>   *mu_isolationR03_hadVetoEt;
   vector<float>   *mu_isolationR05_sumPt;
   vector<float>   *mu_isolationR05_trackerVetoPt;
   vector<float>   *mu_isolationR05_emEt;
   vector<float>   *mu_isolationR05_emVetoEt;
   vector<float>   *mu_isolationR05_hadEt;
   vector<float>   *mu_isolationR05_hadVetoEt;
   vector<float>   *mu_pfIsolationR03_sumChargedHadronPt;
   vector<float>   *mu_pfIsolationR03_sumNeutralHadronEt;
   vector<float>   *mu_pfIsolationR03_sumChargedParticlePt;
   vector<float>   *mu_pfIsolationR03_sumPhotonEt;
   vector<float>   *mu_pfIsolationR03_sumNeutralHadronEtHighThreshold;
   vector<float>   *mu_pfIsolationR03_sumPhotonEtHighThreshold;
   vector<float>   *mu_pfIsolationR03_sumPUPt;
   vector<float>   *mu_pfIsolationR04_sumChargedHadronPt;
   vector<float>   *mu_pfIsolationR04_sumNeutralHadronEt;
   vector<float>   *mu_pfIsolationR04_sumChargedParticlePt;
   vector<float>   *mu_pfIsolationR04_sumPhotonEt;
   vector<float>   *mu_pfIsolationR04_sumNeutralHadronEtHighThreshold;
   vector<float>   *mu_pfIsolationR04_sumPhotonEtHighThreshold;
   vector<float>   *mu_pfIsolationR04_sumPUPt;
   vector<float>   *mu_pfIsoDbCorrected03;
   vector<float>   *mu_pfIsoDbCorrected04;
   vector<float>   *mu_isoTrackerBased03;
   vector<float>   *mu_mc_bestDR;
   vector<int>     *mu_mc_index;
   vector<float>   *mu_mc_ERatio;
   UInt_t          jet_n;
   vector<float>   *jet_px;
   vector<float>   *jet_py;
   vector<float>   *jet_pz;
   vector<float>   *jet_pt;
   vector<float>   *jet_eta;
   vector<float>   *jet_theta;
   vector<float>   *jet_phi;
   vector<float>   *jet_energy;
   vector<float>   *jet_mass;
   vector<float>   *jet_chargedEmEnergyFraction;
   vector<float>   *jet_neutralHadronEnergyFraction;
   vector<float>   *jet_neutralEmEnergyFraction;
   vector<float>   *jet_chargedHadronEnergyFraction;
   vector<float>   *jet_muonEnergyFraction;
   vector<int>     *jet_chargedMultiplicity;
   vector<int>     *jet_neutralMultiplicity;
   vector<int>     *jet_partonFlavour;
   vector<int>     *jet_hadronFlavour;
   vector<float>   *jet_CSVv2;
   vector<float>   *jet_CvsL;
   vector<float>   *jet_CvsB;
   vector<float>   *jet_MVA2BJets;
   vector<float>   *jet_CvsB_DeepJet_charm_tagger;
   vector<float>   *jet_CvsL_DeepJet_charm_tagger;
   vector<float>   *jet_CvsB_DeepCSV_charm_tagger;
   vector<float>   *jet_CvsL_DeepCSV_charm_tagger;
   vector<float>   *jet_DeepJet;
   vector<float>   *jet_DeepCSV;
   vector<int>     *jet_isJetIDLoose_2016;
   vector<int>     *jet_isJetIDTight_2016;
   vector<int>     *jet_isJetIDTightLepVeto_2016;
   vector<int>     *jet_isJetID_2017;
   vector<int>     *jet_isJetIDLepVeto_2017;
   vector<int>     *jet_isJetID_2018;
   vector<int>     *jet_isJetIDLepVeto_2018;
   vector<float>   *jet_Smeared_pt;
   vector<float>   *jet_SmearedJetResUp_pt;
   vector<float>   *jet_SmearedJetResDown_pt;
   vector<float>   *jet_EnUp_pt;
   vector<float>   *jet_EnDown_pt;
   Float_t         BtagSF_Deepcsv_loose;
   Float_t         BtagSF_Deepcsv_Up_loose;
   Float_t         BtagSF_Deepcsv_Down_loose;
   Float_t         BtagSF_Deepcsv_medium;
   Float_t         BtagSF_Deepcsv_Up_medium;
   Float_t         BtagSF_Deepcsv_Down_medium;
   Float_t         BtagSF_Deepcsv_tight;
   Float_t         BtagSF_Deepcsv_Up_tight;
   Float_t         BtagSF_Deepcsv_Down_tight;
   Float_t         MET_nominal_Pt;
   Float_t         MET_nominal_Px;
   Float_t         MET_nominal_Py;
   Float_t         MET_nominal_phi;
   Float_t         MET_nominal_significance;
   Float_t         MET_Pt;
   Float_t         MET_Px;
   Float_t         MET_Py;
   Float_t         MET_phi;
   Float_t         MET_significance;
   Float_t         MET_T1_Pt;
   Float_t         MET_T1_Px;
   Float_t         MET_T1_Py;
   Float_t         MET_T1_phi;
   Float_t         MET_T1_significance;
   Float_t         MET_gen_pt;
   Float_t         MET_gen_phi;
   vector<float>   *MET_Type1Unc_Px;
   vector<float>   *MET_Type1Unc_Py;
   vector<float>   *MET_Type1Unc_Pt;
   Float_t         MET_T1JetEnDown_Pt;
   Float_t         MET_T1JetEnDown_Px;
   Float_t         MET_T1JetEnDown_Py;
   Float_t         MET_T1JetEnDown_phi;
   Float_t         MET_T1JetEnDown_significance;
   Float_t         MET_T1JetEnUp_Pt;
   Float_t         MET_T1JetEnUp_Px;
   Float_t         MET_T1JetEnUp_Py;
   Float_t         MET_T1JetEnUp_phi;
   Float_t         MET_T1JetEnUp_significance;
   Float_t         MET_T1Smear_Pt;
   Float_t         MET_T1Smear_Px;
   Float_t         MET_T1Smear_Py;
   Float_t         MET_T1Smear_phi;
   Float_t         MET_T1Smear_significance;
   Float_t         MET_T1SmearJetEnDown_Pt;
   Float_t         MET_T1SmearJetEnDown_Px;
   Float_t         MET_T1SmearJetEnDown_Py;
   Float_t         MET_T1SmearJetEnDown_phi;
   Float_t         MET_T1SmearJetEnDown_significance;
   Float_t         MET_T1SmearJetEnUp_Pt;
   Float_t         MET_T1SmearJetEnUp_Px;
   Float_t         MET_T1SmearJetEnUp_Py;
   Float_t         MET_T1SmearJetEnUp_phi;
   Float_t         MET_T1SmearJetEnUp_significance;
   Float_t         MET_T1SmearJetResDown_Pt;
   Float_t         MET_T1SmearJetResDown_Px;
   Float_t         MET_T1SmearJetResDown_Py;
   Float_t         MET_T1SmearJetResDown_phi;
   Float_t         MET_T1SmearJetResDown_significance;
   Float_t         MET_T1SmearJetResUp_Pt;
   Float_t         MET_T1SmearJetResUp_Px;
   Float_t         MET_T1SmearJetResUp_Py;
   Float_t         MET_T1SmearJetResUp_phi;
   Float_t         MET_T1SmearJetResUp_significance;
   Float_t         MET_T1Txy_Pt;
   Float_t         MET_T1Txy_Px;
   Float_t         MET_T1Txy_Py;
   Float_t         MET_T1Txy_phi;
   Float_t         MET_T1Txy_significance;
   Float_t         MET_FinalCollection_Pt;
   Float_t         MET_FinalCollection_Px;
   Float_t         MET_FinalCollection_Py;
   Float_t         MET_FinalCollection_phi;
   Float_t         MET_FinalCollection_significance;
   UInt_t          tau_n;
   vector<float>   *tau_px;
   vector<float>   *tau_py;
   vector<float>   *tau_pz;
   vector<float>   *tau_pt;
   vector<float>   *tau_eta;
   vector<float>   *tau_theta;
   vector<float>   *tau_phi;
   vector<float>   *tau_energy;
   vector<float>   *tau_mass;
   vector<float>   *tau_dxy;
   vector<float>   *tau_dxy_error;
   vector<float>   *tau_ptLeadChargedCand;
   vector<float>   *tau_byIsolationMVArun2v1DBoldDMwLTraw;
   vector<float>   *tau_byCombinedIsolationDeltaBetaCorrRaw3Hits;
   vector<float>   *tau_byIsolationMVArun2v1DBnewDMwLTraw;
   vector<float>   *tau_byIsolationMVArun2v1PWoldDMwLTraw;
   vector<float>   *tau_byIsolationMVArun2v1PWnewDMwLTraw;
   vector<float>   *tau_byIsolationMVArun2v1DBdR03oldDMwLTraw;
   vector<float>   *tau_byIsolationMVArun2v1PWdR03oldDMwLTraw;
   vector<float>   *tau_againstElectronMVA6Raw;
   vector<float>   *tau_againstElectronMVA6category;
   vector<float>   *tau_mc_bestDR;
   vector<float>   *tau_mc_ERatio;
   vector<float>   *tau_decayModeFinding;
   vector<float>   *tau_decayModeFindingNewDMs;
   vector<float>   *tau_againstMuonLoose3;
   vector<float>   *tau_againstMuonTight3;
   vector<float>   *tau_byLooseCombinedIsolationDeltaBetaCorr3Hits;
   vector<float>   *tau_byMediumCombinedIsolationDeltaBetaCorr3Hits;
   vector<float>   *tau_byTightCombinedIsolationDeltaBetaCorr3Hits;
   vector<float>   *tau_byVLooseIsolationMVArun2v1DBoldDMwLT;
   vector<float>   *tau_byLooseIsolationMVArun2v1DBoldDMwLT;
   vector<float>   *tau_byMediumIsolationMVArun2v1DBoldDMwLT;
   vector<float>   *tau_byTightIsolationMVArun2v1DBoldDMwLT;
   vector<float>   *tau_byVTightIsolationMVArun2v1DBoldDMwLT;
   vector<float>   *tau_byVVTightIsolationMVArun2v1DBoldDMwLT;
   vector<float>   *tau_byVLooseIsolationMVArun2v1DBnewDMwLT;
   vector<float>   *tau_byLooseIsolationMVArun2v1DBnewDMwLT;
   vector<float>   *tau_byMediumIsolationMVArun2v1DBnewDMwLT;
   vector<float>   *tau_byTightIsolationMVArun2v1DBnewDMwLT;
   vector<float>   *tau_byVTightIsolationMVArun2v1DBnewDMwLT;
   vector<float>   *tau_byVVTightIsolationMVArun2v1DBnewDMwLT;
   vector<float>   *tau_byVLooseIsolationMVArun2v1PWoldDMwLT;
   vector<float>   *tau_byLooseIsolationMVArun2v1PWoldDMwLT;
   vector<float>   *tau_byMediumIsolationMVArun2v1PWoldDMwLT;
   vector<float>   *tau_byTightIsolationMVArun2v1PWoldDMwLT;
   vector<float>   *tau_byVTightIsolationMVArun2v1PWoldDMwLT;
   vector<float>   *tau_byVVTightIsolationMVArun2v1PWoldDMwLT;
   vector<float>   *tau_byVLooseIsolationMVArun2v1PWnewDMwLT;
   vector<float>   *tau_byLooseIsolationMVArun2v1PWnewDMwLT;
   vector<float>   *tau_byMediumIsolationMVArun2v1PWnewDMwLT;
   vector<float>   *tau_byTightIsolationMVArun2v1PWnewDMwLT;
   vector<float>   *tau_byVTightIsolationMVArun2v1PWnewDMwLT;
   vector<float>   *tau_byVVTightIsolationMVArun2v1PWnewDMwLT;
   vector<float>   *tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT;
   vector<float>   *tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT;
   vector<float>   *tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT;
   vector<float>   *tau_byTightIsolationMVArun2v1DBdR03oldDMwLT;
   vector<float>   *tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT;
   vector<float>   *tau_byVVTightIsolationMVArun2v1DBdR03oldDMwLT;
   vector<float>   *tau_byVLooseIsolationMVArun2v1PWdR03oldDMwLT;
   vector<float>   *tau_byLooseIsolationMVArun2v1PWdR03oldDMwLT;
   vector<float>   *tau_byMediumIsolationMVArun2v1PWdR03oldDMwLT;
   vector<float>   *tau_byTightIsolationMVArun2v1PWdR03oldDMwLT;
   vector<float>   *tau_byVTightIsolationMVArun2v1PWdR03oldDMwLT;
   vector<float>   *tau_byVVTightIsolationMVArun2v1PWdR03oldDMwLT;
   vector<float>   *tau_againstElectronVLooseMVA6;
   vector<float>   *tau_againstElectronLooseMVA6;
   vector<float>   *tau_againstElectronMediumMVA6;
   vector<float>   *tau_againstElectronTightMVA6;
   vector<float>   *tau_againstElectronVTightMVA6;
   vector<float>   *tau_chargedIsoPtSum;
   vector<float>   *tau_neutralIsoPtSum;
   vector<float>   *tau_puCorrPtSum;
   vector<float>   *tau_footprintCorrection;
   vector<float>   *tau_neutralIsoPtSumWeight;
   vector<float>   *tau_photonPtSumOutsideSignalCone;
   vector<float>   *tau_byPhotonPtSumOutsideSignalCone;
   vector<float>   *tau_footprintCorrectiondR03;
   vector<float>   *tau_chargedIsoPtSumdR03;
   vector<float>   *tau_neutralIsoPtSumWeightdR03;
   vector<float>   *tau_neutralIsoPtSumdR03;
   vector<float>   *tau_photonPtSumOutsideSignalConedR03;
   vector<float>   *tau_PFChargedHadIso;
   vector<float>   *tau_PFNeutralHadIso;
   vector<float>   *tau_PFPhotonIso;
   vector<float>   *tau_leadChargedParticlePt;
   vector<float>   *tau_trackRefPt;
   vector<float>   *tau_lead_dxy;
   vector<float>   *tau_lead_dz;
   vector<float>   *tau_dxy_Sig;
   vector<float>   *tau_flightLengthSig;
   vector<float>   *tau_ip3d;
   vector<float>   *tau_ip3d_Sig;
   vector<float>   *tau_decayDistX;
   vector<float>   *tau_decayDistY;
   vector<float>   *tau_decayDistZ;
   vector<float>   *tau_decayDistMag;
   vector<float>   *tau_nPhoton;
   vector<float>   *tau_ptWeightedDetaStrip;
   vector<float>   *tau_ptWeightedDphiStrip;
   vector<float>   *tau_ptWeightedDrSignal;
   vector<float>   *tau_ptWeightedDrIsolation;
   vector<float>   *tau_leadingTrackChi2;
   vector<float>   *tau_eRatio;
   vector<float>   *tau_gjAngleDiff;
   vector<unsigned int> *tau_numberOfIsolationChargedHadrCands;
   vector<unsigned int> *tau_numberOfSignalChargedHadrCands;
   vector<unsigned int> *tau_numNeutralHadronsSignalCone;
   vector<unsigned int> *tau_numPhotonsSignalCone;
   vector<unsigned int> *tau_numParticlesSignalCone;
   vector<unsigned int> *tau_numChargedParticlesIsoCone;
   vector<unsigned int> *tau_numNeutralHadronsIsoCone;
   vector<unsigned int> *tau_numPhotonsIsoCone;
   vector<unsigned int> *tau_numParticlesIsoCone;
   vector<int>     *tau_mc_index;
   vector<int>     *tau_decayMode;
   vector<int>     *tau_charge;
   vector<int>     *tau_isPFTau;
   vector<int>     *tau_hasSecondaryVertex;
   vector<int>     *tau_leadChargedHadrAvailable;
   vector<float>   *tau_byIsolationMVArun2017v1DBoldDMwLTraw2017;
   vector<float>   *tau_byIsolationMVArun2017v2DBnewDMwLTraw2017;
   vector<float>   *tau_byIsolationMVArun2017v2DBoldDMdR0p3wLTraw2017;
   vector<float>   *tau_byIsolationMVArun2017v2DBoldDMwLTraw2017;
   vector<float>   *tau_byIsolationMVArun2v1DBnewDMwLTraw2016;
   vector<float>   *tau_byIsolationMVArun2v1DBoldDMwLTraw2016;
   vector<float>   *tau_byVVLooseIsolationMVArun2017v1DBoldDMwLT2017;
   vector<float>   *tau_byVLooseIsolationMVArun2017v1DBoldDMwLT2017;
   vector<float>   *tau_byLooseIsolationMVArun2017v1DBoldDMwLT2017;
   vector<float>   *tau_byMediumIsolationMVArun2017v1DBoldDMwLT2017;
   vector<float>   *tau_byTightIsolationMVArun2017v1DBoldDMwLT2017;
   vector<float>   *tau_byVTightIsolationMVArun2017v1DBoldDMwLT2017;
   vector<float>   *tau_byVVTightIsolationMVArun2017v1DBoldDMwLT2017;
   vector<float>   *tau_byVVLooseIsolationMVArun2017v2DBnewDMwLT2017;
   vector<float>   *tau_byVLooseIsolationMVArun2017v2DBnewDMwLT2017;
   vector<float>   *tau_byLooseIsolationMVArun2017v2DBnewDMwLT2017;
   vector<float>   *tau_byMediumIsolationMVArun2017v2DBnewDMwLT2017;
   vector<float>   *tau_byTightIsolationMVArun2017v2DBnewDMwLT2017;
   vector<float>   *tau_byVTightIsolationMVArun2017v2DBnewDMwLT2017;
   vector<float>   *tau_byVVTightIsolationMVArun2017v2DBnewDMwLT2017;
   vector<float>   *tau_byVVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017;
   vector<float>   *tau_byVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017;
   vector<float>   *tau_byLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017;
   vector<float>   *tau_byMediumIsolationMVArun2017v2DBoldDMdR0p3wLT2017;
   vector<float>   *tau_byTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017;
   vector<float>   *tau_byVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017;
   vector<float>   *tau_byVVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017;
   vector<float>   *tau_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017;
   vector<float>   *tau_byVLooseIsolationMVArun2017v2DBoldDMwLT2017;
   vector<float>   *tau_byLooseIsolationMVArun2017v2DBoldDMwLT2017;
   vector<float>   *tau_byMediumIsolationMVArun2017v2DBoldDMwLT2017;
   vector<float>   *tau_byTightIsolationMVArun2017v2DBoldDMwLT2017;
   vector<float>   *tau_byVTightIsolationMVArun2017v2DBoldDMwLT2017;
   vector<float>   *tau_byVVTightIsolationMVArun2017v2DBoldDMwLT2017;
   vector<float>   *tau_byVLooseIsolationMVArun2v1DBnewDMwLT2016;
   vector<float>   *tau_byLooseIsolationMVArun2v1DBnewDMwLT2016;
   vector<float>   *tau_byMediumIsolationMVArun2v1DBnewDMwLT2016;
   vector<float>   *tau_byTightIsolationMVArun2v1DBnewDMwLT2016;
   vector<float>   *tau_byVTightIsolationMVArun2v1DBnewDMwLT2016;
   vector<float>   *tau_byVVTightIsolationMVArun2v1DBnewDMwLT2016;
   vector<float>   *tau_byVLooseIsolationMVArun2v1DBoldDMwLT2016;
   vector<float>   *tau_byLooseIsolationMVArun2v1DBoldDMwLT2016;
   vector<float>   *tau_byMediumIsolationMVArun2v1DBoldDMwLT2016;
   vector<float>   *tau_byTightIsolationMVArun2v1DBoldDMwLT2016;
   vector<float>   *tau_byVTightIsolationMVArun2v1DBoldDMwLT2016;
   vector<float>   *tau_byVVTightIsolationMVArun2v1DBoldDMwLT2016;
   Int_t           trig_Flag_ecalBadCalibReduced;
   Int_t           trig_Flag_HBHENoiseFilter_accept;
   Int_t           trig_Flag_HBHENoiseIsoFilter_accept;
   Int_t           trig_Flag_CSCTightHaloFilter_accept;
   Int_t           trig_Flag_CSCTightHaloTrkMuUnvetoFilter_accept;
   Int_t           trig_Flag_CSCTightHalo2015Filter_accept;
   Int_t           trig_Flag_globalTightHalo2016Filter_accept;
   Int_t           trig_Flag_globalSuperTightHalo2016Filter_accept;
   Int_t           trig_Flag_HcalStripHaloFilter_accept;
   Int_t           trig_Flag_hcalLaserEventFilter_accept;
   Int_t           trig_Flag_EcalDeadCellTriggerPrimitiveFilter_accept;
   Int_t           trig_Flag_EcalDeadCellBoundaryEnergyFilter_accept;
   Int_t           trig_Flag_ecalBadCalibFilter_accept;
   Int_t           trig_Flag_goodVertices_accept;
   Int_t           trig_Flag_eeBadScFilter_accept;
   Int_t           trig_Flag_ecalLaserCorrFilter_accept;
   Int_t           trig_Flag_trkPOGFilters_accept;
   Int_t           trig_Flag_chargedHadronTrackResolutionFilter_accept;
   Int_t           trig_Flag_muonBadTrackFilter_accept;
   Int_t           trig_Flag_BadChargedCandidateFilter_accept;
   Int_t           trig_Flag_BadPFMuonFilter_accept;
   Int_t           trig_Flag_BadChargedCandidateSummer16Filter_accept;
   Int_t           trig_Flag_BadPFMuonSummer16Filter_accept;
   Int_t           trig_Flag_trkPOG_manystripclus53X_accept;
   Int_t           trig_Flag_trkPOG_toomanystripclus53X_accept;
   Int_t           trig_Flag_trkPOG_logErrorTooManyClusters_accept;
   Int_t           trig_Flag_METFilters_accept;
   Int_t           trig_raw2digi_step_accept;
   Int_t           trig_L1Reco_step_accept;
   Int_t           trig_reconstruction_step_accept;
   Int_t           trig_recosim_step_accept;
   Int_t           trig_eventinterpretaion_step_accept;
   Int_t           trig_HLT_DoubleEle25_CaloIdL_MW_accept;
   Int_t           trig_HLT_DoubleEle27_CaloIdL_MW_accept;
   Int_t           trig_HLT_DoubleEle33_CaloIdL_MW_accept;
   Int_t           trig_HLT_DoubleEle24_eta2p1_WPTight_Gsf_accept;
   Int_t           trig_HLT_Ele27_Ele37_CaloIdL_MW_accept;
   Int_t           trig_HLT_Mu27_Ele37_CaloIdL_MW_accept;
   Int_t           trig_HLT_Mu37_Ele27_CaloIdL_MW_accept;
   Int_t           trig_HLT_Mu37_TkMu27_accept;
   Int_t           trig_HLT_DoubleMu4_3_Bs_accept;
   Int_t           trig_HLT_Ele20_WPTight_Gsf_accept;
   Int_t           trig_HLT_Ele15_WPLoose_Gsf_accept;
   Int_t           trig_HLT_Ele17_WPLoose_Gsf_accept;
   Int_t           trig_HLT_Ele20_WPLoose_Gsf_accept;
   Int_t           trig_HLT_Ele20_eta2p1_WPLoose_Gsf_accept;
   Int_t           trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_accept;
   Int_t           trig_HLT_Ele27_WPTight_Gsf_accept;
   Int_t           trig_HLT_Ele28_WPTight_Gsf_accept;
   Int_t           trig_HLT_Ele30_WPTight_Gsf_accept;
   Int_t           trig_HLT_Ele32_WPTight_Gsf_accept;
   Int_t           trig_HLT_Ele35_WPTight_Gsf_accept;
   Int_t           trig_HLT_Ele35_WPTight_Gsf_L1EGMT_accept;
   Int_t           trig_HLT_Ele38_WPTight_Gsf_accept;
   Int_t           trig_HLT_Ele40_WPTight_Gsf_accept;
   Int_t           trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_accept;
   Int_t           trig_HLT_IsoMu20_accept;
   Int_t           trig_HLT_IsoMu24_accept;
   Int_t           trig_HLT_IsoMu24_eta2p1_accept;
   Int_t           trig_HLT_IsoMu27_accept;
   Int_t           trig_HLT_IsoMu30_accept;
   Int_t           trig_HLT_L1SingleMu18_accept;
   Int_t           trig_HLT_L1SingleMu25_accept;
   Int_t           trig_HLT_L2Mu10_accept;
   Int_t           trig_HLT_L2Mu10_NoVertex_NoBPTX3BX_accept;
   Int_t           trig_HLT_L2Mu10_NoVertex_NoBPTX_accept;
   Int_t           trig_HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX_accept;
   Int_t           trig_HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX_accept;
   Int_t           trig_HLT_L2Mu50_accept;
   Int_t           trig_HLT_L2Mu23NoVtx_2Cha_accept;
   Int_t           trig_HLT_L2Mu23NoVtx_2Cha_CosmicSeed_accept;
   Int_t           trig_HLT_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4_accept;
   Int_t           trig_HLT_DoubleL2Mu30NoVtx_2Cha_Eta2p4_accept;
   Int_t           trig_HLT_DoubleL2Mu50_accept;
   Int_t           trig_HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed_accept;
   Int_t           trig_HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed_NoL2Matched_accept;
   Int_t           trig_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_accept;
   Int_t           trig_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_NoL2Matched_accept;
   Int_t           trig_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4_accept;
   Int_t           trig_HLT_DoubleL2Mu23NoVtx_2Cha_accept;
   Int_t           trig_HLT_DoubleL2Mu23NoVtx_2Cha_NoL2Matched_accept;
   Int_t           trig_HLT_DoubleL2Mu25NoVtx_2Cha_accept;
   Int_t           trig_HLT_DoubleL2Mu25NoVtx_2Cha_NoL2Matched_accept;
   Int_t           trig_HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4_accept;
   Int_t           trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_accept;
   Int_t           trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_accept;
   Int_t           trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_accept;
   Int_t           trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_accept;
   Int_t           trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_accept;
   Int_t           trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8_accept;
   Int_t           trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_accept;
   Int_t           trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8_accept;
   Int_t           trig_HLT_Mu25_TkMu0_Onia_accept;
   Int_t           trig_HLT_Mu30_TkMu0_Psi_accept;
   Int_t           trig_HLT_Mu20_TkMu0_Phi_accept;
   Int_t           trig_HLT_Mu25_TkMu0_Phi_accept;
   Int_t           trig_HLT_Mu12_accept;
   Int_t           trig_HLT_Mu15_accept;
   Int_t           trig_HLT_Mu20_accept;
   Int_t           trig_HLT_Mu27_accept;
   Int_t           trig_HLT_Mu50_accept;
   Int_t           trig_HLT_Mu55_accept;
   Int_t           trig_HLT_OldMu100_accept;
   Int_t           trig_HLT_TkMu100_accept;
   Int_t           trig_HLT_Photon300_NoHE_accept;
   Int_t           trig_HLT_Mu8_TrkIsoVVL_accept;
   Int_t           trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept;
   Int_t           trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept;
   Int_t           trig_HLT_Mu17_TrkIsoVVL_accept;
   Int_t           trig_HLT_Mu19_TrkIsoVVL_accept;
   Int_t           trig_HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL_accept;
   Int_t           trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept;
   Int_t           trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_accept;
   Int_t           trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept;
   Int_t           trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept;
   Int_t           trig_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept;
   Int_t           trig_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept;
   Int_t           trig_HLT_Photon20_accept;
   Int_t           trig_HLT_Photon33_accept;
   Int_t           trig_HLT_Photon50_accept;
   Int_t           trig_HLT_Photon75_accept;
   Int_t           trig_HLT_Photon90_accept;
   Int_t           trig_HLT_Photon120_accept;
   Int_t           trig_HLT_Photon150_accept;
   Int_t           trig_HLT_Photon175_accept;
   Int_t           trig_HLT_Photon200_accept;
   Int_t           trig_HLT_Photon100EB_TightID_TightIso_accept;
   Int_t           trig_HLT_Photon110EB_TightID_TightIso_accept;
   Int_t           trig_HLT_Photon120EB_TightID_TightIso_accept;
   Int_t           trig_HLT_Photon100EBHE10_accept;
   Int_t           trig_HLT_Photon100EEHE10_accept;
   Int_t           trig_HLT_Photon100EE_TightID_TightIso_accept;
   Int_t           trig_HLT_Photon50_R9Id90_HE10_IsoM_accept;
   Int_t           trig_HLT_Photon75_R9Id90_HE10_IsoM_accept;
   Int_t           trig_HLT_Photon90_R9Id90_HE10_IsoM_accept;
   Int_t           trig_HLT_Photon120_R9Id90_HE10_IsoM_accept;
   Int_t           trig_HLT_Photon165_R9Id90_HE10_IsoM_accept;
   Int_t           trig_HLT_Photon35_TwoProngs35_accept;
   Int_t           trig_HLT_IsoMu24_TwoProngs35_accept;
   Int_t           trig_HLT_DoubleMu43NoFiltersNoVtx_accept;
   Int_t           trig_HLT_DoubleMu48NoFiltersNoVtx_accept;
   Int_t           trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_accept;
   Int_t           trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_accept;
   Int_t           trig_HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL_accept;
   Int_t           trig_HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL_accept;
   Int_t           trig_HLT_DoubleMu33NoFiltersNoVtxDisplaced_accept;
   Int_t           trig_HLT_DoubleMu40NoFiltersNoVtxDisplaced_accept;
   Int_t           trig_HLT_DoubleMu20_7_Mass0to30_L1_DM4_accept;
   Int_t           trig_HLT_DoubleMu20_7_Mass0to30_L1_DM4EG_accept;
   Int_t           trig_HLT_DoubleMu20_7_Mass0to30_Photon23_accept;
   Int_t           trig_HLT_Dimuon14_Phi_Barrel_Seagulls_accept;
   Int_t           trig_HLT_Dimuon24_Phi_noCorrL1_accept;
   Int_t           trig_HLT_DoubleIsoMu20_eta2p1_accept;
   Int_t           trig_HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx_accept;
   Int_t           trig_HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx_accept;
   Int_t           trig_HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx_accept;
   Int_t           trig_HLT_Mu8_accept;
   Int_t           trig_HLT_Mu17_accept;
   Int_t           trig_HLT_Mu19_accept;
   Int_t           trig_HLT_Mu17_Photon30_IsoCaloId_accept;
   Int_t           trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_accept;
   Int_t           trig_HLT_Ele135_CaloIdVT_GsfTrkIdT_accept;
   Int_t           trig_HLT_Ele145_CaloIdVT_GsfTrkIdT_accept;
   Int_t           trig_HLT_Ele200_CaloIdVT_GsfTrkIdT_accept;
   Int_t           trig_HLT_Ele250_CaloIdVT_GsfTrkIdT_accept;
   Int_t           trig_HLT_Ele300_CaloIdVT_GsfTrkIdT_accept;
   Int_t           trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_accept;
   Int_t           trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_accept;
   Int_t           trig_HLT_SinglePhoton10_Eta3p1ForPPRef_accept;
   Int_t           trig_HLT_SinglePhoton20_Eta3p1ForPPRef_accept;
   Int_t           trig_HLT_SinglePhoton30_Eta3p1ForPPRef_accept;
   Int_t           trig_HLT_Photon20_HoverELoose_accept;
   Int_t           trig_HLT_Photon30_HoverELoose_accept;
   Int_t           trig_HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142_accept;
   Int_t           trig_AlCa_RPCMuonNormalisation_accept;
   Int_t           trig_MC_DoubleEle5_CaloIdL_MW_accept;
   Int_t           trig_MC_Ele5_WPTight_Gsf_accept;
   Int_t           trig_MC_Ele15_Ele10_CaloIdL_TrackIdL_IsoVL_DZ_accept;
   Int_t           trig_MC_IsoMu_accept;
   Int_t           trig_MC_DoubleMu_TrkIsoVVL_DZ_accept;
   Int_t           trig_MC_DoubleMuNoFiltersNoVtx_accept;
   Int_t           trig_HLT_Mu18_Mu9_SameSign_accept;
   Int_t           trig_HLT_Mu18_Mu9_SameSign_DZ_accept;
   Int_t           trig_HLT_Mu18_Mu9_accept;
   Int_t           trig_HLT_Mu18_Mu9_DZ_accept;
   Int_t           trig_HLT_Mu20_Mu10_SameSign_accept;
   Int_t           trig_HLT_Mu20_Mu10_SameSign_DZ_accept;
   Int_t           trig_HLT_Mu20_Mu10_accept;
   Int_t           trig_HLT_Mu20_Mu10_DZ_accept;
   Int_t           trig_HLT_Mu23_Mu12_SameSign_accept;
   Int_t           trig_HLT_Mu23_Mu12_SameSign_DZ_accept;
   Int_t           trig_HLT_Mu23_Mu12_accept;
   Int_t           trig_HLT_Mu23_Mu12_DZ_accept;
   Int_t           trig_HLT_TrkMu6NoFiltersNoVtx_accept;
   Int_t           trig_HLT_TrkMu16NoFiltersNoVtx_accept;
   Int_t           trig_HLT_DoubleTrkMu_16_6_NoFiltersNoVtx_accept;

   // List of branches
   TBranch        *b_ev_event;   //!
   TBranch        *b_ev_run;   //!
   TBranch        *b_ev_luminosityBlock;   //!
   TBranch        *b_ev_time;   //!
   TBranch        *b_ev_time_unixTime;   //!
   TBranch        *b_ev_time_microsecondOffset;   //!
   TBranch        *b_ev_fixedGridRhoAll;   //!
   TBranch        *b_ev_fixedGridRhoFastjetAll;   //!
   TBranch        *b_ev_fixedGridRhoFastjetAllCalo;   //!
   TBranch        *b_ev_fixedGridRhoFastjetCentralCalo;   //!
   TBranch        *b_ev_fixedGridRhoFastjetCentralChargedPileUp;   //!
   TBranch        *b_ev_fixedGridRhoFastjetCentralNeutral;   //!
   TBranch        *b_LHE_Pt;   //!
   TBranch        *b_LHE_Eta;   //!
   TBranch        *b_LHE_Phi;   //!
   TBranch        *b_LHE_E;   //!
   TBranch        *b_LHE_pdgid;   //!
   TBranch        *b_LHE_status;   //!
   TBranch        *b_LHE_weight_nominal;   //!
   TBranch        *b_LHE_weight_sys;   //!
   TBranch        *b_mc_n;   //!
   TBranch        *b_mc_weight;   //!
   TBranch        *b_mc_w_sign;   //!
   TBranch        *b_mc_id_first;   //!
   TBranch        *b_mc_id_second;   //!
   TBranch        *b_mc_x_first;   //!
   TBranch        *b_mc_x_second;   //!
   TBranch        *b_mc_xPDF_first;   //!
   TBranch        *b_mc_xPDF_second;   //!
   TBranch        *b_mc_scalePDF;   //!
   TBranch        *b_mc_index;   //!
   TBranch        *b_mc_pdgId;   //!
   TBranch        *b_mc_charge;   //!
   TBranch        *b_mc_status;   //!
   TBranch        *b_mc_status_flags;   //!
   TBranch        *b_mc_status_tau_flags;   //!
   TBranch        *b_mc_tau_charge;   //!
   TBranch        *b_mc_tau_pdgId;   //!
   TBranch        *b_mc_tau_decay;   //!
   TBranch        *b_mc_tau_had_status;   //!
   TBranch        *b_mc_tau_had_charge;   //!
   TBranch        *b_mc_tau_had_pdgId;   //!
   TBranch        *b_mc_mass;   //!
   TBranch        *b_mc_px;   //!
   TBranch        *b_mc_py;   //!
   TBranch        *b_mc_pz;   //!
   TBranch        *b_mc_pt;   //!
   TBranch        *b_mc_eta;   //!
   TBranch        *b_mc_phi;   //!
   TBranch        *b_mc_energy;   //!
   TBranch        *b_mc_tau_pt;   //!
   TBranch        *b_mc_tau_eta;   //!
   TBranch        *b_mc_tau_phi;   //!
   TBranch        *b_mc_tau_energy;   //!
   TBranch        *b_mc_tau_had_pt;   //!
   TBranch        *b_mc_tau_had_eta;   //!
   TBranch        *b_mc_tau_had_phi;   //!
   TBranch        *b_mc_tau_had_energy;   //!
   TBranch        *b_mc_numberOfDaughters;   //!
   TBranch        *b_mc_numberOfMothers;   //!
   TBranch        *b_mc_mother_index;   //!
   TBranch        *b_mc_mother_pdgId;   //!
   TBranch        *b_mc_mother_px;   //!
   TBranch        *b_mc_mother_py;   //!
   TBranch        *b_mc_mother_pz;   //!
   TBranch        *b_mc_mother_pt;   //!
   TBranch        *b_mc_mother_eta;   //!
   TBranch        *b_mc_mother_phi;   //!
   TBranch        *b_mc_mother_energy;   //!
   TBranch        *b_mc_mother_mass;   //!
   TBranch        *b_mc_trueNumInteractions;   //!
   TBranch        *b_mc_PU_NumInteractions;   //!
   TBranch        *b_genjet_pt;   //!
   TBranch        *b_genjet_eta;   //!
   TBranch        *b_genjet_phi;   //!
   TBranch        *b_genjet_energy;   //!
   TBranch        *b_gen_weight_sys;   //!
   TBranch        *b_pv_n;   //!
   TBranch        *b_pv_x;   //!
   TBranch        *b_pv_y;   //!
   TBranch        *b_pv_z;   //!
   TBranch        *b_pv_ndof;   //!
   TBranch        *b_pv_normalizedChi2;   //!
   TBranch        *b_pv_isValid;   //!
   TBranch        *b_pv_isFake;   //!
   TBranch        *b_gsf_n;   //!
   TBranch        *b_gsf_classification;   //!
   TBranch        *b_gsf_ecalTrkEnergyPostCorr;   //!
   TBranch        *b_gsf_ecalEnergyPostCorr;   //!
   TBranch        *b_gsf_energy;   //!
   TBranch        *b_gsf_p;   //!
   TBranch        *b_gsf_pt;   //!
   TBranch        *b_gsf_et;   //!
   TBranch        *b_gsf_scE1x5;   //!
   TBranch        *b_gsf_scE5x5;   //!
   TBranch        *b_gsf_scE2x5Max;   //!
   TBranch        *b_gsf_full5x5_e5x5;   //!
   TBranch        *b_gsf_full5x5_e1x5;   //!
   TBranch        *b_gsf_full5x5_e2x5Max;   //!
   TBranch        *b_gsf_full5x5_sigmaIetaIeta;   //!
   TBranch        *b_gsf_full5x5_hcalOverEcal;   //!
   TBranch        *b_gsf_eta;   //!
   TBranch        *b_gsf_phi;   //!
   TBranch        *b_gsf_theta;   //!
   TBranch        *b_gsf_px;   //!
   TBranch        *b_gsf_py;   //!
   TBranch        *b_gsf_pz;   //!
   TBranch        *b_gsf_caloEnergy;   //!
   TBranch        *b_gsf_deltaEtaSuperClusterTrackAtVtx;   //!
   TBranch        *b_gsf_deltaPhiSuperClusterTrackAtVtx;   //!
   TBranch        *b_gsf_hadronicOverEm;   //!
   TBranch        *b_gsf_hcalDepth1OverEcal;   //!
   TBranch        *b_gsf_hcalDepth2OverEcal;   //!
   TBranch        *b_gsf_dr03TkSumPt;   //!
   TBranch        *b_gsf_heepTrkPtIso;   //!
   TBranch        *b_gsf_dr03EcalRecHitSumEt;   //!
   TBranch        *b_gsf_dr03HcalDepth1TowerSumEt;   //!
   TBranch        *b_gsf_dr03HcalDepth2TowerSumEt;   //!
   TBranch        *b_gsf_charge;   //!
   TBranch        *b_gsf_sigmaIetaIeta;   //!
   TBranch        *b_gsf_ecaldrivenSeed;   //!
   TBranch        *b_gsf_trackerdrivenSeed;   //!
   TBranch        *b_gsf_isEB;   //!
   TBranch        *b_gsf_isEE;   //!
   TBranch        *b_gsf_passConversionVeto;   //!
   TBranch        *b_gsf_VID_cutBasedElectronID_Fall17_94X_V1_loose;   //!
   TBranch        *b_gsf_VID_cutBasedElectronID_Fall17_94X_V1_medium;   //!
   TBranch        *b_gsf_VID_cutBasedElectronID_Fall17_94X_V1_tight;   //!
   TBranch        *b_gsf_VID_cutBasedElectronID_Fall17_94X_V1_veto;   //!
   TBranch        *b_gsf_VID_cutBasedElectronID_Summer16_80X_V1_loose;   //!
   TBranch        *b_gsf_VID_cutBasedElectronID_Summer16_80X_V1_medium;   //!
   TBranch        *b_gsf_VID_cutBasedElectronID_Summer16_80X_V1_tight;   //!
   TBranch        *b_gsf_VID_cutBasedElectronID_Summer16_80X_V1_veto;   //!
   TBranch        *b_gsf_VID_heepElectronID_HEEPV70;   //!
   TBranch        *b_gsf_VID_mvaEleID_Fall17_iso_V1_wp80;   //!
   TBranch        *b_gsf_VID_mvaEleID_Fall17_iso_V1_wp90;   //!
   TBranch        *b_gsf_VID_mvaEleID_Fall17_iso_V1_wpLoose;   //!
   TBranch        *b_gsf_VID_mvaEleID_Fall17_noIso_V1_wp80;   //!
   TBranch        *b_gsf_VID_mvaEleID_Fall17_noIso_V1_wp90;   //!
   TBranch        *b_gsf_VID_mvaEleID_Fall17_noIso_V1_wpLoose;   //!
   TBranch        *b_gsf_VID_mvaEleID_Spring16_GeneralPurpose_V1_wp80;   //!
   TBranch        *b_gsf_VID_mvaEleID_Spring16_GeneralPurpose_V1_wp90;   //!
   TBranch        *b_gsf_VID_mvaEleID_Spring16_HZZ_V1_wpLoose;   //!
   TBranch        *b_gsf_deltaEtaSeedClusterTrackAtCalo;   //!
   TBranch        *b_gsf_deltaPhiSeedClusterTrackAtCalo;   //!
   TBranch        *b_gsf_ecalEnergy;   //!
   TBranch        *b_gsf_eSuperClusterOverP;   //!
   TBranch        *b_gsf_dxy;   //!
   TBranch        *b_gsf_dxy_beamSpot;   //!
   TBranch        *b_gsf_dxy_firstPVtx;   //!
   TBranch        *b_gsf_dxyError;   //!
   TBranch        *b_gsf_dz;   //!
   TBranch        *b_gsf_dz_beamSpot;   //!
   TBranch        *b_gsf_dz_firstPVtx;   //!
   TBranch        *b_gsf_dzError;   //!
   TBranch        *b_gsf_vz;   //!
   TBranch        *b_gsf_numberOfValidHits;   //!
   TBranch        *b_gsf_nLostInnerHits;   //!
   TBranch        *b_gsf_nLostOuterHits;   //!
   TBranch        *b_gsf_convFlags;   //!
   TBranch        *b_gsf_convDist;   //!
   TBranch        *b_gsf_convDcot;   //!
   TBranch        *b_gsf_convRadius;   //!
   TBranch        *b_gsf_fBrem;   //!
   TBranch        *b_gsf_e1x5;   //!
   TBranch        *b_gsf_e2x5Max;   //!
   TBranch        *b_gsf_e5x5;   //!
   TBranch        *b_gsf_r9;   //!
   TBranch        *b_gsf_deltaEtaSeedClusterTrackAtVtx;   //!
   TBranch        *b_gsf_effArea;   //!
   TBranch        *b_gsf_sumChargedHadronPt;   //!
   TBranch        *b_gsf_sumNeutralHadronEt;   //!
   TBranch        *b_gsf_sumPhotonEt;   //!
   TBranch        *b_gsf_sumPUPt;   //!
   TBranch        *b_gsf_ooEmooP;   //!
   TBranch        *b_gsf_hitsinfo;   //!
   TBranch        *b_gsf_pixelMatch_dPhi1;   //!
   TBranch        *b_gsf_pixelMatch_dPhi2;   //!
   TBranch        *b_gsf_pixelMatch_dRz1;   //!
   TBranch        *b_gsf_pixelMatch_dRz2;   //!
   TBranch        *b_gsf_pixelMatch_subDetector1;   //!
   TBranch        *b_gsf_pixelMatch_subDetector2;   //!
   TBranch        *b_gsf_mc_bestDR;   //!
   TBranch        *b_gsf_mc_index;   //!
   TBranch        *b_gsf_mc_ERatio;   //!
   TBranch        *b_gsf_sc_energy;   //!
   TBranch        *b_gsf_sc_seed_eta;   //!
   TBranch        *b_gsf_sc_eta;   //!
   TBranch        *b_gsf_sc_etacorr;   //!
   TBranch        *b_gsf_sc_theta;   //!
   TBranch        *b_gsf_sc_thetacorr;   //!
   TBranch        *b_gsf_sc_et;   //!
   TBranch        *b_gsf_sc_phi;   //!
   TBranch        *b_gsf_sc_px;   //!
   TBranch        *b_gsf_sc_py;   //!
   TBranch        *b_gsf_sc_pz;   //!
   TBranch        *b_gsf_sc_x;   //!
   TBranch        *b_gsf_sc_y;   //!
   TBranch        *b_gsf_sc_z;   //!
   TBranch        *b_gsf_sc_phiWidth;   //!
   TBranch        *b_gsf_sc_etaWidth;   //!
   TBranch        *b_gsf_sc_seed_rawId;   //!
   TBranch        *b_gsf_sc_seed_ieta;   //!
   TBranch        *b_gsf_sc_seed_iphi;   //!
   TBranch        *b_gsf_sc_seed_kHasSwitchToGain6;   //!
   TBranch        *b_gsf_sc_seed_kHasSwitchToGain1;   //!
   TBranch        *b_gsf_swissCross;   //!
   TBranch        *b_gsf_sc_rawEnergy;   //!
   TBranch        *b_gsf_sc_preshowerEnergy;   //!
   TBranch        *b_gsf_sc_lazyTools_e2x5Right;   //!
   TBranch        *b_gsf_sc_lazyTools_e2x5Left;   //!
   TBranch        *b_gsf_sc_lazyTools_e2x5Top;   //!
   TBranch        *b_gsf_sc_lazyTools_e2x5Bottom;   //!
   TBranch        *b_gsf_sc_lazyTools_eMax;   //!
   TBranch        *b_gsf_sc_lazyTools_e2nd;   //!
   TBranch        *b_gsf_sc_lazyTools_eRight;   //!
   TBranch        *b_gsf_sc_lazyTools_eLeft;   //!
   TBranch        *b_gsf_sc_lazyTools_eTop;   //!
   TBranch        *b_gsf_sc_lazyTools_eBottom;   //!
   TBranch        *b_gsf_sc_lazyTools_e2x2;   //!
   TBranch        *b_gsf_sc_lazyTools_e3x3;   //!
   TBranch        *b_gsf_sc_lazyTools_e4x4;   //!
   TBranch        *b_gsf_sc_lazyTools_e5x5;   //!
   TBranch        *b_gsf_sc_lazyTools_e1x3;   //!
   TBranch        *b_gsf_sc_lazyTools_e3x1;   //!
   TBranch        *b_gsf_sc_lazyTools_e1x5;   //!
   TBranch        *b_gsf_sc_lazyTools_e5x1;   //!
   TBranch        *b_gsf_sc_lazyTools_eshitsixix;   //!
   TBranch        *b_gsf_sc_lazyTools_eshitsiyiy;   //!
   TBranch        *b_gsf_sc_lazyTools_eseffsixix;   //!
   TBranch        *b_gsf_sc_lazyTools_eseffsiyiy;   //!
   TBranch        *b_gsf_sc_lazyTools_eseffsirir;   //!
   TBranch        *b_gsf_sc_lazyTools_BasicClusterSeedTime;   //!
   TBranch        *b_EHits_isSaturated;   //!
   TBranch        *b_EBHits_rawId;   //!
   TBranch        *b_EBHits_iRechit;   //!
   TBranch        *b_EBHits_ieta;   //!
   TBranch        *b_EBHits_iphi;   //!
   TBranch        *b_EBHits_RecoFlag;   //!
   TBranch        *b_EBHits_kSaturated;   //!
   TBranch        *b_EBHits_kLeadingEdgeRecovered;   //!
   TBranch        *b_EBHits_kNeighboursRecovered;   //!
   TBranch        *b_EBHits_kWeird;   //!
   TBranch        *b_EBHits_energy;   //!
   TBranch        *b_EEHits_rawId;   //!
   TBranch        *b_EEHits_iRechit;   //!
   TBranch        *b_EEHits_ieta;   //!
   TBranch        *b_EEHits_iphi;   //!
   TBranch        *b_EEHits_RecoFlag;   //!
   TBranch        *b_EEHits_kSaturated;   //!
   TBranch        *b_EEHits_kLeadingEdgeRecovered;   //!
   TBranch        *b_EEHits_kNeighboursRecovered;   //!
   TBranch        *b_EEHits_kWeird;   //!
   TBranch        *b_EEHits_energy;   //!
   TBranch        *b_mu_n;   //!
   TBranch        *b_mu_gt_qoverp;   //!
   TBranch        *b_mu_gt_charge;   //!
   TBranch        *b_mu_gt_pt;   //!
   TBranch        *b_mu_gt_eta;   //!
   TBranch        *b_mu_gt_phi;   //!
   TBranch        *b_mu_gt_p;   //!
   TBranch        *b_mu_gt_px;   //!
   TBranch        *b_mu_gt_py;   //!
   TBranch        *b_mu_gt_pz;   //!
   TBranch        *b_mu_gt_theta;   //!
   TBranch        *b_mu_gt_lambda;   //!
   TBranch        *b_mu_gt_d0;   //!
   TBranch        *b_mu_gt_dz;   //!
   TBranch        *b_mu_gt_dz_beamspot;   //!
   TBranch        *b_mu_gt_dz_firstPVtx;   //!
   TBranch        *b_mu_gt_dxy;   //!
   TBranch        *b_mu_gt_dxy_beamspot;   //!
   TBranch        *b_mu_gt_dxy_firstPVtx;   //!
   TBranch        *b_mu_gt_dsz;   //!
   TBranch        *b_mu_gt_vx;   //!
   TBranch        *b_mu_gt_vy;   //!
   TBranch        *b_mu_gt_vz;   //!
   TBranch        *b_mu_gt_qoverpError;   //!
   TBranch        *b_mu_gt_ptError;   //!
   TBranch        *b_mu_gt_thetaError;   //!
   TBranch        *b_mu_gt_lambdaError;   //!
   TBranch        *b_mu_gt_phiError;   //!
   TBranch        *b_mu_gt_dxyError;   //!
   TBranch        *b_mu_gt_d0Error;   //!
   TBranch        *b_mu_gt_dszError;   //!
   TBranch        *b_mu_gt_dzError;   //!
   TBranch        *b_mu_gt_etaError;   //!
   TBranch        *b_mu_gt_chi2;   //!
   TBranch        *b_mu_gt_ndof;   //!
   TBranch        *b_mu_gt_normalizedChi2;   //!
   TBranch        *b_mu_ot_qoverp;   //!
   TBranch        *b_mu_ot_charge;   //!
   TBranch        *b_mu_ot_pt;   //!
   TBranch        *b_mu_ot_eta;   //!
   TBranch        *b_mu_ot_phi;   //!
   TBranch        *b_mu_ot_p;   //!
   TBranch        *b_mu_ot_px;   //!
   TBranch        *b_mu_ot_py;   //!
   TBranch        *b_mu_ot_pz;   //!
   TBranch        *b_mu_ot_theta;   //!
   TBranch        *b_mu_ot_lambda;   //!
   TBranch        *b_mu_ot_d0;   //!
   TBranch        *b_mu_ot_dz;   //!
   TBranch        *b_mu_ot_dz_beamspot;   //!
   TBranch        *b_mu_ot_dz_firstPVtx;   //!
   TBranch        *b_mu_ot_dxy;   //!
   TBranch        *b_mu_ot_dxy_beamspot;   //!
   TBranch        *b_mu_ot_dxy_firstPVtx;   //!
   TBranch        *b_mu_ot_dsz;   //!
   TBranch        *b_mu_ot_vx;   //!
   TBranch        *b_mu_ot_vy;   //!
   TBranch        *b_mu_ot_vz;   //!
   TBranch        *b_mu_ot_qoverpError;   //!
   TBranch        *b_mu_ot_ptError;   //!
   TBranch        *b_mu_ot_thetaError;   //!
   TBranch        *b_mu_ot_lambdaError;   //!
   TBranch        *b_mu_ot_phiError;   //!
   TBranch        *b_mu_ot_dxyError;   //!
   TBranch        *b_mu_ot_d0Error;   //!
   TBranch        *b_mu_ot_dszError;   //!
   TBranch        *b_mu_ot_dzError;   //!
   TBranch        *b_mu_ot_etaError;   //!
   TBranch        *b_mu_ot_chi2;   //!
   TBranch        *b_mu_ot_ndof;   //!
   TBranch        *b_mu_ot_normalizedChi2;   //!
   TBranch        *b_mu_it_qoverp;   //!
   TBranch        *b_mu_it_charge;   //!
   TBranch        *b_mu_it_pt;   //!
   TBranch        *b_mu_it_eta;   //!
   TBranch        *b_mu_it_phi;   //!
   TBranch        *b_mu_it_p;   //!
   TBranch        *b_mu_it_px;   //!
   TBranch        *b_mu_it_py;   //!
   TBranch        *b_mu_it_pz;   //!
   TBranch        *b_mu_it_theta;   //!
   TBranch        *b_mu_it_lambda;   //!
   TBranch        *b_mu_it_d0;   //!
   TBranch        *b_mu_it_dz;   //!
   TBranch        *b_mu_it_dz_beamspot;   //!
   TBranch        *b_mu_it_dz_firstPVtx;   //!
   TBranch        *b_mu_it_dxy;   //!
   TBranch        *b_mu_it_dxy_beamspot;   //!
   TBranch        *b_mu_it_dxy_firstPVtx;   //!
   TBranch        *b_mu_it_dsz;   //!
   TBranch        *b_mu_it_vx;   //!
   TBranch        *b_mu_it_vy;   //!
   TBranch        *b_mu_it_vz;   //!
   TBranch        *b_mu_it_qoverpError;   //!
   TBranch        *b_mu_it_ptError;   //!
   TBranch        *b_mu_it_thetaError;   //!
   TBranch        *b_mu_it_lambdaError;   //!
   TBranch        *b_mu_it_phiError;   //!
   TBranch        *b_mu_it_dxyError;   //!
   TBranch        *b_mu_it_d0Error;   //!
   TBranch        *b_mu_it_dszError;   //!
   TBranch        *b_mu_it_dzError;   //!
   TBranch        *b_mu_it_etaError;   //!
   TBranch        *b_mu_it_chi2;   //!
   TBranch        *b_mu_it_ndof;   //!
   TBranch        *b_mu_it_normalizedChi2;   //!
   TBranch        *b_mu_ibt_qoverp;   //!
   TBranch        *b_mu_ibt_charge;   //!
   TBranch        *b_mu_ibt_pt;   //!
   TBranch        *b_mu_ibt_eta;   //!
   TBranch        *b_mu_ibt_phi;   //!
   TBranch        *b_mu_ibt_p;   //!
   TBranch        *b_mu_ibt_px;   //!
   TBranch        *b_mu_ibt_py;   //!
   TBranch        *b_mu_ibt_pz;   //!
   TBranch        *b_mu_ibt_theta;   //!
   TBranch        *b_mu_ibt_lambda;   //!
   TBranch        *b_mu_ibt_d0;   //!
   TBranch        *b_mu_ibt_dz;   //!
   TBranch        *b_mu_ibt_dz_beamspot;   //!
   TBranch        *b_mu_ibt_dz_firstPVtx;   //!
   TBranch        *b_mu_ibt_dxy;   //!
   TBranch        *b_mu_ibt_dxy_beamspot;   //!
   TBranch        *b_mu_ibt_dxy_firstPVtx;   //!
   TBranch        *b_mu_ibt_dsz;   //!
   TBranch        *b_mu_ibt_vx;   //!
   TBranch        *b_mu_ibt_vy;   //!
   TBranch        *b_mu_ibt_vz;   //!
   TBranch        *b_mu_ibt_qoverpError;   //!
   TBranch        *b_mu_ibt_ptError;   //!
   TBranch        *b_mu_ibt_thetaError;   //!
   TBranch        *b_mu_ibt_lambdaError;   //!
   TBranch        *b_mu_ibt_phiError;   //!
   TBranch        *b_mu_ibt_dxyError;   //!
   TBranch        *b_mu_ibt_d0Error;   //!
   TBranch        *b_mu_ibt_dszError;   //!
   TBranch        *b_mu_ibt_dzError;   //!
   TBranch        *b_mu_ibt_etaError;   //!
   TBranch        *b_mu_ibt_chi2;   //!
   TBranch        *b_mu_ibt_ndof;   //!
   TBranch        *b_mu_ibt_normalizedChi2;   //!
   TBranch        *b_mu_isGlobalMuon;   //!
   TBranch        *b_mu_isStandAloneMuon;   //!
   TBranch        *b_mu_isTrackerMuon;   //!
   TBranch        *b_mu_isPFMuon;   //!
   TBranch        *b_mu_isPFIsolationValid;   //!
   TBranch        *b_mu_isGoodMuonTMLastStationLoose;   //!
   TBranch        *b_mu_isGoodMuonTMLastStationTight;   //!
   TBranch        *b_mu_isGoodMuonTM2DCompatibilityLoose;   //!
   TBranch        *b_mu_isGoodMuonTM2DCompatibilityTight;   //!
   TBranch        *b_mu_isGoodMuonTMOneStationLoose;   //!
   TBranch        *b_mu_isGoodMuonTMOneStationTight;   //!
   TBranch        *b_mu_isGoodMuonTMLastStationOptimizedLowPtLoose;   //!
   TBranch        *b_mu_isGoodMuonTMLastStationOptimizedLowPtTight;   //!
   TBranch        *b_mu_isTightMuon;   //!
   TBranch        *b_mu_isMediumMuon;   //!
   TBranch        *b_mu_isLooseMuon;   //!
   TBranch        *b_mu_isSoftMuon;   //!
   TBranch        *b_mu_isHighPtMuonOld;   //!
   TBranch        *b_mu_isHighPtMuon;   //!
   TBranch        *b_mu_isTrackerHighPtMuon;   //!
   TBranch        *b_mu_CutBasedIdLoose;   //!
   TBranch        *b_mu_CutBasedIdMedium;   //!
   TBranch        *b_mu_CutBasedIdMediumPrompt;   //!
   TBranch        *b_mu_CutBasedIdTight;   //!
   TBranch        *b_mu_CutBasedIdGlobalHighPt;   //!
   TBranch        *b_mu_CutBasedIdTrkHighPt;   //!
   TBranch        *b_mu_PFIsoVeryLoose;   //!
   TBranch        *b_mu_PFIsoLoose;   //!
   TBranch        *b_mu_PFIsoMedium;   //!
   TBranch        *b_mu_PFIsoTight;   //!
   TBranch        *b_mu_PFIsoVeryTight;   //!
   TBranch        *b_mu_PFIsoVeryVeryTight;   //!
   TBranch        *b_mu_TkIsoLoose;   //!
   TBranch        *b_mu_TkIsoTight;   //!
   TBranch        *b_mu_SoftCutBasedId;   //!
   TBranch        *b_mu_SoftMvaId;   //!
   TBranch        *b_mu_MvaLoose;   //!
   TBranch        *b_mu_MvaMedium;   //!
   TBranch        *b_mu_MvaTight;   //!
   TBranch        *b_mu_MiniIsoLoose;   //!
   TBranch        *b_mu_MiniIsoMedium;   //!
   TBranch        *b_mu_MiniIsoTight;   //!
   TBranch        *b_mu_MiniIsoVeryTight;   //!
   TBranch        *b_mu_TriggerIdLoose;   //!
   TBranch        *b_mu_InTimeMuon;   //!
   TBranch        *b_mu_MultiIsoLoose;   //!
   TBranch        *b_mu_MultiIsoMedium;   //!
   TBranch        *b_mu_numberOfMatchedStations;   //!
   TBranch        *b_mu_numberOfValidPixelHits;   //!
   TBranch        *b_mu_trackerLayersWithMeasurement;   //!
   TBranch        *b_mu_numberOfValidMuonHits;   //!
   TBranch        *b_mu_pixelLayersWithMeasurement;   //!
   TBranch        *b_mu_innerTrack_validFraction;   //!
   TBranch        *b_mu_combinedQuality_trkKink;   //!
   TBranch        *b_mu_combinedQuality_chi2LocalPosition;   //!
   TBranch        *b_mu_segmentCompatibility;   //!
   TBranch        *b_mu_dB;   //!
   TBranch        *b_mu_pt_default;   //!
   TBranch        *b_mu_isolationR03_sumPt;   //!
   TBranch        *b_mu_isolationR03_trackerVetoPt;   //!
   TBranch        *b_mu_isolationR03_emEt;   //!
   TBranch        *b_mu_isolationR03_emVetoEt;   //!
   TBranch        *b_mu_isolationR03_hadEt;   //!
   TBranch        *b_mu_isolationR03_hadVetoEt;   //!
   TBranch        *b_mu_isolationR05_sumPt;   //!
   TBranch        *b_mu_isolationR05_trackerVetoPt;   //!
   TBranch        *b_mu_isolationR05_emEt;   //!
   TBranch        *b_mu_isolationR05_emVetoEt;   //!
   TBranch        *b_mu_isolationR05_hadEt;   //!
   TBranch        *b_mu_isolationR05_hadVetoEt;   //!
   TBranch        *b_mu_pfIsolationR03_sumChargedHadronPt;   //!
   TBranch        *b_mu_pfIsolationR03_sumNeutralHadronEt;   //!
   TBranch        *b_mu_pfIsolationR03_sumChargedParticlePt;   //!
   TBranch        *b_mu_pfIsolationR03_sumPhotonEt;   //!
   TBranch        *b_mu_pfIsolationR03_sumNeutralHadronEtHighThreshold;   //!
   TBranch        *b_mu_pfIsolationR03_sumPhotonEtHighThreshold;   //!
   TBranch        *b_mu_pfIsolationR03_sumPUPt;   //!
   TBranch        *b_mu_pfIsolationR04_sumChargedHadronPt;   //!
   TBranch        *b_mu_pfIsolationR04_sumNeutralHadronEt;   //!
   TBranch        *b_mu_pfIsolationR04_sumChargedParticlePt;   //!
   TBranch        *b_mu_pfIsolationR04_sumPhotonEt;   //!
   TBranch        *b_mu_pfIsolationR04_sumNeutralHadronEtHighThreshold;   //!
   TBranch        *b_mu_pfIsolationR04_sumPhotonEtHighThreshold;   //!
   TBranch        *b_mu_pfIsolationR04_sumPUPt;   //!
   TBranch        *b_mu_pfIsoDbCorrected03;   //!
   TBranch        *b_mu_pfIsoDbCorrected04;   //!
   TBranch        *b_mu_isoTrackerBased03;   //!
   TBranch        *b_mu_mc_bestDR;   //!
   TBranch        *b_mu_mc_index;   //!
   TBranch        *b_mu_mc_ERatio;   //!
   TBranch        *b_jet_n;   //!
   TBranch        *b_jet_px;   //!
   TBranch        *b_jet_py;   //!
   TBranch        *b_jet_pz;   //!
   TBranch        *b_jet_pt;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_theta;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_jet_energy;   //!
   TBranch        *b_jet_mass;   //!
   TBranch        *b_jet_chargedEmEnergyFraction;   //!
   TBranch        *b_jet_neutralHadronEnergyFraction;   //!
   TBranch        *b_jet_neutralEmEnergyFraction;   //!
   TBranch        *b_jet_chargedHadronEnergyFraction;   //!
   TBranch        *b_jet_muonEnergyFraction;   //!
   TBranch        *b_jet_chargedMultiplicity;   //!
   TBranch        *b_jet_neutralMultiplicity;   //!
   TBranch        *b_jet_partonFlavour;   //!
   TBranch        *b_jet_hadronFlavour;   //!
   TBranch        *b_jet_CSVv2;   //!
   TBranch        *b_jet_CvsL;   //!
   TBranch        *b_jet_CvsB;   //!
   TBranch        *b_jet_MVA2BJets;   //!
   TBranch        *b_jet_CvsB_DeepJet_charm_tagger;   //!
   TBranch        *b_jet_CvsL_DeepJet_charm_tagger;   //!
   TBranch        *b_jet_CvsB_DeepCSV_charm_tagger;   //!
   TBranch        *b_jet_CvsL_DeepCSV_charm_tagger;   //!
   TBranch        *b_jet_DeepJet;   //!
   TBranch        *b_jet_DeepCSV;   //!
   TBranch        *b_jet_isJetIDLoose_2016;   //!
   TBranch        *b_jet_isJetIDTight_2016;   //!
   TBranch        *b_jet_isJetIDTightLepVeto_2016;   //!
   TBranch        *b_jet_isJetID_2017;   //!
   TBranch        *b_jet_isJetIDLepVeto_2017;   //!
   TBranch        *b_jet_isJetID_2018;   //!
   TBranch        *b_jet_isJetIDLepVeto_2018;   //!
   TBranch        *b_jet_Smeared_pt;   //!
   TBranch        *b_jet_SmearedJetResUp_pt;   //!
   TBranch        *b_jet_SmearedJetResDown_pt;   //!
   TBranch        *b_jet_EnUp_pt;   //!
   TBranch        *b_jet_EnDown_pt;   //!
   TBranch        *b_BtagSF_Deepcsv_loose;   //!
   TBranch        *b_BtagSF_Deepcsv_Up_loose;   //!
   TBranch        *b_BtagSF_Deepcsv_Down_loose;   //!
   TBranch        *b_BtagSF_Deepcsv_medium;   //!
   TBranch        *b_BtagSF_Deepcsv_Up_medium;   //!
   TBranch        *b_BtagSF_Deepcsv_Down_medium;   //!
   TBranch        *b_BtagSF_Deepcsv_tight;   //!
   TBranch        *b_BtagSF_Deepcsv_Up_tight;   //!
   TBranch        *b_BtagSF_Deepcsv_Down_tight;   //!
   TBranch        *b_MET_nominal_Pt;   //!
   TBranch        *b_MET_nominal_Px;   //!
   TBranch        *b_MET_nominal_Py;   //!
   TBranch        *b_MET_nominal_phi;   //!
   TBranch        *b_MET_nominal_significance;   //!
   TBranch        *b_MET_Pt;   //!
   TBranch        *b_MET_Px;   //!
   TBranch        *b_MET_Py;   //!
   TBranch        *b_MET_phi;   //!
   TBranch        *b_MET_significance;   //!
   TBranch        *b_MET_T1_Pt;   //!
   TBranch        *b_MET_T1_Px;   //!
   TBranch        *b_MET_T1_Py;   //!
   TBranch        *b_MET_T1_phi;   //!
   TBranch        *b_MET_T1_significance;   //!
   TBranch        *b_MET_gen_pt;   //!
   TBranch        *b_MET_gen_phi;   //!
   TBranch        *b_MET_Type1Unc_Px;   //!
   TBranch        *b_MET_Type1Unc_Py;   //!
   TBranch        *b_MET_Type1Unc_Pt;   //!
   TBranch        *b_MET_T1JetEnDown_Pt;   //!
   TBranch        *b_MET_T1JetEnDown_Px;   //!
   TBranch        *b_MET_T1JetEnDown_Py;   //!
   TBranch        *b_MET_T1JetEnDown_phi;   //!
   TBranch        *b_MET_T1JetEnDown_significance;   //!
   TBranch        *b_MET_T1JetEnUp_Pt;   //!
   TBranch        *b_MET_T1JetEnUp_Px;   //!
   TBranch        *b_MET_T1JetEnUp_Py;   //!
   TBranch        *b_MET_T1JetEnUp_phi;   //!
   TBranch        *b_MET_T1JetEnUp_significance;   //!
   TBranch        *b_MET_T1Smear_Pt;   //!
   TBranch        *b_MET_T1Smear_Px;   //!
   TBranch        *b_MET_T1Smear_Py;   //!
   TBranch        *b_MET_T1Smear_phi;   //!
   TBranch        *b_MET_T1Smear_significance;   //!
   TBranch        *b_MET_T1SmearJetEnDown_Pt;   //!
   TBranch        *b_MET_T1SmearJetEnDown_Px;   //!
   TBranch        *b_MET_T1SmearJetEnDown_Py;   //!
   TBranch        *b_MET_T1SmearJetEnDown_phi;   //!
   TBranch        *b_MET_T1SmearJetEnDown_significance;   //!
   TBranch        *b_MET_T1SmearJetEnUp_Pt;   //!
   TBranch        *b_MET_T1SmearJetEnUp_Px;   //!
   TBranch        *b_MET_T1SmearJetEnUp_Py;   //!
   TBranch        *b_MET_T1SmearJetEnUp_phi;   //!
   TBranch        *b_MET_T1SmearJetEnUp_significance;   //!
   TBranch        *b_MET_T1SmearJetResDown_Pt;   //!
   TBranch        *b_MET_T1SmearJetResDown_Px;   //!
   TBranch        *b_MET_T1SmearJetResDown_Py;   //!
   TBranch        *b_MET_T1SmearJetResDown_phi;   //!
   TBranch        *b_MET_T1SmearJetResDown_significance;   //!
   TBranch        *b_MET_T1SmearJetResUp_Pt;   //!
   TBranch        *b_MET_T1SmearJetResUp_Px;   //!
   TBranch        *b_MET_T1SmearJetResUp_Py;   //!
   TBranch        *b_MET_T1SmearJetResUp_phi;   //!
   TBranch        *b_MET_T1SmearJetResUp_significance;   //!
   TBranch        *b_MET_T1Txy_Pt;   //!
   TBranch        *b_MET_T1Txy_Px;   //!
   TBranch        *b_MET_T1Txy_Py;   //!
   TBranch        *b_MET_T1Txy_phi;   //!
   TBranch        *b_MET_T1Txy_significance;   //!
   TBranch        *b_MET_FinalCollection_Pt;   //!
   TBranch        *b_MET_FinalCollection_Px;   //!
   TBranch        *b_MET_FinalCollection_Py;   //!
   TBranch        *b_MET_FinalCollection_phi;   //!
   TBranch        *b_MET_FinalCollection_significance;   //!
   TBranch        *b_tau_n;   //!
   TBranch        *b_tau_px;   //!
   TBranch        *b_tau_py;   //!
   TBranch        *b_tau_pz;   //!
   TBranch        *b_tau_pt;   //!
   TBranch        *b_tau_eta;   //!
   TBranch        *b_tau_theta;   //!
   TBranch        *b_tau_phi;   //!
   TBranch        *b_tau_energy;   //!
   TBranch        *b_tau_mass;   //!
   TBranch        *b_tau_dxy;   //!
   TBranch        *b_tau_dxy_error;   //!
   TBranch        *b_tau_ptLeadChargedCand;   //!
   TBranch        *b_tau_byIsolationMVArun2v1DBoldDMwLTraw;   //!
   TBranch        *b_tau_byCombinedIsolationDeltaBetaCorrRaw3Hits;   //!
   TBranch        *b_tau_byIsolationMVArun2v1DBnewDMwLTraw;   //!
   TBranch        *b_tau_byIsolationMVArun2v1PWoldDMwLTraw;   //!
   TBranch        *b_tau_byIsolationMVArun2v1PWnewDMwLTraw;   //!
   TBranch        *b_tau_byIsolationMVArun2v1DBdR03oldDMwLTraw;   //!
   TBranch        *b_tau_byIsolationMVArun2v1PWdR03oldDMwLTraw;   //!
   TBranch        *b_tau_againstElectronMVA6Raw;   //!
   TBranch        *b_tau_againstElectronMVA6category;   //!
   TBranch        *b_tau_mc_bestDR;   //!
   TBranch        *b_tau_mc_ERatio;   //!
   TBranch        *b_tau_decayModeFinding;   //!
   TBranch        *b_tau_decayModeFindingNewDMs;   //!
   TBranch        *b_tau_againstMuonLoose3;   //!
   TBranch        *b_tau_againstMuonTight3;   //!
   TBranch        *b_tau_byLooseCombinedIsolationDeltaBetaCorr3Hits;   //!
   TBranch        *b_tau_byMediumCombinedIsolationDeltaBetaCorr3Hits;   //!
   TBranch        *b_tau_byTightCombinedIsolationDeltaBetaCorr3Hits;   //!
   TBranch        *b_tau_byVLooseIsolationMVArun2v1DBoldDMwLT;   //!
   TBranch        *b_tau_byLooseIsolationMVArun2v1DBoldDMwLT;   //!
   TBranch        *b_tau_byMediumIsolationMVArun2v1DBoldDMwLT;   //!
   TBranch        *b_tau_byTightIsolationMVArun2v1DBoldDMwLT;   //!
   TBranch        *b_tau_byVTightIsolationMVArun2v1DBoldDMwLT;   //!
   TBranch        *b_tau_byVVTightIsolationMVArun2v1DBoldDMwLT;   //!
   TBranch        *b_tau_byVLooseIsolationMVArun2v1DBnewDMwLT;   //!
   TBranch        *b_tau_byLooseIsolationMVArun2v1DBnewDMwLT;   //!
   TBranch        *b_tau_byMediumIsolationMVArun2v1DBnewDMwLT;   //!
   TBranch        *b_tau_byTightIsolationMVArun2v1DBnewDMwLT;   //!
   TBranch        *b_tau_byVTightIsolationMVArun2v1DBnewDMwLT;   //!
   TBranch        *b_tau_byVVTightIsolationMVArun2v1DBnewDMwLT;   //!
   TBranch        *b_tau_byVLooseIsolationMVArun2v1PWoldDMwLT;   //!
   TBranch        *b_tau_byLooseIsolationMVArun2v1PWoldDMwLT;   //!
   TBranch        *b_tau_byMediumIsolationMVArun2v1PWoldDMwLT;   //!
   TBranch        *b_tau_byTightIsolationMVArun2v1PWoldDMwLT;   //!
   TBranch        *b_tau_byVTightIsolationMVArun2v1PWoldDMwLT;   //!
   TBranch        *b_tau_byVVTightIsolationMVArun2v1PWoldDMwLT;   //!
   TBranch        *b_tau_byVLooseIsolationMVArun2v1PWnewDMwLT;   //!
   TBranch        *b_tau_byLooseIsolationMVArun2v1PWnewDMwLT;   //!
   TBranch        *b_tau_byMediumIsolationMVArun2v1PWnewDMwLT;   //!
   TBranch        *b_tau_byTightIsolationMVArun2v1PWnewDMwLT;   //!
   TBranch        *b_tau_byVTightIsolationMVArun2v1PWnewDMwLT;   //!
   TBranch        *b_tau_byVVTightIsolationMVArun2v1PWnewDMwLT;   //!
   TBranch        *b_tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT;   //!
   TBranch        *b_tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT;   //!
   TBranch        *b_tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT;   //!
   TBranch        *b_tau_byTightIsolationMVArun2v1DBdR03oldDMwLT;   //!
   TBranch        *b_tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT;   //!
   TBranch        *b_tau_byVVTightIsolationMVArun2v1DBdR03oldDMwLT;   //!
   TBranch        *b_tau_byVLooseIsolationMVArun2v1PWdR03oldDMwLT;   //!
   TBranch        *b_tau_byLooseIsolationMVArun2v1PWdR03oldDMwLT;   //!
   TBranch        *b_tau_byMediumIsolationMVArun2v1PWdR03oldDMwLT;   //!
   TBranch        *b_tau_byTightIsolationMVArun2v1PWdR03oldDMwLT;   //!
   TBranch        *b_tau_byVTightIsolationMVArun2v1PWdR03oldDMwLT;   //!
   TBranch        *b_tau_byVVTightIsolationMVArun2v1PWdR03oldDMwLT;   //!
   TBranch        *b_tau_againstElectronVLooseMVA6;   //!
   TBranch        *b_tau_againstElectronLooseMVA6;   //!
   TBranch        *b_tau_againstElectronMediumMVA6;   //!
   TBranch        *b_tau_againstElectronTightMVA6;   //!
   TBranch        *b_tau_againstElectronVTightMVA6;   //!
   TBranch        *b_tau_chargedIsoPtSum;   //!
   TBranch        *b_tau_neutralIsoPtSum;   //!
   TBranch        *b_tau_puCorrPtSum;   //!
   TBranch        *b_tau_footprintCorrection;   //!
   TBranch        *b_tau_neutralIsoPtSumWeight;   //!
   TBranch        *b_tau_photonPtSumOutsideSignalCone;   //!
   TBranch        *b_tau_byPhotonPtSumOutsideSignalCone;   //!
   TBranch        *b_tau_footprintCorrectiondR03;   //!
   TBranch        *b_tau_chargedIsoPtSumdR03;   //!
   TBranch        *b_tau_neutralIsoPtSumWeightdR03;   //!
   TBranch        *b_tau_neutralIsoPtSumdR03;   //!
   TBranch        *b_tau_photonPtSumOutsideSignalConedR03;   //!
   TBranch        *b_tau_PFChargedHadIso;   //!
   TBranch        *b_tau_PFNeutralHadIso;   //!
   TBranch        *b_tau_PFPhotonIso;   //!
   TBranch        *b_tau_leadChargedParticlePt;   //!
   TBranch        *b_tau_trackRefPt;   //!
   TBranch        *b_tau_lead_dxy;   //!
   TBranch        *b_tau_lead_dz;   //!
   TBranch        *b_tau_dxy_Sig;   //!
   TBranch        *b_tau_flightLengthSig;   //!
   TBranch        *b_tau_ip3d;   //!
   TBranch        *b_tau_ip3d_Sig;   //!
   TBranch        *b_tau_decayDistX;   //!
   TBranch        *b_tau_decayDistY;   //!
   TBranch        *b_tau_decayDistZ;   //!
   TBranch        *b_tau_decayDistMag;   //!
   TBranch        *b_tau_nPhoton;   //!
   TBranch        *b_tau_ptWeightedDetaStrip;   //!
   TBranch        *b_tau_ptWeightedDphiStrip;   //!
   TBranch        *b_tau_ptWeightedDrSignal;   //!
   TBranch        *b_tau_ptWeightedDrIsolation;   //!
   TBranch        *b_tau_leadingTrackChi2;   //!
   TBranch        *b_tau_eRatio;   //!
   TBranch        *b_tau_gjAngleDiff;   //!
   TBranch        *b_tau_numberOfIsolationChargedHadrCands;   //!
   TBranch        *b_tau_numberOfSignalChargedHadrCands;   //!
   TBranch        *b_tau_numNeutralHadronsSignalCone;   //!
   TBranch        *b_tau_numPhotonsSignalCone;   //!
   TBranch        *b_tau_numParticlesSignalCone;   //!
   TBranch        *b_tau_numChargedParticlesIsoCone;   //!
   TBranch        *b_tau_numNeutralHadronsIsoCone;   //!
   TBranch        *b_tau_numPhotonsIsoCone;   //!
   TBranch        *b_tau_numParticlesIsoCone;   //!
   TBranch        *b_tau_mc_index;   //!
   TBranch        *b_tau_decayMode;   //!
   TBranch        *b_tau_charge;   //!
   TBranch        *b_tau_isPFTau;   //!
   TBranch        *b_tau_hasSecondaryVertex;   //!
   TBranch        *b_tau_leadChargedHadrAvailable;   //!
   TBranch        *b_tau_byIsolationMVArun2017v1DBoldDMwLTraw2017;   //!
   TBranch        *b_tau_byIsolationMVArun2017v2DBnewDMwLTraw2017;   //!
   TBranch        *b_tau_byIsolationMVArun2017v2DBoldDMdR0p3wLTraw2017;   //!
   TBranch        *b_tau_byIsolationMVArun2017v2DBoldDMwLTraw2017;   //!
   TBranch        *b_tau_byIsolationMVArun2v1DBnewDMwLTraw2016;   //!
   TBranch        *b_tau_byIsolationMVArun2v1DBoldDMwLTraw2016;   //!
   TBranch        *b_tau_byVVLooseIsolationMVArun2017v1DBoldDMwLT2017;   //!
   TBranch        *b_tau_byVLooseIsolationMVArun2017v1DBoldDMwLT2017;   //!
   TBranch        *b_tau_byLooseIsolationMVArun2017v1DBoldDMwLT2017;   //!
   TBranch        *b_tau_byMediumIsolationMVArun2017v1DBoldDMwLT2017;   //!
   TBranch        *b_tau_byTightIsolationMVArun2017v1DBoldDMwLT2017;   //!
   TBranch        *b_tau_byVTightIsolationMVArun2017v1DBoldDMwLT2017;   //!
   TBranch        *b_tau_byVVTightIsolationMVArun2017v1DBoldDMwLT2017;   //!
   TBranch        *b_tau_byVVLooseIsolationMVArun2017v2DBnewDMwLT2017;   //!
   TBranch        *b_tau_byVLooseIsolationMVArun2017v2DBnewDMwLT2017;   //!
   TBranch        *b_tau_byLooseIsolationMVArun2017v2DBnewDMwLT2017;   //!
   TBranch        *b_tau_byMediumIsolationMVArun2017v2DBnewDMwLT2017;   //!
   TBranch        *b_tau_byTightIsolationMVArun2017v2DBnewDMwLT2017;   //!
   TBranch        *b_tau_byVTightIsolationMVArun2017v2DBnewDMwLT2017;   //!
   TBranch        *b_tau_byVVTightIsolationMVArun2017v2DBnewDMwLT2017;   //!
   TBranch        *b_tau_byVVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017;   //!
   TBranch        *b_tau_byVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017;   //!
   TBranch        *b_tau_byLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017;   //!
   TBranch        *b_tau_byMediumIsolationMVArun2017v2DBoldDMdR0p3wLT2017;   //!
   TBranch        *b_tau_byTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017;   //!
   TBranch        *b_tau_byVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017;   //!
   TBranch        *b_tau_byVVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017;   //!
   TBranch        *b_tau_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017;   //!
   TBranch        *b_tau_byVLooseIsolationMVArun2017v2DBoldDMwLT2017;   //!
   TBranch        *b_tau_byLooseIsolationMVArun2017v2DBoldDMwLT2017;   //!
   TBranch        *b_tau_byMediumIsolationMVArun2017v2DBoldDMwLT2017;   //!
   TBranch        *b_tau_byTightIsolationMVArun2017v2DBoldDMwLT2017;   //!
   TBranch        *b_tau_byVTightIsolationMVArun2017v2DBoldDMwLT2017;   //!
   TBranch        *b_tau_byVVTightIsolationMVArun2017v2DBoldDMwLT2017;   //!
   TBranch        *b_tau_byVLooseIsolationMVArun2v1DBnewDMwLT2016;   //!
   TBranch        *b_tau_byLooseIsolationMVArun2v1DBnewDMwLT2016;   //!
   TBranch        *b_tau_byMediumIsolationMVArun2v1DBnewDMwLT2016;   //!
   TBranch        *b_tau_byTightIsolationMVArun2v1DBnewDMwLT2016;   //!
   TBranch        *b_tau_byVTightIsolationMVArun2v1DBnewDMwLT2016;   //!
   TBranch        *b_tau_byVVTightIsolationMVArun2v1DBnewDMwLT2016;   //!
   TBranch        *b_tau_byVLooseIsolationMVArun2v1DBoldDMwLT2016;   //!
   TBranch        *b_tau_byLooseIsolationMVArun2v1DBoldDMwLT2016;   //!
   TBranch        *b_tau_byMediumIsolationMVArun2v1DBoldDMwLT2016;   //!
   TBranch        *b_tau_byTightIsolationMVArun2v1DBoldDMwLT2016;   //!
   TBranch        *b_tau_byVTightIsolationMVArun2v1DBoldDMwLT2016;   //!
   TBranch        *b_tau_byVVTightIsolationMVArun2v1DBoldDMwLT2016;   //!
   TBranch        *b_trig_Flag_ecalBadCalibReduced;   //!
   TBranch        *b_trig_Flag_HBHENoiseFilter_accept;   //!
   TBranch        *b_trig_Flag_HBHENoiseIsoFilter_accept;   //!
   TBranch        *b_trig_Flag_CSCTightHaloFilter_accept;   //!
   TBranch        *b_trig_Flag_CSCTightHaloTrkMuUnvetoFilter_accept;   //!
   TBranch        *b_trig_Flag_CSCTightHalo2015Filter_accept;   //!
   TBranch        *b_trig_Flag_globalTightHalo2016Filter_accept;   //!
   TBranch        *b_trig_Flag_globalSuperTightHalo2016Filter_accept;   //!
   TBranch        *b_trig_Flag_HcalStripHaloFilter_accept;   //!
   TBranch        *b_trig_Flag_hcalLaserEventFilter_accept;   //!
   TBranch        *b_trig_Flag_EcalDeadCellTriggerPrimitiveFilter_accept;   //!
   TBranch        *b_trig_Flag_EcalDeadCellBoundaryEnergyFilter_accept;   //!
   TBranch        *b_trig_Flag_ecalBadCalibFilter_accept;   //!
   TBranch        *b_trig_Flag_goodVertices_accept;   //!
   TBranch        *b_trig_Flag_eeBadScFilter_accept;   //!
   TBranch        *b_trig_Flag_ecalLaserCorrFilter_accept;   //!
   TBranch        *b_trig_Flag_trkPOGFilters_accept;   //!
   TBranch        *b_trig_Flag_chargedHadronTrackResolutionFilter_accept;   //!
   TBranch        *b_trig_Flag_muonBadTrackFilter_accept;   //!
   TBranch        *b_trig_Flag_BadChargedCandidateFilter_accept;   //!
   TBranch        *b_trig_Flag_BadPFMuonFilter_accept;   //!
   TBranch        *b_trig_Flag_BadChargedCandidateSummer16Filter_accept;   //!
   TBranch        *b_trig_Flag_BadPFMuonSummer16Filter_accept;   //!
   TBranch        *b_trig_Flag_trkPOG_manystripclus53X_accept;   //!
   TBranch        *b_trig_Flag_trkPOG_toomanystripclus53X_accept;   //!
   TBranch        *b_trig_Flag_trkPOG_logErrorTooManyClusters_accept;   //!
   TBranch        *b_trig_Flag_METFilters_accept;   //!
   TBranch        *b_trig_raw2digi_step_accept;   //!
   TBranch        *b_trig_L1Reco_step_accept;   //!
   TBranch        *b_trig_reconstruction_step_accept;   //!
   TBranch        *b_trig_recosim_step_accept;   //!
   TBranch        *b_trig_eventinterpretaion_step_accept;   //!
   TBranch        *b_trig_HLT_DoubleEle25_CaloIdL_MW_accept;   //!
   TBranch        *b_trig_HLT_DoubleEle27_CaloIdL_MW_accept;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_MW_accept;   //!
   TBranch        *b_trig_HLT_DoubleEle24_eta2p1_WPTight_Gsf_accept;   //!
   TBranch        *b_trig_HLT_Ele27_Ele37_CaloIdL_MW_accept;   //!
   TBranch        *b_trig_HLT_Mu27_Ele37_CaloIdL_MW_accept;   //!
   TBranch        *b_trig_HLT_Mu37_Ele27_CaloIdL_MW_accept;   //!
   TBranch        *b_trig_HLT_Mu37_TkMu27_accept;   //!
   TBranch        *b_trig_HLT_DoubleMu4_3_Bs_accept;   //!
   TBranch        *b_trig_HLT_Ele20_WPTight_Gsf_accept;   //!
   TBranch        *b_trig_HLT_Ele15_WPLoose_Gsf_accept;   //!
   TBranch        *b_trig_HLT_Ele17_WPLoose_Gsf_accept;   //!
   TBranch        *b_trig_HLT_Ele20_WPLoose_Gsf_accept;   //!
   TBranch        *b_trig_HLT_Ele20_eta2p1_WPLoose_Gsf_accept;   //!
   TBranch        *b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_accept;   //!
   TBranch        *b_trig_HLT_Ele27_WPTight_Gsf_accept;   //!
   TBranch        *b_trig_HLT_Ele28_WPTight_Gsf_accept;   //!
   TBranch        *b_trig_HLT_Ele30_WPTight_Gsf_accept;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_accept;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_accept;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_accept;   //!
   TBranch        *b_trig_HLT_Ele38_WPTight_Gsf_accept;   //!
   TBranch        *b_trig_HLT_Ele40_WPTight_Gsf_accept;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_accept;   //!
   TBranch        *b_trig_HLT_IsoMu20_accept;   //!
   TBranch        *b_trig_HLT_IsoMu24_accept;   //!
   TBranch        *b_trig_HLT_IsoMu24_eta2p1_accept;   //!
   TBranch        *b_trig_HLT_IsoMu27_accept;   //!
   TBranch        *b_trig_HLT_IsoMu30_accept;   //!
   TBranch        *b_trig_HLT_L1SingleMu18_accept;   //!
   TBranch        *b_trig_HLT_L1SingleMu25_accept;   //!
   TBranch        *b_trig_HLT_L2Mu10_accept;   //!
   TBranch        *b_trig_HLT_L2Mu10_NoVertex_NoBPTX3BX_accept;   //!
   TBranch        *b_trig_HLT_L2Mu10_NoVertex_NoBPTX_accept;   //!
   TBranch        *b_trig_HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX_accept;   //!
   TBranch        *b_trig_HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX_accept;   //!
   TBranch        *b_trig_HLT_L2Mu50_accept;   //!
   TBranch        *b_trig_HLT_L2Mu23NoVtx_2Cha_accept;   //!
   TBranch        *b_trig_HLT_L2Mu23NoVtx_2Cha_CosmicSeed_accept;   //!
   TBranch        *b_trig_HLT_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4_accept;   //!
   TBranch        *b_trig_HLT_DoubleL2Mu30NoVtx_2Cha_Eta2p4_accept;   //!
   TBranch        *b_trig_HLT_DoubleL2Mu50_accept;   //!
   TBranch        *b_trig_HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed_accept;   //!
   TBranch        *b_trig_HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed_NoL2Matched_accept;   //!
   TBranch        *b_trig_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_accept;   //!
   TBranch        *b_trig_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_NoL2Matched_accept;   //!
   TBranch        *b_trig_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4_accept;   //!
   TBranch        *b_trig_HLT_DoubleL2Mu23NoVtx_2Cha_accept;   //!
   TBranch        *b_trig_HLT_DoubleL2Mu23NoVtx_2Cha_NoL2Matched_accept;   //!
   TBranch        *b_trig_HLT_DoubleL2Mu25NoVtx_2Cha_accept;   //!
   TBranch        *b_trig_HLT_DoubleL2Mu25NoVtx_2Cha_NoL2Matched_accept;   //!
   TBranch        *b_trig_HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4_accept;   //!
   TBranch        *b_trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_accept;   //!
   TBranch        *b_trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_accept;   //!
   TBranch        *b_trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_accept;   //!
   TBranch        *b_trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_accept;   //!
   TBranch        *b_trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_accept;   //!
   TBranch        *b_trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8_accept;   //!
   TBranch        *b_trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_accept;   //!
   TBranch        *b_trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8_accept;   //!
   TBranch        *b_trig_HLT_Mu25_TkMu0_Onia_accept;   //!
   TBranch        *b_trig_HLT_Mu30_TkMu0_Psi_accept;   //!
   TBranch        *b_trig_HLT_Mu20_TkMu0_Phi_accept;   //!
   TBranch        *b_trig_HLT_Mu25_TkMu0_Phi_accept;   //!
   TBranch        *b_trig_HLT_Mu12_accept;   //!
   TBranch        *b_trig_HLT_Mu15_accept;   //!
   TBranch        *b_trig_HLT_Mu20_accept;   //!
   TBranch        *b_trig_HLT_Mu27_accept;   //!
   TBranch        *b_trig_HLT_Mu50_accept;   //!
   TBranch        *b_trig_HLT_Mu55_accept;   //!
   TBranch        *b_trig_HLT_OldMu100_accept;   //!
   TBranch        *b_trig_HLT_TkMu100_accept;   //!
   TBranch        *b_trig_HLT_Photon300_NoHE_accept;   //!
   TBranch        *b_trig_HLT_Mu8_TrkIsoVVL_accept;   //!
   TBranch        *b_trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept;   //!
   TBranch        *b_trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept;   //!
   TBranch        *b_trig_HLT_Mu17_TrkIsoVVL_accept;   //!
   TBranch        *b_trig_HLT_Mu19_TrkIsoVVL_accept;   //!
   TBranch        *b_trig_HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL_accept;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_accept;   //!
   TBranch        *b_trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept;   //!
   TBranch        *b_trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept;   //!
   TBranch        *b_trig_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept;   //!
   TBranch        *b_trig_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept;   //!
   TBranch        *b_trig_HLT_Photon20_accept;   //!
   TBranch        *b_trig_HLT_Photon33_accept;   //!
   TBranch        *b_trig_HLT_Photon50_accept;   //!
   TBranch        *b_trig_HLT_Photon75_accept;   //!
   TBranch        *b_trig_HLT_Photon90_accept;   //!
   TBranch        *b_trig_HLT_Photon120_accept;   //!
   TBranch        *b_trig_HLT_Photon150_accept;   //!
   TBranch        *b_trig_HLT_Photon175_accept;   //!
   TBranch        *b_trig_HLT_Photon200_accept;   //!
   TBranch        *b_trig_HLT_Photon100EB_TightID_TightIso_accept;   //!
   TBranch        *b_trig_HLT_Photon110EB_TightID_TightIso_accept;   //!
   TBranch        *b_trig_HLT_Photon120EB_TightID_TightIso_accept;   //!
   TBranch        *b_trig_HLT_Photon100EBHE10_accept;   //!
   TBranch        *b_trig_HLT_Photon100EEHE10_accept;   //!
   TBranch        *b_trig_HLT_Photon100EE_TightID_TightIso_accept;   //!
   TBranch        *b_trig_HLT_Photon50_R9Id90_HE10_IsoM_accept;   //!
   TBranch        *b_trig_HLT_Photon75_R9Id90_HE10_IsoM_accept;   //!
   TBranch        *b_trig_HLT_Photon90_R9Id90_HE10_IsoM_accept;   //!
   TBranch        *b_trig_HLT_Photon120_R9Id90_HE10_IsoM_accept;   //!
   TBranch        *b_trig_HLT_Photon165_R9Id90_HE10_IsoM_accept;   //!
   TBranch        *b_trig_HLT_Photon35_TwoProngs35_accept;   //!
   TBranch        *b_trig_HLT_IsoMu24_TwoProngs35_accept;   //!
   TBranch        *b_trig_HLT_DoubleMu43NoFiltersNoVtx_accept;   //!
   TBranch        *b_trig_HLT_DoubleMu48NoFiltersNoVtx_accept;   //!
   TBranch        *b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_accept;   //!
   TBranch        *b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_accept;   //!
   TBranch        *b_trig_HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL_accept;   //!
   TBranch        *b_trig_HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL_accept;   //!
   TBranch        *b_trig_HLT_DoubleMu33NoFiltersNoVtxDisplaced_accept;   //!
   TBranch        *b_trig_HLT_DoubleMu40NoFiltersNoVtxDisplaced_accept;   //!
   TBranch        *b_trig_HLT_DoubleMu20_7_Mass0to30_L1_DM4_accept;   //!
   TBranch        *b_trig_HLT_DoubleMu20_7_Mass0to30_L1_DM4EG_accept;   //!
   TBranch        *b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_accept;   //!
   TBranch        *b_trig_HLT_Dimuon14_Phi_Barrel_Seagulls_accept;   //!
   TBranch        *b_trig_HLT_Dimuon24_Phi_noCorrL1_accept;   //!
   TBranch        *b_trig_HLT_DoubleIsoMu20_eta2p1_accept;   //!
   TBranch        *b_trig_HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx_accept;   //!
   TBranch        *b_trig_HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx_accept;   //!
   TBranch        *b_trig_HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx_accept;   //!
   TBranch        *b_trig_HLT_Mu8_accept;   //!
   TBranch        *b_trig_HLT_Mu17_accept;   //!
   TBranch        *b_trig_HLT_Mu19_accept;   //!
   TBranch        *b_trig_HLT_Mu17_Photon30_IsoCaloId_accept;   //!
   TBranch        *b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_accept;   //!
   TBranch        *b_trig_HLT_Ele135_CaloIdVT_GsfTrkIdT_accept;   //!
   TBranch        *b_trig_HLT_Ele145_CaloIdVT_GsfTrkIdT_accept;   //!
   TBranch        *b_trig_HLT_Ele200_CaloIdVT_GsfTrkIdT_accept;   //!
   TBranch        *b_trig_HLT_Ele250_CaloIdVT_GsfTrkIdT_accept;   //!
   TBranch        *b_trig_HLT_Ele300_CaloIdVT_GsfTrkIdT_accept;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_accept;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_accept;   //!
   TBranch        *b_trig_HLT_SinglePhoton10_Eta3p1ForPPRef_accept;   //!
   TBranch        *b_trig_HLT_SinglePhoton20_Eta3p1ForPPRef_accept;   //!
   TBranch        *b_trig_HLT_SinglePhoton30_Eta3p1ForPPRef_accept;   //!
   TBranch        *b_trig_HLT_Photon20_HoverELoose_accept;   //!
   TBranch        *b_trig_HLT_Photon30_HoverELoose_accept;   //!
   TBranch        *b_trig_HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142_accept;   //!
   TBranch        *b_trig_AlCa_RPCMuonNormalisation_accept;   //!
   TBranch        *b_trig_MC_DoubleEle5_CaloIdL_MW_accept;   //!
   TBranch        *b_trig_MC_Ele5_WPTight_Gsf_accept;   //!
   TBranch        *b_trig_MC_Ele15_Ele10_CaloIdL_TrackIdL_IsoVL_DZ_accept;   //!
   TBranch        *b_trig_MC_IsoMu_accept;   //!
   TBranch        *b_trig_MC_DoubleMu_TrkIsoVVL_DZ_accept;   //!
   TBranch        *b_trig_MC_DoubleMuNoFiltersNoVtx_accept;   //!
   TBranch        *b_trig_HLT_Mu18_Mu9_SameSign_accept;   //!
   TBranch        *b_trig_HLT_Mu18_Mu9_SameSign_DZ_accept;   //!
   TBranch        *b_trig_HLT_Mu18_Mu9_accept;   //!
   TBranch        *b_trig_HLT_Mu18_Mu9_DZ_accept;   //!
   TBranch        *b_trig_HLT_Mu20_Mu10_SameSign_accept;   //!
   TBranch        *b_trig_HLT_Mu20_Mu10_SameSign_DZ_accept;   //!
   TBranch        *b_trig_HLT_Mu20_Mu10_accept;   //!
   TBranch        *b_trig_HLT_Mu20_Mu10_DZ_accept;   //!
   TBranch        *b_trig_HLT_Mu23_Mu12_SameSign_accept;   //!
   TBranch        *b_trig_HLT_Mu23_Mu12_SameSign_DZ_accept;   //!
   TBranch        *b_trig_HLT_Mu23_Mu12_accept;   //!
   TBranch        *b_trig_HLT_Mu23_Mu12_DZ_accept;   //!
   TBranch        *b_trig_HLT_TrkMu6NoFiltersNoVtx_accept;   //!
   TBranch        *b_trig_HLT_TrkMu16NoFiltersNoVtx_accept;   //!
   TBranch        *b_trig_HLT_DoubleTrkMu_16_6_NoFiltersNoVtx_accept;   //!



		IIHEAnalysis(TTree *tree);
		~IIHEAnalysis();
		void     Init(TTree *tree);
		Int_t GetEntry(int entry);
		Long64_t GetEntries();
		TTree* GetTree();

};

IIHEAnalysis::IIHEAnalysis(TTree *tree)
{
	Init(tree);
}



IIHEAnalysis::~IIHEAnalysis()
{
}



void IIHEAnalysis::Init(TTree *tree)
{
	// The Init() function is called when the selector needs to initialize
	// a new tree or chain. Typically here the branch addresses and branch
	// pointers of the tree will be set.
	// It is normally not necessary to make changes to the generated
	// code, but the routine can be extended by the user if needed.
	// Init() will be called many times when running on PROOF
	// (once per file to be processed).

   LHE_Pt = 0;
   LHE_Eta = 0;
   LHE_Phi = 0;
   LHE_E = 0;
   LHE_pdgid = 0;
   LHE_status = 0;
   LHE_weight_sys = 0;
   mc_index = 0;
   mc_pdgId = 0;
   mc_charge = 0;
   mc_status = 0;
   mc_status_flags = 0;
   mc_status_tau_flags = 0;
   mc_tau_charge = 0;
   mc_tau_pdgId = 0;
   mc_tau_decay = 0;
   mc_tau_had_status = 0;
   mc_tau_had_charge = 0;
   mc_tau_had_pdgId = 0;
   mc_mass = 0;
   mc_px = 0;
   mc_py = 0;
   mc_pz = 0;
   mc_pt = 0;
   mc_eta = 0;
   mc_phi = 0;
   mc_energy = 0;
   mc_tau_pt = 0;
   mc_tau_eta = 0;
   mc_tau_phi = 0;
   mc_tau_energy = 0;
   mc_tau_had_pt = 0;
   mc_tau_had_eta = 0;
   mc_tau_had_phi = 0;
   mc_tau_had_energy = 0;
   mc_numberOfDaughters = 0;
   mc_numberOfMothers = 0;
   mc_mother_index = 0;
   mc_mother_pdgId = 0;
   mc_mother_px = 0;
   mc_mother_py = 0;
   mc_mother_pz = 0;
   mc_mother_pt = 0;
   mc_mother_eta = 0;
   mc_mother_phi = 0;
   mc_mother_energy = 0;
   mc_mother_mass = 0;
   genjet_pt = 0;
   genjet_eta = 0;
   genjet_phi = 0;
   genjet_energy = 0;
   gen_weight_sys = 0;
   pv_x = 0;
   pv_y = 0;
   pv_z = 0;
   pv_ndof = 0;
   pv_normalizedChi2 = 0;
   pv_isValid = 0;
   pv_isFake = 0;
   gsf_classification = 0;
   gsf_ecalTrkEnergyPostCorr = 0;
   gsf_ecalEnergyPostCorr = 0;
   gsf_energy = 0;
   gsf_p = 0;
   gsf_pt = 0;
   gsf_et = 0;
   gsf_scE1x5 = 0;
   gsf_scE5x5 = 0;
   gsf_scE2x5Max = 0;
   gsf_full5x5_e5x5 = 0;
   gsf_full5x5_e1x5 = 0;
   gsf_full5x5_e2x5Max = 0;
   gsf_full5x5_sigmaIetaIeta = 0;
   gsf_full5x5_hcalOverEcal = 0;
   gsf_eta = 0;
   gsf_phi = 0;
   gsf_theta = 0;
   gsf_px = 0;
   gsf_py = 0;
   gsf_pz = 0;
   gsf_caloEnergy = 0;
   gsf_deltaEtaSuperClusterTrackAtVtx = 0;
   gsf_deltaPhiSuperClusterTrackAtVtx = 0;
   gsf_hadronicOverEm = 0;
   gsf_hcalDepth1OverEcal = 0;
   gsf_hcalDepth2OverEcal = 0;
   gsf_dr03TkSumPt = 0;
   gsf_heepTrkPtIso = 0;
   gsf_dr03EcalRecHitSumEt = 0;
   gsf_dr03HcalDepth1TowerSumEt = 0;
   gsf_dr03HcalDepth2TowerSumEt = 0;
   gsf_charge = 0;
   gsf_sigmaIetaIeta = 0;
   gsf_ecaldrivenSeed = 0;
   gsf_trackerdrivenSeed = 0;
   gsf_isEB = 0;
   gsf_isEE = 0;
   gsf_passConversionVeto = 0;
   gsf_VID_cutBasedElectronID_Fall17_94X_V1_loose = 0;
   gsf_VID_cutBasedElectronID_Fall17_94X_V1_medium = 0;
   gsf_VID_cutBasedElectronID_Fall17_94X_V1_tight = 0;
   gsf_VID_cutBasedElectronID_Fall17_94X_V1_veto = 0;
   gsf_VID_cutBasedElectronID_Summer16_80X_V1_loose = 0;
   gsf_VID_cutBasedElectronID_Summer16_80X_V1_medium = 0;
   gsf_VID_cutBasedElectronID_Summer16_80X_V1_tight = 0;
   gsf_VID_cutBasedElectronID_Summer16_80X_V1_veto = 0;
   gsf_VID_heepElectronID_HEEPV70 = 0;
   gsf_VID_mvaEleID_Fall17_iso_V1_wp80 = 0;
   gsf_VID_mvaEleID_Fall17_iso_V1_wp90 = 0;
   gsf_VID_mvaEleID_Fall17_iso_V1_wpLoose = 0;
   gsf_VID_mvaEleID_Fall17_noIso_V1_wp80 = 0;
   gsf_VID_mvaEleID_Fall17_noIso_V1_wp90 = 0;
   gsf_VID_mvaEleID_Fall17_noIso_V1_wpLoose = 0;
   gsf_VID_mvaEleID_Spring16_GeneralPurpose_V1_wp80 = 0;
   gsf_VID_mvaEleID_Spring16_GeneralPurpose_V1_wp90 = 0;
   gsf_VID_mvaEleID_Spring16_HZZ_V1_wpLoose = 0;
   gsf_deltaEtaSeedClusterTrackAtCalo = 0;
   gsf_deltaPhiSeedClusterTrackAtCalo = 0;
   gsf_ecalEnergy = 0;
   gsf_eSuperClusterOverP = 0;
   gsf_dxy = 0;
   gsf_dxy_beamSpot = 0;
   gsf_dxy_firstPVtx = 0;
   gsf_dxyError = 0;
   gsf_dz = 0;
   gsf_dz_beamSpot = 0;
   gsf_dz_firstPVtx = 0;
   gsf_dzError = 0;
   gsf_vz = 0;
   gsf_numberOfValidHits = 0;
   gsf_nLostInnerHits = 0;
   gsf_nLostOuterHits = 0;
   gsf_convFlags = 0;
   gsf_convDist = 0;
   gsf_convDcot = 0;
   gsf_convRadius = 0;
   gsf_fBrem = 0;
   gsf_e1x5 = 0;
   gsf_e2x5Max = 0;
   gsf_e5x5 = 0;
   gsf_r9 = 0;
   gsf_deltaEtaSeedClusterTrackAtVtx = 0;
   gsf_effArea = 0;
   gsf_sumChargedHadronPt = 0;
   gsf_sumNeutralHadronEt = 0;
   gsf_sumPhotonEt = 0;
   gsf_sumPUPt = 0;
   gsf_ooEmooP = 0;
   gsf_hitsinfo = 0;
   gsf_pixelMatch_dPhi1 = 0;
   gsf_pixelMatch_dPhi2 = 0;
   gsf_pixelMatch_dRz1 = 0;
   gsf_pixelMatch_dRz2 = 0;
   gsf_pixelMatch_subDetector1 = 0;
   gsf_pixelMatch_subDetector2 = 0;
   gsf_mc_bestDR = 0;
   gsf_mc_index = 0;
   gsf_mc_ERatio = 0;
   gsf_sc_energy = 0;
   gsf_sc_seed_eta = 0;
   gsf_sc_eta = 0;
   gsf_sc_etacorr = 0;
   gsf_sc_theta = 0;
   gsf_sc_thetacorr = 0;
   gsf_sc_et = 0;
   gsf_sc_phi = 0;
   gsf_sc_px = 0;
   gsf_sc_py = 0;
   gsf_sc_pz = 0;
   gsf_sc_x = 0;
   gsf_sc_y = 0;
   gsf_sc_z = 0;
   gsf_sc_phiWidth = 0;
   gsf_sc_etaWidth = 0;
   gsf_sc_seed_rawId = 0;
   gsf_sc_seed_ieta = 0;
   gsf_sc_seed_iphi = 0;
   gsf_sc_seed_kHasSwitchToGain6 = 0;
   gsf_sc_seed_kHasSwitchToGain1 = 0;
   gsf_swissCross = 0;
   gsf_sc_rawEnergy = 0;
   gsf_sc_preshowerEnergy = 0;
   gsf_sc_lazyTools_e2x5Right = 0;
   gsf_sc_lazyTools_e2x5Left = 0;
   gsf_sc_lazyTools_e2x5Top = 0;
   gsf_sc_lazyTools_e2x5Bottom = 0;
   gsf_sc_lazyTools_eMax = 0;
   gsf_sc_lazyTools_e2nd = 0;
   gsf_sc_lazyTools_eRight = 0;
   gsf_sc_lazyTools_eLeft = 0;
   gsf_sc_lazyTools_eTop = 0;
   gsf_sc_lazyTools_eBottom = 0;
   gsf_sc_lazyTools_e2x2 = 0;
   gsf_sc_lazyTools_e3x3 = 0;
   gsf_sc_lazyTools_e4x4 = 0;
   gsf_sc_lazyTools_e5x5 = 0;
   gsf_sc_lazyTools_e1x3 = 0;
   gsf_sc_lazyTools_e3x1 = 0;
   gsf_sc_lazyTools_e1x5 = 0;
   gsf_sc_lazyTools_e5x1 = 0;
   gsf_sc_lazyTools_eshitsixix = 0;
   gsf_sc_lazyTools_eshitsiyiy = 0;
   gsf_sc_lazyTools_eseffsixix = 0;
   gsf_sc_lazyTools_eseffsiyiy = 0;
   gsf_sc_lazyTools_eseffsirir = 0;
   gsf_sc_lazyTools_BasicClusterSeedTime = 0;
   EBHits_rawId = 0;
   EBHits_iRechit = 0;
   EBHits_ieta = 0;
   EBHits_iphi = 0;
   EBHits_RecoFlag = 0;
   EBHits_kSaturated = 0;
   EBHits_kLeadingEdgeRecovered = 0;
   EBHits_kNeighboursRecovered = 0;
   EBHits_kWeird = 0;
   EBHits_energy = 0;
   EEHits_rawId = 0;
   EEHits_iRechit = 0;
   EEHits_ieta = 0;
   EEHits_iphi = 0;
   EEHits_RecoFlag = 0;
   EEHits_kSaturated = 0;
   EEHits_kLeadingEdgeRecovered = 0;
   EEHits_kNeighboursRecovered = 0;
   EEHits_kWeird = 0;
   EEHits_energy = 0;
   mu_gt_qoverp = 0;
   mu_gt_charge = 0;
   mu_gt_pt = 0;
   mu_gt_eta = 0;
   mu_gt_phi = 0;
   mu_gt_p = 0;
   mu_gt_px = 0;
   mu_gt_py = 0;
   mu_gt_pz = 0;
   mu_gt_theta = 0;
   mu_gt_lambda = 0;
   mu_gt_d0 = 0;
   mu_gt_dz = 0;
   mu_gt_dz_beamspot = 0;
   mu_gt_dz_firstPVtx = 0;
   mu_gt_dxy = 0;
   mu_gt_dxy_beamspot = 0;
   mu_gt_dxy_firstPVtx = 0;
   mu_gt_dsz = 0;
   mu_gt_vx = 0;
   mu_gt_vy = 0;
   mu_gt_vz = 0;
   mu_gt_qoverpError = 0;
   mu_gt_ptError = 0;
   mu_gt_thetaError = 0;
   mu_gt_lambdaError = 0;
   mu_gt_phiError = 0;
   mu_gt_dxyError = 0;
   mu_gt_d0Error = 0;
   mu_gt_dszError = 0;
   mu_gt_dzError = 0;
   mu_gt_etaError = 0;
   mu_gt_chi2 = 0;
   mu_gt_ndof = 0;
   mu_gt_normalizedChi2 = 0;
   mu_ot_qoverp = 0;
   mu_ot_charge = 0;
   mu_ot_pt = 0;
   mu_ot_eta = 0;
   mu_ot_phi = 0;
   mu_ot_p = 0;
   mu_ot_px = 0;
   mu_ot_py = 0;
   mu_ot_pz = 0;
   mu_ot_theta = 0;
   mu_ot_lambda = 0;
   mu_ot_d0 = 0;
   mu_ot_dz = 0;
   mu_ot_dz_beamspot = 0;
   mu_ot_dz_firstPVtx = 0;
   mu_ot_dxy = 0;
   mu_ot_dxy_beamspot = 0;
   mu_ot_dxy_firstPVtx = 0;
   mu_ot_dsz = 0;
   mu_ot_vx = 0;
   mu_ot_vy = 0;
   mu_ot_vz = 0;
   mu_ot_qoverpError = 0;
   mu_ot_ptError = 0;
   mu_ot_thetaError = 0;
   mu_ot_lambdaError = 0;
   mu_ot_phiError = 0;
   mu_ot_dxyError = 0;
   mu_ot_d0Error = 0;
   mu_ot_dszError = 0;
   mu_ot_dzError = 0;
   mu_ot_etaError = 0;
   mu_ot_chi2 = 0;
   mu_ot_ndof = 0;
   mu_ot_normalizedChi2 = 0;
   mu_it_qoverp = 0;
   mu_it_charge = 0;
   mu_it_pt = 0;
   mu_it_eta = 0;
   mu_it_phi = 0;
   mu_it_p = 0;
   mu_it_px = 0;
   mu_it_py = 0;
   mu_it_pz = 0;
   mu_it_theta = 0;
   mu_it_lambda = 0;
   mu_it_d0 = 0;
   mu_it_dz = 0;
   mu_it_dz_beamspot = 0;
   mu_it_dz_firstPVtx = 0;
   mu_it_dxy = 0;
   mu_it_dxy_beamspot = 0;
   mu_it_dxy_firstPVtx = 0;
   mu_it_dsz = 0;
   mu_it_vx = 0;
   mu_it_vy = 0;
   mu_it_vz = 0;
   mu_it_qoverpError = 0;
   mu_it_ptError = 0;
   mu_it_thetaError = 0;
   mu_it_lambdaError = 0;
   mu_it_phiError = 0;
   mu_it_dxyError = 0;
   mu_it_d0Error = 0;
   mu_it_dszError = 0;
   mu_it_dzError = 0;
   mu_it_etaError = 0;
   mu_it_chi2 = 0;
   mu_it_ndof = 0;
   mu_it_normalizedChi2 = 0;
   mu_ibt_qoverp = 0;
   mu_ibt_charge = 0;
   mu_ibt_pt = 0;
   mu_ibt_eta = 0;
   mu_ibt_phi = 0;
   mu_ibt_p = 0;
   mu_ibt_px = 0;
   mu_ibt_py = 0;
   mu_ibt_pz = 0;
   mu_ibt_theta = 0;
   mu_ibt_lambda = 0;
   mu_ibt_d0 = 0;
   mu_ibt_dz = 0;
   mu_ibt_dz_beamspot = 0;
   mu_ibt_dz_firstPVtx = 0;
   mu_ibt_dxy = 0;
   mu_ibt_dxy_beamspot = 0;
   mu_ibt_dxy_firstPVtx = 0;
   mu_ibt_dsz = 0;
   mu_ibt_vx = 0;
   mu_ibt_vy = 0;
   mu_ibt_vz = 0;
   mu_ibt_qoverpError = 0;
   mu_ibt_ptError = 0;
   mu_ibt_thetaError = 0;
   mu_ibt_lambdaError = 0;
   mu_ibt_phiError = 0;
   mu_ibt_dxyError = 0;
   mu_ibt_d0Error = 0;
   mu_ibt_dszError = 0;
   mu_ibt_dzError = 0;
   mu_ibt_etaError = 0;
   mu_ibt_chi2 = 0;
   mu_ibt_ndof = 0;
   mu_ibt_normalizedChi2 = 0;
   mu_isGlobalMuon = 0;
   mu_isStandAloneMuon = 0;
   mu_isTrackerMuon = 0;
   mu_isPFMuon = 0;
   mu_isPFIsolationValid = 0;
   mu_isGoodMuonTMLastStationLoose = 0;
   mu_isGoodMuonTMLastStationTight = 0;
   mu_isGoodMuonTM2DCompatibilityLoose = 0;
   mu_isGoodMuonTM2DCompatibilityTight = 0;
   mu_isGoodMuonTMOneStationLoose = 0;
   mu_isGoodMuonTMOneStationTight = 0;
   mu_isGoodMuonTMLastStationOptimizedLowPtLoose = 0;
   mu_isGoodMuonTMLastStationOptimizedLowPtTight = 0;
   mu_isTightMuon = 0;
   mu_isMediumMuon = 0;
   mu_isLooseMuon = 0;
   mu_isSoftMuon = 0;
   mu_isHighPtMuonOld = 0;
   mu_isHighPtMuon = 0;
   mu_isTrackerHighPtMuon = 0;
   mu_CutBasedIdLoose = 0;
   mu_CutBasedIdMedium = 0;
   mu_CutBasedIdMediumPrompt = 0;
   mu_CutBasedIdTight = 0;
   mu_CutBasedIdGlobalHighPt = 0;
   mu_CutBasedIdTrkHighPt = 0;
   mu_PFIsoVeryLoose = 0;
   mu_PFIsoLoose = 0;
   mu_PFIsoMedium = 0;
   mu_PFIsoTight = 0;
   mu_PFIsoVeryTight = 0;
   mu_PFIsoVeryVeryTight = 0;
   mu_TkIsoLoose = 0;
   mu_TkIsoTight = 0;
   mu_SoftCutBasedId = 0;
   mu_SoftMvaId = 0;
   mu_MvaLoose = 0;
   mu_MvaMedium = 0;
   mu_MvaTight = 0;
   mu_MiniIsoLoose = 0;
   mu_MiniIsoMedium = 0;
   mu_MiniIsoTight = 0;
   mu_MiniIsoVeryTight = 0;
   mu_TriggerIdLoose = 0;
   mu_InTimeMuon = 0;
   mu_MultiIsoLoose = 0;
   mu_MultiIsoMedium = 0;
   mu_numberOfMatchedStations = 0;
   mu_numberOfValidPixelHits = 0;
   mu_trackerLayersWithMeasurement = 0;
   mu_numberOfValidMuonHits = 0;
   mu_pixelLayersWithMeasurement = 0;
   mu_innerTrack_validFraction = 0;
   mu_combinedQuality_trkKink = 0;
   mu_combinedQuality_chi2LocalPosition = 0;
   mu_segmentCompatibility = 0;
   mu_dB = 0;
   mu_pt_default = 0;
   mu_isolationR03_sumPt = 0;
   mu_isolationR03_trackerVetoPt = 0;
   mu_isolationR03_emEt = 0;
   mu_isolationR03_emVetoEt = 0;
   mu_isolationR03_hadEt = 0;
   mu_isolationR03_hadVetoEt = 0;
   mu_isolationR05_sumPt = 0;
   mu_isolationR05_trackerVetoPt = 0;
   mu_isolationR05_emEt = 0;
   mu_isolationR05_emVetoEt = 0;
   mu_isolationR05_hadEt = 0;
   mu_isolationR05_hadVetoEt = 0;
   mu_pfIsolationR03_sumChargedHadronPt = 0;
   mu_pfIsolationR03_sumNeutralHadronEt = 0;
   mu_pfIsolationR03_sumChargedParticlePt = 0;
   mu_pfIsolationR03_sumPhotonEt = 0;
   mu_pfIsolationR03_sumNeutralHadronEtHighThreshold = 0;
   mu_pfIsolationR03_sumPhotonEtHighThreshold = 0;
   mu_pfIsolationR03_sumPUPt = 0;
   mu_pfIsolationR04_sumChargedHadronPt = 0;
   mu_pfIsolationR04_sumNeutralHadronEt = 0;
   mu_pfIsolationR04_sumChargedParticlePt = 0;
   mu_pfIsolationR04_sumPhotonEt = 0;
   mu_pfIsolationR04_sumNeutralHadronEtHighThreshold = 0;
   mu_pfIsolationR04_sumPhotonEtHighThreshold = 0;
   mu_pfIsolationR04_sumPUPt = 0;
   mu_pfIsoDbCorrected03 = 0;
   mu_pfIsoDbCorrected04 = 0;
   mu_isoTrackerBased03 = 0;
   mu_mc_bestDR = 0;
   mu_mc_index = 0;
   mu_mc_ERatio = 0;
   jet_px = 0;
   jet_py = 0;
   jet_pz = 0;
   jet_pt = 0;
   jet_eta = 0;
   jet_theta = 0;
   jet_phi = 0;
   jet_energy = 0;
   jet_mass = 0;
   jet_chargedEmEnergyFraction = 0;
   jet_neutralHadronEnergyFraction = 0;
   jet_neutralEmEnergyFraction = 0;
   jet_chargedHadronEnergyFraction = 0;
   jet_muonEnergyFraction = 0;
   jet_chargedMultiplicity = 0;
   jet_neutralMultiplicity = 0;
   jet_partonFlavour = 0;
   jet_hadronFlavour = 0;
   jet_CSVv2 = 0;
   jet_CvsL = 0;
   jet_CvsB = 0;
   jet_MVA2BJets = 0;
   jet_CvsB_DeepJet_charm_tagger = 0;
   jet_CvsL_DeepJet_charm_tagger = 0;
   jet_CvsB_DeepCSV_charm_tagger = 0;
   jet_CvsL_DeepCSV_charm_tagger = 0;
   jet_DeepJet = 0;
   jet_DeepCSV = 0;
   jet_isJetIDLoose_2016 = 0;
   jet_isJetIDTight_2016 = 0;
   jet_isJetIDTightLepVeto_2016 = 0;
   jet_isJetID_2017 = 0;
   jet_isJetIDLepVeto_2017 = 0;
   jet_isJetID_2018 = 0;
   jet_isJetIDLepVeto_2018 = 0;
   jet_Smeared_pt = 0;
   jet_SmearedJetResUp_pt = 0;
   jet_SmearedJetResDown_pt = 0;
   jet_EnUp_pt = 0;
   jet_EnDown_pt = 0;
   MET_Type1Unc_Px = 0;
   MET_Type1Unc_Py = 0;
   MET_Type1Unc_Pt = 0;
   tau_px = 0;
   tau_py = 0;
   tau_pz = 0;
   tau_pt = 0;
   tau_eta = 0;
   tau_theta = 0;
   tau_phi = 0;
   tau_energy = 0;
   tau_mass = 0;
   tau_dxy = 0;
   tau_dxy_error = 0;
   tau_ptLeadChargedCand = 0;
   tau_byIsolationMVArun2v1DBoldDMwLTraw = 0;
   tau_byCombinedIsolationDeltaBetaCorrRaw3Hits = 0;
   tau_byIsolationMVArun2v1DBnewDMwLTraw = 0;
   tau_byIsolationMVArun2v1PWoldDMwLTraw = 0;
   tau_byIsolationMVArun2v1PWnewDMwLTraw = 0;
   tau_byIsolationMVArun2v1DBdR03oldDMwLTraw = 0;
   tau_byIsolationMVArun2v1PWdR03oldDMwLTraw = 0;
   tau_againstElectronMVA6Raw = 0;
   tau_againstElectronMVA6category = 0;
   tau_mc_bestDR = 0;
   tau_mc_ERatio = 0;
   tau_decayModeFinding = 0;
   tau_decayModeFindingNewDMs = 0;
   tau_againstMuonLoose3 = 0;
   tau_againstMuonTight3 = 0;
   tau_byLooseCombinedIsolationDeltaBetaCorr3Hits = 0;
   tau_byMediumCombinedIsolationDeltaBetaCorr3Hits = 0;
   tau_byTightCombinedIsolationDeltaBetaCorr3Hits = 0;
   tau_byVLooseIsolationMVArun2v1DBoldDMwLT = 0;
   tau_byLooseIsolationMVArun2v1DBoldDMwLT = 0;
   tau_byMediumIsolationMVArun2v1DBoldDMwLT = 0;
   tau_byTightIsolationMVArun2v1DBoldDMwLT = 0;
   tau_byVTightIsolationMVArun2v1DBoldDMwLT = 0;
   tau_byVVTightIsolationMVArun2v1DBoldDMwLT = 0;
   tau_byVLooseIsolationMVArun2v1DBnewDMwLT = 0;
   tau_byLooseIsolationMVArun2v1DBnewDMwLT = 0;
   tau_byMediumIsolationMVArun2v1DBnewDMwLT = 0;
   tau_byTightIsolationMVArun2v1DBnewDMwLT = 0;
   tau_byVTightIsolationMVArun2v1DBnewDMwLT = 0;
   tau_byVVTightIsolationMVArun2v1DBnewDMwLT = 0;
   tau_byVLooseIsolationMVArun2v1PWoldDMwLT = 0;
   tau_byLooseIsolationMVArun2v1PWoldDMwLT = 0;
   tau_byMediumIsolationMVArun2v1PWoldDMwLT = 0;
   tau_byTightIsolationMVArun2v1PWoldDMwLT = 0;
   tau_byVTightIsolationMVArun2v1PWoldDMwLT = 0;
   tau_byVVTightIsolationMVArun2v1PWoldDMwLT = 0;
   tau_byVLooseIsolationMVArun2v1PWnewDMwLT = 0;
   tau_byLooseIsolationMVArun2v1PWnewDMwLT = 0;
   tau_byMediumIsolationMVArun2v1PWnewDMwLT = 0;
   tau_byTightIsolationMVArun2v1PWnewDMwLT = 0;
   tau_byVTightIsolationMVArun2v1PWnewDMwLT = 0;
   tau_byVVTightIsolationMVArun2v1PWnewDMwLT = 0;
   tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT = 0;
   tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT = 0;
   tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT = 0;
   tau_byTightIsolationMVArun2v1DBdR03oldDMwLT = 0;
   tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT = 0;
   tau_byVVTightIsolationMVArun2v1DBdR03oldDMwLT = 0;
   tau_byVLooseIsolationMVArun2v1PWdR03oldDMwLT = 0;
   tau_byLooseIsolationMVArun2v1PWdR03oldDMwLT = 0;
   tau_byMediumIsolationMVArun2v1PWdR03oldDMwLT = 0;
   tau_byTightIsolationMVArun2v1PWdR03oldDMwLT = 0;
   tau_byVTightIsolationMVArun2v1PWdR03oldDMwLT = 0;
   tau_byVVTightIsolationMVArun2v1PWdR03oldDMwLT = 0;
   tau_againstElectronVLooseMVA6 = 0;
   tau_againstElectronLooseMVA6 = 0;
   tau_againstElectronMediumMVA6 = 0;
   tau_againstElectronTightMVA6 = 0;
   tau_againstElectronVTightMVA6 = 0;
   tau_chargedIsoPtSum = 0;
   tau_neutralIsoPtSum = 0;
   tau_puCorrPtSum = 0;
   tau_footprintCorrection = 0;
   tau_neutralIsoPtSumWeight = 0;
   tau_photonPtSumOutsideSignalCone = 0;
   tau_byPhotonPtSumOutsideSignalCone = 0;
   tau_footprintCorrectiondR03 = 0;
   tau_chargedIsoPtSumdR03 = 0;
   tau_neutralIsoPtSumWeightdR03 = 0;
   tau_neutralIsoPtSumdR03 = 0;
   tau_photonPtSumOutsideSignalConedR03 = 0;
   tau_PFChargedHadIso = 0;
   tau_PFNeutralHadIso = 0;
   tau_PFPhotonIso = 0;
   tau_leadChargedParticlePt = 0;
   tau_trackRefPt = 0;
   tau_lead_dxy = 0;
   tau_lead_dz = 0;
   tau_dxy_Sig = 0;
   tau_flightLengthSig = 0;
   tau_ip3d = 0;
   tau_ip3d_Sig = 0;
   tau_decayDistX = 0;
   tau_decayDistY = 0;
   tau_decayDistZ = 0;
   tau_decayDistMag = 0;
   tau_nPhoton = 0;
   tau_ptWeightedDetaStrip = 0;
   tau_ptWeightedDphiStrip = 0;
   tau_ptWeightedDrSignal = 0;
   tau_ptWeightedDrIsolation = 0;
   tau_leadingTrackChi2 = 0;
   tau_eRatio = 0;
   tau_gjAngleDiff = 0;
   tau_numberOfIsolationChargedHadrCands = 0;
   tau_numberOfSignalChargedHadrCands = 0;
   tau_numNeutralHadronsSignalCone = 0;
   tau_numPhotonsSignalCone = 0;
   tau_numParticlesSignalCone = 0;
   tau_numChargedParticlesIsoCone = 0;
   tau_numNeutralHadronsIsoCone = 0;
   tau_numPhotonsIsoCone = 0;
   tau_numParticlesIsoCone = 0;
   tau_mc_index = 0;
   tau_decayMode = 0;
   tau_charge = 0;
   tau_isPFTau = 0;
   tau_hasSecondaryVertex = 0;
   tau_leadChargedHadrAvailable = 0;
   tau_byIsolationMVArun2017v1DBoldDMwLTraw2017 = 0;
   tau_byIsolationMVArun2017v2DBnewDMwLTraw2017 = 0;
   tau_byIsolationMVArun2017v2DBoldDMdR0p3wLTraw2017 = 0;
   tau_byIsolationMVArun2017v2DBoldDMwLTraw2017 = 0;
   tau_byIsolationMVArun2v1DBnewDMwLTraw2016 = 0;
   tau_byIsolationMVArun2v1DBoldDMwLTraw2016 = 0;
   tau_byVVLooseIsolationMVArun2017v1DBoldDMwLT2017 = 0;
   tau_byVLooseIsolationMVArun2017v1DBoldDMwLT2017 = 0;
   tau_byLooseIsolationMVArun2017v1DBoldDMwLT2017 = 0;
   tau_byMediumIsolationMVArun2017v1DBoldDMwLT2017 = 0;
   tau_byTightIsolationMVArun2017v1DBoldDMwLT2017 = 0;
   tau_byVTightIsolationMVArun2017v1DBoldDMwLT2017 = 0;
   tau_byVVTightIsolationMVArun2017v1DBoldDMwLT2017 = 0;
   tau_byVVLooseIsolationMVArun2017v2DBnewDMwLT2017 = 0;
   tau_byVLooseIsolationMVArun2017v2DBnewDMwLT2017 = 0;
   tau_byLooseIsolationMVArun2017v2DBnewDMwLT2017 = 0;
   tau_byMediumIsolationMVArun2017v2DBnewDMwLT2017 = 0;
   tau_byTightIsolationMVArun2017v2DBnewDMwLT2017 = 0;
   tau_byVTightIsolationMVArun2017v2DBnewDMwLT2017 = 0;
   tau_byVVTightIsolationMVArun2017v2DBnewDMwLT2017 = 0;
   tau_byVVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017 = 0;
   tau_byVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017 = 0;
   tau_byLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017 = 0;
   tau_byMediumIsolationMVArun2017v2DBoldDMdR0p3wLT2017 = 0;
   tau_byTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017 = 0;
   tau_byVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017 = 0;
   tau_byVVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017 = 0;
   tau_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017 = 0;
   tau_byVLooseIsolationMVArun2017v2DBoldDMwLT2017 = 0;
   tau_byLooseIsolationMVArun2017v2DBoldDMwLT2017 = 0;
   tau_byMediumIsolationMVArun2017v2DBoldDMwLT2017 = 0;
   tau_byTightIsolationMVArun2017v2DBoldDMwLT2017 = 0;
   tau_byVTightIsolationMVArun2017v2DBoldDMwLT2017 = 0;
   tau_byVVTightIsolationMVArun2017v2DBoldDMwLT2017 = 0;
   tau_byVLooseIsolationMVArun2v1DBnewDMwLT2016 = 0;
   tau_byLooseIsolationMVArun2v1DBnewDMwLT2016 = 0;
   tau_byMediumIsolationMVArun2v1DBnewDMwLT2016 = 0;
   tau_byTightIsolationMVArun2v1DBnewDMwLT2016 = 0;
   tau_byVTightIsolationMVArun2v1DBnewDMwLT2016 = 0;
   tau_byVVTightIsolationMVArun2v1DBnewDMwLT2016 = 0;
   tau_byVLooseIsolationMVArun2v1DBoldDMwLT2016 = 0;
   tau_byLooseIsolationMVArun2v1DBoldDMwLT2016 = 0;
   tau_byMediumIsolationMVArun2v1DBoldDMwLT2016 = 0;
   tau_byTightIsolationMVArun2v1DBoldDMwLT2016 = 0;
   tau_byVTightIsolationMVArun2v1DBoldDMwLT2016 = 0;
   tau_byVVTightIsolationMVArun2v1DBoldDMwLT2016 = 0;


	// Set object pointer


	//


	// Set branch addresses and branch pointers

	if (!tree) return;
	_tree = tree;
	_tree->SetMakeClass(1);

	//



   _tree->SetBranchAddress("ev_event", &ev_event, &b_ev_event);
   _tree->SetBranchAddress("ev_run", &ev_run, &b_ev_run);
   _tree->SetBranchAddress("ev_luminosityBlock", &ev_luminosityBlock, &b_ev_luminosityBlock);
   _tree->SetBranchAddress("ev_time", &ev_time, &b_ev_time);
   _tree->SetBranchAddress("ev_time_unixTime", &ev_time_unixTime, &b_ev_time_unixTime);
   _tree->SetBranchAddress("ev_time_microsecondOffset", &ev_time_microsecondOffset, &b_ev_time_microsecondOffset);
   _tree->SetBranchAddress("ev_fixedGridRhoAll", &ev_fixedGridRhoAll, &b_ev_fixedGridRhoAll);
   _tree->SetBranchAddress("ev_fixedGridRhoFastjetAll", &ev_fixedGridRhoFastjetAll, &b_ev_fixedGridRhoFastjetAll);
   _tree->SetBranchAddress("ev_fixedGridRhoFastjetAllCalo", &ev_fixedGridRhoFastjetAllCalo, &b_ev_fixedGridRhoFastjetAllCalo);
   _tree->SetBranchAddress("ev_fixedGridRhoFastjetCentralCalo", &ev_fixedGridRhoFastjetCentralCalo, &b_ev_fixedGridRhoFastjetCentralCalo);
   _tree->SetBranchAddress("ev_fixedGridRhoFastjetCentralChargedPileUp", &ev_fixedGridRhoFastjetCentralChargedPileUp, &b_ev_fixedGridRhoFastjetCentralChargedPileUp);
   _tree->SetBranchAddress("ev_fixedGridRhoFastjetCentralNeutral", &ev_fixedGridRhoFastjetCentralNeutral, &b_ev_fixedGridRhoFastjetCentralNeutral);





   if(_tree->GetListOfBranches()->FindObject("LHE_Pt")) {

   _tree->SetBranchAddress("LHE_Pt", &LHE_Pt, &b_LHE_Pt);
   _tree->SetBranchAddress("LHE_Eta", &LHE_Eta, &b_LHE_Eta);
   _tree->SetBranchAddress("LHE_Phi", &LHE_Phi, &b_LHE_Phi);
   _tree->SetBranchAddress("LHE_E", &LHE_E, &b_LHE_E);
   _tree->SetBranchAddress("LHE_pdgid", &LHE_pdgid, &b_LHE_pdgid);
   _tree->SetBranchAddress("LHE_status", &LHE_status, &b_LHE_status);
   _tree->SetBranchAddress("LHE_weight_nominal", &LHE_weight_nominal, &b_LHE_weight_nominal);
   _tree->SetBranchAddress("LHE_weight_sys", &LHE_weight_sys, &b_LHE_weight_sys);
   _tree->SetBranchAddress("mc_n", &mc_n, &b_mc_n);
   _tree->SetBranchAddress("mc_weight", &mc_weight, &b_mc_weight);
   _tree->SetBranchAddress("mc_w_sign", &mc_w_sign, &b_mc_w_sign);
   _tree->SetBranchAddress("mc_id_first", &mc_id_first, &b_mc_id_first);
   _tree->SetBranchAddress("mc_id_second", &mc_id_second, &b_mc_id_second);
   _tree->SetBranchAddress("mc_x_first", &mc_x_first, &b_mc_x_first);
   _tree->SetBranchAddress("mc_x_second", &mc_x_second, &b_mc_x_second);
   _tree->SetBranchAddress("mc_xPDF_first", &mc_xPDF_first, &b_mc_xPDF_first);
   _tree->SetBranchAddress("mc_xPDF_second", &mc_xPDF_second, &b_mc_xPDF_second);
   _tree->SetBranchAddress("mc_scalePDF", &mc_scalePDF, &b_mc_scalePDF);
   _tree->SetBranchAddress("mc_index", &mc_index, &b_mc_index);
   _tree->SetBranchAddress("mc_pdgId", &mc_pdgId, &b_mc_pdgId);
   _tree->SetBranchAddress("mc_charge", &mc_charge, &b_mc_charge);
   _tree->SetBranchAddress("mc_status", &mc_status, &b_mc_status);
   _tree->SetBranchAddress("mc_status_flags", &mc_status_flags, &b_mc_status_flags);
   _tree->SetBranchAddress("mc_status_tau_flags", &mc_status_tau_flags, &b_mc_status_tau_flags);
   _tree->SetBranchAddress("mc_tau_charge", &mc_tau_charge, &b_mc_tau_charge);
   _tree->SetBranchAddress("mc_tau_pdgId", &mc_tau_pdgId, &b_mc_tau_pdgId);
   _tree->SetBranchAddress("mc_tau_decay", &mc_tau_decay, &b_mc_tau_decay);
   _tree->SetBranchAddress("mc_tau_had_status", &mc_tau_had_status, &b_mc_tau_had_status);
   _tree->SetBranchAddress("mc_tau_had_charge", &mc_tau_had_charge, &b_mc_tau_had_charge);
   _tree->SetBranchAddress("mc_tau_had_pdgId", &mc_tau_had_pdgId, &b_mc_tau_had_pdgId);
   _tree->SetBranchAddress("mc_mass", &mc_mass, &b_mc_mass);
   _tree->SetBranchAddress("mc_px", &mc_px, &b_mc_px);
   _tree->SetBranchAddress("mc_py", &mc_py, &b_mc_py);
   _tree->SetBranchAddress("mc_pz", &mc_pz, &b_mc_pz);
   _tree->SetBranchAddress("mc_pt", &mc_pt, &b_mc_pt);
   _tree->SetBranchAddress("mc_eta", &mc_eta, &b_mc_eta);
   _tree->SetBranchAddress("mc_phi", &mc_phi, &b_mc_phi);
   _tree->SetBranchAddress("mc_energy", &mc_energy, &b_mc_energy);
   _tree->SetBranchAddress("mc_tau_pt", &mc_tau_pt, &b_mc_tau_pt);
   _tree->SetBranchAddress("mc_tau_eta", &mc_tau_eta, &b_mc_tau_eta);
   _tree->SetBranchAddress("mc_tau_phi", &mc_tau_phi, &b_mc_tau_phi);
   _tree->SetBranchAddress("mc_tau_energy", &mc_tau_energy, &b_mc_tau_energy);
   _tree->SetBranchAddress("mc_tau_had_pt", &mc_tau_had_pt, &b_mc_tau_had_pt);
   _tree->SetBranchAddress("mc_tau_had_eta", &mc_tau_had_eta, &b_mc_tau_had_eta);
   _tree->SetBranchAddress("mc_tau_had_phi", &mc_tau_had_phi, &b_mc_tau_had_phi);
   _tree->SetBranchAddress("mc_tau_had_energy", &mc_tau_had_energy, &b_mc_tau_had_energy);
   _tree->SetBranchAddress("mc_numberOfDaughters", &mc_numberOfDaughters, &b_mc_numberOfDaughters);
   _tree->SetBranchAddress("mc_numberOfMothers", &mc_numberOfMothers, &b_mc_numberOfMothers);
   _tree->SetBranchAddress("mc_mother_index", &mc_mother_index, &b_mc_mother_index);
   _tree->SetBranchAddress("mc_mother_pdgId", &mc_mother_pdgId, &b_mc_mother_pdgId);
   _tree->SetBranchAddress("mc_mother_px", &mc_mother_px, &b_mc_mother_px);
   _tree->SetBranchAddress("mc_mother_py", &mc_mother_py, &b_mc_mother_py);
   _tree->SetBranchAddress("mc_mother_pz", &mc_mother_pz, &b_mc_mother_pz);
   _tree->SetBranchAddress("mc_mother_pt", &mc_mother_pt, &b_mc_mother_pt);
   _tree->SetBranchAddress("mc_mother_eta", &mc_mother_eta, &b_mc_mother_eta);
   _tree->SetBranchAddress("mc_mother_phi", &mc_mother_phi, &b_mc_mother_phi);
   _tree->SetBranchAddress("mc_mother_energy", &mc_mother_energy, &b_mc_mother_energy);
   _tree->SetBranchAddress("mc_mother_mass", &mc_mother_mass, &b_mc_mother_mass);
   _tree->SetBranchAddress("mc_trueNumInteractions", &mc_trueNumInteractions, &b_mc_trueNumInteractions);
   _tree->SetBranchAddress("mc_PU_NumInteractions", &mc_PU_NumInteractions, &b_mc_PU_NumInteractions);
   _tree->SetBranchAddress("genjet_pt", &genjet_pt, &b_genjet_pt);
   _tree->SetBranchAddress("genjet_eta", &genjet_eta, &b_genjet_eta);
   _tree->SetBranchAddress("genjet_phi", &genjet_phi, &b_genjet_phi);
   _tree->SetBranchAddress("genjet_energy", &genjet_energy, &b_genjet_energy);
   _tree->SetBranchAddress("gen_weight_sys", &gen_weight_sys, &b_gen_weight_sys);




   }


   _tree->SetBranchAddress("pv_n", &pv_n, &b_pv_n);
   _tree->SetBranchAddress("pv_x", &pv_x, &b_pv_x);
   _tree->SetBranchAddress("pv_y", &pv_y, &b_pv_y);
   _tree->SetBranchAddress("pv_z", &pv_z, &b_pv_z);
   _tree->SetBranchAddress("pv_ndof", &pv_ndof, &b_pv_ndof);
   _tree->SetBranchAddress("pv_normalizedChi2", &pv_normalizedChi2, &b_pv_normalizedChi2);
   _tree->SetBranchAddress("pv_isValid", &pv_isValid, &b_pv_isValid);
   _tree->SetBranchAddress("pv_isFake", &pv_isFake, &b_pv_isFake);
   _tree->SetBranchAddress("gsf_n", &gsf_n, &b_gsf_n);
   _tree->SetBranchAddress("gsf_classification", &gsf_classification, &b_gsf_classification);
   _tree->SetBranchAddress("gsf_ecalTrkEnergyPostCorr", &gsf_ecalTrkEnergyPostCorr, &b_gsf_ecalTrkEnergyPostCorr);
   _tree->SetBranchAddress("gsf_ecalEnergyPostCorr", &gsf_ecalEnergyPostCorr, &b_gsf_ecalEnergyPostCorr);
   _tree->SetBranchAddress("gsf_energy", &gsf_energy, &b_gsf_energy);
   _tree->SetBranchAddress("gsf_p", &gsf_p, &b_gsf_p);
   _tree->SetBranchAddress("gsf_pt", &gsf_pt, &b_gsf_pt);
   _tree->SetBranchAddress("gsf_et", &gsf_et, &b_gsf_et);
   _tree->SetBranchAddress("gsf_scE1x5", &gsf_scE1x5, &b_gsf_scE1x5);
   _tree->SetBranchAddress("gsf_scE5x5", &gsf_scE5x5, &b_gsf_scE5x5);
   _tree->SetBranchAddress("gsf_scE2x5Max", &gsf_scE2x5Max, &b_gsf_scE2x5Max);
   _tree->SetBranchAddress("gsf_full5x5_e5x5", &gsf_full5x5_e5x5, &b_gsf_full5x5_e5x5);
   _tree->SetBranchAddress("gsf_full5x5_e1x5", &gsf_full5x5_e1x5, &b_gsf_full5x5_e1x5);
   _tree->SetBranchAddress("gsf_full5x5_e2x5Max", &gsf_full5x5_e2x5Max, &b_gsf_full5x5_e2x5Max);
   _tree->SetBranchAddress("gsf_full5x5_sigmaIetaIeta", &gsf_full5x5_sigmaIetaIeta, &b_gsf_full5x5_sigmaIetaIeta);
   _tree->SetBranchAddress("gsf_full5x5_hcalOverEcal", &gsf_full5x5_hcalOverEcal, &b_gsf_full5x5_hcalOverEcal);
   _tree->SetBranchAddress("gsf_eta", &gsf_eta, &b_gsf_eta);
   _tree->SetBranchAddress("gsf_phi", &gsf_phi, &b_gsf_phi);
   _tree->SetBranchAddress("gsf_theta", &gsf_theta, &b_gsf_theta);
   _tree->SetBranchAddress("gsf_px", &gsf_px, &b_gsf_px);
   _tree->SetBranchAddress("gsf_py", &gsf_py, &b_gsf_py);
   _tree->SetBranchAddress("gsf_pz", &gsf_pz, &b_gsf_pz);
   _tree->SetBranchAddress("gsf_caloEnergy", &gsf_caloEnergy, &b_gsf_caloEnergy);
   _tree->SetBranchAddress("gsf_deltaEtaSuperClusterTrackAtVtx", &gsf_deltaEtaSuperClusterTrackAtVtx, &b_gsf_deltaEtaSuperClusterTrackAtVtx);
   _tree->SetBranchAddress("gsf_deltaPhiSuperClusterTrackAtVtx", &gsf_deltaPhiSuperClusterTrackAtVtx, &b_gsf_deltaPhiSuperClusterTrackAtVtx);
   _tree->SetBranchAddress("gsf_hadronicOverEm", &gsf_hadronicOverEm, &b_gsf_hadronicOverEm);
   _tree->SetBranchAddress("gsf_hcalDepth1OverEcal", &gsf_hcalDepth1OverEcal, &b_gsf_hcalDepth1OverEcal);
   _tree->SetBranchAddress("gsf_hcalDepth2OverEcal", &gsf_hcalDepth2OverEcal, &b_gsf_hcalDepth2OverEcal);
   _tree->SetBranchAddress("gsf_dr03TkSumPt", &gsf_dr03TkSumPt, &b_gsf_dr03TkSumPt);
   _tree->SetBranchAddress("gsf_heepTrkPtIso", &gsf_heepTrkPtIso, &b_gsf_heepTrkPtIso);
   _tree->SetBranchAddress("gsf_dr03EcalRecHitSumEt", &gsf_dr03EcalRecHitSumEt, &b_gsf_dr03EcalRecHitSumEt);
   _tree->SetBranchAddress("gsf_dr03HcalDepth1TowerSumEt", &gsf_dr03HcalDepth1TowerSumEt, &b_gsf_dr03HcalDepth1TowerSumEt);
   _tree->SetBranchAddress("gsf_dr03HcalDepth2TowerSumEt", &gsf_dr03HcalDepth2TowerSumEt, &b_gsf_dr03HcalDepth2TowerSumEt);
   _tree->SetBranchAddress("gsf_charge", &gsf_charge, &b_gsf_charge);
   _tree->SetBranchAddress("gsf_sigmaIetaIeta", &gsf_sigmaIetaIeta, &b_gsf_sigmaIetaIeta);
   _tree->SetBranchAddress("gsf_ecaldrivenSeed", &gsf_ecaldrivenSeed, &b_gsf_ecaldrivenSeed);
   _tree->SetBranchAddress("gsf_trackerdrivenSeed", &gsf_trackerdrivenSeed, &b_gsf_trackerdrivenSeed);
   _tree->SetBranchAddress("gsf_isEB", &gsf_isEB, &b_gsf_isEB);
   _tree->SetBranchAddress("gsf_isEE", &gsf_isEE, &b_gsf_isEE);
   _tree->SetBranchAddress("gsf_passConversionVeto", &gsf_passConversionVeto, &b_gsf_passConversionVeto);
   _tree->SetBranchAddress("gsf_VID_cutBasedElectronID_Fall17_94X_V1_loose", &gsf_VID_cutBasedElectronID_Fall17_94X_V1_loose, &b_gsf_VID_cutBasedElectronID_Fall17_94X_V1_loose);
   _tree->SetBranchAddress("gsf_VID_cutBasedElectronID_Fall17_94X_V1_medium", &gsf_VID_cutBasedElectronID_Fall17_94X_V1_medium, &b_gsf_VID_cutBasedElectronID_Fall17_94X_V1_medium);
   _tree->SetBranchAddress("gsf_VID_cutBasedElectronID_Fall17_94X_V1_tight", &gsf_VID_cutBasedElectronID_Fall17_94X_V1_tight, &b_gsf_VID_cutBasedElectronID_Fall17_94X_V1_tight);
   _tree->SetBranchAddress("gsf_VID_cutBasedElectronID_Fall17_94X_V1_veto", &gsf_VID_cutBasedElectronID_Fall17_94X_V1_veto, &b_gsf_VID_cutBasedElectronID_Fall17_94X_V1_veto);
   _tree->SetBranchAddress("gsf_VID_cutBasedElectronID_Summer16_80X_V1_loose", &gsf_VID_cutBasedElectronID_Summer16_80X_V1_loose, &b_gsf_VID_cutBasedElectronID_Summer16_80X_V1_loose);
   _tree->SetBranchAddress("gsf_VID_cutBasedElectronID_Summer16_80X_V1_medium", &gsf_VID_cutBasedElectronID_Summer16_80X_V1_medium, &b_gsf_VID_cutBasedElectronID_Summer16_80X_V1_medium);
   _tree->SetBranchAddress("gsf_VID_cutBasedElectronID_Summer16_80X_V1_tight", &gsf_VID_cutBasedElectronID_Summer16_80X_V1_tight, &b_gsf_VID_cutBasedElectronID_Summer16_80X_V1_tight);
   _tree->SetBranchAddress("gsf_VID_cutBasedElectronID_Summer16_80X_V1_veto", &gsf_VID_cutBasedElectronID_Summer16_80X_V1_veto, &b_gsf_VID_cutBasedElectronID_Summer16_80X_V1_veto);
   _tree->SetBranchAddress("gsf_VID_heepElectronID_HEEPV70", &gsf_VID_heepElectronID_HEEPV70, &b_gsf_VID_heepElectronID_HEEPV70);
   _tree->SetBranchAddress("gsf_VID_mvaEleID_Fall17_iso_V1_wp80", &gsf_VID_mvaEleID_Fall17_iso_V1_wp80, &b_gsf_VID_mvaEleID_Fall17_iso_V1_wp80);
   _tree->SetBranchAddress("gsf_VID_mvaEleID_Fall17_iso_V1_wp90", &gsf_VID_mvaEleID_Fall17_iso_V1_wp90, &b_gsf_VID_mvaEleID_Fall17_iso_V1_wp90);
   _tree->SetBranchAddress("gsf_VID_mvaEleID_Fall17_iso_V1_wpLoose", &gsf_VID_mvaEleID_Fall17_iso_V1_wpLoose, &b_gsf_VID_mvaEleID_Fall17_iso_V1_wpLoose);
   _tree->SetBranchAddress("gsf_VID_mvaEleID_Fall17_noIso_V1_wp80", &gsf_VID_mvaEleID_Fall17_noIso_V1_wp80, &b_gsf_VID_mvaEleID_Fall17_noIso_V1_wp80);
   _tree->SetBranchAddress("gsf_VID_mvaEleID_Fall17_noIso_V1_wp90", &gsf_VID_mvaEleID_Fall17_noIso_V1_wp90, &b_gsf_VID_mvaEleID_Fall17_noIso_V1_wp90);
   _tree->SetBranchAddress("gsf_VID_mvaEleID_Fall17_noIso_V1_wpLoose", &gsf_VID_mvaEleID_Fall17_noIso_V1_wpLoose, &b_gsf_VID_mvaEleID_Fall17_noIso_V1_wpLoose);
   _tree->SetBranchAddress("gsf_VID_mvaEleID_Spring16_GeneralPurpose_V1_wp80", &gsf_VID_mvaEleID_Spring16_GeneralPurpose_V1_wp80, &b_gsf_VID_mvaEleID_Spring16_GeneralPurpose_V1_wp80);
   _tree->SetBranchAddress("gsf_VID_mvaEleID_Spring16_GeneralPurpose_V1_wp90", &gsf_VID_mvaEleID_Spring16_GeneralPurpose_V1_wp90, &b_gsf_VID_mvaEleID_Spring16_GeneralPurpose_V1_wp90);
   _tree->SetBranchAddress("gsf_VID_mvaEleID_Spring16_HZZ_V1_wpLoose", &gsf_VID_mvaEleID_Spring16_HZZ_V1_wpLoose, &b_gsf_VID_mvaEleID_Spring16_HZZ_V1_wpLoose);
   _tree->SetBranchAddress("gsf_deltaEtaSeedClusterTrackAtCalo", &gsf_deltaEtaSeedClusterTrackAtCalo, &b_gsf_deltaEtaSeedClusterTrackAtCalo);
   _tree->SetBranchAddress("gsf_deltaPhiSeedClusterTrackAtCalo", &gsf_deltaPhiSeedClusterTrackAtCalo, &b_gsf_deltaPhiSeedClusterTrackAtCalo);
   _tree->SetBranchAddress("gsf_ecalEnergy", &gsf_ecalEnergy, &b_gsf_ecalEnergy);
   _tree->SetBranchAddress("gsf_eSuperClusterOverP", &gsf_eSuperClusterOverP, &b_gsf_eSuperClusterOverP);
   _tree->SetBranchAddress("gsf_dxy", &gsf_dxy, &b_gsf_dxy);
   _tree->SetBranchAddress("gsf_dxy_beamSpot", &gsf_dxy_beamSpot, &b_gsf_dxy_beamSpot);
   _tree->SetBranchAddress("gsf_dxy_firstPVtx", &gsf_dxy_firstPVtx, &b_gsf_dxy_firstPVtx);
   _tree->SetBranchAddress("gsf_dxyError", &gsf_dxyError, &b_gsf_dxyError);
   _tree->SetBranchAddress("gsf_dz", &gsf_dz, &b_gsf_dz);
   _tree->SetBranchAddress("gsf_dz_beamSpot", &gsf_dz_beamSpot, &b_gsf_dz_beamSpot);
   _tree->SetBranchAddress("gsf_dz_firstPVtx", &gsf_dz_firstPVtx, &b_gsf_dz_firstPVtx);
   _tree->SetBranchAddress("gsf_dzError", &gsf_dzError, &b_gsf_dzError);
   _tree->SetBranchAddress("gsf_vz", &gsf_vz, &b_gsf_vz);
   _tree->SetBranchAddress("gsf_numberOfValidHits", &gsf_numberOfValidHits, &b_gsf_numberOfValidHits);
   _tree->SetBranchAddress("gsf_nLostInnerHits", &gsf_nLostInnerHits, &b_gsf_nLostInnerHits);
   _tree->SetBranchAddress("gsf_nLostOuterHits", &gsf_nLostOuterHits, &b_gsf_nLostOuterHits);
   _tree->SetBranchAddress("gsf_convFlags", &gsf_convFlags, &b_gsf_convFlags);
   _tree->SetBranchAddress("gsf_convDist", &gsf_convDist, &b_gsf_convDist);
   _tree->SetBranchAddress("gsf_convDcot", &gsf_convDcot, &b_gsf_convDcot);
   _tree->SetBranchAddress("gsf_convRadius", &gsf_convRadius, &b_gsf_convRadius);
   _tree->SetBranchAddress("gsf_fBrem", &gsf_fBrem, &b_gsf_fBrem);
   _tree->SetBranchAddress("gsf_e1x5", &gsf_e1x5, &b_gsf_e1x5);
   _tree->SetBranchAddress("gsf_e2x5Max", &gsf_e2x5Max, &b_gsf_e2x5Max);
   _tree->SetBranchAddress("gsf_e5x5", &gsf_e5x5, &b_gsf_e5x5);
   _tree->SetBranchAddress("gsf_r9", &gsf_r9, &b_gsf_r9);
   _tree->SetBranchAddress("gsf_deltaEtaSeedClusterTrackAtVtx", &gsf_deltaEtaSeedClusterTrackAtVtx, &b_gsf_deltaEtaSeedClusterTrackAtVtx);
   _tree->SetBranchAddress("gsf_effArea", &gsf_effArea, &b_gsf_effArea);
   _tree->SetBranchAddress("gsf_sumChargedHadronPt", &gsf_sumChargedHadronPt, &b_gsf_sumChargedHadronPt);
   _tree->SetBranchAddress("gsf_sumNeutralHadronEt", &gsf_sumNeutralHadronEt, &b_gsf_sumNeutralHadronEt);
   _tree->SetBranchAddress("gsf_sumPhotonEt", &gsf_sumPhotonEt, &b_gsf_sumPhotonEt);
   _tree->SetBranchAddress("gsf_sumPUPt", &gsf_sumPUPt, &b_gsf_sumPUPt);
   _tree->SetBranchAddress("gsf_ooEmooP", &gsf_ooEmooP, &b_gsf_ooEmooP);
   _tree->SetBranchAddress("gsf_hitsinfo", &gsf_hitsinfo, &b_gsf_hitsinfo);
   _tree->SetBranchAddress("gsf_pixelMatch_dPhi1", &gsf_pixelMatch_dPhi1, &b_gsf_pixelMatch_dPhi1);
   _tree->SetBranchAddress("gsf_pixelMatch_dPhi2", &gsf_pixelMatch_dPhi2, &b_gsf_pixelMatch_dPhi2);
   _tree->SetBranchAddress("gsf_pixelMatch_dRz1", &gsf_pixelMatch_dRz1, &b_gsf_pixelMatch_dRz1);
   _tree->SetBranchAddress("gsf_pixelMatch_dRz2", &gsf_pixelMatch_dRz2, &b_gsf_pixelMatch_dRz2);
   _tree->SetBranchAddress("gsf_pixelMatch_subDetector1", &gsf_pixelMatch_subDetector1, &b_gsf_pixelMatch_subDetector1);
   _tree->SetBranchAddress("gsf_pixelMatch_subDetector2", &gsf_pixelMatch_subDetector2, &b_gsf_pixelMatch_subDetector2);
   _tree->SetBranchAddress("gsf_sc_energy", &gsf_sc_energy, &b_gsf_sc_energy);
   _tree->SetBranchAddress("gsf_sc_seed_eta", &gsf_sc_seed_eta, &b_gsf_sc_seed_eta);
   _tree->SetBranchAddress("gsf_sc_eta", &gsf_sc_eta, &b_gsf_sc_eta);
   _tree->SetBranchAddress("gsf_sc_etacorr", &gsf_sc_etacorr, &b_gsf_sc_etacorr);
   _tree->SetBranchAddress("gsf_sc_theta", &gsf_sc_theta, &b_gsf_sc_theta);
   _tree->SetBranchAddress("gsf_sc_thetacorr", &gsf_sc_thetacorr, &b_gsf_sc_thetacorr);
   _tree->SetBranchAddress("gsf_sc_et", &gsf_sc_et, &b_gsf_sc_et);
   _tree->SetBranchAddress("gsf_sc_phi", &gsf_sc_phi, &b_gsf_sc_phi);
   _tree->SetBranchAddress("gsf_sc_px", &gsf_sc_px, &b_gsf_sc_px);
   _tree->SetBranchAddress("gsf_sc_py", &gsf_sc_py, &b_gsf_sc_py);
   _tree->SetBranchAddress("gsf_sc_pz", &gsf_sc_pz, &b_gsf_sc_pz);
   _tree->SetBranchAddress("gsf_sc_x", &gsf_sc_x, &b_gsf_sc_x);
   _tree->SetBranchAddress("gsf_sc_y", &gsf_sc_y, &b_gsf_sc_y);
   _tree->SetBranchAddress("gsf_sc_z", &gsf_sc_z, &b_gsf_sc_z);
   _tree->SetBranchAddress("gsf_sc_phiWidth", &gsf_sc_phiWidth, &b_gsf_sc_phiWidth);
   _tree->SetBranchAddress("gsf_sc_etaWidth", &gsf_sc_etaWidth, &b_gsf_sc_etaWidth);
   _tree->SetBranchAddress("gsf_sc_seed_rawId", &gsf_sc_seed_rawId, &b_gsf_sc_seed_rawId);
   _tree->SetBranchAddress("gsf_sc_seed_ieta", &gsf_sc_seed_ieta, &b_gsf_sc_seed_ieta);
   _tree->SetBranchAddress("gsf_sc_seed_iphi", &gsf_sc_seed_iphi, &b_gsf_sc_seed_iphi);
   _tree->SetBranchAddress("gsf_sc_seed_kHasSwitchToGain6", &gsf_sc_seed_kHasSwitchToGain6, &b_gsf_sc_seed_kHasSwitchToGain6);
   _tree->SetBranchAddress("gsf_sc_seed_kHasSwitchToGain1", &gsf_sc_seed_kHasSwitchToGain1, &b_gsf_sc_seed_kHasSwitchToGain1);
   _tree->SetBranchAddress("gsf_swissCross", &gsf_swissCross, &b_gsf_swissCross);
   _tree->SetBranchAddress("gsf_sc_rawEnergy", &gsf_sc_rawEnergy, &b_gsf_sc_rawEnergy);
   _tree->SetBranchAddress("gsf_sc_preshowerEnergy", &gsf_sc_preshowerEnergy, &b_gsf_sc_preshowerEnergy);
   _tree->SetBranchAddress("gsf_sc_lazyTools_e2x5Right", &gsf_sc_lazyTools_e2x5Right, &b_gsf_sc_lazyTools_e2x5Right);
   _tree->SetBranchAddress("gsf_sc_lazyTools_e2x5Left", &gsf_sc_lazyTools_e2x5Left, &b_gsf_sc_lazyTools_e2x5Left);
   _tree->SetBranchAddress("gsf_sc_lazyTools_e2x5Top", &gsf_sc_lazyTools_e2x5Top, &b_gsf_sc_lazyTools_e2x5Top);
   _tree->SetBranchAddress("gsf_sc_lazyTools_e2x5Bottom", &gsf_sc_lazyTools_e2x5Bottom, &b_gsf_sc_lazyTools_e2x5Bottom);
   _tree->SetBranchAddress("gsf_sc_lazyTools_eMax", &gsf_sc_lazyTools_eMax, &b_gsf_sc_lazyTools_eMax);
   _tree->SetBranchAddress("gsf_sc_lazyTools_e2nd", &gsf_sc_lazyTools_e2nd, &b_gsf_sc_lazyTools_e2nd);
   _tree->SetBranchAddress("gsf_sc_lazyTools_eRight", &gsf_sc_lazyTools_eRight, &b_gsf_sc_lazyTools_eRight);
   _tree->SetBranchAddress("gsf_sc_lazyTools_eLeft", &gsf_sc_lazyTools_eLeft, &b_gsf_sc_lazyTools_eLeft);
   _tree->SetBranchAddress("gsf_sc_lazyTools_eTop", &gsf_sc_lazyTools_eTop, &b_gsf_sc_lazyTools_eTop);
   _tree->SetBranchAddress("gsf_sc_lazyTools_eBottom", &gsf_sc_lazyTools_eBottom, &b_gsf_sc_lazyTools_eBottom);
   _tree->SetBranchAddress("gsf_sc_lazyTools_e2x2", &gsf_sc_lazyTools_e2x2, &b_gsf_sc_lazyTools_e2x2);
   _tree->SetBranchAddress("gsf_sc_lazyTools_e3x3", &gsf_sc_lazyTools_e3x3, &b_gsf_sc_lazyTools_e3x3);
   _tree->SetBranchAddress("gsf_sc_lazyTools_e4x4", &gsf_sc_lazyTools_e4x4, &b_gsf_sc_lazyTools_e4x4);
   _tree->SetBranchAddress("gsf_sc_lazyTools_e5x5", &gsf_sc_lazyTools_e5x5, &b_gsf_sc_lazyTools_e5x5);
   _tree->SetBranchAddress("gsf_sc_lazyTools_e1x3", &gsf_sc_lazyTools_e1x3, &b_gsf_sc_lazyTools_e1x3);
   _tree->SetBranchAddress("gsf_sc_lazyTools_e3x1", &gsf_sc_lazyTools_e3x1, &b_gsf_sc_lazyTools_e3x1);
   _tree->SetBranchAddress("gsf_sc_lazyTools_e1x5", &gsf_sc_lazyTools_e1x5, &b_gsf_sc_lazyTools_e1x5);
   _tree->SetBranchAddress("gsf_sc_lazyTools_e5x1", &gsf_sc_lazyTools_e5x1, &b_gsf_sc_lazyTools_e5x1);
   _tree->SetBranchAddress("gsf_sc_lazyTools_eshitsixix", &gsf_sc_lazyTools_eshitsixix, &b_gsf_sc_lazyTools_eshitsixix);
   _tree->SetBranchAddress("gsf_sc_lazyTools_eshitsiyiy", &gsf_sc_lazyTools_eshitsiyiy, &b_gsf_sc_lazyTools_eshitsiyiy);
   _tree->SetBranchAddress("gsf_sc_lazyTools_eseffsixix", &gsf_sc_lazyTools_eseffsixix, &b_gsf_sc_lazyTools_eseffsixix);
   _tree->SetBranchAddress("gsf_sc_lazyTools_eseffsiyiy", &gsf_sc_lazyTools_eseffsiyiy, &b_gsf_sc_lazyTools_eseffsiyiy);
   _tree->SetBranchAddress("gsf_sc_lazyTools_eseffsirir", &gsf_sc_lazyTools_eseffsirir, &b_gsf_sc_lazyTools_eseffsirir);
   _tree->SetBranchAddress("gsf_sc_lazyTools_BasicClusterSeedTime", &gsf_sc_lazyTools_BasicClusterSeedTime, &b_gsf_sc_lazyTools_BasicClusterSeedTime);
   _tree->SetBranchAddress("EHits_isSaturated", &EHits_isSaturated, &b_EHits_isSaturated);
   _tree->SetBranchAddress("EBHits_rawId", &EBHits_rawId, &b_EBHits_rawId);
   _tree->SetBranchAddress("EBHits_iRechit", &EBHits_iRechit, &b_EBHits_iRechit);
   _tree->SetBranchAddress("EBHits_ieta", &EBHits_ieta, &b_EBHits_ieta);
   _tree->SetBranchAddress("EBHits_iphi", &EBHits_iphi, &b_EBHits_iphi);
   _tree->SetBranchAddress("EBHits_RecoFlag", &EBHits_RecoFlag, &b_EBHits_RecoFlag);
   _tree->SetBranchAddress("EBHits_kSaturated", &EBHits_kSaturated, &b_EBHits_kSaturated);
   _tree->SetBranchAddress("EBHits_kLeadingEdgeRecovered", &EBHits_kLeadingEdgeRecovered, &b_EBHits_kLeadingEdgeRecovered);
   _tree->SetBranchAddress("EBHits_kNeighboursRecovered", &EBHits_kNeighboursRecovered, &b_EBHits_kNeighboursRecovered);
   _tree->SetBranchAddress("EBHits_kWeird", &EBHits_kWeird, &b_EBHits_kWeird);
   _tree->SetBranchAddress("EBHits_energy", &EBHits_energy, &b_EBHits_energy);
   _tree->SetBranchAddress("EEHits_rawId", &EEHits_rawId, &b_EEHits_rawId);
   _tree->SetBranchAddress("EEHits_iRechit", &EEHits_iRechit, &b_EEHits_iRechit);
   _tree->SetBranchAddress("EEHits_ieta", &EEHits_ieta, &b_EEHits_ieta);
   _tree->SetBranchAddress("EEHits_iphi", &EEHits_iphi, &b_EEHits_iphi);
   _tree->SetBranchAddress("EEHits_RecoFlag", &EEHits_RecoFlag, &b_EEHits_RecoFlag);
   _tree->SetBranchAddress("EEHits_kSaturated", &EEHits_kSaturated, &b_EEHits_kSaturated);
   _tree->SetBranchAddress("EEHits_kLeadingEdgeRecovered", &EEHits_kLeadingEdgeRecovered, &b_EEHits_kLeadingEdgeRecovered);
   _tree->SetBranchAddress("EEHits_kNeighboursRecovered", &EEHits_kNeighboursRecovered, &b_EEHits_kNeighboursRecovered);
   _tree->SetBranchAddress("EEHits_kWeird", &EEHits_kWeird, &b_EEHits_kWeird);
   _tree->SetBranchAddress("EEHits_energy", &EEHits_energy, &b_EEHits_energy);






   if(_tree->GetListOfBranches()->FindObject("gsf_mc_bestDR")) {

   _tree->SetBranchAddress("gsf_mc_bestDR", &gsf_mc_bestDR, &b_gsf_mc_bestDR);
   _tree->SetBranchAddress("gsf_mc_index", &gsf_mc_index, &b_gsf_mc_index);
   _tree->SetBranchAddress("gsf_mc_ERatio", &gsf_mc_ERatio, &b_gsf_mc_ERatio);

   }


   _tree->SetBranchAddress("mu_n", &mu_n, &b_mu_n);
   _tree->SetBranchAddress("mu_gt_qoverp", &mu_gt_qoverp, &b_mu_gt_qoverp);
   _tree->SetBranchAddress("mu_gt_charge", &mu_gt_charge, &b_mu_gt_charge);
   _tree->SetBranchAddress("mu_gt_pt", &mu_gt_pt, &b_mu_gt_pt);
   _tree->SetBranchAddress("mu_gt_eta", &mu_gt_eta, &b_mu_gt_eta);
   _tree->SetBranchAddress("mu_gt_phi", &mu_gt_phi, &b_mu_gt_phi);
   _tree->SetBranchAddress("mu_gt_p", &mu_gt_p, &b_mu_gt_p);
   _tree->SetBranchAddress("mu_gt_px", &mu_gt_px, &b_mu_gt_px);
   _tree->SetBranchAddress("mu_gt_py", &mu_gt_py, &b_mu_gt_py);
   _tree->SetBranchAddress("mu_gt_pz", &mu_gt_pz, &b_mu_gt_pz);
   _tree->SetBranchAddress("mu_gt_theta", &mu_gt_theta, &b_mu_gt_theta);
   _tree->SetBranchAddress("mu_gt_lambda", &mu_gt_lambda, &b_mu_gt_lambda);
   _tree->SetBranchAddress("mu_gt_d0", &mu_gt_d0, &b_mu_gt_d0);
   _tree->SetBranchAddress("mu_gt_dz", &mu_gt_dz, &b_mu_gt_dz);
   _tree->SetBranchAddress("mu_gt_dz_beamspot", &mu_gt_dz_beamspot, &b_mu_gt_dz_beamspot);
   _tree->SetBranchAddress("mu_gt_dz_firstPVtx", &mu_gt_dz_firstPVtx, &b_mu_gt_dz_firstPVtx);
   _tree->SetBranchAddress("mu_gt_dxy", &mu_gt_dxy, &b_mu_gt_dxy);
   _tree->SetBranchAddress("mu_gt_dxy_beamspot", &mu_gt_dxy_beamspot, &b_mu_gt_dxy_beamspot);
   _tree->SetBranchAddress("mu_gt_dxy_firstPVtx", &mu_gt_dxy_firstPVtx, &b_mu_gt_dxy_firstPVtx);
   _tree->SetBranchAddress("mu_gt_dsz", &mu_gt_dsz, &b_mu_gt_dsz);
   _tree->SetBranchAddress("mu_gt_vx", &mu_gt_vx, &b_mu_gt_vx);
   _tree->SetBranchAddress("mu_gt_vy", &mu_gt_vy, &b_mu_gt_vy);
   _tree->SetBranchAddress("mu_gt_vz", &mu_gt_vz, &b_mu_gt_vz);
   _tree->SetBranchAddress("mu_gt_qoverpError", &mu_gt_qoverpError, &b_mu_gt_qoverpError);
   _tree->SetBranchAddress("mu_gt_ptError", &mu_gt_ptError, &b_mu_gt_ptError);
   _tree->SetBranchAddress("mu_gt_thetaError", &mu_gt_thetaError, &b_mu_gt_thetaError);
   _tree->SetBranchAddress("mu_gt_lambdaError", &mu_gt_lambdaError, &b_mu_gt_lambdaError);
   _tree->SetBranchAddress("mu_gt_phiError", &mu_gt_phiError, &b_mu_gt_phiError);
   _tree->SetBranchAddress("mu_gt_dxyError", &mu_gt_dxyError, &b_mu_gt_dxyError);
   _tree->SetBranchAddress("mu_gt_d0Error", &mu_gt_d0Error, &b_mu_gt_d0Error);
   _tree->SetBranchAddress("mu_gt_dszError", &mu_gt_dszError, &b_mu_gt_dszError);
   _tree->SetBranchAddress("mu_gt_dzError", &mu_gt_dzError, &b_mu_gt_dzError);
   _tree->SetBranchAddress("mu_gt_etaError", &mu_gt_etaError, &b_mu_gt_etaError);
   _tree->SetBranchAddress("mu_gt_chi2", &mu_gt_chi2, &b_mu_gt_chi2);
   _tree->SetBranchAddress("mu_gt_ndof", &mu_gt_ndof, &b_mu_gt_ndof);
   _tree->SetBranchAddress("mu_gt_normalizedChi2", &mu_gt_normalizedChi2, &b_mu_gt_normalizedChi2);
   _tree->SetBranchAddress("mu_ot_qoverp", &mu_ot_qoverp, &b_mu_ot_qoverp);
   _tree->SetBranchAddress("mu_ot_charge", &mu_ot_charge, &b_mu_ot_charge);
   _tree->SetBranchAddress("mu_ot_pt", &mu_ot_pt, &b_mu_ot_pt);
   _tree->SetBranchAddress("mu_ot_eta", &mu_ot_eta, &b_mu_ot_eta);
   _tree->SetBranchAddress("mu_ot_phi", &mu_ot_phi, &b_mu_ot_phi);
   _tree->SetBranchAddress("mu_ot_p", &mu_ot_p, &b_mu_ot_p);
   _tree->SetBranchAddress("mu_ot_px", &mu_ot_px, &b_mu_ot_px);
   _tree->SetBranchAddress("mu_ot_py", &mu_ot_py, &b_mu_ot_py);
   _tree->SetBranchAddress("mu_ot_pz", &mu_ot_pz, &b_mu_ot_pz);
   _tree->SetBranchAddress("mu_ot_theta", &mu_ot_theta, &b_mu_ot_theta);
   _tree->SetBranchAddress("mu_ot_lambda", &mu_ot_lambda, &b_mu_ot_lambda);
   _tree->SetBranchAddress("mu_ot_d0", &mu_ot_d0, &b_mu_ot_d0);
   _tree->SetBranchAddress("mu_ot_dz", &mu_ot_dz, &b_mu_ot_dz);
   _tree->SetBranchAddress("mu_ot_dz_beamspot", &mu_ot_dz_beamspot, &b_mu_ot_dz_beamspot);
   _tree->SetBranchAddress("mu_ot_dz_firstPVtx", &mu_ot_dz_firstPVtx, &b_mu_ot_dz_firstPVtx);
   _tree->SetBranchAddress("mu_ot_dxy", &mu_ot_dxy, &b_mu_ot_dxy);
   _tree->SetBranchAddress("mu_ot_dxy_beamspot", &mu_ot_dxy_beamspot, &b_mu_ot_dxy_beamspot);
   _tree->SetBranchAddress("mu_ot_dxy_firstPVtx", &mu_ot_dxy_firstPVtx, &b_mu_ot_dxy_firstPVtx);
   _tree->SetBranchAddress("mu_ot_dsz", &mu_ot_dsz, &b_mu_ot_dsz);
   _tree->SetBranchAddress("mu_ot_vx", &mu_ot_vx, &b_mu_ot_vx);
   _tree->SetBranchAddress("mu_ot_vy", &mu_ot_vy, &b_mu_ot_vy);
   _tree->SetBranchAddress("mu_ot_vz", &mu_ot_vz, &b_mu_ot_vz);
   _tree->SetBranchAddress("mu_ot_qoverpError", &mu_ot_qoverpError, &b_mu_ot_qoverpError);
   _tree->SetBranchAddress("mu_ot_ptError", &mu_ot_ptError, &b_mu_ot_ptError);
   _tree->SetBranchAddress("mu_ot_thetaError", &mu_ot_thetaError, &b_mu_ot_thetaError);
   _tree->SetBranchAddress("mu_ot_lambdaError", &mu_ot_lambdaError, &b_mu_ot_lambdaError);
   _tree->SetBranchAddress("mu_ot_phiError", &mu_ot_phiError, &b_mu_ot_phiError);
   _tree->SetBranchAddress("mu_ot_dxyError", &mu_ot_dxyError, &b_mu_ot_dxyError);
   _tree->SetBranchAddress("mu_ot_d0Error", &mu_ot_d0Error, &b_mu_ot_d0Error);
   _tree->SetBranchAddress("mu_ot_dszError", &mu_ot_dszError, &b_mu_ot_dszError);
   _tree->SetBranchAddress("mu_ot_dzError", &mu_ot_dzError, &b_mu_ot_dzError);
   _tree->SetBranchAddress("mu_ot_etaError", &mu_ot_etaError, &b_mu_ot_etaError);
   _tree->SetBranchAddress("mu_ot_chi2", &mu_ot_chi2, &b_mu_ot_chi2);
   _tree->SetBranchAddress("mu_ot_ndof", &mu_ot_ndof, &b_mu_ot_ndof);
   _tree->SetBranchAddress("mu_ot_normalizedChi2", &mu_ot_normalizedChi2, &b_mu_ot_normalizedChi2);
   _tree->SetBranchAddress("mu_it_qoverp", &mu_it_qoverp, &b_mu_it_qoverp);
   _tree->SetBranchAddress("mu_it_charge", &mu_it_charge, &b_mu_it_charge);
   _tree->SetBranchAddress("mu_it_pt", &mu_it_pt, &b_mu_it_pt);
   _tree->SetBranchAddress("mu_it_eta", &mu_it_eta, &b_mu_it_eta);
   _tree->SetBranchAddress("mu_it_phi", &mu_it_phi, &b_mu_it_phi);
   _tree->SetBranchAddress("mu_it_p", &mu_it_p, &b_mu_it_p);
   _tree->SetBranchAddress("mu_it_px", &mu_it_px, &b_mu_it_px);
   _tree->SetBranchAddress("mu_it_py", &mu_it_py, &b_mu_it_py);
   _tree->SetBranchAddress("mu_it_pz", &mu_it_pz, &b_mu_it_pz);
   _tree->SetBranchAddress("mu_it_theta", &mu_it_theta, &b_mu_it_theta);
   _tree->SetBranchAddress("mu_it_lambda", &mu_it_lambda, &b_mu_it_lambda);
   _tree->SetBranchAddress("mu_it_d0", &mu_it_d0, &b_mu_it_d0);
   _tree->SetBranchAddress("mu_it_dz", &mu_it_dz, &b_mu_it_dz);
   _tree->SetBranchAddress("mu_it_dz_beamspot", &mu_it_dz_beamspot, &b_mu_it_dz_beamspot);
   _tree->SetBranchAddress("mu_it_dz_firstPVtx", &mu_it_dz_firstPVtx, &b_mu_it_dz_firstPVtx);
   _tree->SetBranchAddress("mu_it_dxy", &mu_it_dxy, &b_mu_it_dxy);
   _tree->SetBranchAddress("mu_it_dxy_beamspot", &mu_it_dxy_beamspot, &b_mu_it_dxy_beamspot);
   _tree->SetBranchAddress("mu_it_dxy_firstPVtx", &mu_it_dxy_firstPVtx, &b_mu_it_dxy_firstPVtx);
   _tree->SetBranchAddress("mu_it_dsz", &mu_it_dsz, &b_mu_it_dsz);
   _tree->SetBranchAddress("mu_it_vx", &mu_it_vx, &b_mu_it_vx);
   _tree->SetBranchAddress("mu_it_vy", &mu_it_vy, &b_mu_it_vy);
   _tree->SetBranchAddress("mu_it_vz", &mu_it_vz, &b_mu_it_vz);
   _tree->SetBranchAddress("mu_it_qoverpError", &mu_it_qoverpError, &b_mu_it_qoverpError);
   _tree->SetBranchAddress("mu_it_ptError", &mu_it_ptError, &b_mu_it_ptError);
   _tree->SetBranchAddress("mu_it_thetaError", &mu_it_thetaError, &b_mu_it_thetaError);
   _tree->SetBranchAddress("mu_it_lambdaError", &mu_it_lambdaError, &b_mu_it_lambdaError);
   _tree->SetBranchAddress("mu_it_phiError", &mu_it_phiError, &b_mu_it_phiError);
   _tree->SetBranchAddress("mu_it_dxyError", &mu_it_dxyError, &b_mu_it_dxyError);
   _tree->SetBranchAddress("mu_it_d0Error", &mu_it_d0Error, &b_mu_it_d0Error);
   _tree->SetBranchAddress("mu_it_dszError", &mu_it_dszError, &b_mu_it_dszError);
   _tree->SetBranchAddress("mu_it_dzError", &mu_it_dzError, &b_mu_it_dzError);
   _tree->SetBranchAddress("mu_it_etaError", &mu_it_etaError, &b_mu_it_etaError);
   _tree->SetBranchAddress("mu_it_chi2", &mu_it_chi2, &b_mu_it_chi2);
   _tree->SetBranchAddress("mu_it_ndof", &mu_it_ndof, &b_mu_it_ndof);
   _tree->SetBranchAddress("mu_it_normalizedChi2", &mu_it_normalizedChi2, &b_mu_it_normalizedChi2);
   _tree->SetBranchAddress("mu_ibt_qoverp", &mu_ibt_qoverp, &b_mu_ibt_qoverp);
   _tree->SetBranchAddress("mu_ibt_charge", &mu_ibt_charge, &b_mu_ibt_charge);
   _tree->SetBranchAddress("mu_ibt_pt", &mu_ibt_pt, &b_mu_ibt_pt);
   _tree->SetBranchAddress("mu_ibt_eta", &mu_ibt_eta, &b_mu_ibt_eta);
   _tree->SetBranchAddress("mu_ibt_phi", &mu_ibt_phi, &b_mu_ibt_phi);
   _tree->SetBranchAddress("mu_ibt_p", &mu_ibt_p, &b_mu_ibt_p);
   _tree->SetBranchAddress("mu_ibt_px", &mu_ibt_px, &b_mu_ibt_px);
   _tree->SetBranchAddress("mu_ibt_py", &mu_ibt_py, &b_mu_ibt_py);
   _tree->SetBranchAddress("mu_ibt_pz", &mu_ibt_pz, &b_mu_ibt_pz);
   _tree->SetBranchAddress("mu_ibt_theta", &mu_ibt_theta, &b_mu_ibt_theta);
   _tree->SetBranchAddress("mu_ibt_lambda", &mu_ibt_lambda, &b_mu_ibt_lambda);
   _tree->SetBranchAddress("mu_ibt_d0", &mu_ibt_d0, &b_mu_ibt_d0);
   _tree->SetBranchAddress("mu_ibt_dz", &mu_ibt_dz, &b_mu_ibt_dz);
   _tree->SetBranchAddress("mu_ibt_dz_beamspot", &mu_ibt_dz_beamspot, &b_mu_ibt_dz_beamspot);
   _tree->SetBranchAddress("mu_ibt_dz_firstPVtx", &mu_ibt_dz_firstPVtx, &b_mu_ibt_dz_firstPVtx);
   _tree->SetBranchAddress("mu_ibt_dxy", &mu_ibt_dxy, &b_mu_ibt_dxy);
   _tree->SetBranchAddress("mu_ibt_dxy_beamspot", &mu_ibt_dxy_beamspot, &b_mu_ibt_dxy_beamspot);
   _tree->SetBranchAddress("mu_ibt_dxy_firstPVtx", &mu_ibt_dxy_firstPVtx, &b_mu_ibt_dxy_firstPVtx);
   _tree->SetBranchAddress("mu_ibt_dsz", &mu_ibt_dsz, &b_mu_ibt_dsz);
   _tree->SetBranchAddress("mu_ibt_vx", &mu_ibt_vx, &b_mu_ibt_vx);
   _tree->SetBranchAddress("mu_ibt_vy", &mu_ibt_vy, &b_mu_ibt_vy);
   _tree->SetBranchAddress("mu_ibt_vz", &mu_ibt_vz, &b_mu_ibt_vz);
   _tree->SetBranchAddress("mu_ibt_qoverpError", &mu_ibt_qoverpError, &b_mu_ibt_qoverpError);
   _tree->SetBranchAddress("mu_ibt_ptError", &mu_ibt_ptError, &b_mu_ibt_ptError);
   _tree->SetBranchAddress("mu_ibt_thetaError", &mu_ibt_thetaError, &b_mu_ibt_thetaError);
   _tree->SetBranchAddress("mu_ibt_lambdaError", &mu_ibt_lambdaError, &b_mu_ibt_lambdaError);
   _tree->SetBranchAddress("mu_ibt_phiError", &mu_ibt_phiError, &b_mu_ibt_phiError);
   _tree->SetBranchAddress("mu_ibt_dxyError", &mu_ibt_dxyError, &b_mu_ibt_dxyError);
   _tree->SetBranchAddress("mu_ibt_d0Error", &mu_ibt_d0Error, &b_mu_ibt_d0Error);
   _tree->SetBranchAddress("mu_ibt_dszError", &mu_ibt_dszError, &b_mu_ibt_dszError);
   _tree->SetBranchAddress("mu_ibt_dzError", &mu_ibt_dzError, &b_mu_ibt_dzError);
   _tree->SetBranchAddress("mu_ibt_etaError", &mu_ibt_etaError, &b_mu_ibt_etaError);
   _tree->SetBranchAddress("mu_ibt_chi2", &mu_ibt_chi2, &b_mu_ibt_chi2);
   _tree->SetBranchAddress("mu_ibt_ndof", &mu_ibt_ndof, &b_mu_ibt_ndof);
   _tree->SetBranchAddress("mu_ibt_normalizedChi2", &mu_ibt_normalizedChi2, &b_mu_ibt_normalizedChi2);
   _tree->SetBranchAddress("mu_isGlobalMuon", &mu_isGlobalMuon, &b_mu_isGlobalMuon);
   _tree->SetBranchAddress("mu_isStandAloneMuon", &mu_isStandAloneMuon, &b_mu_isStandAloneMuon);
   _tree->SetBranchAddress("mu_isTrackerMuon", &mu_isTrackerMuon, &b_mu_isTrackerMuon);
   _tree->SetBranchAddress("mu_isPFMuon", &mu_isPFMuon, &b_mu_isPFMuon);
   _tree->SetBranchAddress("mu_isPFIsolationValid", &mu_isPFIsolationValid, &b_mu_isPFIsolationValid);
   _tree->SetBranchAddress("mu_isGoodMuonTMLastStationLoose", &mu_isGoodMuonTMLastStationLoose, &b_mu_isGoodMuonTMLastStationLoose);
   _tree->SetBranchAddress("mu_isGoodMuonTMLastStationTight", &mu_isGoodMuonTMLastStationTight, &b_mu_isGoodMuonTMLastStationTight);
   _tree->SetBranchAddress("mu_isGoodMuonTM2DCompatibilityLoose", &mu_isGoodMuonTM2DCompatibilityLoose, &b_mu_isGoodMuonTM2DCompatibilityLoose);
   _tree->SetBranchAddress("mu_isGoodMuonTM2DCompatibilityTight", &mu_isGoodMuonTM2DCompatibilityTight, &b_mu_isGoodMuonTM2DCompatibilityTight);
   _tree->SetBranchAddress("mu_isGoodMuonTMOneStationLoose", &mu_isGoodMuonTMOneStationLoose, &b_mu_isGoodMuonTMOneStationLoose);
   _tree->SetBranchAddress("mu_isGoodMuonTMOneStationTight", &mu_isGoodMuonTMOneStationTight, &b_mu_isGoodMuonTMOneStationTight);
   _tree->SetBranchAddress("mu_isGoodMuonTMLastStationOptimizedLowPtLoose", &mu_isGoodMuonTMLastStationOptimizedLowPtLoose, &b_mu_isGoodMuonTMLastStationOptimizedLowPtLoose);
   _tree->SetBranchAddress("mu_isGoodMuonTMLastStationOptimizedLowPtTight", &mu_isGoodMuonTMLastStationOptimizedLowPtTight, &b_mu_isGoodMuonTMLastStationOptimizedLowPtTight);
   _tree->SetBranchAddress("mu_isTightMuon", &mu_isTightMuon, &b_mu_isTightMuon);
   _tree->SetBranchAddress("mu_isMediumMuon", &mu_isMediumMuon, &b_mu_isMediumMuon);
   _tree->SetBranchAddress("mu_isLooseMuon", &mu_isLooseMuon, &b_mu_isLooseMuon);
   _tree->SetBranchAddress("mu_isSoftMuon", &mu_isSoftMuon, &b_mu_isSoftMuon);
   _tree->SetBranchAddress("mu_isHighPtMuonOld", &mu_isHighPtMuonOld, &b_mu_isHighPtMuonOld);
   _tree->SetBranchAddress("mu_isHighPtMuon", &mu_isHighPtMuon, &b_mu_isHighPtMuon);
   _tree->SetBranchAddress("mu_isTrackerHighPtMuon", &mu_isTrackerHighPtMuon, &b_mu_isTrackerHighPtMuon);
   _tree->SetBranchAddress("mu_CutBasedIdLoose", &mu_CutBasedIdLoose, &b_mu_CutBasedIdLoose);
   _tree->SetBranchAddress("mu_CutBasedIdMedium", &mu_CutBasedIdMedium, &b_mu_CutBasedIdMedium);
   _tree->SetBranchAddress("mu_CutBasedIdMediumPrompt", &mu_CutBasedIdMediumPrompt, &b_mu_CutBasedIdMediumPrompt);
   _tree->SetBranchAddress("mu_CutBasedIdTight", &mu_CutBasedIdTight, &b_mu_CutBasedIdTight);
   _tree->SetBranchAddress("mu_CutBasedIdGlobalHighPt", &mu_CutBasedIdGlobalHighPt, &b_mu_CutBasedIdGlobalHighPt);
   _tree->SetBranchAddress("mu_CutBasedIdTrkHighPt", &mu_CutBasedIdTrkHighPt, &b_mu_CutBasedIdTrkHighPt);
   _tree->SetBranchAddress("mu_PFIsoVeryLoose", &mu_PFIsoVeryLoose, &b_mu_PFIsoVeryLoose);
   _tree->SetBranchAddress("mu_PFIsoLoose", &mu_PFIsoLoose, &b_mu_PFIsoLoose);
   _tree->SetBranchAddress("mu_PFIsoMedium", &mu_PFIsoMedium, &b_mu_PFIsoMedium);
   _tree->SetBranchAddress("mu_PFIsoTight", &mu_PFIsoTight, &b_mu_PFIsoTight);
   _tree->SetBranchAddress("mu_PFIsoVeryTight", &mu_PFIsoVeryTight, &b_mu_PFIsoVeryTight);
   _tree->SetBranchAddress("mu_PFIsoVeryVeryTight", &mu_PFIsoVeryVeryTight, &b_mu_PFIsoVeryVeryTight);
   _tree->SetBranchAddress("mu_TkIsoLoose", &mu_TkIsoLoose, &b_mu_TkIsoLoose);
   _tree->SetBranchAddress("mu_TkIsoTight", &mu_TkIsoTight, &b_mu_TkIsoTight);
   _tree->SetBranchAddress("mu_SoftCutBasedId", &mu_SoftCutBasedId, &b_mu_SoftCutBasedId);
   _tree->SetBranchAddress("mu_SoftMvaId", &mu_SoftMvaId, &b_mu_SoftMvaId);
   _tree->SetBranchAddress("mu_MvaLoose", &mu_MvaLoose, &b_mu_MvaLoose);
   _tree->SetBranchAddress("mu_MvaMedium", &mu_MvaMedium, &b_mu_MvaMedium);
   _tree->SetBranchAddress("mu_MvaTight", &mu_MvaTight, &b_mu_MvaTight);
   _tree->SetBranchAddress("mu_MiniIsoLoose", &mu_MiniIsoLoose, &b_mu_MiniIsoLoose);
   _tree->SetBranchAddress("mu_MiniIsoMedium", &mu_MiniIsoMedium, &b_mu_MiniIsoMedium);
   _tree->SetBranchAddress("mu_MiniIsoTight", &mu_MiniIsoTight, &b_mu_MiniIsoTight);
   _tree->SetBranchAddress("mu_MiniIsoVeryTight", &mu_MiniIsoVeryTight, &b_mu_MiniIsoVeryTight);
   _tree->SetBranchAddress("mu_TriggerIdLoose", &mu_TriggerIdLoose, &b_mu_TriggerIdLoose);
   _tree->SetBranchAddress("mu_InTimeMuon", &mu_InTimeMuon, &b_mu_InTimeMuon);
   _tree->SetBranchAddress("mu_MultiIsoLoose", &mu_MultiIsoLoose, &b_mu_MultiIsoLoose);
   _tree->SetBranchAddress("mu_MultiIsoMedium", &mu_MultiIsoMedium, &b_mu_MultiIsoMedium);
   _tree->SetBranchAddress("mu_numberOfMatchedStations", &mu_numberOfMatchedStations, &b_mu_numberOfMatchedStations);
   _tree->SetBranchAddress("mu_numberOfValidPixelHits", &mu_numberOfValidPixelHits, &b_mu_numberOfValidPixelHits);
   _tree->SetBranchAddress("mu_trackerLayersWithMeasurement", &mu_trackerLayersWithMeasurement, &b_mu_trackerLayersWithMeasurement);
   _tree->SetBranchAddress("mu_numberOfValidMuonHits", &mu_numberOfValidMuonHits, &b_mu_numberOfValidMuonHits);
   _tree->SetBranchAddress("mu_pixelLayersWithMeasurement", &mu_pixelLayersWithMeasurement, &b_mu_pixelLayersWithMeasurement);
   _tree->SetBranchAddress("mu_innerTrack_validFraction", &mu_innerTrack_validFraction, &b_mu_innerTrack_validFraction);
   _tree->SetBranchAddress("mu_combinedQuality_trkKink", &mu_combinedQuality_trkKink, &b_mu_combinedQuality_trkKink);
   _tree->SetBranchAddress("mu_combinedQuality_chi2LocalPosition", &mu_combinedQuality_chi2LocalPosition, &b_mu_combinedQuality_chi2LocalPosition);
   _tree->SetBranchAddress("mu_segmentCompatibility", &mu_segmentCompatibility, &b_mu_segmentCompatibility);
   _tree->SetBranchAddress("mu_dB", &mu_dB, &b_mu_dB);
   _tree->SetBranchAddress("mu_pt_default", &mu_pt_default, &b_mu_pt_default);
   _tree->SetBranchAddress("mu_isolationR03_sumPt", &mu_isolationR03_sumPt, &b_mu_isolationR03_sumPt);
   _tree->SetBranchAddress("mu_isolationR03_trackerVetoPt", &mu_isolationR03_trackerVetoPt, &b_mu_isolationR03_trackerVetoPt);
   _tree->SetBranchAddress("mu_isolationR03_emEt", &mu_isolationR03_emEt, &b_mu_isolationR03_emEt);
   _tree->SetBranchAddress("mu_isolationR03_emVetoEt", &mu_isolationR03_emVetoEt, &b_mu_isolationR03_emVetoEt);
   _tree->SetBranchAddress("mu_isolationR03_hadEt", &mu_isolationR03_hadEt, &b_mu_isolationR03_hadEt);
   _tree->SetBranchAddress("mu_isolationR03_hadVetoEt", &mu_isolationR03_hadVetoEt, &b_mu_isolationR03_hadVetoEt);
   _tree->SetBranchAddress("mu_isolationR05_sumPt", &mu_isolationR05_sumPt, &b_mu_isolationR05_sumPt);
   _tree->SetBranchAddress("mu_isolationR05_trackerVetoPt", &mu_isolationR05_trackerVetoPt, &b_mu_isolationR05_trackerVetoPt);
   _tree->SetBranchAddress("mu_isolationR05_emEt", &mu_isolationR05_emEt, &b_mu_isolationR05_emEt);
   _tree->SetBranchAddress("mu_isolationR05_emVetoEt", &mu_isolationR05_emVetoEt, &b_mu_isolationR05_emVetoEt);
   _tree->SetBranchAddress("mu_isolationR05_hadEt", &mu_isolationR05_hadEt, &b_mu_isolationR05_hadEt);
   _tree->SetBranchAddress("mu_isolationR05_hadVetoEt", &mu_isolationR05_hadVetoEt, &b_mu_isolationR05_hadVetoEt);
   _tree->SetBranchAddress("mu_pfIsolationR03_sumChargedHadronPt", &mu_pfIsolationR03_sumChargedHadronPt, &b_mu_pfIsolationR03_sumChargedHadronPt);
   _tree->SetBranchAddress("mu_pfIsolationR03_sumNeutralHadronEt", &mu_pfIsolationR03_sumNeutralHadronEt, &b_mu_pfIsolationR03_sumNeutralHadronEt);
   _tree->SetBranchAddress("mu_pfIsolationR03_sumChargedParticlePt", &mu_pfIsolationR03_sumChargedParticlePt, &b_mu_pfIsolationR03_sumChargedParticlePt);
   _tree->SetBranchAddress("mu_pfIsolationR03_sumPhotonEt", &mu_pfIsolationR03_sumPhotonEt, &b_mu_pfIsolationR03_sumPhotonEt);
   _tree->SetBranchAddress("mu_pfIsolationR03_sumNeutralHadronEtHighThreshold", &mu_pfIsolationR03_sumNeutralHadronEtHighThreshold, &b_mu_pfIsolationR03_sumNeutralHadronEtHighThreshold);
   _tree->SetBranchAddress("mu_pfIsolationR03_sumPhotonEtHighThreshold", &mu_pfIsolationR03_sumPhotonEtHighThreshold, &b_mu_pfIsolationR03_sumPhotonEtHighThreshold);
   _tree->SetBranchAddress("mu_pfIsolationR03_sumPUPt", &mu_pfIsolationR03_sumPUPt, &b_mu_pfIsolationR03_sumPUPt);
   _tree->SetBranchAddress("mu_pfIsolationR04_sumChargedHadronPt", &mu_pfIsolationR04_sumChargedHadronPt, &b_mu_pfIsolationR04_sumChargedHadronPt);
   _tree->SetBranchAddress("mu_pfIsolationR04_sumNeutralHadronEt", &mu_pfIsolationR04_sumNeutralHadronEt, &b_mu_pfIsolationR04_sumNeutralHadronEt);
   _tree->SetBranchAddress("mu_pfIsolationR04_sumChargedParticlePt", &mu_pfIsolationR04_sumChargedParticlePt, &b_mu_pfIsolationR04_sumChargedParticlePt);
   _tree->SetBranchAddress("mu_pfIsolationR04_sumPhotonEt", &mu_pfIsolationR04_sumPhotonEt, &b_mu_pfIsolationR04_sumPhotonEt);
   _tree->SetBranchAddress("mu_pfIsolationR04_sumNeutralHadronEtHighThreshold", &mu_pfIsolationR04_sumNeutralHadronEtHighThreshold, &b_mu_pfIsolationR04_sumNeutralHadronEtHighThreshold);
   _tree->SetBranchAddress("mu_pfIsolationR04_sumPhotonEtHighThreshold", &mu_pfIsolationR04_sumPhotonEtHighThreshold, &b_mu_pfIsolationR04_sumPhotonEtHighThreshold);
   _tree->SetBranchAddress("mu_pfIsolationR04_sumPUPt", &mu_pfIsolationR04_sumPUPt, &b_mu_pfIsolationR04_sumPUPt);
   _tree->SetBranchAddress("mu_pfIsoDbCorrected03", &mu_pfIsoDbCorrected03, &b_mu_pfIsoDbCorrected03);
   _tree->SetBranchAddress("mu_pfIsoDbCorrected04", &mu_pfIsoDbCorrected04, &b_mu_pfIsoDbCorrected04);
   _tree->SetBranchAddress("mu_isoTrackerBased03", &mu_isoTrackerBased03, &b_mu_isoTrackerBased03);



if(_tree->GetListOfBranches()->FindObject("mu_mc_bestDR")) {


   _tree->SetBranchAddress("mu_mc_bestDR", &mu_mc_bestDR, &b_mu_mc_bestDR);
   _tree->SetBranchAddress("mu_mc_index", &mu_mc_index, &b_mu_mc_index);
   _tree->SetBranchAddress("mu_mc_ERatio", &mu_mc_ERatio, &b_mu_mc_ERatio);

}


   _tree->SetBranchAddress("jet_n", &jet_n, &b_jet_n);
   _tree->SetBranchAddress("jet_px", &jet_px, &b_jet_px);
   _tree->SetBranchAddress("jet_py", &jet_py, &b_jet_py);
   _tree->SetBranchAddress("jet_pz", &jet_pz, &b_jet_pz);
   _tree->SetBranchAddress("jet_pt", &jet_pt, &b_jet_pt);
   _tree->SetBranchAddress("jet_eta", &jet_eta, &b_jet_eta);
   _tree->SetBranchAddress("jet_theta", &jet_theta, &b_jet_theta);
   _tree->SetBranchAddress("jet_phi", &jet_phi, &b_jet_phi);
   _tree->SetBranchAddress("jet_energy", &jet_energy, &b_jet_energy);
   _tree->SetBranchAddress("jet_mass", &jet_mass, &b_jet_mass);
   _tree->SetBranchAddress("jet_chargedEmEnergyFraction", &jet_chargedEmEnergyFraction, &b_jet_chargedEmEnergyFraction);
   _tree->SetBranchAddress("jet_neutralHadronEnergyFraction", &jet_neutralHadronEnergyFraction, &b_jet_neutralHadronEnergyFraction);
   _tree->SetBranchAddress("jet_neutralEmEnergyFraction", &jet_neutralEmEnergyFraction, &b_jet_neutralEmEnergyFraction);
   _tree->SetBranchAddress("jet_chargedHadronEnergyFraction", &jet_chargedHadronEnergyFraction, &b_jet_chargedHadronEnergyFraction);
   _tree->SetBranchAddress("jet_muonEnergyFraction", &jet_muonEnergyFraction, &b_jet_muonEnergyFraction);
   _tree->SetBranchAddress("jet_chargedMultiplicity", &jet_chargedMultiplicity, &b_jet_chargedMultiplicity);
   _tree->SetBranchAddress("jet_neutralMultiplicity", &jet_neutralMultiplicity, &b_jet_neutralMultiplicity);
   _tree->SetBranchAddress("jet_partonFlavour", &jet_partonFlavour, &b_jet_partonFlavour);
   _tree->SetBranchAddress("jet_hadronFlavour", &jet_hadronFlavour, &b_jet_hadronFlavour);
   _tree->SetBranchAddress("jet_CSVv2", &jet_CSVv2, &b_jet_CSVv2);
   _tree->SetBranchAddress("jet_CvsL", &jet_CvsL, &b_jet_CvsL);
   _tree->SetBranchAddress("jet_CvsB", &jet_CvsB, &b_jet_CvsB);
   _tree->SetBranchAddress("jet_MVA2BJets", &jet_MVA2BJets, &b_jet_MVA2BJets);
   _tree->SetBranchAddress("jet_CvsB_DeepJet_charm_tagger", &jet_CvsB_DeepJet_charm_tagger, &b_jet_CvsB_DeepJet_charm_tagger);
   _tree->SetBranchAddress("jet_CvsL_DeepJet_charm_tagger", &jet_CvsL_DeepJet_charm_tagger, &b_jet_CvsL_DeepJet_charm_tagger);
   _tree->SetBranchAddress("jet_CvsB_DeepCSV_charm_tagger", &jet_CvsB_DeepCSV_charm_tagger, &b_jet_CvsB_DeepCSV_charm_tagger);
   _tree->SetBranchAddress("jet_CvsL_DeepCSV_charm_tagger", &jet_CvsL_DeepCSV_charm_tagger, &b_jet_CvsL_DeepCSV_charm_tagger);
   _tree->SetBranchAddress("jet_DeepJet", &jet_DeepJet, &b_jet_DeepJet);
   _tree->SetBranchAddress("jet_DeepCSV", &jet_DeepCSV, &b_jet_DeepCSV);
   _tree->SetBranchAddress("jet_isJetIDLoose_2016", &jet_isJetIDLoose_2016, &b_jet_isJetIDLoose_2016);
   _tree->SetBranchAddress("jet_isJetIDTight_2016", &jet_isJetIDTight_2016, &b_jet_isJetIDTight_2016);
   _tree->SetBranchAddress("jet_isJetIDTightLepVeto_2016", &jet_isJetIDTightLepVeto_2016, &b_jet_isJetIDTightLepVeto_2016);
   _tree->SetBranchAddress("jet_isJetID_2017", &jet_isJetID_2017, &b_jet_isJetID_2017);
   _tree->SetBranchAddress("jet_isJetIDLepVeto_2017", &jet_isJetIDLepVeto_2017, &b_jet_isJetIDLepVeto_2017);
   _tree->SetBranchAddress("jet_isJetID_2018", &jet_isJetID_2018, &b_jet_isJetID_2018);
   _tree->SetBranchAddress("jet_isJetIDLepVeto_2018", &jet_isJetIDLepVeto_2018, &b_jet_isJetIDLepVeto_2018);


if(_tree->GetListOfBranches()->FindObject("jet_Smeared_pt")) {


   _tree->SetBranchAddress("jet_Smeared_pt", &jet_Smeared_pt, &b_jet_Smeared_pt);
   _tree->SetBranchAddress("jet_SmearedJetResUp_pt", &jet_SmearedJetResUp_pt, &b_jet_SmearedJetResUp_pt);
   _tree->SetBranchAddress("jet_SmearedJetResDown_pt", &jet_SmearedJetResDown_pt, &b_jet_SmearedJetResDown_pt);
   _tree->SetBranchAddress("jet_EnUp_pt", &jet_EnUp_pt, &b_jet_EnUp_pt);
   _tree->SetBranchAddress("jet_EnDown_pt", &jet_EnDown_pt, &b_jet_EnDown_pt);
   _tree->SetBranchAddress("BtagSF_Deepcsv_loose", &BtagSF_Deepcsv_loose, &b_BtagSF_Deepcsv_loose);
   _tree->SetBranchAddress("BtagSF_Deepcsv_Up_loose", &BtagSF_Deepcsv_Up_loose, &b_BtagSF_Deepcsv_Up_loose);
   _tree->SetBranchAddress("BtagSF_Deepcsv_Down_loose", &BtagSF_Deepcsv_Down_loose, &b_BtagSF_Deepcsv_Down_loose);
   _tree->SetBranchAddress("BtagSF_Deepcsv_medium", &BtagSF_Deepcsv_medium, &b_BtagSF_Deepcsv_medium);
   _tree->SetBranchAddress("BtagSF_Deepcsv_Up_medium", &BtagSF_Deepcsv_Up_medium, &b_BtagSF_Deepcsv_Up_medium);
   _tree->SetBranchAddress("BtagSF_Deepcsv_Down_medium", &BtagSF_Deepcsv_Down_medium, &b_BtagSF_Deepcsv_Down_medium);
   _tree->SetBranchAddress("BtagSF_Deepcsv_tight", &BtagSF_Deepcsv_tight, &b_BtagSF_Deepcsv_tight);
   _tree->SetBranchAddress("BtagSF_Deepcsv_Up_tight", &BtagSF_Deepcsv_Up_tight, &b_BtagSF_Deepcsv_Up_tight);
   _tree->SetBranchAddress("BtagSF_Deepcsv_Down_tight", &BtagSF_Deepcsv_Down_tight, &b_BtagSF_Deepcsv_Down_tight);

 }




   _tree->SetBranchAddress("MET_nominal_Pt", &MET_nominal_Pt, &b_MET_nominal_Pt);
   _tree->SetBranchAddress("MET_nominal_Px", &MET_nominal_Px, &b_MET_nominal_Px);
   _tree->SetBranchAddress("MET_nominal_Py", &MET_nominal_Py, &b_MET_nominal_Py);
   _tree->SetBranchAddress("MET_nominal_phi", &MET_nominal_phi, &b_MET_nominal_phi);
   _tree->SetBranchAddress("MET_nominal_significance", &MET_nominal_significance, &b_MET_nominal_significance);
   _tree->SetBranchAddress("MET_Pt", &MET_Pt, &b_MET_Pt);
   _tree->SetBranchAddress("MET_Px", &MET_Px, &b_MET_Px);
   _tree->SetBranchAddress("MET_Py", &MET_Py, &b_MET_Py);
   _tree->SetBranchAddress("MET_phi", &MET_phi, &b_MET_phi);
   _tree->SetBranchAddress("MET_significance", &MET_significance, &b_MET_significance);
   _tree->SetBranchAddress("MET_T1_Pt", &MET_T1_Pt, &b_MET_T1_Pt);
   _tree->SetBranchAddress("MET_T1_Px", &MET_T1_Px, &b_MET_T1_Px);
   _tree->SetBranchAddress("MET_T1_Py", &MET_T1_Py, &b_MET_T1_Py);
   _tree->SetBranchAddress("MET_T1_phi", &MET_T1_phi, &b_MET_T1_phi);
   _tree->SetBranchAddress("MET_T1_significance", &MET_T1_significance, &b_MET_T1_significance);


   if(_tree->GetListOfBranches()->FindObject("MET_gen_pt")){
   _tree->SetBranchAddress("MET_gen_pt", &MET_gen_pt, &b_MET_gen_pt);
   _tree->SetBranchAddress("MET_gen_phi", &MET_gen_phi, &b_MET_gen_phi);
   }


   if(_tree->GetListOfBranches()->FindObject("MET_Type1Unc_Px")){
   _tree->SetBranchAddress("MET_Type1Unc_Px", &MET_Type1Unc_Px, &b_MET_Type1Unc_Px);
   _tree->SetBranchAddress("MET_Type1Unc_Py", &MET_Type1Unc_Py, &b_MET_Type1Unc_Py);
   _tree->SetBranchAddress("MET_Type1Unc_Pt", &MET_Type1Unc_Pt, &b_MET_Type1Unc_Pt);
   _tree->SetBranchAddress("MET_T1JetEnDown_Pt", &MET_T1JetEnDown_Pt, &b_MET_T1JetEnDown_Pt);
   _tree->SetBranchAddress("MET_T1JetEnDown_Px", &MET_T1JetEnDown_Px, &b_MET_T1JetEnDown_Px);
   _tree->SetBranchAddress("MET_T1JetEnDown_Py", &MET_T1JetEnDown_Py, &b_MET_T1JetEnDown_Py);
   _tree->SetBranchAddress("MET_T1JetEnDown_phi", &MET_T1JetEnDown_phi, &b_MET_T1JetEnDown_phi);
   _tree->SetBranchAddress("MET_T1JetEnDown_significance", &MET_T1JetEnDown_significance, &b_MET_T1JetEnDown_significance);
   _tree->SetBranchAddress("MET_T1JetEnUp_Pt", &MET_T1JetEnUp_Pt, &b_MET_T1JetEnUp_Pt);
   _tree->SetBranchAddress("MET_T1JetEnUp_Px", &MET_T1JetEnUp_Px, &b_MET_T1JetEnUp_Px);
   _tree->SetBranchAddress("MET_T1JetEnUp_Py", &MET_T1JetEnUp_Py, &b_MET_T1JetEnUp_Py);
   _tree->SetBranchAddress("MET_T1JetEnUp_phi", &MET_T1JetEnUp_phi, &b_MET_T1JetEnUp_phi);
   _tree->SetBranchAddress("MET_T1JetEnUp_significance", &MET_T1JetEnUp_significance, &b_MET_T1JetEnUp_significance);
   _tree->SetBranchAddress("MET_T1Smear_Pt", &MET_T1Smear_Pt, &b_MET_T1Smear_Pt);
   _tree->SetBranchAddress("MET_T1Smear_Px", &MET_T1Smear_Px, &b_MET_T1Smear_Px);
   _tree->SetBranchAddress("MET_T1Smear_Py", &MET_T1Smear_Py, &b_MET_T1Smear_Py);
   _tree->SetBranchAddress("MET_T1Smear_phi", &MET_T1Smear_phi, &b_MET_T1Smear_phi);
   _tree->SetBranchAddress("MET_T1Smear_significance", &MET_T1Smear_significance, &b_MET_T1Smear_significance);
   _tree->SetBranchAddress("MET_T1SmearJetEnDown_Pt", &MET_T1SmearJetEnDown_Pt, &b_MET_T1SmearJetEnDown_Pt);
   _tree->SetBranchAddress("MET_T1SmearJetEnDown_Px", &MET_T1SmearJetEnDown_Px, &b_MET_T1SmearJetEnDown_Px);
   _tree->SetBranchAddress("MET_T1SmearJetEnDown_Py", &MET_T1SmearJetEnDown_Py, &b_MET_T1SmearJetEnDown_Py);
   _tree->SetBranchAddress("MET_T1SmearJetEnDown_phi", &MET_T1SmearJetEnDown_phi, &b_MET_T1SmearJetEnDown_phi);
   _tree->SetBranchAddress("MET_T1SmearJetEnDown_significance", &MET_T1SmearJetEnDown_significance, &b_MET_T1SmearJetEnDown_significance);
   _tree->SetBranchAddress("MET_T1SmearJetEnUp_Pt", &MET_T1SmearJetEnUp_Pt, &b_MET_T1SmearJetEnUp_Pt);
   _tree->SetBranchAddress("MET_T1SmearJetEnUp_Px", &MET_T1SmearJetEnUp_Px, &b_MET_T1SmearJetEnUp_Px);
   _tree->SetBranchAddress("MET_T1SmearJetEnUp_Py", &MET_T1SmearJetEnUp_Py, &b_MET_T1SmearJetEnUp_Py);
   _tree->SetBranchAddress("MET_T1SmearJetEnUp_phi", &MET_T1SmearJetEnUp_phi, &b_MET_T1SmearJetEnUp_phi);
   _tree->SetBranchAddress("MET_T1SmearJetEnUp_significance", &MET_T1SmearJetEnUp_significance, &b_MET_T1SmearJetEnUp_significance);
   _tree->SetBranchAddress("MET_T1SmearJetResDown_Pt", &MET_T1SmearJetResDown_Pt, &b_MET_T1SmearJetResDown_Pt);
   _tree->SetBranchAddress("MET_T1SmearJetResDown_Px", &MET_T1SmearJetResDown_Px, &b_MET_T1SmearJetResDown_Px);
   _tree->SetBranchAddress("MET_T1SmearJetResDown_Py", &MET_T1SmearJetResDown_Py, &b_MET_T1SmearJetResDown_Py);
   _tree->SetBranchAddress("MET_T1SmearJetResDown_phi", &MET_T1SmearJetResDown_phi, &b_MET_T1SmearJetResDown_phi);
   _tree->SetBranchAddress("MET_T1SmearJetResDown_significance", &MET_T1SmearJetResDown_significance, &b_MET_T1SmearJetResDown_significance);
   _tree->SetBranchAddress("MET_T1SmearJetResUp_Pt", &MET_T1SmearJetResUp_Pt, &b_MET_T1SmearJetResUp_Pt);
   _tree->SetBranchAddress("MET_T1SmearJetResUp_Px", &MET_T1SmearJetResUp_Px, &b_MET_T1SmearJetResUp_Px);
   _tree->SetBranchAddress("MET_T1SmearJetResUp_Py", &MET_T1SmearJetResUp_Py, &b_MET_T1SmearJetResUp_Py);
   _tree->SetBranchAddress("MET_T1SmearJetResUp_phi", &MET_T1SmearJetResUp_phi, &b_MET_T1SmearJetResUp_phi);
   _tree->SetBranchAddress("MET_T1SmearJetResUp_significance", &MET_T1SmearJetResUp_significance, &b_MET_T1SmearJetResUp_significance);

   }


   _tree->SetBranchAddress("MET_T1Txy_Pt", &MET_T1Txy_Pt, &b_MET_T1Txy_Pt);
   _tree->SetBranchAddress("MET_T1Txy_Px", &MET_T1Txy_Px, &b_MET_T1Txy_Px);
   _tree->SetBranchAddress("MET_T1Txy_Py", &MET_T1Txy_Py, &b_MET_T1Txy_Py);
   _tree->SetBranchAddress("MET_T1Txy_phi", &MET_T1Txy_phi, &b_MET_T1Txy_phi);
   _tree->SetBranchAddress("MET_T1Txy_significance", &MET_T1Txy_significance, &b_MET_T1Txy_significance);
   _tree->SetBranchAddress("MET_FinalCollection_Pt", &MET_FinalCollection_Pt, &b_MET_FinalCollection_Pt);
   _tree->SetBranchAddress("MET_FinalCollection_Px", &MET_FinalCollection_Px, &b_MET_FinalCollection_Px);
   _tree->SetBranchAddress("MET_FinalCollection_Py", &MET_FinalCollection_Py, &b_MET_FinalCollection_Py);
   _tree->SetBranchAddress("MET_FinalCollection_phi", &MET_FinalCollection_phi, &b_MET_FinalCollection_phi);
   _tree->SetBranchAddress("MET_FinalCollection_significance", &MET_FinalCollection_significance, &b_MET_FinalCollection_significance);

   _tree->SetBranchAddress("tau_n", &tau_n, &b_tau_n);
   _tree->SetBranchAddress("tau_px", &tau_px, &b_tau_px);
   _tree->SetBranchAddress("tau_py", &tau_py, &b_tau_py);
   _tree->SetBranchAddress("tau_pz", &tau_pz, &b_tau_pz);
   _tree->SetBranchAddress("tau_pt", &tau_pt, &b_tau_pt);
   _tree->SetBranchAddress("tau_eta", &tau_eta, &b_tau_eta);
   _tree->SetBranchAddress("tau_theta", &tau_theta, &b_tau_theta);
   _tree->SetBranchAddress("tau_phi", &tau_phi, &b_tau_phi);
   _tree->SetBranchAddress("tau_energy", &tau_energy, &b_tau_energy);
   _tree->SetBranchAddress("tau_mass", &tau_mass, &b_tau_mass);
   _tree->SetBranchAddress("tau_dxy", &tau_dxy, &b_tau_dxy);
   _tree->SetBranchAddress("tau_dxy_error", &tau_dxy_error, &b_tau_dxy_error);
   _tree->SetBranchAddress("tau_ptLeadChargedCand", &tau_ptLeadChargedCand, &b_tau_ptLeadChargedCand);
   _tree->SetBranchAddress("tau_byIsolationMVArun2v1DBoldDMwLTraw", &tau_byIsolationMVArun2v1DBoldDMwLTraw, &b_tau_byIsolationMVArun2v1DBoldDMwLTraw);
   _tree->SetBranchAddress("tau_byCombinedIsolationDeltaBetaCorrRaw3Hits", &tau_byCombinedIsolationDeltaBetaCorrRaw3Hits, &b_tau_byCombinedIsolationDeltaBetaCorrRaw3Hits);
   _tree->SetBranchAddress("tau_byIsolationMVArun2v1DBnewDMwLTraw", &tau_byIsolationMVArun2v1DBnewDMwLTraw, &b_tau_byIsolationMVArun2v1DBnewDMwLTraw);
   _tree->SetBranchAddress("tau_byIsolationMVArun2v1PWoldDMwLTraw", &tau_byIsolationMVArun2v1PWoldDMwLTraw, &b_tau_byIsolationMVArun2v1PWoldDMwLTraw);
   _tree->SetBranchAddress("tau_byIsolationMVArun2v1PWnewDMwLTraw", &tau_byIsolationMVArun2v1PWnewDMwLTraw, &b_tau_byIsolationMVArun2v1PWnewDMwLTraw);
   _tree->SetBranchAddress("tau_byIsolationMVArun2v1DBdR03oldDMwLTraw", &tau_byIsolationMVArun2v1DBdR03oldDMwLTraw, &b_tau_byIsolationMVArun2v1DBdR03oldDMwLTraw);
   _tree->SetBranchAddress("tau_byIsolationMVArun2v1PWdR03oldDMwLTraw", &tau_byIsolationMVArun2v1PWdR03oldDMwLTraw, &b_tau_byIsolationMVArun2v1PWdR03oldDMwLTraw);
   _tree->SetBranchAddress("tau_againstElectronMVA6Raw", &tau_againstElectronMVA6Raw, &b_tau_againstElectronMVA6Raw);
   _tree->SetBranchAddress("tau_againstElectronMVA6category", &tau_againstElectronMVA6category, &b_tau_againstElectronMVA6category);

   if(_tree->GetListOfBranches()->FindObject("tau_mc_bestDR")){
   _tree->SetBranchAddress("tau_mc_bestDR", &tau_mc_bestDR, &b_tau_mc_bestDR);
   _tree->SetBranchAddress("tau_mc_ERatio", &tau_mc_ERatio, &b_tau_mc_ERatio);
   _tree->SetBranchAddress("tau_mc_index", &tau_mc_index, &b_tau_mc_index);

   }

   _tree->SetBranchAddress("tau_decayModeFinding", &tau_decayModeFinding, &b_tau_decayModeFinding);
   _tree->SetBranchAddress("tau_decayModeFindingNewDMs", &tau_decayModeFindingNewDMs, &b_tau_decayModeFindingNewDMs);
   _tree->SetBranchAddress("tau_againstMuonLoose3", &tau_againstMuonLoose3, &b_tau_againstMuonLoose3);
   _tree->SetBranchAddress("tau_againstMuonTight3", &tau_againstMuonTight3, &b_tau_againstMuonTight3);
   _tree->SetBranchAddress("tau_byLooseCombinedIsolationDeltaBetaCorr3Hits", &tau_byLooseCombinedIsolationDeltaBetaCorr3Hits, &b_tau_byLooseCombinedIsolationDeltaBetaCorr3Hits);
   _tree->SetBranchAddress("tau_byMediumCombinedIsolationDeltaBetaCorr3Hits", &tau_byMediumCombinedIsolationDeltaBetaCorr3Hits, &b_tau_byMediumCombinedIsolationDeltaBetaCorr3Hits);
   _tree->SetBranchAddress("tau_byTightCombinedIsolationDeltaBetaCorr3Hits", &tau_byTightCombinedIsolationDeltaBetaCorr3Hits, &b_tau_byTightCombinedIsolationDeltaBetaCorr3Hits);
   _tree->SetBranchAddress("tau_byVLooseIsolationMVArun2v1DBoldDMwLT", &tau_byVLooseIsolationMVArun2v1DBoldDMwLT, &b_tau_byVLooseIsolationMVArun2v1DBoldDMwLT);
   _tree->SetBranchAddress("tau_byLooseIsolationMVArun2v1DBoldDMwLT", &tau_byLooseIsolationMVArun2v1DBoldDMwLT, &b_tau_byLooseIsolationMVArun2v1DBoldDMwLT);
   _tree->SetBranchAddress("tau_byMediumIsolationMVArun2v1DBoldDMwLT", &tau_byMediumIsolationMVArun2v1DBoldDMwLT, &b_tau_byMediumIsolationMVArun2v1DBoldDMwLT);
   _tree->SetBranchAddress("tau_byTightIsolationMVArun2v1DBoldDMwLT", &tau_byTightIsolationMVArun2v1DBoldDMwLT, &b_tau_byTightIsolationMVArun2v1DBoldDMwLT);
   _tree->SetBranchAddress("tau_byVTightIsolationMVArun2v1DBoldDMwLT", &tau_byVTightIsolationMVArun2v1DBoldDMwLT, &b_tau_byVTightIsolationMVArun2v1DBoldDMwLT);
   _tree->SetBranchAddress("tau_byVVTightIsolationMVArun2v1DBoldDMwLT", &tau_byVVTightIsolationMVArun2v1DBoldDMwLT, &b_tau_byVVTightIsolationMVArun2v1DBoldDMwLT);
   _tree->SetBranchAddress("tau_byVLooseIsolationMVArun2v1DBnewDMwLT", &tau_byVLooseIsolationMVArun2v1DBnewDMwLT, &b_tau_byVLooseIsolationMVArun2v1DBnewDMwLT);
   _tree->SetBranchAddress("tau_byLooseIsolationMVArun2v1DBnewDMwLT", &tau_byLooseIsolationMVArun2v1DBnewDMwLT, &b_tau_byLooseIsolationMVArun2v1DBnewDMwLT);
   _tree->SetBranchAddress("tau_byMediumIsolationMVArun2v1DBnewDMwLT", &tau_byMediumIsolationMVArun2v1DBnewDMwLT, &b_tau_byMediumIsolationMVArun2v1DBnewDMwLT);
   _tree->SetBranchAddress("tau_byTightIsolationMVArun2v1DBnewDMwLT", &tau_byTightIsolationMVArun2v1DBnewDMwLT, &b_tau_byTightIsolationMVArun2v1DBnewDMwLT);
   _tree->SetBranchAddress("tau_byVTightIsolationMVArun2v1DBnewDMwLT", &tau_byVTightIsolationMVArun2v1DBnewDMwLT, &b_tau_byVTightIsolationMVArun2v1DBnewDMwLT);
   _tree->SetBranchAddress("tau_byVVTightIsolationMVArun2v1DBnewDMwLT", &tau_byVVTightIsolationMVArun2v1DBnewDMwLT, &b_tau_byVVTightIsolationMVArun2v1DBnewDMwLT);
   _tree->SetBranchAddress("tau_byVLooseIsolationMVArun2v1PWoldDMwLT", &tau_byVLooseIsolationMVArun2v1PWoldDMwLT, &b_tau_byVLooseIsolationMVArun2v1PWoldDMwLT);
   _tree->SetBranchAddress("tau_byLooseIsolationMVArun2v1PWoldDMwLT", &tau_byLooseIsolationMVArun2v1PWoldDMwLT, &b_tau_byLooseIsolationMVArun2v1PWoldDMwLT);
   _tree->SetBranchAddress("tau_byMediumIsolationMVArun2v1PWoldDMwLT", &tau_byMediumIsolationMVArun2v1PWoldDMwLT, &b_tau_byMediumIsolationMVArun2v1PWoldDMwLT);
   _tree->SetBranchAddress("tau_byTightIsolationMVArun2v1PWoldDMwLT", &tau_byTightIsolationMVArun2v1PWoldDMwLT, &b_tau_byTightIsolationMVArun2v1PWoldDMwLT);
   _tree->SetBranchAddress("tau_byVTightIsolationMVArun2v1PWoldDMwLT", &tau_byVTightIsolationMVArun2v1PWoldDMwLT, &b_tau_byVTightIsolationMVArun2v1PWoldDMwLT);
   _tree->SetBranchAddress("tau_byVVTightIsolationMVArun2v1PWoldDMwLT", &tau_byVVTightIsolationMVArun2v1PWoldDMwLT, &b_tau_byVVTightIsolationMVArun2v1PWoldDMwLT);
   _tree->SetBranchAddress("tau_byVLooseIsolationMVArun2v1PWnewDMwLT", &tau_byVLooseIsolationMVArun2v1PWnewDMwLT, &b_tau_byVLooseIsolationMVArun2v1PWnewDMwLT);
   _tree->SetBranchAddress("tau_byLooseIsolationMVArun2v1PWnewDMwLT", &tau_byLooseIsolationMVArun2v1PWnewDMwLT, &b_tau_byLooseIsolationMVArun2v1PWnewDMwLT);
   _tree->SetBranchAddress("tau_byMediumIsolationMVArun2v1PWnewDMwLT", &tau_byMediumIsolationMVArun2v1PWnewDMwLT, &b_tau_byMediumIsolationMVArun2v1PWnewDMwLT);
   _tree->SetBranchAddress("tau_byTightIsolationMVArun2v1PWnewDMwLT", &tau_byTightIsolationMVArun2v1PWnewDMwLT, &b_tau_byTightIsolationMVArun2v1PWnewDMwLT);
   _tree->SetBranchAddress("tau_byVTightIsolationMVArun2v1PWnewDMwLT", &tau_byVTightIsolationMVArun2v1PWnewDMwLT, &b_tau_byVTightIsolationMVArun2v1PWnewDMwLT);
   _tree->SetBranchAddress("tau_byVVTightIsolationMVArun2v1PWnewDMwLT", &tau_byVVTightIsolationMVArun2v1PWnewDMwLT, &b_tau_byVVTightIsolationMVArun2v1PWnewDMwLT);
   _tree->SetBranchAddress("tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT", &tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT, &b_tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT);
   _tree->SetBranchAddress("tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT", &tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT, &b_tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT);
   _tree->SetBranchAddress("tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT", &tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT, &b_tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT);
   _tree->SetBranchAddress("tau_byTightIsolationMVArun2v1DBdR03oldDMwLT", &tau_byTightIsolationMVArun2v1DBdR03oldDMwLT, &b_tau_byTightIsolationMVArun2v1DBdR03oldDMwLT);
   _tree->SetBranchAddress("tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT", &tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT, &b_tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT);
   _tree->SetBranchAddress("tau_byVVTightIsolationMVArun2v1DBdR03oldDMwLT", &tau_byVVTightIsolationMVArun2v1DBdR03oldDMwLT, &b_tau_byVVTightIsolationMVArun2v1DBdR03oldDMwLT);
   _tree->SetBranchAddress("tau_byVLooseIsolationMVArun2v1PWdR03oldDMwLT", &tau_byVLooseIsolationMVArun2v1PWdR03oldDMwLT, &b_tau_byVLooseIsolationMVArun2v1PWdR03oldDMwLT);
   _tree->SetBranchAddress("tau_byLooseIsolationMVArun2v1PWdR03oldDMwLT", &tau_byLooseIsolationMVArun2v1PWdR03oldDMwLT, &b_tau_byLooseIsolationMVArun2v1PWdR03oldDMwLT);
   _tree->SetBranchAddress("tau_byMediumIsolationMVArun2v1PWdR03oldDMwLT", &tau_byMediumIsolationMVArun2v1PWdR03oldDMwLT, &b_tau_byMediumIsolationMVArun2v1PWdR03oldDMwLT);
   _tree->SetBranchAddress("tau_byTightIsolationMVArun2v1PWdR03oldDMwLT", &tau_byTightIsolationMVArun2v1PWdR03oldDMwLT, &b_tau_byTightIsolationMVArun2v1PWdR03oldDMwLT);
   _tree->SetBranchAddress("tau_byVTightIsolationMVArun2v1PWdR03oldDMwLT", &tau_byVTightIsolationMVArun2v1PWdR03oldDMwLT, &b_tau_byVTightIsolationMVArun2v1PWdR03oldDMwLT);
   _tree->SetBranchAddress("tau_byVVTightIsolationMVArun2v1PWdR03oldDMwLT", &tau_byVVTightIsolationMVArun2v1PWdR03oldDMwLT, &b_tau_byVVTightIsolationMVArun2v1PWdR03oldDMwLT);
   _tree->SetBranchAddress("tau_againstElectronVLooseMVA6", &tau_againstElectronVLooseMVA6, &b_tau_againstElectronVLooseMVA6);
   _tree->SetBranchAddress("tau_againstElectronLooseMVA6", &tau_againstElectronLooseMVA6, &b_tau_againstElectronLooseMVA6);
   _tree->SetBranchAddress("tau_againstElectronMediumMVA6", &tau_againstElectronMediumMVA6, &b_tau_againstElectronMediumMVA6);
   _tree->SetBranchAddress("tau_againstElectronTightMVA6", &tau_againstElectronTightMVA6, &b_tau_againstElectronTightMVA6);
   _tree->SetBranchAddress("tau_againstElectronVTightMVA6", &tau_againstElectronVTightMVA6, &b_tau_againstElectronVTightMVA6);
   _tree->SetBranchAddress("tau_chargedIsoPtSum", &tau_chargedIsoPtSum, &b_tau_chargedIsoPtSum);
   _tree->SetBranchAddress("tau_neutralIsoPtSum", &tau_neutralIsoPtSum, &b_tau_neutralIsoPtSum);
   _tree->SetBranchAddress("tau_puCorrPtSum", &tau_puCorrPtSum, &b_tau_puCorrPtSum);
   _tree->SetBranchAddress("tau_footprintCorrection", &tau_footprintCorrection, &b_tau_footprintCorrection);
   _tree->SetBranchAddress("tau_neutralIsoPtSumWeight", &tau_neutralIsoPtSumWeight, &b_tau_neutralIsoPtSumWeight);
   _tree->SetBranchAddress("tau_photonPtSumOutsideSignalCone", &tau_photonPtSumOutsideSignalCone, &b_tau_photonPtSumOutsideSignalCone);
   _tree->SetBranchAddress("tau_byPhotonPtSumOutsideSignalCone", &tau_byPhotonPtSumOutsideSignalCone, &b_tau_byPhotonPtSumOutsideSignalCone);
   _tree->SetBranchAddress("tau_footprintCorrectiondR03", &tau_footprintCorrectiondR03, &b_tau_footprintCorrectiondR03);
   _tree->SetBranchAddress("tau_chargedIsoPtSumdR03", &tau_chargedIsoPtSumdR03, &b_tau_chargedIsoPtSumdR03);
   _tree->SetBranchAddress("tau_neutralIsoPtSumWeightdR03", &tau_neutralIsoPtSumWeightdR03, &b_tau_neutralIsoPtSumWeightdR03);
   _tree->SetBranchAddress("tau_neutralIsoPtSumdR03", &tau_neutralIsoPtSumdR03, &b_tau_neutralIsoPtSumdR03);
   _tree->SetBranchAddress("tau_photonPtSumOutsideSignalConedR03", &tau_photonPtSumOutsideSignalConedR03, &b_tau_photonPtSumOutsideSignalConedR03);
   _tree->SetBranchAddress("tau_PFChargedHadIso", &tau_PFChargedHadIso, &b_tau_PFChargedHadIso);
   _tree->SetBranchAddress("tau_PFNeutralHadIso", &tau_PFNeutralHadIso, &b_tau_PFNeutralHadIso);
   _tree->SetBranchAddress("tau_PFPhotonIso", &tau_PFPhotonIso, &b_tau_PFPhotonIso);
   _tree->SetBranchAddress("tau_leadChargedParticlePt", &tau_leadChargedParticlePt, &b_tau_leadChargedParticlePt);
   _tree->SetBranchAddress("tau_trackRefPt", &tau_trackRefPt, &b_tau_trackRefPt);
   _tree->SetBranchAddress("tau_lead_dxy", &tau_lead_dxy, &b_tau_lead_dxy);
   _tree->SetBranchAddress("tau_lead_dz", &tau_lead_dz, &b_tau_lead_dz);
   _tree->SetBranchAddress("tau_dxy_Sig", &tau_dxy_Sig, &b_tau_dxy_Sig);
   _tree->SetBranchAddress("tau_flightLengthSig", &tau_flightLengthSig, &b_tau_flightLengthSig);
   _tree->SetBranchAddress("tau_ip3d", &tau_ip3d, &b_tau_ip3d);
   _tree->SetBranchAddress("tau_ip3d_Sig", &tau_ip3d_Sig, &b_tau_ip3d_Sig);
   _tree->SetBranchAddress("tau_decayDistX", &tau_decayDistX, &b_tau_decayDistX);
   _tree->SetBranchAddress("tau_decayDistY", &tau_decayDistY, &b_tau_decayDistY);
   _tree->SetBranchAddress("tau_decayDistZ", &tau_decayDistZ, &b_tau_decayDistZ);
   _tree->SetBranchAddress("tau_decayDistMag", &tau_decayDistMag, &b_tau_decayDistMag);
   _tree->SetBranchAddress("tau_nPhoton", &tau_nPhoton, &b_tau_nPhoton);
   _tree->SetBranchAddress("tau_ptWeightedDetaStrip", &tau_ptWeightedDetaStrip, &b_tau_ptWeightedDetaStrip);
   _tree->SetBranchAddress("tau_ptWeightedDphiStrip", &tau_ptWeightedDphiStrip, &b_tau_ptWeightedDphiStrip);
   _tree->SetBranchAddress("tau_ptWeightedDrSignal", &tau_ptWeightedDrSignal, &b_tau_ptWeightedDrSignal);
   _tree->SetBranchAddress("tau_ptWeightedDrIsolation", &tau_ptWeightedDrIsolation, &b_tau_ptWeightedDrIsolation);
   _tree->SetBranchAddress("tau_leadingTrackChi2", &tau_leadingTrackChi2, &b_tau_leadingTrackChi2);
   _tree->SetBranchAddress("tau_eRatio", &tau_eRatio, &b_tau_eRatio);
   _tree->SetBranchAddress("tau_gjAngleDiff", &tau_gjAngleDiff, &b_tau_gjAngleDiff);
   _tree->SetBranchAddress("tau_numberOfIsolationChargedHadrCands", &tau_numberOfIsolationChargedHadrCands, &b_tau_numberOfIsolationChargedHadrCands);
   _tree->SetBranchAddress("tau_numberOfSignalChargedHadrCands", &tau_numberOfSignalChargedHadrCands, &b_tau_numberOfSignalChargedHadrCands);
   _tree->SetBranchAddress("tau_numNeutralHadronsSignalCone", &tau_numNeutralHadronsSignalCone, &b_tau_numNeutralHadronsSignalCone);
   _tree->SetBranchAddress("tau_numPhotonsSignalCone", &tau_numPhotonsSignalCone, &b_tau_numPhotonsSignalCone);
   _tree->SetBranchAddress("tau_numParticlesSignalCone", &tau_numParticlesSignalCone, &b_tau_numParticlesSignalCone);
   _tree->SetBranchAddress("tau_numChargedParticlesIsoCone", &tau_numChargedParticlesIsoCone, &b_tau_numChargedParticlesIsoCone);
   _tree->SetBranchAddress("tau_numNeutralHadronsIsoCone", &tau_numNeutralHadronsIsoCone, &b_tau_numNeutralHadronsIsoCone);
   _tree->SetBranchAddress("tau_numPhotonsIsoCone", &tau_numPhotonsIsoCone, &b_tau_numPhotonsIsoCone);
   _tree->SetBranchAddress("tau_numParticlesIsoCone", &tau_numParticlesIsoCone, &b_tau_numParticlesIsoCone);

   _tree->SetBranchAddress("tau_decayMode", &tau_decayMode, &b_tau_decayMode);
   _tree->SetBranchAddress("tau_charge", &tau_charge, &b_tau_charge);
   _tree->SetBranchAddress("tau_isPFTau", &tau_isPFTau, &b_tau_isPFTau);
   _tree->SetBranchAddress("tau_hasSecondaryVertex", &tau_hasSecondaryVertex, &b_tau_hasSecondaryVertex);
   _tree->SetBranchAddress("tau_leadChargedHadrAvailable", &tau_leadChargedHadrAvailable, &b_tau_leadChargedHadrAvailable);
   _tree->SetBranchAddress("tau_byIsolationMVArun2017v1DBoldDMwLTraw2017", &tau_byIsolationMVArun2017v1DBoldDMwLTraw2017, &b_tau_byIsolationMVArun2017v1DBoldDMwLTraw2017);
   _tree->SetBranchAddress("tau_byIsolationMVArun2017v2DBnewDMwLTraw2017", &tau_byIsolationMVArun2017v2DBnewDMwLTraw2017, &b_tau_byIsolationMVArun2017v2DBnewDMwLTraw2017);
   _tree->SetBranchAddress("tau_byIsolationMVArun2017v2DBoldDMdR0p3wLTraw2017", &tau_byIsolationMVArun2017v2DBoldDMdR0p3wLTraw2017, &b_tau_byIsolationMVArun2017v2DBoldDMdR0p3wLTraw2017);
   _tree->SetBranchAddress("tau_byIsolationMVArun2017v2DBoldDMwLTraw2017", &tau_byIsolationMVArun2017v2DBoldDMwLTraw2017, &b_tau_byIsolationMVArun2017v2DBoldDMwLTraw2017);
   _tree->SetBranchAddress("tau_byIsolationMVArun2v1DBnewDMwLTraw2016", &tau_byIsolationMVArun2v1DBnewDMwLTraw2016, &b_tau_byIsolationMVArun2v1DBnewDMwLTraw2016);
   _tree->SetBranchAddress("tau_byIsolationMVArun2v1DBoldDMwLTraw2016", &tau_byIsolationMVArun2v1DBoldDMwLTraw2016, &b_tau_byIsolationMVArun2v1DBoldDMwLTraw2016);
   _tree->SetBranchAddress("tau_byVVLooseIsolationMVArun2017v1DBoldDMwLT2017", &tau_byVVLooseIsolationMVArun2017v1DBoldDMwLT2017, &b_tau_byVVLooseIsolationMVArun2017v1DBoldDMwLT2017);
   _tree->SetBranchAddress("tau_byVLooseIsolationMVArun2017v1DBoldDMwLT2017", &tau_byVLooseIsolationMVArun2017v1DBoldDMwLT2017, &b_tau_byVLooseIsolationMVArun2017v1DBoldDMwLT2017);
   _tree->SetBranchAddress("tau_byLooseIsolationMVArun2017v1DBoldDMwLT2017", &tau_byLooseIsolationMVArun2017v1DBoldDMwLT2017, &b_tau_byLooseIsolationMVArun2017v1DBoldDMwLT2017);
   _tree->SetBranchAddress("tau_byMediumIsolationMVArun2017v1DBoldDMwLT2017", &tau_byMediumIsolationMVArun2017v1DBoldDMwLT2017, &b_tau_byMediumIsolationMVArun2017v1DBoldDMwLT2017);
   _tree->SetBranchAddress("tau_byTightIsolationMVArun2017v1DBoldDMwLT2017", &tau_byTightIsolationMVArun2017v1DBoldDMwLT2017, &b_tau_byTightIsolationMVArun2017v1DBoldDMwLT2017);
   _tree->SetBranchAddress("tau_byVTightIsolationMVArun2017v1DBoldDMwLT2017", &tau_byVTightIsolationMVArun2017v1DBoldDMwLT2017, &b_tau_byVTightIsolationMVArun2017v1DBoldDMwLT2017);
   _tree->SetBranchAddress("tau_byVVTightIsolationMVArun2017v1DBoldDMwLT2017", &tau_byVVTightIsolationMVArun2017v1DBoldDMwLT2017, &b_tau_byVVTightIsolationMVArun2017v1DBoldDMwLT2017);
   _tree->SetBranchAddress("tau_byVVLooseIsolationMVArun2017v2DBnewDMwLT2017", &tau_byVVLooseIsolationMVArun2017v2DBnewDMwLT2017, &b_tau_byVVLooseIsolationMVArun2017v2DBnewDMwLT2017);
   _tree->SetBranchAddress("tau_byVLooseIsolationMVArun2017v2DBnewDMwLT2017", &tau_byVLooseIsolationMVArun2017v2DBnewDMwLT2017, &b_tau_byVLooseIsolationMVArun2017v2DBnewDMwLT2017);
   _tree->SetBranchAddress("tau_byLooseIsolationMVArun2017v2DBnewDMwLT2017", &tau_byLooseIsolationMVArun2017v2DBnewDMwLT2017, &b_tau_byLooseIsolationMVArun2017v2DBnewDMwLT2017);
   _tree->SetBranchAddress("tau_byMediumIsolationMVArun2017v2DBnewDMwLT2017", &tau_byMediumIsolationMVArun2017v2DBnewDMwLT2017, &b_tau_byMediumIsolationMVArun2017v2DBnewDMwLT2017);
   _tree->SetBranchAddress("tau_byTightIsolationMVArun2017v2DBnewDMwLT2017", &tau_byTightIsolationMVArun2017v2DBnewDMwLT2017, &b_tau_byTightIsolationMVArun2017v2DBnewDMwLT2017);
   _tree->SetBranchAddress("tau_byVTightIsolationMVArun2017v2DBnewDMwLT2017", &tau_byVTightIsolationMVArun2017v2DBnewDMwLT2017, &b_tau_byVTightIsolationMVArun2017v2DBnewDMwLT2017);
   _tree->SetBranchAddress("tau_byVVTightIsolationMVArun2017v2DBnewDMwLT2017", &tau_byVVTightIsolationMVArun2017v2DBnewDMwLT2017, &b_tau_byVVTightIsolationMVArun2017v2DBnewDMwLT2017);
   _tree->SetBranchAddress("tau_byVVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017", &tau_byVVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017, &b_tau_byVVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017);
   _tree->SetBranchAddress("tau_byVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017", &tau_byVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017, &b_tau_byVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017);
   _tree->SetBranchAddress("tau_byLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017", &tau_byLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017, &b_tau_byLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017);
   _tree->SetBranchAddress("tau_byMediumIsolationMVArun2017v2DBoldDMdR0p3wLT2017", &tau_byMediumIsolationMVArun2017v2DBoldDMdR0p3wLT2017, &b_tau_byMediumIsolationMVArun2017v2DBoldDMdR0p3wLT2017);
   _tree->SetBranchAddress("tau_byTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017", &tau_byTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017, &b_tau_byTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017);
   _tree->SetBranchAddress("tau_byVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017", &tau_byVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017, &b_tau_byVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017);
   _tree->SetBranchAddress("tau_byVVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017", &tau_byVVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017, &b_tau_byVVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017);
   _tree->SetBranchAddress("tau_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017", &tau_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017, &b_tau_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017);
   _tree->SetBranchAddress("tau_byVLooseIsolationMVArun2017v2DBoldDMwLT2017", &tau_byVLooseIsolationMVArun2017v2DBoldDMwLT2017, &b_tau_byVLooseIsolationMVArun2017v2DBoldDMwLT2017);
   _tree->SetBranchAddress("tau_byLooseIsolationMVArun2017v2DBoldDMwLT2017", &tau_byLooseIsolationMVArun2017v2DBoldDMwLT2017, &b_tau_byLooseIsolationMVArun2017v2DBoldDMwLT2017);
   _tree->SetBranchAddress("tau_byMediumIsolationMVArun2017v2DBoldDMwLT2017", &tau_byMediumIsolationMVArun2017v2DBoldDMwLT2017, &b_tau_byMediumIsolationMVArun2017v2DBoldDMwLT2017);
   _tree->SetBranchAddress("tau_byTightIsolationMVArun2017v2DBoldDMwLT2017", &tau_byTightIsolationMVArun2017v2DBoldDMwLT2017, &b_tau_byTightIsolationMVArun2017v2DBoldDMwLT2017);
   _tree->SetBranchAddress("tau_byVTightIsolationMVArun2017v2DBoldDMwLT2017", &tau_byVTightIsolationMVArun2017v2DBoldDMwLT2017, &b_tau_byVTightIsolationMVArun2017v2DBoldDMwLT2017);
   _tree->SetBranchAddress("tau_byVVTightIsolationMVArun2017v2DBoldDMwLT2017", &tau_byVVTightIsolationMVArun2017v2DBoldDMwLT2017, &b_tau_byVVTightIsolationMVArun2017v2DBoldDMwLT2017);
   _tree->SetBranchAddress("tau_byVLooseIsolationMVArun2v1DBnewDMwLT2016", &tau_byVLooseIsolationMVArun2v1DBnewDMwLT2016, &b_tau_byVLooseIsolationMVArun2v1DBnewDMwLT2016);
   _tree->SetBranchAddress("tau_byLooseIsolationMVArun2v1DBnewDMwLT2016", &tau_byLooseIsolationMVArun2v1DBnewDMwLT2016, &b_tau_byLooseIsolationMVArun2v1DBnewDMwLT2016);
   _tree->SetBranchAddress("tau_byMediumIsolationMVArun2v1DBnewDMwLT2016", &tau_byMediumIsolationMVArun2v1DBnewDMwLT2016, &b_tau_byMediumIsolationMVArun2v1DBnewDMwLT2016);
   _tree->SetBranchAddress("tau_byTightIsolationMVArun2v1DBnewDMwLT2016", &tau_byTightIsolationMVArun2v1DBnewDMwLT2016, &b_tau_byTightIsolationMVArun2v1DBnewDMwLT2016);
   _tree->SetBranchAddress("tau_byVTightIsolationMVArun2v1DBnewDMwLT2016", &tau_byVTightIsolationMVArun2v1DBnewDMwLT2016, &b_tau_byVTightIsolationMVArun2v1DBnewDMwLT2016);
   _tree->SetBranchAddress("tau_byVVTightIsolationMVArun2v1DBnewDMwLT2016", &tau_byVVTightIsolationMVArun2v1DBnewDMwLT2016, &b_tau_byVVTightIsolationMVArun2v1DBnewDMwLT2016);
   _tree->SetBranchAddress("tau_byVLooseIsolationMVArun2v1DBoldDMwLT2016", &tau_byVLooseIsolationMVArun2v1DBoldDMwLT2016, &b_tau_byVLooseIsolationMVArun2v1DBoldDMwLT2016);
   _tree->SetBranchAddress("tau_byLooseIsolationMVArun2v1DBoldDMwLT2016", &tau_byLooseIsolationMVArun2v1DBoldDMwLT2016, &b_tau_byLooseIsolationMVArun2v1DBoldDMwLT2016);
   _tree->SetBranchAddress("tau_byMediumIsolationMVArun2v1DBoldDMwLT2016", &tau_byMediumIsolationMVArun2v1DBoldDMwLT2016, &b_tau_byMediumIsolationMVArun2v1DBoldDMwLT2016);
   _tree->SetBranchAddress("tau_byTightIsolationMVArun2v1DBoldDMwLT2016", &tau_byTightIsolationMVArun2v1DBoldDMwLT2016, &b_tau_byTightIsolationMVArun2v1DBoldDMwLT2016);
   _tree->SetBranchAddress("tau_byVTightIsolationMVArun2v1DBoldDMwLT2016", &tau_byVTightIsolationMVArun2v1DBoldDMwLT2016, &b_tau_byVTightIsolationMVArun2v1DBoldDMwLT2016);
   _tree->SetBranchAddress("tau_byVVTightIsolationMVArun2v1DBoldDMwLT2016", &tau_byVVTightIsolationMVArun2v1DBoldDMwLT2016, &b_tau_byVVTightIsolationMVArun2v1DBoldDMwLT2016);
   _tree->SetBranchAddress("trig_Flag_ecalBadCalibReduced", &trig_Flag_ecalBadCalibReduced, &b_trig_Flag_ecalBadCalibReduced);
   _tree->SetBranchAddress("trig_Flag_HBHENoiseFilter_accept", &trig_Flag_HBHENoiseFilter_accept, &b_trig_Flag_HBHENoiseFilter_accept);
   _tree->SetBranchAddress("trig_Flag_HBHENoiseIsoFilter_accept", &trig_Flag_HBHENoiseIsoFilter_accept, &b_trig_Flag_HBHENoiseIsoFilter_accept);
   _tree->SetBranchAddress("trig_Flag_CSCTightHaloFilter_accept", &trig_Flag_CSCTightHaloFilter_accept, &b_trig_Flag_CSCTightHaloFilter_accept);
   _tree->SetBranchAddress("trig_Flag_CSCTightHaloTrkMuUnvetoFilter_accept", &trig_Flag_CSCTightHaloTrkMuUnvetoFilter_accept, &b_trig_Flag_CSCTightHaloTrkMuUnvetoFilter_accept);
   _tree->SetBranchAddress("trig_Flag_CSCTightHalo2015Filter_accept", &trig_Flag_CSCTightHalo2015Filter_accept, &b_trig_Flag_CSCTightHalo2015Filter_accept);
   _tree->SetBranchAddress("trig_Flag_globalTightHalo2016Filter_accept", &trig_Flag_globalTightHalo2016Filter_accept, &b_trig_Flag_globalTightHalo2016Filter_accept);
   _tree->SetBranchAddress("trig_Flag_globalSuperTightHalo2016Filter_accept", &trig_Flag_globalSuperTightHalo2016Filter_accept, &b_trig_Flag_globalSuperTightHalo2016Filter_accept);
   _tree->SetBranchAddress("trig_Flag_HcalStripHaloFilter_accept", &trig_Flag_HcalStripHaloFilter_accept, &b_trig_Flag_HcalStripHaloFilter_accept);
   _tree->SetBranchAddress("trig_Flag_hcalLaserEventFilter_accept", &trig_Flag_hcalLaserEventFilter_accept, &b_trig_Flag_hcalLaserEventFilter_accept);
   _tree->SetBranchAddress("trig_Flag_EcalDeadCellTriggerPrimitiveFilter_accept", &trig_Flag_EcalDeadCellTriggerPrimitiveFilter_accept, &b_trig_Flag_EcalDeadCellTriggerPrimitiveFilter_accept);
   _tree->SetBranchAddress("trig_Flag_EcalDeadCellBoundaryEnergyFilter_accept", &trig_Flag_EcalDeadCellBoundaryEnergyFilter_accept, &b_trig_Flag_EcalDeadCellBoundaryEnergyFilter_accept);
   _tree->SetBranchAddress("trig_Flag_ecalBadCalibFilter_accept", &trig_Flag_ecalBadCalibFilter_accept, &b_trig_Flag_ecalBadCalibFilter_accept);
   _tree->SetBranchAddress("trig_Flag_goodVertices_accept", &trig_Flag_goodVertices_accept, &b_trig_Flag_goodVertices_accept);
   _tree->SetBranchAddress("trig_Flag_eeBadScFilter_accept", &trig_Flag_eeBadScFilter_accept, &b_trig_Flag_eeBadScFilter_accept);
   _tree->SetBranchAddress("trig_Flag_ecalLaserCorrFilter_accept", &trig_Flag_ecalLaserCorrFilter_accept, &b_trig_Flag_ecalLaserCorrFilter_accept);
   _tree->SetBranchAddress("trig_Flag_trkPOGFilters_accept", &trig_Flag_trkPOGFilters_accept, &b_trig_Flag_trkPOGFilters_accept);
   _tree->SetBranchAddress("trig_Flag_chargedHadronTrackResolutionFilter_accept", &trig_Flag_chargedHadronTrackResolutionFilter_accept, &b_trig_Flag_chargedHadronTrackResolutionFilter_accept);
   _tree->SetBranchAddress("trig_Flag_muonBadTrackFilter_accept", &trig_Flag_muonBadTrackFilter_accept, &b_trig_Flag_muonBadTrackFilter_accept);
   _tree->SetBranchAddress("trig_Flag_BadChargedCandidateFilter_accept", &trig_Flag_BadChargedCandidateFilter_accept, &b_trig_Flag_BadChargedCandidateFilter_accept);
   _tree->SetBranchAddress("trig_Flag_BadPFMuonFilter_accept", &trig_Flag_BadPFMuonFilter_accept, &b_trig_Flag_BadPFMuonFilter_accept);
   _tree->SetBranchAddress("trig_Flag_BadChargedCandidateSummer16Filter_accept", &trig_Flag_BadChargedCandidateSummer16Filter_accept, &b_trig_Flag_BadChargedCandidateSummer16Filter_accept);
   _tree->SetBranchAddress("trig_Flag_BadPFMuonSummer16Filter_accept", &trig_Flag_BadPFMuonSummer16Filter_accept, &b_trig_Flag_BadPFMuonSummer16Filter_accept);
   _tree->SetBranchAddress("trig_Flag_trkPOG_manystripclus53X_accept", &trig_Flag_trkPOG_manystripclus53X_accept, &b_trig_Flag_trkPOG_manystripclus53X_accept);
   _tree->SetBranchAddress("trig_Flag_trkPOG_toomanystripclus53X_accept", &trig_Flag_trkPOG_toomanystripclus53X_accept, &b_trig_Flag_trkPOG_toomanystripclus53X_accept);
   _tree->SetBranchAddress("trig_Flag_trkPOG_logErrorTooManyClusters_accept", &trig_Flag_trkPOG_logErrorTooManyClusters_accept, &b_trig_Flag_trkPOG_logErrorTooManyClusters_accept);
   _tree->SetBranchAddress("trig_Flag_METFilters_accept", &trig_Flag_METFilters_accept, &b_trig_Flag_METFilters_accept);
   //   _tree->SetBranchAddress("trig_raw2digi_step_accept", &trig_raw2digi_step_accept, &b_trig_raw2digi_step_accept);
   /* _tree->SetBranchAddress("trig_L1Reco_step_accept", &trig_L1Reco_step_accept, &b_trig_L1Reco_step_accept);
   _tree->SetBranchAddress("trig_reconstruction_step_accept", &trig_reconstruction_step_accept, &b_trig_reconstruction_step_accept);
   _tree->SetBranchAddress("trig_recosim_step_accept", &trig_recosim_step_accept, &b_trig_recosim_step_accept);
   _tree->SetBranchAddress("trig_eventinterpretaion_step_accept", &trig_eventinterpretaion_step_accept, &b_trig_eventinterpretaion_step_accept);
   _tree->SetBranchAddress("trig_HLT_DoubleEle25_CaloIdL_MW_accept", &trig_HLT_DoubleEle25_CaloIdL_MW_accept, &b_trig_HLT_DoubleEle25_CaloIdL_MW_accept);
   _tree->SetBranchAddress("trig_HLT_DoubleEle27_CaloIdL_MW_accept", &trig_HLT_DoubleEle27_CaloIdL_MW_accept, &b_trig_HLT_DoubleEle27_CaloIdL_MW_accept);
   _tree->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_MW_accept", &trig_HLT_DoubleEle33_CaloIdL_MW_accept, &b_trig_HLT_DoubleEle33_CaloIdL_MW_accept);
   _tree->SetBranchAddress("trig_HLT_DoubleEle24_eta2p1_WPTight_Gsf_accept", &trig_HLT_DoubleEle24_eta2p1_WPTight_Gsf_accept, &b_trig_HLT_DoubleEle24_eta2p1_WPTight_Gsf_accept);
   _tree->SetBranchAddress("trig_HLT_Ele27_Ele37_CaloIdL_MW_accept", &trig_HLT_Ele27_Ele37_CaloIdL_MW_accept, &b_trig_HLT_Ele27_Ele37_CaloIdL_MW_accept);
   _tree->SetBranchAddress("trig_HLT_Mu27_Ele37_CaloIdL_MW_accept", &trig_HLT_Mu27_Ele37_CaloIdL_MW_accept, &b_trig_HLT_Mu27_Ele37_CaloIdL_MW_accept);
   _tree->SetBranchAddress("trig_HLT_Mu37_Ele27_CaloIdL_MW_accept", &trig_HLT_Mu37_Ele27_CaloIdL_MW_accept, &b_trig_HLT_Mu37_Ele27_CaloIdL_MW_accept);
   _tree->SetBranchAddress("trig_HLT_Mu37_TkMu27_accept", &trig_HLT_Mu37_TkMu27_accept, &b_trig_HLT_Mu37_TkMu27_accept);
   _tree->SetBranchAddress("trig_HLT_DoubleMu4_3_Bs_accept", &trig_HLT_DoubleMu4_3_Bs_accept, &b_trig_HLT_DoubleMu4_3_Bs_accept);
   _tree->SetBranchAddress("trig_HLT_Ele20_WPTight_Gsf_accept", &trig_HLT_Ele20_WPTight_Gsf_accept, &b_trig_HLT_Ele20_WPTight_Gsf_accept);
   _tree->SetBranchAddress("trig_HLT_Ele15_WPLoose_Gsf_accept", &trig_HLT_Ele15_WPLoose_Gsf_accept, &b_trig_HLT_Ele15_WPLoose_Gsf_accept);
   _tree->SetBranchAddress("trig_HLT_Ele17_WPLoose_Gsf_accept", &trig_HLT_Ele17_WPLoose_Gsf_accept, &b_trig_HLT_Ele17_WPLoose_Gsf_accept);
   _tree->SetBranchAddress("trig_HLT_Ele20_WPLoose_Gsf_accept", &trig_HLT_Ele20_WPLoose_Gsf_accept, &b_trig_HLT_Ele20_WPLoose_Gsf_accept);
   _tree->SetBranchAddress("trig_HLT_Ele20_eta2p1_WPLoose_Gsf_accept", &trig_HLT_Ele20_eta2p1_WPLoose_Gsf_accept, &b_trig_HLT_Ele20_eta2p1_WPLoose_Gsf_accept);
   _tree->SetBranchAddress("trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_accept", &trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_accept, &b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_accept);
   _tree->SetBranchAddress("trig_HLT_Ele27_WPTight_Gsf_accept", &trig_HLT_Ele27_WPTight_Gsf_accept, &b_trig_HLT_Ele27_WPTight_Gsf_accept);
   _tree->SetBranchAddress("trig_HLT_Ele28_WPTight_Gsf_accept", &trig_HLT_Ele28_WPTight_Gsf_accept, &b_trig_HLT_Ele28_WPTight_Gsf_accept);
   _tree->SetBranchAddress("trig_HLT_Ele30_WPTight_Gsf_accept", &trig_HLT_Ele30_WPTight_Gsf_accept, &b_trig_HLT_Ele30_WPTight_Gsf_accept);
   */
   _tree->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_accept", &trig_HLT_Ele32_WPTight_Gsf_accept, &b_trig_HLT_Ele32_WPTight_Gsf_accept);
   _tree->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_accept", &trig_HLT_Ele35_WPTight_Gsf_accept, &b_trig_HLT_Ele35_WPTight_Gsf_accept);
   _tree->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_L1EGMT_accept", &trig_HLT_Ele35_WPTight_Gsf_L1EGMT_accept, &b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_accept);
   _tree->SetBranchAddress("trig_HLT_Ele38_WPTight_Gsf_accept", &trig_HLT_Ele38_WPTight_Gsf_accept, &b_trig_HLT_Ele38_WPTight_Gsf_accept);
   _tree->SetBranchAddress("trig_HLT_Ele40_WPTight_Gsf_accept", &trig_HLT_Ele40_WPTight_Gsf_accept, &b_trig_HLT_Ele40_WPTight_Gsf_accept);
   _tree->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_accept", &trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_accept, &b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_accept);
   /*
   _tree->SetBranchAddress("trig_HLT_IsoMu20_accept", &trig_HLT_IsoMu20_accept, &b_trig_HLT_IsoMu20_accept);
   _tree->SetBranchAddress("trig_HLT_IsoMu24_accept", &trig_HLT_IsoMu24_accept, &b_trig_HLT_IsoMu24_accept);
   _tree->SetBranchAddress("trig_HLT_IsoMu24_eta2p1_accept", &trig_HLT_IsoMu24_eta2p1_accept, &b_trig_HLT_IsoMu24_eta2p1_accept);
   _tree->SetBranchAddress("trig_HLT_IsoMu27_accept", &trig_HLT_IsoMu27_accept, &b_trig_HLT_IsoMu27_accept);
   _tree->SetBranchAddress("trig_HLT_IsoMu30_accept", &trig_HLT_IsoMu30_accept, &b_trig_HLT_IsoMu30_accept);
   _tree->SetBranchAddress("trig_HLT_L1SingleMu18_accept", &trig_HLT_L1SingleMu18_accept, &b_trig_HLT_L1SingleMu18_accept);
   _tree->SetBranchAddress("trig_HLT_L1SingleMu25_accept", &trig_HLT_L1SingleMu25_accept, &b_trig_HLT_L1SingleMu25_accept);
   _tree->SetBranchAddress("trig_HLT_L2Mu10_accept", &trig_HLT_L2Mu10_accept, &b_trig_HLT_L2Mu10_accept);
   _tree->SetBranchAddress("trig_HLT_L2Mu10_NoVertex_NoBPTX3BX_accept", &trig_HLT_L2Mu10_NoVertex_NoBPTX3BX_accept, &b_trig_HLT_L2Mu10_NoVertex_NoBPTX3BX_accept);
   _tree->SetBranchAddress("trig_HLT_L2Mu10_NoVertex_NoBPTX_accept", &trig_HLT_L2Mu10_NoVertex_NoBPTX_accept, &b_trig_HLT_L2Mu10_NoVertex_NoBPTX_accept);
   _tree->SetBranchAddress("trig_HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX_accept", &trig_HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX_accept, &b_trig_HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX_accept);
   _tree->SetBranchAddress("trig_HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX_accept", &trig_HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX_accept, &b_trig_HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX_accept);
   _tree->SetBranchAddress("trig_HLT_L2Mu50_accept", &trig_HLT_L2Mu50_accept, &b_trig_HLT_L2Mu50_accept);
   _tree->SetBranchAddress("trig_HLT_L2Mu23NoVtx_2Cha_accept", &trig_HLT_L2Mu23NoVtx_2Cha_accept, &b_trig_HLT_L2Mu23NoVtx_2Cha_accept);
   _tree->SetBranchAddress("trig_HLT_L2Mu23NoVtx_2Cha_CosmicSeed_accept", &trig_HLT_L2Mu23NoVtx_2Cha_CosmicSeed_accept, &b_trig_HLT_L2Mu23NoVtx_2Cha_CosmicSeed_accept);
   _tree->SetBranchAddress("trig_HLT_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4_accept", &trig_HLT_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4_accept, &b_trig_HLT_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4_accept);
   _tree->SetBranchAddress("trig_HLT_DoubleL2Mu30NoVtx_2Cha_Eta2p4_accept", &trig_HLT_DoubleL2Mu30NoVtx_2Cha_Eta2p4_accept, &b_trig_HLT_DoubleL2Mu30NoVtx_2Cha_Eta2p4_accept);
   _tree->SetBranchAddress("trig_HLT_DoubleL2Mu50_accept", &trig_HLT_DoubleL2Mu50_accept, &b_trig_HLT_DoubleL2Mu50_accept);
   _tree->SetBranchAddress("trig_HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed_accept", &trig_HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed_accept, &b_trig_HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed_accept);
   _tree->SetBranchAddress("trig_HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed_NoL2Matched_accept", &trig_HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed_NoL2Matched_accept, &b_trig_HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed_NoL2Matched_accept);
   _tree->SetBranchAddress("trig_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_accept", &trig_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_accept, &b_trig_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_accept);
   _tree->SetBranchAddress("trig_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_NoL2Matched_accept", &trig_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_NoL2Matched_accept, &b_trig_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_NoL2Matched_accept);
   _tree->SetBranchAddress("trig_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4_accept", &trig_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4_accept, &b_trig_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4_accept);
   _tree->SetBranchAddress("trig_HLT_DoubleL2Mu23NoVtx_2Cha_accept", &trig_HLT_DoubleL2Mu23NoVtx_2Cha_accept, &b_trig_HLT_DoubleL2Mu23NoVtx_2Cha_accept);
   _tree->SetBranchAddress("trig_HLT_DoubleL2Mu23NoVtx_2Cha_NoL2Matched_accept", &trig_HLT_DoubleL2Mu23NoVtx_2Cha_NoL2Matched_accept, &b_trig_HLT_DoubleL2Mu23NoVtx_2Cha_NoL2Matched_accept);
   _tree->SetBranchAddress("trig_HLT_DoubleL2Mu25NoVtx_2Cha_accept", &trig_HLT_DoubleL2Mu25NoVtx_2Cha_accept, &b_trig_HLT_DoubleL2Mu25NoVtx_2Cha_accept);
   _tree->SetBranchAddress("trig_HLT_DoubleL2Mu25NoVtx_2Cha_NoL2Matched_accept", &trig_HLT_DoubleL2Mu25NoVtx_2Cha_NoL2Matched_accept, &b_trig_HLT_DoubleL2Mu25NoVtx_2Cha_NoL2Matched_accept);
   _tree->SetBranchAddress("trig_HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4_accept", &trig_HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4_accept, &b_trig_HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4_accept);
   _tree->SetBranchAddress("trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_accept", &trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_accept, &b_trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_accept);
   _tree->SetBranchAddress("trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_accept", &trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_accept, &b_trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_accept);
   _tree->SetBranchAddress("trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_accept", &trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_accept, &b_trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_accept);
   _tree->SetBranchAddress("trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_accept", &trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_accept, &b_trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_accept);
   _tree->SetBranchAddress("trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_accept", &trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_accept, &b_trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_accept);
   _tree->SetBranchAddress("trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8_accept", &trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8_accept, &b_trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8_accept);
   _tree->SetBranchAddress("trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_accept", &trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_accept, &b_trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_accept);
   _tree->SetBranchAddress("trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8_accept", &trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8_accept, &b_trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8_accept);
   _tree->SetBranchAddress("trig_HLT_Mu25_TkMu0_Onia_accept", &trig_HLT_Mu25_TkMu0_Onia_accept, &b_trig_HLT_Mu25_TkMu0_Onia_accept);
   _tree->SetBranchAddress("trig_HLT_Mu30_TkMu0_Psi_accept", &trig_HLT_Mu30_TkMu0_Psi_accept, &b_trig_HLT_Mu30_TkMu0_Psi_accept);
   _tree->SetBranchAddress("trig_HLT_Mu20_TkMu0_Phi_accept", &trig_HLT_Mu20_TkMu0_Phi_accept, &b_trig_HLT_Mu20_TkMu0_Phi_accept);
   _tree->SetBranchAddress("trig_HLT_Mu25_TkMu0_Phi_accept", &trig_HLT_Mu25_TkMu0_Phi_accept, &b_trig_HLT_Mu25_TkMu0_Phi_accept);
   _tree->SetBranchAddress("trig_HLT_Mu12_accept", &trig_HLT_Mu12_accept, &b_trig_HLT_Mu12_accept);
   _tree->SetBranchAddress("trig_HLT_Mu15_accept", &trig_HLT_Mu15_accept, &b_trig_HLT_Mu15_accept);
   _tree->SetBranchAddress("trig_HLT_Mu20_accept", &trig_HLT_Mu20_accept, &b_trig_HLT_Mu20_accept);
   _tree->SetBranchAddress("trig_HLT_Mu27_accept", &trig_HLT_Mu27_accept, &b_trig_HLT_Mu27_accept);
   _tree->SetBranchAddress("trig_HLT_Mu50_accept", &trig_HLT_Mu50_accept, &b_trig_HLT_Mu50_accept);
   _tree->SetBranchAddress("trig_HLT_Mu55_accept", &trig_HLT_Mu55_accept, &b_trig_HLT_Mu55_accept);
   _tree->SetBranchAddress("trig_HLT_OldMu100_accept", &trig_HLT_OldMu100_accept, &b_trig_HLT_OldMu100_accept);
   _tree->SetBranchAddress("trig_HLT_TkMu100_accept", &trig_HLT_TkMu100_accept, &b_trig_HLT_TkMu100_accept);
   _tree->SetBranchAddress("trig_HLT_Photon300_NoHE_accept", &trig_HLT_Photon300_NoHE_accept, &b_trig_HLT_Photon300_NoHE_accept);
   _tree->SetBranchAddress("trig_HLT_Mu8_TrkIsoVVL_accept", &trig_HLT_Mu8_TrkIsoVVL_accept, &b_trig_HLT_Mu8_TrkIsoVVL_accept);
   _tree->SetBranchAddress("trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept", &trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept, &b_trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept);
   _tree->SetBranchAddress("trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept", &trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept, &b_trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept);
   _tree->SetBranchAddress("trig_HLT_Mu17_TrkIsoVVL_accept", &trig_HLT_Mu17_TrkIsoVVL_accept, &b_trig_HLT_Mu17_TrkIsoVVL_accept);
   _tree->SetBranchAddress("trig_HLT_Mu19_TrkIsoVVL_accept", &trig_HLT_Mu19_TrkIsoVVL_accept, &b_trig_HLT_Mu19_TrkIsoVVL_accept);
   _tree->SetBranchAddress("trig_HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL_accept", &trig_HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL_accept, &b_trig_HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL_accept);
   _tree->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept);
   _tree->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_accept", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_accept, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_accept);
   _tree->SetBranchAddress("trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept", &trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept, &b_trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept);
   _tree->SetBranchAddress("trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept", &trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept, &b_trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept);
   _tree->SetBranchAddress("trig_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept", &trig_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept, &b_trig_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept);
   _tree->SetBranchAddress("trig_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept", &trig_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept, &b_trig_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept);
   _tree->SetBranchAddress("trig_HLT_Photon20_accept", &trig_HLT_Photon20_accept, &b_trig_HLT_Photon20_accept);
   _tree->SetBranchAddress("trig_HLT_Photon33_accept", &trig_HLT_Photon33_accept, &b_trig_HLT_Photon33_accept);
   _tree->SetBranchAddress("trig_HLT_Photon50_accept", &trig_HLT_Photon50_accept, &b_trig_HLT_Photon50_accept);
   _tree->SetBranchAddress("trig_HLT_Photon75_accept", &trig_HLT_Photon75_accept, &b_trig_HLT_Photon75_accept);
   _tree->SetBranchAddress("trig_HLT_Photon90_accept", &trig_HLT_Photon90_accept, &b_trig_HLT_Photon90_accept);
   _tree->SetBranchAddress("trig_HLT_Photon120_accept", &trig_HLT_Photon120_accept, &b_trig_HLT_Photon120_accept);
   _tree->SetBranchAddress("trig_HLT_Photon150_accept", &trig_HLT_Photon150_accept, &b_trig_HLT_Photon150_accept);
   _tree->SetBranchAddress("trig_HLT_Photon175_accept", &trig_HLT_Photon175_accept, &b_trig_HLT_Photon175_accept);
   */
   _tree->SetBranchAddress("trig_HLT_Photon200_accept", &trig_HLT_Photon200_accept, &b_trig_HLT_Photon200_accept);
   /* _tree->SetBranchAddress("trig_HLT_Photon100EB_TightID_TightIso_accept", &trig_HLT_Photon100EB_TightID_TightIso_accept, &b_trig_HLT_Photon100EB_TightID_TightIso_accept);
   _tree->SetBranchAddress("trig_HLT_Photon110EB_TightID_TightIso_accept", &trig_HLT_Photon110EB_TightID_TightIso_accept, &b_trig_HLT_Photon110EB_TightID_TightIso_accept);
   _tree->SetBranchAddress("trig_HLT_Photon120EB_TightID_TightIso_accept", &trig_HLT_Photon120EB_TightID_TightIso_accept, &b_trig_HLT_Photon120EB_TightID_TightIso_accept);
   _tree->SetBranchAddress("trig_HLT_Photon100EBHE10_accept", &trig_HLT_Photon100EBHE10_accept, &b_trig_HLT_Photon100EBHE10_accept);
   _tree->SetBranchAddress("trig_HLT_Photon100EEHE10_accept", &trig_HLT_Photon100EEHE10_accept, &b_trig_HLT_Photon100EEHE10_accept);
   _tree->SetBranchAddress("trig_HLT_Photon100EE_TightID_TightIso_accept", &trig_HLT_Photon100EE_TightID_TightIso_accept, &b_trig_HLT_Photon100EE_TightID_TightIso_accept);
   _tree->SetBranchAddress("trig_HLT_Photon50_R9Id90_HE10_IsoM_accept", &trig_HLT_Photon50_R9Id90_HE10_IsoM_accept, &b_trig_HLT_Photon50_R9Id90_HE10_IsoM_accept);
   _tree->SetBranchAddress("trig_HLT_Photon75_R9Id90_HE10_IsoM_accept", &trig_HLT_Photon75_R9Id90_HE10_IsoM_accept, &b_trig_HLT_Photon75_R9Id90_HE10_IsoM_accept);
   _tree->SetBranchAddress("trig_HLT_Photon90_R9Id90_HE10_IsoM_accept", &trig_HLT_Photon90_R9Id90_HE10_IsoM_accept, &b_trig_HLT_Photon90_R9Id90_HE10_IsoM_accept);
   _tree->SetBranchAddress("trig_HLT_Photon120_R9Id90_HE10_IsoM_accept", &trig_HLT_Photon120_R9Id90_HE10_IsoM_accept, &b_trig_HLT_Photon120_R9Id90_HE10_IsoM_accept);
   _tree->SetBranchAddress("trig_HLT_Photon165_R9Id90_HE10_IsoM_accept", &trig_HLT_Photon165_R9Id90_HE10_IsoM_accept, &b_trig_HLT_Photon165_R9Id90_HE10_IsoM_accept);
   _tree->SetBranchAddress("trig_HLT_Photon35_TwoProngs35_accept", &trig_HLT_Photon35_TwoProngs35_accept, &b_trig_HLT_Photon35_TwoProngs35_accept);
   _tree->SetBranchAddress("trig_HLT_IsoMu24_TwoProngs35_accept", &trig_HLT_IsoMu24_TwoProngs35_accept, &b_trig_HLT_IsoMu24_TwoProngs35_accept);
   _tree->SetBranchAddress("trig_HLT_DoubleMu43NoFiltersNoVtx_accept", &trig_HLT_DoubleMu43NoFiltersNoVtx_accept, &b_trig_HLT_DoubleMu43NoFiltersNoVtx_accept);
   _tree->SetBranchAddress("trig_HLT_DoubleMu48NoFiltersNoVtx_accept", &trig_HLT_DoubleMu48NoFiltersNoVtx_accept, &b_trig_HLT_DoubleMu48NoFiltersNoVtx_accept);
   _tree->SetBranchAddress("trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_accept", &trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_accept, &b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_accept);
   _tree->SetBranchAddress("trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_accept", &trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_accept, &b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_accept);
   _tree->SetBranchAddress("trig_HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL_accept", &trig_HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL_accept, &b_trig_HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL_accept);
   _tree->SetBranchAddress("trig_HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL_accept", &trig_HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL_accept, &b_trig_HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL_accept);
   _tree->SetBranchAddress("trig_HLT_DoubleMu33NoFiltersNoVtxDisplaced_accept", &trig_HLT_DoubleMu33NoFiltersNoVtxDisplaced_accept, &b_trig_HLT_DoubleMu33NoFiltersNoVtxDisplaced_accept);
   _tree->SetBranchAddress("trig_HLT_DoubleMu40NoFiltersNoVtxDisplaced_accept", &trig_HLT_DoubleMu40NoFiltersNoVtxDisplaced_accept, &b_trig_HLT_DoubleMu40NoFiltersNoVtxDisplaced_accept);
   _tree->SetBranchAddress("trig_HLT_DoubleMu20_7_Mass0to30_L1_DM4_accept", &trig_HLT_DoubleMu20_7_Mass0to30_L1_DM4_accept, &b_trig_HLT_DoubleMu20_7_Mass0to30_L1_DM4_accept);
   _tree->SetBranchAddress("trig_HLT_DoubleMu20_7_Mass0to30_L1_DM4EG_accept", &trig_HLT_DoubleMu20_7_Mass0to30_L1_DM4EG_accept, &b_trig_HLT_DoubleMu20_7_Mass0to30_L1_DM4EG_accept);
   _tree->SetBranchAddress("trig_HLT_DoubleMu20_7_Mass0to30_Photon23_accept", &trig_HLT_DoubleMu20_7_Mass0to30_Photon23_accept, &b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_accept);
   _tree->SetBranchAddress("trig_HLT_Dimuon14_Phi_Barrel_Seagulls_accept", &trig_HLT_Dimuon14_Phi_Barrel_Seagulls_accept, &b_trig_HLT_Dimuon14_Phi_Barrel_Seagulls_accept);
   _tree->SetBranchAddress("trig_HLT_Dimuon24_Phi_noCorrL1_accept", &trig_HLT_Dimuon24_Phi_noCorrL1_accept, &b_trig_HLT_Dimuon24_Phi_noCorrL1_accept);
   _tree->SetBranchAddress("trig_HLT_DoubleIsoMu20_eta2p1_accept", &trig_HLT_DoubleIsoMu20_eta2p1_accept, &b_trig_HLT_DoubleIsoMu20_eta2p1_accept);
   _tree->SetBranchAddress("trig_HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx_accept", &trig_HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx_accept, &b_trig_HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx_accept);
   _tree->SetBranchAddress("trig_HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx_accept", &trig_HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx_accept, &b_trig_HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx_accept);
   _tree->SetBranchAddress("trig_HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx_accept", &trig_HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx_accept, &b_trig_HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx_accept);
   _tree->SetBranchAddress("trig_HLT_Mu8_accept", &trig_HLT_Mu8_accept, &b_trig_HLT_Mu8_accept);
   _tree->SetBranchAddress("trig_HLT_Mu17_accept", &trig_HLT_Mu17_accept, &b_trig_HLT_Mu17_accept);
   _tree->SetBranchAddress("trig_HLT_Mu19_accept", &trig_HLT_Mu19_accept, &b_trig_HLT_Mu19_accept);
   _tree->SetBranchAddress("trig_HLT_Mu17_Photon30_IsoCaloId_accept", &trig_HLT_Mu17_Photon30_IsoCaloId_accept, &b_trig_HLT_Mu17_Photon30_IsoCaloId_accept);
   */
   _tree->SetBranchAddress("trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_accept", &trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_accept, &b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_accept);

   /*
   _tree->SetBranchAddress("trig_HLT_Ele135_CaloIdVT_GsfTrkIdT_accept", &trig_HLT_Ele135_CaloIdVT_GsfTrkIdT_accept, &b_trig_HLT_Ele135_CaloIdVT_GsfTrkIdT_accept);
   _tree->SetBranchAddress("trig_HLT_Ele145_CaloIdVT_GsfTrkIdT_accept", &trig_HLT_Ele145_CaloIdVT_GsfTrkIdT_accept, &b_trig_HLT_Ele145_CaloIdVT_GsfTrkIdT_accept);
   _tree->SetBranchAddress("trig_HLT_Ele200_CaloIdVT_GsfTrkIdT_accept", &trig_HLT_Ele200_CaloIdVT_GsfTrkIdT_accept, &b_trig_HLT_Ele200_CaloIdVT_GsfTrkIdT_accept);
   _tree->SetBranchAddress("trig_HLT_Ele250_CaloIdVT_GsfTrkIdT_accept", &trig_HLT_Ele250_CaloIdVT_GsfTrkIdT_accept, &b_trig_HLT_Ele250_CaloIdVT_GsfTrkIdT_accept);
   _tree->SetBranchAddress("trig_HLT_Ele300_CaloIdVT_GsfTrkIdT_accept", &trig_HLT_Ele300_CaloIdVT_GsfTrkIdT_accept, &b_trig_HLT_Ele300_CaloIdVT_GsfTrkIdT_accept);
   _tree->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_accept", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_accept, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_accept);
   _tree->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_accept", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_accept, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_accept);
   _tree->SetBranchAddress("trig_HLT_SinglePhoton10_Eta3p1ForPPRef_accept", &trig_HLT_SinglePhoton10_Eta3p1ForPPRef_accept, &b_trig_HLT_SinglePhoton10_Eta3p1ForPPRef_accept);
   _tree->SetBranchAddress("trig_HLT_SinglePhoton20_Eta3p1ForPPRef_accept", &trig_HLT_SinglePhoton20_Eta3p1ForPPRef_accept, &b_trig_HLT_SinglePhoton20_Eta3p1ForPPRef_accept);
   _tree->SetBranchAddress("trig_HLT_SinglePhoton30_Eta3p1ForPPRef_accept", &trig_HLT_SinglePhoton30_Eta3p1ForPPRef_accept, &b_trig_HLT_SinglePhoton30_Eta3p1ForPPRef_accept);
   _tree->SetBranchAddress("trig_HLT_Photon20_HoverELoose_accept", &trig_HLT_Photon20_HoverELoose_accept, &b_trig_HLT_Photon20_HoverELoose_accept);
   _tree->SetBranchAddress("trig_HLT_Photon30_HoverELoose_accept", &trig_HLT_Photon30_HoverELoose_accept, &b_trig_HLT_Photon30_HoverELoose_accept);
   _tree->SetBranchAddress("trig_HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142_accept", &trig_HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142_accept, &b_trig_HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142_accept);
   _tree->SetBranchAddress("trig_AlCa_RPCMuonNormalisation_accept", &trig_AlCa_RPCMuonNormalisation_accept, &b_trig_AlCa_RPCMuonNormalisation_accept);
   _tree->SetBranchAddress("trig_MC_DoubleEle5_CaloIdL_MW_accept", &trig_MC_DoubleEle5_CaloIdL_MW_accept, &b_trig_MC_DoubleEle5_CaloIdL_MW_accept);
   _tree->SetBranchAddress("trig_MC_Ele5_WPTight_Gsf_accept", &trig_MC_Ele5_WPTight_Gsf_accept, &b_trig_MC_Ele5_WPTight_Gsf_accept);
   _tree->SetBranchAddress("trig_MC_Ele15_Ele10_CaloIdL_TrackIdL_IsoVL_DZ_accept", &trig_MC_Ele15_Ele10_CaloIdL_TrackIdL_IsoVL_DZ_accept, &b_trig_MC_Ele15_Ele10_CaloIdL_TrackIdL_IsoVL_DZ_accept);
   _tree->SetBranchAddress("trig_MC_IsoMu_accept", &trig_MC_IsoMu_accept, &b_trig_MC_IsoMu_accept);
   _tree->SetBranchAddress("trig_MC_DoubleMu_TrkIsoVVL_DZ_accept", &trig_MC_DoubleMu_TrkIsoVVL_DZ_accept, &b_trig_MC_DoubleMu_TrkIsoVVL_DZ_accept);
   _tree->SetBranchAddress("trig_MC_DoubleMuNoFiltersNoVtx_accept", &trig_MC_DoubleMuNoFiltersNoVtx_accept, &b_trig_MC_DoubleMuNoFiltersNoVtx_accept);
   _tree->SetBranchAddress("trig_HLT_Mu18_Mu9_SameSign_accept", &trig_HLT_Mu18_Mu9_SameSign_accept, &b_trig_HLT_Mu18_Mu9_SameSign_accept);
   _tree->SetBranchAddress("trig_HLT_Mu18_Mu9_SameSign_DZ_accept", &trig_HLT_Mu18_Mu9_SameSign_DZ_accept, &b_trig_HLT_Mu18_Mu9_SameSign_DZ_accept);
   _tree->SetBranchAddress("trig_HLT_Mu18_Mu9_accept", &trig_HLT_Mu18_Mu9_accept, &b_trig_HLT_Mu18_Mu9_accept);
   _tree->SetBranchAddress("trig_HLT_Mu18_Mu9_DZ_accept", &trig_HLT_Mu18_Mu9_DZ_accept, &b_trig_HLT_Mu18_Mu9_DZ_accept);
   _tree->SetBranchAddress("trig_HLT_Mu20_Mu10_SameSign_accept", &trig_HLT_Mu20_Mu10_SameSign_accept, &b_trig_HLT_Mu20_Mu10_SameSign_accept);
   _tree->SetBranchAddress("trig_HLT_Mu20_Mu10_SameSign_DZ_accept", &trig_HLT_Mu20_Mu10_SameSign_DZ_accept, &b_trig_HLT_Mu20_Mu10_SameSign_DZ_accept);
   _tree->SetBranchAddress("trig_HLT_Mu20_Mu10_accept", &trig_HLT_Mu20_Mu10_accept, &b_trig_HLT_Mu20_Mu10_accept);
   _tree->SetBranchAddress("trig_HLT_Mu20_Mu10_DZ_accept", &trig_HLT_Mu20_Mu10_DZ_accept, &b_trig_HLT_Mu20_Mu10_DZ_accept);
   _tree->SetBranchAddress("trig_HLT_Mu23_Mu12_SameSign_accept", &trig_HLT_Mu23_Mu12_SameSign_accept, &b_trig_HLT_Mu23_Mu12_SameSign_accept);
   _tree->SetBranchAddress("trig_HLT_Mu23_Mu12_SameSign_DZ_accept", &trig_HLT_Mu23_Mu12_SameSign_DZ_accept, &b_trig_HLT_Mu23_Mu12_SameSign_DZ_accept);
   _tree->SetBranchAddress("trig_HLT_Mu23_Mu12_accept", &trig_HLT_Mu23_Mu12_accept, &b_trig_HLT_Mu23_Mu12_accept);
   _tree->SetBranchAddress("trig_HLT_Mu23_Mu12_DZ_accept", &trig_HLT_Mu23_Mu12_DZ_accept, &b_trig_HLT_Mu23_Mu12_DZ_accept);
   _tree->SetBranchAddress("trig_HLT_TrkMu6NoFiltersNoVtx_accept", &trig_HLT_TrkMu6NoFiltersNoVtx_accept, &b_trig_HLT_TrkMu6NoFiltersNoVtx_accept);
   _tree->SetBranchAddress("trig_HLT_TrkMu16NoFiltersNoVtx_accept", &trig_HLT_TrkMu16NoFiltersNoVtx_accept, &b_trig_HLT_TrkMu16NoFiltersNoVtx_accept);
   _tree->SetBranchAddress("trig_HLT_DoubleTrkMu_16_6_NoFiltersNoVtx_accept", &trig_HLT_DoubleTrkMu_16_6_NoFiltersNoVtx_accept, &b_trig_HLT_DoubleTrkMu_16_6_NoFiltersNoVtx_accept);

   */

}

Int_t IIHEAnalysis::GetEntry(int entry)
{
	return _tree->GetEntry(entry);
}

Long64_t IIHEAnalysis::GetEntries()
{
	return _tree->GetEntries();
}

TTree* IIHEAnalysis::GetTree()
{
	return _tree;
}

#endif // #ifdef IIHEAnalysis_cxx
