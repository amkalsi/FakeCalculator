//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Apr 30 14:07:12 2019 by ROOT version 6.06/01
// from TTree IIHEAnalysis/IIHEAnalysis
// found on file: outfile_956.root
//////////////////////////////////////////////////////////

#ifndef Nochecl_h
#define Nochecl_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"

class Nochecl {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
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
   vector<int>     *tau_decayModeFinding;
   vector<int>     *tau_decayModeFindingNewDMs;
   vector<int>     *tau_againstMuonLoose3;
   vector<int>     *tau_againstMuonTight3;
   vector<int>     *tau_byLooseCombinedIsolationDeltaBetaCorr3Hits;
   vector<int>     *tau_byMediumCombinedIsolationDeltaBetaCorr3Hits;
   vector<int>     *tau_byTightCombinedIsolationDeltaBetaCorr3Hits;
   vector<int>     *tau_byVLooseIsolationMVArun2v1DBoldDMwLT;
   vector<int>     *tau_byLooseIsolationMVArun2v1DBoldDMwLT;
   vector<int>     *tau_byMediumIsolationMVArun2v1DBoldDMwLT;
   vector<int>     *tau_byTightIsolationMVArun2v1DBoldDMwLT;
   vector<int>     *tau_byVTightIsolationMVArun2v1DBoldDMwLT;
   vector<int>     *tau_byVVTightIsolationMVArun2v1DBoldDMwLT;
   vector<int>     *tau_byVLooseIsolationMVArun2v1DBnewDMwLT;
   vector<int>     *tau_byLooseIsolationMVArun2v1DBnewDMwLT;
   vector<int>     *tau_byMediumIsolationMVArun2v1DBnewDMwLT;
   vector<int>     *tau_byTightIsolationMVArun2v1DBnewDMwLT;
   vector<int>     *tau_byVTightIsolationMVArun2v1DBnewDMwLT;
   vector<int>     *tau_byVVTightIsolationMVArun2v1DBnewDMwLT;
   vector<int>     *tau_byVLooseIsolationMVArun2v1PWoldDMwLT;
   vector<int>     *tau_byLooseIsolationMVArun2v1PWoldDMwLT;
   vector<int>     *tau_byMediumIsolationMVArun2v1PWoldDMwLT;
   vector<int>     *tau_byTightIsolationMVArun2v1PWoldDMwLT;
   vector<int>     *tau_byVTightIsolationMVArun2v1PWoldDMwLT;
   vector<int>     *tau_byVVTightIsolationMVArun2v1PWoldDMwLT;
   vector<int>     *tau_byVLooseIsolationMVArun2v1PWnewDMwLT;
   vector<int>     *tau_byLooseIsolationMVArun2v1PWnewDMwLT;
   vector<int>     *tau_byMediumIsolationMVArun2v1PWnewDMwLT;
   vector<int>     *tau_byTightIsolationMVArun2v1PWnewDMwLT;
   vector<int>     *tau_byVTightIsolationMVArun2v1PWnewDMwLT;
   vector<int>     *tau_byVVTightIsolationMVArun2v1PWnewDMwLT;
   vector<int>     *tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT;
   vector<int>     *tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT;
   vector<int>     *tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT;
   vector<int>     *tau_byTightIsolationMVArun2v1DBdR03oldDMwLT;
   vector<int>     *tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT;
   vector<int>     *tau_byVVTightIsolationMVArun2v1DBdR03oldDMwLT;
   vector<int>     *tau_byVLooseIsolationMVArun2v1PWdR03oldDMwLT;
   vector<int>     *tau_byLooseIsolationMVArun2v1PWdR03oldDMwLT;
   vector<int>     *tau_byMediumIsolationMVArun2v1PWdR03oldDMwLT;
   vector<int>     *tau_byTightIsolationMVArun2v1PWdR03oldDMwLT;
   vector<int>     *tau_byVTightIsolationMVArun2v1PWdR03oldDMwLT;
   vector<int>     *tau_byVVTightIsolationMVArun2v1PWdR03oldDMwLT;
   vector<int>     *tau_againstElectronVLooseMVA6;
   vector<int>     *tau_againstElectronLooseMVA6;
   vector<int>     *tau_againstElectronMediumMVA6;
   vector<int>     *tau_againstElectronTightMVA6;
   vector<int>     *tau_againstElectronVTightMVA6;
   vector<int>     *tau_chargedIsoPtSum;
   vector<int>     *tau_neutralIsoPtSum;
   vector<int>     *tau_puCorrPtSum;
   vector<int>     *tau_footprintCorrection;
   vector<int>     *tau_neutralIsoPtSumWeight;
   vector<int>     *tau_photonPtSumOutsideSignalCone;
   vector<int>     *tau_byPhotonPtSumOutsideSignalCone;
   vector<int>     *tau_footprintCorrectiondR03;
   vector<int>     *tau_chargedIsoPtSumdR03;
   vector<int>     *tau_neutralIsoPtSumWeightdR03;
   vector<int>     *tau_neutralIsoPtSumdR03;
   vector<int>     *tau_photonPtSumOutsideSignalConedR03;
   vector<int>     *tau_PFChargedHadIso;
   vector<int>     *tau_PFNeutralHadIso;
   vector<int>     *tau_PFPhotonIso;
   vector<int>     *tau_leadChargedParticlePt;
   vector<int>     *tau_trackRefPt;
   vector<int>     *tau_lead_dxy;
   vector<int>     *tau_lead_dz;
   vector<int>     *tau_dxy_Sig;
   vector<int>     *tau_flightLengthSig;
   vector<int>     *tau_ip3d;
   vector<int>     *tau_ip3d_Sig;
   vector<int>     *tau_decayDistX;
   vector<int>     *tau_decayDistY;
   vector<int>     *tau_decayDistZ;
   vector<int>     *tau_decayDistMag;
   vector<int>     *tau_nPhoton;
   vector<int>     *tau_ptWeightedDetaStrip;
   vector<int>     *tau_ptWeightedDphiStrip;
   vector<int>     *tau_ptWeightedDrSignal;
   vector<int>     *tau_ptWeightedDrIsolation;
   vector<int>     *tau_leadingTrackChi2;
   vector<int>     *tau_eRatio;
   vector<int>     *tau_gjAngleDiff;
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
   vector<int>     *tau_byVVLooseIsolationMVArun2017v1DBoldDMwLT2017;
   vector<int>     *tau_byVLooseIsolationMVArun2017v1DBoldDMwLT2017;
   vector<int>     *tau_byLooseIsolationMVArun2017v1DBoldDMwLT2017;
   vector<int>     *tau_byMediumIsolationMVArun2017v1DBoldDMwLT2017;
   vector<int>     *tau_byTightIsolationMVArun2017v1DBoldDMwLT2017;
   vector<int>     *tau_byVTightIsolationMVArun2017v1DBoldDMwLT2017;
   vector<int>     *tau_byVVTightIsolationMVArun2017v1DBoldDMwLT2017;
   vector<int>     *tau_byVVLooseIsolationMVArun2017v2DBnewDMwLT2017;
   vector<int>     *tau_byVLooseIsolationMVArun2017v2DBnewDMwLT2017;
   vector<int>     *tau_byLooseIsolationMVArun2017v2DBnewDMwLT2017;
   vector<int>     *tau_byMediumIsolationMVArun2017v2DBnewDMwLT2017;
   vector<int>     *tau_byTightIsolationMVArun2017v2DBnewDMwLT2017;
   vector<int>     *tau_byVTightIsolationMVArun2017v2DBnewDMwLT2017;
   vector<int>     *tau_byVVTightIsolationMVArun2017v2DBnewDMwLT2017;
   vector<int>     *tau_byVVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017;
   vector<int>     *tau_byVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017;
   vector<int>     *tau_byLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017;
   vector<int>     *tau_byMediumIsolationMVArun2017v2DBoldDMdR0p3wLT2017;
   vector<int>     *tau_byTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017;
   vector<int>     *tau_byVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017;
   vector<int>     *tau_byVVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017;
   vector<int>     *tau_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017;
   vector<int>     *tau_byVLooseIsolationMVArun2017v2DBoldDMwLT2017;
   vector<int>     *tau_byLooseIsolationMVArun2017v2DBoldDMwLT2017;
   vector<int>     *tau_byMediumIsolationMVArun2017v2DBoldDMwLT2017;
   vector<int>     *tau_byTightIsolationMVArun2017v2DBoldDMwLT2017;
   vector<int>     *tau_byVTightIsolationMVArun2017v2DBoldDMwLT2017;
   vector<int>     *tau_byVVTightIsolationMVArun2017v2DBoldDMwLT2017;
   vector<int>     *tau_byVLooseIsolationMVArun2v1DBnewDMwLT2016;
   vector<int>     *tau_byLooseIsolationMVArun2v1DBnewDMwLT2016;
   vector<int>     *tau_byMediumIsolationMVArun2v1DBnewDMwLT2016;
   vector<int>     *tau_byTightIsolationMVArun2v1DBnewDMwLT2016;
   vector<int>     *tau_byVTightIsolationMVArun2v1DBnewDMwLT2016;
   vector<int>     *tau_byVVTightIsolationMVArun2v1DBnewDMwLT2016;
   vector<int>     *tau_byVLooseIsolationMVArun2v1DBoldDMwLT2016;
   vector<int>     *tau_byLooseIsolationMVArun2v1DBoldDMwLT2016;
   vector<int>     *tau_byMediumIsolationMVArun2v1DBoldDMwLT2016;
   vector<int>     *tau_byTightIsolationMVArun2v1DBoldDMwLT2016;
   vector<int>     *tau_byVTightIsolationMVArun2v1DBoldDMwLT2016;
   vector<int>     *tau_byVVTightIsolationMVArun2v1DBoldDMwLT2016;
   Int_t           trig_Flag_ecalBadCalibReduced;
   Int_t           trig_raw2digi_step_accept;
   Int_t           trig_L1Reco_step_accept;
   Int_t           trig_reconstruction_step_accept;
   Int_t           trig_EXOMONOPOLEPath_accept;
   Int_t           trig_ZElectronPathGsf_accept;
   Int_t           trig_pathALCARECOEcalESAlign_accept;
   Int_t           trig_pathALCARECOEcalUncalWElectron_accept;
   Int_t           trig_pathALCARECOHcalCalIsoTrkFilter_accept;
   Int_t           trig_pathALCARECOHcalCalIterativePhiSym_accept;
   Int_t           trig_pathALCARECOEcalUncalZElectron_accept;
   Int_t           trig_pathALCARECOEcalUncalZSCElectron_accept;
   Int_t           trig_eventinterpretaion_step_accept;
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
   Int_t           trig_HLT_DoubleEle25_CaloIdL_MW_accept;
   vector<float>   *trig_HLT_DoubleEle25_CaloIdL_MW_hltL1sSingleAndDoubleEGNonIsoOrWithEG26WithJetAndTau_eta;
   vector<float>   *trig_HLT_DoubleEle25_CaloIdL_MW_hltL1sSingleAndDoubleEGNonIsoOrWithEG26WithJetAndTau_phi;
   vector<float>   *trig_HLT_DoubleEle25_CaloIdL_MW_hltL1sSingleAndDoubleEGNonIsoOrWithEG26WithJetAndTau_et;
   vector<float>   *trig_HLT_DoubleEle25_CaloIdL_MW_hltEGL1SingleAndDoubleEGNonIsoOrWithEG26WithJetAndTauFilter_eta;
   vector<float>   *trig_HLT_DoubleEle25_CaloIdL_MW_hltEGL1SingleAndDoubleEGNonIsoOrWithEG26WithJetAndTauFilter_phi;
   vector<float>   *trig_HLT_DoubleEle25_CaloIdL_MW_hltEGL1SingleAndDoubleEGNonIsoOrWithEG26WithJetAndTauFilter_et;
   vector<float>   *trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25EtFilter_eta;
   vector<float>   *trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25EtFilter_phi;
   vector<float>   *trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25EtFilter_et;
   vector<float>   *trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25HEFilter_eta;
   vector<float>   *trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25HEFilter_phi;
   vector<float>   *trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25HEFilter_et;
   vector<float>   *trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25CaloIdLClusterShapeFilter_eta;
   vector<float>   *trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25CaloIdLClusterShapeFilter_phi;
   vector<float>   *trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25CaloIdLClusterShapeFilter_et;
   vector<float>   *trig_HLT_DoubleEle25_CaloIdL_MW_hltEle25CaloIdLPixelMatchFilter_eta;
   vector<float>   *trig_HLT_DoubleEle25_CaloIdL_MW_hltEle25CaloIdLPixelMatchFilter_phi;
   vector<float>   *trig_HLT_DoubleEle25_CaloIdL_MW_hltEle25CaloIdLPixelMatchFilter_et;
   vector<float>   *trig_HLT_DoubleEle25_CaloIdL_MW_hltEle25CaloIdLMWPMS2Filter_eta;
   vector<float>   *trig_HLT_DoubleEle25_CaloIdL_MW_hltEle25CaloIdLMWPMS2Filter_phi;
   vector<float>   *trig_HLT_DoubleEle25_CaloIdL_MW_hltEle25CaloIdLMWPMS2Filter_et;
   vector<float>   *trig_HLT_DoubleEle25_CaloIdL_MW_hltEgammaCandidatesWrapperUnseeded_eta;
   vector<float>   *trig_HLT_DoubleEle25_CaloIdL_MW_hltEgammaCandidatesWrapperUnseeded_phi;
   vector<float>   *trig_HLT_DoubleEle25_CaloIdL_MW_hltEgammaCandidatesWrapperUnseeded_et;
   vector<float>   *trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25EtUnseededFilter_eta;
   vector<float>   *trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25EtUnseededFilter_phi;
   vector<float>   *trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25EtUnseededFilter_et;
   vector<float>   *trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25HEUnseededFilter_eta;
   vector<float>   *trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25HEUnseededFilter_phi;
   vector<float>   *trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25HEUnseededFilter_et;
   vector<float>   *trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25CaloIdLClusterShapeUnseededFilter_eta;
   vector<float>   *trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25CaloIdLClusterShapeUnseededFilter_phi;
   vector<float>   *trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25CaloIdLClusterShapeUnseededFilter_et;
   vector<float>   *trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEle25CaloIdLPixelMatchUnseededFilter_eta;
   vector<float>   *trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEle25CaloIdLPixelMatchUnseededFilter_phi;
   vector<float>   *trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEle25CaloIdLPixelMatchUnseededFilter_et;
   vector<float>   *trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEle25CaloIdLMWPMS2UnseededFilter_eta;
   vector<float>   *trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEle25CaloIdLMWPMS2UnseededFilter_phi;
   vector<float>   *trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEle25CaloIdLMWPMS2UnseededFilter_et;
   Int_t           trig_HLT_DoubleEle27_CaloIdL_MW_accept;
   Int_t           trig_HLT_DoubleEle33_CaloIdL_MW_accept;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_MW_hltL1sSingleAndDoubleEGNonIsoOrWithEG26WithJetAndTau_eta;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_MW_hltL1sSingleAndDoubleEGNonIsoOrWithEG26WithJetAndTau_phi;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_MW_hltL1sSingleAndDoubleEGNonIsoOrWithEG26WithJetAndTau_et;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_MW_hltEGL1SingleAndDoubleEGNonIsoOrWithEG26WithJetAndTauFilter_eta;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_MW_hltEGL1SingleAndDoubleEGNonIsoOrWithEG26WithJetAndTauFilter_phi;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_MW_hltEGL1SingleAndDoubleEGNonIsoOrWithEG26WithJetAndTauFilter_et;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33EtFilter_eta;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33EtFilter_phi;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33EtFilter_et;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33HEFilter_eta;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33HEFilter_phi;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33HEFilter_et;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33CaloIdLClusterShapeFilter_eta;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33CaloIdLClusterShapeFilter_phi;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33CaloIdLClusterShapeFilter_et;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_MW_hltEle33CaloIdLPixelMatchFilter_eta;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_MW_hltEle33CaloIdLPixelMatchFilter_phi;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_MW_hltEle33CaloIdLPixelMatchFilter_et;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_MW_hltEle33CaloIdLMWPMS2Filter_eta;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_MW_hltEle33CaloIdLMWPMS2Filter_phi;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_MW_hltEle33CaloIdLMWPMS2Filter_et;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_MW_hltEgammaCandidatesWrapperUnseeded_eta;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_MW_hltEgammaCandidatesWrapperUnseeded_phi;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_MW_hltEgammaCandidatesWrapperUnseeded_et;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33EtUnseededFilter_eta;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33EtUnseededFilter_phi;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33EtUnseededFilter_et;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33HEUnseededFilter_eta;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33HEUnseededFilter_phi;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33HEUnseededFilter_et;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33CaloIdLClusterShapeUnseededFilter_eta;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33CaloIdLClusterShapeUnseededFilter_phi;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33CaloIdLClusterShapeUnseededFilter_et;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEle33CaloIdLPixelMatchUnseededFilter_eta;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEle33CaloIdLPixelMatchUnseededFilter_phi;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEle33CaloIdLPixelMatchUnseededFilter_et;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEle33CaloIdLMWPMS2UnseededFilter_eta;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEle33CaloIdLMWPMS2UnseededFilter_phi;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEle33CaloIdLMWPMS2UnseededFilter_et;
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
   vector<float>   *trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltL1sSingleAndDoubleEGor_eta;
   vector<float>   *trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltL1sSingleAndDoubleEGor_phi;
   vector<float>   *trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltL1sSingleAndDoubleEGor_et;
   vector<float>   *trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEGL1SingleAndDoubleEGOrSingleFilter_eta;
   vector<float>   *trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEGL1SingleAndDoubleEGOrSingleFilter_phi;
   vector<float>   *trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEGL1SingleAndDoubleEGOrSingleFilter_et;
   vector<float>   *trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEG27L1SingleAndDoubleEGEtFilter_eta;
   vector<float>   *trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEG27L1SingleAndDoubleEGEtFilter_phi;
   vector<float>   *trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEG27L1SingleAndDoubleEGEtFilter_et;
   vector<float>   *trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightClusterShapeFilter_eta;
   vector<float>   *trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightClusterShapeFilter_phi;
   vector<float>   *trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightClusterShapeFilter_et;
   vector<float>   *trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightHEFilter_eta;
   vector<float>   *trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightHEFilter_phi;
   vector<float>   *trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightHEFilter_et;
   vector<float>   *trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightEcalIsoFilter_eta;
   vector<float>   *trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightEcalIsoFilter_phi;
   vector<float>   *trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightEcalIsoFilter_et;
   vector<float>   *trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightHcalIsoFilter_eta;
   vector<float>   *trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightHcalIsoFilter_phi;
   vector<float>   *trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightHcalIsoFilter_et;
   vector<float>   *trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEG27L1SingleAndDoubleEGEtFilter_eta;
   vector<float>   *trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEG27L1SingleAndDoubleEGEtFilter_phi;
   vector<float>   *trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEG27L1SingleAndDoubleEGEtFilter_et;
   vector<float>   *trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightClusterShapeFilter_eta;
   vector<float>   *trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightClusterShapeFilter_phi;
   vector<float>   *trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightClusterShapeFilter_et;
   vector<float>   *trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightHEFilter_eta;
   vector<float>   *trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightHEFilter_phi;
   vector<float>   *trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightHEFilter_et;
   vector<float>   *trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightEcalIsoFilter_eta;
   vector<float>   *trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightEcalIsoFilter_phi;
   vector<float>   *trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightEcalIsoFilter_et;
   vector<float>   *trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightHcalIsoFilter_eta;
   vector<float>   *trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightHcalIsoFilter_phi;
   vector<float>   *trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightHcalIsoFilter_et;
   Int_t           trig_HLT_Ele27_WPTight_Gsf_accept;
   Int_t           trig_HLT_Ele32_WPTight_Gsf_accept;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_hltL1sSingleEGor_eta;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_hltL1sSingleEGor_phi;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_hltL1sSingleEGor_et;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_hltEGL1SingleEGOrFilter_eta;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_hltEGL1SingleEGOrFilter_phi;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_hltEGL1SingleEGOrFilter_et;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_hltEG32L1SingleEGOrEtFilter_eta;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_hltEG32L1SingleEGOrEtFilter_phi;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_hltEG32L1SingleEGOrEtFilter_et;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightClusterShapeFilter_eta;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightClusterShapeFilter_phi;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightClusterShapeFilter_et;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightHEFilter_eta;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightHEFilter_phi;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightHEFilter_et;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightEcalIsoFilter_eta;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightEcalIsoFilter_phi;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightEcalIsoFilter_et;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightHcalIsoFilter_eta;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightHcalIsoFilter_phi;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightHcalIsoFilter_et;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightPixelMatchFilter_eta;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightPixelMatchFilter_phi;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightPixelMatchFilter_et;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightPMS2Filter_eta;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightPMS2Filter_phi;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightPMS2Filter_et;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfOneOEMinusOneOPFilter_eta;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfOneOEMinusOneOPFilter_phi;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfOneOEMinusOneOPFilter_et;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfMissingHitsFilter_eta;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfMissingHitsFilter_phi;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfMissingHitsFilter_et;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfDetaFilter_eta;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfDetaFilter_phi;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfDetaFilter_et;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfDphiFilter_eta;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfDphiFilter_phi;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfDphiFilter_et;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfTrackIsoFilter_eta;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfTrackIsoFilter_phi;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfTrackIsoFilter_et;
   Int_t           trig_HLT_Ele35_WPTight_Gsf_accept;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_hltL1sSingleEGor_eta;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_hltL1sSingleEGor_phi;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_hltL1sSingleEGor_et;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_hltEGL1SingleEGOrFilter_eta;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_hltEGL1SingleEGOrFilter_phi;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_hltEGL1SingleEGOrFilter_et;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_hltEG35L1SingleEGOrEtFilter_eta;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_hltEG35L1SingleEGOrEtFilter_phi;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_hltEG35L1SingleEGOrEtFilter_et;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightClusterShapeFilter_eta;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightClusterShapeFilter_phi;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightClusterShapeFilter_et;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightHEFilter_eta;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightHEFilter_phi;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightHEFilter_et;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightEcalIsoFilter_eta;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightEcalIsoFilter_phi;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightEcalIsoFilter_et;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightHcalIsoFilter_eta;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightHcalIsoFilter_phi;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightHcalIsoFilter_et;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightPixelMatchFilter_eta;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightPixelMatchFilter_phi;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightPixelMatchFilter_et;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightPMS2Filter_eta;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightPMS2Filter_phi;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightPMS2Filter_et;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfOneOEMinusOneOPFilter_eta;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfOneOEMinusOneOPFilter_phi;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfOneOEMinusOneOPFilter_et;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfMissingHitsFilter_eta;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfMissingHitsFilter_phi;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfMissingHitsFilter_et;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfDetaFilter_eta;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfDetaFilter_phi;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfDetaFilter_et;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfDphiFilter_eta;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfDphiFilter_phi;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfDphiFilter_et;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfTrackIsoFilter_eta;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfTrackIsoFilter_phi;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfTrackIsoFilter_et;
   Int_t           trig_HLT_Ele35_WPTight_Gsf_L1EGMT_accept;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltL1sAlCaSingleEle_eta;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltL1sAlCaSingleEle_phi;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltL1sAlCaSingleEle_et;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTHEFilter_eta;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTHEFilter_phi;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTHEFilter_et;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTEcalIsoFilter_eta;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTEcalIsoFilter_phi;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTEcalIsoFilter_et;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTHcalIsoFilter_eta;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTHcalIsoFilter_phi;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTHcalIsoFilter_et;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTPixelMatchFilter_eta;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTPixelMatchFilter_phi;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTPixelMatchFilter_et;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTPMS2Filter_eta;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTPMS2Filter_phi;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTPMS2Filter_et;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTDetaFilter_eta;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTDetaFilter_phi;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTDetaFilter_et;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTDphiFilter_eta;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTDphiFilter_phi;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTDphiFilter_et;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTTrackIsoFilter_eta;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTTrackIsoFilter_phi;
   vector<float>   *trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTTrackIsoFilter_et;
   Int_t           trig_HLT_Ele38_WPTight_Gsf_accept;
   Int_t           trig_HLT_Ele40_WPTight_Gsf_accept;
   Int_t           trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_accept;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltL1sSingleAndDoubleEGor_eta;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltL1sSingleAndDoubleEGor_phi;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltL1sSingleAndDoubleEGor_et;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEGL1SingleAndDoubleEGOrSingleFilter_eta;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEGL1SingleAndDoubleEGOrSingleFilter_phi;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEGL1SingleAndDoubleEGOrSingleFilter_et;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEG32L1SingleAndDoubleEGEtFilter_eta;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEG32L1SingleAndDoubleEGEtFilter_phi;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEG32L1SingleAndDoubleEGEtFilter_et;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightClusterShapeFilter_eta;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightClusterShapeFilter_phi;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightClusterShapeFilter_et;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightHEFilter_eta;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightHEFilter_phi;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightHEFilter_et;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightEcalIsoFilter_eta;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightEcalIsoFilter_phi;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightEcalIsoFilter_et;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightHcalIsoFilter_eta;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightHcalIsoFilter_phi;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightHcalIsoFilter_et;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightPixelMatchFilter_eta;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightPixelMatchFilter_phi;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightPixelMatchFilter_et;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightPMS2Filter_eta;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightPMS2Filter_phi;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightPMS2Filter_et;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfOneOEMinusOneOPFilter_eta;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfOneOEMinusOneOPFilter_phi;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfOneOEMinusOneOPFilter_et;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfMissingHitsFilter_eta;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfMissingHitsFilter_phi;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfMissingHitsFilter_et;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfDetaFilter_eta;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfDetaFilter_phi;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfDetaFilter_et;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfDphiFilter_eta;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfDphiFilter_phi;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfDphiFilter_et;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfTrackIsoFilter_eta;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfTrackIsoFilter_phi;
   vector<float>   *trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfTrackIsoFilter_et;
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
   Int_t           trig_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_accept;
   Int_t           trig_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4_accept;
   Int_t           trig_HLT_DoubleL2Mu23NoVtx_2Cha_accept;
   Int_t           trig_HLT_DoubleL2Mu25NoVtx_2Cha_accept;
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
   Int_t           trig_HLT_Mu30_TkMu0_Onia_accept;
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
   vector<float>   *trig_HLT_Photon300_NoHE_hltL1sSingleEGNonIsoOrWithJetAndTau_eta;
   vector<float>   *trig_HLT_Photon300_NoHE_hltL1sSingleEGNonIsoOrWithJetAndTau_phi;
   vector<float>   *trig_HLT_Photon300_NoHE_hltL1sSingleEGNonIsoOrWithJetAndTau_et;
   vector<float>   *trig_HLT_Photon300_NoHE_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_eta;
   vector<float>   *trig_HLT_Photon300_NoHE_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_phi;
   vector<float>   *trig_HLT_Photon300_NoHE_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_et;
   vector<float>   *trig_HLT_Photon300_NoHE_hltEG300erEtFilter_eta;
   vector<float>   *trig_HLT_Photon300_NoHE_hltEG300erEtFilter_phi;
   vector<float>   *trig_HLT_Photon300_NoHE_hltEG300erEtFilter_et;
   Int_t           trig_HLT_Mu8_TrkIsoVVL_accept;
   Int_t           trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept;
   Int_t           trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept;
   Int_t           trig_HLT_Mu17_TrkIsoVVL_accept;
   Int_t           trig_HLT_Mu19_TrkIsoVVL_accept;
   Int_t           trig_HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL_accept;
   Int_t           trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltL1sSingleAndDoubleEG_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltL1sSingleAndDoubleEG_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltL1sSingleAndDoubleEG_et;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEGL1SingleAndDoubleEGOrPairFilter_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEGL1SingleAndDoubleEGOrPairFilter_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEGL1SingleAndDoubleEGOrPairFilter_et;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg1Filter_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg1Filter_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg1Filter_et;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg2Filter_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg2Filter_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg2Filter_et;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg1Filter_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg1Filter_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg1Filter_et;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg2Filter_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg2Filter_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg2Filter_et;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg1Filter_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg1Filter_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg1Filter_et;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg2Filter_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg2Filter_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg2Filter_et;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg1Filter_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg1Filter_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg1Filter_et;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg2Filter_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg2Filter_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg2Filter_et;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg1Filter_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg1Filter_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg1Filter_et;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg2Filter_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg2Filter_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg2Filter_et;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg1Filter_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg1Filter_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg1Filter_et;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg2Filter_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg2Filter_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg2Filter_et;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg1Filter_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg1Filter_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg1Filter_et;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg2Filter_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg2Filter_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg2Filter_et;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg1Filter_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg1Filter_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg1Filter_et;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg2Filter_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg2Filter_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg2Filter_et;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg1Filter_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg1Filter_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg1Filter_et;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg2Filter_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg2Filter_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg2Filter_et;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_et;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_et;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter_et;
   Int_t           trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_accept;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltL1sSingleAndDoubleEG_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltL1sSingleAndDoubleEG_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltL1sSingleAndDoubleEG_et;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEGL1SingleAndDoubleEGOrPairFilter_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEGL1SingleAndDoubleEGOrPairFilter_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEGL1SingleAndDoubleEGOrPairFilter_et;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg1Filter_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg1Filter_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg1Filter_et;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg2Filter_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg2Filter_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg2Filter_et;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg1Filter_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg1Filter_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg1Filter_et;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg2Filter_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg2Filter_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg2Filter_et;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg1Filter_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg1Filter_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg1Filter_et;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg2Filter_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg2Filter_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg2Filter_et;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg1Filter_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg1Filter_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg1Filter_et;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg2Filter_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg2Filter_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg2Filter_et;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg1Filter_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg1Filter_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg1Filter_et;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg2Filter_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg2Filter_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg2Filter_et;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg1Filter_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg1Filter_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg1Filter_et;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg2Filter_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg2Filter_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg2Filter_et;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg1Filter_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg1Filter_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg1Filter_et;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg2Filter_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg2Filter_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg2Filter_et;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg1Filter_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg1Filter_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg1Filter_et;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg2Filter_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg2Filter_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg2Filter_et;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg1Filter_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg1Filter_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg1Filter_et;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg2Filter_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg2Filter_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg2Filter_et;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_et;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_et;
   Int_t           trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept;
   Int_t           trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept;
   Int_t           trig_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept;
   Int_t           trig_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept;
   Int_t           trig_HLT_Photon20_accept;
   vector<float>   *trig_HLT_Photon20_hltL1sSingleEG15_eta;
   vector<float>   *trig_HLT_Photon20_hltL1sSingleEG15_phi;
   vector<float>   *trig_HLT_Photon20_hltL1sSingleEG15_et;
   vector<float>   *trig_HLT_Photon20_hltEGL1SingleEG15Filter_eta;
   vector<float>   *trig_HLT_Photon20_hltEGL1SingleEG15Filter_phi;
   vector<float>   *trig_HLT_Photon20_hltEGL1SingleEG15Filter_et;
   vector<float>   *trig_HLT_Photon20_hltEG20L1EG15EtFilter_eta;
   vector<float>   *trig_HLT_Photon20_hltEG20L1EG15EtFilter_phi;
   vector<float>   *trig_HLT_Photon20_hltEG20L1EG15EtFilter_et;
   vector<float>   *trig_HLT_Photon20_hltEG20L1EG10OR15HEFilter_eta;
   vector<float>   *trig_HLT_Photon20_hltEG20L1EG10OR15HEFilter_phi;
   vector<float>   *trig_HLT_Photon20_hltEG20L1EG10OR15HEFilter_et;
   Int_t           trig_HLT_Photon33_accept;
   vector<float>   *trig_HLT_Photon33_hltL1sSingleEG26_eta;
   vector<float>   *trig_HLT_Photon33_hltL1sSingleEG26_phi;
   vector<float>   *trig_HLT_Photon33_hltL1sSingleEG26_et;
   vector<float>   *trig_HLT_Photon33_hltEGL1SingleEG26Filter_eta;
   vector<float>   *trig_HLT_Photon33_hltEGL1SingleEG26Filter_phi;
   vector<float>   *trig_HLT_Photon33_hltEGL1SingleEG26Filter_et;
   vector<float>   *trig_HLT_Photon33_hltEG33L1EG26EtFilter_eta;
   vector<float>   *trig_HLT_Photon33_hltEG33L1EG26EtFilter_phi;
   vector<float>   *trig_HLT_Photon33_hltEG33L1EG26EtFilter_et;
   vector<float>   *trig_HLT_Photon33_hltEG33L1EG26HEFilter_eta;
   vector<float>   *trig_HLT_Photon33_hltEG33L1EG26HEFilter_phi;
   vector<float>   *trig_HLT_Photon33_hltEG33L1EG26HEFilter_et;
   Int_t           trig_HLT_Photon50_accept;
   vector<float>   *trig_HLT_Photon50_hltL1sSingleEG34to45_eta;
   vector<float>   *trig_HLT_Photon50_hltL1sSingleEG34to45_phi;
   vector<float>   *trig_HLT_Photon50_hltL1sSingleEG34to45_et;
   vector<float>   *trig_HLT_Photon50_hltEGL1SingleEG34to45Filter_eta;
   vector<float>   *trig_HLT_Photon50_hltEGL1SingleEG34to45Filter_phi;
   vector<float>   *trig_HLT_Photon50_hltEGL1SingleEG34to45Filter_et;
   vector<float>   *trig_HLT_Photon50_hltEG50EtFilter_eta;
   vector<float>   *trig_HLT_Photon50_hltEG50EtFilter_phi;
   vector<float>   *trig_HLT_Photon50_hltEG50EtFilter_et;
   vector<float>   *trig_HLT_Photon50_hltEG50HEFilter_eta;
   vector<float>   *trig_HLT_Photon50_hltEG50HEFilter_phi;
   vector<float>   *trig_HLT_Photon50_hltEG50HEFilter_et;
   Int_t           trig_HLT_Photon75_accept;
   vector<float>   *trig_HLT_Photon75_hltL1sSingleEG34to50_eta;
   vector<float>   *trig_HLT_Photon75_hltL1sSingleEG34to50_phi;
   vector<float>   *trig_HLT_Photon75_hltL1sSingleEG34to50_et;
   vector<float>   *trig_HLT_Photon75_hltEGL1SingleEG34ORL1SingleEG40Filter_eta;
   vector<float>   *trig_HLT_Photon75_hltEGL1SingleEG34ORL1SingleEG40Filter_phi;
   vector<float>   *trig_HLT_Photon75_hltEGL1SingleEG34ORL1SingleEG40Filter_et;
   vector<float>   *trig_HLT_Photon75_hltEG75EtFilter_eta;
   vector<float>   *trig_HLT_Photon75_hltEG75EtFilter_phi;
   vector<float>   *trig_HLT_Photon75_hltEG75EtFilter_et;
   vector<float>   *trig_HLT_Photon75_hltEG75HEFilter_eta;
   vector<float>   *trig_HLT_Photon75_hltEG75HEFilter_phi;
   vector<float>   *trig_HLT_Photon75_hltEG75HEFilter_et;
   Int_t           trig_HLT_Photon90_accept;
   vector<float>   *trig_HLT_Photon90_hltL1sSingleEG34to50_eta;
   vector<float>   *trig_HLT_Photon90_hltL1sSingleEG34to50_phi;
   vector<float>   *trig_HLT_Photon90_hltL1sSingleEG34to50_et;
   vector<float>   *trig_HLT_Photon90_hltEGL1SingleEG34ORL1SingleEG40Filter_eta;
   vector<float>   *trig_HLT_Photon90_hltEGL1SingleEG34ORL1SingleEG40Filter_phi;
   vector<float>   *trig_HLT_Photon90_hltEGL1SingleEG34ORL1SingleEG40Filter_et;
   vector<float>   *trig_HLT_Photon90_hltEG90EtFilter_eta;
   vector<float>   *trig_HLT_Photon90_hltEG90EtFilter_phi;
   vector<float>   *trig_HLT_Photon90_hltEG90EtFilter_et;
   vector<float>   *trig_HLT_Photon90_hltEG90HEFilter_eta;
   vector<float>   *trig_HLT_Photon90_hltEG90HEFilter_phi;
   vector<float>   *trig_HLT_Photon90_hltEG90HEFilter_et;
   Int_t           trig_HLT_Photon120_accept;
   vector<float>   *trig_HLT_Photon120_hltL1sSingleEG34to50_eta;
   vector<float>   *trig_HLT_Photon120_hltL1sSingleEG34to50_phi;
   vector<float>   *trig_HLT_Photon120_hltL1sSingleEG34to50_et;
   vector<float>   *trig_HLT_Photon120_hltEGL1SingleEG34ORL1SingleEG40Filter_eta;
   vector<float>   *trig_HLT_Photon120_hltEGL1SingleEG34ORL1SingleEG40Filter_phi;
   vector<float>   *trig_HLT_Photon120_hltEGL1SingleEG34ORL1SingleEG40Filter_et;
   vector<float>   *trig_HLT_Photon120_hltEG120EtFilter_eta;
   vector<float>   *trig_HLT_Photon120_hltEG120EtFilter_phi;
   vector<float>   *trig_HLT_Photon120_hltEG120EtFilter_et;
   vector<float>   *trig_HLT_Photon120_hltEG120HEFilter_eta;
   vector<float>   *trig_HLT_Photon120_hltEG120HEFilter_phi;
   vector<float>   *trig_HLT_Photon120_hltEG120HEFilter_et;
   Int_t           trig_HLT_Photon150_accept;
   vector<float>   *trig_HLT_Photon150_hltL1sSingleEGNonIsoOrWithJetAndTau_eta;
   vector<float>   *trig_HLT_Photon150_hltL1sSingleEGNonIsoOrWithJetAndTau_phi;
   vector<float>   *trig_HLT_Photon150_hltL1sSingleEGNonIsoOrWithJetAndTau_et;
   vector<float>   *trig_HLT_Photon150_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_eta;
   vector<float>   *trig_HLT_Photon150_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_phi;
   vector<float>   *trig_HLT_Photon150_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_et;
   vector<float>   *trig_HLT_Photon150_hltEG150EtFilter_eta;
   vector<float>   *trig_HLT_Photon150_hltEG150EtFilter_phi;
   vector<float>   *trig_HLT_Photon150_hltEG150EtFilter_et;
   vector<float>   *trig_HLT_Photon150_hltEG150HEFilter_eta;
   vector<float>   *trig_HLT_Photon150_hltEG150HEFilter_phi;
   vector<float>   *trig_HLT_Photon150_hltEG150HEFilter_et;
   Int_t           trig_HLT_Photon175_accept;
   vector<float>   *trig_HLT_Photon175_hltL1sSingleEGNonIsoOrWithJetAndTau_eta;
   vector<float>   *trig_HLT_Photon175_hltL1sSingleEGNonIsoOrWithJetAndTau_phi;
   vector<float>   *trig_HLT_Photon175_hltL1sSingleEGNonIsoOrWithJetAndTau_et;
   vector<float>   *trig_HLT_Photon175_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_eta;
   vector<float>   *trig_HLT_Photon175_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_phi;
   vector<float>   *trig_HLT_Photon175_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_et;
   vector<float>   *trig_HLT_Photon175_hltEG175EtFilter_eta;
   vector<float>   *trig_HLT_Photon175_hltEG175EtFilter_phi;
   vector<float>   *trig_HLT_Photon175_hltEG175EtFilter_et;
   vector<float>   *trig_HLT_Photon175_hltEG175HEFilter_eta;
   vector<float>   *trig_HLT_Photon175_hltEG175HEFilter_phi;
   vector<float>   *trig_HLT_Photon175_hltEG175HEFilter_et;
   Int_t           trig_HLT_Photon200_accept;
   vector<float>   *trig_HLT_Photon200_hltL1sSingleEGNonIsoOrWithJetAndTau_eta;
   vector<float>   *trig_HLT_Photon200_hltL1sSingleEGNonIsoOrWithJetAndTau_phi;
   vector<float>   *trig_HLT_Photon200_hltL1sSingleEGNonIsoOrWithJetAndTau_et;
   vector<float>   *trig_HLT_Photon200_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_eta;
   vector<float>   *trig_HLT_Photon200_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_phi;
   vector<float>   *trig_HLT_Photon200_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_et;
   vector<float>   *trig_HLT_Photon200_hltEG200EtFilter_eta;
   vector<float>   *trig_HLT_Photon200_hltEG200EtFilter_phi;
   vector<float>   *trig_HLT_Photon200_hltEG200EtFilter_et;
   vector<float>   *trig_HLT_Photon200_hltEG200HEFilter_eta;
   vector<float>   *trig_HLT_Photon200_hltEG200HEFilter_phi;
   vector<float>   *trig_HLT_Photon200_hltEG200HEFilter_et;
   Int_t           trig_HLT_Photon50_R9Id90_HE10_IsoM_accept;
   vector<float>   *trig_HLT_Photon50_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_eta;
   vector<float>   *trig_HLT_Photon50_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_phi;
   vector<float>   *trig_HLT_Photon50_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_et;
   vector<float>   *trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_eta;
   vector<float>   *trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_phi;
   vector<float>   *trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_et;
   vector<float>   *trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMEtFilter_eta;
   vector<float>   *trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMEtFilter_phi;
   vector<float>   *trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMEtFilter_et;
   vector<float>   *trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMHEFilter_eta;
   vector<float>   *trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMHEFilter_phi;
   vector<float>   *trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMHEFilter_et;
   vector<float>   *trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMR9Filter_eta;
   vector<float>   *trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMR9Filter_phi;
   vector<float>   *trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMR9Filter_et;
   vector<float>   *trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMEcalIsoFilter_eta;
   vector<float>   *trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMEcalIsoFilter_phi;
   vector<float>   *trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMEcalIsoFilter_et;
   vector<float>   *trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMHcalIsoFilter_eta;
   vector<float>   *trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMHcalIsoFilter_phi;
   vector<float>   *trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMHcalIsoFilter_et;
   vector<float>   *trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMTrackIsoFilter_eta;
   vector<float>   *trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMTrackIsoFilter_phi;
   vector<float>   *trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMTrackIsoFilter_et;
   Int_t           trig_HLT_Photon75_R9Id90_HE10_IsoM_accept;
   vector<float>   *trig_HLT_Photon75_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_eta;
   vector<float>   *trig_HLT_Photon75_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_phi;
   vector<float>   *trig_HLT_Photon75_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_et;
   vector<float>   *trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_eta;
   vector<float>   *trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_phi;
   vector<float>   *trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_et;
   vector<float>   *trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMEtFilter_eta;
   vector<float>   *trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMEtFilter_phi;
   vector<float>   *trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMEtFilter_et;
   vector<float>   *trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMHEFilter_eta;
   vector<float>   *trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMHEFilter_phi;
   vector<float>   *trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMHEFilter_et;
   vector<float>   *trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMR9Filter_eta;
   vector<float>   *trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMR9Filter_phi;
   vector<float>   *trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMR9Filter_et;
   vector<float>   *trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMEcalIsoFilter_eta;
   vector<float>   *trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMEcalIsoFilter_phi;
   vector<float>   *trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMEcalIsoFilter_et;
   vector<float>   *trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMHcalIsoFilter_eta;
   vector<float>   *trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMHcalIsoFilter_phi;
   vector<float>   *trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMHcalIsoFilter_et;
   vector<float>   *trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMTrackIsoFilter_eta;
   vector<float>   *trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMTrackIsoFilter_phi;
   vector<float>   *trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMTrackIsoFilter_et;
   Int_t           trig_HLT_Photon90_R9Id90_HE10_IsoM_accept;
   vector<float>   *trig_HLT_Photon90_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_eta;
   vector<float>   *trig_HLT_Photon90_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_phi;
   vector<float>   *trig_HLT_Photon90_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_et;
   vector<float>   *trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_eta;
   vector<float>   *trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_phi;
   vector<float>   *trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_et;
   vector<float>   *trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMEtFilter_eta;
   vector<float>   *trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMEtFilter_phi;
   vector<float>   *trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMEtFilter_et;
   vector<float>   *trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMHEFilter_eta;
   vector<float>   *trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMHEFilter_phi;
   vector<float>   *trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMHEFilter_et;
   vector<float>   *trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMR9Filter_eta;
   vector<float>   *trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMR9Filter_phi;
   vector<float>   *trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMR9Filter_et;
   vector<float>   *trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMEcalIsoFilter_eta;
   vector<float>   *trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMEcalIsoFilter_phi;
   vector<float>   *trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMEcalIsoFilter_et;
   vector<float>   *trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMHcalIsoFilter_eta;
   vector<float>   *trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMHcalIsoFilter_phi;
   vector<float>   *trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMHcalIsoFilter_et;
   vector<float>   *trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMTrackIsoFilter_eta;
   vector<float>   *trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMTrackIsoFilter_phi;
   vector<float>   *trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMTrackIsoFilter_et;
   Int_t           trig_HLT_Photon120_R9Id90_HE10_IsoM_accept;
   vector<float>   *trig_HLT_Photon120_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_eta;
   vector<float>   *trig_HLT_Photon120_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_phi;
   vector<float>   *trig_HLT_Photon120_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_et;
   vector<float>   *trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_eta;
   vector<float>   *trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_phi;
   vector<float>   *trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_et;
   vector<float>   *trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMEtFilter_eta;
   vector<float>   *trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMEtFilter_phi;
   vector<float>   *trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMEtFilter_et;
   vector<float>   *trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMHEFilter_eta;
   vector<float>   *trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMHEFilter_phi;
   vector<float>   *trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMHEFilter_et;
   vector<float>   *trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMR9Filter_eta;
   vector<float>   *trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMR9Filter_phi;
   vector<float>   *trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMR9Filter_et;
   vector<float>   *trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMEcalIsoFilter_eta;
   vector<float>   *trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMEcalIsoFilter_phi;
   vector<float>   *trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMEcalIsoFilter_et;
   vector<float>   *trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMHcalIsoFilter_eta;
   vector<float>   *trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMHcalIsoFilter_phi;
   vector<float>   *trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMHcalIsoFilter_et;
   vector<float>   *trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMTrackIsoFilter_eta;
   vector<float>   *trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMTrackIsoFilter_phi;
   vector<float>   *trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMTrackIsoFilter_et;
   Int_t           trig_HLT_Photon165_R9Id90_HE10_IsoM_accept;
   vector<float>   *trig_HLT_Photon165_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_eta;
   vector<float>   *trig_HLT_Photon165_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_phi;
   vector<float>   *trig_HLT_Photon165_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_et;
   vector<float>   *trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_eta;
   vector<float>   *trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_phi;
   vector<float>   *trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_et;
   vector<float>   *trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMEtFilter_eta;
   vector<float>   *trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMEtFilter_phi;
   vector<float>   *trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMEtFilter_et;
   vector<float>   *trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMHEFilter_eta;
   vector<float>   *trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMHEFilter_phi;
   vector<float>   *trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMHEFilter_et;
   vector<float>   *trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMR9Filter_eta;
   vector<float>   *trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMR9Filter_phi;
   vector<float>   *trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMR9Filter_et;
   vector<float>   *trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMEcalIsoFilter_eta;
   vector<float>   *trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMEcalIsoFilter_phi;
   vector<float>   *trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMEcalIsoFilter_et;
   vector<float>   *trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMHcalIsoFilter_eta;
   vector<float>   *trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMHcalIsoFilter_phi;
   vector<float>   *trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMHcalIsoFilter_et;
   vector<float>   *trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMTrackIsoFilter_eta;
   vector<float>   *trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMTrackIsoFilter_phi;
   vector<float>   *trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMTrackIsoFilter_et;
   Int_t           trig_HLT_DoubleMu43NoFiltersNoVtx_accept;
   Int_t           trig_HLT_DoubleMu48NoFiltersNoVtx_accept;
   Int_t           trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_accept;
   vector<float>   *trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL1sMu5EG23IorMu7EG23IorMu20EG17IorMu23EG10_eta;
   vector<float>   *trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL1sMu5EG23IorMu7EG23IorMu20EG17IorMu23EG10_phi;
   vector<float>   *trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL1sMu5EG23IorMu7EG23IorMu20EG17IorMu23EG10_et;
   vector<float>   *trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL1fL1sMu5EG20orMu20EG15L1Filtered5_eta;
   vector<float>   *trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL1fL1sMu5EG20orMu20EG15L1Filtered5_phi;
   vector<float>   *trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL1fL1sMu5EG20orMu20EG15L1Filtered5_et;
   vector<float>   *trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_-hltL2fL1sMu5EG20orMu20EG15L1f5L2NoVtxFiltered16_eta;
   vector<float>   *trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_-hltL2fL1sMu5EG20orMu20EG15L1f5L2NoVtxFiltered16_phi;
   vector<float>   *trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_-hltL2fL1sMu5EG20orMu20EG15L1f5L2NoVtxFiltered16_et;
   vector<float>   *trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL3fL1sMu5EG20orMu20EG15L1f5L2NVf16L3NoFiltersNoVtxFiltered43_eta;
   vector<float>   *trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL3fL1sMu5EG20orMu20EG15L1f5L2NVf16L3NoFiltersNoVtxFiltered43_phi;
   vector<float>   *trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL3fL1sMu5EG20orMu20EG15L1f5L2NVf16L3NoFiltersNoVtxFiltered43_et;
   vector<float>   *trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltEGL1Mu5EG20Filter_eta;
   vector<float>   *trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltEGL1Mu5EG20Filter_phi;
   vector<float>   *trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltEGL1Mu5EG20Filter_et;
   vector<float>   *trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLEtFilter_eta;
   vector<float>   *trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLEtFilter_phi;
   vector<float>   *trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLEtFilter_et;
   vector<float>   *trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLClusterShapeFilter_eta;
   vector<float>   *trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLClusterShapeFilter_phi;
   vector<float>   *trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLClusterShapeFilter_et;
   vector<float>   *trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLHEFilter_eta;
   vector<float>   *trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLHEFilter_phi;
   vector<float>   *trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLHEFilter_et;
   Int_t           trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_accept;
   vector<float>   *trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL1sMu5EG23IorMu7EG23IorMu20EG17IorMu23EG10_eta;
   vector<float>   *trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL1sMu5EG23IorMu7EG23IorMu20EG17IorMu23EG10_phi;
   vector<float>   *trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL1sMu5EG23IorMu7EG23IorMu20EG17IorMu23EG10_et;
   vector<float>   *trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL1fL1sMu5EG20orMu20EG15L1Filtered5_eta;
   vector<float>   *trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL1fL1sMu5EG20orMu20EG15L1Filtered5_phi;
   vector<float>   *trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL1fL1sMu5EG20orMu20EG15L1Filtered5_et;
   vector<float>   *trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_-hltL2fL1sMu5EG20orMu20EG15L1f5L2NoVtxFiltered16_eta;
   vector<float>   *trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_-hltL2fL1sMu5EG20orMu20EG15L1f5L2NoVtxFiltered16_phi;
   vector<float>   *trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_-hltL2fL1sMu5EG20orMu20EG15L1f5L2NoVtxFiltered16_et;
   vector<float>   *trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL3fL1sMu5EG20orMu20EG15L1f5L2NVf16L3NoFiltersNoVtxFiltered48_eta;
   vector<float>   *trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL3fL1sMu5EG20orMu20EG15L1f5L2NVf16L3NoFiltersNoVtxFiltered48_phi;
   vector<float>   *trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL3fL1sMu5EG20orMu20EG15L1f5L2NVf16L3NoFiltersNoVtxFiltered48_et;
   vector<float>   *trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltEGL1Mu5EG20Filter_eta;
   vector<float>   *trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltEGL1Mu5EG20Filter_phi;
   vector<float>   *trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltEGL1Mu5EG20Filter_et;
   vector<float>   *trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLEtFilter_eta;
   vector<float>   *trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLEtFilter_phi;
   vector<float>   *trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLEtFilter_et;
   vector<float>   *trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLClusterShapeFilter_eta;
   vector<float>   *trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLClusterShapeFilter_phi;
   vector<float>   *trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLClusterShapeFilter_et;
   vector<float>   *trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLHEFilter_eta;
   vector<float>   *trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLHEFilter_phi;
   vector<float>   *trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLHEFilter_et;
   Int_t           trig_HLT_DoubleMu20_7_Mass0to30_L1_DM4_accept;
   Int_t           trig_HLT_DoubleMu20_7_Mass0to30_L1_DM4EG_accept;
   Int_t           trig_HLT_DoubleMu20_7_Mass0to30_Photon23_accept;
   vector<float>   *trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltL1sDoubleMu4OSEG12ORDoubleMu5OSEG12_eta;
   vector<float>   *trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltL1sDoubleMu4OSEG12ORDoubleMu5OSEG12_phi;
   vector<float>   *trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltL1sDoubleMu4OSEG12ORDoubleMu5OSEG12_et;
   vector<float>   *trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltL1fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1Filtered0_eta;
   vector<float>   *trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltL1fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1Filtered0_phi;
   vector<float>   *trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltL1fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1Filtered0_et;
   vector<float>   *trig_HLT_DoubleMu20_7_Mass0to30_Photon23_-hltL2fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1f0L2PreFiltered0_eta;
   vector<float>   *trig_HLT_DoubleMu20_7_Mass0to30_Photon23_-hltL2fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1f0L2PreFiltered0_phi;
   vector<float>   *trig_HLT_DoubleMu20_7_Mass0to30_Photon23_-hltL2fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1f0L2PreFiltered0_et;
   vector<float>   *trig_HLT_DoubleMu20_7_Mass0to30_Photon23_-hltL1fForIterL3L1fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1Filtered0_eta;
   vector<float>   *trig_HLT_DoubleMu20_7_Mass0to30_Photon23_-hltL1fForIterL3L1fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1Filtered0_phi;
   vector<float>   *trig_HLT_DoubleMu20_7_Mass0to30_Photon23_-hltL1fForIterL3L1fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1Filtered0_et;
   vector<float>   *trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltDoubleMu207Mass0to30Photon23L3Filtered_eta;
   vector<float>   *trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltDoubleMu207Mass0to30Photon23L3Filtered_phi;
   vector<float>   *trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltDoubleMu207Mass0to30Photon23L3Filtered_et;
   vector<float>   *trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEGDoubleMu4OSEG12ORDoubleMu5OSEG12Filter_eta;
   vector<float>   *trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEGDoubleMu4OSEG12ORDoubleMu5OSEG12Filter_phi;
   vector<float>   *trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEGDoubleMu4OSEG12ORDoubleMu5OSEG12Filter_et;
   vector<float>   *trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEG23EtFilter_eta;
   vector<float>   *trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEG23EtFilter_phi;
   vector<float>   *trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEG23EtFilter_et;
   vector<float>   *trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEG23HEFilter_eta;
   vector<float>   *trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEG23HEFilter_phi;
   vector<float>   *trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEG23HEFilter_et;
   Int_t           trig_HLT_Dimuon14_Phi_Barrel_Seagulls_accept;
   Int_t           trig_HLT_Dimuon24_Phi_noCorrL1_accept;
   Int_t           trig_HLT_DoubleIsoMu20_eta2p1_accept;
   Int_t           trig_HLT_DoubleIsoMu24_eta2p1_accept;
   Int_t           trig_HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx_accept;
   Int_t           trig_HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx_accept;
   Int_t           trig_HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx_accept;
   Int_t           trig_HLT_Mu8_accept;
   Int_t           trig_HLT_Mu17_accept;
   Int_t           trig_HLT_Mu19_accept;
   Int_t           trig_HLT_Mu17_Photon30_IsoCaloId_accept;
   vector<float>   *trig_HLT_Mu17_Photon30_IsoCaloId_hltL1sMu5EG23IorMu5IsoEG20IorMu7EG23IorMu7IsoEG20IorMuIso7EG23_eta;
   vector<float>   *trig_HLT_Mu17_Photon30_IsoCaloId_hltL1sMu5EG23IorMu5IsoEG20IorMu7EG23IorMu7IsoEG20IorMuIso7EG23_phi;
   vector<float>   *trig_HLT_Mu17_Photon30_IsoCaloId_hltL1sMu5EG23IorMu5IsoEG20IorMu7EG23IorMu7IsoEG20IorMuIso7EG23_et;
   vector<float>   *trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdMuonlegL1Filtered7_eta;
   vector<float>   *trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdMuonlegL1Filtered7_phi;
   vector<float>   *trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdMuonlegL1Filtered7_et;
   vector<float>   *trig_HLT_Mu17_Photon30_IsoCaloId_-hltMu17Photon30IsoCaloIdMuonlegL2Filtered10_eta;
   vector<float>   *trig_HLT_Mu17_Photon30_IsoCaloId_-hltMu17Photon30IsoCaloIdMuonlegL2Filtered10_phi;
   vector<float>   *trig_HLT_Mu17_Photon30_IsoCaloId_-hltMu17Photon30IsoCaloIdMuonlegL2Filtered10_et;
   vector<float>   *trig_HLT_Mu17_Photon30_IsoCaloId_-hltL1fForIterL3Mu17Photon30IsoCaloIdMuonlegL1Filtered7_eta;
   vector<float>   *trig_HLT_Mu17_Photon30_IsoCaloId_-hltL1fForIterL3Mu17Photon30IsoCaloIdMuonlegL1Filtered7_phi;
   vector<float>   *trig_HLT_Mu17_Photon30_IsoCaloId_-hltL1fForIterL3Mu17Photon30IsoCaloIdMuonlegL1Filtered7_et;
   vector<float>   *trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdMuonlegL3Filtered17Q_eta;
   vector<float>   *trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdMuonlegL3Filtered17Q_phi;
   vector<float>   *trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdMuonlegL3Filtered17Q_et;
   vector<float>   *trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegFilter_eta;
   vector<float>   *trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegFilter_phi;
   vector<float>   *trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegFilter_et;
   vector<float>   *trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegEtFilter_eta;
   vector<float>   *trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegEtFilter_phi;
   vector<float>   *trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegEtFilter_et;
   vector<float>   *trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegClusterShapeFilter_eta;
   vector<float>   *trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegClusterShapeFilter_phi;
   vector<float>   *trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegClusterShapeFilter_et;
   vector<float>   *trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegHEFilter_eta;
   vector<float>   *trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegHEFilter_phi;
   vector<float>   *trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegHEFilter_et;
   vector<float>   *trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegEcalIsoFilter_eta;
   vector<float>   *trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegEcalIsoFilter_phi;
   vector<float>   *trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegEcalIsoFilter_et;
   vector<float>   *trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegTrackIsoFilter_eta;
   vector<float>   *trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegTrackIsoFilter_phi;
   vector<float>   *trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegTrackIsoFilter_et;
   Int_t           trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_accept;
   vector<float>   *trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltL1sSingleEGNonIsoOrWithJetAndTau_eta;
   vector<float>   *trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltL1sSingleEGNonIsoOrWithJetAndTau_phi;
   vector<float>   *trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltL1sSingleEGNonIsoOrWithJetAndTau_et;
   vector<float>   *trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_eta;
   vector<float>   *trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_phi;
   vector<float>   *trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_et;
   vector<float>   *trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115EtFilter_eta;
   vector<float>   *trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115EtFilter_phi;
   vector<float>   *trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115EtFilter_et;
   vector<float>   *trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115CaloIdVTClusterShapeFilter_eta;
   vector<float>   *trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115CaloIdVTClusterShapeFilter_phi;
   vector<float>   *trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115CaloIdVTClusterShapeFilter_et;
   vector<float>   *trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115CaloIdVTHEFilter_eta;
   vector<float>   *trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115CaloIdVTHEFilter_phi;
   vector<float>   *trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115CaloIdVTHEFilter_et;
   vector<float>   *trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTPixelMatchFilter_eta;
   vector<float>   *trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTPixelMatchFilter_phi;
   vector<float>   *trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTPixelMatchFilter_et;
   vector<float>   *trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTGsfTrkIdTGsfDetaFilter_eta;
   vector<float>   *trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTGsfTrkIdTGsfDetaFilter_phi;
   vector<float>   *trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTGsfTrkIdTGsfDetaFilter_et;
   vector<float>   *trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTGsfTrkIdTGsfDphiFilter_eta;
   vector<float>   *trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTGsfTrkIdTGsfDphiFilter_phi;
   vector<float>   *trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTGsfTrkIdTGsfDphiFilter_et;
   Int_t           trig_HLT_Ele135_CaloIdVT_GsfTrkIdT_accept;
   Int_t           trig_HLT_Ele145_CaloIdVT_GsfTrkIdT_accept;
   Int_t           trig_HLT_Ele200_CaloIdVT_GsfTrkIdT_accept;
   Int_t           trig_HLT_Ele250_CaloIdVT_GsfTrkIdT_accept;
   Int_t           trig_HLT_Ele300_CaloIdVT_GsfTrkIdT_accept;
   Int_t           trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_accept;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltL1sSingleEGNonIsoOrWithJetAndTauNoPS_eta;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltL1sSingleEGNonIsoOrWithJetAndTauNoPS_phi;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltL1sSingleEGNonIsoOrWithJetAndTauNoPS_et;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEGL1SingleEGNonIsoOrWithJetAndTauNoPSFilter_eta;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEGL1SingleEGNonIsoOrWithJetAndTauNoPSFilter_phi;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEGL1SingleEGNonIsoOrWithJetAndTauNoPSFilter_et;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60EtFilter_eta;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60EtFilter_phi;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60EtFilter_et;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60HEFilter_eta;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60HEFilter_phi;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60HEFilter_et;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLR9IdFilter_eta;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLR9IdFilter_phi;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLR9IdFilter_et;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLClusterShapeFilter_eta;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLClusterShapeFilter_phi;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLClusterShapeFilter_et;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLEcalPFClusterIsoFilter_eta;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLEcalPFClusterIsoFilter_phi;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLEcalPFClusterIsoFilter_et;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLHcalPFClusterIsoFilter_eta;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLHcalPFClusterIsoFilter_phi;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLHcalPFClusterIsoFilter_et;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLHollowTrackIsoFilter_eta;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLHollowTrackIsoFilter_phi;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLHollowTrackIsoFilter_et;
   Int_t           trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_accept;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltL1sSingleEGNonIsoOrWithJetAndTauNoPS_eta;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltL1sSingleEGNonIsoOrWithJetAndTauNoPS_phi;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltL1sSingleEGNonIsoOrWithJetAndTauNoPS_et;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEGL1SingleEGNonIsoOrWithJetAndTauNoPSFilter_eta;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEGL1SingleEGNonIsoOrWithJetAndTauNoPSFilter_phi;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEGL1SingleEGNonIsoOrWithJetAndTauNoPSFilter_et;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60EtFilter_eta;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60EtFilter_phi;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60EtFilter_et;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60HEFilter_eta;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60HEFilter_phi;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60HEFilter_et;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLR9IdFilter_eta;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLR9IdFilter_phi;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLR9IdFilter_et;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLClusterShapeFilter_eta;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLClusterShapeFilter_phi;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLClusterShapeFilter_et;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLEcalPFClusterIsoFilter_eta;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLEcalPFClusterIsoFilter_phi;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLEcalPFClusterIsoFilter_et;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLHcalPFClusterIsoFilter_eta;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLHcalPFClusterIsoFilter_phi;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLHcalPFClusterIsoFilter_et;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLHollowTrackIsoFilter_eta;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLHollowTrackIsoFilter_phi;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLHollowTrackIsoFilter_et;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLDisplacedIdFilter_eta;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLDisplacedIdFilter_phi;
   vector<float>   *trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLDisplacedIdFilter_et;
   Int_t           trig_HLT_SinglePhoton10_Eta3p1ForPPRef_accept;
   vector<float>   *trig_HLT_SinglePhoton10_Eta3p1ForPPRef_hltL1sZeroBias_eta;
   vector<float>   *trig_HLT_SinglePhoton10_Eta3p1ForPPRef_hltL1sZeroBias_phi;
   vector<float>   *trig_HLT_SinglePhoton10_Eta3p1ForPPRef_hltL1sZeroBias_et;
   vector<float>   *trig_HLT_SinglePhoton10_Eta3p1ForPPRef_hltHIPhoton10Eta3p1_eta;
   vector<float>   *trig_HLT_SinglePhoton10_Eta3p1ForPPRef_hltHIPhoton10Eta3p1_phi;
   vector<float>   *trig_HLT_SinglePhoton10_Eta3p1ForPPRef_hltHIPhoton10Eta3p1_et;
   Int_t           trig_HLT_SinglePhoton20_Eta3p1ForPPRef_accept;
   vector<float>   *trig_HLT_SinglePhoton20_Eta3p1ForPPRef_hltL1sSingleEG10_eta;
   vector<float>   *trig_HLT_SinglePhoton20_Eta3p1ForPPRef_hltL1sSingleEG10_phi;
   vector<float>   *trig_HLT_SinglePhoton20_Eta3p1ForPPRef_hltL1sSingleEG10_et;
   vector<float>   *trig_HLT_SinglePhoton20_Eta3p1ForPPRef_hltHIPhoton20Eta3p1_eta;
   vector<float>   *trig_HLT_SinglePhoton20_Eta3p1ForPPRef_hltHIPhoton20Eta3p1_phi;
   vector<float>   *trig_HLT_SinglePhoton20_Eta3p1ForPPRef_hltHIPhoton20Eta3p1_et;
   Int_t           trig_HLT_SinglePhoton30_Eta3p1ForPPRef_accept;
   vector<float>   *trig_HLT_SinglePhoton30_Eta3p1ForPPRef_hltL1sSingleEG10_eta;
   vector<float>   *trig_HLT_SinglePhoton30_Eta3p1ForPPRef_hltL1sSingleEG10_phi;
   vector<float>   *trig_HLT_SinglePhoton30_Eta3p1ForPPRef_hltL1sSingleEG10_et;
   vector<float>   *trig_HLT_SinglePhoton30_Eta3p1ForPPRef_hltHIPhoton30Eta3p1_eta;
   vector<float>   *trig_HLT_SinglePhoton30_Eta3p1ForPPRef_hltHIPhoton30Eta3p1_phi;
   vector<float>   *trig_HLT_SinglePhoton30_Eta3p1ForPPRef_hltHIPhoton30Eta3p1_et;
   Int_t           trig_HLT_Photon20_HoverELoose_accept;
   vector<float>   *trig_HLT_Photon20_HoverELoose_hltL1sSingleEG10_eta;
   vector<float>   *trig_HLT_Photon20_HoverELoose_hltL1sSingleEG10_phi;
   vector<float>   *trig_HLT_Photon20_HoverELoose_hltL1sSingleEG10_et;
   vector<float>   *trig_HLT_Photon20_HoverELoose_hltEGL1SingleEG10Filter_eta;
   vector<float>   *trig_HLT_Photon20_HoverELoose_hltEGL1SingleEG10Filter_phi;
   vector<float>   *trig_HLT_Photon20_HoverELoose_hltEGL1SingleEG10Filter_et;
   vector<float>   *trig_HLT_Photon20_HoverELoose_hltEG20EtFilterLooseHoverE_eta;
   vector<float>   *trig_HLT_Photon20_HoverELoose_hltEG20EtFilterLooseHoverE_phi;
   vector<float>   *trig_HLT_Photon20_HoverELoose_hltEG20EtFilterLooseHoverE_et;
   vector<float>   *trig_HLT_Photon20_HoverELoose_hltEG20HEFilterLooseHoverE_eta;
   vector<float>   *trig_HLT_Photon20_HoverELoose_hltEG20HEFilterLooseHoverE_phi;
   vector<float>   *trig_HLT_Photon20_HoverELoose_hltEG20HEFilterLooseHoverE_et;
   Int_t           trig_HLT_Photon30_HoverELoose_accept;
   vector<float>   *trig_HLT_Photon30_HoverELoose_hltL1sSingleEG10_eta;
   vector<float>   *trig_HLT_Photon30_HoverELoose_hltL1sSingleEG10_phi;
   vector<float>   *trig_HLT_Photon30_HoverELoose_hltL1sSingleEG10_et;
   vector<float>   *trig_HLT_Photon30_HoverELoose_hltEGL1SingleEG10Filter_eta;
   vector<float>   *trig_HLT_Photon30_HoverELoose_hltEGL1SingleEG10Filter_phi;
   vector<float>   *trig_HLT_Photon30_HoverELoose_hltEGL1SingleEG10Filter_et;
   vector<float>   *trig_HLT_Photon30_HoverELoose_hltEG30EtFilterLooseHoverE_eta;
   vector<float>   *trig_HLT_Photon30_HoverELoose_hltEG30EtFilterLooseHoverE_phi;
   vector<float>   *trig_HLT_Photon30_HoverELoose_hltEG30EtFilterLooseHoverE_et;
   vector<float>   *trig_HLT_Photon30_HoverELoose_hltEG30HEFilterLooseHoverE_eta;
   vector<float>   *trig_HLT_Photon30_HoverELoose_hltEG30HEFilterLooseHoverE_phi;
   vector<float>   *trig_HLT_Photon30_HoverELoose_hltEG30HEFilterLooseHoverE_et;
   Int_t           trig_HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142_accept;
   Int_t           trig_AlCa_RPCMuonNormalisation_accept;
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
   TBranch        *b_trig_raw2digi_step_accept;   //!
   TBranch        *b_trig_L1Reco_step_accept;   //!
   TBranch        *b_trig_reconstruction_step_accept;   //!
   TBranch        *b_trig_EXOMONOPOLEPath_accept;   //!
   TBranch        *b_trig_ZElectronPathGsf_accept;   //!
   TBranch        *b_trig_pathALCARECOEcalESAlign_accept;   //!
   TBranch        *b_trig_pathALCARECOEcalUncalWElectron_accept;   //!
   TBranch        *b_trig_pathALCARECOHcalCalIsoTrkFilter_accept;   //!
   TBranch        *b_trig_pathALCARECOHcalCalIterativePhiSym_accept;   //!
   TBranch        *b_trig_pathALCARECOEcalUncalZElectron_accept;   //!
   TBranch        *b_trig_pathALCARECOEcalUncalZSCElectron_accept;   //!
   TBranch        *b_trig_eventinterpretaion_step_accept;   //!
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
   TBranch        *b_trig_HLT_DoubleEle25_CaloIdL_MW_accept;   //!
   TBranch        *b_trig_HLT_DoubleEle25_CaloIdL_MW_hltL1sSingleAndDoubleEGNonIsoOrWithEG26WithJetAndTau_eta;   //!
   TBranch        *b_trig_HLT_DoubleEle25_CaloIdL_MW_hltL1sSingleAndDoubleEGNonIsoOrWithEG26WithJetAndTau_phi;   //!
   TBranch        *b_trig_HLT_DoubleEle25_CaloIdL_MW_hltL1sSingleAndDoubleEGNonIsoOrWithEG26WithJetAndTau_et;   //!
   TBranch        *b_trig_HLT_DoubleEle25_CaloIdL_MW_hltEGL1SingleAndDoubleEGNonIsoOrWithEG26WithJetAndTauFilter_eta;   //!
   TBranch        *b_trig_HLT_DoubleEle25_CaloIdL_MW_hltEGL1SingleAndDoubleEGNonIsoOrWithEG26WithJetAndTauFilter_phi;   //!
   TBranch        *b_trig_HLT_DoubleEle25_CaloIdL_MW_hltEGL1SingleAndDoubleEGNonIsoOrWithEG26WithJetAndTauFilter_et;   //!
   TBranch        *b_trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25EtFilter_eta;   //!
   TBranch        *b_trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25EtFilter_phi;   //!
   TBranch        *b_trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25EtFilter_et;   //!
   TBranch        *b_trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25HEFilter_eta;   //!
   TBranch        *b_trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25HEFilter_phi;   //!
   TBranch        *b_trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25HEFilter_et;   //!
   TBranch        *b_trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25CaloIdLClusterShapeFilter_eta;   //!
   TBranch        *b_trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25CaloIdLClusterShapeFilter_phi;   //!
   TBranch        *b_trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25CaloIdLClusterShapeFilter_et;   //!
   TBranch        *b_trig_HLT_DoubleEle25_CaloIdL_MW_hltEle25CaloIdLPixelMatchFilter_eta;   //!
   TBranch        *b_trig_HLT_DoubleEle25_CaloIdL_MW_hltEle25CaloIdLPixelMatchFilter_phi;   //!
   TBranch        *b_trig_HLT_DoubleEle25_CaloIdL_MW_hltEle25CaloIdLPixelMatchFilter_et;   //!
   TBranch        *b_trig_HLT_DoubleEle25_CaloIdL_MW_hltEle25CaloIdLMWPMS2Filter_eta;   //!
   TBranch        *b_trig_HLT_DoubleEle25_CaloIdL_MW_hltEle25CaloIdLMWPMS2Filter_phi;   //!
   TBranch        *b_trig_HLT_DoubleEle25_CaloIdL_MW_hltEle25CaloIdLMWPMS2Filter_et;   //!
   TBranch        *b_trig_HLT_DoubleEle25_CaloIdL_MW_hltEgammaCandidatesWrapperUnseeded_eta;   //!
   TBranch        *b_trig_HLT_DoubleEle25_CaloIdL_MW_hltEgammaCandidatesWrapperUnseeded_phi;   //!
   TBranch        *b_trig_HLT_DoubleEle25_CaloIdL_MW_hltEgammaCandidatesWrapperUnseeded_et;   //!
   TBranch        *b_trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25EtUnseededFilter_eta;   //!
   TBranch        *b_trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25EtUnseededFilter_phi;   //!
   TBranch        *b_trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25EtUnseededFilter_et;   //!
   TBranch        *b_trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25HEUnseededFilter_eta;   //!
   TBranch        *b_trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25HEUnseededFilter_phi;   //!
   TBranch        *b_trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25HEUnseededFilter_et;   //!
   TBranch        *b_trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25CaloIdLClusterShapeUnseededFilter_eta;   //!
   TBranch        *b_trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25CaloIdLClusterShapeUnseededFilter_phi;   //!
   TBranch        *b_trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25CaloIdLClusterShapeUnseededFilter_et;   //!
   TBranch        *b_trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEle25CaloIdLPixelMatchUnseededFilter_eta;   //!
   TBranch        *b_trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEle25CaloIdLPixelMatchUnseededFilter_phi;   //!
   TBranch        *b_trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEle25CaloIdLPixelMatchUnseededFilter_et;   //!
   TBranch        *b_trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEle25CaloIdLMWPMS2UnseededFilter_eta;   //!
   TBranch        *b_trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEle25CaloIdLMWPMS2UnseededFilter_phi;   //!
   TBranch        *b_trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEle25CaloIdLMWPMS2UnseededFilter_et;   //!
   TBranch        *b_trig_HLT_DoubleEle27_CaloIdL_MW_accept;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_MW_accept;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_MW_hltL1sSingleAndDoubleEGNonIsoOrWithEG26WithJetAndTau_eta;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_MW_hltL1sSingleAndDoubleEGNonIsoOrWithEG26WithJetAndTau_phi;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_MW_hltL1sSingleAndDoubleEGNonIsoOrWithEG26WithJetAndTau_et;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_MW_hltEGL1SingleAndDoubleEGNonIsoOrWithEG26WithJetAndTauFilter_eta;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_MW_hltEGL1SingleAndDoubleEGNonIsoOrWithEG26WithJetAndTauFilter_phi;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_MW_hltEGL1SingleAndDoubleEGNonIsoOrWithEG26WithJetAndTauFilter_et;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33EtFilter_eta;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33EtFilter_phi;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33EtFilter_et;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33HEFilter_eta;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33HEFilter_phi;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33HEFilter_et;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33CaloIdLClusterShapeFilter_eta;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33CaloIdLClusterShapeFilter_phi;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33CaloIdLClusterShapeFilter_et;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_MW_hltEle33CaloIdLPixelMatchFilter_eta;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_MW_hltEle33CaloIdLPixelMatchFilter_phi;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_MW_hltEle33CaloIdLPixelMatchFilter_et;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_MW_hltEle33CaloIdLMWPMS2Filter_eta;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_MW_hltEle33CaloIdLMWPMS2Filter_phi;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_MW_hltEle33CaloIdLMWPMS2Filter_et;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_MW_hltEgammaCandidatesWrapperUnseeded_eta;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_MW_hltEgammaCandidatesWrapperUnseeded_phi;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_MW_hltEgammaCandidatesWrapperUnseeded_et;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33EtUnseededFilter_eta;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33EtUnseededFilter_phi;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33EtUnseededFilter_et;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33HEUnseededFilter_eta;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33HEUnseededFilter_phi;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33HEUnseededFilter_et;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33CaloIdLClusterShapeUnseededFilter_eta;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33CaloIdLClusterShapeUnseededFilter_phi;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33CaloIdLClusterShapeUnseededFilter_et;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEle33CaloIdLPixelMatchUnseededFilter_eta;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEle33CaloIdLPixelMatchUnseededFilter_phi;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEle33CaloIdLPixelMatchUnseededFilter_et;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEle33CaloIdLMWPMS2UnseededFilter_eta;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEle33CaloIdLMWPMS2UnseededFilter_phi;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEle33CaloIdLMWPMS2UnseededFilter_et;   //!
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
   TBranch        *b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltL1sSingleAndDoubleEGor_eta;   //!
   TBranch        *b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltL1sSingleAndDoubleEGor_phi;   //!
   TBranch        *b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltL1sSingleAndDoubleEGor_et;   //!
   TBranch        *b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEGL1SingleAndDoubleEGOrSingleFilter_eta;   //!
   TBranch        *b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEGL1SingleAndDoubleEGOrSingleFilter_phi;   //!
   TBranch        *b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEGL1SingleAndDoubleEGOrSingleFilter_et;   //!
   TBranch        *b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEG27L1SingleAndDoubleEGEtFilter_eta;   //!
   TBranch        *b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEG27L1SingleAndDoubleEGEtFilter_phi;   //!
   TBranch        *b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEG27L1SingleAndDoubleEGEtFilter_et;   //!
   TBranch        *b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightClusterShapeFilter_eta;   //!
   TBranch        *b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightClusterShapeFilter_phi;   //!
   TBranch        *b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightClusterShapeFilter_et;   //!
   TBranch        *b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightHEFilter_eta;   //!
   TBranch        *b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightHEFilter_phi;   //!
   TBranch        *b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightHEFilter_et;   //!
   TBranch        *b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightEcalIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightEcalIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightEcalIsoFilter_et;   //!
   TBranch        *b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightHcalIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightHcalIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightHcalIsoFilter_et;   //!
   TBranch        *b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEG27L1SingleAndDoubleEGEtFilter_eta;   //!
   TBranch        *b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEG27L1SingleAndDoubleEGEtFilter_phi;   //!
   TBranch        *b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEG27L1SingleAndDoubleEGEtFilter_et;   //!
   TBranch        *b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightClusterShapeFilter_eta;   //!
   TBranch        *b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightClusterShapeFilter_phi;   //!
   TBranch        *b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightClusterShapeFilter_et;   //!
   TBranch        *b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightHEFilter_eta;   //!
   TBranch        *b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightHEFilter_phi;   //!
   TBranch        *b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightHEFilter_et;   //!
   TBranch        *b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightEcalIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightEcalIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightEcalIsoFilter_et;   //!
   TBranch        *b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightHcalIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightHcalIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightHcalIsoFilter_et;   //!
   TBranch        *b_trig_HLT_Ele27_WPTight_Gsf_accept;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_accept;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_hltL1sSingleEGor_eta;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_hltL1sSingleEGor_phi;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_hltL1sSingleEGor_et;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_hltEGL1SingleEGOrFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_hltEGL1SingleEGOrFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_hltEGL1SingleEGOrFilter_et;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_hltEG32L1SingleEGOrEtFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_hltEG32L1SingleEGOrEtFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_hltEG32L1SingleEGOrEtFilter_et;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightClusterShapeFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightClusterShapeFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightClusterShapeFilter_et;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightHEFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightHEFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightHEFilter_et;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightEcalIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightEcalIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightEcalIsoFilter_et;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightHcalIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightHcalIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightHcalIsoFilter_et;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightPixelMatchFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightPixelMatchFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightPixelMatchFilter_et;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightPMS2Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightPMS2Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightPMS2Filter_et;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfOneOEMinusOneOPFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfOneOEMinusOneOPFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfOneOEMinusOneOPFilter_et;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfMissingHitsFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfMissingHitsFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfMissingHitsFilter_et;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfDetaFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfDetaFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfDetaFilter_et;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfDphiFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfDphiFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfDphiFilter_et;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfTrackIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfTrackIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfTrackIsoFilter_et;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_accept;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_hltL1sSingleEGor_eta;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_hltL1sSingleEGor_phi;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_hltL1sSingleEGor_et;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_hltEGL1SingleEGOrFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_hltEGL1SingleEGOrFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_hltEGL1SingleEGOrFilter_et;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_hltEG35L1SingleEGOrEtFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_hltEG35L1SingleEGOrEtFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_hltEG35L1SingleEGOrEtFilter_et;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightClusterShapeFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightClusterShapeFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightClusterShapeFilter_et;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightHEFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightHEFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightHEFilter_et;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightEcalIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightEcalIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightEcalIsoFilter_et;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightHcalIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightHcalIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightHcalIsoFilter_et;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightPixelMatchFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightPixelMatchFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightPixelMatchFilter_et;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightPMS2Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightPMS2Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightPMS2Filter_et;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfOneOEMinusOneOPFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfOneOEMinusOneOPFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfOneOEMinusOneOPFilter_et;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfMissingHitsFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfMissingHitsFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfMissingHitsFilter_et;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfDetaFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfDetaFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfDetaFilter_et;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfDphiFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfDphiFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfDphiFilter_et;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfTrackIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfTrackIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfTrackIsoFilter_et;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_accept;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltL1sAlCaSingleEle_eta;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltL1sAlCaSingleEle_phi;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltL1sAlCaSingleEle_et;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTHEFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTHEFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTHEFilter_et;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTEcalIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTEcalIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTEcalIsoFilter_et;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTHcalIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTHcalIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTHcalIsoFilter_et;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTPixelMatchFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTPixelMatchFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTPixelMatchFilter_et;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTPMS2Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTPMS2Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTPMS2Filter_et;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTDetaFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTDetaFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTDetaFilter_et;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTDphiFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTDphiFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTDphiFilter_et;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTTrackIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTTrackIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTTrackIsoFilter_et;   //!
   TBranch        *b_trig_HLT_Ele38_WPTight_Gsf_accept;   //!
   TBranch        *b_trig_HLT_Ele40_WPTight_Gsf_accept;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_accept;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltL1sSingleAndDoubleEGor_eta;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltL1sSingleAndDoubleEGor_phi;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltL1sSingleAndDoubleEGor_et;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEGL1SingleAndDoubleEGOrSingleFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEGL1SingleAndDoubleEGOrSingleFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEGL1SingleAndDoubleEGOrSingleFilter_et;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEG32L1SingleAndDoubleEGEtFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEG32L1SingleAndDoubleEGEtFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEG32L1SingleAndDoubleEGEtFilter_et;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightClusterShapeFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightClusterShapeFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightClusterShapeFilter_et;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightHEFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightHEFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightHEFilter_et;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightEcalIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightEcalIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightEcalIsoFilter_et;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightHcalIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightHcalIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightHcalIsoFilter_et;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightPixelMatchFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightPixelMatchFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightPixelMatchFilter_et;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightPMS2Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightPMS2Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightPMS2Filter_et;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfOneOEMinusOneOPFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfOneOEMinusOneOPFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfOneOEMinusOneOPFilter_et;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfMissingHitsFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfMissingHitsFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfMissingHitsFilter_et;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfDetaFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfDetaFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfDetaFilter_et;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfDphiFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfDphiFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfDphiFilter_et;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfTrackIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfTrackIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfTrackIsoFilter_et;   //!
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
   TBranch        *b_trig_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_accept;   //!
   TBranch        *b_trig_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4_accept;   //!
   TBranch        *b_trig_HLT_DoubleL2Mu23NoVtx_2Cha_accept;   //!
   TBranch        *b_trig_HLT_DoubleL2Mu25NoVtx_2Cha_accept;   //!
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
   TBranch        *b_trig_HLT_Mu30_TkMu0_Onia_accept;   //!
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
   TBranch        *b_trig_HLT_Photon300_NoHE_hltL1sSingleEGNonIsoOrWithJetAndTau_eta;   //!
   TBranch        *b_trig_HLT_Photon300_NoHE_hltL1sSingleEGNonIsoOrWithJetAndTau_phi;   //!
   TBranch        *b_trig_HLT_Photon300_NoHE_hltL1sSingleEGNonIsoOrWithJetAndTau_et;   //!
   TBranch        *b_trig_HLT_Photon300_NoHE_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon300_NoHE_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon300_NoHE_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_et;   //!
   TBranch        *b_trig_HLT_Photon300_NoHE_hltEG300erEtFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon300_NoHE_hltEG300erEtFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon300_NoHE_hltEG300erEtFilter_et;   //!
   TBranch        *b_trig_HLT_Mu8_TrkIsoVVL_accept;   //!
   TBranch        *b_trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept;   //!
   TBranch        *b_trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept;   //!
   TBranch        *b_trig_HLT_Mu17_TrkIsoVVL_accept;   //!
   TBranch        *b_trig_HLT_Mu19_TrkIsoVVL_accept;   //!
   TBranch        *b_trig_HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL_accept;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltL1sSingleAndDoubleEG_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltL1sSingleAndDoubleEG_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltL1sSingleAndDoubleEG_et;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEGL1SingleAndDoubleEGOrPairFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEGL1SingleAndDoubleEGOrPairFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEGL1SingleAndDoubleEGOrPairFilter_et;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg1Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg1Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg1Filter_et;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg2Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg2Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg2Filter_et;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg1Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg1Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg1Filter_et;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg2Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg2Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg2Filter_et;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg1Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg1Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg1Filter_et;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg2Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg2Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg2Filter_et;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg1Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg1Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg1Filter_et;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg2Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg2Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg2Filter_et;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg1Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg1Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg1Filter_et;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg2Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg2Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg2Filter_et;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg1Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg1Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg1Filter_et;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg2Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg2Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg2Filter_et;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg1Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg1Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg1Filter_et;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg2Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg2Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg2Filter_et;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg1Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg1Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg1Filter_et;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg2Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg2Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg2Filter_et;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg1Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg1Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg1Filter_et;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg2Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg2Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg2Filter_et;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_et;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_et;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter_et;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_accept;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltL1sSingleAndDoubleEG_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltL1sSingleAndDoubleEG_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltL1sSingleAndDoubleEG_et;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEGL1SingleAndDoubleEGOrPairFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEGL1SingleAndDoubleEGOrPairFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEGL1SingleAndDoubleEGOrPairFilter_et;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg1Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg1Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg1Filter_et;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg2Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg2Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg2Filter_et;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg1Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg1Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg1Filter_et;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg2Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg2Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg2Filter_et;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg1Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg1Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg1Filter_et;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg2Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg2Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg2Filter_et;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg1Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg1Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg1Filter_et;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg2Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg2Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg2Filter_et;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg1Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg1Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg1Filter_et;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg2Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg2Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg2Filter_et;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg1Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg1Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg1Filter_et;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg2Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg2Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg2Filter_et;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg1Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg1Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg1Filter_et;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg2Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg2Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg2Filter_et;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg1Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg1Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg1Filter_et;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg2Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg2Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg2Filter_et;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg1Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg1Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg1Filter_et;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg2Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg2Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg2Filter_et;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_et;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_et;   //!
   TBranch        *b_trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept;   //!
   TBranch        *b_trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept;   //!
   TBranch        *b_trig_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept;   //!
   TBranch        *b_trig_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept;   //!
   TBranch        *b_trig_HLT_Photon20_accept;   //!
   TBranch        *b_trig_HLT_Photon20_hltL1sSingleEG15_eta;   //!
   TBranch        *b_trig_HLT_Photon20_hltL1sSingleEG15_phi;   //!
   TBranch        *b_trig_HLT_Photon20_hltL1sSingleEG15_et;   //!
   TBranch        *b_trig_HLT_Photon20_hltEGL1SingleEG15Filter_eta;   //!
   TBranch        *b_trig_HLT_Photon20_hltEGL1SingleEG15Filter_phi;   //!
   TBranch        *b_trig_HLT_Photon20_hltEGL1SingleEG15Filter_et;   //!
   TBranch        *b_trig_HLT_Photon20_hltEG20L1EG15EtFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon20_hltEG20L1EG15EtFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon20_hltEG20L1EG15EtFilter_et;   //!
   TBranch        *b_trig_HLT_Photon20_hltEG20L1EG10OR15HEFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon20_hltEG20L1EG10OR15HEFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon20_hltEG20L1EG10OR15HEFilter_et;   //!
   TBranch        *b_trig_HLT_Photon33_accept;   //!
   TBranch        *b_trig_HLT_Photon33_hltL1sSingleEG26_eta;   //!
   TBranch        *b_trig_HLT_Photon33_hltL1sSingleEG26_phi;   //!
   TBranch        *b_trig_HLT_Photon33_hltL1sSingleEG26_et;   //!
   TBranch        *b_trig_HLT_Photon33_hltEGL1SingleEG26Filter_eta;   //!
   TBranch        *b_trig_HLT_Photon33_hltEGL1SingleEG26Filter_phi;   //!
   TBranch        *b_trig_HLT_Photon33_hltEGL1SingleEG26Filter_et;   //!
   TBranch        *b_trig_HLT_Photon33_hltEG33L1EG26EtFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon33_hltEG33L1EG26EtFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon33_hltEG33L1EG26EtFilter_et;   //!
   TBranch        *b_trig_HLT_Photon33_hltEG33L1EG26HEFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon33_hltEG33L1EG26HEFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon33_hltEG33L1EG26HEFilter_et;   //!
   TBranch        *b_trig_HLT_Photon50_accept;   //!
   TBranch        *b_trig_HLT_Photon50_hltL1sSingleEG34to45_eta;   //!
   TBranch        *b_trig_HLT_Photon50_hltL1sSingleEG34to45_phi;   //!
   TBranch        *b_trig_HLT_Photon50_hltL1sSingleEG34to45_et;   //!
   TBranch        *b_trig_HLT_Photon50_hltEGL1SingleEG34to45Filter_eta;   //!
   TBranch        *b_trig_HLT_Photon50_hltEGL1SingleEG34to45Filter_phi;   //!
   TBranch        *b_trig_HLT_Photon50_hltEGL1SingleEG34to45Filter_et;   //!
   TBranch        *b_trig_HLT_Photon50_hltEG50EtFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon50_hltEG50EtFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon50_hltEG50EtFilter_et;   //!
   TBranch        *b_trig_HLT_Photon50_hltEG50HEFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon50_hltEG50HEFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon50_hltEG50HEFilter_et;   //!
   TBranch        *b_trig_HLT_Photon75_accept;   //!
   TBranch        *b_trig_HLT_Photon75_hltL1sSingleEG34to50_eta;   //!
   TBranch        *b_trig_HLT_Photon75_hltL1sSingleEG34to50_phi;   //!
   TBranch        *b_trig_HLT_Photon75_hltL1sSingleEG34to50_et;   //!
   TBranch        *b_trig_HLT_Photon75_hltEGL1SingleEG34ORL1SingleEG40Filter_eta;   //!
   TBranch        *b_trig_HLT_Photon75_hltEGL1SingleEG34ORL1SingleEG40Filter_phi;   //!
   TBranch        *b_trig_HLT_Photon75_hltEGL1SingleEG34ORL1SingleEG40Filter_et;   //!
   TBranch        *b_trig_HLT_Photon75_hltEG75EtFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon75_hltEG75EtFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon75_hltEG75EtFilter_et;   //!
   TBranch        *b_trig_HLT_Photon75_hltEG75HEFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon75_hltEG75HEFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon75_hltEG75HEFilter_et;   //!
   TBranch        *b_trig_HLT_Photon90_accept;   //!
   TBranch        *b_trig_HLT_Photon90_hltL1sSingleEG34to50_eta;   //!
   TBranch        *b_trig_HLT_Photon90_hltL1sSingleEG34to50_phi;   //!
   TBranch        *b_trig_HLT_Photon90_hltL1sSingleEG34to50_et;   //!
   TBranch        *b_trig_HLT_Photon90_hltEGL1SingleEG34ORL1SingleEG40Filter_eta;   //!
   TBranch        *b_trig_HLT_Photon90_hltEGL1SingleEG34ORL1SingleEG40Filter_phi;   //!
   TBranch        *b_trig_HLT_Photon90_hltEGL1SingleEG34ORL1SingleEG40Filter_et;   //!
   TBranch        *b_trig_HLT_Photon90_hltEG90EtFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon90_hltEG90EtFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon90_hltEG90EtFilter_et;   //!
   TBranch        *b_trig_HLT_Photon90_hltEG90HEFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon90_hltEG90HEFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon90_hltEG90HEFilter_et;   //!
   TBranch        *b_trig_HLT_Photon120_accept;   //!
   TBranch        *b_trig_HLT_Photon120_hltL1sSingleEG34to50_eta;   //!
   TBranch        *b_trig_HLT_Photon120_hltL1sSingleEG34to50_phi;   //!
   TBranch        *b_trig_HLT_Photon120_hltL1sSingleEG34to50_et;   //!
   TBranch        *b_trig_HLT_Photon120_hltEGL1SingleEG34ORL1SingleEG40Filter_eta;   //!
   TBranch        *b_trig_HLT_Photon120_hltEGL1SingleEG34ORL1SingleEG40Filter_phi;   //!
   TBranch        *b_trig_HLT_Photon120_hltEGL1SingleEG34ORL1SingleEG40Filter_et;   //!
   TBranch        *b_trig_HLT_Photon120_hltEG120EtFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon120_hltEG120EtFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon120_hltEG120EtFilter_et;   //!
   TBranch        *b_trig_HLT_Photon120_hltEG120HEFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon120_hltEG120HEFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon120_hltEG120HEFilter_et;   //!
   TBranch        *b_trig_HLT_Photon150_accept;   //!
   TBranch        *b_trig_HLT_Photon150_hltL1sSingleEGNonIsoOrWithJetAndTau_eta;   //!
   TBranch        *b_trig_HLT_Photon150_hltL1sSingleEGNonIsoOrWithJetAndTau_phi;   //!
   TBranch        *b_trig_HLT_Photon150_hltL1sSingleEGNonIsoOrWithJetAndTau_et;   //!
   TBranch        *b_trig_HLT_Photon150_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon150_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon150_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_et;   //!
   TBranch        *b_trig_HLT_Photon150_hltEG150EtFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon150_hltEG150EtFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon150_hltEG150EtFilter_et;   //!
   TBranch        *b_trig_HLT_Photon150_hltEG150HEFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon150_hltEG150HEFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon150_hltEG150HEFilter_et;   //!
   TBranch        *b_trig_HLT_Photon175_accept;   //!
   TBranch        *b_trig_HLT_Photon175_hltL1sSingleEGNonIsoOrWithJetAndTau_eta;   //!
   TBranch        *b_trig_HLT_Photon175_hltL1sSingleEGNonIsoOrWithJetAndTau_phi;   //!
   TBranch        *b_trig_HLT_Photon175_hltL1sSingleEGNonIsoOrWithJetAndTau_et;   //!
   TBranch        *b_trig_HLT_Photon175_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon175_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon175_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_et;   //!
   TBranch        *b_trig_HLT_Photon175_hltEG175EtFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon175_hltEG175EtFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon175_hltEG175EtFilter_et;   //!
   TBranch        *b_trig_HLT_Photon175_hltEG175HEFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon175_hltEG175HEFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon175_hltEG175HEFilter_et;   //!
   TBranch        *b_trig_HLT_Photon200_accept;   //!
   TBranch        *b_trig_HLT_Photon200_hltL1sSingleEGNonIsoOrWithJetAndTau_eta;   //!
   TBranch        *b_trig_HLT_Photon200_hltL1sSingleEGNonIsoOrWithJetAndTau_phi;   //!
   TBranch        *b_trig_HLT_Photon200_hltL1sSingleEGNonIsoOrWithJetAndTau_et;   //!
   TBranch        *b_trig_HLT_Photon200_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon200_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon200_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_et;   //!
   TBranch        *b_trig_HLT_Photon200_hltEG200EtFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon200_hltEG200EtFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon200_hltEG200EtFilter_et;   //!
   TBranch        *b_trig_HLT_Photon200_hltEG200HEFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon200_hltEG200HEFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon200_hltEG200HEFilter_et;   //!
   TBranch        *b_trig_HLT_Photon50_R9Id90_HE10_IsoM_accept;   //!
   TBranch        *b_trig_HLT_Photon50_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_eta;   //!
   TBranch        *b_trig_HLT_Photon50_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_phi;   //!
   TBranch        *b_trig_HLT_Photon50_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_et;   //!
   TBranch        *b_trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_eta;   //!
   TBranch        *b_trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_phi;   //!
   TBranch        *b_trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_et;   //!
   TBranch        *b_trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMEtFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMEtFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMEtFilter_et;   //!
   TBranch        *b_trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMHEFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMHEFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMHEFilter_et;   //!
   TBranch        *b_trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMR9Filter_eta;   //!
   TBranch        *b_trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMR9Filter_phi;   //!
   TBranch        *b_trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMR9Filter_et;   //!
   TBranch        *b_trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMEcalIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMEcalIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMEcalIsoFilter_et;   //!
   TBranch        *b_trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMHcalIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMHcalIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMHcalIsoFilter_et;   //!
   TBranch        *b_trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMTrackIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMTrackIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMTrackIsoFilter_et;   //!
   TBranch        *b_trig_HLT_Photon75_R9Id90_HE10_IsoM_accept;   //!
   TBranch        *b_trig_HLT_Photon75_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_eta;   //!
   TBranch        *b_trig_HLT_Photon75_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_phi;   //!
   TBranch        *b_trig_HLT_Photon75_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_et;   //!
   TBranch        *b_trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_eta;   //!
   TBranch        *b_trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_phi;   //!
   TBranch        *b_trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_et;   //!
   TBranch        *b_trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMEtFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMEtFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMEtFilter_et;   //!
   TBranch        *b_trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMHEFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMHEFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMHEFilter_et;   //!
   TBranch        *b_trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMR9Filter_eta;   //!
   TBranch        *b_trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMR9Filter_phi;   //!
   TBranch        *b_trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMR9Filter_et;   //!
   TBranch        *b_trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMEcalIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMEcalIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMEcalIsoFilter_et;   //!
   TBranch        *b_trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMHcalIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMHcalIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMHcalIsoFilter_et;   //!
   TBranch        *b_trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMTrackIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMTrackIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMTrackIsoFilter_et;   //!
   TBranch        *b_trig_HLT_Photon90_R9Id90_HE10_IsoM_accept;   //!
   TBranch        *b_trig_HLT_Photon90_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_eta;   //!
   TBranch        *b_trig_HLT_Photon90_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_phi;   //!
   TBranch        *b_trig_HLT_Photon90_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_et;   //!
   TBranch        *b_trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_eta;   //!
   TBranch        *b_trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_phi;   //!
   TBranch        *b_trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_et;   //!
   TBranch        *b_trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMEtFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMEtFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMEtFilter_et;   //!
   TBranch        *b_trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMHEFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMHEFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMHEFilter_et;   //!
   TBranch        *b_trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMR9Filter_eta;   //!
   TBranch        *b_trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMR9Filter_phi;   //!
   TBranch        *b_trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMR9Filter_et;   //!
   TBranch        *b_trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMEcalIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMEcalIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMEcalIsoFilter_et;   //!
   TBranch        *b_trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMHcalIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMHcalIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMHcalIsoFilter_et;   //!
   TBranch        *b_trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMTrackIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMTrackIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMTrackIsoFilter_et;   //!
   TBranch        *b_trig_HLT_Photon120_R9Id90_HE10_IsoM_accept;   //!
   TBranch        *b_trig_HLT_Photon120_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_eta;   //!
   TBranch        *b_trig_HLT_Photon120_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_phi;   //!
   TBranch        *b_trig_HLT_Photon120_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_et;   //!
   TBranch        *b_trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_eta;   //!
   TBranch        *b_trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_phi;   //!
   TBranch        *b_trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_et;   //!
   TBranch        *b_trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMEtFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMEtFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMEtFilter_et;   //!
   TBranch        *b_trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMHEFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMHEFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMHEFilter_et;   //!
   TBranch        *b_trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMR9Filter_eta;   //!
   TBranch        *b_trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMR9Filter_phi;   //!
   TBranch        *b_trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMR9Filter_et;   //!
   TBranch        *b_trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMEcalIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMEcalIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMEcalIsoFilter_et;   //!
   TBranch        *b_trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMHcalIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMHcalIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMHcalIsoFilter_et;   //!
   TBranch        *b_trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMTrackIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMTrackIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMTrackIsoFilter_et;   //!
   TBranch        *b_trig_HLT_Photon165_R9Id90_HE10_IsoM_accept;   //!
   TBranch        *b_trig_HLT_Photon165_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_eta;   //!
   TBranch        *b_trig_HLT_Photon165_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_phi;   //!
   TBranch        *b_trig_HLT_Photon165_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_et;   //!
   TBranch        *b_trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_eta;   //!
   TBranch        *b_trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_phi;   //!
   TBranch        *b_trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_et;   //!
   TBranch        *b_trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMEtFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMEtFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMEtFilter_et;   //!
   TBranch        *b_trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMHEFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMHEFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMHEFilter_et;   //!
   TBranch        *b_trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMR9Filter_eta;   //!
   TBranch        *b_trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMR9Filter_phi;   //!
   TBranch        *b_trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMR9Filter_et;   //!
   TBranch        *b_trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMEcalIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMEcalIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMEcalIsoFilter_et;   //!
   TBranch        *b_trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMHcalIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMHcalIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMHcalIsoFilter_et;   //!
   TBranch        *b_trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMTrackIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMTrackIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMTrackIsoFilter_et;   //!
   TBranch        *b_trig_HLT_DoubleMu43NoFiltersNoVtx_accept;   //!
   TBranch        *b_trig_HLT_DoubleMu48NoFiltersNoVtx_accept;   //!
   TBranch        *b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_accept;   //!
   TBranch        *b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL1sMu5EG23IorMu7EG23IorMu20EG17IorMu23EG10_eta;   //!
   TBranch        *b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL1sMu5EG23IorMu7EG23IorMu20EG17IorMu23EG10_phi;   //!
   TBranch        *b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL1sMu5EG23IorMu7EG23IorMu20EG17IorMu23EG10_et;   //!
   TBranch        *b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL1fL1sMu5EG20orMu20EG15L1Filtered5_eta;   //!
   TBranch        *b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL1fL1sMu5EG20orMu20EG15L1Filtered5_phi;   //!
   TBranch        *b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL1fL1sMu5EG20orMu20EG15L1Filtered5_et;   //!
   TBranch        *b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_-hltL2fL1sMu5EG20orMu20EG15L1f5L2NoVtxFiltered16_eta;   //!
   TBranch        *b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_-hltL2fL1sMu5EG20orMu20EG15L1f5L2NoVtxFiltered16_phi;   //!
   TBranch        *b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_-hltL2fL1sMu5EG20orMu20EG15L1f5L2NoVtxFiltered16_et;   //!
   TBranch        *b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL3fL1sMu5EG20orMu20EG15L1f5L2NVf16L3NoFiltersNoVtxFiltered43_eta;   //!
   TBranch        *b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL3fL1sMu5EG20orMu20EG15L1f5L2NVf16L3NoFiltersNoVtxFiltered43_phi;   //!
   TBranch        *b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL3fL1sMu5EG20orMu20EG15L1f5L2NVf16L3NoFiltersNoVtxFiltered43_et;   //!
   TBranch        *b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltEGL1Mu5EG20Filter_eta;   //!
   TBranch        *b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltEGL1Mu5EG20Filter_phi;   //!
   TBranch        *b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltEGL1Mu5EG20Filter_et;   //!
   TBranch        *b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLEtFilter_eta;   //!
   TBranch        *b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLEtFilter_phi;   //!
   TBranch        *b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLEtFilter_et;   //!
   TBranch        *b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLClusterShapeFilter_eta;   //!
   TBranch        *b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLClusterShapeFilter_phi;   //!
   TBranch        *b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLClusterShapeFilter_et;   //!
   TBranch        *b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLHEFilter_eta;   //!
   TBranch        *b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLHEFilter_phi;   //!
   TBranch        *b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLHEFilter_et;   //!
   TBranch        *b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_accept;   //!
   TBranch        *b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL1sMu5EG23IorMu7EG23IorMu20EG17IorMu23EG10_eta;   //!
   TBranch        *b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL1sMu5EG23IorMu7EG23IorMu20EG17IorMu23EG10_phi;   //!
   TBranch        *b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL1sMu5EG23IorMu7EG23IorMu20EG17IorMu23EG10_et;   //!
   TBranch        *b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL1fL1sMu5EG20orMu20EG15L1Filtered5_eta;   //!
   TBranch        *b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL1fL1sMu5EG20orMu20EG15L1Filtered5_phi;   //!
   TBranch        *b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL1fL1sMu5EG20orMu20EG15L1Filtered5_et;   //!
   TBranch        *b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_-hltL2fL1sMu5EG20orMu20EG15L1f5L2NoVtxFiltered16_eta;   //!
   TBranch        *b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_-hltL2fL1sMu5EG20orMu20EG15L1f5L2NoVtxFiltered16_phi;   //!
   TBranch        *b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_-hltL2fL1sMu5EG20orMu20EG15L1f5L2NoVtxFiltered16_et;   //!
   TBranch        *b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL3fL1sMu5EG20orMu20EG15L1f5L2NVf16L3NoFiltersNoVtxFiltered48_eta;   //!
   TBranch        *b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL3fL1sMu5EG20orMu20EG15L1f5L2NVf16L3NoFiltersNoVtxFiltered48_phi;   //!
   TBranch        *b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL3fL1sMu5EG20orMu20EG15L1f5L2NVf16L3NoFiltersNoVtxFiltered48_et;   //!
   TBranch        *b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltEGL1Mu5EG20Filter_eta;   //!
   TBranch        *b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltEGL1Mu5EG20Filter_phi;   //!
   TBranch        *b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltEGL1Mu5EG20Filter_et;   //!
   TBranch        *b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLEtFilter_eta;   //!
   TBranch        *b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLEtFilter_phi;   //!
   TBranch        *b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLEtFilter_et;   //!
   TBranch        *b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLClusterShapeFilter_eta;   //!
   TBranch        *b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLClusterShapeFilter_phi;   //!
   TBranch        *b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLClusterShapeFilter_et;   //!
   TBranch        *b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLHEFilter_eta;   //!
   TBranch        *b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLHEFilter_phi;   //!
   TBranch        *b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLHEFilter_et;   //!
   TBranch        *b_trig_HLT_DoubleMu20_7_Mass0to30_L1_DM4_accept;   //!
   TBranch        *b_trig_HLT_DoubleMu20_7_Mass0to30_L1_DM4EG_accept;   //!
   TBranch        *b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_accept;   //!
   TBranch        *b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltL1sDoubleMu4OSEG12ORDoubleMu5OSEG12_eta;   //!
   TBranch        *b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltL1sDoubleMu4OSEG12ORDoubleMu5OSEG12_phi;   //!
   TBranch        *b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltL1sDoubleMu4OSEG12ORDoubleMu5OSEG12_et;   //!
   TBranch        *b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltL1fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1Filtered0_eta;   //!
   TBranch        *b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltL1fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1Filtered0_phi;   //!
   TBranch        *b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltL1fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1Filtered0_et;   //!
   TBranch        *b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_-hltL2fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1f0L2PreFiltered0_eta;   //!
   TBranch        *b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_-hltL2fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1f0L2PreFiltered0_phi;   //!
   TBranch        *b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_-hltL2fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1f0L2PreFiltered0_et;   //!
   TBranch        *b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_-hltL1fForIterL3L1fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1Filtered0_eta;   //!
   TBranch        *b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_-hltL1fForIterL3L1fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1Filtered0_phi;   //!
   TBranch        *b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_-hltL1fForIterL3L1fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1Filtered0_et;   //!
   TBranch        *b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltDoubleMu207Mass0to30Photon23L3Filtered_eta;   //!
   TBranch        *b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltDoubleMu207Mass0to30Photon23L3Filtered_phi;   //!
   TBranch        *b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltDoubleMu207Mass0to30Photon23L3Filtered_et;   //!
   TBranch        *b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEGDoubleMu4OSEG12ORDoubleMu5OSEG12Filter_eta;   //!
   TBranch        *b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEGDoubleMu4OSEG12ORDoubleMu5OSEG12Filter_phi;   //!
   TBranch        *b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEGDoubleMu4OSEG12ORDoubleMu5OSEG12Filter_et;   //!
   TBranch        *b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEG23EtFilter_eta;   //!
   TBranch        *b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEG23EtFilter_phi;   //!
   TBranch        *b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEG23EtFilter_et;   //!
   TBranch        *b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEG23HEFilter_eta;   //!
   TBranch        *b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEG23HEFilter_phi;   //!
   TBranch        *b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEG23HEFilter_et;   //!
   TBranch        *b_trig_HLT_Dimuon14_Phi_Barrel_Seagulls_accept;   //!
   TBranch        *b_trig_HLT_Dimuon24_Phi_noCorrL1_accept;   //!
   TBranch        *b_trig_HLT_DoubleIsoMu20_eta2p1_accept;   //!
   TBranch        *b_trig_HLT_DoubleIsoMu24_eta2p1_accept;   //!
   TBranch        *b_trig_HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx_accept;   //!
   TBranch        *b_trig_HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx_accept;   //!
   TBranch        *b_trig_HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx_accept;   //!
   TBranch        *b_trig_HLT_Mu8_accept;   //!
   TBranch        *b_trig_HLT_Mu17_accept;   //!
   TBranch        *b_trig_HLT_Mu19_accept;   //!
   TBranch        *b_trig_HLT_Mu17_Photon30_IsoCaloId_accept;   //!
   TBranch        *b_trig_HLT_Mu17_Photon30_IsoCaloId_hltL1sMu5EG23IorMu5IsoEG20IorMu7EG23IorMu7IsoEG20IorMuIso7EG23_eta;   //!
   TBranch        *b_trig_HLT_Mu17_Photon30_IsoCaloId_hltL1sMu5EG23IorMu5IsoEG20IorMu7EG23IorMu7IsoEG20IorMuIso7EG23_phi;   //!
   TBranch        *b_trig_HLT_Mu17_Photon30_IsoCaloId_hltL1sMu5EG23IorMu5IsoEG20IorMu7EG23IorMu7IsoEG20IorMuIso7EG23_et;   //!
   TBranch        *b_trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdMuonlegL1Filtered7_eta;   //!
   TBranch        *b_trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdMuonlegL1Filtered7_phi;   //!
   TBranch        *b_trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdMuonlegL1Filtered7_et;   //!
   TBranch        *b_trig_HLT_Mu17_Photon30_IsoCaloId_-hltMu17Photon30IsoCaloIdMuonlegL2Filtered10_eta;   //!
   TBranch        *b_trig_HLT_Mu17_Photon30_IsoCaloId_-hltMu17Photon30IsoCaloIdMuonlegL2Filtered10_phi;   //!
   TBranch        *b_trig_HLT_Mu17_Photon30_IsoCaloId_-hltMu17Photon30IsoCaloIdMuonlegL2Filtered10_et;   //!
   TBranch        *b_trig_HLT_Mu17_Photon30_IsoCaloId_-hltL1fForIterL3Mu17Photon30IsoCaloIdMuonlegL1Filtered7_eta;   //!
   TBranch        *b_trig_HLT_Mu17_Photon30_IsoCaloId_-hltL1fForIterL3Mu17Photon30IsoCaloIdMuonlegL1Filtered7_phi;   //!
   TBranch        *b_trig_HLT_Mu17_Photon30_IsoCaloId_-hltL1fForIterL3Mu17Photon30IsoCaloIdMuonlegL1Filtered7_et;   //!
   TBranch        *b_trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdMuonlegL3Filtered17Q_eta;   //!
   TBranch        *b_trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdMuonlegL3Filtered17Q_phi;   //!
   TBranch        *b_trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdMuonlegL3Filtered17Q_et;   //!
   TBranch        *b_trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegFilter_eta;   //!
   TBranch        *b_trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegFilter_phi;   //!
   TBranch        *b_trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegFilter_et;   //!
   TBranch        *b_trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegEtFilter_eta;   //!
   TBranch        *b_trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegEtFilter_phi;   //!
   TBranch        *b_trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegEtFilter_et;   //!
   TBranch        *b_trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegClusterShapeFilter_eta;   //!
   TBranch        *b_trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegClusterShapeFilter_phi;   //!
   TBranch        *b_trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegClusterShapeFilter_et;   //!
   TBranch        *b_trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegHEFilter_eta;   //!
   TBranch        *b_trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegHEFilter_phi;   //!
   TBranch        *b_trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegHEFilter_et;   //!
   TBranch        *b_trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegEcalIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegEcalIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegEcalIsoFilter_et;   //!
   TBranch        *b_trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegTrackIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegTrackIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegTrackIsoFilter_et;   //!
   TBranch        *b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_accept;   //!
   TBranch        *b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltL1sSingleEGNonIsoOrWithJetAndTau_eta;   //!
   TBranch        *b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltL1sSingleEGNonIsoOrWithJetAndTau_phi;   //!
   TBranch        *b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltL1sSingleEGNonIsoOrWithJetAndTau_et;   //!
   TBranch        *b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_et;   //!
   TBranch        *b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115EtFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115EtFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115EtFilter_et;   //!
   TBranch        *b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115CaloIdVTClusterShapeFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115CaloIdVTClusterShapeFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115CaloIdVTClusterShapeFilter_et;   //!
   TBranch        *b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115CaloIdVTHEFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115CaloIdVTHEFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115CaloIdVTHEFilter_et;   //!
   TBranch        *b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTPixelMatchFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTPixelMatchFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTPixelMatchFilter_et;   //!
   TBranch        *b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTGsfTrkIdTGsfDetaFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTGsfTrkIdTGsfDetaFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTGsfTrkIdTGsfDetaFilter_et;   //!
   TBranch        *b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTGsfTrkIdTGsfDphiFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTGsfTrkIdTGsfDphiFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTGsfTrkIdTGsfDphiFilter_et;   //!
   TBranch        *b_trig_HLT_Ele135_CaloIdVT_GsfTrkIdT_accept;   //!
   TBranch        *b_trig_HLT_Ele145_CaloIdVT_GsfTrkIdT_accept;   //!
   TBranch        *b_trig_HLT_Ele200_CaloIdVT_GsfTrkIdT_accept;   //!
   TBranch        *b_trig_HLT_Ele250_CaloIdVT_GsfTrkIdT_accept;   //!
   TBranch        *b_trig_HLT_Ele300_CaloIdVT_GsfTrkIdT_accept;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_accept;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltL1sSingleEGNonIsoOrWithJetAndTauNoPS_eta;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltL1sSingleEGNonIsoOrWithJetAndTauNoPS_phi;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltL1sSingleEGNonIsoOrWithJetAndTauNoPS_et;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEGL1SingleEGNonIsoOrWithJetAndTauNoPSFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEGL1SingleEGNonIsoOrWithJetAndTauNoPSFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEGL1SingleEGNonIsoOrWithJetAndTauNoPSFilter_et;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60EtFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60EtFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60EtFilter_et;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60HEFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60HEFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60HEFilter_et;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLR9IdFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLR9IdFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLR9IdFilter_et;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLClusterShapeFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLClusterShapeFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLClusterShapeFilter_et;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLEcalPFClusterIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLEcalPFClusterIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLEcalPFClusterIsoFilter_et;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLHcalPFClusterIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLHcalPFClusterIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLHcalPFClusterIsoFilter_et;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLHollowTrackIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLHollowTrackIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLHollowTrackIsoFilter_et;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_accept;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltL1sSingleEGNonIsoOrWithJetAndTauNoPS_eta;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltL1sSingleEGNonIsoOrWithJetAndTauNoPS_phi;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltL1sSingleEGNonIsoOrWithJetAndTauNoPS_et;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEGL1SingleEGNonIsoOrWithJetAndTauNoPSFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEGL1SingleEGNonIsoOrWithJetAndTauNoPSFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEGL1SingleEGNonIsoOrWithJetAndTauNoPSFilter_et;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60EtFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60EtFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60EtFilter_et;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60HEFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60HEFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60HEFilter_et;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLR9IdFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLR9IdFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLR9IdFilter_et;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLClusterShapeFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLClusterShapeFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLClusterShapeFilter_et;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLEcalPFClusterIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLEcalPFClusterIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLEcalPFClusterIsoFilter_et;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLHcalPFClusterIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLHcalPFClusterIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLHcalPFClusterIsoFilter_et;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLHollowTrackIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLHollowTrackIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLHollowTrackIsoFilter_et;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLDisplacedIdFilter_eta;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLDisplacedIdFilter_phi;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLDisplacedIdFilter_et;   //!
   TBranch        *b_trig_HLT_SinglePhoton10_Eta3p1ForPPRef_accept;   //!
   TBranch        *b_trig_HLT_SinglePhoton10_Eta3p1ForPPRef_hltL1sZeroBias_eta;   //!
   TBranch        *b_trig_HLT_SinglePhoton10_Eta3p1ForPPRef_hltL1sZeroBias_phi;   //!
   TBranch        *b_trig_HLT_SinglePhoton10_Eta3p1ForPPRef_hltL1sZeroBias_et;   //!
   TBranch        *b_trig_HLT_SinglePhoton10_Eta3p1ForPPRef_hltHIPhoton10Eta3p1_eta;   //!
   TBranch        *b_trig_HLT_SinglePhoton10_Eta3p1ForPPRef_hltHIPhoton10Eta3p1_phi;   //!
   TBranch        *b_trig_HLT_SinglePhoton10_Eta3p1ForPPRef_hltHIPhoton10Eta3p1_et;   //!
   TBranch        *b_trig_HLT_SinglePhoton20_Eta3p1ForPPRef_accept;   //!
   TBranch        *b_trig_HLT_SinglePhoton20_Eta3p1ForPPRef_hltL1sSingleEG10_eta;   //!
   TBranch        *b_trig_HLT_SinglePhoton20_Eta3p1ForPPRef_hltL1sSingleEG10_phi;   //!
   TBranch        *b_trig_HLT_SinglePhoton20_Eta3p1ForPPRef_hltL1sSingleEG10_et;   //!
   TBranch        *b_trig_HLT_SinglePhoton20_Eta3p1ForPPRef_hltHIPhoton20Eta3p1_eta;   //!
   TBranch        *b_trig_HLT_SinglePhoton20_Eta3p1ForPPRef_hltHIPhoton20Eta3p1_phi;   //!
   TBranch        *b_trig_HLT_SinglePhoton20_Eta3p1ForPPRef_hltHIPhoton20Eta3p1_et;   //!
   TBranch        *b_trig_HLT_SinglePhoton30_Eta3p1ForPPRef_accept;   //!
   TBranch        *b_trig_HLT_SinglePhoton30_Eta3p1ForPPRef_hltL1sSingleEG10_eta;   //!
   TBranch        *b_trig_HLT_SinglePhoton30_Eta3p1ForPPRef_hltL1sSingleEG10_phi;   //!
   TBranch        *b_trig_HLT_SinglePhoton30_Eta3p1ForPPRef_hltL1sSingleEG10_et;   //!
   TBranch        *b_trig_HLT_SinglePhoton30_Eta3p1ForPPRef_hltHIPhoton30Eta3p1_eta;   //!
   TBranch        *b_trig_HLT_SinglePhoton30_Eta3p1ForPPRef_hltHIPhoton30Eta3p1_phi;   //!
   TBranch        *b_trig_HLT_SinglePhoton30_Eta3p1ForPPRef_hltHIPhoton30Eta3p1_et;   //!
   TBranch        *b_trig_HLT_Photon20_HoverELoose_accept;   //!
   TBranch        *b_trig_HLT_Photon20_HoverELoose_hltL1sSingleEG10_eta;   //!
   TBranch        *b_trig_HLT_Photon20_HoverELoose_hltL1sSingleEG10_phi;   //!
   TBranch        *b_trig_HLT_Photon20_HoverELoose_hltL1sSingleEG10_et;   //!
   TBranch        *b_trig_HLT_Photon20_HoverELoose_hltEGL1SingleEG10Filter_eta;   //!
   TBranch        *b_trig_HLT_Photon20_HoverELoose_hltEGL1SingleEG10Filter_phi;   //!
   TBranch        *b_trig_HLT_Photon20_HoverELoose_hltEGL1SingleEG10Filter_et;   //!
   TBranch        *b_trig_HLT_Photon20_HoverELoose_hltEG20EtFilterLooseHoverE_eta;   //!
   TBranch        *b_trig_HLT_Photon20_HoverELoose_hltEG20EtFilterLooseHoverE_phi;   //!
   TBranch        *b_trig_HLT_Photon20_HoverELoose_hltEG20EtFilterLooseHoverE_et;   //!
   TBranch        *b_trig_HLT_Photon20_HoverELoose_hltEG20HEFilterLooseHoverE_eta;   //!
   TBranch        *b_trig_HLT_Photon20_HoverELoose_hltEG20HEFilterLooseHoverE_phi;   //!
   TBranch        *b_trig_HLT_Photon20_HoverELoose_hltEG20HEFilterLooseHoverE_et;   //!
   TBranch        *b_trig_HLT_Photon30_HoverELoose_accept;   //!
   TBranch        *b_trig_HLT_Photon30_HoverELoose_hltL1sSingleEG10_eta;   //!
   TBranch        *b_trig_HLT_Photon30_HoverELoose_hltL1sSingleEG10_phi;   //!
   TBranch        *b_trig_HLT_Photon30_HoverELoose_hltL1sSingleEG10_et;   //!
   TBranch        *b_trig_HLT_Photon30_HoverELoose_hltEGL1SingleEG10Filter_eta;   //!
   TBranch        *b_trig_HLT_Photon30_HoverELoose_hltEGL1SingleEG10Filter_phi;   //!
   TBranch        *b_trig_HLT_Photon30_HoverELoose_hltEGL1SingleEG10Filter_et;   //!
   TBranch        *b_trig_HLT_Photon30_HoverELoose_hltEG30EtFilterLooseHoverE_eta;   //!
   TBranch        *b_trig_HLT_Photon30_HoverELoose_hltEG30EtFilterLooseHoverE_phi;   //!
   TBranch        *b_trig_HLT_Photon30_HoverELoose_hltEG30EtFilterLooseHoverE_et;   //!
   TBranch        *b_trig_HLT_Photon30_HoverELoose_hltEG30HEFilterLooseHoverE_eta;   //!
   TBranch        *b_trig_HLT_Photon30_HoverELoose_hltEG30HEFilterLooseHoverE_phi;   //!
   TBranch        *b_trig_HLT_Photon30_HoverELoose_hltEG30HEFilterLooseHoverE_et;   //!
   TBranch        *b_trig_HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142_accept;   //!
   TBranch        *b_trig_AlCa_RPCMuonNormalisation_accept;   //!
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

   Nochecl(TTree *tree=0);
   virtual ~Nochecl();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Nochecl_cxx
Nochecl::Nochecl(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("outfile_956.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("outfile_956.root");
      }
      f->GetObject("IIHEAnalysis",tree);

   }
   Init(tree);
}

Nochecl::~Nochecl()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Nochecl::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Nochecl::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Nochecl::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
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
   trig_HLT_DoubleEle25_CaloIdL_MW_hltL1sSingleAndDoubleEGNonIsoOrWithEG26WithJetAndTau_eta = 0;
   trig_HLT_DoubleEle25_CaloIdL_MW_hltL1sSingleAndDoubleEGNonIsoOrWithEG26WithJetAndTau_phi = 0;
   trig_HLT_DoubleEle25_CaloIdL_MW_hltL1sSingleAndDoubleEGNonIsoOrWithEG26WithJetAndTau_et = 0;
   trig_HLT_DoubleEle25_CaloIdL_MW_hltEGL1SingleAndDoubleEGNonIsoOrWithEG26WithJetAndTauFilter_eta = 0;
   trig_HLT_DoubleEle25_CaloIdL_MW_hltEGL1SingleAndDoubleEGNonIsoOrWithEG26WithJetAndTauFilter_phi = 0;
   trig_HLT_DoubleEle25_CaloIdL_MW_hltEGL1SingleAndDoubleEGNonIsoOrWithEG26WithJetAndTauFilter_et = 0;
   trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25EtFilter_eta = 0;
   trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25EtFilter_phi = 0;
   trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25EtFilter_et = 0;
   trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25HEFilter_eta = 0;
   trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25HEFilter_phi = 0;
   trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25HEFilter_et = 0;
   trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25CaloIdLClusterShapeFilter_eta = 0;
   trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25CaloIdLClusterShapeFilter_phi = 0;
   trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25CaloIdLClusterShapeFilter_et = 0;
   trig_HLT_DoubleEle25_CaloIdL_MW_hltEle25CaloIdLPixelMatchFilter_eta = 0;
   trig_HLT_DoubleEle25_CaloIdL_MW_hltEle25CaloIdLPixelMatchFilter_phi = 0;
   trig_HLT_DoubleEle25_CaloIdL_MW_hltEle25CaloIdLPixelMatchFilter_et = 0;
   trig_HLT_DoubleEle25_CaloIdL_MW_hltEle25CaloIdLMWPMS2Filter_eta = 0;
   trig_HLT_DoubleEle25_CaloIdL_MW_hltEle25CaloIdLMWPMS2Filter_phi = 0;
   trig_HLT_DoubleEle25_CaloIdL_MW_hltEle25CaloIdLMWPMS2Filter_et = 0;
   trig_HLT_DoubleEle25_CaloIdL_MW_hltEgammaCandidatesWrapperUnseeded_eta = 0;
   trig_HLT_DoubleEle25_CaloIdL_MW_hltEgammaCandidatesWrapperUnseeded_phi = 0;
   trig_HLT_DoubleEle25_CaloIdL_MW_hltEgammaCandidatesWrapperUnseeded_et = 0;
   trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25EtUnseededFilter_eta = 0;
   trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25EtUnseededFilter_phi = 0;
   trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25EtUnseededFilter_et = 0;
   trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25HEUnseededFilter_eta = 0;
   trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25HEUnseededFilter_phi = 0;
   trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25HEUnseededFilter_et = 0;
   trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25CaloIdLClusterShapeUnseededFilter_eta = 0;
   trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25CaloIdLClusterShapeUnseededFilter_phi = 0;
   trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25CaloIdLClusterShapeUnseededFilter_et = 0;
   trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEle25CaloIdLPixelMatchUnseededFilter_eta = 0;
   trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEle25CaloIdLPixelMatchUnseededFilter_phi = 0;
   trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEle25CaloIdLPixelMatchUnseededFilter_et = 0;
   trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEle25CaloIdLMWPMS2UnseededFilter_eta = 0;
   trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEle25CaloIdLMWPMS2UnseededFilter_phi = 0;
   trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEle25CaloIdLMWPMS2UnseededFilter_et = 0;
   trig_HLT_DoubleEle33_CaloIdL_MW_hltL1sSingleAndDoubleEGNonIsoOrWithEG26WithJetAndTau_eta = 0;
   trig_HLT_DoubleEle33_CaloIdL_MW_hltL1sSingleAndDoubleEGNonIsoOrWithEG26WithJetAndTau_phi = 0;
   trig_HLT_DoubleEle33_CaloIdL_MW_hltL1sSingleAndDoubleEGNonIsoOrWithEG26WithJetAndTau_et = 0;
   trig_HLT_DoubleEle33_CaloIdL_MW_hltEGL1SingleAndDoubleEGNonIsoOrWithEG26WithJetAndTauFilter_eta = 0;
   trig_HLT_DoubleEle33_CaloIdL_MW_hltEGL1SingleAndDoubleEGNonIsoOrWithEG26WithJetAndTauFilter_phi = 0;
   trig_HLT_DoubleEle33_CaloIdL_MW_hltEGL1SingleAndDoubleEGNonIsoOrWithEG26WithJetAndTauFilter_et = 0;
   trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33EtFilter_eta = 0;
   trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33EtFilter_phi = 0;
   trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33EtFilter_et = 0;
   trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33HEFilter_eta = 0;
   trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33HEFilter_phi = 0;
   trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33HEFilter_et = 0;
   trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33CaloIdLClusterShapeFilter_eta = 0;
   trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33CaloIdLClusterShapeFilter_phi = 0;
   trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33CaloIdLClusterShapeFilter_et = 0;
   trig_HLT_DoubleEle33_CaloIdL_MW_hltEle33CaloIdLPixelMatchFilter_eta = 0;
   trig_HLT_DoubleEle33_CaloIdL_MW_hltEle33CaloIdLPixelMatchFilter_phi = 0;
   trig_HLT_DoubleEle33_CaloIdL_MW_hltEle33CaloIdLPixelMatchFilter_et = 0;
   trig_HLT_DoubleEle33_CaloIdL_MW_hltEle33CaloIdLMWPMS2Filter_eta = 0;
   trig_HLT_DoubleEle33_CaloIdL_MW_hltEle33CaloIdLMWPMS2Filter_phi = 0;
   trig_HLT_DoubleEle33_CaloIdL_MW_hltEle33CaloIdLMWPMS2Filter_et = 0;
   trig_HLT_DoubleEle33_CaloIdL_MW_hltEgammaCandidatesWrapperUnseeded_eta = 0;
   trig_HLT_DoubleEle33_CaloIdL_MW_hltEgammaCandidatesWrapperUnseeded_phi = 0;
   trig_HLT_DoubleEle33_CaloIdL_MW_hltEgammaCandidatesWrapperUnseeded_et = 0;
   trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33EtUnseededFilter_eta = 0;
   trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33EtUnseededFilter_phi = 0;
   trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33EtUnseededFilter_et = 0;
   trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33HEUnseededFilter_eta = 0;
   trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33HEUnseededFilter_phi = 0;
   trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33HEUnseededFilter_et = 0;
   trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33CaloIdLClusterShapeUnseededFilter_eta = 0;
   trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33CaloIdLClusterShapeUnseededFilter_phi = 0;
   trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33CaloIdLClusterShapeUnseededFilter_et = 0;
   trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEle33CaloIdLPixelMatchUnseededFilter_eta = 0;
   trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEle33CaloIdLPixelMatchUnseededFilter_phi = 0;
   trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEle33CaloIdLPixelMatchUnseededFilter_et = 0;
   trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEle33CaloIdLMWPMS2UnseededFilter_eta = 0;
   trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEle33CaloIdLMWPMS2UnseededFilter_phi = 0;
   trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEle33CaloIdLMWPMS2UnseededFilter_et = 0;
   trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltL1sSingleAndDoubleEGor_eta = 0;
   trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltL1sSingleAndDoubleEGor_phi = 0;
   trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltL1sSingleAndDoubleEGor_et = 0;
   trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEGL1SingleAndDoubleEGOrSingleFilter_eta = 0;
   trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEGL1SingleAndDoubleEGOrSingleFilter_phi = 0;
   trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEGL1SingleAndDoubleEGOrSingleFilter_et = 0;
   trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEG27L1SingleAndDoubleEGEtFilter_eta = 0;
   trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEG27L1SingleAndDoubleEGEtFilter_phi = 0;
   trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEG27L1SingleAndDoubleEGEtFilter_et = 0;
   trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightClusterShapeFilter_eta = 0;
   trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightClusterShapeFilter_phi = 0;
   trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightClusterShapeFilter_et = 0;
   trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightHEFilter_eta = 0;
   trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightHEFilter_phi = 0;
   trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightHEFilter_et = 0;
   trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightEcalIsoFilter_eta = 0;
   trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightEcalIsoFilter_phi = 0;
   trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightEcalIsoFilter_et = 0;
   trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightHcalIsoFilter_eta = 0;
   trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightHcalIsoFilter_phi = 0;
   trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightHcalIsoFilter_et = 0;
   trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEG27L1SingleAndDoubleEGEtFilter_eta = 0;
   trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEG27L1SingleAndDoubleEGEtFilter_phi = 0;
   trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEG27L1SingleAndDoubleEGEtFilter_et = 0;
   trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightClusterShapeFilter_eta = 0;
   trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightClusterShapeFilter_phi = 0;
   trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightClusterShapeFilter_et = 0;
   trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightHEFilter_eta = 0;
   trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightHEFilter_phi = 0;
   trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightHEFilter_et = 0;
   trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightEcalIsoFilter_eta = 0;
   trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightEcalIsoFilter_phi = 0;
   trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightEcalIsoFilter_et = 0;
   trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightHcalIsoFilter_eta = 0;
   trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightHcalIsoFilter_phi = 0;
   trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightHcalIsoFilter_et = 0;
   trig_HLT_Ele32_WPTight_Gsf_hltL1sSingleEGor_eta = 0;
   trig_HLT_Ele32_WPTight_Gsf_hltL1sSingleEGor_phi = 0;
   trig_HLT_Ele32_WPTight_Gsf_hltL1sSingleEGor_et = 0;
   trig_HLT_Ele32_WPTight_Gsf_hltEGL1SingleEGOrFilter_eta = 0;
   trig_HLT_Ele32_WPTight_Gsf_hltEGL1SingleEGOrFilter_phi = 0;
   trig_HLT_Ele32_WPTight_Gsf_hltEGL1SingleEGOrFilter_et = 0;
   trig_HLT_Ele32_WPTight_Gsf_hltEG32L1SingleEGOrEtFilter_eta = 0;
   trig_HLT_Ele32_WPTight_Gsf_hltEG32L1SingleEGOrEtFilter_phi = 0;
   trig_HLT_Ele32_WPTight_Gsf_hltEG32L1SingleEGOrEtFilter_et = 0;
   trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightClusterShapeFilter_eta = 0;
   trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightClusterShapeFilter_phi = 0;
   trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightClusterShapeFilter_et = 0;
   trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightHEFilter_eta = 0;
   trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightHEFilter_phi = 0;
   trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightHEFilter_et = 0;
   trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightEcalIsoFilter_eta = 0;
   trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightEcalIsoFilter_phi = 0;
   trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightEcalIsoFilter_et = 0;
   trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightHcalIsoFilter_eta = 0;
   trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightHcalIsoFilter_phi = 0;
   trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightHcalIsoFilter_et = 0;
   trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightPixelMatchFilter_eta = 0;
   trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightPixelMatchFilter_phi = 0;
   trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightPixelMatchFilter_et = 0;
   trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightPMS2Filter_eta = 0;
   trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightPMS2Filter_phi = 0;
   trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightPMS2Filter_et = 0;
   trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfOneOEMinusOneOPFilter_eta = 0;
   trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfOneOEMinusOneOPFilter_phi = 0;
   trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfOneOEMinusOneOPFilter_et = 0;
   trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfMissingHitsFilter_eta = 0;
   trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfMissingHitsFilter_phi = 0;
   trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfMissingHitsFilter_et = 0;
   trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfDetaFilter_eta = 0;
   trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfDetaFilter_phi = 0;
   trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfDetaFilter_et = 0;
   trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfDphiFilter_eta = 0;
   trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfDphiFilter_phi = 0;
   trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfDphiFilter_et = 0;
   trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfTrackIsoFilter_eta = 0;
   trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfTrackIsoFilter_phi = 0;
   trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfTrackIsoFilter_et = 0;
   trig_HLT_Ele35_WPTight_Gsf_hltL1sSingleEGor_eta = 0;
   trig_HLT_Ele35_WPTight_Gsf_hltL1sSingleEGor_phi = 0;
   trig_HLT_Ele35_WPTight_Gsf_hltL1sSingleEGor_et = 0;
   trig_HLT_Ele35_WPTight_Gsf_hltEGL1SingleEGOrFilter_eta = 0;
   trig_HLT_Ele35_WPTight_Gsf_hltEGL1SingleEGOrFilter_phi = 0;
   trig_HLT_Ele35_WPTight_Gsf_hltEGL1SingleEGOrFilter_et = 0;
   trig_HLT_Ele35_WPTight_Gsf_hltEG35L1SingleEGOrEtFilter_eta = 0;
   trig_HLT_Ele35_WPTight_Gsf_hltEG35L1SingleEGOrEtFilter_phi = 0;
   trig_HLT_Ele35_WPTight_Gsf_hltEG35L1SingleEGOrEtFilter_et = 0;
   trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightClusterShapeFilter_eta = 0;
   trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightClusterShapeFilter_phi = 0;
   trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightClusterShapeFilter_et = 0;
   trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightHEFilter_eta = 0;
   trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightHEFilter_phi = 0;
   trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightHEFilter_et = 0;
   trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightEcalIsoFilter_eta = 0;
   trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightEcalIsoFilter_phi = 0;
   trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightEcalIsoFilter_et = 0;
   trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightHcalIsoFilter_eta = 0;
   trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightHcalIsoFilter_phi = 0;
   trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightHcalIsoFilter_et = 0;
   trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightPixelMatchFilter_eta = 0;
   trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightPixelMatchFilter_phi = 0;
   trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightPixelMatchFilter_et = 0;
   trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightPMS2Filter_eta = 0;
   trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightPMS2Filter_phi = 0;
   trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightPMS2Filter_et = 0;
   trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfOneOEMinusOneOPFilter_eta = 0;
   trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfOneOEMinusOneOPFilter_phi = 0;
   trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfOneOEMinusOneOPFilter_et = 0;
   trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfMissingHitsFilter_eta = 0;
   trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfMissingHitsFilter_phi = 0;
   trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfMissingHitsFilter_et = 0;
   trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfDetaFilter_eta = 0;
   trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfDetaFilter_phi = 0;
   trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfDetaFilter_et = 0;
   trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfDphiFilter_eta = 0;
   trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfDphiFilter_phi = 0;
   trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfDphiFilter_et = 0;
   trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfTrackIsoFilter_eta = 0;
   trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfTrackIsoFilter_phi = 0;
   trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfTrackIsoFilter_et = 0;
   trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltL1sAlCaSingleEle_eta = 0;
   trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltL1sAlCaSingleEle_phi = 0;
   trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltL1sAlCaSingleEle_et = 0;
   trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTHEFilter_eta = 0;
   trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTHEFilter_phi = 0;
   trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTHEFilter_et = 0;
   trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTEcalIsoFilter_eta = 0;
   trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTEcalIsoFilter_phi = 0;
   trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTEcalIsoFilter_et = 0;
   trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTHcalIsoFilter_eta = 0;
   trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTHcalIsoFilter_phi = 0;
   trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTHcalIsoFilter_et = 0;
   trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTPixelMatchFilter_eta = 0;
   trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTPixelMatchFilter_phi = 0;
   trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTPixelMatchFilter_et = 0;
   trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTPMS2Filter_eta = 0;
   trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTPMS2Filter_phi = 0;
   trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTPMS2Filter_et = 0;
   trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTDetaFilter_eta = 0;
   trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTDetaFilter_phi = 0;
   trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTDetaFilter_et = 0;
   trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTDphiFilter_eta = 0;
   trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTDphiFilter_phi = 0;
   trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTDphiFilter_et = 0;
   trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTTrackIsoFilter_eta = 0;
   trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTTrackIsoFilter_phi = 0;
   trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTTrackIsoFilter_et = 0;
   trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltL1sSingleAndDoubleEGor_eta = 0;
   trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltL1sSingleAndDoubleEGor_phi = 0;
   trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltL1sSingleAndDoubleEGor_et = 0;
   trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEGL1SingleAndDoubleEGOrSingleFilter_eta = 0;
   trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEGL1SingleAndDoubleEGOrSingleFilter_phi = 0;
   trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEGL1SingleAndDoubleEGOrSingleFilter_et = 0;
   trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEG32L1SingleAndDoubleEGEtFilter_eta = 0;
   trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEG32L1SingleAndDoubleEGEtFilter_phi = 0;
   trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEG32L1SingleAndDoubleEGEtFilter_et = 0;
   trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightClusterShapeFilter_eta = 0;
   trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightClusterShapeFilter_phi = 0;
   trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightClusterShapeFilter_et = 0;
   trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightHEFilter_eta = 0;
   trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightHEFilter_phi = 0;
   trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightHEFilter_et = 0;
   trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightEcalIsoFilter_eta = 0;
   trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightEcalIsoFilter_phi = 0;
   trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightEcalIsoFilter_et = 0;
   trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightHcalIsoFilter_eta = 0;
   trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightHcalIsoFilter_phi = 0;
   trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightHcalIsoFilter_et = 0;
   trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightPixelMatchFilter_eta = 0;
   trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightPixelMatchFilter_phi = 0;
   trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightPixelMatchFilter_et = 0;
   trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightPMS2Filter_eta = 0;
   trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightPMS2Filter_phi = 0;
   trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightPMS2Filter_et = 0;
   trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfOneOEMinusOneOPFilter_eta = 0;
   trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfOneOEMinusOneOPFilter_phi = 0;
   trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfOneOEMinusOneOPFilter_et = 0;
   trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfMissingHitsFilter_eta = 0;
   trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfMissingHitsFilter_phi = 0;
   trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfMissingHitsFilter_et = 0;
   trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfDetaFilter_eta = 0;
   trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfDetaFilter_phi = 0;
   trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfDetaFilter_et = 0;
   trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfDphiFilter_eta = 0;
   trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfDphiFilter_phi = 0;
   trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfDphiFilter_et = 0;
   trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfTrackIsoFilter_eta = 0;
   trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfTrackIsoFilter_phi = 0;
   trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfTrackIsoFilter_et = 0;
   trig_HLT_Photon300_NoHE_hltL1sSingleEGNonIsoOrWithJetAndTau_eta = 0;
   trig_HLT_Photon300_NoHE_hltL1sSingleEGNonIsoOrWithJetAndTau_phi = 0;
   trig_HLT_Photon300_NoHE_hltL1sSingleEGNonIsoOrWithJetAndTau_et = 0;
   trig_HLT_Photon300_NoHE_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_eta = 0;
   trig_HLT_Photon300_NoHE_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_phi = 0;
   trig_HLT_Photon300_NoHE_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_et = 0;
   trig_HLT_Photon300_NoHE_hltEG300erEtFilter_eta = 0;
   trig_HLT_Photon300_NoHE_hltEG300erEtFilter_phi = 0;
   trig_HLT_Photon300_NoHE_hltEG300erEtFilter_et = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltL1sSingleAndDoubleEG_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltL1sSingleAndDoubleEG_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltL1sSingleAndDoubleEG_et = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEGL1SingleAndDoubleEGOrPairFilter_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEGL1SingleAndDoubleEGOrPairFilter_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEGL1SingleAndDoubleEGOrPairFilter_et = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg1Filter_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg1Filter_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg1Filter_et = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg2Filter_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg2Filter_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg2Filter_et = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg1Filter_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg1Filter_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg1Filter_et = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg2Filter_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg2Filter_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg2Filter_et = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg1Filter_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg1Filter_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg1Filter_et = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg2Filter_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg2Filter_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg2Filter_et = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg1Filter_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg1Filter_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg1Filter_et = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg2Filter_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg2Filter_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg2Filter_et = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg1Filter_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg1Filter_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg1Filter_et = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg2Filter_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg2Filter_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg2Filter_et = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg1Filter_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg1Filter_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg1Filter_et = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg2Filter_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg2Filter_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg2Filter_et = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg1Filter_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg1Filter_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg1Filter_et = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg2Filter_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg2Filter_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg2Filter_et = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg1Filter_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg1Filter_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg1Filter_et = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg2Filter_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg2Filter_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg2Filter_et = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg1Filter_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg1Filter_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg1Filter_et = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg2Filter_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg2Filter_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg2Filter_et = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_et = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_et = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter_et = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltL1sSingleAndDoubleEG_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltL1sSingleAndDoubleEG_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltL1sSingleAndDoubleEG_et = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEGL1SingleAndDoubleEGOrPairFilter_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEGL1SingleAndDoubleEGOrPairFilter_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEGL1SingleAndDoubleEGOrPairFilter_et = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg1Filter_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg1Filter_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg1Filter_et = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg2Filter_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg2Filter_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg2Filter_et = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg1Filter_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg1Filter_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg1Filter_et = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg2Filter_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg2Filter_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg2Filter_et = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg1Filter_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg1Filter_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg1Filter_et = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg2Filter_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg2Filter_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg2Filter_et = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg1Filter_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg1Filter_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg1Filter_et = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg2Filter_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg2Filter_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg2Filter_et = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg1Filter_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg1Filter_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg1Filter_et = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg2Filter_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg2Filter_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg2Filter_et = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg1Filter_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg1Filter_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg1Filter_et = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg2Filter_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg2Filter_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg2Filter_et = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg1Filter_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg1Filter_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg1Filter_et = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg2Filter_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg2Filter_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg2Filter_et = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg1Filter_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg1Filter_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg1Filter_et = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg2Filter_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg2Filter_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg2Filter_et = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg1Filter_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg1Filter_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg1Filter_et = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg2Filter_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg2Filter_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg2Filter_et = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_et = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_et = 0;
   trig_HLT_Photon20_hltL1sSingleEG15_eta = 0;
   trig_HLT_Photon20_hltL1sSingleEG15_phi = 0;
   trig_HLT_Photon20_hltL1sSingleEG15_et = 0;
   trig_HLT_Photon20_hltEGL1SingleEG15Filter_eta = 0;
   trig_HLT_Photon20_hltEGL1SingleEG15Filter_phi = 0;
   trig_HLT_Photon20_hltEGL1SingleEG15Filter_et = 0;
   trig_HLT_Photon20_hltEG20L1EG15EtFilter_eta = 0;
   trig_HLT_Photon20_hltEG20L1EG15EtFilter_phi = 0;
   trig_HLT_Photon20_hltEG20L1EG15EtFilter_et = 0;
   trig_HLT_Photon20_hltEG20L1EG10OR15HEFilter_eta = 0;
   trig_HLT_Photon20_hltEG20L1EG10OR15HEFilter_phi = 0;
   trig_HLT_Photon20_hltEG20L1EG10OR15HEFilter_et = 0;
   trig_HLT_Photon33_hltL1sSingleEG26_eta = 0;
   trig_HLT_Photon33_hltL1sSingleEG26_phi = 0;
   trig_HLT_Photon33_hltL1sSingleEG26_et = 0;
   trig_HLT_Photon33_hltEGL1SingleEG26Filter_eta = 0;
   trig_HLT_Photon33_hltEGL1SingleEG26Filter_phi = 0;
   trig_HLT_Photon33_hltEGL1SingleEG26Filter_et = 0;
   trig_HLT_Photon33_hltEG33L1EG26EtFilter_eta = 0;
   trig_HLT_Photon33_hltEG33L1EG26EtFilter_phi = 0;
   trig_HLT_Photon33_hltEG33L1EG26EtFilter_et = 0;
   trig_HLT_Photon33_hltEG33L1EG26HEFilter_eta = 0;
   trig_HLT_Photon33_hltEG33L1EG26HEFilter_phi = 0;
   trig_HLT_Photon33_hltEG33L1EG26HEFilter_et = 0;
   trig_HLT_Photon50_hltL1sSingleEG34to45_eta = 0;
   trig_HLT_Photon50_hltL1sSingleEG34to45_phi = 0;
   trig_HLT_Photon50_hltL1sSingleEG34to45_et = 0;
   trig_HLT_Photon50_hltEGL1SingleEG34to45Filter_eta = 0;
   trig_HLT_Photon50_hltEGL1SingleEG34to45Filter_phi = 0;
   trig_HLT_Photon50_hltEGL1SingleEG34to45Filter_et = 0;
   trig_HLT_Photon50_hltEG50EtFilter_eta = 0;
   trig_HLT_Photon50_hltEG50EtFilter_phi = 0;
   trig_HLT_Photon50_hltEG50EtFilter_et = 0;
   trig_HLT_Photon50_hltEG50HEFilter_eta = 0;
   trig_HLT_Photon50_hltEG50HEFilter_phi = 0;
   trig_HLT_Photon50_hltEG50HEFilter_et = 0;
   trig_HLT_Photon75_hltL1sSingleEG34to50_eta = 0;
   trig_HLT_Photon75_hltL1sSingleEG34to50_phi = 0;
   trig_HLT_Photon75_hltL1sSingleEG34to50_et = 0;
   trig_HLT_Photon75_hltEGL1SingleEG34ORL1SingleEG40Filter_eta = 0;
   trig_HLT_Photon75_hltEGL1SingleEG34ORL1SingleEG40Filter_phi = 0;
   trig_HLT_Photon75_hltEGL1SingleEG34ORL1SingleEG40Filter_et = 0;
   trig_HLT_Photon75_hltEG75EtFilter_eta = 0;
   trig_HLT_Photon75_hltEG75EtFilter_phi = 0;
   trig_HLT_Photon75_hltEG75EtFilter_et = 0;
   trig_HLT_Photon75_hltEG75HEFilter_eta = 0;
   trig_HLT_Photon75_hltEG75HEFilter_phi = 0;
   trig_HLT_Photon75_hltEG75HEFilter_et = 0;
   trig_HLT_Photon90_hltL1sSingleEG34to50_eta = 0;
   trig_HLT_Photon90_hltL1sSingleEG34to50_phi = 0;
   trig_HLT_Photon90_hltL1sSingleEG34to50_et = 0;
   trig_HLT_Photon90_hltEGL1SingleEG34ORL1SingleEG40Filter_eta = 0;
   trig_HLT_Photon90_hltEGL1SingleEG34ORL1SingleEG40Filter_phi = 0;
   trig_HLT_Photon90_hltEGL1SingleEG34ORL1SingleEG40Filter_et = 0;
   trig_HLT_Photon90_hltEG90EtFilter_eta = 0;
   trig_HLT_Photon90_hltEG90EtFilter_phi = 0;
   trig_HLT_Photon90_hltEG90EtFilter_et = 0;
   trig_HLT_Photon90_hltEG90HEFilter_eta = 0;
   trig_HLT_Photon90_hltEG90HEFilter_phi = 0;
   trig_HLT_Photon90_hltEG90HEFilter_et = 0;
   trig_HLT_Photon120_hltL1sSingleEG34to50_eta = 0;
   trig_HLT_Photon120_hltL1sSingleEG34to50_phi = 0;
   trig_HLT_Photon120_hltL1sSingleEG34to50_et = 0;
   trig_HLT_Photon120_hltEGL1SingleEG34ORL1SingleEG40Filter_eta = 0;
   trig_HLT_Photon120_hltEGL1SingleEG34ORL1SingleEG40Filter_phi = 0;
   trig_HLT_Photon120_hltEGL1SingleEG34ORL1SingleEG40Filter_et = 0;
   trig_HLT_Photon120_hltEG120EtFilter_eta = 0;
   trig_HLT_Photon120_hltEG120EtFilter_phi = 0;
   trig_HLT_Photon120_hltEG120EtFilter_et = 0;
   trig_HLT_Photon120_hltEG120HEFilter_eta = 0;
   trig_HLT_Photon120_hltEG120HEFilter_phi = 0;
   trig_HLT_Photon120_hltEG120HEFilter_et = 0;
   trig_HLT_Photon150_hltL1sSingleEGNonIsoOrWithJetAndTau_eta = 0;
   trig_HLT_Photon150_hltL1sSingleEGNonIsoOrWithJetAndTau_phi = 0;
   trig_HLT_Photon150_hltL1sSingleEGNonIsoOrWithJetAndTau_et = 0;
   trig_HLT_Photon150_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_eta = 0;
   trig_HLT_Photon150_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_phi = 0;
   trig_HLT_Photon150_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_et = 0;
   trig_HLT_Photon150_hltEG150EtFilter_eta = 0;
   trig_HLT_Photon150_hltEG150EtFilter_phi = 0;
   trig_HLT_Photon150_hltEG150EtFilter_et = 0;
   trig_HLT_Photon150_hltEG150HEFilter_eta = 0;
   trig_HLT_Photon150_hltEG150HEFilter_phi = 0;
   trig_HLT_Photon150_hltEG150HEFilter_et = 0;
   trig_HLT_Photon175_hltL1sSingleEGNonIsoOrWithJetAndTau_eta = 0;
   trig_HLT_Photon175_hltL1sSingleEGNonIsoOrWithJetAndTau_phi = 0;
   trig_HLT_Photon175_hltL1sSingleEGNonIsoOrWithJetAndTau_et = 0;
   trig_HLT_Photon175_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_eta = 0;
   trig_HLT_Photon175_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_phi = 0;
   trig_HLT_Photon175_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_et = 0;
   trig_HLT_Photon175_hltEG175EtFilter_eta = 0;
   trig_HLT_Photon175_hltEG175EtFilter_phi = 0;
   trig_HLT_Photon175_hltEG175EtFilter_et = 0;
   trig_HLT_Photon175_hltEG175HEFilter_eta = 0;
   trig_HLT_Photon175_hltEG175HEFilter_phi = 0;
   trig_HLT_Photon175_hltEG175HEFilter_et = 0;
   trig_HLT_Photon200_hltL1sSingleEGNonIsoOrWithJetAndTau_eta = 0;
   trig_HLT_Photon200_hltL1sSingleEGNonIsoOrWithJetAndTau_phi = 0;
   trig_HLT_Photon200_hltL1sSingleEGNonIsoOrWithJetAndTau_et = 0;
   trig_HLT_Photon200_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_eta = 0;
   trig_HLT_Photon200_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_phi = 0;
   trig_HLT_Photon200_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_et = 0;
   trig_HLT_Photon200_hltEG200EtFilter_eta = 0;
   trig_HLT_Photon200_hltEG200EtFilter_phi = 0;
   trig_HLT_Photon200_hltEG200EtFilter_et = 0;
   trig_HLT_Photon200_hltEG200HEFilter_eta = 0;
   trig_HLT_Photon200_hltEG200HEFilter_phi = 0;
   trig_HLT_Photon200_hltEG200HEFilter_et = 0;
   trig_HLT_Photon50_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_eta = 0;
   trig_HLT_Photon50_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_phi = 0;
   trig_HLT_Photon50_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_et = 0;
   trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_eta = 0;
   trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_phi = 0;
   trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_et = 0;
   trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMEtFilter_eta = 0;
   trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMEtFilter_phi = 0;
   trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMEtFilter_et = 0;
   trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMHEFilter_eta = 0;
   trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMHEFilter_phi = 0;
   trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMHEFilter_et = 0;
   trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMR9Filter_eta = 0;
   trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMR9Filter_phi = 0;
   trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMR9Filter_et = 0;
   trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMEcalIsoFilter_eta = 0;
   trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMEcalIsoFilter_phi = 0;
   trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMEcalIsoFilter_et = 0;
   trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMHcalIsoFilter_eta = 0;
   trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMHcalIsoFilter_phi = 0;
   trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMHcalIsoFilter_et = 0;
   trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMTrackIsoFilter_eta = 0;
   trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMTrackIsoFilter_phi = 0;
   trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMTrackIsoFilter_et = 0;
   trig_HLT_Photon75_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_eta = 0;
   trig_HLT_Photon75_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_phi = 0;
   trig_HLT_Photon75_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_et = 0;
   trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_eta = 0;
   trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_phi = 0;
   trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_et = 0;
   trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMEtFilter_eta = 0;
   trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMEtFilter_phi = 0;
   trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMEtFilter_et = 0;
   trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMHEFilter_eta = 0;
   trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMHEFilter_phi = 0;
   trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMHEFilter_et = 0;
   trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMR9Filter_eta = 0;
   trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMR9Filter_phi = 0;
   trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMR9Filter_et = 0;
   trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMEcalIsoFilter_eta = 0;
   trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMEcalIsoFilter_phi = 0;
   trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMEcalIsoFilter_et = 0;
   trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMHcalIsoFilter_eta = 0;
   trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMHcalIsoFilter_phi = 0;
   trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMHcalIsoFilter_et = 0;
   trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMTrackIsoFilter_eta = 0;
   trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMTrackIsoFilter_phi = 0;
   trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMTrackIsoFilter_et = 0;
   trig_HLT_Photon90_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_eta = 0;
   trig_HLT_Photon90_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_phi = 0;
   trig_HLT_Photon90_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_et = 0;
   trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_eta = 0;
   trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_phi = 0;
   trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_et = 0;
   trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMEtFilter_eta = 0;
   trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMEtFilter_phi = 0;
   trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMEtFilter_et = 0;
   trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMHEFilter_eta = 0;
   trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMHEFilter_phi = 0;
   trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMHEFilter_et = 0;
   trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMR9Filter_eta = 0;
   trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMR9Filter_phi = 0;
   trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMR9Filter_et = 0;
   trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMEcalIsoFilter_eta = 0;
   trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMEcalIsoFilter_phi = 0;
   trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMEcalIsoFilter_et = 0;
   trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMHcalIsoFilter_eta = 0;
   trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMHcalIsoFilter_phi = 0;
   trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMHcalIsoFilter_et = 0;
   trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMTrackIsoFilter_eta = 0;
   trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMTrackIsoFilter_phi = 0;
   trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMTrackIsoFilter_et = 0;
   trig_HLT_Photon120_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_eta = 0;
   trig_HLT_Photon120_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_phi = 0;
   trig_HLT_Photon120_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_et = 0;
   trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_eta = 0;
   trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_phi = 0;
   trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_et = 0;
   trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMEtFilter_eta = 0;
   trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMEtFilter_phi = 0;
   trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMEtFilter_et = 0;
   trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMHEFilter_eta = 0;
   trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMHEFilter_phi = 0;
   trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMHEFilter_et = 0;
   trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMR9Filter_eta = 0;
   trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMR9Filter_phi = 0;
   trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMR9Filter_et = 0;
   trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMEcalIsoFilter_eta = 0;
   trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMEcalIsoFilter_phi = 0;
   trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMEcalIsoFilter_et = 0;
   trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMHcalIsoFilter_eta = 0;
   trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMHcalIsoFilter_phi = 0;
   trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMHcalIsoFilter_et = 0;
   trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMTrackIsoFilter_eta = 0;
   trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMTrackIsoFilter_phi = 0;
   trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMTrackIsoFilter_et = 0;
   trig_HLT_Photon165_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_eta = 0;
   trig_HLT_Photon165_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_phi = 0;
   trig_HLT_Photon165_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_et = 0;
   trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_eta = 0;
   trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_phi = 0;
   trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_et = 0;
   trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMEtFilter_eta = 0;
   trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMEtFilter_phi = 0;
   trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMEtFilter_et = 0;
   trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMHEFilter_eta = 0;
   trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMHEFilter_phi = 0;
   trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMHEFilter_et = 0;
   trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMR9Filter_eta = 0;
   trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMR9Filter_phi = 0;
   trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMR9Filter_et = 0;
   trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMEcalIsoFilter_eta = 0;
   trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMEcalIsoFilter_phi = 0;
   trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMEcalIsoFilter_et = 0;
   trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMHcalIsoFilter_eta = 0;
   trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMHcalIsoFilter_phi = 0;
   trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMHcalIsoFilter_et = 0;
   trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMTrackIsoFilter_eta = 0;
   trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMTrackIsoFilter_phi = 0;
   trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMTrackIsoFilter_et = 0;
   trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL1sMu5EG23IorMu7EG23IorMu20EG17IorMu23EG10_eta = 0;
   trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL1sMu5EG23IorMu7EG23IorMu20EG17IorMu23EG10_phi = 0;
   trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL1sMu5EG23IorMu7EG23IorMu20EG17IorMu23EG10_et = 0;
   trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL1fL1sMu5EG20orMu20EG15L1Filtered5_eta = 0;
   trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL1fL1sMu5EG20orMu20EG15L1Filtered5_phi = 0;
   trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL1fL1sMu5EG20orMu20EG15L1Filtered5_et = 0;
   trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_-hltL2fL1sMu5EG20orMu20EG15L1f5L2NoVtxFiltered16_eta = 0;
   trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_-hltL2fL1sMu5EG20orMu20EG15L1f5L2NoVtxFiltered16_phi = 0;
   trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_-hltL2fL1sMu5EG20orMu20EG15L1f5L2NoVtxFiltered16_et = 0;
   trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL3fL1sMu5EG20orMu20EG15L1f5L2NVf16L3NoFiltersNoVtxFiltered43_eta = 0;
   trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL3fL1sMu5EG20orMu20EG15L1f5L2NVf16L3NoFiltersNoVtxFiltered43_phi = 0;
   trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL3fL1sMu5EG20orMu20EG15L1f5L2NVf16L3NoFiltersNoVtxFiltered43_et = 0;
   trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltEGL1Mu5EG20Filter_eta = 0;
   trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltEGL1Mu5EG20Filter_phi = 0;
   trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltEGL1Mu5EG20Filter_et = 0;
   trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLEtFilter_eta = 0;
   trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLEtFilter_phi = 0;
   trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLEtFilter_et = 0;
   trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLClusterShapeFilter_eta = 0;
   trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLClusterShapeFilter_phi = 0;
   trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLClusterShapeFilter_et = 0;
   trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLHEFilter_eta = 0;
   trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLHEFilter_phi = 0;
   trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLHEFilter_et = 0;
   trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL1sMu5EG23IorMu7EG23IorMu20EG17IorMu23EG10_eta = 0;
   trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL1sMu5EG23IorMu7EG23IorMu20EG17IorMu23EG10_phi = 0;
   trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL1sMu5EG23IorMu7EG23IorMu20EG17IorMu23EG10_et = 0;
   trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL1fL1sMu5EG20orMu20EG15L1Filtered5_eta = 0;
   trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL1fL1sMu5EG20orMu20EG15L1Filtered5_phi = 0;
   trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL1fL1sMu5EG20orMu20EG15L1Filtered5_et = 0;
   trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_-hltL2fL1sMu5EG20orMu20EG15L1f5L2NoVtxFiltered16_eta = 0;
   trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_-hltL2fL1sMu5EG20orMu20EG15L1f5L2NoVtxFiltered16_phi = 0;
   trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_-hltL2fL1sMu5EG20orMu20EG15L1f5L2NoVtxFiltered16_et = 0;
   trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL3fL1sMu5EG20orMu20EG15L1f5L2NVf16L3NoFiltersNoVtxFiltered48_eta = 0;
   trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL3fL1sMu5EG20orMu20EG15L1f5L2NVf16L3NoFiltersNoVtxFiltered48_phi = 0;
   trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL3fL1sMu5EG20orMu20EG15L1f5L2NVf16L3NoFiltersNoVtxFiltered48_et = 0;
   trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltEGL1Mu5EG20Filter_eta = 0;
   trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltEGL1Mu5EG20Filter_phi = 0;
   trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltEGL1Mu5EG20Filter_et = 0;
   trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLEtFilter_eta = 0;
   trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLEtFilter_phi = 0;
   trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLEtFilter_et = 0;
   trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLClusterShapeFilter_eta = 0;
   trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLClusterShapeFilter_phi = 0;
   trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLClusterShapeFilter_et = 0;
   trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLHEFilter_eta = 0;
   trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLHEFilter_phi = 0;
   trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLHEFilter_et = 0;
   trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltL1sDoubleMu4OSEG12ORDoubleMu5OSEG12_eta = 0;
   trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltL1sDoubleMu4OSEG12ORDoubleMu5OSEG12_phi = 0;
   trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltL1sDoubleMu4OSEG12ORDoubleMu5OSEG12_et = 0;
   trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltL1fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1Filtered0_eta = 0;
   trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltL1fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1Filtered0_phi = 0;
   trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltL1fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1Filtered0_et = 0;
   trig_HLT_DoubleMu20_7_Mass0to30_Photon23_-hltL2fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1f0L2PreFiltered0_eta = 0;
   trig_HLT_DoubleMu20_7_Mass0to30_Photon23_-hltL2fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1f0L2PreFiltered0_phi = 0;
   trig_HLT_DoubleMu20_7_Mass0to30_Photon23_-hltL2fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1f0L2PreFiltered0_et = 0;
   trig_HLT_DoubleMu20_7_Mass0to30_Photon23_-hltL1fForIterL3L1fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1Filtered0_eta = 0;
   trig_HLT_DoubleMu20_7_Mass0to30_Photon23_-hltL1fForIterL3L1fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1Filtered0_phi = 0;
   trig_HLT_DoubleMu20_7_Mass0to30_Photon23_-hltL1fForIterL3L1fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1Filtered0_et = 0;
   trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltDoubleMu207Mass0to30Photon23L3Filtered_eta = 0;
   trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltDoubleMu207Mass0to30Photon23L3Filtered_phi = 0;
   trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltDoubleMu207Mass0to30Photon23L3Filtered_et = 0;
   trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEGDoubleMu4OSEG12ORDoubleMu5OSEG12Filter_eta = 0;
   trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEGDoubleMu4OSEG12ORDoubleMu5OSEG12Filter_phi = 0;
   trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEGDoubleMu4OSEG12ORDoubleMu5OSEG12Filter_et = 0;
   trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEG23EtFilter_eta = 0;
   trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEG23EtFilter_phi = 0;
   trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEG23EtFilter_et = 0;
   trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEG23HEFilter_eta = 0;
   trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEG23HEFilter_phi = 0;
   trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEG23HEFilter_et = 0;
   trig_HLT_Mu17_Photon30_IsoCaloId_hltL1sMu5EG23IorMu5IsoEG20IorMu7EG23IorMu7IsoEG20IorMuIso7EG23_eta = 0;
   trig_HLT_Mu17_Photon30_IsoCaloId_hltL1sMu5EG23IorMu5IsoEG20IorMu7EG23IorMu7IsoEG20IorMuIso7EG23_phi = 0;
   trig_HLT_Mu17_Photon30_IsoCaloId_hltL1sMu5EG23IorMu5IsoEG20IorMu7EG23IorMu7IsoEG20IorMuIso7EG23_et = 0;
   trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdMuonlegL1Filtered7_eta = 0;
   trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdMuonlegL1Filtered7_phi = 0;
   trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdMuonlegL1Filtered7_et = 0;
   trig_HLT_Mu17_Photon30_IsoCaloId_-hltMu17Photon30IsoCaloIdMuonlegL2Filtered10_eta = 0;
   trig_HLT_Mu17_Photon30_IsoCaloId_-hltMu17Photon30IsoCaloIdMuonlegL2Filtered10_phi = 0;
   trig_HLT_Mu17_Photon30_IsoCaloId_-hltMu17Photon30IsoCaloIdMuonlegL2Filtered10_et = 0;
   trig_HLT_Mu17_Photon30_IsoCaloId_-hltL1fForIterL3Mu17Photon30IsoCaloIdMuonlegL1Filtered7_eta = 0;
   trig_HLT_Mu17_Photon30_IsoCaloId_-hltL1fForIterL3Mu17Photon30IsoCaloIdMuonlegL1Filtered7_phi = 0;
   trig_HLT_Mu17_Photon30_IsoCaloId_-hltL1fForIterL3Mu17Photon30IsoCaloIdMuonlegL1Filtered7_et = 0;
   trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdMuonlegL3Filtered17Q_eta = 0;
   trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdMuonlegL3Filtered17Q_phi = 0;
   trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdMuonlegL3Filtered17Q_et = 0;
   trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegFilter_eta = 0;
   trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegFilter_phi = 0;
   trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegFilter_et = 0;
   trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegEtFilter_eta = 0;
   trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegEtFilter_phi = 0;
   trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegEtFilter_et = 0;
   trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegClusterShapeFilter_eta = 0;
   trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegClusterShapeFilter_phi = 0;
   trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegClusterShapeFilter_et = 0;
   trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegHEFilter_eta = 0;
   trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegHEFilter_phi = 0;
   trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegHEFilter_et = 0;
   trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegEcalIsoFilter_eta = 0;
   trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegEcalIsoFilter_phi = 0;
   trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegEcalIsoFilter_et = 0;
   trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegTrackIsoFilter_eta = 0;
   trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegTrackIsoFilter_phi = 0;
   trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegTrackIsoFilter_et = 0;
   trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltL1sSingleEGNonIsoOrWithJetAndTau_eta = 0;
   trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltL1sSingleEGNonIsoOrWithJetAndTau_phi = 0;
   trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltL1sSingleEGNonIsoOrWithJetAndTau_et = 0;
   trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_eta = 0;
   trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_phi = 0;
   trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_et = 0;
   trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115EtFilter_eta = 0;
   trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115EtFilter_phi = 0;
   trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115EtFilter_et = 0;
   trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115CaloIdVTClusterShapeFilter_eta = 0;
   trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115CaloIdVTClusterShapeFilter_phi = 0;
   trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115CaloIdVTClusterShapeFilter_et = 0;
   trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115CaloIdVTHEFilter_eta = 0;
   trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115CaloIdVTHEFilter_phi = 0;
   trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115CaloIdVTHEFilter_et = 0;
   trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTPixelMatchFilter_eta = 0;
   trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTPixelMatchFilter_phi = 0;
   trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTPixelMatchFilter_et = 0;
   trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTGsfTrkIdTGsfDetaFilter_eta = 0;
   trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTGsfTrkIdTGsfDetaFilter_phi = 0;
   trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTGsfTrkIdTGsfDetaFilter_et = 0;
   trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTGsfTrkIdTGsfDphiFilter_eta = 0;
   trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTGsfTrkIdTGsfDphiFilter_phi = 0;
   trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTGsfTrkIdTGsfDphiFilter_et = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltL1sSingleEGNonIsoOrWithJetAndTauNoPS_eta = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltL1sSingleEGNonIsoOrWithJetAndTauNoPS_phi = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltL1sSingleEGNonIsoOrWithJetAndTauNoPS_et = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEGL1SingleEGNonIsoOrWithJetAndTauNoPSFilter_eta = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEGL1SingleEGNonIsoOrWithJetAndTauNoPSFilter_phi = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEGL1SingleEGNonIsoOrWithJetAndTauNoPSFilter_et = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60EtFilter_eta = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60EtFilter_phi = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60EtFilter_et = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60HEFilter_eta = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60HEFilter_phi = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60HEFilter_et = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLR9IdFilter_eta = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLR9IdFilter_phi = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLR9IdFilter_et = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLClusterShapeFilter_eta = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLClusterShapeFilter_phi = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLClusterShapeFilter_et = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLEcalPFClusterIsoFilter_eta = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLEcalPFClusterIsoFilter_phi = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLEcalPFClusterIsoFilter_et = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLHcalPFClusterIsoFilter_eta = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLHcalPFClusterIsoFilter_phi = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLHcalPFClusterIsoFilter_et = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLHollowTrackIsoFilter_eta = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLHollowTrackIsoFilter_phi = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLHollowTrackIsoFilter_et = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltL1sSingleEGNonIsoOrWithJetAndTauNoPS_eta = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltL1sSingleEGNonIsoOrWithJetAndTauNoPS_phi = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltL1sSingleEGNonIsoOrWithJetAndTauNoPS_et = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEGL1SingleEGNonIsoOrWithJetAndTauNoPSFilter_eta = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEGL1SingleEGNonIsoOrWithJetAndTauNoPSFilter_phi = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEGL1SingleEGNonIsoOrWithJetAndTauNoPSFilter_et = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60EtFilter_eta = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60EtFilter_phi = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60EtFilter_et = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60HEFilter_eta = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60HEFilter_phi = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60HEFilter_et = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLR9IdFilter_eta = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLR9IdFilter_phi = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLR9IdFilter_et = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLClusterShapeFilter_eta = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLClusterShapeFilter_phi = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLClusterShapeFilter_et = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLEcalPFClusterIsoFilter_eta = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLEcalPFClusterIsoFilter_phi = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLEcalPFClusterIsoFilter_et = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLHcalPFClusterIsoFilter_eta = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLHcalPFClusterIsoFilter_phi = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLHcalPFClusterIsoFilter_et = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLHollowTrackIsoFilter_eta = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLHollowTrackIsoFilter_phi = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLHollowTrackIsoFilter_et = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLDisplacedIdFilter_eta = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLDisplacedIdFilter_phi = 0;
   trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLDisplacedIdFilter_et = 0;
   trig_HLT_SinglePhoton10_Eta3p1ForPPRef_hltL1sZeroBias_eta = 0;
   trig_HLT_SinglePhoton10_Eta3p1ForPPRef_hltL1sZeroBias_phi = 0;
   trig_HLT_SinglePhoton10_Eta3p1ForPPRef_hltL1sZeroBias_et = 0;
   trig_HLT_SinglePhoton10_Eta3p1ForPPRef_hltHIPhoton10Eta3p1_eta = 0;
   trig_HLT_SinglePhoton10_Eta3p1ForPPRef_hltHIPhoton10Eta3p1_phi = 0;
   trig_HLT_SinglePhoton10_Eta3p1ForPPRef_hltHIPhoton10Eta3p1_et = 0;
   trig_HLT_SinglePhoton20_Eta3p1ForPPRef_hltL1sSingleEG10_eta = 0;
   trig_HLT_SinglePhoton20_Eta3p1ForPPRef_hltL1sSingleEG10_phi = 0;
   trig_HLT_SinglePhoton20_Eta3p1ForPPRef_hltL1sSingleEG10_et = 0;
   trig_HLT_SinglePhoton20_Eta3p1ForPPRef_hltHIPhoton20Eta3p1_eta = 0;
   trig_HLT_SinglePhoton20_Eta3p1ForPPRef_hltHIPhoton20Eta3p1_phi = 0;
   trig_HLT_SinglePhoton20_Eta3p1ForPPRef_hltHIPhoton20Eta3p1_et = 0;
   trig_HLT_SinglePhoton30_Eta3p1ForPPRef_hltL1sSingleEG10_eta = 0;
   trig_HLT_SinglePhoton30_Eta3p1ForPPRef_hltL1sSingleEG10_phi = 0;
   trig_HLT_SinglePhoton30_Eta3p1ForPPRef_hltL1sSingleEG10_et = 0;
   trig_HLT_SinglePhoton30_Eta3p1ForPPRef_hltHIPhoton30Eta3p1_eta = 0;
   trig_HLT_SinglePhoton30_Eta3p1ForPPRef_hltHIPhoton30Eta3p1_phi = 0;
   trig_HLT_SinglePhoton30_Eta3p1ForPPRef_hltHIPhoton30Eta3p1_et = 0;
   trig_HLT_Photon20_HoverELoose_hltL1sSingleEG10_eta = 0;
   trig_HLT_Photon20_HoverELoose_hltL1sSingleEG10_phi = 0;
   trig_HLT_Photon20_HoverELoose_hltL1sSingleEG10_et = 0;
   trig_HLT_Photon20_HoverELoose_hltEGL1SingleEG10Filter_eta = 0;
   trig_HLT_Photon20_HoverELoose_hltEGL1SingleEG10Filter_phi = 0;
   trig_HLT_Photon20_HoverELoose_hltEGL1SingleEG10Filter_et = 0;
   trig_HLT_Photon20_HoverELoose_hltEG20EtFilterLooseHoverE_eta = 0;
   trig_HLT_Photon20_HoverELoose_hltEG20EtFilterLooseHoverE_phi = 0;
   trig_HLT_Photon20_HoverELoose_hltEG20EtFilterLooseHoverE_et = 0;
   trig_HLT_Photon20_HoverELoose_hltEG20HEFilterLooseHoverE_eta = 0;
   trig_HLT_Photon20_HoverELoose_hltEG20HEFilterLooseHoverE_phi = 0;
   trig_HLT_Photon20_HoverELoose_hltEG20HEFilterLooseHoverE_et = 0;
   trig_HLT_Photon30_HoverELoose_hltL1sSingleEG10_eta = 0;
   trig_HLT_Photon30_HoverELoose_hltL1sSingleEG10_phi = 0;
   trig_HLT_Photon30_HoverELoose_hltL1sSingleEG10_et = 0;
   trig_HLT_Photon30_HoverELoose_hltEGL1SingleEG10Filter_eta = 0;
   trig_HLT_Photon30_HoverELoose_hltEGL1SingleEG10Filter_phi = 0;
   trig_HLT_Photon30_HoverELoose_hltEGL1SingleEG10Filter_et = 0;
   trig_HLT_Photon30_HoverELoose_hltEG30EtFilterLooseHoverE_eta = 0;
   trig_HLT_Photon30_HoverELoose_hltEG30EtFilterLooseHoverE_phi = 0;
   trig_HLT_Photon30_HoverELoose_hltEG30EtFilterLooseHoverE_et = 0;
   trig_HLT_Photon30_HoverELoose_hltEG30HEFilterLooseHoverE_eta = 0;
   trig_HLT_Photon30_HoverELoose_hltEG30HEFilterLooseHoverE_phi = 0;
   trig_HLT_Photon30_HoverELoose_hltEG30HEFilterLooseHoverE_et = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("ev_event", &ev_event, &b_ev_event);
   fChain->SetBranchAddress("ev_run", &ev_run, &b_ev_run);
   fChain->SetBranchAddress("ev_luminosityBlock", &ev_luminosityBlock, &b_ev_luminosityBlock);
   fChain->SetBranchAddress("ev_time", &ev_time, &b_ev_time);
   fChain->SetBranchAddress("ev_time_unixTime", &ev_time_unixTime, &b_ev_time_unixTime);
   fChain->SetBranchAddress("ev_time_microsecondOffset", &ev_time_microsecondOffset, &b_ev_time_microsecondOffset);
   fChain->SetBranchAddress("ev_fixedGridRhoAll", &ev_fixedGridRhoAll, &b_ev_fixedGridRhoAll);
   fChain->SetBranchAddress("ev_fixedGridRhoFastjetAll", &ev_fixedGridRhoFastjetAll, &b_ev_fixedGridRhoFastjetAll);
   fChain->SetBranchAddress("ev_fixedGridRhoFastjetAllCalo", &ev_fixedGridRhoFastjetAllCalo, &b_ev_fixedGridRhoFastjetAllCalo);
   fChain->SetBranchAddress("ev_fixedGridRhoFastjetCentralCalo", &ev_fixedGridRhoFastjetCentralCalo, &b_ev_fixedGridRhoFastjetCentralCalo);
   fChain->SetBranchAddress("ev_fixedGridRhoFastjetCentralChargedPileUp", &ev_fixedGridRhoFastjetCentralChargedPileUp, &b_ev_fixedGridRhoFastjetCentralChargedPileUp);
   fChain->SetBranchAddress("ev_fixedGridRhoFastjetCentralNeutral", &ev_fixedGridRhoFastjetCentralNeutral, &b_ev_fixedGridRhoFastjetCentralNeutral);
   fChain->SetBranchAddress("pv_n", &pv_n, &b_pv_n);
   fChain->SetBranchAddress("pv_x", &pv_x, &b_pv_x);
   fChain->SetBranchAddress("pv_y", &pv_y, &b_pv_y);
   fChain->SetBranchAddress("pv_z", &pv_z, &b_pv_z);
   fChain->SetBranchAddress("pv_ndof", &pv_ndof, &b_pv_ndof);
   fChain->SetBranchAddress("pv_normalizedChi2", &pv_normalizedChi2, &b_pv_normalizedChi2);
   fChain->SetBranchAddress("pv_isValid", &pv_isValid, &b_pv_isValid);
   fChain->SetBranchAddress("pv_isFake", &pv_isFake, &b_pv_isFake);
   fChain->SetBranchAddress("gsf_n", &gsf_n, &b_gsf_n);
   fChain->SetBranchAddress("gsf_classification", &gsf_classification, &b_gsf_classification);
   fChain->SetBranchAddress("gsf_ecalTrkEnergyPostCorr", &gsf_ecalTrkEnergyPostCorr, &b_gsf_ecalTrkEnergyPostCorr);
   fChain->SetBranchAddress("gsf_ecalEnergyPostCorr", &gsf_ecalEnergyPostCorr, &b_gsf_ecalEnergyPostCorr);
   fChain->SetBranchAddress("gsf_energy", &gsf_energy, &b_gsf_energy);
   fChain->SetBranchAddress("gsf_p", &gsf_p, &b_gsf_p);
   fChain->SetBranchAddress("gsf_pt", &gsf_pt, &b_gsf_pt);
   fChain->SetBranchAddress("gsf_et", &gsf_et, &b_gsf_et);
   fChain->SetBranchAddress("gsf_scE1x5", &gsf_scE1x5, &b_gsf_scE1x5);
   fChain->SetBranchAddress("gsf_scE5x5", &gsf_scE5x5, &b_gsf_scE5x5);
   fChain->SetBranchAddress("gsf_scE2x5Max", &gsf_scE2x5Max, &b_gsf_scE2x5Max);
   fChain->SetBranchAddress("gsf_full5x5_e5x5", &gsf_full5x5_e5x5, &b_gsf_full5x5_e5x5);
   fChain->SetBranchAddress("gsf_full5x5_e1x5", &gsf_full5x5_e1x5, &b_gsf_full5x5_e1x5);
   fChain->SetBranchAddress("gsf_full5x5_e2x5Max", &gsf_full5x5_e2x5Max, &b_gsf_full5x5_e2x5Max);
   fChain->SetBranchAddress("gsf_full5x5_sigmaIetaIeta", &gsf_full5x5_sigmaIetaIeta, &b_gsf_full5x5_sigmaIetaIeta);
   fChain->SetBranchAddress("gsf_full5x5_hcalOverEcal", &gsf_full5x5_hcalOverEcal, &b_gsf_full5x5_hcalOverEcal);
   fChain->SetBranchAddress("gsf_eta", &gsf_eta, &b_gsf_eta);
   fChain->SetBranchAddress("gsf_phi", &gsf_phi, &b_gsf_phi);
   fChain->SetBranchAddress("gsf_theta", &gsf_theta, &b_gsf_theta);
   fChain->SetBranchAddress("gsf_px", &gsf_px, &b_gsf_px);
   fChain->SetBranchAddress("gsf_py", &gsf_py, &b_gsf_py);
   fChain->SetBranchAddress("gsf_pz", &gsf_pz, &b_gsf_pz);
   fChain->SetBranchAddress("gsf_caloEnergy", &gsf_caloEnergy, &b_gsf_caloEnergy);
   fChain->SetBranchAddress("gsf_deltaEtaSuperClusterTrackAtVtx", &gsf_deltaEtaSuperClusterTrackAtVtx, &b_gsf_deltaEtaSuperClusterTrackAtVtx);
   fChain->SetBranchAddress("gsf_deltaPhiSuperClusterTrackAtVtx", &gsf_deltaPhiSuperClusterTrackAtVtx, &b_gsf_deltaPhiSuperClusterTrackAtVtx);
   fChain->SetBranchAddress("gsf_hadronicOverEm", &gsf_hadronicOverEm, &b_gsf_hadronicOverEm);
   fChain->SetBranchAddress("gsf_hcalDepth1OverEcal", &gsf_hcalDepth1OverEcal, &b_gsf_hcalDepth1OverEcal);
   fChain->SetBranchAddress("gsf_hcalDepth2OverEcal", &gsf_hcalDepth2OverEcal, &b_gsf_hcalDepth2OverEcal);
   fChain->SetBranchAddress("gsf_dr03TkSumPt", &gsf_dr03TkSumPt, &b_gsf_dr03TkSumPt);
   fChain->SetBranchAddress("gsf_heepTrkPtIso", &gsf_heepTrkPtIso, &b_gsf_heepTrkPtIso);
   fChain->SetBranchAddress("gsf_dr03EcalRecHitSumEt", &gsf_dr03EcalRecHitSumEt, &b_gsf_dr03EcalRecHitSumEt);
   fChain->SetBranchAddress("gsf_dr03HcalDepth1TowerSumEt", &gsf_dr03HcalDepth1TowerSumEt, &b_gsf_dr03HcalDepth1TowerSumEt);
   fChain->SetBranchAddress("gsf_dr03HcalDepth2TowerSumEt", &gsf_dr03HcalDepth2TowerSumEt, &b_gsf_dr03HcalDepth2TowerSumEt);
   fChain->SetBranchAddress("gsf_charge", &gsf_charge, &b_gsf_charge);
   fChain->SetBranchAddress("gsf_sigmaIetaIeta", &gsf_sigmaIetaIeta, &b_gsf_sigmaIetaIeta);
   fChain->SetBranchAddress("gsf_ecaldrivenSeed", &gsf_ecaldrivenSeed, &b_gsf_ecaldrivenSeed);
   fChain->SetBranchAddress("gsf_trackerdrivenSeed", &gsf_trackerdrivenSeed, &b_gsf_trackerdrivenSeed);
   fChain->SetBranchAddress("gsf_isEB", &gsf_isEB, &b_gsf_isEB);
   fChain->SetBranchAddress("gsf_isEE", &gsf_isEE, &b_gsf_isEE);
   fChain->SetBranchAddress("gsf_passConversionVeto", &gsf_passConversionVeto, &b_gsf_passConversionVeto);
   fChain->SetBranchAddress("gsf_VID_cutBasedElectronID_Fall17_94X_V1_loose", &gsf_VID_cutBasedElectronID_Fall17_94X_V1_loose, &b_gsf_VID_cutBasedElectronID_Fall17_94X_V1_loose);
   fChain->SetBranchAddress("gsf_VID_cutBasedElectronID_Fall17_94X_V1_medium", &gsf_VID_cutBasedElectronID_Fall17_94X_V1_medium, &b_gsf_VID_cutBasedElectronID_Fall17_94X_V1_medium);
   fChain->SetBranchAddress("gsf_VID_cutBasedElectronID_Fall17_94X_V1_tight", &gsf_VID_cutBasedElectronID_Fall17_94X_V1_tight, &b_gsf_VID_cutBasedElectronID_Fall17_94X_V1_tight);
   fChain->SetBranchAddress("gsf_VID_cutBasedElectronID_Fall17_94X_V1_veto", &gsf_VID_cutBasedElectronID_Fall17_94X_V1_veto, &b_gsf_VID_cutBasedElectronID_Fall17_94X_V1_veto);
   fChain->SetBranchAddress("gsf_VID_cutBasedElectronID_Summer16_80X_V1_loose", &gsf_VID_cutBasedElectronID_Summer16_80X_V1_loose, &b_gsf_VID_cutBasedElectronID_Summer16_80X_V1_loose);
   fChain->SetBranchAddress("gsf_VID_cutBasedElectronID_Summer16_80X_V1_medium", &gsf_VID_cutBasedElectronID_Summer16_80X_V1_medium, &b_gsf_VID_cutBasedElectronID_Summer16_80X_V1_medium);
   fChain->SetBranchAddress("gsf_VID_cutBasedElectronID_Summer16_80X_V1_tight", &gsf_VID_cutBasedElectronID_Summer16_80X_V1_tight, &b_gsf_VID_cutBasedElectronID_Summer16_80X_V1_tight);
   fChain->SetBranchAddress("gsf_VID_cutBasedElectronID_Summer16_80X_V1_veto", &gsf_VID_cutBasedElectronID_Summer16_80X_V1_veto, &b_gsf_VID_cutBasedElectronID_Summer16_80X_V1_veto);
   fChain->SetBranchAddress("gsf_VID_heepElectronID_HEEPV70", &gsf_VID_heepElectronID_HEEPV70, &b_gsf_VID_heepElectronID_HEEPV70);
   fChain->SetBranchAddress("gsf_VID_mvaEleID_Fall17_iso_V1_wp80", &gsf_VID_mvaEleID_Fall17_iso_V1_wp80, &b_gsf_VID_mvaEleID_Fall17_iso_V1_wp80);
   fChain->SetBranchAddress("gsf_VID_mvaEleID_Fall17_iso_V1_wp90", &gsf_VID_mvaEleID_Fall17_iso_V1_wp90, &b_gsf_VID_mvaEleID_Fall17_iso_V1_wp90);
   fChain->SetBranchAddress("gsf_VID_mvaEleID_Fall17_iso_V1_wpLoose", &gsf_VID_mvaEleID_Fall17_iso_V1_wpLoose, &b_gsf_VID_mvaEleID_Fall17_iso_V1_wpLoose);
   fChain->SetBranchAddress("gsf_VID_mvaEleID_Fall17_noIso_V1_wp80", &gsf_VID_mvaEleID_Fall17_noIso_V1_wp80, &b_gsf_VID_mvaEleID_Fall17_noIso_V1_wp80);
   fChain->SetBranchAddress("gsf_VID_mvaEleID_Fall17_noIso_V1_wp90", &gsf_VID_mvaEleID_Fall17_noIso_V1_wp90, &b_gsf_VID_mvaEleID_Fall17_noIso_V1_wp90);
   fChain->SetBranchAddress("gsf_VID_mvaEleID_Fall17_noIso_V1_wpLoose", &gsf_VID_mvaEleID_Fall17_noIso_V1_wpLoose, &b_gsf_VID_mvaEleID_Fall17_noIso_V1_wpLoose);
   fChain->SetBranchAddress("gsf_VID_mvaEleID_Spring16_GeneralPurpose_V1_wp80", &gsf_VID_mvaEleID_Spring16_GeneralPurpose_V1_wp80, &b_gsf_VID_mvaEleID_Spring16_GeneralPurpose_V1_wp80);
   fChain->SetBranchAddress("gsf_VID_mvaEleID_Spring16_GeneralPurpose_V1_wp90", &gsf_VID_mvaEleID_Spring16_GeneralPurpose_V1_wp90, &b_gsf_VID_mvaEleID_Spring16_GeneralPurpose_V1_wp90);
   fChain->SetBranchAddress("gsf_VID_mvaEleID_Spring16_HZZ_V1_wpLoose", &gsf_VID_mvaEleID_Spring16_HZZ_V1_wpLoose, &b_gsf_VID_mvaEleID_Spring16_HZZ_V1_wpLoose);
   fChain->SetBranchAddress("gsf_deltaEtaSeedClusterTrackAtCalo", &gsf_deltaEtaSeedClusterTrackAtCalo, &b_gsf_deltaEtaSeedClusterTrackAtCalo);
   fChain->SetBranchAddress("gsf_deltaPhiSeedClusterTrackAtCalo", &gsf_deltaPhiSeedClusterTrackAtCalo, &b_gsf_deltaPhiSeedClusterTrackAtCalo);
   fChain->SetBranchAddress("gsf_ecalEnergy", &gsf_ecalEnergy, &b_gsf_ecalEnergy);
   fChain->SetBranchAddress("gsf_eSuperClusterOverP", &gsf_eSuperClusterOverP, &b_gsf_eSuperClusterOverP);
   fChain->SetBranchAddress("gsf_dxy", &gsf_dxy, &b_gsf_dxy);
   fChain->SetBranchAddress("gsf_dxy_beamSpot", &gsf_dxy_beamSpot, &b_gsf_dxy_beamSpot);
   fChain->SetBranchAddress("gsf_dxy_firstPVtx", &gsf_dxy_firstPVtx, &b_gsf_dxy_firstPVtx);
   fChain->SetBranchAddress("gsf_dxyError", &gsf_dxyError, &b_gsf_dxyError);
   fChain->SetBranchAddress("gsf_dz", &gsf_dz, &b_gsf_dz);
   fChain->SetBranchAddress("gsf_dz_beamSpot", &gsf_dz_beamSpot, &b_gsf_dz_beamSpot);
   fChain->SetBranchAddress("gsf_dz_firstPVtx", &gsf_dz_firstPVtx, &b_gsf_dz_firstPVtx);
   fChain->SetBranchAddress("gsf_dzError", &gsf_dzError, &b_gsf_dzError);
   fChain->SetBranchAddress("gsf_vz", &gsf_vz, &b_gsf_vz);
   fChain->SetBranchAddress("gsf_numberOfValidHits", &gsf_numberOfValidHits, &b_gsf_numberOfValidHits);
   fChain->SetBranchAddress("gsf_nLostInnerHits", &gsf_nLostInnerHits, &b_gsf_nLostInnerHits);
   fChain->SetBranchAddress("gsf_nLostOuterHits", &gsf_nLostOuterHits, &b_gsf_nLostOuterHits);
   fChain->SetBranchAddress("gsf_convFlags", &gsf_convFlags, &b_gsf_convFlags);
   fChain->SetBranchAddress("gsf_convDist", &gsf_convDist, &b_gsf_convDist);
   fChain->SetBranchAddress("gsf_convDcot", &gsf_convDcot, &b_gsf_convDcot);
   fChain->SetBranchAddress("gsf_convRadius", &gsf_convRadius, &b_gsf_convRadius);
   fChain->SetBranchAddress("gsf_fBrem", &gsf_fBrem, &b_gsf_fBrem);
   fChain->SetBranchAddress("gsf_e1x5", &gsf_e1x5, &b_gsf_e1x5);
   fChain->SetBranchAddress("gsf_e2x5Max", &gsf_e2x5Max, &b_gsf_e2x5Max);
   fChain->SetBranchAddress("gsf_e5x5", &gsf_e5x5, &b_gsf_e5x5);
   fChain->SetBranchAddress("gsf_r9", &gsf_r9, &b_gsf_r9);
   fChain->SetBranchAddress("gsf_deltaEtaSeedClusterTrackAtVtx", &gsf_deltaEtaSeedClusterTrackAtVtx, &b_gsf_deltaEtaSeedClusterTrackAtVtx);
   fChain->SetBranchAddress("gsf_effArea", &gsf_effArea, &b_gsf_effArea);
   fChain->SetBranchAddress("gsf_sumChargedHadronPt", &gsf_sumChargedHadronPt, &b_gsf_sumChargedHadronPt);
   fChain->SetBranchAddress("gsf_sumNeutralHadronEt", &gsf_sumNeutralHadronEt, &b_gsf_sumNeutralHadronEt);
   fChain->SetBranchAddress("gsf_sumPhotonEt", &gsf_sumPhotonEt, &b_gsf_sumPhotonEt);
   fChain->SetBranchAddress("gsf_sumPUPt", &gsf_sumPUPt, &b_gsf_sumPUPt);
   fChain->SetBranchAddress("gsf_ooEmooP", &gsf_ooEmooP, &b_gsf_ooEmooP);
   fChain->SetBranchAddress("gsf_hitsinfo", &gsf_hitsinfo, &b_gsf_hitsinfo);
   fChain->SetBranchAddress("gsf_pixelMatch_dPhi1", &gsf_pixelMatch_dPhi1, &b_gsf_pixelMatch_dPhi1);
   fChain->SetBranchAddress("gsf_pixelMatch_dPhi2", &gsf_pixelMatch_dPhi2, &b_gsf_pixelMatch_dPhi2);
   fChain->SetBranchAddress("gsf_pixelMatch_dRz1", &gsf_pixelMatch_dRz1, &b_gsf_pixelMatch_dRz1);
   fChain->SetBranchAddress("gsf_pixelMatch_dRz2", &gsf_pixelMatch_dRz2, &b_gsf_pixelMatch_dRz2);
   fChain->SetBranchAddress("gsf_pixelMatch_subDetector1", &gsf_pixelMatch_subDetector1, &b_gsf_pixelMatch_subDetector1);
   fChain->SetBranchAddress("gsf_pixelMatch_subDetector2", &gsf_pixelMatch_subDetector2, &b_gsf_pixelMatch_subDetector2);
   fChain->SetBranchAddress("gsf_mc_bestDR", &gsf_mc_bestDR, &b_gsf_mc_bestDR);
   fChain->SetBranchAddress("gsf_mc_index", &gsf_mc_index, &b_gsf_mc_index);
   fChain->SetBranchAddress("gsf_mc_ERatio", &gsf_mc_ERatio, &b_gsf_mc_ERatio);
   fChain->SetBranchAddress("gsf_sc_energy", &gsf_sc_energy, &b_gsf_sc_energy);
   fChain->SetBranchAddress("gsf_sc_seed_eta", &gsf_sc_seed_eta, &b_gsf_sc_seed_eta);
   fChain->SetBranchAddress("gsf_sc_eta", &gsf_sc_eta, &b_gsf_sc_eta);
   fChain->SetBranchAddress("gsf_sc_etacorr", &gsf_sc_etacorr, &b_gsf_sc_etacorr);
   fChain->SetBranchAddress("gsf_sc_theta", &gsf_sc_theta, &b_gsf_sc_theta);
   fChain->SetBranchAddress("gsf_sc_thetacorr", &gsf_sc_thetacorr, &b_gsf_sc_thetacorr);
   fChain->SetBranchAddress("gsf_sc_et", &gsf_sc_et, &b_gsf_sc_et);
   fChain->SetBranchAddress("gsf_sc_phi", &gsf_sc_phi, &b_gsf_sc_phi);
   fChain->SetBranchAddress("gsf_sc_px", &gsf_sc_px, &b_gsf_sc_px);
   fChain->SetBranchAddress("gsf_sc_py", &gsf_sc_py, &b_gsf_sc_py);
   fChain->SetBranchAddress("gsf_sc_pz", &gsf_sc_pz, &b_gsf_sc_pz);
   fChain->SetBranchAddress("gsf_sc_x", &gsf_sc_x, &b_gsf_sc_x);
   fChain->SetBranchAddress("gsf_sc_y", &gsf_sc_y, &b_gsf_sc_y);
   fChain->SetBranchAddress("gsf_sc_z", &gsf_sc_z, &b_gsf_sc_z);
   fChain->SetBranchAddress("gsf_sc_phiWidth", &gsf_sc_phiWidth, &b_gsf_sc_phiWidth);
   fChain->SetBranchAddress("gsf_sc_etaWidth", &gsf_sc_etaWidth, &b_gsf_sc_etaWidth);
   fChain->SetBranchAddress("gsf_sc_seed_rawId", &gsf_sc_seed_rawId, &b_gsf_sc_seed_rawId);
   fChain->SetBranchAddress("gsf_sc_seed_ieta", &gsf_sc_seed_ieta, &b_gsf_sc_seed_ieta);
   fChain->SetBranchAddress("gsf_sc_seed_iphi", &gsf_sc_seed_iphi, &b_gsf_sc_seed_iphi);
   fChain->SetBranchAddress("gsf_sc_seed_kHasSwitchToGain6", &gsf_sc_seed_kHasSwitchToGain6, &b_gsf_sc_seed_kHasSwitchToGain6);
   fChain->SetBranchAddress("gsf_sc_seed_kHasSwitchToGain1", &gsf_sc_seed_kHasSwitchToGain1, &b_gsf_sc_seed_kHasSwitchToGain1);
   fChain->SetBranchAddress("gsf_swissCross", &gsf_swissCross, &b_gsf_swissCross);
   fChain->SetBranchAddress("gsf_sc_rawEnergy", &gsf_sc_rawEnergy, &b_gsf_sc_rawEnergy);
   fChain->SetBranchAddress("gsf_sc_preshowerEnergy", &gsf_sc_preshowerEnergy, &b_gsf_sc_preshowerEnergy);
   fChain->SetBranchAddress("gsf_sc_lazyTools_e2x5Right", &gsf_sc_lazyTools_e2x5Right, &b_gsf_sc_lazyTools_e2x5Right);
   fChain->SetBranchAddress("gsf_sc_lazyTools_e2x5Left", &gsf_sc_lazyTools_e2x5Left, &b_gsf_sc_lazyTools_e2x5Left);
   fChain->SetBranchAddress("gsf_sc_lazyTools_e2x5Top", &gsf_sc_lazyTools_e2x5Top, &b_gsf_sc_lazyTools_e2x5Top);
   fChain->SetBranchAddress("gsf_sc_lazyTools_e2x5Bottom", &gsf_sc_lazyTools_e2x5Bottom, &b_gsf_sc_lazyTools_e2x5Bottom);
   fChain->SetBranchAddress("gsf_sc_lazyTools_eMax", &gsf_sc_lazyTools_eMax, &b_gsf_sc_lazyTools_eMax);
   fChain->SetBranchAddress("gsf_sc_lazyTools_e2nd", &gsf_sc_lazyTools_e2nd, &b_gsf_sc_lazyTools_e2nd);
   fChain->SetBranchAddress("gsf_sc_lazyTools_eRight", &gsf_sc_lazyTools_eRight, &b_gsf_sc_lazyTools_eRight);
   fChain->SetBranchAddress("gsf_sc_lazyTools_eLeft", &gsf_sc_lazyTools_eLeft, &b_gsf_sc_lazyTools_eLeft);
   fChain->SetBranchAddress("gsf_sc_lazyTools_eTop", &gsf_sc_lazyTools_eTop, &b_gsf_sc_lazyTools_eTop);
   fChain->SetBranchAddress("gsf_sc_lazyTools_eBottom", &gsf_sc_lazyTools_eBottom, &b_gsf_sc_lazyTools_eBottom);
   fChain->SetBranchAddress("gsf_sc_lazyTools_e2x2", &gsf_sc_lazyTools_e2x2, &b_gsf_sc_lazyTools_e2x2);
   fChain->SetBranchAddress("gsf_sc_lazyTools_e3x3", &gsf_sc_lazyTools_e3x3, &b_gsf_sc_lazyTools_e3x3);
   fChain->SetBranchAddress("gsf_sc_lazyTools_e4x4", &gsf_sc_lazyTools_e4x4, &b_gsf_sc_lazyTools_e4x4);
   fChain->SetBranchAddress("gsf_sc_lazyTools_e5x5", &gsf_sc_lazyTools_e5x5, &b_gsf_sc_lazyTools_e5x5);
   fChain->SetBranchAddress("gsf_sc_lazyTools_e1x3", &gsf_sc_lazyTools_e1x3, &b_gsf_sc_lazyTools_e1x3);
   fChain->SetBranchAddress("gsf_sc_lazyTools_e3x1", &gsf_sc_lazyTools_e3x1, &b_gsf_sc_lazyTools_e3x1);
   fChain->SetBranchAddress("gsf_sc_lazyTools_e1x5", &gsf_sc_lazyTools_e1x5, &b_gsf_sc_lazyTools_e1x5);
   fChain->SetBranchAddress("gsf_sc_lazyTools_e5x1", &gsf_sc_lazyTools_e5x1, &b_gsf_sc_lazyTools_e5x1);
   fChain->SetBranchAddress("gsf_sc_lazyTools_eshitsixix", &gsf_sc_lazyTools_eshitsixix, &b_gsf_sc_lazyTools_eshitsixix);
   fChain->SetBranchAddress("gsf_sc_lazyTools_eshitsiyiy", &gsf_sc_lazyTools_eshitsiyiy, &b_gsf_sc_lazyTools_eshitsiyiy);
   fChain->SetBranchAddress("gsf_sc_lazyTools_eseffsixix", &gsf_sc_lazyTools_eseffsixix, &b_gsf_sc_lazyTools_eseffsixix);
   fChain->SetBranchAddress("gsf_sc_lazyTools_eseffsiyiy", &gsf_sc_lazyTools_eseffsiyiy, &b_gsf_sc_lazyTools_eseffsiyiy);
   fChain->SetBranchAddress("gsf_sc_lazyTools_eseffsirir", &gsf_sc_lazyTools_eseffsirir, &b_gsf_sc_lazyTools_eseffsirir);
   fChain->SetBranchAddress("gsf_sc_lazyTools_BasicClusterSeedTime", &gsf_sc_lazyTools_BasicClusterSeedTime, &b_gsf_sc_lazyTools_BasicClusterSeedTime);
   fChain->SetBranchAddress("EHits_isSaturated", &EHits_isSaturated, &b_EHits_isSaturated);
   fChain->SetBranchAddress("EBHits_rawId", &EBHits_rawId, &b_EBHits_rawId);
   fChain->SetBranchAddress("EBHits_iRechit", &EBHits_iRechit, &b_EBHits_iRechit);
   fChain->SetBranchAddress("EBHits_ieta", &EBHits_ieta, &b_EBHits_ieta);
   fChain->SetBranchAddress("EBHits_iphi", &EBHits_iphi, &b_EBHits_iphi);
   fChain->SetBranchAddress("EBHits_RecoFlag", &EBHits_RecoFlag, &b_EBHits_RecoFlag);
   fChain->SetBranchAddress("EBHits_kSaturated", &EBHits_kSaturated, &b_EBHits_kSaturated);
   fChain->SetBranchAddress("EBHits_kLeadingEdgeRecovered", &EBHits_kLeadingEdgeRecovered, &b_EBHits_kLeadingEdgeRecovered);
   fChain->SetBranchAddress("EBHits_kNeighboursRecovered", &EBHits_kNeighboursRecovered, &b_EBHits_kNeighboursRecovered);
   fChain->SetBranchAddress("EBHits_kWeird", &EBHits_kWeird, &b_EBHits_kWeird);
   fChain->SetBranchAddress("EBHits_energy", &EBHits_energy, &b_EBHits_energy);
   fChain->SetBranchAddress("EEHits_rawId", &EEHits_rawId, &b_EEHits_rawId);
   fChain->SetBranchAddress("EEHits_iRechit", &EEHits_iRechit, &b_EEHits_iRechit);
   fChain->SetBranchAddress("EEHits_ieta", &EEHits_ieta, &b_EEHits_ieta);
   fChain->SetBranchAddress("EEHits_iphi", &EEHits_iphi, &b_EEHits_iphi);
   fChain->SetBranchAddress("EEHits_RecoFlag", &EEHits_RecoFlag, &b_EEHits_RecoFlag);
   fChain->SetBranchAddress("EEHits_kSaturated", &EEHits_kSaturated, &b_EEHits_kSaturated);
   fChain->SetBranchAddress("EEHits_kLeadingEdgeRecovered", &EEHits_kLeadingEdgeRecovered, &b_EEHits_kLeadingEdgeRecovered);
   fChain->SetBranchAddress("EEHits_kNeighboursRecovered", &EEHits_kNeighboursRecovered, &b_EEHits_kNeighboursRecovered);
   fChain->SetBranchAddress("EEHits_kWeird", &EEHits_kWeird, &b_EEHits_kWeird);
   fChain->SetBranchAddress("EEHits_energy", &EEHits_energy, &b_EEHits_energy);
   fChain->SetBranchAddress("mu_n", &mu_n, &b_mu_n);
   fChain->SetBranchAddress("mu_gt_qoverp", &mu_gt_qoverp, &b_mu_gt_qoverp);
   fChain->SetBranchAddress("mu_gt_charge", &mu_gt_charge, &b_mu_gt_charge);
   fChain->SetBranchAddress("mu_gt_pt", &mu_gt_pt, &b_mu_gt_pt);
   fChain->SetBranchAddress("mu_gt_eta", &mu_gt_eta, &b_mu_gt_eta);
   fChain->SetBranchAddress("mu_gt_phi", &mu_gt_phi, &b_mu_gt_phi);
   fChain->SetBranchAddress("mu_gt_p", &mu_gt_p, &b_mu_gt_p);
   fChain->SetBranchAddress("mu_gt_px", &mu_gt_px, &b_mu_gt_px);
   fChain->SetBranchAddress("mu_gt_py", &mu_gt_py, &b_mu_gt_py);
   fChain->SetBranchAddress("mu_gt_pz", &mu_gt_pz, &b_mu_gt_pz);
   fChain->SetBranchAddress("mu_gt_theta", &mu_gt_theta, &b_mu_gt_theta);
   fChain->SetBranchAddress("mu_gt_lambda", &mu_gt_lambda, &b_mu_gt_lambda);
   fChain->SetBranchAddress("mu_gt_d0", &mu_gt_d0, &b_mu_gt_d0);
   fChain->SetBranchAddress("mu_gt_dz", &mu_gt_dz, &b_mu_gt_dz);
   fChain->SetBranchAddress("mu_gt_dz_beamspot", &mu_gt_dz_beamspot, &b_mu_gt_dz_beamspot);
   fChain->SetBranchAddress("mu_gt_dz_firstPVtx", &mu_gt_dz_firstPVtx, &b_mu_gt_dz_firstPVtx);
   fChain->SetBranchAddress("mu_gt_dxy", &mu_gt_dxy, &b_mu_gt_dxy);
   fChain->SetBranchAddress("mu_gt_dxy_beamspot", &mu_gt_dxy_beamspot, &b_mu_gt_dxy_beamspot);
   fChain->SetBranchAddress("mu_gt_dxy_firstPVtx", &mu_gt_dxy_firstPVtx, &b_mu_gt_dxy_firstPVtx);
   fChain->SetBranchAddress("mu_gt_dsz", &mu_gt_dsz, &b_mu_gt_dsz);
   fChain->SetBranchAddress("mu_gt_vx", &mu_gt_vx, &b_mu_gt_vx);
   fChain->SetBranchAddress("mu_gt_vy", &mu_gt_vy, &b_mu_gt_vy);
   fChain->SetBranchAddress("mu_gt_vz", &mu_gt_vz, &b_mu_gt_vz);
   fChain->SetBranchAddress("mu_gt_qoverpError", &mu_gt_qoverpError, &b_mu_gt_qoverpError);
   fChain->SetBranchAddress("mu_gt_ptError", &mu_gt_ptError, &b_mu_gt_ptError);
   fChain->SetBranchAddress("mu_gt_thetaError", &mu_gt_thetaError, &b_mu_gt_thetaError);
   fChain->SetBranchAddress("mu_gt_lambdaError", &mu_gt_lambdaError, &b_mu_gt_lambdaError);
   fChain->SetBranchAddress("mu_gt_phiError", &mu_gt_phiError, &b_mu_gt_phiError);
   fChain->SetBranchAddress("mu_gt_dxyError", &mu_gt_dxyError, &b_mu_gt_dxyError);
   fChain->SetBranchAddress("mu_gt_d0Error", &mu_gt_d0Error, &b_mu_gt_d0Error);
   fChain->SetBranchAddress("mu_gt_dszError", &mu_gt_dszError, &b_mu_gt_dszError);
   fChain->SetBranchAddress("mu_gt_dzError", &mu_gt_dzError, &b_mu_gt_dzError);
   fChain->SetBranchAddress("mu_gt_etaError", &mu_gt_etaError, &b_mu_gt_etaError);
   fChain->SetBranchAddress("mu_gt_chi2", &mu_gt_chi2, &b_mu_gt_chi2);
   fChain->SetBranchAddress("mu_gt_ndof", &mu_gt_ndof, &b_mu_gt_ndof);
   fChain->SetBranchAddress("mu_gt_normalizedChi2", &mu_gt_normalizedChi2, &b_mu_gt_normalizedChi2);
   fChain->SetBranchAddress("mu_ot_qoverp", &mu_ot_qoverp, &b_mu_ot_qoverp);
   fChain->SetBranchAddress("mu_ot_charge", &mu_ot_charge, &b_mu_ot_charge);
   fChain->SetBranchAddress("mu_ot_pt", &mu_ot_pt, &b_mu_ot_pt);
   fChain->SetBranchAddress("mu_ot_eta", &mu_ot_eta, &b_mu_ot_eta);
   fChain->SetBranchAddress("mu_ot_phi", &mu_ot_phi, &b_mu_ot_phi);
   fChain->SetBranchAddress("mu_ot_p", &mu_ot_p, &b_mu_ot_p);
   fChain->SetBranchAddress("mu_ot_px", &mu_ot_px, &b_mu_ot_px);
   fChain->SetBranchAddress("mu_ot_py", &mu_ot_py, &b_mu_ot_py);
   fChain->SetBranchAddress("mu_ot_pz", &mu_ot_pz, &b_mu_ot_pz);
   fChain->SetBranchAddress("mu_ot_theta", &mu_ot_theta, &b_mu_ot_theta);
   fChain->SetBranchAddress("mu_ot_lambda", &mu_ot_lambda, &b_mu_ot_lambda);
   fChain->SetBranchAddress("mu_ot_d0", &mu_ot_d0, &b_mu_ot_d0);
   fChain->SetBranchAddress("mu_ot_dz", &mu_ot_dz, &b_mu_ot_dz);
   fChain->SetBranchAddress("mu_ot_dz_beamspot", &mu_ot_dz_beamspot, &b_mu_ot_dz_beamspot);
   fChain->SetBranchAddress("mu_ot_dz_firstPVtx", &mu_ot_dz_firstPVtx, &b_mu_ot_dz_firstPVtx);
   fChain->SetBranchAddress("mu_ot_dxy", &mu_ot_dxy, &b_mu_ot_dxy);
   fChain->SetBranchAddress("mu_ot_dxy_beamspot", &mu_ot_dxy_beamspot, &b_mu_ot_dxy_beamspot);
   fChain->SetBranchAddress("mu_ot_dxy_firstPVtx", &mu_ot_dxy_firstPVtx, &b_mu_ot_dxy_firstPVtx);
   fChain->SetBranchAddress("mu_ot_dsz", &mu_ot_dsz, &b_mu_ot_dsz);
   fChain->SetBranchAddress("mu_ot_vx", &mu_ot_vx, &b_mu_ot_vx);
   fChain->SetBranchAddress("mu_ot_vy", &mu_ot_vy, &b_mu_ot_vy);
   fChain->SetBranchAddress("mu_ot_vz", &mu_ot_vz, &b_mu_ot_vz);
   fChain->SetBranchAddress("mu_ot_qoverpError", &mu_ot_qoverpError, &b_mu_ot_qoverpError);
   fChain->SetBranchAddress("mu_ot_ptError", &mu_ot_ptError, &b_mu_ot_ptError);
   fChain->SetBranchAddress("mu_ot_thetaError", &mu_ot_thetaError, &b_mu_ot_thetaError);
   fChain->SetBranchAddress("mu_ot_lambdaError", &mu_ot_lambdaError, &b_mu_ot_lambdaError);
   fChain->SetBranchAddress("mu_ot_phiError", &mu_ot_phiError, &b_mu_ot_phiError);
   fChain->SetBranchAddress("mu_ot_dxyError", &mu_ot_dxyError, &b_mu_ot_dxyError);
   fChain->SetBranchAddress("mu_ot_d0Error", &mu_ot_d0Error, &b_mu_ot_d0Error);
   fChain->SetBranchAddress("mu_ot_dszError", &mu_ot_dszError, &b_mu_ot_dszError);
   fChain->SetBranchAddress("mu_ot_dzError", &mu_ot_dzError, &b_mu_ot_dzError);
   fChain->SetBranchAddress("mu_ot_etaError", &mu_ot_etaError, &b_mu_ot_etaError);
   fChain->SetBranchAddress("mu_ot_chi2", &mu_ot_chi2, &b_mu_ot_chi2);
   fChain->SetBranchAddress("mu_ot_ndof", &mu_ot_ndof, &b_mu_ot_ndof);
   fChain->SetBranchAddress("mu_ot_normalizedChi2", &mu_ot_normalizedChi2, &b_mu_ot_normalizedChi2);
   fChain->SetBranchAddress("mu_it_qoverp", &mu_it_qoverp, &b_mu_it_qoverp);
   fChain->SetBranchAddress("mu_it_charge", &mu_it_charge, &b_mu_it_charge);
   fChain->SetBranchAddress("mu_it_pt", &mu_it_pt, &b_mu_it_pt);
   fChain->SetBranchAddress("mu_it_eta", &mu_it_eta, &b_mu_it_eta);
   fChain->SetBranchAddress("mu_it_phi", &mu_it_phi, &b_mu_it_phi);
   fChain->SetBranchAddress("mu_it_p", &mu_it_p, &b_mu_it_p);
   fChain->SetBranchAddress("mu_it_px", &mu_it_px, &b_mu_it_px);
   fChain->SetBranchAddress("mu_it_py", &mu_it_py, &b_mu_it_py);
   fChain->SetBranchAddress("mu_it_pz", &mu_it_pz, &b_mu_it_pz);
   fChain->SetBranchAddress("mu_it_theta", &mu_it_theta, &b_mu_it_theta);
   fChain->SetBranchAddress("mu_it_lambda", &mu_it_lambda, &b_mu_it_lambda);
   fChain->SetBranchAddress("mu_it_d0", &mu_it_d0, &b_mu_it_d0);
   fChain->SetBranchAddress("mu_it_dz", &mu_it_dz, &b_mu_it_dz);
   fChain->SetBranchAddress("mu_it_dz_beamspot", &mu_it_dz_beamspot, &b_mu_it_dz_beamspot);
   fChain->SetBranchAddress("mu_it_dz_firstPVtx", &mu_it_dz_firstPVtx, &b_mu_it_dz_firstPVtx);
   fChain->SetBranchAddress("mu_it_dxy", &mu_it_dxy, &b_mu_it_dxy);
   fChain->SetBranchAddress("mu_it_dxy_beamspot", &mu_it_dxy_beamspot, &b_mu_it_dxy_beamspot);
   fChain->SetBranchAddress("mu_it_dxy_firstPVtx", &mu_it_dxy_firstPVtx, &b_mu_it_dxy_firstPVtx);
   fChain->SetBranchAddress("mu_it_dsz", &mu_it_dsz, &b_mu_it_dsz);
   fChain->SetBranchAddress("mu_it_vx", &mu_it_vx, &b_mu_it_vx);
   fChain->SetBranchAddress("mu_it_vy", &mu_it_vy, &b_mu_it_vy);
   fChain->SetBranchAddress("mu_it_vz", &mu_it_vz, &b_mu_it_vz);
   fChain->SetBranchAddress("mu_it_qoverpError", &mu_it_qoverpError, &b_mu_it_qoverpError);
   fChain->SetBranchAddress("mu_it_ptError", &mu_it_ptError, &b_mu_it_ptError);
   fChain->SetBranchAddress("mu_it_thetaError", &mu_it_thetaError, &b_mu_it_thetaError);
   fChain->SetBranchAddress("mu_it_lambdaError", &mu_it_lambdaError, &b_mu_it_lambdaError);
   fChain->SetBranchAddress("mu_it_phiError", &mu_it_phiError, &b_mu_it_phiError);
   fChain->SetBranchAddress("mu_it_dxyError", &mu_it_dxyError, &b_mu_it_dxyError);
   fChain->SetBranchAddress("mu_it_d0Error", &mu_it_d0Error, &b_mu_it_d0Error);
   fChain->SetBranchAddress("mu_it_dszError", &mu_it_dszError, &b_mu_it_dszError);
   fChain->SetBranchAddress("mu_it_dzError", &mu_it_dzError, &b_mu_it_dzError);
   fChain->SetBranchAddress("mu_it_etaError", &mu_it_etaError, &b_mu_it_etaError);
   fChain->SetBranchAddress("mu_it_chi2", &mu_it_chi2, &b_mu_it_chi2);
   fChain->SetBranchAddress("mu_it_ndof", &mu_it_ndof, &b_mu_it_ndof);
   fChain->SetBranchAddress("mu_it_normalizedChi2", &mu_it_normalizedChi2, &b_mu_it_normalizedChi2);
   fChain->SetBranchAddress("mu_ibt_qoverp", &mu_ibt_qoverp, &b_mu_ibt_qoverp);
   fChain->SetBranchAddress("mu_ibt_charge", &mu_ibt_charge, &b_mu_ibt_charge);
   fChain->SetBranchAddress("mu_ibt_pt", &mu_ibt_pt, &b_mu_ibt_pt);
   fChain->SetBranchAddress("mu_ibt_eta", &mu_ibt_eta, &b_mu_ibt_eta);
   fChain->SetBranchAddress("mu_ibt_phi", &mu_ibt_phi, &b_mu_ibt_phi);
   fChain->SetBranchAddress("mu_ibt_p", &mu_ibt_p, &b_mu_ibt_p);
   fChain->SetBranchAddress("mu_ibt_px", &mu_ibt_px, &b_mu_ibt_px);
   fChain->SetBranchAddress("mu_ibt_py", &mu_ibt_py, &b_mu_ibt_py);
   fChain->SetBranchAddress("mu_ibt_pz", &mu_ibt_pz, &b_mu_ibt_pz);
   fChain->SetBranchAddress("mu_ibt_theta", &mu_ibt_theta, &b_mu_ibt_theta);
   fChain->SetBranchAddress("mu_ibt_lambda", &mu_ibt_lambda, &b_mu_ibt_lambda);
   fChain->SetBranchAddress("mu_ibt_d0", &mu_ibt_d0, &b_mu_ibt_d0);
   fChain->SetBranchAddress("mu_ibt_dz", &mu_ibt_dz, &b_mu_ibt_dz);
   fChain->SetBranchAddress("mu_ibt_dz_beamspot", &mu_ibt_dz_beamspot, &b_mu_ibt_dz_beamspot);
   fChain->SetBranchAddress("mu_ibt_dz_firstPVtx", &mu_ibt_dz_firstPVtx, &b_mu_ibt_dz_firstPVtx);
   fChain->SetBranchAddress("mu_ibt_dxy", &mu_ibt_dxy, &b_mu_ibt_dxy);
   fChain->SetBranchAddress("mu_ibt_dxy_beamspot", &mu_ibt_dxy_beamspot, &b_mu_ibt_dxy_beamspot);
   fChain->SetBranchAddress("mu_ibt_dxy_firstPVtx", &mu_ibt_dxy_firstPVtx, &b_mu_ibt_dxy_firstPVtx);
   fChain->SetBranchAddress("mu_ibt_dsz", &mu_ibt_dsz, &b_mu_ibt_dsz);
   fChain->SetBranchAddress("mu_ibt_vx", &mu_ibt_vx, &b_mu_ibt_vx);
   fChain->SetBranchAddress("mu_ibt_vy", &mu_ibt_vy, &b_mu_ibt_vy);
   fChain->SetBranchAddress("mu_ibt_vz", &mu_ibt_vz, &b_mu_ibt_vz);
   fChain->SetBranchAddress("mu_ibt_qoverpError", &mu_ibt_qoverpError, &b_mu_ibt_qoverpError);
   fChain->SetBranchAddress("mu_ibt_ptError", &mu_ibt_ptError, &b_mu_ibt_ptError);
   fChain->SetBranchAddress("mu_ibt_thetaError", &mu_ibt_thetaError, &b_mu_ibt_thetaError);
   fChain->SetBranchAddress("mu_ibt_lambdaError", &mu_ibt_lambdaError, &b_mu_ibt_lambdaError);
   fChain->SetBranchAddress("mu_ibt_phiError", &mu_ibt_phiError, &b_mu_ibt_phiError);
   fChain->SetBranchAddress("mu_ibt_dxyError", &mu_ibt_dxyError, &b_mu_ibt_dxyError);
   fChain->SetBranchAddress("mu_ibt_d0Error", &mu_ibt_d0Error, &b_mu_ibt_d0Error);
   fChain->SetBranchAddress("mu_ibt_dszError", &mu_ibt_dszError, &b_mu_ibt_dszError);
   fChain->SetBranchAddress("mu_ibt_dzError", &mu_ibt_dzError, &b_mu_ibt_dzError);
   fChain->SetBranchAddress("mu_ibt_etaError", &mu_ibt_etaError, &b_mu_ibt_etaError);
   fChain->SetBranchAddress("mu_ibt_chi2", &mu_ibt_chi2, &b_mu_ibt_chi2);
   fChain->SetBranchAddress("mu_ibt_ndof", &mu_ibt_ndof, &b_mu_ibt_ndof);
   fChain->SetBranchAddress("mu_ibt_normalizedChi2", &mu_ibt_normalizedChi2, &b_mu_ibt_normalizedChi2);
   fChain->SetBranchAddress("mu_isGlobalMuon", &mu_isGlobalMuon, &b_mu_isGlobalMuon);
   fChain->SetBranchAddress("mu_isStandAloneMuon", &mu_isStandAloneMuon, &b_mu_isStandAloneMuon);
   fChain->SetBranchAddress("mu_isTrackerMuon", &mu_isTrackerMuon, &b_mu_isTrackerMuon);
   fChain->SetBranchAddress("mu_isPFMuon", &mu_isPFMuon, &b_mu_isPFMuon);
   fChain->SetBranchAddress("mu_isPFIsolationValid", &mu_isPFIsolationValid, &b_mu_isPFIsolationValid);
   fChain->SetBranchAddress("mu_isGoodMuonTMLastStationLoose", &mu_isGoodMuonTMLastStationLoose, &b_mu_isGoodMuonTMLastStationLoose);
   fChain->SetBranchAddress("mu_isGoodMuonTMLastStationTight", &mu_isGoodMuonTMLastStationTight, &b_mu_isGoodMuonTMLastStationTight);
   fChain->SetBranchAddress("mu_isGoodMuonTM2DCompatibilityLoose", &mu_isGoodMuonTM2DCompatibilityLoose, &b_mu_isGoodMuonTM2DCompatibilityLoose);
   fChain->SetBranchAddress("mu_isGoodMuonTM2DCompatibilityTight", &mu_isGoodMuonTM2DCompatibilityTight, &b_mu_isGoodMuonTM2DCompatibilityTight);
   fChain->SetBranchAddress("mu_isGoodMuonTMOneStationLoose", &mu_isGoodMuonTMOneStationLoose, &b_mu_isGoodMuonTMOneStationLoose);
   fChain->SetBranchAddress("mu_isGoodMuonTMOneStationTight", &mu_isGoodMuonTMOneStationTight, &b_mu_isGoodMuonTMOneStationTight);
   fChain->SetBranchAddress("mu_isGoodMuonTMLastStationOptimizedLowPtLoose", &mu_isGoodMuonTMLastStationOptimizedLowPtLoose, &b_mu_isGoodMuonTMLastStationOptimizedLowPtLoose);
   fChain->SetBranchAddress("mu_isGoodMuonTMLastStationOptimizedLowPtTight", &mu_isGoodMuonTMLastStationOptimizedLowPtTight, &b_mu_isGoodMuonTMLastStationOptimizedLowPtTight);
   fChain->SetBranchAddress("mu_isTightMuon", &mu_isTightMuon, &b_mu_isTightMuon);
   fChain->SetBranchAddress("mu_isMediumMuon", &mu_isMediumMuon, &b_mu_isMediumMuon);
   fChain->SetBranchAddress("mu_isLooseMuon", &mu_isLooseMuon, &b_mu_isLooseMuon);
   fChain->SetBranchAddress("mu_isSoftMuon", &mu_isSoftMuon, &b_mu_isSoftMuon);
   fChain->SetBranchAddress("mu_isHighPtMuonOld", &mu_isHighPtMuonOld, &b_mu_isHighPtMuonOld);
   fChain->SetBranchAddress("mu_isHighPtMuon", &mu_isHighPtMuon, &b_mu_isHighPtMuon);
   fChain->SetBranchAddress("mu_isTrackerHighPtMuon", &mu_isTrackerHighPtMuon, &b_mu_isTrackerHighPtMuon);
   fChain->SetBranchAddress("mu_CutBasedIdLoose", &mu_CutBasedIdLoose, &b_mu_CutBasedIdLoose);
   fChain->SetBranchAddress("mu_CutBasedIdMedium", &mu_CutBasedIdMedium, &b_mu_CutBasedIdMedium);
   fChain->SetBranchAddress("mu_CutBasedIdMediumPrompt", &mu_CutBasedIdMediumPrompt, &b_mu_CutBasedIdMediumPrompt);
   fChain->SetBranchAddress("mu_CutBasedIdTight", &mu_CutBasedIdTight, &b_mu_CutBasedIdTight);
   fChain->SetBranchAddress("mu_CutBasedIdGlobalHighPt", &mu_CutBasedIdGlobalHighPt, &b_mu_CutBasedIdGlobalHighPt);
   fChain->SetBranchAddress("mu_CutBasedIdTrkHighPt", &mu_CutBasedIdTrkHighPt, &b_mu_CutBasedIdTrkHighPt);
   fChain->SetBranchAddress("mu_PFIsoVeryLoose", &mu_PFIsoVeryLoose, &b_mu_PFIsoVeryLoose);
   fChain->SetBranchAddress("mu_PFIsoLoose", &mu_PFIsoLoose, &b_mu_PFIsoLoose);
   fChain->SetBranchAddress("mu_PFIsoMedium", &mu_PFIsoMedium, &b_mu_PFIsoMedium);
   fChain->SetBranchAddress("mu_PFIsoTight", &mu_PFIsoTight, &b_mu_PFIsoTight);
   fChain->SetBranchAddress("mu_PFIsoVeryTight", &mu_PFIsoVeryTight, &b_mu_PFIsoVeryTight);
   fChain->SetBranchAddress("mu_PFIsoVeryVeryTight", &mu_PFIsoVeryVeryTight, &b_mu_PFIsoVeryVeryTight);
   fChain->SetBranchAddress("mu_TkIsoLoose", &mu_TkIsoLoose, &b_mu_TkIsoLoose);
   fChain->SetBranchAddress("mu_TkIsoTight", &mu_TkIsoTight, &b_mu_TkIsoTight);
   fChain->SetBranchAddress("mu_SoftCutBasedId", &mu_SoftCutBasedId, &b_mu_SoftCutBasedId);
   fChain->SetBranchAddress("mu_SoftMvaId", &mu_SoftMvaId, &b_mu_SoftMvaId);
   fChain->SetBranchAddress("mu_MvaLoose", &mu_MvaLoose, &b_mu_MvaLoose);
   fChain->SetBranchAddress("mu_MvaMedium", &mu_MvaMedium, &b_mu_MvaMedium);
   fChain->SetBranchAddress("mu_MvaTight", &mu_MvaTight, &b_mu_MvaTight);
   fChain->SetBranchAddress("mu_MiniIsoLoose", &mu_MiniIsoLoose, &b_mu_MiniIsoLoose);
   fChain->SetBranchAddress("mu_MiniIsoMedium", &mu_MiniIsoMedium, &b_mu_MiniIsoMedium);
   fChain->SetBranchAddress("mu_MiniIsoTight", &mu_MiniIsoTight, &b_mu_MiniIsoTight);
   fChain->SetBranchAddress("mu_MiniIsoVeryTight", &mu_MiniIsoVeryTight, &b_mu_MiniIsoVeryTight);
   fChain->SetBranchAddress("mu_TriggerIdLoose", &mu_TriggerIdLoose, &b_mu_TriggerIdLoose);
   fChain->SetBranchAddress("mu_InTimeMuon", &mu_InTimeMuon, &b_mu_InTimeMuon);
   fChain->SetBranchAddress("mu_MultiIsoLoose", &mu_MultiIsoLoose, &b_mu_MultiIsoLoose);
   fChain->SetBranchAddress("mu_MultiIsoMedium", &mu_MultiIsoMedium, &b_mu_MultiIsoMedium);
   fChain->SetBranchAddress("mu_numberOfMatchedStations", &mu_numberOfMatchedStations, &b_mu_numberOfMatchedStations);
   fChain->SetBranchAddress("mu_numberOfValidPixelHits", &mu_numberOfValidPixelHits, &b_mu_numberOfValidPixelHits);
   fChain->SetBranchAddress("mu_trackerLayersWithMeasurement", &mu_trackerLayersWithMeasurement, &b_mu_trackerLayersWithMeasurement);
   fChain->SetBranchAddress("mu_numberOfValidMuonHits", &mu_numberOfValidMuonHits, &b_mu_numberOfValidMuonHits);
   fChain->SetBranchAddress("mu_pixelLayersWithMeasurement", &mu_pixelLayersWithMeasurement, &b_mu_pixelLayersWithMeasurement);
   fChain->SetBranchAddress("mu_innerTrack_validFraction", &mu_innerTrack_validFraction, &b_mu_innerTrack_validFraction);
   fChain->SetBranchAddress("mu_combinedQuality_trkKink", &mu_combinedQuality_trkKink, &b_mu_combinedQuality_trkKink);
   fChain->SetBranchAddress("mu_combinedQuality_chi2LocalPosition", &mu_combinedQuality_chi2LocalPosition, &b_mu_combinedQuality_chi2LocalPosition);
   fChain->SetBranchAddress("mu_segmentCompatibility", &mu_segmentCompatibility, &b_mu_segmentCompatibility);
   fChain->SetBranchAddress("mu_dB", &mu_dB, &b_mu_dB);
   fChain->SetBranchAddress("mu_pt_default", &mu_pt_default, &b_mu_pt_default);
   fChain->SetBranchAddress("mu_isolationR03_sumPt", &mu_isolationR03_sumPt, &b_mu_isolationR03_sumPt);
   fChain->SetBranchAddress("mu_isolationR03_trackerVetoPt", &mu_isolationR03_trackerVetoPt, &b_mu_isolationR03_trackerVetoPt);
   fChain->SetBranchAddress("mu_isolationR03_emEt", &mu_isolationR03_emEt, &b_mu_isolationR03_emEt);
   fChain->SetBranchAddress("mu_isolationR03_emVetoEt", &mu_isolationR03_emVetoEt, &b_mu_isolationR03_emVetoEt);
   fChain->SetBranchAddress("mu_isolationR03_hadEt", &mu_isolationR03_hadEt, &b_mu_isolationR03_hadEt);
   fChain->SetBranchAddress("mu_isolationR03_hadVetoEt", &mu_isolationR03_hadVetoEt, &b_mu_isolationR03_hadVetoEt);
   fChain->SetBranchAddress("mu_isolationR05_sumPt", &mu_isolationR05_sumPt, &b_mu_isolationR05_sumPt);
   fChain->SetBranchAddress("mu_isolationR05_trackerVetoPt", &mu_isolationR05_trackerVetoPt, &b_mu_isolationR05_trackerVetoPt);
   fChain->SetBranchAddress("mu_isolationR05_emEt", &mu_isolationR05_emEt, &b_mu_isolationR05_emEt);
   fChain->SetBranchAddress("mu_isolationR05_emVetoEt", &mu_isolationR05_emVetoEt, &b_mu_isolationR05_emVetoEt);
   fChain->SetBranchAddress("mu_isolationR05_hadEt", &mu_isolationR05_hadEt, &b_mu_isolationR05_hadEt);
   fChain->SetBranchAddress("mu_isolationR05_hadVetoEt", &mu_isolationR05_hadVetoEt, &b_mu_isolationR05_hadVetoEt);
   fChain->SetBranchAddress("mu_pfIsolationR03_sumChargedHadronPt", &mu_pfIsolationR03_sumChargedHadronPt, &b_mu_pfIsolationR03_sumChargedHadronPt);
   fChain->SetBranchAddress("mu_pfIsolationR03_sumNeutralHadronEt", &mu_pfIsolationR03_sumNeutralHadronEt, &b_mu_pfIsolationR03_sumNeutralHadronEt);
   fChain->SetBranchAddress("mu_pfIsolationR03_sumChargedParticlePt", &mu_pfIsolationR03_sumChargedParticlePt, &b_mu_pfIsolationR03_sumChargedParticlePt);
   fChain->SetBranchAddress("mu_pfIsolationR03_sumPhotonEt", &mu_pfIsolationR03_sumPhotonEt, &b_mu_pfIsolationR03_sumPhotonEt);
   fChain->SetBranchAddress("mu_pfIsolationR03_sumNeutralHadronEtHighThreshold", &mu_pfIsolationR03_sumNeutralHadronEtHighThreshold, &b_mu_pfIsolationR03_sumNeutralHadronEtHighThreshold);
   fChain->SetBranchAddress("mu_pfIsolationR03_sumPhotonEtHighThreshold", &mu_pfIsolationR03_sumPhotonEtHighThreshold, &b_mu_pfIsolationR03_sumPhotonEtHighThreshold);
   fChain->SetBranchAddress("mu_pfIsolationR03_sumPUPt", &mu_pfIsolationR03_sumPUPt, &b_mu_pfIsolationR03_sumPUPt);
   fChain->SetBranchAddress("mu_pfIsolationR04_sumChargedHadronPt", &mu_pfIsolationR04_sumChargedHadronPt, &b_mu_pfIsolationR04_sumChargedHadronPt);
   fChain->SetBranchAddress("mu_pfIsolationR04_sumNeutralHadronEt", &mu_pfIsolationR04_sumNeutralHadronEt, &b_mu_pfIsolationR04_sumNeutralHadronEt);
   fChain->SetBranchAddress("mu_pfIsolationR04_sumChargedParticlePt", &mu_pfIsolationR04_sumChargedParticlePt, &b_mu_pfIsolationR04_sumChargedParticlePt);
   fChain->SetBranchAddress("mu_pfIsolationR04_sumPhotonEt", &mu_pfIsolationR04_sumPhotonEt, &b_mu_pfIsolationR04_sumPhotonEt);
   fChain->SetBranchAddress("mu_pfIsolationR04_sumNeutralHadronEtHighThreshold", &mu_pfIsolationR04_sumNeutralHadronEtHighThreshold, &b_mu_pfIsolationR04_sumNeutralHadronEtHighThreshold);
   fChain->SetBranchAddress("mu_pfIsolationR04_sumPhotonEtHighThreshold", &mu_pfIsolationR04_sumPhotonEtHighThreshold, &b_mu_pfIsolationR04_sumPhotonEtHighThreshold);
   fChain->SetBranchAddress("mu_pfIsolationR04_sumPUPt", &mu_pfIsolationR04_sumPUPt, &b_mu_pfIsolationR04_sumPUPt);
   fChain->SetBranchAddress("mu_pfIsoDbCorrected03", &mu_pfIsoDbCorrected03, &b_mu_pfIsoDbCorrected03);
   fChain->SetBranchAddress("mu_pfIsoDbCorrected04", &mu_pfIsoDbCorrected04, &b_mu_pfIsoDbCorrected04);
   fChain->SetBranchAddress("mu_isoTrackerBased03", &mu_isoTrackerBased03, &b_mu_isoTrackerBased03);
   fChain->SetBranchAddress("jet_n", &jet_n, &b_jet_n);
   fChain->SetBranchAddress("jet_px", &jet_px, &b_jet_px);
   fChain->SetBranchAddress("jet_py", &jet_py, &b_jet_py);
   fChain->SetBranchAddress("jet_pz", &jet_pz, &b_jet_pz);
   fChain->SetBranchAddress("jet_pt", &jet_pt, &b_jet_pt);
   fChain->SetBranchAddress("jet_eta", &jet_eta, &b_jet_eta);
   fChain->SetBranchAddress("jet_theta", &jet_theta, &b_jet_theta);
   fChain->SetBranchAddress("jet_phi", &jet_phi, &b_jet_phi);
   fChain->SetBranchAddress("jet_energy", &jet_energy, &b_jet_energy);
   fChain->SetBranchAddress("jet_mass", &jet_mass, &b_jet_mass);
   fChain->SetBranchAddress("jet_chargedEmEnergyFraction", &jet_chargedEmEnergyFraction, &b_jet_chargedEmEnergyFraction);
   fChain->SetBranchAddress("jet_neutralHadronEnergyFraction", &jet_neutralHadronEnergyFraction, &b_jet_neutralHadronEnergyFraction);
   fChain->SetBranchAddress("jet_neutralEmEnergyFraction", &jet_neutralEmEnergyFraction, &b_jet_neutralEmEnergyFraction);
   fChain->SetBranchAddress("jet_chargedHadronEnergyFraction", &jet_chargedHadronEnergyFraction, &b_jet_chargedHadronEnergyFraction);
   fChain->SetBranchAddress("jet_muonEnergyFraction", &jet_muonEnergyFraction, &b_jet_muonEnergyFraction);
   fChain->SetBranchAddress("jet_chargedMultiplicity", &jet_chargedMultiplicity, &b_jet_chargedMultiplicity);
   fChain->SetBranchAddress("jet_neutralMultiplicity", &jet_neutralMultiplicity, &b_jet_neutralMultiplicity);
   fChain->SetBranchAddress("jet_partonFlavour", &jet_partonFlavour, &b_jet_partonFlavour);
   fChain->SetBranchAddress("jet_hadronFlavour", &jet_hadronFlavour, &b_jet_hadronFlavour);
   fChain->SetBranchAddress("jet_CSVv2", &jet_CSVv2, &b_jet_CSVv2);
   fChain->SetBranchAddress("jet_CvsL", &jet_CvsL, &b_jet_CvsL);
   fChain->SetBranchAddress("jet_CvsB", &jet_CvsB, &b_jet_CvsB);
   fChain->SetBranchAddress("jet_MVA2BJets", &jet_MVA2BJets, &b_jet_MVA2BJets);
   fChain->SetBranchAddress("jet_CvsB_DeepJet_charm_tagger", &jet_CvsB_DeepJet_charm_tagger, &b_jet_CvsB_DeepJet_charm_tagger);
   fChain->SetBranchAddress("jet_CvsL_DeepJet_charm_tagger", &jet_CvsL_DeepJet_charm_tagger, &b_jet_CvsL_DeepJet_charm_tagger);
   fChain->SetBranchAddress("jet_CvsB_DeepCSV_charm_tagger", &jet_CvsB_DeepCSV_charm_tagger, &b_jet_CvsB_DeepCSV_charm_tagger);
   fChain->SetBranchAddress("jet_CvsL_DeepCSV_charm_tagger", &jet_CvsL_DeepCSV_charm_tagger, &b_jet_CvsL_DeepCSV_charm_tagger);
   fChain->SetBranchAddress("jet_DeepJet", &jet_DeepJet, &b_jet_DeepJet);
   fChain->SetBranchAddress("jet_DeepCSV", &jet_DeepCSV, &b_jet_DeepCSV);
   fChain->SetBranchAddress("jet_isJetIDLoose_2016", &jet_isJetIDLoose_2016, &b_jet_isJetIDLoose_2016);
   fChain->SetBranchAddress("jet_isJetIDTight_2016", &jet_isJetIDTight_2016, &b_jet_isJetIDTight_2016);
   fChain->SetBranchAddress("jet_isJetIDTightLepVeto_2016", &jet_isJetIDTightLepVeto_2016, &b_jet_isJetIDTightLepVeto_2016);
   fChain->SetBranchAddress("jet_isJetID_2017", &jet_isJetID_2017, &b_jet_isJetID_2017);
   fChain->SetBranchAddress("jet_isJetIDLepVeto_2017", &jet_isJetIDLepVeto_2017, &b_jet_isJetIDLepVeto_2017);
   fChain->SetBranchAddress("jet_isJetID_2018", &jet_isJetID_2018, &b_jet_isJetID_2018);
   fChain->SetBranchAddress("jet_isJetIDLepVeto_2018", &jet_isJetIDLepVeto_2018, &b_jet_isJetIDLepVeto_2018);
   fChain->SetBranchAddress("MET_nominal_Pt", &MET_nominal_Pt, &b_MET_nominal_Pt);
   fChain->SetBranchAddress("MET_nominal_Px", &MET_nominal_Px, &b_MET_nominal_Px);
   fChain->SetBranchAddress("MET_nominal_Py", &MET_nominal_Py, &b_MET_nominal_Py);
   fChain->SetBranchAddress("MET_nominal_phi", &MET_nominal_phi, &b_MET_nominal_phi);
   fChain->SetBranchAddress("MET_nominal_significance", &MET_nominal_significance, &b_MET_nominal_significance);
   fChain->SetBranchAddress("MET_Pt", &MET_Pt, &b_MET_Pt);
   fChain->SetBranchAddress("MET_Px", &MET_Px, &b_MET_Px);
   fChain->SetBranchAddress("MET_Py", &MET_Py, &b_MET_Py);
   fChain->SetBranchAddress("MET_phi", &MET_phi, &b_MET_phi);
   fChain->SetBranchAddress("MET_significance", &MET_significance, &b_MET_significance);
   fChain->SetBranchAddress("MET_T1_Pt", &MET_T1_Pt, &b_MET_T1_Pt);
   fChain->SetBranchAddress("MET_T1_Px", &MET_T1_Px, &b_MET_T1_Px);
   fChain->SetBranchAddress("MET_T1_Py", &MET_T1_Py, &b_MET_T1_Py);
   fChain->SetBranchAddress("MET_T1_phi", &MET_T1_phi, &b_MET_T1_phi);
   fChain->SetBranchAddress("MET_T1_significance", &MET_T1_significance, &b_MET_T1_significance);
   fChain->SetBranchAddress("MET_T1Txy_Pt", &MET_T1Txy_Pt, &b_MET_T1Txy_Pt);
   fChain->SetBranchAddress("MET_T1Txy_Px", &MET_T1Txy_Px, &b_MET_T1Txy_Px);
   fChain->SetBranchAddress("MET_T1Txy_Py", &MET_T1Txy_Py, &b_MET_T1Txy_Py);
   fChain->SetBranchAddress("MET_T1Txy_phi", &MET_T1Txy_phi, &b_MET_T1Txy_phi);
   fChain->SetBranchAddress("MET_T1Txy_significance", &MET_T1Txy_significance, &b_MET_T1Txy_significance);
   fChain->SetBranchAddress("MET_FinalCollection_Pt", &MET_FinalCollection_Pt, &b_MET_FinalCollection_Pt);
   fChain->SetBranchAddress("MET_FinalCollection_Px", &MET_FinalCollection_Px, &b_MET_FinalCollection_Px);
   fChain->SetBranchAddress("MET_FinalCollection_Py", &MET_FinalCollection_Py, &b_MET_FinalCollection_Py);
   fChain->SetBranchAddress("MET_FinalCollection_phi", &MET_FinalCollection_phi, &b_MET_FinalCollection_phi);
   fChain->SetBranchAddress("MET_FinalCollection_significance", &MET_FinalCollection_significance, &b_MET_FinalCollection_significance);
   fChain->SetBranchAddress("tau_n", &tau_n, &b_tau_n);
   fChain->SetBranchAddress("tau_px", &tau_px, &b_tau_px);
   fChain->SetBranchAddress("tau_py", &tau_py, &b_tau_py);
   fChain->SetBranchAddress("tau_pz", &tau_pz, &b_tau_pz);
   fChain->SetBranchAddress("tau_pt", &tau_pt, &b_tau_pt);
   fChain->SetBranchAddress("tau_eta", &tau_eta, &b_tau_eta);
   fChain->SetBranchAddress("tau_theta", &tau_theta, &b_tau_theta);
   fChain->SetBranchAddress("tau_phi", &tau_phi, &b_tau_phi);
   fChain->SetBranchAddress("tau_energy", &tau_energy, &b_tau_energy);
   fChain->SetBranchAddress("tau_mass", &tau_mass, &b_tau_mass);
   fChain->SetBranchAddress("tau_dxy", &tau_dxy, &b_tau_dxy);
   fChain->SetBranchAddress("tau_dxy_error", &tau_dxy_error, &b_tau_dxy_error);
   fChain->SetBranchAddress("tau_ptLeadChargedCand", &tau_ptLeadChargedCand, &b_tau_ptLeadChargedCand);
   fChain->SetBranchAddress("tau_byIsolationMVArun2v1DBoldDMwLTraw", &tau_byIsolationMVArun2v1DBoldDMwLTraw, &b_tau_byIsolationMVArun2v1DBoldDMwLTraw);
   fChain->SetBranchAddress("tau_byCombinedIsolationDeltaBetaCorrRaw3Hits", &tau_byCombinedIsolationDeltaBetaCorrRaw3Hits, &b_tau_byCombinedIsolationDeltaBetaCorrRaw3Hits);
   fChain->SetBranchAddress("tau_byIsolationMVArun2v1DBnewDMwLTraw", &tau_byIsolationMVArun2v1DBnewDMwLTraw, &b_tau_byIsolationMVArun2v1DBnewDMwLTraw);
   fChain->SetBranchAddress("tau_byIsolationMVArun2v1PWoldDMwLTraw", &tau_byIsolationMVArun2v1PWoldDMwLTraw, &b_tau_byIsolationMVArun2v1PWoldDMwLTraw);
   fChain->SetBranchAddress("tau_byIsolationMVArun2v1PWnewDMwLTraw", &tau_byIsolationMVArun2v1PWnewDMwLTraw, &b_tau_byIsolationMVArun2v1PWnewDMwLTraw);
   fChain->SetBranchAddress("tau_byIsolationMVArun2v1DBdR03oldDMwLTraw", &tau_byIsolationMVArun2v1DBdR03oldDMwLTraw, &b_tau_byIsolationMVArun2v1DBdR03oldDMwLTraw);
   fChain->SetBranchAddress("tau_byIsolationMVArun2v1PWdR03oldDMwLTraw", &tau_byIsolationMVArun2v1PWdR03oldDMwLTraw, &b_tau_byIsolationMVArun2v1PWdR03oldDMwLTraw);
   fChain->SetBranchAddress("tau_againstElectronMVA6Raw", &tau_againstElectronMVA6Raw, &b_tau_againstElectronMVA6Raw);
   fChain->SetBranchAddress("tau_againstElectronMVA6category", &tau_againstElectronMVA6category, &b_tau_againstElectronMVA6category);
   fChain->SetBranchAddress("tau_mc_bestDR", &tau_mc_bestDR, &b_tau_mc_bestDR);
   fChain->SetBranchAddress("tau_mc_ERatio", &tau_mc_ERatio, &b_tau_mc_ERatio);
   fChain->SetBranchAddress("tau_decayModeFinding", &tau_decayModeFinding, &b_tau_decayModeFinding);
   fChain->SetBranchAddress("tau_decayModeFindingNewDMs", &tau_decayModeFindingNewDMs, &b_tau_decayModeFindingNewDMs);
   fChain->SetBranchAddress("tau_againstMuonLoose3", &tau_againstMuonLoose3, &b_tau_againstMuonLoose3);
   fChain->SetBranchAddress("tau_againstMuonTight3", &tau_againstMuonTight3, &b_tau_againstMuonTight3);
   fChain->SetBranchAddress("tau_byLooseCombinedIsolationDeltaBetaCorr3Hits", &tau_byLooseCombinedIsolationDeltaBetaCorr3Hits, &b_tau_byLooseCombinedIsolationDeltaBetaCorr3Hits);
   fChain->SetBranchAddress("tau_byMediumCombinedIsolationDeltaBetaCorr3Hits", &tau_byMediumCombinedIsolationDeltaBetaCorr3Hits, &b_tau_byMediumCombinedIsolationDeltaBetaCorr3Hits);
   fChain->SetBranchAddress("tau_byTightCombinedIsolationDeltaBetaCorr3Hits", &tau_byTightCombinedIsolationDeltaBetaCorr3Hits, &b_tau_byTightCombinedIsolationDeltaBetaCorr3Hits);
   fChain->SetBranchAddress("tau_byVLooseIsolationMVArun2v1DBoldDMwLT", &tau_byVLooseIsolationMVArun2v1DBoldDMwLT, &b_tau_byVLooseIsolationMVArun2v1DBoldDMwLT);
   fChain->SetBranchAddress("tau_byLooseIsolationMVArun2v1DBoldDMwLT", &tau_byLooseIsolationMVArun2v1DBoldDMwLT, &b_tau_byLooseIsolationMVArun2v1DBoldDMwLT);
   fChain->SetBranchAddress("tau_byMediumIsolationMVArun2v1DBoldDMwLT", &tau_byMediumIsolationMVArun2v1DBoldDMwLT, &b_tau_byMediumIsolationMVArun2v1DBoldDMwLT);
   fChain->SetBranchAddress("tau_byTightIsolationMVArun2v1DBoldDMwLT", &tau_byTightIsolationMVArun2v1DBoldDMwLT, &b_tau_byTightIsolationMVArun2v1DBoldDMwLT);
   fChain->SetBranchAddress("tau_byVTightIsolationMVArun2v1DBoldDMwLT", &tau_byVTightIsolationMVArun2v1DBoldDMwLT, &b_tau_byVTightIsolationMVArun2v1DBoldDMwLT);
   fChain->SetBranchAddress("tau_byVVTightIsolationMVArun2v1DBoldDMwLT", &tau_byVVTightIsolationMVArun2v1DBoldDMwLT, &b_tau_byVVTightIsolationMVArun2v1DBoldDMwLT);
   fChain->SetBranchAddress("tau_byVLooseIsolationMVArun2v1DBnewDMwLT", &tau_byVLooseIsolationMVArun2v1DBnewDMwLT, &b_tau_byVLooseIsolationMVArun2v1DBnewDMwLT);
   fChain->SetBranchAddress("tau_byLooseIsolationMVArun2v1DBnewDMwLT", &tau_byLooseIsolationMVArun2v1DBnewDMwLT, &b_tau_byLooseIsolationMVArun2v1DBnewDMwLT);
   fChain->SetBranchAddress("tau_byMediumIsolationMVArun2v1DBnewDMwLT", &tau_byMediumIsolationMVArun2v1DBnewDMwLT, &b_tau_byMediumIsolationMVArun2v1DBnewDMwLT);
   fChain->SetBranchAddress("tau_byTightIsolationMVArun2v1DBnewDMwLT", &tau_byTightIsolationMVArun2v1DBnewDMwLT, &b_tau_byTightIsolationMVArun2v1DBnewDMwLT);
   fChain->SetBranchAddress("tau_byVTightIsolationMVArun2v1DBnewDMwLT", &tau_byVTightIsolationMVArun2v1DBnewDMwLT, &b_tau_byVTightIsolationMVArun2v1DBnewDMwLT);
   fChain->SetBranchAddress("tau_byVVTightIsolationMVArun2v1DBnewDMwLT", &tau_byVVTightIsolationMVArun2v1DBnewDMwLT, &b_tau_byVVTightIsolationMVArun2v1DBnewDMwLT);
   fChain->SetBranchAddress("tau_byVLooseIsolationMVArun2v1PWoldDMwLT", &tau_byVLooseIsolationMVArun2v1PWoldDMwLT, &b_tau_byVLooseIsolationMVArun2v1PWoldDMwLT);
   fChain->SetBranchAddress("tau_byLooseIsolationMVArun2v1PWoldDMwLT", &tau_byLooseIsolationMVArun2v1PWoldDMwLT, &b_tau_byLooseIsolationMVArun2v1PWoldDMwLT);
   fChain->SetBranchAddress("tau_byMediumIsolationMVArun2v1PWoldDMwLT", &tau_byMediumIsolationMVArun2v1PWoldDMwLT, &b_tau_byMediumIsolationMVArun2v1PWoldDMwLT);
   fChain->SetBranchAddress("tau_byTightIsolationMVArun2v1PWoldDMwLT", &tau_byTightIsolationMVArun2v1PWoldDMwLT, &b_tau_byTightIsolationMVArun2v1PWoldDMwLT);
   fChain->SetBranchAddress("tau_byVTightIsolationMVArun2v1PWoldDMwLT", &tau_byVTightIsolationMVArun2v1PWoldDMwLT, &b_tau_byVTightIsolationMVArun2v1PWoldDMwLT);
   fChain->SetBranchAddress("tau_byVVTightIsolationMVArun2v1PWoldDMwLT", &tau_byVVTightIsolationMVArun2v1PWoldDMwLT, &b_tau_byVVTightIsolationMVArun2v1PWoldDMwLT);
   fChain->SetBranchAddress("tau_byVLooseIsolationMVArun2v1PWnewDMwLT", &tau_byVLooseIsolationMVArun2v1PWnewDMwLT, &b_tau_byVLooseIsolationMVArun2v1PWnewDMwLT);
   fChain->SetBranchAddress("tau_byLooseIsolationMVArun2v1PWnewDMwLT", &tau_byLooseIsolationMVArun2v1PWnewDMwLT, &b_tau_byLooseIsolationMVArun2v1PWnewDMwLT);
   fChain->SetBranchAddress("tau_byMediumIsolationMVArun2v1PWnewDMwLT", &tau_byMediumIsolationMVArun2v1PWnewDMwLT, &b_tau_byMediumIsolationMVArun2v1PWnewDMwLT);
   fChain->SetBranchAddress("tau_byTightIsolationMVArun2v1PWnewDMwLT", &tau_byTightIsolationMVArun2v1PWnewDMwLT, &b_tau_byTightIsolationMVArun2v1PWnewDMwLT);
   fChain->SetBranchAddress("tau_byVTightIsolationMVArun2v1PWnewDMwLT", &tau_byVTightIsolationMVArun2v1PWnewDMwLT, &b_tau_byVTightIsolationMVArun2v1PWnewDMwLT);
   fChain->SetBranchAddress("tau_byVVTightIsolationMVArun2v1PWnewDMwLT", &tau_byVVTightIsolationMVArun2v1PWnewDMwLT, &b_tau_byVVTightIsolationMVArun2v1PWnewDMwLT);
   fChain->SetBranchAddress("tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT", &tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT, &b_tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT);
   fChain->SetBranchAddress("tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT", &tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT, &b_tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT);
   fChain->SetBranchAddress("tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT", &tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT, &b_tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT);
   fChain->SetBranchAddress("tau_byTightIsolationMVArun2v1DBdR03oldDMwLT", &tau_byTightIsolationMVArun2v1DBdR03oldDMwLT, &b_tau_byTightIsolationMVArun2v1DBdR03oldDMwLT);
   fChain->SetBranchAddress("tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT", &tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT, &b_tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT);
   fChain->SetBranchAddress("tau_byVVTightIsolationMVArun2v1DBdR03oldDMwLT", &tau_byVVTightIsolationMVArun2v1DBdR03oldDMwLT, &b_tau_byVVTightIsolationMVArun2v1DBdR03oldDMwLT);
   fChain->SetBranchAddress("tau_byVLooseIsolationMVArun2v1PWdR03oldDMwLT", &tau_byVLooseIsolationMVArun2v1PWdR03oldDMwLT, &b_tau_byVLooseIsolationMVArun2v1PWdR03oldDMwLT);
   fChain->SetBranchAddress("tau_byLooseIsolationMVArun2v1PWdR03oldDMwLT", &tau_byLooseIsolationMVArun2v1PWdR03oldDMwLT, &b_tau_byLooseIsolationMVArun2v1PWdR03oldDMwLT);
   fChain->SetBranchAddress("tau_byMediumIsolationMVArun2v1PWdR03oldDMwLT", &tau_byMediumIsolationMVArun2v1PWdR03oldDMwLT, &b_tau_byMediumIsolationMVArun2v1PWdR03oldDMwLT);
   fChain->SetBranchAddress("tau_byTightIsolationMVArun2v1PWdR03oldDMwLT", &tau_byTightIsolationMVArun2v1PWdR03oldDMwLT, &b_tau_byTightIsolationMVArun2v1PWdR03oldDMwLT);
   fChain->SetBranchAddress("tau_byVTightIsolationMVArun2v1PWdR03oldDMwLT", &tau_byVTightIsolationMVArun2v1PWdR03oldDMwLT, &b_tau_byVTightIsolationMVArun2v1PWdR03oldDMwLT);
   fChain->SetBranchAddress("tau_byVVTightIsolationMVArun2v1PWdR03oldDMwLT", &tau_byVVTightIsolationMVArun2v1PWdR03oldDMwLT, &b_tau_byVVTightIsolationMVArun2v1PWdR03oldDMwLT);
   fChain->SetBranchAddress("tau_againstElectronVLooseMVA6", &tau_againstElectronVLooseMVA6, &b_tau_againstElectronVLooseMVA6);
   fChain->SetBranchAddress("tau_againstElectronLooseMVA6", &tau_againstElectronLooseMVA6, &b_tau_againstElectronLooseMVA6);
   fChain->SetBranchAddress("tau_againstElectronMediumMVA6", &tau_againstElectronMediumMVA6, &b_tau_againstElectronMediumMVA6);
   fChain->SetBranchAddress("tau_againstElectronTightMVA6", &tau_againstElectronTightMVA6, &b_tau_againstElectronTightMVA6);
   fChain->SetBranchAddress("tau_againstElectronVTightMVA6", &tau_againstElectronVTightMVA6, &b_tau_againstElectronVTightMVA6);
   fChain->SetBranchAddress("tau_chargedIsoPtSum", &tau_chargedIsoPtSum, &b_tau_chargedIsoPtSum);
   fChain->SetBranchAddress("tau_neutralIsoPtSum", &tau_neutralIsoPtSum, &b_tau_neutralIsoPtSum);
   fChain->SetBranchAddress("tau_puCorrPtSum", &tau_puCorrPtSum, &b_tau_puCorrPtSum);
   fChain->SetBranchAddress("tau_footprintCorrection", &tau_footprintCorrection, &b_tau_footprintCorrection);
   fChain->SetBranchAddress("tau_neutralIsoPtSumWeight", &tau_neutralIsoPtSumWeight, &b_tau_neutralIsoPtSumWeight);
   fChain->SetBranchAddress("tau_photonPtSumOutsideSignalCone", &tau_photonPtSumOutsideSignalCone, &b_tau_photonPtSumOutsideSignalCone);
   fChain->SetBranchAddress("tau_byPhotonPtSumOutsideSignalCone", &tau_byPhotonPtSumOutsideSignalCone, &b_tau_byPhotonPtSumOutsideSignalCone);
   fChain->SetBranchAddress("tau_footprintCorrectiondR03", &tau_footprintCorrectiondR03, &b_tau_footprintCorrectiondR03);
   fChain->SetBranchAddress("tau_chargedIsoPtSumdR03", &tau_chargedIsoPtSumdR03, &b_tau_chargedIsoPtSumdR03);
   fChain->SetBranchAddress("tau_neutralIsoPtSumWeightdR03", &tau_neutralIsoPtSumWeightdR03, &b_tau_neutralIsoPtSumWeightdR03);
   fChain->SetBranchAddress("tau_neutralIsoPtSumdR03", &tau_neutralIsoPtSumdR03, &b_tau_neutralIsoPtSumdR03);
   fChain->SetBranchAddress("tau_photonPtSumOutsideSignalConedR03", &tau_photonPtSumOutsideSignalConedR03, &b_tau_photonPtSumOutsideSignalConedR03);
   fChain->SetBranchAddress("tau_PFChargedHadIso", &tau_PFChargedHadIso, &b_tau_PFChargedHadIso);
   fChain->SetBranchAddress("tau_PFNeutralHadIso", &tau_PFNeutralHadIso, &b_tau_PFNeutralHadIso);
   fChain->SetBranchAddress("tau_PFPhotonIso", &tau_PFPhotonIso, &b_tau_PFPhotonIso);
   fChain->SetBranchAddress("tau_leadChargedParticlePt", &tau_leadChargedParticlePt, &b_tau_leadChargedParticlePt);
   fChain->SetBranchAddress("tau_trackRefPt", &tau_trackRefPt, &b_tau_trackRefPt);
   fChain->SetBranchAddress("tau_lead_dxy", &tau_lead_dxy, &b_tau_lead_dxy);
   fChain->SetBranchAddress("tau_lead_dz", &tau_lead_dz, &b_tau_lead_dz);
   fChain->SetBranchAddress("tau_dxy_Sig", &tau_dxy_Sig, &b_tau_dxy_Sig);
   fChain->SetBranchAddress("tau_flightLengthSig", &tau_flightLengthSig, &b_tau_flightLengthSig);
   fChain->SetBranchAddress("tau_ip3d", &tau_ip3d, &b_tau_ip3d);
   fChain->SetBranchAddress("tau_ip3d_Sig", &tau_ip3d_Sig, &b_tau_ip3d_Sig);
   fChain->SetBranchAddress("tau_decayDistX", &tau_decayDistX, &b_tau_decayDistX);
   fChain->SetBranchAddress("tau_decayDistY", &tau_decayDistY, &b_tau_decayDistY);
   fChain->SetBranchAddress("tau_decayDistZ", &tau_decayDistZ, &b_tau_decayDistZ);
   fChain->SetBranchAddress("tau_decayDistMag", &tau_decayDistMag, &b_tau_decayDistMag);
   fChain->SetBranchAddress("tau_nPhoton", &tau_nPhoton, &b_tau_nPhoton);
   fChain->SetBranchAddress("tau_ptWeightedDetaStrip", &tau_ptWeightedDetaStrip, &b_tau_ptWeightedDetaStrip);
   fChain->SetBranchAddress("tau_ptWeightedDphiStrip", &tau_ptWeightedDphiStrip, &b_tau_ptWeightedDphiStrip);
   fChain->SetBranchAddress("tau_ptWeightedDrSignal", &tau_ptWeightedDrSignal, &b_tau_ptWeightedDrSignal);
   fChain->SetBranchAddress("tau_ptWeightedDrIsolation", &tau_ptWeightedDrIsolation, &b_tau_ptWeightedDrIsolation);
   fChain->SetBranchAddress("tau_leadingTrackChi2", &tau_leadingTrackChi2, &b_tau_leadingTrackChi2);
   fChain->SetBranchAddress("tau_eRatio", &tau_eRatio, &b_tau_eRatio);
   fChain->SetBranchAddress("tau_gjAngleDiff", &tau_gjAngleDiff, &b_tau_gjAngleDiff);
   fChain->SetBranchAddress("tau_numberOfIsolationChargedHadrCands", &tau_numberOfIsolationChargedHadrCands, &b_tau_numberOfIsolationChargedHadrCands);
   fChain->SetBranchAddress("tau_numberOfSignalChargedHadrCands", &tau_numberOfSignalChargedHadrCands, &b_tau_numberOfSignalChargedHadrCands);
   fChain->SetBranchAddress("tau_numNeutralHadronsSignalCone", &tau_numNeutralHadronsSignalCone, &b_tau_numNeutralHadronsSignalCone);
   fChain->SetBranchAddress("tau_numPhotonsSignalCone", &tau_numPhotonsSignalCone, &b_tau_numPhotonsSignalCone);
   fChain->SetBranchAddress("tau_numParticlesSignalCone", &tau_numParticlesSignalCone, &b_tau_numParticlesSignalCone);
   fChain->SetBranchAddress("tau_numChargedParticlesIsoCone", &tau_numChargedParticlesIsoCone, &b_tau_numChargedParticlesIsoCone);
   fChain->SetBranchAddress("tau_numNeutralHadronsIsoCone", &tau_numNeutralHadronsIsoCone, &b_tau_numNeutralHadronsIsoCone);
   fChain->SetBranchAddress("tau_numPhotonsIsoCone", &tau_numPhotonsIsoCone, &b_tau_numPhotonsIsoCone);
   fChain->SetBranchAddress("tau_numParticlesIsoCone", &tau_numParticlesIsoCone, &b_tau_numParticlesIsoCone);
   fChain->SetBranchAddress("tau_mc_index", &tau_mc_index, &b_tau_mc_index);
   fChain->SetBranchAddress("tau_decayMode", &tau_decayMode, &b_tau_decayMode);
   fChain->SetBranchAddress("tau_charge", &tau_charge, &b_tau_charge);
   fChain->SetBranchAddress("tau_isPFTau", &tau_isPFTau, &b_tau_isPFTau);
   fChain->SetBranchAddress("tau_hasSecondaryVertex", &tau_hasSecondaryVertex, &b_tau_hasSecondaryVertex);
   fChain->SetBranchAddress("tau_leadChargedHadrAvailable", &tau_leadChargedHadrAvailable, &b_tau_leadChargedHadrAvailable);
   fChain->SetBranchAddress("tau_byIsolationMVArun2017v1DBoldDMwLTraw2017", &tau_byIsolationMVArun2017v1DBoldDMwLTraw2017, &b_tau_byIsolationMVArun2017v1DBoldDMwLTraw2017);
   fChain->SetBranchAddress("tau_byIsolationMVArun2017v2DBnewDMwLTraw2017", &tau_byIsolationMVArun2017v2DBnewDMwLTraw2017, &b_tau_byIsolationMVArun2017v2DBnewDMwLTraw2017);
   fChain->SetBranchAddress("tau_byIsolationMVArun2017v2DBoldDMdR0p3wLTraw2017", &tau_byIsolationMVArun2017v2DBoldDMdR0p3wLTraw2017, &b_tau_byIsolationMVArun2017v2DBoldDMdR0p3wLTraw2017);
   fChain->SetBranchAddress("tau_byIsolationMVArun2017v2DBoldDMwLTraw2017", &tau_byIsolationMVArun2017v2DBoldDMwLTraw2017, &b_tau_byIsolationMVArun2017v2DBoldDMwLTraw2017);
   fChain->SetBranchAddress("tau_byIsolationMVArun2v1DBnewDMwLTraw2016", &tau_byIsolationMVArun2v1DBnewDMwLTraw2016, &b_tau_byIsolationMVArun2v1DBnewDMwLTraw2016);
   fChain->SetBranchAddress("tau_byIsolationMVArun2v1DBoldDMwLTraw2016", &tau_byIsolationMVArun2v1DBoldDMwLTraw2016, &b_tau_byIsolationMVArun2v1DBoldDMwLTraw2016);
   fChain->SetBranchAddress("tau_byVVLooseIsolationMVArun2017v1DBoldDMwLT2017", &tau_byVVLooseIsolationMVArun2017v1DBoldDMwLT2017, &b_tau_byVVLooseIsolationMVArun2017v1DBoldDMwLT2017);
   fChain->SetBranchAddress("tau_byVLooseIsolationMVArun2017v1DBoldDMwLT2017", &tau_byVLooseIsolationMVArun2017v1DBoldDMwLT2017, &b_tau_byVLooseIsolationMVArun2017v1DBoldDMwLT2017);
   fChain->SetBranchAddress("tau_byLooseIsolationMVArun2017v1DBoldDMwLT2017", &tau_byLooseIsolationMVArun2017v1DBoldDMwLT2017, &b_tau_byLooseIsolationMVArun2017v1DBoldDMwLT2017);
   fChain->SetBranchAddress("tau_byMediumIsolationMVArun2017v1DBoldDMwLT2017", &tau_byMediumIsolationMVArun2017v1DBoldDMwLT2017, &b_tau_byMediumIsolationMVArun2017v1DBoldDMwLT2017);
   fChain->SetBranchAddress("tau_byTightIsolationMVArun2017v1DBoldDMwLT2017", &tau_byTightIsolationMVArun2017v1DBoldDMwLT2017, &b_tau_byTightIsolationMVArun2017v1DBoldDMwLT2017);
   fChain->SetBranchAddress("tau_byVTightIsolationMVArun2017v1DBoldDMwLT2017", &tau_byVTightIsolationMVArun2017v1DBoldDMwLT2017, &b_tau_byVTightIsolationMVArun2017v1DBoldDMwLT2017);
   fChain->SetBranchAddress("tau_byVVTightIsolationMVArun2017v1DBoldDMwLT2017", &tau_byVVTightIsolationMVArun2017v1DBoldDMwLT2017, &b_tau_byVVTightIsolationMVArun2017v1DBoldDMwLT2017);
   fChain->SetBranchAddress("tau_byVVLooseIsolationMVArun2017v2DBnewDMwLT2017", &tau_byVVLooseIsolationMVArun2017v2DBnewDMwLT2017, &b_tau_byVVLooseIsolationMVArun2017v2DBnewDMwLT2017);
   fChain->SetBranchAddress("tau_byVLooseIsolationMVArun2017v2DBnewDMwLT2017", &tau_byVLooseIsolationMVArun2017v2DBnewDMwLT2017, &b_tau_byVLooseIsolationMVArun2017v2DBnewDMwLT2017);
   fChain->SetBranchAddress("tau_byLooseIsolationMVArun2017v2DBnewDMwLT2017", &tau_byLooseIsolationMVArun2017v2DBnewDMwLT2017, &b_tau_byLooseIsolationMVArun2017v2DBnewDMwLT2017);
   fChain->SetBranchAddress("tau_byMediumIsolationMVArun2017v2DBnewDMwLT2017", &tau_byMediumIsolationMVArun2017v2DBnewDMwLT2017, &b_tau_byMediumIsolationMVArun2017v2DBnewDMwLT2017);
   fChain->SetBranchAddress("tau_byTightIsolationMVArun2017v2DBnewDMwLT2017", &tau_byTightIsolationMVArun2017v2DBnewDMwLT2017, &b_tau_byTightIsolationMVArun2017v2DBnewDMwLT2017);
   fChain->SetBranchAddress("tau_byVTightIsolationMVArun2017v2DBnewDMwLT2017", &tau_byVTightIsolationMVArun2017v2DBnewDMwLT2017, &b_tau_byVTightIsolationMVArun2017v2DBnewDMwLT2017);
   fChain->SetBranchAddress("tau_byVVTightIsolationMVArun2017v2DBnewDMwLT2017", &tau_byVVTightIsolationMVArun2017v2DBnewDMwLT2017, &b_tau_byVVTightIsolationMVArun2017v2DBnewDMwLT2017);
   fChain->SetBranchAddress("tau_byVVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017", &tau_byVVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017, &b_tau_byVVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017);
   fChain->SetBranchAddress("tau_byVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017", &tau_byVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017, &b_tau_byVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017);
   fChain->SetBranchAddress("tau_byLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017", &tau_byLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017, &b_tau_byLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017);
   fChain->SetBranchAddress("tau_byMediumIsolationMVArun2017v2DBoldDMdR0p3wLT2017", &tau_byMediumIsolationMVArun2017v2DBoldDMdR0p3wLT2017, &b_tau_byMediumIsolationMVArun2017v2DBoldDMdR0p3wLT2017);
   fChain->SetBranchAddress("tau_byTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017", &tau_byTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017, &b_tau_byTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017);
   fChain->SetBranchAddress("tau_byVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017", &tau_byVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017, &b_tau_byVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017);
   fChain->SetBranchAddress("tau_byVVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017", &tau_byVVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017, &b_tau_byVVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017);
   fChain->SetBranchAddress("tau_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017", &tau_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017, &b_tau_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017);
   fChain->SetBranchAddress("tau_byVLooseIsolationMVArun2017v2DBoldDMwLT2017", &tau_byVLooseIsolationMVArun2017v2DBoldDMwLT2017, &b_tau_byVLooseIsolationMVArun2017v2DBoldDMwLT2017);
   fChain->SetBranchAddress("tau_byLooseIsolationMVArun2017v2DBoldDMwLT2017", &tau_byLooseIsolationMVArun2017v2DBoldDMwLT2017, &b_tau_byLooseIsolationMVArun2017v2DBoldDMwLT2017);
   fChain->SetBranchAddress("tau_byMediumIsolationMVArun2017v2DBoldDMwLT2017", &tau_byMediumIsolationMVArun2017v2DBoldDMwLT2017, &b_tau_byMediumIsolationMVArun2017v2DBoldDMwLT2017);
   fChain->SetBranchAddress("tau_byTightIsolationMVArun2017v2DBoldDMwLT2017", &tau_byTightIsolationMVArun2017v2DBoldDMwLT2017, &b_tau_byTightIsolationMVArun2017v2DBoldDMwLT2017);
   fChain->SetBranchAddress("tau_byVTightIsolationMVArun2017v2DBoldDMwLT2017", &tau_byVTightIsolationMVArun2017v2DBoldDMwLT2017, &b_tau_byVTightIsolationMVArun2017v2DBoldDMwLT2017);
   fChain->SetBranchAddress("tau_byVVTightIsolationMVArun2017v2DBoldDMwLT2017", &tau_byVVTightIsolationMVArun2017v2DBoldDMwLT2017, &b_tau_byVVTightIsolationMVArun2017v2DBoldDMwLT2017);
   fChain->SetBranchAddress("tau_byVLooseIsolationMVArun2v1DBnewDMwLT2016", &tau_byVLooseIsolationMVArun2v1DBnewDMwLT2016, &b_tau_byVLooseIsolationMVArun2v1DBnewDMwLT2016);
   fChain->SetBranchAddress("tau_byLooseIsolationMVArun2v1DBnewDMwLT2016", &tau_byLooseIsolationMVArun2v1DBnewDMwLT2016, &b_tau_byLooseIsolationMVArun2v1DBnewDMwLT2016);
   fChain->SetBranchAddress("tau_byMediumIsolationMVArun2v1DBnewDMwLT2016", &tau_byMediumIsolationMVArun2v1DBnewDMwLT2016, &b_tau_byMediumIsolationMVArun2v1DBnewDMwLT2016);
   fChain->SetBranchAddress("tau_byTightIsolationMVArun2v1DBnewDMwLT2016", &tau_byTightIsolationMVArun2v1DBnewDMwLT2016, &b_tau_byTightIsolationMVArun2v1DBnewDMwLT2016);
   fChain->SetBranchAddress("tau_byVTightIsolationMVArun2v1DBnewDMwLT2016", &tau_byVTightIsolationMVArun2v1DBnewDMwLT2016, &b_tau_byVTightIsolationMVArun2v1DBnewDMwLT2016);
   fChain->SetBranchAddress("tau_byVVTightIsolationMVArun2v1DBnewDMwLT2016", &tau_byVVTightIsolationMVArun2v1DBnewDMwLT2016, &b_tau_byVVTightIsolationMVArun2v1DBnewDMwLT2016);
   fChain->SetBranchAddress("tau_byVLooseIsolationMVArun2v1DBoldDMwLT2016", &tau_byVLooseIsolationMVArun2v1DBoldDMwLT2016, &b_tau_byVLooseIsolationMVArun2v1DBoldDMwLT2016);
   fChain->SetBranchAddress("tau_byLooseIsolationMVArun2v1DBoldDMwLT2016", &tau_byLooseIsolationMVArun2v1DBoldDMwLT2016, &b_tau_byLooseIsolationMVArun2v1DBoldDMwLT2016);
   fChain->SetBranchAddress("tau_byMediumIsolationMVArun2v1DBoldDMwLT2016", &tau_byMediumIsolationMVArun2v1DBoldDMwLT2016, &b_tau_byMediumIsolationMVArun2v1DBoldDMwLT2016);
   fChain->SetBranchAddress("tau_byTightIsolationMVArun2v1DBoldDMwLT2016", &tau_byTightIsolationMVArun2v1DBoldDMwLT2016, &b_tau_byTightIsolationMVArun2v1DBoldDMwLT2016);
   fChain->SetBranchAddress("tau_byVTightIsolationMVArun2v1DBoldDMwLT2016", &tau_byVTightIsolationMVArun2v1DBoldDMwLT2016, &b_tau_byVTightIsolationMVArun2v1DBoldDMwLT2016);
   fChain->SetBranchAddress("tau_byVVTightIsolationMVArun2v1DBoldDMwLT2016", &tau_byVVTightIsolationMVArun2v1DBoldDMwLT2016, &b_tau_byVVTightIsolationMVArun2v1DBoldDMwLT2016);
   fChain->SetBranchAddress("trig_Flag_ecalBadCalibReduced", &trig_Flag_ecalBadCalibReduced, &b_trig_Flag_ecalBadCalibReduced);
   fChain->SetBranchAddress("trig_raw2digi_step_accept", &trig_raw2digi_step_accept, &b_trig_raw2digi_step_accept);
   fChain->SetBranchAddress("trig_L1Reco_step_accept", &trig_L1Reco_step_accept, &b_trig_L1Reco_step_accept);
   fChain->SetBranchAddress("trig_reconstruction_step_accept", &trig_reconstruction_step_accept, &b_trig_reconstruction_step_accept);
   fChain->SetBranchAddress("trig_EXOMONOPOLEPath_accept", &trig_EXOMONOPOLEPath_accept, &b_trig_EXOMONOPOLEPath_accept);
   fChain->SetBranchAddress("trig_ZElectronPathGsf_accept", &trig_ZElectronPathGsf_accept, &b_trig_ZElectronPathGsf_accept);
   fChain->SetBranchAddress("trig_pathALCARECOEcalESAlign_accept", &trig_pathALCARECOEcalESAlign_accept, &b_trig_pathALCARECOEcalESAlign_accept);
   fChain->SetBranchAddress("trig_pathALCARECOEcalUncalWElectron_accept", &trig_pathALCARECOEcalUncalWElectron_accept, &b_trig_pathALCARECOEcalUncalWElectron_accept);
   fChain->SetBranchAddress("trig_pathALCARECOHcalCalIsoTrkFilter_accept", &trig_pathALCARECOHcalCalIsoTrkFilter_accept, &b_trig_pathALCARECOHcalCalIsoTrkFilter_accept);
   fChain->SetBranchAddress("trig_pathALCARECOHcalCalIterativePhiSym_accept", &trig_pathALCARECOHcalCalIterativePhiSym_accept, &b_trig_pathALCARECOHcalCalIterativePhiSym_accept);
   fChain->SetBranchAddress("trig_pathALCARECOEcalUncalZElectron_accept", &trig_pathALCARECOEcalUncalZElectron_accept, &b_trig_pathALCARECOEcalUncalZElectron_accept);
   fChain->SetBranchAddress("trig_pathALCARECOEcalUncalZSCElectron_accept", &trig_pathALCARECOEcalUncalZSCElectron_accept, &b_trig_pathALCARECOEcalUncalZSCElectron_accept);
   fChain->SetBranchAddress("trig_eventinterpretaion_step_accept", &trig_eventinterpretaion_step_accept, &b_trig_eventinterpretaion_step_accept);
   fChain->SetBranchAddress("trig_Flag_HBHENoiseFilter_accept", &trig_Flag_HBHENoiseFilter_accept, &b_trig_Flag_HBHENoiseFilter_accept);
   fChain->SetBranchAddress("trig_Flag_HBHENoiseIsoFilter_accept", &trig_Flag_HBHENoiseIsoFilter_accept, &b_trig_Flag_HBHENoiseIsoFilter_accept);
   fChain->SetBranchAddress("trig_Flag_CSCTightHaloFilter_accept", &trig_Flag_CSCTightHaloFilter_accept, &b_trig_Flag_CSCTightHaloFilter_accept);
   fChain->SetBranchAddress("trig_Flag_CSCTightHaloTrkMuUnvetoFilter_accept", &trig_Flag_CSCTightHaloTrkMuUnvetoFilter_accept, &b_trig_Flag_CSCTightHaloTrkMuUnvetoFilter_accept);
   fChain->SetBranchAddress("trig_Flag_CSCTightHalo2015Filter_accept", &trig_Flag_CSCTightHalo2015Filter_accept, &b_trig_Flag_CSCTightHalo2015Filter_accept);
   fChain->SetBranchAddress("trig_Flag_globalTightHalo2016Filter_accept", &trig_Flag_globalTightHalo2016Filter_accept, &b_trig_Flag_globalTightHalo2016Filter_accept);
   fChain->SetBranchAddress("trig_Flag_globalSuperTightHalo2016Filter_accept", &trig_Flag_globalSuperTightHalo2016Filter_accept, &b_trig_Flag_globalSuperTightHalo2016Filter_accept);
   fChain->SetBranchAddress("trig_Flag_HcalStripHaloFilter_accept", &trig_Flag_HcalStripHaloFilter_accept, &b_trig_Flag_HcalStripHaloFilter_accept);
   fChain->SetBranchAddress("trig_Flag_hcalLaserEventFilter_accept", &trig_Flag_hcalLaserEventFilter_accept, &b_trig_Flag_hcalLaserEventFilter_accept);
   fChain->SetBranchAddress("trig_Flag_EcalDeadCellTriggerPrimitiveFilter_accept", &trig_Flag_EcalDeadCellTriggerPrimitiveFilter_accept, &b_trig_Flag_EcalDeadCellTriggerPrimitiveFilter_accept);
   fChain->SetBranchAddress("trig_Flag_EcalDeadCellBoundaryEnergyFilter_accept", &trig_Flag_EcalDeadCellBoundaryEnergyFilter_accept, &b_trig_Flag_EcalDeadCellBoundaryEnergyFilter_accept);
   fChain->SetBranchAddress("trig_Flag_ecalBadCalibFilter_accept", &trig_Flag_ecalBadCalibFilter_accept, &b_trig_Flag_ecalBadCalibFilter_accept);
   fChain->SetBranchAddress("trig_Flag_goodVertices_accept", &trig_Flag_goodVertices_accept, &b_trig_Flag_goodVertices_accept);
   fChain->SetBranchAddress("trig_Flag_eeBadScFilter_accept", &trig_Flag_eeBadScFilter_accept, &b_trig_Flag_eeBadScFilter_accept);
   fChain->SetBranchAddress("trig_Flag_ecalLaserCorrFilter_accept", &trig_Flag_ecalLaserCorrFilter_accept, &b_trig_Flag_ecalLaserCorrFilter_accept);
   fChain->SetBranchAddress("trig_Flag_trkPOGFilters_accept", &trig_Flag_trkPOGFilters_accept, &b_trig_Flag_trkPOGFilters_accept);
   fChain->SetBranchAddress("trig_Flag_chargedHadronTrackResolutionFilter_accept", &trig_Flag_chargedHadronTrackResolutionFilter_accept, &b_trig_Flag_chargedHadronTrackResolutionFilter_accept);
   fChain->SetBranchAddress("trig_Flag_muonBadTrackFilter_accept", &trig_Flag_muonBadTrackFilter_accept, &b_trig_Flag_muonBadTrackFilter_accept);
   fChain->SetBranchAddress("trig_Flag_BadChargedCandidateFilter_accept", &trig_Flag_BadChargedCandidateFilter_accept, &b_trig_Flag_BadChargedCandidateFilter_accept);
   fChain->SetBranchAddress("trig_Flag_BadPFMuonFilter_accept", &trig_Flag_BadPFMuonFilter_accept, &b_trig_Flag_BadPFMuonFilter_accept);
   fChain->SetBranchAddress("trig_Flag_BadChargedCandidateSummer16Filter_accept", &trig_Flag_BadChargedCandidateSummer16Filter_accept, &b_trig_Flag_BadChargedCandidateSummer16Filter_accept);
   fChain->SetBranchAddress("trig_Flag_BadPFMuonSummer16Filter_accept", &trig_Flag_BadPFMuonSummer16Filter_accept, &b_trig_Flag_BadPFMuonSummer16Filter_accept);
   fChain->SetBranchAddress("trig_Flag_trkPOG_manystripclus53X_accept", &trig_Flag_trkPOG_manystripclus53X_accept, &b_trig_Flag_trkPOG_manystripclus53X_accept);
   fChain->SetBranchAddress("trig_Flag_trkPOG_toomanystripclus53X_accept", &trig_Flag_trkPOG_toomanystripclus53X_accept, &b_trig_Flag_trkPOG_toomanystripclus53X_accept);
   fChain->SetBranchAddress("trig_Flag_trkPOG_logErrorTooManyClusters_accept", &trig_Flag_trkPOG_logErrorTooManyClusters_accept, &b_trig_Flag_trkPOG_logErrorTooManyClusters_accept);
   fChain->SetBranchAddress("trig_Flag_METFilters_accept", &trig_Flag_METFilters_accept, &b_trig_Flag_METFilters_accept);
   fChain->SetBranchAddress("trig_HLT_DoubleEle25_CaloIdL_MW_accept", &trig_HLT_DoubleEle25_CaloIdL_MW_accept, &b_trig_HLT_DoubleEle25_CaloIdL_MW_accept);
   fChain->SetBranchAddress("trig_HLT_DoubleEle25_CaloIdL_MW_hltL1sSingleAndDoubleEGNonIsoOrWithEG26WithJetAndTau_eta", &trig_HLT_DoubleEle25_CaloIdL_MW_hltL1sSingleAndDoubleEGNonIsoOrWithEG26WithJetAndTau_eta, &b_trig_HLT_DoubleEle25_CaloIdL_MW_hltL1sSingleAndDoubleEGNonIsoOrWithEG26WithJetAndTau_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleEle25_CaloIdL_MW_hltL1sSingleAndDoubleEGNonIsoOrWithEG26WithJetAndTau_phi", &trig_HLT_DoubleEle25_CaloIdL_MW_hltL1sSingleAndDoubleEGNonIsoOrWithEG26WithJetAndTau_phi, &b_trig_HLT_DoubleEle25_CaloIdL_MW_hltL1sSingleAndDoubleEGNonIsoOrWithEG26WithJetAndTau_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleEle25_CaloIdL_MW_hltL1sSingleAndDoubleEGNonIsoOrWithEG26WithJetAndTau_et", &trig_HLT_DoubleEle25_CaloIdL_MW_hltL1sSingleAndDoubleEGNonIsoOrWithEG26WithJetAndTau_et, &b_trig_HLT_DoubleEle25_CaloIdL_MW_hltL1sSingleAndDoubleEGNonIsoOrWithEG26WithJetAndTau_et);
   fChain->SetBranchAddress("trig_HLT_DoubleEle25_CaloIdL_MW_hltEGL1SingleAndDoubleEGNonIsoOrWithEG26WithJetAndTauFilter_eta", &trig_HLT_DoubleEle25_CaloIdL_MW_hltEGL1SingleAndDoubleEGNonIsoOrWithEG26WithJetAndTauFilter_eta, &b_trig_HLT_DoubleEle25_CaloIdL_MW_hltEGL1SingleAndDoubleEGNonIsoOrWithEG26WithJetAndTauFilter_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleEle25_CaloIdL_MW_hltEGL1SingleAndDoubleEGNonIsoOrWithEG26WithJetAndTauFilter_phi", &trig_HLT_DoubleEle25_CaloIdL_MW_hltEGL1SingleAndDoubleEGNonIsoOrWithEG26WithJetAndTauFilter_phi, &b_trig_HLT_DoubleEle25_CaloIdL_MW_hltEGL1SingleAndDoubleEGNonIsoOrWithEG26WithJetAndTauFilter_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleEle25_CaloIdL_MW_hltEGL1SingleAndDoubleEGNonIsoOrWithEG26WithJetAndTauFilter_et", &trig_HLT_DoubleEle25_CaloIdL_MW_hltEGL1SingleAndDoubleEGNonIsoOrWithEG26WithJetAndTauFilter_et, &b_trig_HLT_DoubleEle25_CaloIdL_MW_hltEGL1SingleAndDoubleEGNonIsoOrWithEG26WithJetAndTauFilter_et);
   fChain->SetBranchAddress("trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25EtFilter_eta", &trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25EtFilter_eta, &b_trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25EtFilter_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25EtFilter_phi", &trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25EtFilter_phi, &b_trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25EtFilter_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25EtFilter_et", &trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25EtFilter_et, &b_trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25EtFilter_et);
   fChain->SetBranchAddress("trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25HEFilter_eta", &trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25HEFilter_eta, &b_trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25HEFilter_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25HEFilter_phi", &trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25HEFilter_phi, &b_trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25HEFilter_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25HEFilter_et", &trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25HEFilter_et, &b_trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25HEFilter_et);
   fChain->SetBranchAddress("trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25CaloIdLClusterShapeFilter_eta", &trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25CaloIdLClusterShapeFilter_eta, &b_trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25CaloIdLClusterShapeFilter_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25CaloIdLClusterShapeFilter_phi", &trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25CaloIdLClusterShapeFilter_phi, &b_trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25CaloIdLClusterShapeFilter_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25CaloIdLClusterShapeFilter_et", &trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25CaloIdLClusterShapeFilter_et, &b_trig_HLT_DoubleEle25_CaloIdL_MW_hltEG25CaloIdLClusterShapeFilter_et);
   fChain->SetBranchAddress("trig_HLT_DoubleEle25_CaloIdL_MW_hltEle25CaloIdLPixelMatchFilter_eta", &trig_HLT_DoubleEle25_CaloIdL_MW_hltEle25CaloIdLPixelMatchFilter_eta, &b_trig_HLT_DoubleEle25_CaloIdL_MW_hltEle25CaloIdLPixelMatchFilter_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleEle25_CaloIdL_MW_hltEle25CaloIdLPixelMatchFilter_phi", &trig_HLT_DoubleEle25_CaloIdL_MW_hltEle25CaloIdLPixelMatchFilter_phi, &b_trig_HLT_DoubleEle25_CaloIdL_MW_hltEle25CaloIdLPixelMatchFilter_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleEle25_CaloIdL_MW_hltEle25CaloIdLPixelMatchFilter_et", &trig_HLT_DoubleEle25_CaloIdL_MW_hltEle25CaloIdLPixelMatchFilter_et, &b_trig_HLT_DoubleEle25_CaloIdL_MW_hltEle25CaloIdLPixelMatchFilter_et);
   fChain->SetBranchAddress("trig_HLT_DoubleEle25_CaloIdL_MW_hltEle25CaloIdLMWPMS2Filter_eta", &trig_HLT_DoubleEle25_CaloIdL_MW_hltEle25CaloIdLMWPMS2Filter_eta, &b_trig_HLT_DoubleEle25_CaloIdL_MW_hltEle25CaloIdLMWPMS2Filter_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleEle25_CaloIdL_MW_hltEle25CaloIdLMWPMS2Filter_phi", &trig_HLT_DoubleEle25_CaloIdL_MW_hltEle25CaloIdLMWPMS2Filter_phi, &b_trig_HLT_DoubleEle25_CaloIdL_MW_hltEle25CaloIdLMWPMS2Filter_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleEle25_CaloIdL_MW_hltEle25CaloIdLMWPMS2Filter_et", &trig_HLT_DoubleEle25_CaloIdL_MW_hltEle25CaloIdLMWPMS2Filter_et, &b_trig_HLT_DoubleEle25_CaloIdL_MW_hltEle25CaloIdLMWPMS2Filter_et);
   fChain->SetBranchAddress("trig_HLT_DoubleEle25_CaloIdL_MW_hltEgammaCandidatesWrapperUnseeded_eta", &trig_HLT_DoubleEle25_CaloIdL_MW_hltEgammaCandidatesWrapperUnseeded_eta, &b_trig_HLT_DoubleEle25_CaloIdL_MW_hltEgammaCandidatesWrapperUnseeded_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleEle25_CaloIdL_MW_hltEgammaCandidatesWrapperUnseeded_phi", &trig_HLT_DoubleEle25_CaloIdL_MW_hltEgammaCandidatesWrapperUnseeded_phi, &b_trig_HLT_DoubleEle25_CaloIdL_MW_hltEgammaCandidatesWrapperUnseeded_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleEle25_CaloIdL_MW_hltEgammaCandidatesWrapperUnseeded_et", &trig_HLT_DoubleEle25_CaloIdL_MW_hltEgammaCandidatesWrapperUnseeded_et, &b_trig_HLT_DoubleEle25_CaloIdL_MW_hltEgammaCandidatesWrapperUnseeded_et);
   fChain->SetBranchAddress("trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25EtUnseededFilter_eta", &trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25EtUnseededFilter_eta, &b_trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25EtUnseededFilter_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25EtUnseededFilter_phi", &trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25EtUnseededFilter_phi, &b_trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25EtUnseededFilter_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25EtUnseededFilter_et", &trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25EtUnseededFilter_et, &b_trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25EtUnseededFilter_et);
   fChain->SetBranchAddress("trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25HEUnseededFilter_eta", &trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25HEUnseededFilter_eta, &b_trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25HEUnseededFilter_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25HEUnseededFilter_phi", &trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25HEUnseededFilter_phi, &b_trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25HEUnseededFilter_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25HEUnseededFilter_et", &trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25HEUnseededFilter_et, &b_trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25HEUnseededFilter_et);
   fChain->SetBranchAddress("trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25CaloIdLClusterShapeUnseededFilter_eta", &trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25CaloIdLClusterShapeUnseededFilter_eta, &b_trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25CaloIdLClusterShapeUnseededFilter_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25CaloIdLClusterShapeUnseededFilter_phi", &trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25CaloIdLClusterShapeUnseededFilter_phi, &b_trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25CaloIdLClusterShapeUnseededFilter_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25CaloIdLClusterShapeUnseededFilter_et", &trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25CaloIdLClusterShapeUnseededFilter_et, &b_trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEG25CaloIdLClusterShapeUnseededFilter_et);
   fChain->SetBranchAddress("trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEle25CaloIdLPixelMatchUnseededFilter_eta", &trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEle25CaloIdLPixelMatchUnseededFilter_eta, &b_trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEle25CaloIdLPixelMatchUnseededFilter_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEle25CaloIdLPixelMatchUnseededFilter_phi", &trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEle25CaloIdLPixelMatchUnseededFilter_phi, &b_trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEle25CaloIdLPixelMatchUnseededFilter_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEle25CaloIdLPixelMatchUnseededFilter_et", &trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEle25CaloIdLPixelMatchUnseededFilter_et, &b_trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEle25CaloIdLPixelMatchUnseededFilter_et);
   fChain->SetBranchAddress("trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEle25CaloIdLMWPMS2UnseededFilter_eta", &trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEle25CaloIdLMWPMS2UnseededFilter_eta, &b_trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEle25CaloIdLMWPMS2UnseededFilter_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEle25CaloIdLMWPMS2UnseededFilter_phi", &trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEle25CaloIdLMWPMS2UnseededFilter_phi, &b_trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEle25CaloIdLMWPMS2UnseededFilter_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEle25CaloIdLMWPMS2UnseededFilter_et", &trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEle25CaloIdLMWPMS2UnseededFilter_et, &b_trig_HLT_DoubleEle25_CaloIdL_MW_hltDiEle25CaloIdLMWPMS2UnseededFilter_et);
   fChain->SetBranchAddress("trig_HLT_DoubleEle27_CaloIdL_MW_accept", &trig_HLT_DoubleEle27_CaloIdL_MW_accept, &b_trig_HLT_DoubleEle27_CaloIdL_MW_accept);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_MW_accept", &trig_HLT_DoubleEle33_CaloIdL_MW_accept, &b_trig_HLT_DoubleEle33_CaloIdL_MW_accept);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_MW_hltL1sSingleAndDoubleEGNonIsoOrWithEG26WithJetAndTau_eta", &trig_HLT_DoubleEle33_CaloIdL_MW_hltL1sSingleAndDoubleEGNonIsoOrWithEG26WithJetAndTau_eta, &b_trig_HLT_DoubleEle33_CaloIdL_MW_hltL1sSingleAndDoubleEGNonIsoOrWithEG26WithJetAndTau_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_MW_hltL1sSingleAndDoubleEGNonIsoOrWithEG26WithJetAndTau_phi", &trig_HLT_DoubleEle33_CaloIdL_MW_hltL1sSingleAndDoubleEGNonIsoOrWithEG26WithJetAndTau_phi, &b_trig_HLT_DoubleEle33_CaloIdL_MW_hltL1sSingleAndDoubleEGNonIsoOrWithEG26WithJetAndTau_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_MW_hltL1sSingleAndDoubleEGNonIsoOrWithEG26WithJetAndTau_et", &trig_HLT_DoubleEle33_CaloIdL_MW_hltL1sSingleAndDoubleEGNonIsoOrWithEG26WithJetAndTau_et, &b_trig_HLT_DoubleEle33_CaloIdL_MW_hltL1sSingleAndDoubleEGNonIsoOrWithEG26WithJetAndTau_et);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_MW_hltEGL1SingleAndDoubleEGNonIsoOrWithEG26WithJetAndTauFilter_eta", &trig_HLT_DoubleEle33_CaloIdL_MW_hltEGL1SingleAndDoubleEGNonIsoOrWithEG26WithJetAndTauFilter_eta, &b_trig_HLT_DoubleEle33_CaloIdL_MW_hltEGL1SingleAndDoubleEGNonIsoOrWithEG26WithJetAndTauFilter_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_MW_hltEGL1SingleAndDoubleEGNonIsoOrWithEG26WithJetAndTauFilter_phi", &trig_HLT_DoubleEle33_CaloIdL_MW_hltEGL1SingleAndDoubleEGNonIsoOrWithEG26WithJetAndTauFilter_phi, &b_trig_HLT_DoubleEle33_CaloIdL_MW_hltEGL1SingleAndDoubleEGNonIsoOrWithEG26WithJetAndTauFilter_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_MW_hltEGL1SingleAndDoubleEGNonIsoOrWithEG26WithJetAndTauFilter_et", &trig_HLT_DoubleEle33_CaloIdL_MW_hltEGL1SingleAndDoubleEGNonIsoOrWithEG26WithJetAndTauFilter_et, &b_trig_HLT_DoubleEle33_CaloIdL_MW_hltEGL1SingleAndDoubleEGNonIsoOrWithEG26WithJetAndTauFilter_et);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33EtFilter_eta", &trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33EtFilter_eta, &b_trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33EtFilter_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33EtFilter_phi", &trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33EtFilter_phi, &b_trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33EtFilter_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33EtFilter_et", &trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33EtFilter_et, &b_trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33EtFilter_et);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33HEFilter_eta", &trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33HEFilter_eta, &b_trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33HEFilter_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33HEFilter_phi", &trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33HEFilter_phi, &b_trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33HEFilter_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33HEFilter_et", &trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33HEFilter_et, &b_trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33HEFilter_et);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33CaloIdLClusterShapeFilter_eta", &trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33CaloIdLClusterShapeFilter_eta, &b_trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33CaloIdLClusterShapeFilter_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33CaloIdLClusterShapeFilter_phi", &trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33CaloIdLClusterShapeFilter_phi, &b_trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33CaloIdLClusterShapeFilter_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33CaloIdLClusterShapeFilter_et", &trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33CaloIdLClusterShapeFilter_et, &b_trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33CaloIdLClusterShapeFilter_et);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_MW_hltEle33CaloIdLPixelMatchFilter_eta", &trig_HLT_DoubleEle33_CaloIdL_MW_hltEle33CaloIdLPixelMatchFilter_eta, &b_trig_HLT_DoubleEle33_CaloIdL_MW_hltEle33CaloIdLPixelMatchFilter_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_MW_hltEle33CaloIdLPixelMatchFilter_phi", &trig_HLT_DoubleEle33_CaloIdL_MW_hltEle33CaloIdLPixelMatchFilter_phi, &b_trig_HLT_DoubleEle33_CaloIdL_MW_hltEle33CaloIdLPixelMatchFilter_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_MW_hltEle33CaloIdLPixelMatchFilter_et", &trig_HLT_DoubleEle33_CaloIdL_MW_hltEle33CaloIdLPixelMatchFilter_et, &b_trig_HLT_DoubleEle33_CaloIdL_MW_hltEle33CaloIdLPixelMatchFilter_et);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_MW_hltEle33CaloIdLMWPMS2Filter_eta", &trig_HLT_DoubleEle33_CaloIdL_MW_hltEle33CaloIdLMWPMS2Filter_eta, &b_trig_HLT_DoubleEle33_CaloIdL_MW_hltEle33CaloIdLMWPMS2Filter_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_MW_hltEle33CaloIdLMWPMS2Filter_phi", &trig_HLT_DoubleEle33_CaloIdL_MW_hltEle33CaloIdLMWPMS2Filter_phi, &b_trig_HLT_DoubleEle33_CaloIdL_MW_hltEle33CaloIdLMWPMS2Filter_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_MW_hltEle33CaloIdLMWPMS2Filter_et", &trig_HLT_DoubleEle33_CaloIdL_MW_hltEle33CaloIdLMWPMS2Filter_et, &b_trig_HLT_DoubleEle33_CaloIdL_MW_hltEle33CaloIdLMWPMS2Filter_et);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_MW_hltEgammaCandidatesWrapperUnseeded_eta", &trig_HLT_DoubleEle33_CaloIdL_MW_hltEgammaCandidatesWrapperUnseeded_eta, &b_trig_HLT_DoubleEle33_CaloIdL_MW_hltEgammaCandidatesWrapperUnseeded_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_MW_hltEgammaCandidatesWrapperUnseeded_phi", &trig_HLT_DoubleEle33_CaloIdL_MW_hltEgammaCandidatesWrapperUnseeded_phi, &b_trig_HLT_DoubleEle33_CaloIdL_MW_hltEgammaCandidatesWrapperUnseeded_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_MW_hltEgammaCandidatesWrapperUnseeded_et", &trig_HLT_DoubleEle33_CaloIdL_MW_hltEgammaCandidatesWrapperUnseeded_et, &b_trig_HLT_DoubleEle33_CaloIdL_MW_hltEgammaCandidatesWrapperUnseeded_et);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33EtUnseededFilter_eta", &trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33EtUnseededFilter_eta, &b_trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33EtUnseededFilter_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33EtUnseededFilter_phi", &trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33EtUnseededFilter_phi, &b_trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33EtUnseededFilter_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33EtUnseededFilter_et", &trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33EtUnseededFilter_et, &b_trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33EtUnseededFilter_et);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33HEUnseededFilter_eta", &trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33HEUnseededFilter_eta, &b_trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33HEUnseededFilter_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33HEUnseededFilter_phi", &trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33HEUnseededFilter_phi, &b_trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33HEUnseededFilter_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33HEUnseededFilter_et", &trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33HEUnseededFilter_et, &b_trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33HEUnseededFilter_et);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33CaloIdLClusterShapeUnseededFilter_eta", &trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33CaloIdLClusterShapeUnseededFilter_eta, &b_trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33CaloIdLClusterShapeUnseededFilter_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33CaloIdLClusterShapeUnseededFilter_phi", &trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33CaloIdLClusterShapeUnseededFilter_phi, &b_trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33CaloIdLClusterShapeUnseededFilter_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33CaloIdLClusterShapeUnseededFilter_et", &trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33CaloIdLClusterShapeUnseededFilter_et, &b_trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEG33CaloIdLClusterShapeUnseededFilter_et);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEle33CaloIdLPixelMatchUnseededFilter_eta", &trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEle33CaloIdLPixelMatchUnseededFilter_eta, &b_trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEle33CaloIdLPixelMatchUnseededFilter_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEle33CaloIdLPixelMatchUnseededFilter_phi", &trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEle33CaloIdLPixelMatchUnseededFilter_phi, &b_trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEle33CaloIdLPixelMatchUnseededFilter_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEle33CaloIdLPixelMatchUnseededFilter_et", &trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEle33CaloIdLPixelMatchUnseededFilter_et, &b_trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEle33CaloIdLPixelMatchUnseededFilter_et);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEle33CaloIdLMWPMS2UnseededFilter_eta", &trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEle33CaloIdLMWPMS2UnseededFilter_eta, &b_trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEle33CaloIdLMWPMS2UnseededFilter_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEle33CaloIdLMWPMS2UnseededFilter_phi", &trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEle33CaloIdLMWPMS2UnseededFilter_phi, &b_trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEle33CaloIdLMWPMS2UnseededFilter_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEle33CaloIdLMWPMS2UnseededFilter_et", &trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEle33CaloIdLMWPMS2UnseededFilter_et, &b_trig_HLT_DoubleEle33_CaloIdL_MW_hltDiEle33CaloIdLMWPMS2UnseededFilter_et);
   fChain->SetBranchAddress("trig_HLT_DoubleEle24_eta2p1_WPTight_Gsf_accept", &trig_HLT_DoubleEle24_eta2p1_WPTight_Gsf_accept, &b_trig_HLT_DoubleEle24_eta2p1_WPTight_Gsf_accept);
   fChain->SetBranchAddress("trig_HLT_Ele27_Ele37_CaloIdL_MW_accept", &trig_HLT_Ele27_Ele37_CaloIdL_MW_accept, &b_trig_HLT_Ele27_Ele37_CaloIdL_MW_accept);
   fChain->SetBranchAddress("trig_HLT_Mu27_Ele37_CaloIdL_MW_accept", &trig_HLT_Mu27_Ele37_CaloIdL_MW_accept, &b_trig_HLT_Mu27_Ele37_CaloIdL_MW_accept);
   fChain->SetBranchAddress("trig_HLT_Mu37_Ele27_CaloIdL_MW_accept", &trig_HLT_Mu37_Ele27_CaloIdL_MW_accept, &b_trig_HLT_Mu37_Ele27_CaloIdL_MW_accept);
   fChain->SetBranchAddress("trig_HLT_Mu37_TkMu27_accept", &trig_HLT_Mu37_TkMu27_accept, &b_trig_HLT_Mu37_TkMu27_accept);
   fChain->SetBranchAddress("trig_HLT_DoubleMu4_3_Bs_accept", &trig_HLT_DoubleMu4_3_Bs_accept, &b_trig_HLT_DoubleMu4_3_Bs_accept);
   fChain->SetBranchAddress("trig_HLT_Ele20_WPTight_Gsf_accept", &trig_HLT_Ele20_WPTight_Gsf_accept, &b_trig_HLT_Ele20_WPTight_Gsf_accept);
   fChain->SetBranchAddress("trig_HLT_Ele15_WPLoose_Gsf_accept", &trig_HLT_Ele15_WPLoose_Gsf_accept, &b_trig_HLT_Ele15_WPLoose_Gsf_accept);
   fChain->SetBranchAddress("trig_HLT_Ele17_WPLoose_Gsf_accept", &trig_HLT_Ele17_WPLoose_Gsf_accept, &b_trig_HLT_Ele17_WPLoose_Gsf_accept);
   fChain->SetBranchAddress("trig_HLT_Ele20_WPLoose_Gsf_accept", &trig_HLT_Ele20_WPLoose_Gsf_accept, &b_trig_HLT_Ele20_WPLoose_Gsf_accept);
   fChain->SetBranchAddress("trig_HLT_Ele20_eta2p1_WPLoose_Gsf_accept", &trig_HLT_Ele20_eta2p1_WPLoose_Gsf_accept, &b_trig_HLT_Ele20_eta2p1_WPLoose_Gsf_accept);
   fChain->SetBranchAddress("trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_accept", &trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_accept, &b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_accept);
   fChain->SetBranchAddress("trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltL1sSingleAndDoubleEGor_eta", &trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltL1sSingleAndDoubleEGor_eta, &b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltL1sSingleAndDoubleEGor_eta);
   fChain->SetBranchAddress("trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltL1sSingleAndDoubleEGor_phi", &trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltL1sSingleAndDoubleEGor_phi, &b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltL1sSingleAndDoubleEGor_phi);
   fChain->SetBranchAddress("trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltL1sSingleAndDoubleEGor_et", &trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltL1sSingleAndDoubleEGor_et, &b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltL1sSingleAndDoubleEGor_et);
   fChain->SetBranchAddress("trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEGL1SingleAndDoubleEGOrSingleFilter_eta", &trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEGL1SingleAndDoubleEGOrSingleFilter_eta, &b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEGL1SingleAndDoubleEGOrSingleFilter_eta);
   fChain->SetBranchAddress("trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEGL1SingleAndDoubleEGOrSingleFilter_phi", &trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEGL1SingleAndDoubleEGOrSingleFilter_phi, &b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEGL1SingleAndDoubleEGOrSingleFilter_phi);
   fChain->SetBranchAddress("trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEGL1SingleAndDoubleEGOrSingleFilter_et", &trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEGL1SingleAndDoubleEGOrSingleFilter_et, &b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEGL1SingleAndDoubleEGOrSingleFilter_et);
   fChain->SetBranchAddress("trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEG27L1SingleAndDoubleEGEtFilter_eta", &trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEG27L1SingleAndDoubleEGEtFilter_eta, &b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEG27L1SingleAndDoubleEGEtFilter_eta);
   fChain->SetBranchAddress("trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEG27L1SingleAndDoubleEGEtFilter_phi", &trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEG27L1SingleAndDoubleEGEtFilter_phi, &b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEG27L1SingleAndDoubleEGEtFilter_phi);
   fChain->SetBranchAddress("trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEG27L1SingleAndDoubleEGEtFilter_et", &trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEG27L1SingleAndDoubleEGEtFilter_et, &b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEG27L1SingleAndDoubleEGEtFilter_et);
   fChain->SetBranchAddress("trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightClusterShapeFilter_eta", &trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightClusterShapeFilter_eta, &b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightClusterShapeFilter_eta);
   fChain->SetBranchAddress("trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightClusterShapeFilter_phi", &trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightClusterShapeFilter_phi, &b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightClusterShapeFilter_phi);
   fChain->SetBranchAddress("trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightClusterShapeFilter_et", &trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightClusterShapeFilter_et, &b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightClusterShapeFilter_et);
   fChain->SetBranchAddress("trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightHEFilter_eta", &trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightHEFilter_eta, &b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightHEFilter_eta);
   fChain->SetBranchAddress("trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightHEFilter_phi", &trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightHEFilter_phi, &b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightHEFilter_phi);
   fChain->SetBranchAddress("trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightHEFilter_et", &trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightHEFilter_et, &b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightHEFilter_et);
   fChain->SetBranchAddress("trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightEcalIsoFilter_eta", &trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightEcalIsoFilter_eta, &b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightEcalIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightEcalIsoFilter_phi", &trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightEcalIsoFilter_phi, &b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightEcalIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightEcalIsoFilter_et", &trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightEcalIsoFilter_et, &b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightEcalIsoFilter_et);
   fChain->SetBranchAddress("trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightHcalIsoFilter_eta", &trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightHcalIsoFilter_eta, &b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightHcalIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightHcalIsoFilter_phi", &trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightHcalIsoFilter_phi, &b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightHcalIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightHcalIsoFilter_et", &trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightHcalIsoFilter_et, &b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltEle27L1DoubleEGWPTightHcalIsoFilter_et);
   fChain->SetBranchAddress("trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEG27L1SingleAndDoubleEGEtFilter_eta", &trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEG27L1SingleAndDoubleEGEtFilter_eta, &b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEG27L1SingleAndDoubleEGEtFilter_eta);
   fChain->SetBranchAddress("trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEG27L1SingleAndDoubleEGEtFilter_phi", &trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEG27L1SingleAndDoubleEGEtFilter_phi, &b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEG27L1SingleAndDoubleEGEtFilter_phi);
   fChain->SetBranchAddress("trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEG27L1SingleAndDoubleEGEtFilter_et", &trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEG27L1SingleAndDoubleEGEtFilter_et, &b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEG27L1SingleAndDoubleEGEtFilter_et);
   fChain->SetBranchAddress("trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightClusterShapeFilter_eta", &trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightClusterShapeFilter_eta, &b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightClusterShapeFilter_eta);
   fChain->SetBranchAddress("trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightClusterShapeFilter_phi", &trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightClusterShapeFilter_phi, &b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightClusterShapeFilter_phi);
   fChain->SetBranchAddress("trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightClusterShapeFilter_et", &trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightClusterShapeFilter_et, &b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightClusterShapeFilter_et);
   fChain->SetBranchAddress("trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightHEFilter_eta", &trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightHEFilter_eta, &b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightHEFilter_eta);
   fChain->SetBranchAddress("trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightHEFilter_phi", &trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightHEFilter_phi, &b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightHEFilter_phi);
   fChain->SetBranchAddress("trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightHEFilter_et", &trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightHEFilter_et, &b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightHEFilter_et);
   fChain->SetBranchAddress("trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightEcalIsoFilter_eta", &trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightEcalIsoFilter_eta, &b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightEcalIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightEcalIsoFilter_phi", &trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightEcalIsoFilter_phi, &b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightEcalIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightEcalIsoFilter_et", &trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightEcalIsoFilter_et, &b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightEcalIsoFilter_et);
   fChain->SetBranchAddress("trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightHcalIsoFilter_eta", &trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightHcalIsoFilter_eta, &b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightHcalIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightHcalIsoFilter_phi", &trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightHcalIsoFilter_phi, &b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightHcalIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightHcalIsoFilter_et", &trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightHcalIsoFilter_et, &b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_hltDiEle27L1DoubleEGWPTightHcalIsoFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele27_WPTight_Gsf_accept", &trig_HLT_Ele27_WPTight_Gsf_accept, &b_trig_HLT_Ele27_WPTight_Gsf_accept);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_accept", &trig_HLT_Ele32_WPTight_Gsf_accept, &b_trig_HLT_Ele32_WPTight_Gsf_accept);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_hltL1sSingleEGor_eta", &trig_HLT_Ele32_WPTight_Gsf_hltL1sSingleEGor_eta, &b_trig_HLT_Ele32_WPTight_Gsf_hltL1sSingleEGor_eta);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_hltL1sSingleEGor_phi", &trig_HLT_Ele32_WPTight_Gsf_hltL1sSingleEGor_phi, &b_trig_HLT_Ele32_WPTight_Gsf_hltL1sSingleEGor_phi);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_hltL1sSingleEGor_et", &trig_HLT_Ele32_WPTight_Gsf_hltL1sSingleEGor_et, &b_trig_HLT_Ele32_WPTight_Gsf_hltL1sSingleEGor_et);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_hltEGL1SingleEGOrFilter_eta", &trig_HLT_Ele32_WPTight_Gsf_hltEGL1SingleEGOrFilter_eta, &b_trig_HLT_Ele32_WPTight_Gsf_hltEGL1SingleEGOrFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_hltEGL1SingleEGOrFilter_phi", &trig_HLT_Ele32_WPTight_Gsf_hltEGL1SingleEGOrFilter_phi, &b_trig_HLT_Ele32_WPTight_Gsf_hltEGL1SingleEGOrFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_hltEGL1SingleEGOrFilter_et", &trig_HLT_Ele32_WPTight_Gsf_hltEGL1SingleEGOrFilter_et, &b_trig_HLT_Ele32_WPTight_Gsf_hltEGL1SingleEGOrFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_hltEG32L1SingleEGOrEtFilter_eta", &trig_HLT_Ele32_WPTight_Gsf_hltEG32L1SingleEGOrEtFilter_eta, &b_trig_HLT_Ele32_WPTight_Gsf_hltEG32L1SingleEGOrEtFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_hltEG32L1SingleEGOrEtFilter_phi", &trig_HLT_Ele32_WPTight_Gsf_hltEG32L1SingleEGOrEtFilter_phi, &b_trig_HLT_Ele32_WPTight_Gsf_hltEG32L1SingleEGOrEtFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_hltEG32L1SingleEGOrEtFilter_et", &trig_HLT_Ele32_WPTight_Gsf_hltEG32L1SingleEGOrEtFilter_et, &b_trig_HLT_Ele32_WPTight_Gsf_hltEG32L1SingleEGOrEtFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightClusterShapeFilter_eta", &trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightClusterShapeFilter_eta, &b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightClusterShapeFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightClusterShapeFilter_phi", &trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightClusterShapeFilter_phi, &b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightClusterShapeFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightClusterShapeFilter_et", &trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightClusterShapeFilter_et, &b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightClusterShapeFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightHEFilter_eta", &trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightHEFilter_eta, &b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightHEFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightHEFilter_phi", &trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightHEFilter_phi, &b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightHEFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightHEFilter_et", &trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightHEFilter_et, &b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightHEFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightEcalIsoFilter_eta", &trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightEcalIsoFilter_eta, &b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightEcalIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightEcalIsoFilter_phi", &trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightEcalIsoFilter_phi, &b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightEcalIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightEcalIsoFilter_et", &trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightEcalIsoFilter_et, &b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightEcalIsoFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightHcalIsoFilter_eta", &trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightHcalIsoFilter_eta, &b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightHcalIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightHcalIsoFilter_phi", &trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightHcalIsoFilter_phi, &b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightHcalIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightHcalIsoFilter_et", &trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightHcalIsoFilter_et, &b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightHcalIsoFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightPixelMatchFilter_eta", &trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightPixelMatchFilter_eta, &b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightPixelMatchFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightPixelMatchFilter_phi", &trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightPixelMatchFilter_phi, &b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightPixelMatchFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightPixelMatchFilter_et", &trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightPixelMatchFilter_et, &b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightPixelMatchFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightPMS2Filter_eta", &trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightPMS2Filter_eta, &b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightPMS2Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightPMS2Filter_phi", &trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightPMS2Filter_phi, &b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightPMS2Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightPMS2Filter_et", &trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightPMS2Filter_et, &b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightPMS2Filter_et);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfOneOEMinusOneOPFilter_eta", &trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfOneOEMinusOneOPFilter_eta, &b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfOneOEMinusOneOPFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfOneOEMinusOneOPFilter_phi", &trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfOneOEMinusOneOPFilter_phi, &b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfOneOEMinusOneOPFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfOneOEMinusOneOPFilter_et", &trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfOneOEMinusOneOPFilter_et, &b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfOneOEMinusOneOPFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfMissingHitsFilter_eta", &trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfMissingHitsFilter_eta, &b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfMissingHitsFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfMissingHitsFilter_phi", &trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfMissingHitsFilter_phi, &b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfMissingHitsFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfMissingHitsFilter_et", &trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfMissingHitsFilter_et, &b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfMissingHitsFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfDetaFilter_eta", &trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfDetaFilter_eta, &b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfDetaFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfDetaFilter_phi", &trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfDetaFilter_phi, &b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfDetaFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfDetaFilter_et", &trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfDetaFilter_et, &b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfDetaFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfDphiFilter_eta", &trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfDphiFilter_eta, &b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfDphiFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfDphiFilter_phi", &trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfDphiFilter_phi, &b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfDphiFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfDphiFilter_et", &trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfDphiFilter_et, &b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfDphiFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfTrackIsoFilter_eta", &trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfTrackIsoFilter_eta, &b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfTrackIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfTrackIsoFilter_phi", &trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfTrackIsoFilter_phi, &b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfTrackIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfTrackIsoFilter_et", &trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfTrackIsoFilter_et, &b_trig_HLT_Ele32_WPTight_Gsf_hltEle32WPTightGsfTrackIsoFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_accept", &trig_HLT_Ele35_WPTight_Gsf_accept, &b_trig_HLT_Ele35_WPTight_Gsf_accept);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_hltL1sSingleEGor_eta", &trig_HLT_Ele35_WPTight_Gsf_hltL1sSingleEGor_eta, &b_trig_HLT_Ele35_WPTight_Gsf_hltL1sSingleEGor_eta);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_hltL1sSingleEGor_phi", &trig_HLT_Ele35_WPTight_Gsf_hltL1sSingleEGor_phi, &b_trig_HLT_Ele35_WPTight_Gsf_hltL1sSingleEGor_phi);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_hltL1sSingleEGor_et", &trig_HLT_Ele35_WPTight_Gsf_hltL1sSingleEGor_et, &b_trig_HLT_Ele35_WPTight_Gsf_hltL1sSingleEGor_et);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_hltEGL1SingleEGOrFilter_eta", &trig_HLT_Ele35_WPTight_Gsf_hltEGL1SingleEGOrFilter_eta, &b_trig_HLT_Ele35_WPTight_Gsf_hltEGL1SingleEGOrFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_hltEGL1SingleEGOrFilter_phi", &trig_HLT_Ele35_WPTight_Gsf_hltEGL1SingleEGOrFilter_phi, &b_trig_HLT_Ele35_WPTight_Gsf_hltEGL1SingleEGOrFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_hltEGL1SingleEGOrFilter_et", &trig_HLT_Ele35_WPTight_Gsf_hltEGL1SingleEGOrFilter_et, &b_trig_HLT_Ele35_WPTight_Gsf_hltEGL1SingleEGOrFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_hltEG35L1SingleEGOrEtFilter_eta", &trig_HLT_Ele35_WPTight_Gsf_hltEG35L1SingleEGOrEtFilter_eta, &b_trig_HLT_Ele35_WPTight_Gsf_hltEG35L1SingleEGOrEtFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_hltEG35L1SingleEGOrEtFilter_phi", &trig_HLT_Ele35_WPTight_Gsf_hltEG35L1SingleEGOrEtFilter_phi, &b_trig_HLT_Ele35_WPTight_Gsf_hltEG35L1SingleEGOrEtFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_hltEG35L1SingleEGOrEtFilter_et", &trig_HLT_Ele35_WPTight_Gsf_hltEG35L1SingleEGOrEtFilter_et, &b_trig_HLT_Ele35_WPTight_Gsf_hltEG35L1SingleEGOrEtFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightClusterShapeFilter_eta", &trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightClusterShapeFilter_eta, &b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightClusterShapeFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightClusterShapeFilter_phi", &trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightClusterShapeFilter_phi, &b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightClusterShapeFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightClusterShapeFilter_et", &trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightClusterShapeFilter_et, &b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightClusterShapeFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightHEFilter_eta", &trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightHEFilter_eta, &b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightHEFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightHEFilter_phi", &trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightHEFilter_phi, &b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightHEFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightHEFilter_et", &trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightHEFilter_et, &b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightHEFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightEcalIsoFilter_eta", &trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightEcalIsoFilter_eta, &b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightEcalIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightEcalIsoFilter_phi", &trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightEcalIsoFilter_phi, &b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightEcalIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightEcalIsoFilter_et", &trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightEcalIsoFilter_et, &b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightEcalIsoFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightHcalIsoFilter_eta", &trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightHcalIsoFilter_eta, &b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightHcalIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightHcalIsoFilter_phi", &trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightHcalIsoFilter_phi, &b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightHcalIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightHcalIsoFilter_et", &trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightHcalIsoFilter_et, &b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightHcalIsoFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightPixelMatchFilter_eta", &trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightPixelMatchFilter_eta, &b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightPixelMatchFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightPixelMatchFilter_phi", &trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightPixelMatchFilter_phi, &b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightPixelMatchFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightPixelMatchFilter_et", &trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightPixelMatchFilter_et, &b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightPixelMatchFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightPMS2Filter_eta", &trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightPMS2Filter_eta, &b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightPMS2Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightPMS2Filter_phi", &trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightPMS2Filter_phi, &b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightPMS2Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightPMS2Filter_et", &trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightPMS2Filter_et, &b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightPMS2Filter_et);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfOneOEMinusOneOPFilter_eta", &trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfOneOEMinusOneOPFilter_eta, &b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfOneOEMinusOneOPFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfOneOEMinusOneOPFilter_phi", &trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfOneOEMinusOneOPFilter_phi, &b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfOneOEMinusOneOPFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfOneOEMinusOneOPFilter_et", &trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfOneOEMinusOneOPFilter_et, &b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfOneOEMinusOneOPFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfMissingHitsFilter_eta", &trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfMissingHitsFilter_eta, &b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfMissingHitsFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfMissingHitsFilter_phi", &trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfMissingHitsFilter_phi, &b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfMissingHitsFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfMissingHitsFilter_et", &trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfMissingHitsFilter_et, &b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfMissingHitsFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfDetaFilter_eta", &trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfDetaFilter_eta, &b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfDetaFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfDetaFilter_phi", &trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfDetaFilter_phi, &b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfDetaFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfDetaFilter_et", &trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfDetaFilter_et, &b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfDetaFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfDphiFilter_eta", &trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfDphiFilter_eta, &b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfDphiFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfDphiFilter_phi", &trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfDphiFilter_phi, &b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfDphiFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfDphiFilter_et", &trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfDphiFilter_et, &b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfDphiFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfTrackIsoFilter_eta", &trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfTrackIsoFilter_eta, &b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfTrackIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfTrackIsoFilter_phi", &trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfTrackIsoFilter_phi, &b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfTrackIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfTrackIsoFilter_et", &trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfTrackIsoFilter_et, &b_trig_HLT_Ele35_WPTight_Gsf_hltEle35noerWPTightGsfTrackIsoFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_L1EGMT_accept", &trig_HLT_Ele35_WPTight_Gsf_L1EGMT_accept, &b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_accept);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltL1sAlCaSingleEle_eta", &trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltL1sAlCaSingleEle_eta, &b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltL1sAlCaSingleEle_eta);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltL1sAlCaSingleEle_phi", &trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltL1sAlCaSingleEle_phi, &b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltL1sAlCaSingleEle_phi);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltL1sAlCaSingleEle_et", &trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltL1sAlCaSingleEle_et, &b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltL1sAlCaSingleEle_et);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTHEFilter_eta", &trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTHEFilter_eta, &b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTHEFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTHEFilter_phi", &trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTHEFilter_phi, &b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTHEFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTHEFilter_et", &trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTHEFilter_et, &b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTHEFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTEcalIsoFilter_eta", &trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTEcalIsoFilter_eta, &b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTEcalIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTEcalIsoFilter_phi", &trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTEcalIsoFilter_phi, &b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTEcalIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTEcalIsoFilter_et", &trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTEcalIsoFilter_et, &b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTEcalIsoFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTHcalIsoFilter_eta", &trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTHcalIsoFilter_eta, &b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTHcalIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTHcalIsoFilter_phi", &trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTHcalIsoFilter_phi, &b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTHcalIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTHcalIsoFilter_et", &trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTHcalIsoFilter_et, &b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTHcalIsoFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTPixelMatchFilter_eta", &trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTPixelMatchFilter_eta, &b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTPixelMatchFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTPixelMatchFilter_phi", &trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTPixelMatchFilter_phi, &b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTPixelMatchFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTPixelMatchFilter_et", &trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTPixelMatchFilter_et, &b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTPixelMatchFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTPMS2Filter_eta", &trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTPMS2Filter_eta, &b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTPMS2Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTPMS2Filter_phi", &trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTPMS2Filter_phi, &b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTPMS2Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTPMS2Filter_et", &trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTPMS2Filter_et, &b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTPMS2Filter_et);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTDetaFilter_eta", &trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTDetaFilter_eta, &b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTDetaFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTDetaFilter_phi", &trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTDetaFilter_phi, &b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTDetaFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTDetaFilter_et", &trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTDetaFilter_et, &b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTDetaFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTDphiFilter_eta", &trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTDphiFilter_eta, &b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTDphiFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTDphiFilter_phi", &trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTDphiFilter_phi, &b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTDphiFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTDphiFilter_et", &trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTDphiFilter_et, &b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTDphiFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTTrackIsoFilter_eta", &trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTTrackIsoFilter_eta, &b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTTrackIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTTrackIsoFilter_phi", &trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTTrackIsoFilter_phi, &b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTTrackIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTTrackIsoFilter_et", &trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTTrackIsoFilter_et, &b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_hltSingleEle35WPTightGsfL1EGMTTrackIsoFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele38_WPTight_Gsf_accept", &trig_HLT_Ele38_WPTight_Gsf_accept, &b_trig_HLT_Ele38_WPTight_Gsf_accept);
   fChain->SetBranchAddress("trig_HLT_Ele40_WPTight_Gsf_accept", &trig_HLT_Ele40_WPTight_Gsf_accept, &b_trig_HLT_Ele40_WPTight_Gsf_accept);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_accept", &trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_accept, &b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_accept);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltL1sSingleAndDoubleEGor_eta", &trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltL1sSingleAndDoubleEGor_eta, &b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltL1sSingleAndDoubleEGor_eta);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltL1sSingleAndDoubleEGor_phi", &trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltL1sSingleAndDoubleEGor_phi, &b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltL1sSingleAndDoubleEGor_phi);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltL1sSingleAndDoubleEGor_et", &trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltL1sSingleAndDoubleEGor_et, &b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltL1sSingleAndDoubleEGor_et);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEGL1SingleAndDoubleEGOrSingleFilter_eta", &trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEGL1SingleAndDoubleEGOrSingleFilter_eta, &b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEGL1SingleAndDoubleEGOrSingleFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEGL1SingleAndDoubleEGOrSingleFilter_phi", &trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEGL1SingleAndDoubleEGOrSingleFilter_phi, &b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEGL1SingleAndDoubleEGOrSingleFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEGL1SingleAndDoubleEGOrSingleFilter_et", &trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEGL1SingleAndDoubleEGOrSingleFilter_et, &b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEGL1SingleAndDoubleEGOrSingleFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEG32L1SingleAndDoubleEGEtFilter_eta", &trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEG32L1SingleAndDoubleEGEtFilter_eta, &b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEG32L1SingleAndDoubleEGEtFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEG32L1SingleAndDoubleEGEtFilter_phi", &trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEG32L1SingleAndDoubleEGEtFilter_phi, &b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEG32L1SingleAndDoubleEGEtFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEG32L1SingleAndDoubleEGEtFilter_et", &trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEG32L1SingleAndDoubleEGEtFilter_et, &b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEG32L1SingleAndDoubleEGEtFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightClusterShapeFilter_eta", &trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightClusterShapeFilter_eta, &b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightClusterShapeFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightClusterShapeFilter_phi", &trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightClusterShapeFilter_phi, &b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightClusterShapeFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightClusterShapeFilter_et", &trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightClusterShapeFilter_et, &b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightClusterShapeFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightHEFilter_eta", &trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightHEFilter_eta, &b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightHEFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightHEFilter_phi", &trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightHEFilter_phi, &b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightHEFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightHEFilter_et", &trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightHEFilter_et, &b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightHEFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightEcalIsoFilter_eta", &trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightEcalIsoFilter_eta, &b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightEcalIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightEcalIsoFilter_phi", &trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightEcalIsoFilter_phi, &b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightEcalIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightEcalIsoFilter_et", &trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightEcalIsoFilter_et, &b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightEcalIsoFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightHcalIsoFilter_eta", &trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightHcalIsoFilter_eta, &b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightHcalIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightHcalIsoFilter_phi", &trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightHcalIsoFilter_phi, &b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightHcalIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightHcalIsoFilter_et", &trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightHcalIsoFilter_et, &b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightHcalIsoFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightPixelMatchFilter_eta", &trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightPixelMatchFilter_eta, &b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightPixelMatchFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightPixelMatchFilter_phi", &trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightPixelMatchFilter_phi, &b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightPixelMatchFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightPixelMatchFilter_et", &trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightPixelMatchFilter_et, &b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightPixelMatchFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightPMS2Filter_eta", &trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightPMS2Filter_eta, &b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightPMS2Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightPMS2Filter_phi", &trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightPMS2Filter_phi, &b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightPMS2Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightPMS2Filter_et", &trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightPMS2Filter_et, &b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightPMS2Filter_et);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfOneOEMinusOneOPFilter_eta", &trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfOneOEMinusOneOPFilter_eta, &b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfOneOEMinusOneOPFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfOneOEMinusOneOPFilter_phi", &trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfOneOEMinusOneOPFilter_phi, &b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfOneOEMinusOneOPFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfOneOEMinusOneOPFilter_et", &trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfOneOEMinusOneOPFilter_et, &b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfOneOEMinusOneOPFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfMissingHitsFilter_eta", &trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfMissingHitsFilter_eta, &b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfMissingHitsFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfMissingHitsFilter_phi", &trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfMissingHitsFilter_phi, &b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfMissingHitsFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfMissingHitsFilter_et", &trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfMissingHitsFilter_et, &b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfMissingHitsFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfDetaFilter_eta", &trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfDetaFilter_eta, &b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfDetaFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfDetaFilter_phi", &trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfDetaFilter_phi, &b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfDetaFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfDetaFilter_et", &trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfDetaFilter_et, &b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfDetaFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfDphiFilter_eta", &trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfDphiFilter_eta, &b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfDphiFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfDphiFilter_phi", &trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfDphiFilter_phi, &b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfDphiFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfDphiFilter_et", &trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfDphiFilter_et, &b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfDphiFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfTrackIsoFilter_eta", &trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfTrackIsoFilter_eta, &b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfTrackIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfTrackIsoFilter_phi", &trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfTrackIsoFilter_phi, &b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfTrackIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfTrackIsoFilter_et", &trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfTrackIsoFilter_et, &b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_hltEle32L1DoubleEGWPTightGsfTrackIsoFilter_et);
   fChain->SetBranchAddress("trig_HLT_IsoMu20_accept", &trig_HLT_IsoMu20_accept, &b_trig_HLT_IsoMu20_accept);
   fChain->SetBranchAddress("trig_HLT_IsoMu24_accept", &trig_HLT_IsoMu24_accept, &b_trig_HLT_IsoMu24_accept);
   fChain->SetBranchAddress("trig_HLT_IsoMu24_eta2p1_accept", &trig_HLT_IsoMu24_eta2p1_accept, &b_trig_HLT_IsoMu24_eta2p1_accept);
   fChain->SetBranchAddress("trig_HLT_IsoMu27_accept", &trig_HLT_IsoMu27_accept, &b_trig_HLT_IsoMu27_accept);
   fChain->SetBranchAddress("trig_HLT_IsoMu30_accept", &trig_HLT_IsoMu30_accept, &b_trig_HLT_IsoMu30_accept);
   fChain->SetBranchAddress("trig_HLT_L1SingleMu18_accept", &trig_HLT_L1SingleMu18_accept, &b_trig_HLT_L1SingleMu18_accept);
   fChain->SetBranchAddress("trig_HLT_L1SingleMu25_accept", &trig_HLT_L1SingleMu25_accept, &b_trig_HLT_L1SingleMu25_accept);
   fChain->SetBranchAddress("trig_HLT_L2Mu10_accept", &trig_HLT_L2Mu10_accept, &b_trig_HLT_L2Mu10_accept);
   fChain->SetBranchAddress("trig_HLT_L2Mu10_NoVertex_NoBPTX3BX_accept", &trig_HLT_L2Mu10_NoVertex_NoBPTX3BX_accept, &b_trig_HLT_L2Mu10_NoVertex_NoBPTX3BX_accept);
   fChain->SetBranchAddress("trig_HLT_L2Mu10_NoVertex_NoBPTX_accept", &trig_HLT_L2Mu10_NoVertex_NoBPTX_accept, &b_trig_HLT_L2Mu10_NoVertex_NoBPTX_accept);
   fChain->SetBranchAddress("trig_HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX_accept", &trig_HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX_accept, &b_trig_HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX_accept);
   fChain->SetBranchAddress("trig_HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX_accept", &trig_HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX_accept, &b_trig_HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX_accept);
   fChain->SetBranchAddress("trig_HLT_L2Mu50_accept", &trig_HLT_L2Mu50_accept, &b_trig_HLT_L2Mu50_accept);
   fChain->SetBranchAddress("trig_HLT_L2Mu23NoVtx_2Cha_accept", &trig_HLT_L2Mu23NoVtx_2Cha_accept, &b_trig_HLT_L2Mu23NoVtx_2Cha_accept);
   fChain->SetBranchAddress("trig_HLT_L2Mu23NoVtx_2Cha_CosmicSeed_accept", &trig_HLT_L2Mu23NoVtx_2Cha_CosmicSeed_accept, &b_trig_HLT_L2Mu23NoVtx_2Cha_CosmicSeed_accept);
   fChain->SetBranchAddress("trig_HLT_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4_accept", &trig_HLT_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4_accept, &b_trig_HLT_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4_accept);
   fChain->SetBranchAddress("trig_HLT_DoubleL2Mu30NoVtx_2Cha_Eta2p4_accept", &trig_HLT_DoubleL2Mu30NoVtx_2Cha_Eta2p4_accept, &b_trig_HLT_DoubleL2Mu30NoVtx_2Cha_Eta2p4_accept);
   fChain->SetBranchAddress("trig_HLT_DoubleL2Mu50_accept", &trig_HLT_DoubleL2Mu50_accept, &b_trig_HLT_DoubleL2Mu50_accept);
   fChain->SetBranchAddress("trig_HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed_accept", &trig_HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed_accept, &b_trig_HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed_accept);
   fChain->SetBranchAddress("trig_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_accept", &trig_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_accept, &b_trig_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_accept);
   fChain->SetBranchAddress("trig_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4_accept", &trig_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4_accept, &b_trig_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4_accept);
   fChain->SetBranchAddress("trig_HLT_DoubleL2Mu23NoVtx_2Cha_accept", &trig_HLT_DoubleL2Mu23NoVtx_2Cha_accept, &b_trig_HLT_DoubleL2Mu23NoVtx_2Cha_accept);
   fChain->SetBranchAddress("trig_HLT_DoubleL2Mu25NoVtx_2Cha_accept", &trig_HLT_DoubleL2Mu25NoVtx_2Cha_accept, &b_trig_HLT_DoubleL2Mu25NoVtx_2Cha_accept);
   fChain->SetBranchAddress("trig_HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4_accept", &trig_HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4_accept, &b_trig_HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4_accept);
   fChain->SetBranchAddress("trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_accept", &trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_accept, &b_trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_accept);
   fChain->SetBranchAddress("trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_accept", &trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_accept, &b_trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_accept);
   fChain->SetBranchAddress("trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_accept", &trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_accept, &b_trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_accept);
   fChain->SetBranchAddress("trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_accept", &trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_accept, &b_trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_accept);
   fChain->SetBranchAddress("trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_accept", &trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_accept, &b_trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_accept);
   fChain->SetBranchAddress("trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8_accept", &trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8_accept, &b_trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8_accept);
   fChain->SetBranchAddress("trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_accept", &trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_accept, &b_trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_accept);
   fChain->SetBranchAddress("trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8_accept", &trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8_accept, &b_trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8_accept);
   fChain->SetBranchAddress("trig_HLT_Mu25_TkMu0_Onia_accept", &trig_HLT_Mu25_TkMu0_Onia_accept, &b_trig_HLT_Mu25_TkMu0_Onia_accept);
   fChain->SetBranchAddress("trig_HLT_Mu30_TkMu0_Onia_accept", &trig_HLT_Mu30_TkMu0_Onia_accept, &b_trig_HLT_Mu30_TkMu0_Onia_accept);
   fChain->SetBranchAddress("trig_HLT_Mu20_TkMu0_Phi_accept", &trig_HLT_Mu20_TkMu0_Phi_accept, &b_trig_HLT_Mu20_TkMu0_Phi_accept);
   fChain->SetBranchAddress("trig_HLT_Mu25_TkMu0_Phi_accept", &trig_HLT_Mu25_TkMu0_Phi_accept, &b_trig_HLT_Mu25_TkMu0_Phi_accept);
   fChain->SetBranchAddress("trig_HLT_Mu12_accept", &trig_HLT_Mu12_accept, &b_trig_HLT_Mu12_accept);
   fChain->SetBranchAddress("trig_HLT_Mu15_accept", &trig_HLT_Mu15_accept, &b_trig_HLT_Mu15_accept);
   fChain->SetBranchAddress("trig_HLT_Mu20_accept", &trig_HLT_Mu20_accept, &b_trig_HLT_Mu20_accept);
   fChain->SetBranchAddress("trig_HLT_Mu27_accept", &trig_HLT_Mu27_accept, &b_trig_HLT_Mu27_accept);
   fChain->SetBranchAddress("trig_HLT_Mu50_accept", &trig_HLT_Mu50_accept, &b_trig_HLT_Mu50_accept);
   fChain->SetBranchAddress("trig_HLT_Mu55_accept", &trig_HLT_Mu55_accept, &b_trig_HLT_Mu55_accept);
   fChain->SetBranchAddress("trig_HLT_OldMu100_accept", &trig_HLT_OldMu100_accept, &b_trig_HLT_OldMu100_accept);
   fChain->SetBranchAddress("trig_HLT_TkMu100_accept", &trig_HLT_TkMu100_accept, &b_trig_HLT_TkMu100_accept);
   fChain->SetBranchAddress("trig_HLT_Photon300_NoHE_accept", &trig_HLT_Photon300_NoHE_accept, &b_trig_HLT_Photon300_NoHE_accept);
   fChain->SetBranchAddress("trig_HLT_Photon300_NoHE_hltL1sSingleEGNonIsoOrWithJetAndTau_eta", &trig_HLT_Photon300_NoHE_hltL1sSingleEGNonIsoOrWithJetAndTau_eta, &b_trig_HLT_Photon300_NoHE_hltL1sSingleEGNonIsoOrWithJetAndTau_eta);
   fChain->SetBranchAddress("trig_HLT_Photon300_NoHE_hltL1sSingleEGNonIsoOrWithJetAndTau_phi", &trig_HLT_Photon300_NoHE_hltL1sSingleEGNonIsoOrWithJetAndTau_phi, &b_trig_HLT_Photon300_NoHE_hltL1sSingleEGNonIsoOrWithJetAndTau_phi);
   fChain->SetBranchAddress("trig_HLT_Photon300_NoHE_hltL1sSingleEGNonIsoOrWithJetAndTau_et", &trig_HLT_Photon300_NoHE_hltL1sSingleEGNonIsoOrWithJetAndTau_et, &b_trig_HLT_Photon300_NoHE_hltL1sSingleEGNonIsoOrWithJetAndTau_et);
   fChain->SetBranchAddress("trig_HLT_Photon300_NoHE_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_eta", &trig_HLT_Photon300_NoHE_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_eta, &b_trig_HLT_Photon300_NoHE_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon300_NoHE_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_phi", &trig_HLT_Photon300_NoHE_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_phi, &b_trig_HLT_Photon300_NoHE_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon300_NoHE_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_et", &trig_HLT_Photon300_NoHE_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_et, &b_trig_HLT_Photon300_NoHE_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon300_NoHE_hltEG300erEtFilter_eta", &trig_HLT_Photon300_NoHE_hltEG300erEtFilter_eta, &b_trig_HLT_Photon300_NoHE_hltEG300erEtFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon300_NoHE_hltEG300erEtFilter_phi", &trig_HLT_Photon300_NoHE_hltEG300erEtFilter_phi, &b_trig_HLT_Photon300_NoHE_hltEG300erEtFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon300_NoHE_hltEG300erEtFilter_et", &trig_HLT_Photon300_NoHE_hltEG300erEtFilter_et, &b_trig_HLT_Photon300_NoHE_hltEG300erEtFilter_et);
   fChain->SetBranchAddress("trig_HLT_Mu8_TrkIsoVVL_accept", &trig_HLT_Mu8_TrkIsoVVL_accept, &b_trig_HLT_Mu8_TrkIsoVVL_accept);
   fChain->SetBranchAddress("trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept", &trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept, &b_trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept);
   fChain->SetBranchAddress("trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept", &trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept, &b_trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept);
   fChain->SetBranchAddress("trig_HLT_Mu17_TrkIsoVVL_accept", &trig_HLT_Mu17_TrkIsoVVL_accept, &b_trig_HLT_Mu17_TrkIsoVVL_accept);
   fChain->SetBranchAddress("trig_HLT_Mu19_TrkIsoVVL_accept", &trig_HLT_Mu19_TrkIsoVVL_accept, &b_trig_HLT_Mu19_TrkIsoVVL_accept);
   fChain->SetBranchAddress("trig_HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL_accept", &trig_HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL_accept, &b_trig_HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL_accept);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltL1sSingleAndDoubleEG_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltL1sSingleAndDoubleEG_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltL1sSingleAndDoubleEG_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltL1sSingleAndDoubleEG_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltL1sSingleAndDoubleEG_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltL1sSingleAndDoubleEG_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltL1sSingleAndDoubleEG_et", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltL1sSingleAndDoubleEG_et, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltL1sSingleAndDoubleEG_et);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEGL1SingleAndDoubleEGOrPairFilter_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEGL1SingleAndDoubleEGOrPairFilter_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEGL1SingleAndDoubleEGOrPairFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEGL1SingleAndDoubleEGOrPairFilter_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEGL1SingleAndDoubleEGOrPairFilter_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEGL1SingleAndDoubleEGOrPairFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEGL1SingleAndDoubleEGOrPairFilter_et", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEGL1SingleAndDoubleEGOrPairFilter_et, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEGL1SingleAndDoubleEGOrPairFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg1Filter_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg1Filter_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg1Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg1Filter_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg1Filter_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg1Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg1Filter_et", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg1Filter_et, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg1Filter_et);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg2Filter_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg2Filter_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg2Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg2Filter_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg2Filter_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg2Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg2Filter_et", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg2Filter_et, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg2Filter_et);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg1Filter_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg1Filter_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg1Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg1Filter_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg1Filter_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg1Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg1Filter_et", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg1Filter_et, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg1Filter_et);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg2Filter_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg2Filter_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg2Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg2Filter_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg2Filter_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg2Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg2Filter_et", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg2Filter_et, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg2Filter_et);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg1Filter_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg1Filter_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg1Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg1Filter_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg1Filter_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg1Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg1Filter_et", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg1Filter_et, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg1Filter_et);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg2Filter_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg2Filter_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg2Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg2Filter_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg2Filter_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg2Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg2Filter_et", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg2Filter_et, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg2Filter_et);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg1Filter_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg1Filter_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg1Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg1Filter_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg1Filter_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg1Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg1Filter_et", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg1Filter_et, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg1Filter_et);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg2Filter_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg2Filter_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg2Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg2Filter_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg2Filter_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg2Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg2Filter_et", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg2Filter_et, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg2Filter_et);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg1Filter_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg1Filter_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg1Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg1Filter_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg1Filter_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg1Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg1Filter_et", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg1Filter_et, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg1Filter_et);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg2Filter_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg2Filter_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg2Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg2Filter_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg2Filter_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg2Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg2Filter_et", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg2Filter_et, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg2Filter_et);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg1Filter_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg1Filter_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg1Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg1Filter_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg1Filter_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg1Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg1Filter_et", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg1Filter_et, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg1Filter_et);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg2Filter_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg2Filter_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg2Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg2Filter_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg2Filter_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg2Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg2Filter_et", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg2Filter_et, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg2Filter_et);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg1Filter_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg1Filter_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg1Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg1Filter_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg1Filter_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg1Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg1Filter_et", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg1Filter_et, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg1Filter_et);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg2Filter_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg2Filter_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg2Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg2Filter_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg2Filter_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg2Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg2Filter_et", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg2Filter_et, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg2Filter_et);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg1Filter_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg1Filter_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg1Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg1Filter_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg1Filter_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg1Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg1Filter_et", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg1Filter_et, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg1Filter_et);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg2Filter_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg2Filter_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg2Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg2Filter_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg2Filter_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg2Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg2Filter_et", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg2Filter_et, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg2Filter_et);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg1Filter_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg1Filter_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg1Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg1Filter_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg1Filter_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg1Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg1Filter_et", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg1Filter_et, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg1Filter_et);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg2Filter_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg2Filter_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg2Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg2Filter_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg2Filter_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg2Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg2Filter_et", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg2Filter_et, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg2Filter_et);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_et", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_et, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_et);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_et", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_et, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_et);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter_et", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter_et, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_accept", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_accept, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_accept);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltL1sSingleAndDoubleEG_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltL1sSingleAndDoubleEG_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltL1sSingleAndDoubleEG_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltL1sSingleAndDoubleEG_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltL1sSingleAndDoubleEG_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltL1sSingleAndDoubleEG_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltL1sSingleAndDoubleEG_et", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltL1sSingleAndDoubleEG_et, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltL1sSingleAndDoubleEG_et);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEGL1SingleAndDoubleEGOrPairFilter_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEGL1SingleAndDoubleEGOrPairFilter_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEGL1SingleAndDoubleEGOrPairFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEGL1SingleAndDoubleEGOrPairFilter_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEGL1SingleAndDoubleEGOrPairFilter_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEGL1SingleAndDoubleEGOrPairFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEGL1SingleAndDoubleEGOrPairFilter_et", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEGL1SingleAndDoubleEGOrPairFilter_et, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEGL1SingleAndDoubleEGOrPairFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg1Filter_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg1Filter_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg1Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg1Filter_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg1Filter_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg1Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg1Filter_et", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg1Filter_et, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg1Filter_et);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg2Filter_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg2Filter_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg2Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg2Filter_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg2Filter_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg2Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg2Filter_et", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg2Filter_et, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg2Filter_et);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg1Filter_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg1Filter_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg1Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg1Filter_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg1Filter_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg1Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg1Filter_et", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg1Filter_et, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg1Filter_et);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg2Filter_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg2Filter_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg2Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg2Filter_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg2Filter_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg2Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg2Filter_et", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg2Filter_et, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLClusterShapeLeg2Filter_et);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg1Filter_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg1Filter_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg1Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg1Filter_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg1Filter_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg1Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg1Filter_et", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg1Filter_et, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg1Filter_et);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg2Filter_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg2Filter_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg2Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg2Filter_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg2Filter_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg2Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg2Filter_et", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg2Filter_et, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHELeg2Filter_et);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg1Filter_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg1Filter_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg1Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg1Filter_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg1Filter_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg1Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg1Filter_et", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg1Filter_et, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg1Filter_et);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg2Filter_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg2Filter_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg2Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg2Filter_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg2Filter_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg2Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg2Filter_et", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg2Filter_et, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLEcalIsoLeg2Filter_et);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg1Filter_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg1Filter_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg1Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg1Filter_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg1Filter_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg1Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg1Filter_et", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg1Filter_et, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg1Filter_et);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg2Filter_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg2Filter_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg2Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg2Filter_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg2Filter_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg2Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg2Filter_et", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg2Filter_et, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLHcalIsoLeg2Filter_et);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg1Filter_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg1Filter_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg1Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg1Filter_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg1Filter_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg1Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg1Filter_et", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg1Filter_et, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg1Filter_et);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg2Filter_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg2Filter_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg2Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg2Filter_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg2Filter_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg2Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg2Filter_et", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg2Filter_et, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLPixelMatchLeg2Filter_et);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg1Filter_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg1Filter_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg1Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg1Filter_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg1Filter_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg1Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg1Filter_et", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg1Filter_et, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg1Filter_et);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg2Filter_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg2Filter_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg2Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg2Filter_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg2Filter_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg2Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg2Filter_et", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg2Filter_et, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLOneOEMinusOneOPLeg2Filter_et);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg1Filter_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg1Filter_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg1Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg1Filter_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg1Filter_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg1Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg1Filter_et", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg1Filter_et, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg1Filter_et);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg2Filter_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg2Filter_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg2Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg2Filter_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg2Filter_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg2Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg2Filter_et", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg2Filter_et, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDetaLeg2Filter_et);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg1Filter_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg1Filter_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg1Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg1Filter_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg1Filter_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg1Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg1Filter_et", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg1Filter_et, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg1Filter_et);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg2Filter_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg2Filter_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg2Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg2Filter_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg2Filter_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg2Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg2Filter_et", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg2Filter_et, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLDphiLeg2Filter_et);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_et", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_et, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_et);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_et", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_et, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_et);
   fChain->SetBranchAddress("trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept", &trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept, &b_trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept);
   fChain->SetBranchAddress("trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept", &trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept, &b_trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept);
   fChain->SetBranchAddress("trig_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept", &trig_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept, &b_trig_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept);
   fChain->SetBranchAddress("trig_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept", &trig_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept, &b_trig_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept);
   fChain->SetBranchAddress("trig_HLT_Photon20_accept", &trig_HLT_Photon20_accept, &b_trig_HLT_Photon20_accept);
   fChain->SetBranchAddress("trig_HLT_Photon20_hltL1sSingleEG15_eta", &trig_HLT_Photon20_hltL1sSingleEG15_eta, &b_trig_HLT_Photon20_hltL1sSingleEG15_eta);
   fChain->SetBranchAddress("trig_HLT_Photon20_hltL1sSingleEG15_phi", &trig_HLT_Photon20_hltL1sSingleEG15_phi, &b_trig_HLT_Photon20_hltL1sSingleEG15_phi);
   fChain->SetBranchAddress("trig_HLT_Photon20_hltL1sSingleEG15_et", &trig_HLT_Photon20_hltL1sSingleEG15_et, &b_trig_HLT_Photon20_hltL1sSingleEG15_et);
   fChain->SetBranchAddress("trig_HLT_Photon20_hltEGL1SingleEG15Filter_eta", &trig_HLT_Photon20_hltEGL1SingleEG15Filter_eta, &b_trig_HLT_Photon20_hltEGL1SingleEG15Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon20_hltEGL1SingleEG15Filter_phi", &trig_HLT_Photon20_hltEGL1SingleEG15Filter_phi, &b_trig_HLT_Photon20_hltEGL1SingleEG15Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon20_hltEGL1SingleEG15Filter_et", &trig_HLT_Photon20_hltEGL1SingleEG15Filter_et, &b_trig_HLT_Photon20_hltEGL1SingleEG15Filter_et);
   fChain->SetBranchAddress("trig_HLT_Photon20_hltEG20L1EG15EtFilter_eta", &trig_HLT_Photon20_hltEG20L1EG15EtFilter_eta, &b_trig_HLT_Photon20_hltEG20L1EG15EtFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon20_hltEG20L1EG15EtFilter_phi", &trig_HLT_Photon20_hltEG20L1EG15EtFilter_phi, &b_trig_HLT_Photon20_hltEG20L1EG15EtFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon20_hltEG20L1EG15EtFilter_et", &trig_HLT_Photon20_hltEG20L1EG15EtFilter_et, &b_trig_HLT_Photon20_hltEG20L1EG15EtFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon20_hltEG20L1EG10OR15HEFilter_eta", &trig_HLT_Photon20_hltEG20L1EG10OR15HEFilter_eta, &b_trig_HLT_Photon20_hltEG20L1EG10OR15HEFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon20_hltEG20L1EG10OR15HEFilter_phi", &trig_HLT_Photon20_hltEG20L1EG10OR15HEFilter_phi, &b_trig_HLT_Photon20_hltEG20L1EG10OR15HEFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon20_hltEG20L1EG10OR15HEFilter_et", &trig_HLT_Photon20_hltEG20L1EG10OR15HEFilter_et, &b_trig_HLT_Photon20_hltEG20L1EG10OR15HEFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon33_accept", &trig_HLT_Photon33_accept, &b_trig_HLT_Photon33_accept);
   fChain->SetBranchAddress("trig_HLT_Photon33_hltL1sSingleEG26_eta", &trig_HLT_Photon33_hltL1sSingleEG26_eta, &b_trig_HLT_Photon33_hltL1sSingleEG26_eta);
   fChain->SetBranchAddress("trig_HLT_Photon33_hltL1sSingleEG26_phi", &trig_HLT_Photon33_hltL1sSingleEG26_phi, &b_trig_HLT_Photon33_hltL1sSingleEG26_phi);
   fChain->SetBranchAddress("trig_HLT_Photon33_hltL1sSingleEG26_et", &trig_HLT_Photon33_hltL1sSingleEG26_et, &b_trig_HLT_Photon33_hltL1sSingleEG26_et);
   fChain->SetBranchAddress("trig_HLT_Photon33_hltEGL1SingleEG26Filter_eta", &trig_HLT_Photon33_hltEGL1SingleEG26Filter_eta, &b_trig_HLT_Photon33_hltEGL1SingleEG26Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon33_hltEGL1SingleEG26Filter_phi", &trig_HLT_Photon33_hltEGL1SingleEG26Filter_phi, &b_trig_HLT_Photon33_hltEGL1SingleEG26Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon33_hltEGL1SingleEG26Filter_et", &trig_HLT_Photon33_hltEGL1SingleEG26Filter_et, &b_trig_HLT_Photon33_hltEGL1SingleEG26Filter_et);
   fChain->SetBranchAddress("trig_HLT_Photon33_hltEG33L1EG26EtFilter_eta", &trig_HLT_Photon33_hltEG33L1EG26EtFilter_eta, &b_trig_HLT_Photon33_hltEG33L1EG26EtFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon33_hltEG33L1EG26EtFilter_phi", &trig_HLT_Photon33_hltEG33L1EG26EtFilter_phi, &b_trig_HLT_Photon33_hltEG33L1EG26EtFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon33_hltEG33L1EG26EtFilter_et", &trig_HLT_Photon33_hltEG33L1EG26EtFilter_et, &b_trig_HLT_Photon33_hltEG33L1EG26EtFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon33_hltEG33L1EG26HEFilter_eta", &trig_HLT_Photon33_hltEG33L1EG26HEFilter_eta, &b_trig_HLT_Photon33_hltEG33L1EG26HEFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon33_hltEG33L1EG26HEFilter_phi", &trig_HLT_Photon33_hltEG33L1EG26HEFilter_phi, &b_trig_HLT_Photon33_hltEG33L1EG26HEFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon33_hltEG33L1EG26HEFilter_et", &trig_HLT_Photon33_hltEG33L1EG26HEFilter_et, &b_trig_HLT_Photon33_hltEG33L1EG26HEFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon50_accept", &trig_HLT_Photon50_accept, &b_trig_HLT_Photon50_accept);
   fChain->SetBranchAddress("trig_HLT_Photon50_hltL1sSingleEG34to45_eta", &trig_HLT_Photon50_hltL1sSingleEG34to45_eta, &b_trig_HLT_Photon50_hltL1sSingleEG34to45_eta);
   fChain->SetBranchAddress("trig_HLT_Photon50_hltL1sSingleEG34to45_phi", &trig_HLT_Photon50_hltL1sSingleEG34to45_phi, &b_trig_HLT_Photon50_hltL1sSingleEG34to45_phi);
   fChain->SetBranchAddress("trig_HLT_Photon50_hltL1sSingleEG34to45_et", &trig_HLT_Photon50_hltL1sSingleEG34to45_et, &b_trig_HLT_Photon50_hltL1sSingleEG34to45_et);
   fChain->SetBranchAddress("trig_HLT_Photon50_hltEGL1SingleEG34to45Filter_eta", &trig_HLT_Photon50_hltEGL1SingleEG34to45Filter_eta, &b_trig_HLT_Photon50_hltEGL1SingleEG34to45Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon50_hltEGL1SingleEG34to45Filter_phi", &trig_HLT_Photon50_hltEGL1SingleEG34to45Filter_phi, &b_trig_HLT_Photon50_hltEGL1SingleEG34to45Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon50_hltEGL1SingleEG34to45Filter_et", &trig_HLT_Photon50_hltEGL1SingleEG34to45Filter_et, &b_trig_HLT_Photon50_hltEGL1SingleEG34to45Filter_et);
   fChain->SetBranchAddress("trig_HLT_Photon50_hltEG50EtFilter_eta", &trig_HLT_Photon50_hltEG50EtFilter_eta, &b_trig_HLT_Photon50_hltEG50EtFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon50_hltEG50EtFilter_phi", &trig_HLT_Photon50_hltEG50EtFilter_phi, &b_trig_HLT_Photon50_hltEG50EtFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon50_hltEG50EtFilter_et", &trig_HLT_Photon50_hltEG50EtFilter_et, &b_trig_HLT_Photon50_hltEG50EtFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon50_hltEG50HEFilter_eta", &trig_HLT_Photon50_hltEG50HEFilter_eta, &b_trig_HLT_Photon50_hltEG50HEFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon50_hltEG50HEFilter_phi", &trig_HLT_Photon50_hltEG50HEFilter_phi, &b_trig_HLT_Photon50_hltEG50HEFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon50_hltEG50HEFilter_et", &trig_HLT_Photon50_hltEG50HEFilter_et, &b_trig_HLT_Photon50_hltEG50HEFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon75_accept", &trig_HLT_Photon75_accept, &b_trig_HLT_Photon75_accept);
   fChain->SetBranchAddress("trig_HLT_Photon75_hltL1sSingleEG34to50_eta", &trig_HLT_Photon75_hltL1sSingleEG34to50_eta, &b_trig_HLT_Photon75_hltL1sSingleEG34to50_eta);
   fChain->SetBranchAddress("trig_HLT_Photon75_hltL1sSingleEG34to50_phi", &trig_HLT_Photon75_hltL1sSingleEG34to50_phi, &b_trig_HLT_Photon75_hltL1sSingleEG34to50_phi);
   fChain->SetBranchAddress("trig_HLT_Photon75_hltL1sSingleEG34to50_et", &trig_HLT_Photon75_hltL1sSingleEG34to50_et, &b_trig_HLT_Photon75_hltL1sSingleEG34to50_et);
   fChain->SetBranchAddress("trig_HLT_Photon75_hltEGL1SingleEG34ORL1SingleEG40Filter_eta", &trig_HLT_Photon75_hltEGL1SingleEG34ORL1SingleEG40Filter_eta, &b_trig_HLT_Photon75_hltEGL1SingleEG34ORL1SingleEG40Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon75_hltEGL1SingleEG34ORL1SingleEG40Filter_phi", &trig_HLT_Photon75_hltEGL1SingleEG34ORL1SingleEG40Filter_phi, &b_trig_HLT_Photon75_hltEGL1SingleEG34ORL1SingleEG40Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon75_hltEGL1SingleEG34ORL1SingleEG40Filter_et", &trig_HLT_Photon75_hltEGL1SingleEG34ORL1SingleEG40Filter_et, &b_trig_HLT_Photon75_hltEGL1SingleEG34ORL1SingleEG40Filter_et);
   fChain->SetBranchAddress("trig_HLT_Photon75_hltEG75EtFilter_eta", &trig_HLT_Photon75_hltEG75EtFilter_eta, &b_trig_HLT_Photon75_hltEG75EtFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon75_hltEG75EtFilter_phi", &trig_HLT_Photon75_hltEG75EtFilter_phi, &b_trig_HLT_Photon75_hltEG75EtFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon75_hltEG75EtFilter_et", &trig_HLT_Photon75_hltEG75EtFilter_et, &b_trig_HLT_Photon75_hltEG75EtFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon75_hltEG75HEFilter_eta", &trig_HLT_Photon75_hltEG75HEFilter_eta, &b_trig_HLT_Photon75_hltEG75HEFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon75_hltEG75HEFilter_phi", &trig_HLT_Photon75_hltEG75HEFilter_phi, &b_trig_HLT_Photon75_hltEG75HEFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon75_hltEG75HEFilter_et", &trig_HLT_Photon75_hltEG75HEFilter_et, &b_trig_HLT_Photon75_hltEG75HEFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon90_accept", &trig_HLT_Photon90_accept, &b_trig_HLT_Photon90_accept);
   fChain->SetBranchAddress("trig_HLT_Photon90_hltL1sSingleEG34to50_eta", &trig_HLT_Photon90_hltL1sSingleEG34to50_eta, &b_trig_HLT_Photon90_hltL1sSingleEG34to50_eta);
   fChain->SetBranchAddress("trig_HLT_Photon90_hltL1sSingleEG34to50_phi", &trig_HLT_Photon90_hltL1sSingleEG34to50_phi, &b_trig_HLT_Photon90_hltL1sSingleEG34to50_phi);
   fChain->SetBranchAddress("trig_HLT_Photon90_hltL1sSingleEG34to50_et", &trig_HLT_Photon90_hltL1sSingleEG34to50_et, &b_trig_HLT_Photon90_hltL1sSingleEG34to50_et);
   fChain->SetBranchAddress("trig_HLT_Photon90_hltEGL1SingleEG34ORL1SingleEG40Filter_eta", &trig_HLT_Photon90_hltEGL1SingleEG34ORL1SingleEG40Filter_eta, &b_trig_HLT_Photon90_hltEGL1SingleEG34ORL1SingleEG40Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon90_hltEGL1SingleEG34ORL1SingleEG40Filter_phi", &trig_HLT_Photon90_hltEGL1SingleEG34ORL1SingleEG40Filter_phi, &b_trig_HLT_Photon90_hltEGL1SingleEG34ORL1SingleEG40Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon90_hltEGL1SingleEG34ORL1SingleEG40Filter_et", &trig_HLT_Photon90_hltEGL1SingleEG34ORL1SingleEG40Filter_et, &b_trig_HLT_Photon90_hltEGL1SingleEG34ORL1SingleEG40Filter_et);
   fChain->SetBranchAddress("trig_HLT_Photon90_hltEG90EtFilter_eta", &trig_HLT_Photon90_hltEG90EtFilter_eta, &b_trig_HLT_Photon90_hltEG90EtFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon90_hltEG90EtFilter_phi", &trig_HLT_Photon90_hltEG90EtFilter_phi, &b_trig_HLT_Photon90_hltEG90EtFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon90_hltEG90EtFilter_et", &trig_HLT_Photon90_hltEG90EtFilter_et, &b_trig_HLT_Photon90_hltEG90EtFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon90_hltEG90HEFilter_eta", &trig_HLT_Photon90_hltEG90HEFilter_eta, &b_trig_HLT_Photon90_hltEG90HEFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon90_hltEG90HEFilter_phi", &trig_HLT_Photon90_hltEG90HEFilter_phi, &b_trig_HLT_Photon90_hltEG90HEFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon90_hltEG90HEFilter_et", &trig_HLT_Photon90_hltEG90HEFilter_et, &b_trig_HLT_Photon90_hltEG90HEFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon120_accept", &trig_HLT_Photon120_accept, &b_trig_HLT_Photon120_accept);
   fChain->SetBranchAddress("trig_HLT_Photon120_hltL1sSingleEG34to50_eta", &trig_HLT_Photon120_hltL1sSingleEG34to50_eta, &b_trig_HLT_Photon120_hltL1sSingleEG34to50_eta);
   fChain->SetBranchAddress("trig_HLT_Photon120_hltL1sSingleEG34to50_phi", &trig_HLT_Photon120_hltL1sSingleEG34to50_phi, &b_trig_HLT_Photon120_hltL1sSingleEG34to50_phi);
   fChain->SetBranchAddress("trig_HLT_Photon120_hltL1sSingleEG34to50_et", &trig_HLT_Photon120_hltL1sSingleEG34to50_et, &b_trig_HLT_Photon120_hltL1sSingleEG34to50_et);
   fChain->SetBranchAddress("trig_HLT_Photon120_hltEGL1SingleEG34ORL1SingleEG40Filter_eta", &trig_HLT_Photon120_hltEGL1SingleEG34ORL1SingleEG40Filter_eta, &b_trig_HLT_Photon120_hltEGL1SingleEG34ORL1SingleEG40Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon120_hltEGL1SingleEG34ORL1SingleEG40Filter_phi", &trig_HLT_Photon120_hltEGL1SingleEG34ORL1SingleEG40Filter_phi, &b_trig_HLT_Photon120_hltEGL1SingleEG34ORL1SingleEG40Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon120_hltEGL1SingleEG34ORL1SingleEG40Filter_et", &trig_HLT_Photon120_hltEGL1SingleEG34ORL1SingleEG40Filter_et, &b_trig_HLT_Photon120_hltEGL1SingleEG34ORL1SingleEG40Filter_et);
   fChain->SetBranchAddress("trig_HLT_Photon120_hltEG120EtFilter_eta", &trig_HLT_Photon120_hltEG120EtFilter_eta, &b_trig_HLT_Photon120_hltEG120EtFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon120_hltEG120EtFilter_phi", &trig_HLT_Photon120_hltEG120EtFilter_phi, &b_trig_HLT_Photon120_hltEG120EtFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon120_hltEG120EtFilter_et", &trig_HLT_Photon120_hltEG120EtFilter_et, &b_trig_HLT_Photon120_hltEG120EtFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon120_hltEG120HEFilter_eta", &trig_HLT_Photon120_hltEG120HEFilter_eta, &b_trig_HLT_Photon120_hltEG120HEFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon120_hltEG120HEFilter_phi", &trig_HLT_Photon120_hltEG120HEFilter_phi, &b_trig_HLT_Photon120_hltEG120HEFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon120_hltEG120HEFilter_et", &trig_HLT_Photon120_hltEG120HEFilter_et, &b_trig_HLT_Photon120_hltEG120HEFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon150_accept", &trig_HLT_Photon150_accept, &b_trig_HLT_Photon150_accept);
   fChain->SetBranchAddress("trig_HLT_Photon150_hltL1sSingleEGNonIsoOrWithJetAndTau_eta", &trig_HLT_Photon150_hltL1sSingleEGNonIsoOrWithJetAndTau_eta, &b_trig_HLT_Photon150_hltL1sSingleEGNonIsoOrWithJetAndTau_eta);
   fChain->SetBranchAddress("trig_HLT_Photon150_hltL1sSingleEGNonIsoOrWithJetAndTau_phi", &trig_HLT_Photon150_hltL1sSingleEGNonIsoOrWithJetAndTau_phi, &b_trig_HLT_Photon150_hltL1sSingleEGNonIsoOrWithJetAndTau_phi);
   fChain->SetBranchAddress("trig_HLT_Photon150_hltL1sSingleEGNonIsoOrWithJetAndTau_et", &trig_HLT_Photon150_hltL1sSingleEGNonIsoOrWithJetAndTau_et, &b_trig_HLT_Photon150_hltL1sSingleEGNonIsoOrWithJetAndTau_et);
   fChain->SetBranchAddress("trig_HLT_Photon150_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_eta", &trig_HLT_Photon150_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_eta, &b_trig_HLT_Photon150_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon150_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_phi", &trig_HLT_Photon150_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_phi, &b_trig_HLT_Photon150_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon150_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_et", &trig_HLT_Photon150_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_et, &b_trig_HLT_Photon150_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon150_hltEG150EtFilter_eta", &trig_HLT_Photon150_hltEG150EtFilter_eta, &b_trig_HLT_Photon150_hltEG150EtFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon150_hltEG150EtFilter_phi", &trig_HLT_Photon150_hltEG150EtFilter_phi, &b_trig_HLT_Photon150_hltEG150EtFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon150_hltEG150EtFilter_et", &trig_HLT_Photon150_hltEG150EtFilter_et, &b_trig_HLT_Photon150_hltEG150EtFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon150_hltEG150HEFilter_eta", &trig_HLT_Photon150_hltEG150HEFilter_eta, &b_trig_HLT_Photon150_hltEG150HEFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon150_hltEG150HEFilter_phi", &trig_HLT_Photon150_hltEG150HEFilter_phi, &b_trig_HLT_Photon150_hltEG150HEFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon150_hltEG150HEFilter_et", &trig_HLT_Photon150_hltEG150HEFilter_et, &b_trig_HLT_Photon150_hltEG150HEFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon175_accept", &trig_HLT_Photon175_accept, &b_trig_HLT_Photon175_accept);
   fChain->SetBranchAddress("trig_HLT_Photon175_hltL1sSingleEGNonIsoOrWithJetAndTau_eta", &trig_HLT_Photon175_hltL1sSingleEGNonIsoOrWithJetAndTau_eta, &b_trig_HLT_Photon175_hltL1sSingleEGNonIsoOrWithJetAndTau_eta);
   fChain->SetBranchAddress("trig_HLT_Photon175_hltL1sSingleEGNonIsoOrWithJetAndTau_phi", &trig_HLT_Photon175_hltL1sSingleEGNonIsoOrWithJetAndTau_phi, &b_trig_HLT_Photon175_hltL1sSingleEGNonIsoOrWithJetAndTau_phi);
   fChain->SetBranchAddress("trig_HLT_Photon175_hltL1sSingleEGNonIsoOrWithJetAndTau_et", &trig_HLT_Photon175_hltL1sSingleEGNonIsoOrWithJetAndTau_et, &b_trig_HLT_Photon175_hltL1sSingleEGNonIsoOrWithJetAndTau_et);
   fChain->SetBranchAddress("trig_HLT_Photon175_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_eta", &trig_HLT_Photon175_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_eta, &b_trig_HLT_Photon175_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon175_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_phi", &trig_HLT_Photon175_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_phi, &b_trig_HLT_Photon175_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon175_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_et", &trig_HLT_Photon175_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_et, &b_trig_HLT_Photon175_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon175_hltEG175EtFilter_eta", &trig_HLT_Photon175_hltEG175EtFilter_eta, &b_trig_HLT_Photon175_hltEG175EtFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon175_hltEG175EtFilter_phi", &trig_HLT_Photon175_hltEG175EtFilter_phi, &b_trig_HLT_Photon175_hltEG175EtFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon175_hltEG175EtFilter_et", &trig_HLT_Photon175_hltEG175EtFilter_et, &b_trig_HLT_Photon175_hltEG175EtFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon175_hltEG175HEFilter_eta", &trig_HLT_Photon175_hltEG175HEFilter_eta, &b_trig_HLT_Photon175_hltEG175HEFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon175_hltEG175HEFilter_phi", &trig_HLT_Photon175_hltEG175HEFilter_phi, &b_trig_HLT_Photon175_hltEG175HEFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon175_hltEG175HEFilter_et", &trig_HLT_Photon175_hltEG175HEFilter_et, &b_trig_HLT_Photon175_hltEG175HEFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon200_accept", &trig_HLT_Photon200_accept, &b_trig_HLT_Photon200_accept);
   fChain->SetBranchAddress("trig_HLT_Photon200_hltL1sSingleEGNonIsoOrWithJetAndTau_eta", &trig_HLT_Photon200_hltL1sSingleEGNonIsoOrWithJetAndTau_eta, &b_trig_HLT_Photon200_hltL1sSingleEGNonIsoOrWithJetAndTau_eta);
   fChain->SetBranchAddress("trig_HLT_Photon200_hltL1sSingleEGNonIsoOrWithJetAndTau_phi", &trig_HLT_Photon200_hltL1sSingleEGNonIsoOrWithJetAndTau_phi, &b_trig_HLT_Photon200_hltL1sSingleEGNonIsoOrWithJetAndTau_phi);
   fChain->SetBranchAddress("trig_HLT_Photon200_hltL1sSingleEGNonIsoOrWithJetAndTau_et", &trig_HLT_Photon200_hltL1sSingleEGNonIsoOrWithJetAndTau_et, &b_trig_HLT_Photon200_hltL1sSingleEGNonIsoOrWithJetAndTau_et);
   fChain->SetBranchAddress("trig_HLT_Photon200_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_eta", &trig_HLT_Photon200_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_eta, &b_trig_HLT_Photon200_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon200_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_phi", &trig_HLT_Photon200_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_phi, &b_trig_HLT_Photon200_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon200_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_et", &trig_HLT_Photon200_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_et, &b_trig_HLT_Photon200_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon200_hltEG200EtFilter_eta", &trig_HLT_Photon200_hltEG200EtFilter_eta, &b_trig_HLT_Photon200_hltEG200EtFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon200_hltEG200EtFilter_phi", &trig_HLT_Photon200_hltEG200EtFilter_phi, &b_trig_HLT_Photon200_hltEG200EtFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon200_hltEG200EtFilter_et", &trig_HLT_Photon200_hltEG200EtFilter_et, &b_trig_HLT_Photon200_hltEG200EtFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon200_hltEG200HEFilter_eta", &trig_HLT_Photon200_hltEG200HEFilter_eta, &b_trig_HLT_Photon200_hltEG200HEFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon200_hltEG200HEFilter_phi", &trig_HLT_Photon200_hltEG200HEFilter_phi, &b_trig_HLT_Photon200_hltEG200HEFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon200_hltEG200HEFilter_et", &trig_HLT_Photon200_hltEG200HEFilter_et, &b_trig_HLT_Photon200_hltEG200HEFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon50_R9Id90_HE10_IsoM_accept", &trig_HLT_Photon50_R9Id90_HE10_IsoM_accept, &b_trig_HLT_Photon50_R9Id90_HE10_IsoM_accept);
   fChain->SetBranchAddress("trig_HLT_Photon50_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_eta", &trig_HLT_Photon50_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_eta, &b_trig_HLT_Photon50_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_eta);
   fChain->SetBranchAddress("trig_HLT_Photon50_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_phi", &trig_HLT_Photon50_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_phi, &b_trig_HLT_Photon50_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_phi);
   fChain->SetBranchAddress("trig_HLT_Photon50_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_et", &trig_HLT_Photon50_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_et, &b_trig_HLT_Photon50_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_et);
   fChain->SetBranchAddress("trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_eta", &trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_eta, &b_trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_phi", &trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_phi, &b_trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_et", &trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_et, &b_trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_et);
   fChain->SetBranchAddress("trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMEtFilter_eta", &trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMEtFilter_eta, &b_trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMEtFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMEtFilter_phi", &trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMEtFilter_phi, &b_trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMEtFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMEtFilter_et", &trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMEtFilter_et, &b_trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMEtFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMHEFilter_eta", &trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMHEFilter_eta, &b_trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMHEFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMHEFilter_phi", &trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMHEFilter_phi, &b_trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMHEFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMHEFilter_et", &trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMHEFilter_et, &b_trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMHEFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMR9Filter_eta", &trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMR9Filter_eta, &b_trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMR9Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMR9Filter_phi", &trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMR9Filter_phi, &b_trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMR9Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMR9Filter_et", &trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMR9Filter_et, &b_trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMR9Filter_et);
   fChain->SetBranchAddress("trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMEcalIsoFilter_eta", &trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMEcalIsoFilter_eta, &b_trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMEcalIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMEcalIsoFilter_phi", &trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMEcalIsoFilter_phi, &b_trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMEcalIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMEcalIsoFilter_et", &trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMEcalIsoFilter_et, &b_trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMEcalIsoFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMHcalIsoFilter_eta", &trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMHcalIsoFilter_eta, &b_trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMHcalIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMHcalIsoFilter_phi", &trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMHcalIsoFilter_phi, &b_trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMHcalIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMHcalIsoFilter_et", &trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMHcalIsoFilter_et, &b_trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMHcalIsoFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMTrackIsoFilter_eta", &trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMTrackIsoFilter_eta, &b_trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMTrackIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMTrackIsoFilter_phi", &trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMTrackIsoFilter_phi, &b_trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMTrackIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMTrackIsoFilter_et", &trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMTrackIsoFilter_et, &b_trig_HLT_Photon50_R9Id90_HE10_IsoM_hltEG50R9Id90HE10IsoMTrackIsoFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon75_R9Id90_HE10_IsoM_accept", &trig_HLT_Photon75_R9Id90_HE10_IsoM_accept, &b_trig_HLT_Photon75_R9Id90_HE10_IsoM_accept);
   fChain->SetBranchAddress("trig_HLT_Photon75_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_eta", &trig_HLT_Photon75_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_eta, &b_trig_HLT_Photon75_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_eta);
   fChain->SetBranchAddress("trig_HLT_Photon75_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_phi", &trig_HLT_Photon75_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_phi, &b_trig_HLT_Photon75_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_phi);
   fChain->SetBranchAddress("trig_HLT_Photon75_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_et", &trig_HLT_Photon75_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_et, &b_trig_HLT_Photon75_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_et);
   fChain->SetBranchAddress("trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_eta", &trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_eta, &b_trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_phi", &trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_phi, &b_trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_et", &trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_et, &b_trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_et);
   fChain->SetBranchAddress("trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMEtFilter_eta", &trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMEtFilter_eta, &b_trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMEtFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMEtFilter_phi", &trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMEtFilter_phi, &b_trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMEtFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMEtFilter_et", &trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMEtFilter_et, &b_trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMEtFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMHEFilter_eta", &trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMHEFilter_eta, &b_trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMHEFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMHEFilter_phi", &trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMHEFilter_phi, &b_trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMHEFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMHEFilter_et", &trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMHEFilter_et, &b_trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMHEFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMR9Filter_eta", &trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMR9Filter_eta, &b_trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMR9Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMR9Filter_phi", &trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMR9Filter_phi, &b_trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMR9Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMR9Filter_et", &trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMR9Filter_et, &b_trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMR9Filter_et);
   fChain->SetBranchAddress("trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMEcalIsoFilter_eta", &trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMEcalIsoFilter_eta, &b_trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMEcalIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMEcalIsoFilter_phi", &trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMEcalIsoFilter_phi, &b_trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMEcalIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMEcalIsoFilter_et", &trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMEcalIsoFilter_et, &b_trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMEcalIsoFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMHcalIsoFilter_eta", &trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMHcalIsoFilter_eta, &b_trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMHcalIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMHcalIsoFilter_phi", &trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMHcalIsoFilter_phi, &b_trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMHcalIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMHcalIsoFilter_et", &trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMHcalIsoFilter_et, &b_trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMHcalIsoFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMTrackIsoFilter_eta", &trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMTrackIsoFilter_eta, &b_trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMTrackIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMTrackIsoFilter_phi", &trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMTrackIsoFilter_phi, &b_trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMTrackIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMTrackIsoFilter_et", &trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMTrackIsoFilter_et, &b_trig_HLT_Photon75_R9Id90_HE10_IsoM_hltEG75R9Id90HE10IsoMTrackIsoFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon90_R9Id90_HE10_IsoM_accept", &trig_HLT_Photon90_R9Id90_HE10_IsoM_accept, &b_trig_HLT_Photon90_R9Id90_HE10_IsoM_accept);
   fChain->SetBranchAddress("trig_HLT_Photon90_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_eta", &trig_HLT_Photon90_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_eta, &b_trig_HLT_Photon90_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_eta);
   fChain->SetBranchAddress("trig_HLT_Photon90_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_phi", &trig_HLT_Photon90_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_phi, &b_trig_HLT_Photon90_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_phi);
   fChain->SetBranchAddress("trig_HLT_Photon90_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_et", &trig_HLT_Photon90_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_et, &b_trig_HLT_Photon90_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_et);
   fChain->SetBranchAddress("trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_eta", &trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_eta, &b_trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_phi", &trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_phi, &b_trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_et", &trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_et, &b_trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_et);
   fChain->SetBranchAddress("trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMEtFilter_eta", &trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMEtFilter_eta, &b_trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMEtFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMEtFilter_phi", &trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMEtFilter_phi, &b_trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMEtFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMEtFilter_et", &trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMEtFilter_et, &b_trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMEtFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMHEFilter_eta", &trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMHEFilter_eta, &b_trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMHEFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMHEFilter_phi", &trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMHEFilter_phi, &b_trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMHEFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMHEFilter_et", &trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMHEFilter_et, &b_trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMHEFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMR9Filter_eta", &trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMR9Filter_eta, &b_trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMR9Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMR9Filter_phi", &trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMR9Filter_phi, &b_trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMR9Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMR9Filter_et", &trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMR9Filter_et, &b_trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMR9Filter_et);
   fChain->SetBranchAddress("trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMEcalIsoFilter_eta", &trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMEcalIsoFilter_eta, &b_trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMEcalIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMEcalIsoFilter_phi", &trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMEcalIsoFilter_phi, &b_trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMEcalIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMEcalIsoFilter_et", &trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMEcalIsoFilter_et, &b_trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMEcalIsoFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMHcalIsoFilter_eta", &trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMHcalIsoFilter_eta, &b_trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMHcalIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMHcalIsoFilter_phi", &trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMHcalIsoFilter_phi, &b_trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMHcalIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMHcalIsoFilter_et", &trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMHcalIsoFilter_et, &b_trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMHcalIsoFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMTrackIsoFilter_eta", &trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMTrackIsoFilter_eta, &b_trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMTrackIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMTrackIsoFilter_phi", &trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMTrackIsoFilter_phi, &b_trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMTrackIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMTrackIsoFilter_et", &trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMTrackIsoFilter_et, &b_trig_HLT_Photon90_R9Id90_HE10_IsoM_hltEG90R9Id90HE10IsoMTrackIsoFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon120_R9Id90_HE10_IsoM_accept", &trig_HLT_Photon120_R9Id90_HE10_IsoM_accept, &b_trig_HLT_Photon120_R9Id90_HE10_IsoM_accept);
   fChain->SetBranchAddress("trig_HLT_Photon120_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_eta", &trig_HLT_Photon120_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_eta, &b_trig_HLT_Photon120_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_eta);
   fChain->SetBranchAddress("trig_HLT_Photon120_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_phi", &trig_HLT_Photon120_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_phi, &b_trig_HLT_Photon120_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_phi);
   fChain->SetBranchAddress("trig_HLT_Photon120_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_et", &trig_HLT_Photon120_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_et, &b_trig_HLT_Photon120_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_et);
   fChain->SetBranchAddress("trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_eta", &trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_eta, &b_trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_phi", &trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_phi, &b_trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_et", &trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_et, &b_trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_et);
   fChain->SetBranchAddress("trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMEtFilter_eta", &trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMEtFilter_eta, &b_trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMEtFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMEtFilter_phi", &trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMEtFilter_phi, &b_trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMEtFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMEtFilter_et", &trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMEtFilter_et, &b_trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMEtFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMHEFilter_eta", &trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMHEFilter_eta, &b_trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMHEFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMHEFilter_phi", &trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMHEFilter_phi, &b_trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMHEFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMHEFilter_et", &trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMHEFilter_et, &b_trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMHEFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMR9Filter_eta", &trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMR9Filter_eta, &b_trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMR9Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMR9Filter_phi", &trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMR9Filter_phi, &b_trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMR9Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMR9Filter_et", &trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMR9Filter_et, &b_trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMR9Filter_et);
   fChain->SetBranchAddress("trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMEcalIsoFilter_eta", &trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMEcalIsoFilter_eta, &b_trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMEcalIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMEcalIsoFilter_phi", &trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMEcalIsoFilter_phi, &b_trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMEcalIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMEcalIsoFilter_et", &trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMEcalIsoFilter_et, &b_trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMEcalIsoFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMHcalIsoFilter_eta", &trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMHcalIsoFilter_eta, &b_trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMHcalIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMHcalIsoFilter_phi", &trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMHcalIsoFilter_phi, &b_trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMHcalIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMHcalIsoFilter_et", &trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMHcalIsoFilter_et, &b_trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMHcalIsoFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMTrackIsoFilter_eta", &trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMTrackIsoFilter_eta, &b_trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMTrackIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMTrackIsoFilter_phi", &trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMTrackIsoFilter_phi, &b_trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMTrackIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMTrackIsoFilter_et", &trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMTrackIsoFilter_et, &b_trig_HLT_Photon120_R9Id90_HE10_IsoM_hltEG120R9Id90HE10IsoMTrackIsoFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon165_R9Id90_HE10_IsoM_accept", &trig_HLT_Photon165_R9Id90_HE10_IsoM_accept, &b_trig_HLT_Photon165_R9Id90_HE10_IsoM_accept);
   fChain->SetBranchAddress("trig_HLT_Photon165_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_eta", &trig_HLT_Photon165_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_eta, &b_trig_HLT_Photon165_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_eta);
   fChain->SetBranchAddress("trig_HLT_Photon165_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_phi", &trig_HLT_Photon165_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_phi, &b_trig_HLT_Photon165_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_phi);
   fChain->SetBranchAddress("trig_HLT_Photon165_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_et", &trig_HLT_Photon165_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_et, &b_trig_HLT_Photon165_R9Id90_HE10_IsoM_hltL1sSingleEG40to50_et);
   fChain->SetBranchAddress("trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_eta", &trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_eta, &b_trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_phi", &trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_phi, &b_trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_et", &trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_et, &b_trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEGL1SingleEG40Filter_et);
   fChain->SetBranchAddress("trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMEtFilter_eta", &trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMEtFilter_eta, &b_trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMEtFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMEtFilter_phi", &trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMEtFilter_phi, &b_trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMEtFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMEtFilter_et", &trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMEtFilter_et, &b_trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMEtFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMHEFilter_eta", &trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMHEFilter_eta, &b_trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMHEFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMHEFilter_phi", &trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMHEFilter_phi, &b_trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMHEFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMHEFilter_et", &trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMHEFilter_et, &b_trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMHEFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMR9Filter_eta", &trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMR9Filter_eta, &b_trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMR9Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMR9Filter_phi", &trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMR9Filter_phi, &b_trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMR9Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMR9Filter_et", &trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMR9Filter_et, &b_trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMR9Filter_et);
   fChain->SetBranchAddress("trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMEcalIsoFilter_eta", &trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMEcalIsoFilter_eta, &b_trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMEcalIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMEcalIsoFilter_phi", &trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMEcalIsoFilter_phi, &b_trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMEcalIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMEcalIsoFilter_et", &trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMEcalIsoFilter_et, &b_trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMEcalIsoFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMHcalIsoFilter_eta", &trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMHcalIsoFilter_eta, &b_trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMHcalIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMHcalIsoFilter_phi", &trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMHcalIsoFilter_phi, &b_trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMHcalIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMHcalIsoFilter_et", &trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMHcalIsoFilter_et, &b_trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMHcalIsoFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMTrackIsoFilter_eta", &trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMTrackIsoFilter_eta, &b_trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMTrackIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMTrackIsoFilter_phi", &trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMTrackIsoFilter_phi, &b_trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMTrackIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMTrackIsoFilter_et", &trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMTrackIsoFilter_et, &b_trig_HLT_Photon165_R9Id90_HE10_IsoM_hltEG165R9Id90HE10IsoMTrackIsoFilter_et);
   fChain->SetBranchAddress("trig_HLT_DoubleMu43NoFiltersNoVtx_accept", &trig_HLT_DoubleMu43NoFiltersNoVtx_accept, &b_trig_HLT_DoubleMu43NoFiltersNoVtx_accept);
   fChain->SetBranchAddress("trig_HLT_DoubleMu48NoFiltersNoVtx_accept", &trig_HLT_DoubleMu48NoFiltersNoVtx_accept, &b_trig_HLT_DoubleMu48NoFiltersNoVtx_accept);
   fChain->SetBranchAddress("trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_accept", &trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_accept, &b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_accept);
   fChain->SetBranchAddress("trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL1sMu5EG23IorMu7EG23IorMu20EG17IorMu23EG10_eta", &trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL1sMu5EG23IorMu7EG23IorMu20EG17IorMu23EG10_eta, &b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL1sMu5EG23IorMu7EG23IorMu20EG17IorMu23EG10_eta);
   fChain->SetBranchAddress("trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL1sMu5EG23IorMu7EG23IorMu20EG17IorMu23EG10_phi", &trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL1sMu5EG23IorMu7EG23IorMu20EG17IorMu23EG10_phi, &b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL1sMu5EG23IorMu7EG23IorMu20EG17IorMu23EG10_phi);
   fChain->SetBranchAddress("trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL1sMu5EG23IorMu7EG23IorMu20EG17IorMu23EG10_et", &trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL1sMu5EG23IorMu7EG23IorMu20EG17IorMu23EG10_et, &b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL1sMu5EG23IorMu7EG23IorMu20EG17IorMu23EG10_et);
   fChain->SetBranchAddress("trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL1fL1sMu5EG20orMu20EG15L1Filtered5_eta", &trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL1fL1sMu5EG20orMu20EG15L1Filtered5_eta, &b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL1fL1sMu5EG20orMu20EG15L1Filtered5_eta);
   fChain->SetBranchAddress("trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL1fL1sMu5EG20orMu20EG15L1Filtered5_phi", &trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL1fL1sMu5EG20orMu20EG15L1Filtered5_phi, &b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL1fL1sMu5EG20orMu20EG15L1Filtered5_phi);
   fChain->SetBranchAddress("trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL1fL1sMu5EG20orMu20EG15L1Filtered5_et", &trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL1fL1sMu5EG20orMu20EG15L1Filtered5_et, &b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL1fL1sMu5EG20orMu20EG15L1Filtered5_et);
   fChain->SetBranchAddress("trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_-hltL2fL1sMu5EG20orMu20EG15L1f5L2NoVtxFiltered16_eta", &trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_-hltL2fL1sMu5EG20orMu20EG15L1f5L2NoVtxFiltered16_eta, &b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_-hltL2fL1sMu5EG20orMu20EG15L1f5L2NoVtxFiltered16_eta);
   fChain->SetBranchAddress("trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_-hltL2fL1sMu5EG20orMu20EG15L1f5L2NoVtxFiltered16_phi", &trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_-hltL2fL1sMu5EG20orMu20EG15L1f5L2NoVtxFiltered16_phi, &b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_-hltL2fL1sMu5EG20orMu20EG15L1f5L2NoVtxFiltered16_phi);
   fChain->SetBranchAddress("trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_-hltL2fL1sMu5EG20orMu20EG15L1f5L2NoVtxFiltered16_et", &trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_-hltL2fL1sMu5EG20orMu20EG15L1f5L2NoVtxFiltered16_et, &b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_-hltL2fL1sMu5EG20orMu20EG15L1f5L2NoVtxFiltered16_et);
   fChain->SetBranchAddress("trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL3fL1sMu5EG20orMu20EG15L1f5L2NVf16L3NoFiltersNoVtxFiltered43_eta", &trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL3fL1sMu5EG20orMu20EG15L1f5L2NVf16L3NoFiltersNoVtxFiltered43_eta, &b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL3fL1sMu5EG20orMu20EG15L1f5L2NVf16L3NoFiltersNoVtxFiltered43_eta);
   fChain->SetBranchAddress("trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL3fL1sMu5EG20orMu20EG15L1f5L2NVf16L3NoFiltersNoVtxFiltered43_phi", &trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL3fL1sMu5EG20orMu20EG15L1f5L2NVf16L3NoFiltersNoVtxFiltered43_phi, &b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL3fL1sMu5EG20orMu20EG15L1f5L2NVf16L3NoFiltersNoVtxFiltered43_phi);
   fChain->SetBranchAddress("trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL3fL1sMu5EG20orMu20EG15L1f5L2NVf16L3NoFiltersNoVtxFiltered43_et", &trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL3fL1sMu5EG20orMu20EG15L1f5L2NVf16L3NoFiltersNoVtxFiltered43_et, &b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltL3fL1sMu5EG20orMu20EG15L1f5L2NVf16L3NoFiltersNoVtxFiltered43_et);
   fChain->SetBranchAddress("trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltEGL1Mu5EG20Filter_eta", &trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltEGL1Mu5EG20Filter_eta, &b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltEGL1Mu5EG20Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltEGL1Mu5EG20Filter_phi", &trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltEGL1Mu5EG20Filter_phi, &b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltEGL1Mu5EG20Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltEGL1Mu5EG20Filter_et", &trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltEGL1Mu5EG20Filter_et, &b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltEGL1Mu5EG20Filter_et);
   fChain->SetBranchAddress("trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLEtFilter_eta", &trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLEtFilter_eta, &b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLEtFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLEtFilter_phi", &trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLEtFilter_phi, &b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLEtFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLEtFilter_et", &trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLEtFilter_et, &b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLEtFilter_et);
   fChain->SetBranchAddress("trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLClusterShapeFilter_eta", &trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLClusterShapeFilter_eta, &b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLClusterShapeFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLClusterShapeFilter_phi", &trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLClusterShapeFilter_phi, &b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLClusterShapeFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLClusterShapeFilter_et", &trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLClusterShapeFilter_et, &b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLClusterShapeFilter_et);
   fChain->SetBranchAddress("trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLHEFilter_eta", &trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLHEFilter_eta, &b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLHEFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLHEFilter_phi", &trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLHEFilter_phi, &b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLHEFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLHEFilter_et", &trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLHEFilter_et, &b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_hltMu43NoFiltersNoVtxPhoton43CaloIdLHEFilter_et);
   fChain->SetBranchAddress("trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_accept", &trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_accept, &b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_accept);
   fChain->SetBranchAddress("trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL1sMu5EG23IorMu7EG23IorMu20EG17IorMu23EG10_eta", &trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL1sMu5EG23IorMu7EG23IorMu20EG17IorMu23EG10_eta, &b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL1sMu5EG23IorMu7EG23IorMu20EG17IorMu23EG10_eta);
   fChain->SetBranchAddress("trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL1sMu5EG23IorMu7EG23IorMu20EG17IorMu23EG10_phi", &trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL1sMu5EG23IorMu7EG23IorMu20EG17IorMu23EG10_phi, &b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL1sMu5EG23IorMu7EG23IorMu20EG17IorMu23EG10_phi);
   fChain->SetBranchAddress("trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL1sMu5EG23IorMu7EG23IorMu20EG17IorMu23EG10_et", &trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL1sMu5EG23IorMu7EG23IorMu20EG17IorMu23EG10_et, &b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL1sMu5EG23IorMu7EG23IorMu20EG17IorMu23EG10_et);
   fChain->SetBranchAddress("trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL1fL1sMu5EG20orMu20EG15L1Filtered5_eta", &trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL1fL1sMu5EG20orMu20EG15L1Filtered5_eta, &b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL1fL1sMu5EG20orMu20EG15L1Filtered5_eta);
   fChain->SetBranchAddress("trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL1fL1sMu5EG20orMu20EG15L1Filtered5_phi", &trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL1fL1sMu5EG20orMu20EG15L1Filtered5_phi, &b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL1fL1sMu5EG20orMu20EG15L1Filtered5_phi);
   fChain->SetBranchAddress("trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL1fL1sMu5EG20orMu20EG15L1Filtered5_et", &trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL1fL1sMu5EG20orMu20EG15L1Filtered5_et, &b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL1fL1sMu5EG20orMu20EG15L1Filtered5_et);
   fChain->SetBranchAddress("trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_-hltL2fL1sMu5EG20orMu20EG15L1f5L2NoVtxFiltered16_eta", &trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_-hltL2fL1sMu5EG20orMu20EG15L1f5L2NoVtxFiltered16_eta, &b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_-hltL2fL1sMu5EG20orMu20EG15L1f5L2NoVtxFiltered16_eta);
   fChain->SetBranchAddress("trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_-hltL2fL1sMu5EG20orMu20EG15L1f5L2NoVtxFiltered16_phi", &trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_-hltL2fL1sMu5EG20orMu20EG15L1f5L2NoVtxFiltered16_phi, &b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_-hltL2fL1sMu5EG20orMu20EG15L1f5L2NoVtxFiltered16_phi);
   fChain->SetBranchAddress("trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_-hltL2fL1sMu5EG20orMu20EG15L1f5L2NoVtxFiltered16_et", &trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_-hltL2fL1sMu5EG20orMu20EG15L1f5L2NoVtxFiltered16_et, &b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_-hltL2fL1sMu5EG20orMu20EG15L1f5L2NoVtxFiltered16_et);
   fChain->SetBranchAddress("trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL3fL1sMu5EG20orMu20EG15L1f5L2NVf16L3NoFiltersNoVtxFiltered48_eta", &trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL3fL1sMu5EG20orMu20EG15L1f5L2NVf16L3NoFiltersNoVtxFiltered48_eta, &b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL3fL1sMu5EG20orMu20EG15L1f5L2NVf16L3NoFiltersNoVtxFiltered48_eta);
   fChain->SetBranchAddress("trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL3fL1sMu5EG20orMu20EG15L1f5L2NVf16L3NoFiltersNoVtxFiltered48_phi", &trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL3fL1sMu5EG20orMu20EG15L1f5L2NVf16L3NoFiltersNoVtxFiltered48_phi, &b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL3fL1sMu5EG20orMu20EG15L1f5L2NVf16L3NoFiltersNoVtxFiltered48_phi);
   fChain->SetBranchAddress("trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL3fL1sMu5EG20orMu20EG15L1f5L2NVf16L3NoFiltersNoVtxFiltered48_et", &trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL3fL1sMu5EG20orMu20EG15L1f5L2NVf16L3NoFiltersNoVtxFiltered48_et, &b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltL3fL1sMu5EG20orMu20EG15L1f5L2NVf16L3NoFiltersNoVtxFiltered48_et);
   fChain->SetBranchAddress("trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltEGL1Mu5EG20Filter_eta", &trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltEGL1Mu5EG20Filter_eta, &b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltEGL1Mu5EG20Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltEGL1Mu5EG20Filter_phi", &trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltEGL1Mu5EG20Filter_phi, &b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltEGL1Mu5EG20Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltEGL1Mu5EG20Filter_et", &trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltEGL1Mu5EG20Filter_et, &b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltEGL1Mu5EG20Filter_et);
   fChain->SetBranchAddress("trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLEtFilter_eta", &trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLEtFilter_eta, &b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLEtFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLEtFilter_phi", &trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLEtFilter_phi, &b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLEtFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLEtFilter_et", &trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLEtFilter_et, &b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLEtFilter_et);
   fChain->SetBranchAddress("trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLClusterShapeFilter_eta", &trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLClusterShapeFilter_eta, &b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLClusterShapeFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLClusterShapeFilter_phi", &trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLClusterShapeFilter_phi, &b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLClusterShapeFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLClusterShapeFilter_et", &trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLClusterShapeFilter_et, &b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLClusterShapeFilter_et);
   fChain->SetBranchAddress("trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLHEFilter_eta", &trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLHEFilter_eta, &b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLHEFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLHEFilter_phi", &trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLHEFilter_phi, &b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLHEFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLHEFilter_et", &trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLHEFilter_et, &b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_hltMu48NoFiltersNoVtxPhoton48CaloIdLHEFilter_et);
   fChain->SetBranchAddress("trig_HLT_DoubleMu20_7_Mass0to30_L1_DM4_accept", &trig_HLT_DoubleMu20_7_Mass0to30_L1_DM4_accept, &b_trig_HLT_DoubleMu20_7_Mass0to30_L1_DM4_accept);
   fChain->SetBranchAddress("trig_HLT_DoubleMu20_7_Mass0to30_L1_DM4EG_accept", &trig_HLT_DoubleMu20_7_Mass0to30_L1_DM4EG_accept, &b_trig_HLT_DoubleMu20_7_Mass0to30_L1_DM4EG_accept);
   fChain->SetBranchAddress("trig_HLT_DoubleMu20_7_Mass0to30_Photon23_accept", &trig_HLT_DoubleMu20_7_Mass0to30_Photon23_accept, &b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_accept);
   fChain->SetBranchAddress("trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltL1sDoubleMu4OSEG12ORDoubleMu5OSEG12_eta", &trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltL1sDoubleMu4OSEG12ORDoubleMu5OSEG12_eta, &b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltL1sDoubleMu4OSEG12ORDoubleMu5OSEG12_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltL1sDoubleMu4OSEG12ORDoubleMu5OSEG12_phi", &trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltL1sDoubleMu4OSEG12ORDoubleMu5OSEG12_phi, &b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltL1sDoubleMu4OSEG12ORDoubleMu5OSEG12_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltL1sDoubleMu4OSEG12ORDoubleMu5OSEG12_et", &trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltL1sDoubleMu4OSEG12ORDoubleMu5OSEG12_et, &b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltL1sDoubleMu4OSEG12ORDoubleMu5OSEG12_et);
   fChain->SetBranchAddress("trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltL1fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1Filtered0_eta", &trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltL1fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1Filtered0_eta, &b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltL1fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1Filtered0_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltL1fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1Filtered0_phi", &trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltL1fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1Filtered0_phi, &b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltL1fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1Filtered0_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltL1fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1Filtered0_et", &trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltL1fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1Filtered0_et, &b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltL1fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1Filtered0_et);
   fChain->SetBranchAddress("trig_HLT_DoubleMu20_7_Mass0to30_Photon23_-hltL2fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1f0L2PreFiltered0_eta", &trig_HLT_DoubleMu20_7_Mass0to30_Photon23_-hltL2fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1f0L2PreFiltered0_eta, &b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_-hltL2fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1f0L2PreFiltered0_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleMu20_7_Mass0to30_Photon23_-hltL2fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1f0L2PreFiltered0_phi", &trig_HLT_DoubleMu20_7_Mass0to30_Photon23_-hltL2fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1f0L2PreFiltered0_phi, &b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_-hltL2fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1f0L2PreFiltered0_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleMu20_7_Mass0to30_Photon23_-hltL2fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1f0L2PreFiltered0_et", &trig_HLT_DoubleMu20_7_Mass0to30_Photon23_-hltL2fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1f0L2PreFiltered0_et, &b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_-hltL2fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1f0L2PreFiltered0_et);
   fChain->SetBranchAddress("trig_HLT_DoubleMu20_7_Mass0to30_Photon23_-hltL1fForIterL3L1fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1Filtered0_eta", &trig_HLT_DoubleMu20_7_Mass0to30_Photon23_-hltL1fForIterL3L1fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1Filtered0_eta, &b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_-hltL1fForIterL3L1fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1Filtered0_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleMu20_7_Mass0to30_Photon23_-hltL1fForIterL3L1fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1Filtered0_phi", &trig_HLT_DoubleMu20_7_Mass0to30_Photon23_-hltL1fForIterL3L1fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1Filtered0_phi, &b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_-hltL1fForIterL3L1fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1Filtered0_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleMu20_7_Mass0to30_Photon23_-hltL1fForIterL3L1fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1Filtered0_et", &trig_HLT_DoubleMu20_7_Mass0to30_Photon23_-hltL1fForIterL3L1fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1Filtered0_et, &b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_-hltL1fForIterL3L1fL1sDoubleMu4OSEG12ORDoubleMu5OSEG12L1Filtered0_et);
   fChain->SetBranchAddress("trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltDoubleMu207Mass0to30Photon23L3Filtered_eta", &trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltDoubleMu207Mass0to30Photon23L3Filtered_eta, &b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltDoubleMu207Mass0to30Photon23L3Filtered_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltDoubleMu207Mass0to30Photon23L3Filtered_phi", &trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltDoubleMu207Mass0to30Photon23L3Filtered_phi, &b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltDoubleMu207Mass0to30Photon23L3Filtered_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltDoubleMu207Mass0to30Photon23L3Filtered_et", &trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltDoubleMu207Mass0to30Photon23L3Filtered_et, &b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltDoubleMu207Mass0to30Photon23L3Filtered_et);
   fChain->SetBranchAddress("trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEGDoubleMu4OSEG12ORDoubleMu5OSEG12Filter_eta", &trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEGDoubleMu4OSEG12ORDoubleMu5OSEG12Filter_eta, &b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEGDoubleMu4OSEG12ORDoubleMu5OSEG12Filter_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEGDoubleMu4OSEG12ORDoubleMu5OSEG12Filter_phi", &trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEGDoubleMu4OSEG12ORDoubleMu5OSEG12Filter_phi, &b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEGDoubleMu4OSEG12ORDoubleMu5OSEG12Filter_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEGDoubleMu4OSEG12ORDoubleMu5OSEG12Filter_et", &trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEGDoubleMu4OSEG12ORDoubleMu5OSEG12Filter_et, &b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEGDoubleMu4OSEG12ORDoubleMu5OSEG12Filter_et);
   fChain->SetBranchAddress("trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEG23EtFilter_eta", &trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEG23EtFilter_eta, &b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEG23EtFilter_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEG23EtFilter_phi", &trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEG23EtFilter_phi, &b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEG23EtFilter_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEG23EtFilter_et", &trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEG23EtFilter_et, &b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEG23EtFilter_et);
   fChain->SetBranchAddress("trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEG23HEFilter_eta", &trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEG23HEFilter_eta, &b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEG23HEFilter_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEG23HEFilter_phi", &trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEG23HEFilter_phi, &b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEG23HEFilter_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEG23HEFilter_et", &trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEG23HEFilter_et, &b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_hltEG23HEFilter_et);
   fChain->SetBranchAddress("trig_HLT_Dimuon14_Phi_Barrel_Seagulls_accept", &trig_HLT_Dimuon14_Phi_Barrel_Seagulls_accept, &b_trig_HLT_Dimuon14_Phi_Barrel_Seagulls_accept);
   fChain->SetBranchAddress("trig_HLT_Dimuon24_Phi_noCorrL1_accept", &trig_HLT_Dimuon24_Phi_noCorrL1_accept, &b_trig_HLT_Dimuon24_Phi_noCorrL1_accept);
   fChain->SetBranchAddress("trig_HLT_DoubleIsoMu20_eta2p1_accept", &trig_HLT_DoubleIsoMu20_eta2p1_accept, &b_trig_HLT_DoubleIsoMu20_eta2p1_accept);
   fChain->SetBranchAddress("trig_HLT_DoubleIsoMu24_eta2p1_accept", &trig_HLT_DoubleIsoMu24_eta2p1_accept, &b_trig_HLT_DoubleIsoMu24_eta2p1_accept);
   fChain->SetBranchAddress("trig_HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx_accept", &trig_HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx_accept, &b_trig_HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx_accept);
   fChain->SetBranchAddress("trig_HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx_accept", &trig_HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx_accept, &b_trig_HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx_accept);
   fChain->SetBranchAddress("trig_HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx_accept", &trig_HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx_accept, &b_trig_HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx_accept);
   fChain->SetBranchAddress("trig_HLT_Mu8_accept", &trig_HLT_Mu8_accept, &b_trig_HLT_Mu8_accept);
   fChain->SetBranchAddress("trig_HLT_Mu17_accept", &trig_HLT_Mu17_accept, &b_trig_HLT_Mu17_accept);
   fChain->SetBranchAddress("trig_HLT_Mu19_accept", &trig_HLT_Mu19_accept, &b_trig_HLT_Mu19_accept);
   fChain->SetBranchAddress("trig_HLT_Mu17_Photon30_IsoCaloId_accept", &trig_HLT_Mu17_Photon30_IsoCaloId_accept, &b_trig_HLT_Mu17_Photon30_IsoCaloId_accept);
   fChain->SetBranchAddress("trig_HLT_Mu17_Photon30_IsoCaloId_hltL1sMu5EG23IorMu5IsoEG20IorMu7EG23IorMu7IsoEG20IorMuIso7EG23_eta", &trig_HLT_Mu17_Photon30_IsoCaloId_hltL1sMu5EG23IorMu5IsoEG20IorMu7EG23IorMu7IsoEG20IorMuIso7EG23_eta, &b_trig_HLT_Mu17_Photon30_IsoCaloId_hltL1sMu5EG23IorMu5IsoEG20IorMu7EG23IorMu7IsoEG20IorMuIso7EG23_eta);
   fChain->SetBranchAddress("trig_HLT_Mu17_Photon30_IsoCaloId_hltL1sMu5EG23IorMu5IsoEG20IorMu7EG23IorMu7IsoEG20IorMuIso7EG23_phi", &trig_HLT_Mu17_Photon30_IsoCaloId_hltL1sMu5EG23IorMu5IsoEG20IorMu7EG23IorMu7IsoEG20IorMuIso7EG23_phi, &b_trig_HLT_Mu17_Photon30_IsoCaloId_hltL1sMu5EG23IorMu5IsoEG20IorMu7EG23IorMu7IsoEG20IorMuIso7EG23_phi);
   fChain->SetBranchAddress("trig_HLT_Mu17_Photon30_IsoCaloId_hltL1sMu5EG23IorMu5IsoEG20IorMu7EG23IorMu7IsoEG20IorMuIso7EG23_et", &trig_HLT_Mu17_Photon30_IsoCaloId_hltL1sMu5EG23IorMu5IsoEG20IorMu7EG23IorMu7IsoEG20IorMuIso7EG23_et, &b_trig_HLT_Mu17_Photon30_IsoCaloId_hltL1sMu5EG23IorMu5IsoEG20IorMu7EG23IorMu7IsoEG20IorMuIso7EG23_et);
   fChain->SetBranchAddress("trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdMuonlegL1Filtered7_eta", &trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdMuonlegL1Filtered7_eta, &b_trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdMuonlegL1Filtered7_eta);
   fChain->SetBranchAddress("trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdMuonlegL1Filtered7_phi", &trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdMuonlegL1Filtered7_phi, &b_trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdMuonlegL1Filtered7_phi);
   fChain->SetBranchAddress("trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdMuonlegL1Filtered7_et", &trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdMuonlegL1Filtered7_et, &b_trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdMuonlegL1Filtered7_et);
   fChain->SetBranchAddress("trig_HLT_Mu17_Photon30_IsoCaloId_-hltMu17Photon30IsoCaloIdMuonlegL2Filtered10_eta", &trig_HLT_Mu17_Photon30_IsoCaloId_-hltMu17Photon30IsoCaloIdMuonlegL2Filtered10_eta, &b_trig_HLT_Mu17_Photon30_IsoCaloId_-hltMu17Photon30IsoCaloIdMuonlegL2Filtered10_eta);
   fChain->SetBranchAddress("trig_HLT_Mu17_Photon30_IsoCaloId_-hltMu17Photon30IsoCaloIdMuonlegL2Filtered10_phi", &trig_HLT_Mu17_Photon30_IsoCaloId_-hltMu17Photon30IsoCaloIdMuonlegL2Filtered10_phi, &b_trig_HLT_Mu17_Photon30_IsoCaloId_-hltMu17Photon30IsoCaloIdMuonlegL2Filtered10_phi);
   fChain->SetBranchAddress("trig_HLT_Mu17_Photon30_IsoCaloId_-hltMu17Photon30IsoCaloIdMuonlegL2Filtered10_et", &trig_HLT_Mu17_Photon30_IsoCaloId_-hltMu17Photon30IsoCaloIdMuonlegL2Filtered10_et, &b_trig_HLT_Mu17_Photon30_IsoCaloId_-hltMu17Photon30IsoCaloIdMuonlegL2Filtered10_et);
   fChain->SetBranchAddress("trig_HLT_Mu17_Photon30_IsoCaloId_-hltL1fForIterL3Mu17Photon30IsoCaloIdMuonlegL1Filtered7_eta", &trig_HLT_Mu17_Photon30_IsoCaloId_-hltL1fForIterL3Mu17Photon30IsoCaloIdMuonlegL1Filtered7_eta, &b_trig_HLT_Mu17_Photon30_IsoCaloId_-hltL1fForIterL3Mu17Photon30IsoCaloIdMuonlegL1Filtered7_eta);
   fChain->SetBranchAddress("trig_HLT_Mu17_Photon30_IsoCaloId_-hltL1fForIterL3Mu17Photon30IsoCaloIdMuonlegL1Filtered7_phi", &trig_HLT_Mu17_Photon30_IsoCaloId_-hltL1fForIterL3Mu17Photon30IsoCaloIdMuonlegL1Filtered7_phi, &b_trig_HLT_Mu17_Photon30_IsoCaloId_-hltL1fForIterL3Mu17Photon30IsoCaloIdMuonlegL1Filtered7_phi);
   fChain->SetBranchAddress("trig_HLT_Mu17_Photon30_IsoCaloId_-hltL1fForIterL3Mu17Photon30IsoCaloIdMuonlegL1Filtered7_et", &trig_HLT_Mu17_Photon30_IsoCaloId_-hltL1fForIterL3Mu17Photon30IsoCaloIdMuonlegL1Filtered7_et, &b_trig_HLT_Mu17_Photon30_IsoCaloId_-hltL1fForIterL3Mu17Photon30IsoCaloIdMuonlegL1Filtered7_et);
   fChain->SetBranchAddress("trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdMuonlegL3Filtered17Q_eta", &trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdMuonlegL3Filtered17Q_eta, &b_trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdMuonlegL3Filtered17Q_eta);
   fChain->SetBranchAddress("trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdMuonlegL3Filtered17Q_phi", &trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdMuonlegL3Filtered17Q_phi, &b_trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdMuonlegL3Filtered17Q_phi);
   fChain->SetBranchAddress("trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdMuonlegL3Filtered17Q_et", &trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdMuonlegL3Filtered17Q_et, &b_trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdMuonlegL3Filtered17Q_et);
   fChain->SetBranchAddress("trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegFilter_eta", &trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegFilter_eta, &b_trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegFilter_phi", &trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegFilter_phi, &b_trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegFilter_et", &trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegFilter_et, &b_trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegFilter_et);
   fChain->SetBranchAddress("trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegEtFilter_eta", &trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegEtFilter_eta, &b_trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegEtFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegEtFilter_phi", &trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegEtFilter_phi, &b_trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegEtFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegEtFilter_et", &trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegEtFilter_et, &b_trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegEtFilter_et);
   fChain->SetBranchAddress("trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegClusterShapeFilter_eta", &trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegClusterShapeFilter_eta, &b_trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegClusterShapeFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegClusterShapeFilter_phi", &trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegClusterShapeFilter_phi, &b_trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegClusterShapeFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegClusterShapeFilter_et", &trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegClusterShapeFilter_et, &b_trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegClusterShapeFilter_et);
   fChain->SetBranchAddress("trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegHEFilter_eta", &trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegHEFilter_eta, &b_trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegHEFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegHEFilter_phi", &trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegHEFilter_phi, &b_trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegHEFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegHEFilter_et", &trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegHEFilter_et, &b_trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegHEFilter_et);
   fChain->SetBranchAddress("trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegEcalIsoFilter_eta", &trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegEcalIsoFilter_eta, &b_trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegEcalIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegEcalIsoFilter_phi", &trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegEcalIsoFilter_phi, &b_trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegEcalIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegEcalIsoFilter_et", &trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegEcalIsoFilter_et, &b_trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegEcalIsoFilter_et);
   fChain->SetBranchAddress("trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegTrackIsoFilter_eta", &trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegTrackIsoFilter_eta, &b_trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegTrackIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegTrackIsoFilter_phi", &trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegTrackIsoFilter_phi, &b_trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegTrackIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegTrackIsoFilter_et", &trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegTrackIsoFilter_et, &b_trig_HLT_Mu17_Photon30_IsoCaloId_hltMu17Photon30IsoCaloIdPhotonlegTrackIsoFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_accept", &trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_accept, &b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_accept);
   fChain->SetBranchAddress("trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltL1sSingleEGNonIsoOrWithJetAndTau_eta", &trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltL1sSingleEGNonIsoOrWithJetAndTau_eta, &b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltL1sSingleEGNonIsoOrWithJetAndTau_eta);
   fChain->SetBranchAddress("trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltL1sSingleEGNonIsoOrWithJetAndTau_phi", &trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltL1sSingleEGNonIsoOrWithJetAndTau_phi, &b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltL1sSingleEGNonIsoOrWithJetAndTau_phi);
   fChain->SetBranchAddress("trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltL1sSingleEGNonIsoOrWithJetAndTau_et", &trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltL1sSingleEGNonIsoOrWithJetAndTau_et, &b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltL1sSingleEGNonIsoOrWithJetAndTau_et);
   fChain->SetBranchAddress("trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_eta", &trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_eta, &b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_phi", &trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_phi, &b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_et", &trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_et, &b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEGL1SingleEGNonIsoOrWithJetAndTauFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115EtFilter_eta", &trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115EtFilter_eta, &b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115EtFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115EtFilter_phi", &trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115EtFilter_phi, &b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115EtFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115EtFilter_et", &trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115EtFilter_et, &b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115EtFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115CaloIdVTClusterShapeFilter_eta", &trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115CaloIdVTClusterShapeFilter_eta, &b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115CaloIdVTClusterShapeFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115CaloIdVTClusterShapeFilter_phi", &trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115CaloIdVTClusterShapeFilter_phi, &b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115CaloIdVTClusterShapeFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115CaloIdVTClusterShapeFilter_et", &trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115CaloIdVTClusterShapeFilter_et, &b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115CaloIdVTClusterShapeFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115CaloIdVTHEFilter_eta", &trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115CaloIdVTHEFilter_eta, &b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115CaloIdVTHEFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115CaloIdVTHEFilter_phi", &trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115CaloIdVTHEFilter_phi, &b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115CaloIdVTHEFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115CaloIdVTHEFilter_et", &trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115CaloIdVTHEFilter_et, &b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEG115CaloIdVTHEFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTPixelMatchFilter_eta", &trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTPixelMatchFilter_eta, &b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTPixelMatchFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTPixelMatchFilter_phi", &trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTPixelMatchFilter_phi, &b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTPixelMatchFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTPixelMatchFilter_et", &trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTPixelMatchFilter_et, &b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTPixelMatchFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTGsfTrkIdTGsfDetaFilter_eta", &trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTGsfTrkIdTGsfDetaFilter_eta, &b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTGsfTrkIdTGsfDetaFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTGsfTrkIdTGsfDetaFilter_phi", &trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTGsfTrkIdTGsfDetaFilter_phi, &b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTGsfTrkIdTGsfDetaFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTGsfTrkIdTGsfDetaFilter_et", &trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTGsfTrkIdTGsfDetaFilter_et, &b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTGsfTrkIdTGsfDetaFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTGsfTrkIdTGsfDphiFilter_eta", &trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTGsfTrkIdTGsfDphiFilter_eta, &b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTGsfTrkIdTGsfDphiFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTGsfTrkIdTGsfDphiFilter_phi", &trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTGsfTrkIdTGsfDphiFilter_phi, &b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTGsfTrkIdTGsfDphiFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTGsfTrkIdTGsfDphiFilter_et", &trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTGsfTrkIdTGsfDphiFilter_et, &b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_hltEle115CaloIdVTGsfTrkIdTGsfDphiFilter_et);
   fChain->SetBranchAddress("trig_HLT_Ele135_CaloIdVT_GsfTrkIdT_accept", &trig_HLT_Ele135_CaloIdVT_GsfTrkIdT_accept, &b_trig_HLT_Ele135_CaloIdVT_GsfTrkIdT_accept);
   fChain->SetBranchAddress("trig_HLT_Ele145_CaloIdVT_GsfTrkIdT_accept", &trig_HLT_Ele145_CaloIdVT_GsfTrkIdT_accept, &b_trig_HLT_Ele145_CaloIdVT_GsfTrkIdT_accept);
   fChain->SetBranchAddress("trig_HLT_Ele200_CaloIdVT_GsfTrkIdT_accept", &trig_HLT_Ele200_CaloIdVT_GsfTrkIdT_accept, &b_trig_HLT_Ele200_CaloIdVT_GsfTrkIdT_accept);
   fChain->SetBranchAddress("trig_HLT_Ele250_CaloIdVT_GsfTrkIdT_accept", &trig_HLT_Ele250_CaloIdVT_GsfTrkIdT_accept, &b_trig_HLT_Ele250_CaloIdVT_GsfTrkIdT_accept);
   fChain->SetBranchAddress("trig_HLT_Ele300_CaloIdVT_GsfTrkIdT_accept", &trig_HLT_Ele300_CaloIdVT_GsfTrkIdT_accept, &b_trig_HLT_Ele300_CaloIdVT_GsfTrkIdT_accept);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_accept", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_accept, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_accept);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltL1sSingleEGNonIsoOrWithJetAndTauNoPS_eta", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltL1sSingleEGNonIsoOrWithJetAndTauNoPS_eta, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltL1sSingleEGNonIsoOrWithJetAndTauNoPS_eta);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltL1sSingleEGNonIsoOrWithJetAndTauNoPS_phi", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltL1sSingleEGNonIsoOrWithJetAndTauNoPS_phi, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltL1sSingleEGNonIsoOrWithJetAndTauNoPS_phi);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltL1sSingleEGNonIsoOrWithJetAndTauNoPS_et", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltL1sSingleEGNonIsoOrWithJetAndTauNoPS_et, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltL1sSingleEGNonIsoOrWithJetAndTauNoPS_et);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEGL1SingleEGNonIsoOrWithJetAndTauNoPSFilter_eta", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEGL1SingleEGNonIsoOrWithJetAndTauNoPSFilter_eta, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEGL1SingleEGNonIsoOrWithJetAndTauNoPSFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEGL1SingleEGNonIsoOrWithJetAndTauNoPSFilter_phi", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEGL1SingleEGNonIsoOrWithJetAndTauNoPSFilter_phi, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEGL1SingleEGNonIsoOrWithJetAndTauNoPSFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEGL1SingleEGNonIsoOrWithJetAndTauNoPSFilter_et", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEGL1SingleEGNonIsoOrWithJetAndTauNoPSFilter_et, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEGL1SingleEGNonIsoOrWithJetAndTauNoPSFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60EtFilter_eta", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60EtFilter_eta, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60EtFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60EtFilter_phi", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60EtFilter_phi, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60EtFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60EtFilter_et", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60EtFilter_et, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60EtFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60HEFilter_eta", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60HEFilter_eta, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60HEFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60HEFilter_phi", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60HEFilter_phi, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60HEFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60HEFilter_et", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60HEFilter_et, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60HEFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLR9IdFilter_eta", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLR9IdFilter_eta, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLR9IdFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLR9IdFilter_phi", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLR9IdFilter_phi, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLR9IdFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLR9IdFilter_et", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLR9IdFilter_et, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLR9IdFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLClusterShapeFilter_eta", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLClusterShapeFilter_eta, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLClusterShapeFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLClusterShapeFilter_phi", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLClusterShapeFilter_phi, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLClusterShapeFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLClusterShapeFilter_et", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLClusterShapeFilter_et, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLClusterShapeFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLEcalPFClusterIsoFilter_eta", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLEcalPFClusterIsoFilter_eta, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLEcalPFClusterIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLEcalPFClusterIsoFilter_phi", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLEcalPFClusterIsoFilter_phi, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLEcalPFClusterIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLEcalPFClusterIsoFilter_et", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLEcalPFClusterIsoFilter_et, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLEcalPFClusterIsoFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLHcalPFClusterIsoFilter_eta", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLHcalPFClusterIsoFilter_eta, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLHcalPFClusterIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLHcalPFClusterIsoFilter_phi", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLHcalPFClusterIsoFilter_phi, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLHcalPFClusterIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLHcalPFClusterIsoFilter_et", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLHcalPFClusterIsoFilter_et, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLHcalPFClusterIsoFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLHollowTrackIsoFilter_eta", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLHollowTrackIsoFilter_eta, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLHollowTrackIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLHollowTrackIsoFilter_phi", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLHollowTrackIsoFilter_phi, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLHollowTrackIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLHollowTrackIsoFilter_et", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLHollowTrackIsoFilter_et, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_hltEG60R9Id90CaloIdLIsoLHollowTrackIsoFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_accept", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_accept, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_accept);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltL1sSingleEGNonIsoOrWithJetAndTauNoPS_eta", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltL1sSingleEGNonIsoOrWithJetAndTauNoPS_eta, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltL1sSingleEGNonIsoOrWithJetAndTauNoPS_eta);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltL1sSingleEGNonIsoOrWithJetAndTauNoPS_phi", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltL1sSingleEGNonIsoOrWithJetAndTauNoPS_phi, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltL1sSingleEGNonIsoOrWithJetAndTauNoPS_phi);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltL1sSingleEGNonIsoOrWithJetAndTauNoPS_et", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltL1sSingleEGNonIsoOrWithJetAndTauNoPS_et, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltL1sSingleEGNonIsoOrWithJetAndTauNoPS_et);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEGL1SingleEGNonIsoOrWithJetAndTauNoPSFilter_eta", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEGL1SingleEGNonIsoOrWithJetAndTauNoPSFilter_eta, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEGL1SingleEGNonIsoOrWithJetAndTauNoPSFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEGL1SingleEGNonIsoOrWithJetAndTauNoPSFilter_phi", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEGL1SingleEGNonIsoOrWithJetAndTauNoPSFilter_phi, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEGL1SingleEGNonIsoOrWithJetAndTauNoPSFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEGL1SingleEGNonIsoOrWithJetAndTauNoPSFilter_et", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEGL1SingleEGNonIsoOrWithJetAndTauNoPSFilter_et, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEGL1SingleEGNonIsoOrWithJetAndTauNoPSFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60EtFilter_eta", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60EtFilter_eta, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60EtFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60EtFilter_phi", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60EtFilter_phi, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60EtFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60EtFilter_et", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60EtFilter_et, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60EtFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60HEFilter_eta", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60HEFilter_eta, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60HEFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60HEFilter_phi", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60HEFilter_phi, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60HEFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60HEFilter_et", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60HEFilter_et, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60HEFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLR9IdFilter_eta", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLR9IdFilter_eta, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLR9IdFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLR9IdFilter_phi", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLR9IdFilter_phi, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLR9IdFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLR9IdFilter_et", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLR9IdFilter_et, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLR9IdFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLClusterShapeFilter_eta", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLClusterShapeFilter_eta, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLClusterShapeFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLClusterShapeFilter_phi", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLClusterShapeFilter_phi, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLClusterShapeFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLClusterShapeFilter_et", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLClusterShapeFilter_et, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLClusterShapeFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLEcalPFClusterIsoFilter_eta", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLEcalPFClusterIsoFilter_eta, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLEcalPFClusterIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLEcalPFClusterIsoFilter_phi", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLEcalPFClusterIsoFilter_phi, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLEcalPFClusterIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLEcalPFClusterIsoFilter_et", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLEcalPFClusterIsoFilter_et, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLEcalPFClusterIsoFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLHcalPFClusterIsoFilter_eta", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLHcalPFClusterIsoFilter_eta, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLHcalPFClusterIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLHcalPFClusterIsoFilter_phi", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLHcalPFClusterIsoFilter_phi, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLHcalPFClusterIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLHcalPFClusterIsoFilter_et", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLHcalPFClusterIsoFilter_et, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLHcalPFClusterIsoFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLHollowTrackIsoFilter_eta", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLHollowTrackIsoFilter_eta, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLHollowTrackIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLHollowTrackIsoFilter_phi", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLHollowTrackIsoFilter_phi, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLHollowTrackIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLHollowTrackIsoFilter_et", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLHollowTrackIsoFilter_et, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLHollowTrackIsoFilter_et);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLDisplacedIdFilter_eta", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLDisplacedIdFilter_eta, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLDisplacedIdFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLDisplacedIdFilter_phi", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLDisplacedIdFilter_phi, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLDisplacedIdFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLDisplacedIdFilter_et", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLDisplacedIdFilter_et, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_hltEG60R9Id90CaloIdLIsoLDisplacedIdFilter_et);
   fChain->SetBranchAddress("trig_HLT_SinglePhoton10_Eta3p1ForPPRef_accept", &trig_HLT_SinglePhoton10_Eta3p1ForPPRef_accept, &b_trig_HLT_SinglePhoton10_Eta3p1ForPPRef_accept);
   fChain->SetBranchAddress("trig_HLT_SinglePhoton10_Eta3p1ForPPRef_hltL1sZeroBias_eta", &trig_HLT_SinglePhoton10_Eta3p1ForPPRef_hltL1sZeroBias_eta, &b_trig_HLT_SinglePhoton10_Eta3p1ForPPRef_hltL1sZeroBias_eta);
   fChain->SetBranchAddress("trig_HLT_SinglePhoton10_Eta3p1ForPPRef_hltL1sZeroBias_phi", &trig_HLT_SinglePhoton10_Eta3p1ForPPRef_hltL1sZeroBias_phi, &b_trig_HLT_SinglePhoton10_Eta3p1ForPPRef_hltL1sZeroBias_phi);
   fChain->SetBranchAddress("trig_HLT_SinglePhoton10_Eta3p1ForPPRef_hltL1sZeroBias_et", &trig_HLT_SinglePhoton10_Eta3p1ForPPRef_hltL1sZeroBias_et, &b_trig_HLT_SinglePhoton10_Eta3p1ForPPRef_hltL1sZeroBias_et);
   fChain->SetBranchAddress("trig_HLT_SinglePhoton10_Eta3p1ForPPRef_hltHIPhoton10Eta3p1_eta", &trig_HLT_SinglePhoton10_Eta3p1ForPPRef_hltHIPhoton10Eta3p1_eta, &b_trig_HLT_SinglePhoton10_Eta3p1ForPPRef_hltHIPhoton10Eta3p1_eta);
   fChain->SetBranchAddress("trig_HLT_SinglePhoton10_Eta3p1ForPPRef_hltHIPhoton10Eta3p1_phi", &trig_HLT_SinglePhoton10_Eta3p1ForPPRef_hltHIPhoton10Eta3p1_phi, &b_trig_HLT_SinglePhoton10_Eta3p1ForPPRef_hltHIPhoton10Eta3p1_phi);
   fChain->SetBranchAddress("trig_HLT_SinglePhoton10_Eta3p1ForPPRef_hltHIPhoton10Eta3p1_et", &trig_HLT_SinglePhoton10_Eta3p1ForPPRef_hltHIPhoton10Eta3p1_et, &b_trig_HLT_SinglePhoton10_Eta3p1ForPPRef_hltHIPhoton10Eta3p1_et);
   fChain->SetBranchAddress("trig_HLT_SinglePhoton20_Eta3p1ForPPRef_accept", &trig_HLT_SinglePhoton20_Eta3p1ForPPRef_accept, &b_trig_HLT_SinglePhoton20_Eta3p1ForPPRef_accept);
   fChain->SetBranchAddress("trig_HLT_SinglePhoton20_Eta3p1ForPPRef_hltL1sSingleEG10_eta", &trig_HLT_SinglePhoton20_Eta3p1ForPPRef_hltL1sSingleEG10_eta, &b_trig_HLT_SinglePhoton20_Eta3p1ForPPRef_hltL1sSingleEG10_eta);
   fChain->SetBranchAddress("trig_HLT_SinglePhoton20_Eta3p1ForPPRef_hltL1sSingleEG10_phi", &trig_HLT_SinglePhoton20_Eta3p1ForPPRef_hltL1sSingleEG10_phi, &b_trig_HLT_SinglePhoton20_Eta3p1ForPPRef_hltL1sSingleEG10_phi);
   fChain->SetBranchAddress("trig_HLT_SinglePhoton20_Eta3p1ForPPRef_hltL1sSingleEG10_et", &trig_HLT_SinglePhoton20_Eta3p1ForPPRef_hltL1sSingleEG10_et, &b_trig_HLT_SinglePhoton20_Eta3p1ForPPRef_hltL1sSingleEG10_et);
   fChain->SetBranchAddress("trig_HLT_SinglePhoton20_Eta3p1ForPPRef_hltHIPhoton20Eta3p1_eta", &trig_HLT_SinglePhoton20_Eta3p1ForPPRef_hltHIPhoton20Eta3p1_eta, &b_trig_HLT_SinglePhoton20_Eta3p1ForPPRef_hltHIPhoton20Eta3p1_eta);
   fChain->SetBranchAddress("trig_HLT_SinglePhoton20_Eta3p1ForPPRef_hltHIPhoton20Eta3p1_phi", &trig_HLT_SinglePhoton20_Eta3p1ForPPRef_hltHIPhoton20Eta3p1_phi, &b_trig_HLT_SinglePhoton20_Eta3p1ForPPRef_hltHIPhoton20Eta3p1_phi);
   fChain->SetBranchAddress("trig_HLT_SinglePhoton20_Eta3p1ForPPRef_hltHIPhoton20Eta3p1_et", &trig_HLT_SinglePhoton20_Eta3p1ForPPRef_hltHIPhoton20Eta3p1_et, &b_trig_HLT_SinglePhoton20_Eta3p1ForPPRef_hltHIPhoton20Eta3p1_et);
   fChain->SetBranchAddress("trig_HLT_SinglePhoton30_Eta3p1ForPPRef_accept", &trig_HLT_SinglePhoton30_Eta3p1ForPPRef_accept, &b_trig_HLT_SinglePhoton30_Eta3p1ForPPRef_accept);
   fChain->SetBranchAddress("trig_HLT_SinglePhoton30_Eta3p1ForPPRef_hltL1sSingleEG10_eta", &trig_HLT_SinglePhoton30_Eta3p1ForPPRef_hltL1sSingleEG10_eta, &b_trig_HLT_SinglePhoton30_Eta3p1ForPPRef_hltL1sSingleEG10_eta);
   fChain->SetBranchAddress("trig_HLT_SinglePhoton30_Eta3p1ForPPRef_hltL1sSingleEG10_phi", &trig_HLT_SinglePhoton30_Eta3p1ForPPRef_hltL1sSingleEG10_phi, &b_trig_HLT_SinglePhoton30_Eta3p1ForPPRef_hltL1sSingleEG10_phi);
   fChain->SetBranchAddress("trig_HLT_SinglePhoton30_Eta3p1ForPPRef_hltL1sSingleEG10_et", &trig_HLT_SinglePhoton30_Eta3p1ForPPRef_hltL1sSingleEG10_et, &b_trig_HLT_SinglePhoton30_Eta3p1ForPPRef_hltL1sSingleEG10_et);
   fChain->SetBranchAddress("trig_HLT_SinglePhoton30_Eta3p1ForPPRef_hltHIPhoton30Eta3p1_eta", &trig_HLT_SinglePhoton30_Eta3p1ForPPRef_hltHIPhoton30Eta3p1_eta, &b_trig_HLT_SinglePhoton30_Eta3p1ForPPRef_hltHIPhoton30Eta3p1_eta);
   fChain->SetBranchAddress("trig_HLT_SinglePhoton30_Eta3p1ForPPRef_hltHIPhoton30Eta3p1_phi", &trig_HLT_SinglePhoton30_Eta3p1ForPPRef_hltHIPhoton30Eta3p1_phi, &b_trig_HLT_SinglePhoton30_Eta3p1ForPPRef_hltHIPhoton30Eta3p1_phi);
   fChain->SetBranchAddress("trig_HLT_SinglePhoton30_Eta3p1ForPPRef_hltHIPhoton30Eta3p1_et", &trig_HLT_SinglePhoton30_Eta3p1ForPPRef_hltHIPhoton30Eta3p1_et, &b_trig_HLT_SinglePhoton30_Eta3p1ForPPRef_hltHIPhoton30Eta3p1_et);
   fChain->SetBranchAddress("trig_HLT_Photon20_HoverELoose_accept", &trig_HLT_Photon20_HoverELoose_accept, &b_trig_HLT_Photon20_HoverELoose_accept);
   fChain->SetBranchAddress("trig_HLT_Photon20_HoverELoose_hltL1sSingleEG10_eta", &trig_HLT_Photon20_HoverELoose_hltL1sSingleEG10_eta, &b_trig_HLT_Photon20_HoverELoose_hltL1sSingleEG10_eta);
   fChain->SetBranchAddress("trig_HLT_Photon20_HoverELoose_hltL1sSingleEG10_phi", &trig_HLT_Photon20_HoverELoose_hltL1sSingleEG10_phi, &b_trig_HLT_Photon20_HoverELoose_hltL1sSingleEG10_phi);
   fChain->SetBranchAddress("trig_HLT_Photon20_HoverELoose_hltL1sSingleEG10_et", &trig_HLT_Photon20_HoverELoose_hltL1sSingleEG10_et, &b_trig_HLT_Photon20_HoverELoose_hltL1sSingleEG10_et);
   fChain->SetBranchAddress("trig_HLT_Photon20_HoverELoose_hltEGL1SingleEG10Filter_eta", &trig_HLT_Photon20_HoverELoose_hltEGL1SingleEG10Filter_eta, &b_trig_HLT_Photon20_HoverELoose_hltEGL1SingleEG10Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon20_HoverELoose_hltEGL1SingleEG10Filter_phi", &trig_HLT_Photon20_HoverELoose_hltEGL1SingleEG10Filter_phi, &b_trig_HLT_Photon20_HoverELoose_hltEGL1SingleEG10Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon20_HoverELoose_hltEGL1SingleEG10Filter_et", &trig_HLT_Photon20_HoverELoose_hltEGL1SingleEG10Filter_et, &b_trig_HLT_Photon20_HoverELoose_hltEGL1SingleEG10Filter_et);
   fChain->SetBranchAddress("trig_HLT_Photon20_HoverELoose_hltEG20EtFilterLooseHoverE_eta", &trig_HLT_Photon20_HoverELoose_hltEG20EtFilterLooseHoverE_eta, &b_trig_HLT_Photon20_HoverELoose_hltEG20EtFilterLooseHoverE_eta);
   fChain->SetBranchAddress("trig_HLT_Photon20_HoverELoose_hltEG20EtFilterLooseHoverE_phi", &trig_HLT_Photon20_HoverELoose_hltEG20EtFilterLooseHoverE_phi, &b_trig_HLT_Photon20_HoverELoose_hltEG20EtFilterLooseHoverE_phi);
   fChain->SetBranchAddress("trig_HLT_Photon20_HoverELoose_hltEG20EtFilterLooseHoverE_et", &trig_HLT_Photon20_HoverELoose_hltEG20EtFilterLooseHoverE_et, &b_trig_HLT_Photon20_HoverELoose_hltEG20EtFilterLooseHoverE_et);
   fChain->SetBranchAddress("trig_HLT_Photon20_HoverELoose_hltEG20HEFilterLooseHoverE_eta", &trig_HLT_Photon20_HoverELoose_hltEG20HEFilterLooseHoverE_eta, &b_trig_HLT_Photon20_HoverELoose_hltEG20HEFilterLooseHoverE_eta);
   fChain->SetBranchAddress("trig_HLT_Photon20_HoverELoose_hltEG20HEFilterLooseHoverE_phi", &trig_HLT_Photon20_HoverELoose_hltEG20HEFilterLooseHoverE_phi, &b_trig_HLT_Photon20_HoverELoose_hltEG20HEFilterLooseHoverE_phi);
   fChain->SetBranchAddress("trig_HLT_Photon20_HoverELoose_hltEG20HEFilterLooseHoverE_et", &trig_HLT_Photon20_HoverELoose_hltEG20HEFilterLooseHoverE_et, &b_trig_HLT_Photon20_HoverELoose_hltEG20HEFilterLooseHoverE_et);
   fChain->SetBranchAddress("trig_HLT_Photon30_HoverELoose_accept", &trig_HLT_Photon30_HoverELoose_accept, &b_trig_HLT_Photon30_HoverELoose_accept);
   fChain->SetBranchAddress("trig_HLT_Photon30_HoverELoose_hltL1sSingleEG10_eta", &trig_HLT_Photon30_HoverELoose_hltL1sSingleEG10_eta, &b_trig_HLT_Photon30_HoverELoose_hltL1sSingleEG10_eta);
   fChain->SetBranchAddress("trig_HLT_Photon30_HoverELoose_hltL1sSingleEG10_phi", &trig_HLT_Photon30_HoverELoose_hltL1sSingleEG10_phi, &b_trig_HLT_Photon30_HoverELoose_hltL1sSingleEG10_phi);
   fChain->SetBranchAddress("trig_HLT_Photon30_HoverELoose_hltL1sSingleEG10_et", &trig_HLT_Photon30_HoverELoose_hltL1sSingleEG10_et, &b_trig_HLT_Photon30_HoverELoose_hltL1sSingleEG10_et);
   fChain->SetBranchAddress("trig_HLT_Photon30_HoverELoose_hltEGL1SingleEG10Filter_eta", &trig_HLT_Photon30_HoverELoose_hltEGL1SingleEG10Filter_eta, &b_trig_HLT_Photon30_HoverELoose_hltEGL1SingleEG10Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Photon30_HoverELoose_hltEGL1SingleEG10Filter_phi", &trig_HLT_Photon30_HoverELoose_hltEGL1SingleEG10Filter_phi, &b_trig_HLT_Photon30_HoverELoose_hltEGL1SingleEG10Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Photon30_HoverELoose_hltEGL1SingleEG10Filter_et", &trig_HLT_Photon30_HoverELoose_hltEGL1SingleEG10Filter_et, &b_trig_HLT_Photon30_HoverELoose_hltEGL1SingleEG10Filter_et);
   fChain->SetBranchAddress("trig_HLT_Photon30_HoverELoose_hltEG30EtFilterLooseHoverE_eta", &trig_HLT_Photon30_HoverELoose_hltEG30EtFilterLooseHoverE_eta, &b_trig_HLT_Photon30_HoverELoose_hltEG30EtFilterLooseHoverE_eta);
   fChain->SetBranchAddress("trig_HLT_Photon30_HoverELoose_hltEG30EtFilterLooseHoverE_phi", &trig_HLT_Photon30_HoverELoose_hltEG30EtFilterLooseHoverE_phi, &b_trig_HLT_Photon30_HoverELoose_hltEG30EtFilterLooseHoverE_phi);
   fChain->SetBranchAddress("trig_HLT_Photon30_HoverELoose_hltEG30EtFilterLooseHoverE_et", &trig_HLT_Photon30_HoverELoose_hltEG30EtFilterLooseHoverE_et, &b_trig_HLT_Photon30_HoverELoose_hltEG30EtFilterLooseHoverE_et);
   fChain->SetBranchAddress("trig_HLT_Photon30_HoverELoose_hltEG30HEFilterLooseHoverE_eta", &trig_HLT_Photon30_HoverELoose_hltEG30HEFilterLooseHoverE_eta, &b_trig_HLT_Photon30_HoverELoose_hltEG30HEFilterLooseHoverE_eta);
   fChain->SetBranchAddress("trig_HLT_Photon30_HoverELoose_hltEG30HEFilterLooseHoverE_phi", &trig_HLT_Photon30_HoverELoose_hltEG30HEFilterLooseHoverE_phi, &b_trig_HLT_Photon30_HoverELoose_hltEG30HEFilterLooseHoverE_phi);
   fChain->SetBranchAddress("trig_HLT_Photon30_HoverELoose_hltEG30HEFilterLooseHoverE_et", &trig_HLT_Photon30_HoverELoose_hltEG30HEFilterLooseHoverE_et, &b_trig_HLT_Photon30_HoverELoose_hltEG30HEFilterLooseHoverE_et);
   fChain->SetBranchAddress("trig_HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142_accept", &trig_HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142_accept, &b_trig_HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142_accept);
   fChain->SetBranchAddress("trig_AlCa_RPCMuonNormalisation_accept", &trig_AlCa_RPCMuonNormalisation_accept, &b_trig_AlCa_RPCMuonNormalisation_accept);
   fChain->SetBranchAddress("trig_HLT_Mu18_Mu9_SameSign_accept", &trig_HLT_Mu18_Mu9_SameSign_accept, &b_trig_HLT_Mu18_Mu9_SameSign_accept);
   fChain->SetBranchAddress("trig_HLT_Mu18_Mu9_SameSign_DZ_accept", &trig_HLT_Mu18_Mu9_SameSign_DZ_accept, &b_trig_HLT_Mu18_Mu9_SameSign_DZ_accept);
   fChain->SetBranchAddress("trig_HLT_Mu18_Mu9_accept", &trig_HLT_Mu18_Mu9_accept, &b_trig_HLT_Mu18_Mu9_accept);
   fChain->SetBranchAddress("trig_HLT_Mu18_Mu9_DZ_accept", &trig_HLT_Mu18_Mu9_DZ_accept, &b_trig_HLT_Mu18_Mu9_DZ_accept);
   fChain->SetBranchAddress("trig_HLT_Mu20_Mu10_SameSign_accept", &trig_HLT_Mu20_Mu10_SameSign_accept, &b_trig_HLT_Mu20_Mu10_SameSign_accept);
   fChain->SetBranchAddress("trig_HLT_Mu20_Mu10_SameSign_DZ_accept", &trig_HLT_Mu20_Mu10_SameSign_DZ_accept, &b_trig_HLT_Mu20_Mu10_SameSign_DZ_accept);
   fChain->SetBranchAddress("trig_HLT_Mu20_Mu10_accept", &trig_HLT_Mu20_Mu10_accept, &b_trig_HLT_Mu20_Mu10_accept);
   fChain->SetBranchAddress("trig_HLT_Mu20_Mu10_DZ_accept", &trig_HLT_Mu20_Mu10_DZ_accept, &b_trig_HLT_Mu20_Mu10_DZ_accept);
   fChain->SetBranchAddress("trig_HLT_Mu23_Mu12_SameSign_accept", &trig_HLT_Mu23_Mu12_SameSign_accept, &b_trig_HLT_Mu23_Mu12_SameSign_accept);
   fChain->SetBranchAddress("trig_HLT_Mu23_Mu12_SameSign_DZ_accept", &trig_HLT_Mu23_Mu12_SameSign_DZ_accept, &b_trig_HLT_Mu23_Mu12_SameSign_DZ_accept);
   fChain->SetBranchAddress("trig_HLT_Mu23_Mu12_accept", &trig_HLT_Mu23_Mu12_accept, &b_trig_HLT_Mu23_Mu12_accept);
   fChain->SetBranchAddress("trig_HLT_Mu23_Mu12_DZ_accept", &trig_HLT_Mu23_Mu12_DZ_accept, &b_trig_HLT_Mu23_Mu12_DZ_accept);
   Notify();
}

Bool_t Nochecl::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Nochecl::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Nochecl::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Nochecl_cxx
