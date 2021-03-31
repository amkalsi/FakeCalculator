// change done 7 July 2018 -- add >= 1 bjets requirement
#include "LowMTRegion_2018.h"
using namespace std;

LowMTRegion_2018::LowMTRegion_2018(const edm::ParameterSet& iConfig)

{

	isoTau = iConfig.getParameter<string>("isolationTau");
	InputFile = iConfig.getParameter<string>("InputFile");
	DYSample = iConfig.getParameter<string>("DYSample");                                                                                   

	MuonIDIsoSF = iConfig.getParameter<string>("MuonIDIsoSF");
	MuonTriggerSF = iConfig.getParameter<string>("MuonTriggerSF");
	FakeRateSSfile = iConfig.getParameter<string>("FakeRateSSfile");
	FakeRateDYJetsfile = iConfig.getParameter<string>("FakeRateDYJetsfile");
	size_t npos = -1;
	if( DYSample.find("Egamma") != npos || DYSample.find("SE") != npos || DYSample.find("SingleElectron") != npos || DYSample.find("singlephoton") != npos ) { isdata = true; } else { isdata = false;}
	std::cout<<"file:"<<DYSample<<"\tis data \t"<<isdata<<std::endl;
	muscale=1.00;


	isOS = iConfig.getParameter<bool>("isOS");
	isSS = iConfig.getParameter<bool>("isSS");

}


LowMTRegion_2018::~LowMTRegion_2018()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
	void
LowMTRegion_2018::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

	using namespace edm;
	nEventsRaw = nEventsStored = mc_nEventsWeighted = nEventsiihe = 0;
	std::cout<<"datafile::"<<InputFile<<std::endl;
	TFile *file_in=TFile::Open(InputFile.c_str(),"READ");


	/// tree to get META data ///
	if(!isdata){
		TTreeReader evt_reader("meta", file_in);
		TTreeReaderValue<float> evt_nEventsRaw(evt_reader, "nEventsRaw");
		TTreeReaderValue<float> evt_mc_nEventsWeighted(evt_reader, "mc_nEventsWeighted");
		while (evt_reader.Next()) {
			nEventsRaw = nEventsRaw+ (*evt_nEventsRaw);
			mc_nEventsWeighted = mc_nEventsWeighted + (*evt_mc_nEventsWeighted);
		}
	} else {
		TTreeReader evt_reader("meta", file_in);

		TTreeReaderValue<float> evt_nEventsRaw(evt_reader, "nEventsRaw");
		while (evt_reader.Next()) {
			nEventsRaw = nEventsRaw + (*evt_nEventsRaw);
		}

	}

	h_Fill_RawEvents->SetBinContent(1, nEventsRaw);
	h_Fill_WeightedEvents->SetBinContent(1,mc_nEventsWeighted);

	TTree* treePtr = (TTree*) file_in->Get("IIHEAnalysis"); 
	tree = new IIHEAnalysis (treePtr);   

	std::cout << "entries in tree : "<< tree->GetEntries() << std::endl;
	for (int iEntry = 0; iEntry < tree->GetEntries(); iEntry++)
	{     
		tree->GetEntry(iEntry);
		h_Events_Before_Skim->Fill(1.);
		// gen filters
		reject_event=false;
		size_t npos = -1;
		Mass=0.;



		if(isdata ){
			if( !( tree->trig_HLT_Ele32_WPTight_Gsf_accept || tree->trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_accept )) continue;
		}
		if(!isdata){
			if( !( tree->trig_HLT_Ele32_WPTight_Gsf_accept || tree->trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_accept )) continue;
		}

		if(isdata) {h_Events_After_Trigger->Fill(1.);}
		if(!isdata) {h_Events_After_Trigger->Fill(1.,tree->mc_w_sign);}


		if(!tree->trig_Flag_goodVertices_accept) continue;
		if(!tree->trig_Flag_globalSuperTightHalo2016Filter_accept) continue;
		if(!tree->trig_Flag_HBHENoiseFilter_accept) continue;
		if(!tree->trig_Flag_HBHENoiseIsoFilter_accept) continue;
		if(!tree->trig_Flag_EcalDeadCellTriggerPrimitiveFilter_accept) continue;
		if(!tree->trig_Flag_BadPFMuonFilter_accept) continue;
		// new filter -- updated in miniAOD
		if(!tree->trig_Flag_ecalBadCalibReduced) continue;

		//		if(!tree->trig_Flag_BadChargedCandidateFilter_accept) continue;
		//		if(isdata) if(!tree->trig_Flag_eeBadScFilter_accept) continue;

		if(isdata) h_Events_After_MetFilters->Fill(1.);
		if(!isdata) h_Events_After_MetFilters->Fill(1., tree->mc_w_sign);


		// loop for OS and SS events
		// events filters

		sumweight =0;
		weighthis=1.;

		if(!isdata) {pu_weight =PU_2018_Rereco::MC_pileup_weight(tree->mc_trueNumInteractions, "MC", "Data_2018AtoD"); } else {pu_weight = 1;}

		if(!isdata) weighthis = tree->mc_w_sign*pu_weight; 
		if(isdata) { weighthis = 1.;}

		if(!isdata &&  ( DYSample.find("TT") != npos || DYSample.find("TTinc") != npos || DYSample.find("TT_2l2nu") != npos || DYSample.find("TT_had") != npos || DYSample.find("TT_semilep") != npos )) {

			double weight_top1 = 1.0;
			//double weight_top2 = 1.0;
			int tops=0;

			double tmp_t1 , tmp_t2;
			tmp_t1=0;
			tmp_t2=0;
			for (unsigned int iLHE = 0; iLHE < tree->LHE_Pt->size(); ++iLHE) {

				if ( tree->LHE_pdgid->at(iLHE) == 6) { tmp_t1 = exp(0.0615-0.0005*tree->LHE_Pt->at(iLHE)); }

				if ( tree->LHE_pdgid->at(iLHE) == -6) { tmp_t2 = exp(0.0615-0.0005*tree->LHE_Pt->at(iLHE)); }

				tops++;
			}
			if(tops == 2) weight_top1 = sqrt(tmp_t1 * tmp_t2);


			weighthis = weighthis* weight_top1;


		}
//		std::cout<<"debug1:"<< std::endl;
		vector<unsigned int> ele_indexes, tau_indexes;
		ele_indexes.clear();
		tau_indexes.clear();

		for (unsigned int dau1index = 0; dau1index < tree->gsf_caloEnergy->size(); dau1index++){
			float ET1 = tree->gsf_caloEnergy->at(dau1index)*sin(2.*atan(exp(-1.*tree->gsf_eta->at(dau1index)))) ;


			TLorentzVector DauEle ;
			DauEle.SetPtEtaPhiM(ET1, tree->gsf_eta->at(dau1index), tree->gsf_phi->at(dau1index),m_el);
			if(!( ET1 > 50.  && fabs(tree->gsf_sc_eta->at(dau1index)) < 2.5) ) continue;
			if(!(tree->gsf_VID_heepElectronID_HEEPV70->at(dau1index))) continue;

			TVector2 METcorr;
			METcorr.SetMagPhi(tree->MET_nominal_Pt,tree->MET_nominal_phi);
			double mt_cut = mTCalculation(METcorr.Px(), METcorr.Py(), DauEle.Px(), DauEle.Py(), DauEle.Pt());
			if(mt_cut > 120.) continue;
			ele_indexes.push_back(dau1index);

		}

		if(ele_indexes.size() == 0 ) continue;

  //              std::cout<<"debug1 part1:"<< std::endl;

		if(isdata) h_Fill_ElePass->Fill(1.);
		if(!isdata) h_Fill_ElePass->Fill(1., tree->mc_w_sign);

		for (unsigned int dau2index = 0; dau2index < tree->tau_pt->size(); dau2index++){

			if(!(tree->tau_pt->at(dau2index) > 30.)) continue;
			if(!(fabs(tree->tau_eta->at(dau2index)) < 2.3)) continue;
			if(!(tree->tau_decayModeFinding->at(dau2index) > 0.5)) continue;
			if(!(tree->tau_byVLooseIsolationMVArun2017v2DBnewDMwLT2017->at(dau2index) > 0.5)) continue;
			if(!(tree->tau_againstMuonLoose3->at(dau2index) > 0.5)) continue;
			if(!(tree->tau_againstElectronTightMVA6->at(dau2index) > 0.5)) continue;

			TLorentzVector DauTau ;
			DauTau.SetPxPyPzE(tree->tau_px->at(dau2index), tree->tau_py->at(dau2index), tree->tau_pz->at(dau2index),tree->tau_energy->at(dau2index));
			tau_indexes.push_back(dau2index);
		}

		if(tau_indexes.size() == 0 ) continue;

    //            std::cout<<"debug1 part2:"<< std::endl;

		if(isdata) h_Fill_TauPass->Fill(1.);
		if(!isdata) h_Fill_TauPass->Fill(1., tree->mc_w_sign);

		int muon_counter = 0;
		for(unsigned int imuon=0; imuon < tree->mu_ibt_pt->size(); imuon++) {
			if(!(tree->mu_ibt_pt->at(imuon) > 30.)) continue;
			if(!(fabs(tree->mu_ibt_eta->at(imuon)) < 2.4)) continue;
			if(!(tree->mu_isHighPtMuon->at(imuon) == 1)) continue;//mu_isHighPtMuon->at(imuon))) continue;
			if(!(tree->mu_isoTrackerBased03->at(imuon) < 0.15)) continue;
			bool matched_not_to_lepton=false;
			bool matched_not_to_tau=false;
			TLorentzVector DauMu;
			DauMu.SetPtEtaPhiM(tree->mu_ibt_pt->at(imuon), tree->mu_ibt_eta->at(imuon), tree->mu_ibt_phi->at(imuon), m_mu);
			// electron loop
			for (unsigned int dau1index = 0; dau1index < ele_indexes.size(); dau1index++){
				float ET1 = tree->gsf_caloEnergy->at(ele_indexes.at(dau1index))*sin(2.*atan(exp(-1.*tree->gsf_eta->at(ele_indexes.at(dau1index)))));

				TLorentzVector DauEle ;
				DauEle.SetPtEtaPhiM(ET1, tree->gsf_eta->at(ele_indexes.at(dau1index)), tree->gsf_phi->at(ele_indexes.at(dau1index)),m_el);
				if(DauMu.DeltaR(DauEle) > 0.5){
					matched_not_to_lepton = true;
				}
			}
			// tau lloop
			for (unsigned int dau2index = 0; dau2index < tau_indexes.size(); dau2index++){
				TLorentzVector DauTau ;
				DauTau.SetPxPyPzE(tree->tau_px->at(tau_indexes.at(dau2index)), tree->tau_py->at(tau_indexes.at(dau2index)), tree->tau_pz->at(tau_indexes.at(dau2index)),tree->tau_energy->at(tau_indexes.at(dau2index)));
				if(DauMu.DeltaR(DauTau) > 0.5) {
					matched_not_to_tau = true;

				}

			}
			// DR matching
			if(matched_not_to_lepton && matched_not_to_tau) muon_counter++;
		}
		if(muon_counter != 0 ) continue;
//		std::cout<<"debug2:"<< std::endl;


		if(isdata) h_Fill_MuonExtra->Fill(1.);
		if(!isdata) h_Fill_MuonExtra->Fill(1., tree->mc_w_sign);

		TLorentzVector FirstObj(0.,0.,0.,0), SecondObj(0.,0.,0.,0);

		TLorentzVector FirstObj_FF(0.,0.,0.,0), SecondObj_FF(0.,0.,0.,0);
		first_index_FF= -1;
		second_index_FF = -1;
		bool found_one_ele = false;

		for (unsigned int dau1index = 0; dau1index < tree->gsf_caloEnergy->size(); dau1index++){
			if(found_one_ele) continue;

			float ET1 = tree->gsf_caloEnergy->at(dau1index)*sin(2.*atan(exp(-1.*tree->gsf_eta->at(dau1index)))) ;

			TLorentzVector DauEle ;
			DauEle.SetPtEtaPhiM(ET1, tree->gsf_eta->at(dau1index), tree->gsf_phi->at(dau1index),m_el);
			if(!( ET1 > 50.  && fabs(tree->gsf_sc_eta->at(dau1index)) < 2.5) ) continue;
			if(!(tree->gsf_VID_heepElectronID_HEEPV70->at(dau1index))) continue;
			TVector2 METcorr;
			METcorr.SetMagPhi(tree->MET_nominal_Pt,tree->MET_nominal_phi);
			double mt_cut = mTCalculation(METcorr.Px(), METcorr.Py(), DauEle.Px(), DauEle.Py(), DauEle.Pt());
			if(mt_cut > 120.) continue;



			for (unsigned int dau2index = 0; dau2index < tree->tau_pt->size(); dau2index++){

				if(!(tree->tau_pt->at(dau2index) > 30.)) continue;
				if(!(fabs(tree->tau_eta->at(dau2index)) < 2.3)) continue;


				if(!(tree->tau_decayModeFinding->at(dau2index) > 0.5)) continue;
				if(!(tree->tau_byVLooseIsolationMVArun2017v2DBnewDMwLT2017->at(dau2index) > 0.5)) continue;
				if(!(tree->tau_againstMuonLoose3->at(dau2index) > 0.5)) continue;
				if(!(tree->tau_againstElectronTightMVA6->at(dau2index) > 0.5)) continue;

				TLorentzVector DauTau ;
				DauTau.SetPxPyPzE(tree->tau_px->at(dau2index), tree->tau_py->at(dau2index), tree->tau_pz->at(dau2index),tree->tau_energy->at(dau2index));

				if(!(DauEle.DeltaR(DauTau) > 0.5 ))  continue;

				bool applySF;
				applySF= false;
				int genindex=-1;
				unsigned int matchgen_obj_type;
				if(!isdata) { if( (MatchingToGenTaus(DauTau, genindex) && genindex!=-1) || matchedToGenObjetcs(DauTau, matchgen_obj_type)) {applySF = true;} }
				if(isdata) {applySF= true;}
				if(!applySF) continue;


				bool matched_to_reco_jet=false;
				TLorentzVector jet_p4(0.,0.,0.,0.);
				for (unsigned int ijet = 0; ijet < tree->jet_pt->size(); ijet++){
					if( !(tree->jet_pt->at(ijet) > 0.) ) continue;
					if(!(fabs(tree->jet_eta->at(ijet)) < 2.3)) continue;
					if(!(tree->jet_isJetID_2018->at(ijet))) continue;
					TLorentzVector jet_p4_tmp;
					jet_p4_tmp.SetPtEtaPhiE(tree->jet_pt->at(ijet), tree->jet_eta->at(ijet), tree->jet_phi->at(ijet), tree->jet_energy->at(ijet));
					if(!(DauTau.DeltaR(jet_p4_tmp) < 0.2)) continue;
					matched_to_reco_jet=true;
					jet_p4=jet_p4_tmp;
					break;

				}

				if(!(matched_to_reco_jet)) continue;
				int j_dm = -1, k_eta = -1;
				if (tree->tau_decayMode->at(dau2index) == 0) {
					j_dm = 0;
				}
				else if (tree->tau_decayMode->at(dau2index) == 1 || tree->tau_decayMode->at(dau2index) == 2) {
					j_dm = 1;
				}
				else if (tree->tau_decayMode->at(dau2index) == 10 || tree->tau_decayMode->at(dau2index) == 11) {
					j_dm = 2;
				}

				if (fabs(tree->tau_eta->at(dau2index)) < 1.46) {
					k_eta = 0;
				}
				else if (fabs(tree->tau_eta->at(dau2index)) > 1.558) {
					k_eta = 1;
				}
				else {
					continue;
				}

				double weight_sd = 1.;
				double fake_wt = 1.;
				//				if( (!isdata) )  {  if( fabs(tree->gsf_sc_eta->at(dau1index)) < 1.4442) {weight_sd = eTrigger_barrel->Eval(DauEle.Pt()); } else { weight_sd = eTrigger_endcap->Eval(DauEle.Pt()); }}



				//				if(!isdata) { fake_wt =  MutoTauFR(DauTau);}
//std::cout<<"debug3:"<< std::endl;

				double final_weight = weight_sd*fake_wt*weighthis;
				found_one_ele = true;

				if(tree->gsf_charge->at(dau1index) * tree->tau_charge->at(dau2index) == -1 ) {
					if (tree->tau_byTightIsolationMVArun2017v2DBnewDMwLT2017->at(dau2index) > 0.5) {
						hh[0][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index), jet_p4.Pt(), final_weight);
						h1[0][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index), final_weight);
						h2[0][j_dm][k_eta]->Fill(jet_p4.Pt(), final_weight);
						h3[0][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index)/jet_p4.Pt(), final_weight);

						if(tree->tau_pt->at(dau2index) < 150.) { h4[0][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index), tree->tau_pt->at(dau2index)/jet_p4.Pt(), final_weight); }
						if(tree->tau_pt->at(dau2index)  >= 150. ) { h5[0][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index), tree->tau_pt->at(dau2index)/jet_p4.Pt(), final_weight); }

						//						if(tree->tau_pt->at(dau2index) >= 30. && tree->tau_pt->at(dau2index)  < 50.) { h4[0][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index)/jet_p4.Pt(), final_weight); }
						//						if(tree->tau_pt->at(dau2index)  >= 50. && tree->tau_pt->at(dau2index)  < 70.) { h5[0][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index)/jet_p4.Pt(), final_weight); }
						//				if(tree->tau_pt->at(dau2index)  >= 70. && tree->tau_pt->at(dau2index)  < 100.) { h6[0][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index)/jet_p4.Pt(), final_weight); }
						//				if(tree->tau_pt->at(dau2index)  >= 100. ) { h7[0][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index)/jet_p4.Pt(), final_weight); }


					}
					if ((tree->tau_byTightIsolationMVArun2017v2DBnewDMwLT2017->at(dau2index) < 0.5) && (tree->tau_byVLooseIsolationMVArun2017v2DBnewDMwLT2017->at(dau2index) > 0.5)) {
						hh[1][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index), jet_p4.Pt(), final_weight);
						h1[1][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index), final_weight);
						h2[1][j_dm][k_eta]->Fill(jet_p4.Pt(), final_weight);
						h3[1][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index)/jet_p4.Pt(), final_weight);

						if(tree->tau_pt->at(dau2index) < 150. ) { h4[1][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index), tree->tau_pt->at(dau2index)/jet_p4.Pt(), final_weight); }
						if(tree->tau_pt->at(dau2index)  >= 150. ) { h5[1][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index), tree->tau_pt->at(dau2index)/jet_p4.Pt(), final_weight); }

						//						if(tree->tau_pt->at(dau2index) >= 30. && tree->tau_pt->at(dau2index)  < 50.) { h4[1][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index)/jet_p4.Pt(), final_weight); }
						//						if(tree->tau_pt->at(dau2index)  >= 50. && tree->tau_pt->at(dau2index)  < 70.) { h5[1][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index)/jet_p4.Pt(), final_weight); }
						//				if(tree->tau_pt->at(dau2index)  >= 70. && tree->tau_pt->at(dau2index)  < 100.) { h6[1][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index)/jet_p4.Pt(), final_weight); }
						//				if(tree->tau_pt->at(dau2index)  >= 100. ) { h7[1][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index)/jet_p4.Pt(), final_weight); }


					}

				} else if (tree->gsf_charge->at(dau1index) * tree->tau_charge->at(dau2index) == 1) {


					if (tree->tau_byTightIsolationMVArun2017v2DBnewDMwLT2017->at(dau2index) > 0.5) {
						hh[2][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index), jet_p4.Pt(), final_weight);
						h1[2][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index), final_weight);
						h2[2][j_dm][k_eta]->Fill(jet_p4.Pt(), final_weight);
						h3[2][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index)/jet_p4.Pt(), final_weight);

						if(tree->tau_pt->at(dau2index)   < 150.) { h4[2][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index), tree->tau_pt->at(dau2index)/jet_p4.Pt(), final_weight); }
						if(tree->tau_pt->at(dau2index)  >= 150. ) { h5[2][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index), tree->tau_pt->at(dau2index)/jet_p4.Pt(), final_weight); }

						//						if(tree->tau_pt->at(dau2index) >= 30. && tree->tau_pt->at(dau2index)  < 50.) { h4[2][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index)/jet_p4.Pt(), final_weight); }
						//						if(tree->tau_pt->at(dau2index)  >= 50. && tree->tau_pt->at(dau2index)  < 70.) { h5[2][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index)/jet_p4.Pt(), final_weight); }
						//				if(tree->tau_pt->at(dau2index)  >= 70. && tree->tau_pt->at(dau2index)  < 100.) { h6[2][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index)/jet_p4.Pt(), final_weight); }
						//				if(tree->tau_pt->at(dau2index)  >= 100. ) { h7[2][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index)/jet_p4.Pt(), final_weight); }


					}
					if ((tree->tau_byTightIsolationMVArun2017v2DBnewDMwLT2017->at(dau2index) < 0.5) && (tree->tau_byVLooseIsolationMVArun2017v2DBnewDMwLT2017->at(dau2index) > 0.5)) {
						hh[3][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index), jet_p4.Pt(), final_weight);
						h1[3][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index), final_weight);
						h2[3][j_dm][k_eta]->Fill(jet_p4.Pt(), final_weight);
						h3[3][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index)/jet_p4.Pt(), final_weight);

						if(tree->tau_pt->at(dau2index)  < 150.) { h4[3][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index), tree->tau_pt->at(dau2index)/jet_p4.Pt(), final_weight); }
						if(tree->tau_pt->at(dau2index)  >= 150. ) { h5[3][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index), tree->tau_pt->at(dau2index)/jet_p4.Pt(), final_weight); }

						//				if(tree->tau_pt->at(dau2index) >= 30. && tree->tau_pt->at(dau2index)  < 50.) { h4[3][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index)/jet_p4.Pt(), final_weight); }
						//				if(tree->tau_pt->at(dau2index)  >= 50. && tree->tau_pt->at(dau2index)  < 70.) { h5[3][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index)/jet_p4.Pt(), final_weight); }
						//				if(tree->tau_pt->at(dau2index)  >= 70. && tree->tau_pt->at(dau2index)  < 100.) { h6[3][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index)/jet_p4.Pt(), final_weight); }
						//				if(tree->tau_pt->at(dau2index)  >= 100. ) { h7[3][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index)/jet_p4.Pt(), final_weight); }


					}


				}


			}
		}


	}
	delete tree;

	file_in->Close();



#ifdef THIS_IS_AN_EVENT_EXAMPLE
	Handle<ExampleData> pIn;
	iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
	ESHandle<SetupData> pSetup;
	iSetup.get<SetupRecord>().get(pSetup);
#endif

}
// ------------ method called once each job just before starting event loop  ------------


	void 
LowMTRegion_2018::beginJob()
{
	//	fIn = TFile::Open(InputFile.c_str());
	std::cout<<"InputFile::"<<InputFile<< std::endl;
	//	file_db1.open(InputFile);                                                                                       
  float xbins_left[] = {0, 30, 40, 50, 60, 70, 80, 100, 150};
        float ybins_left[]= {0, 0.5, 0.6, 0.7, 0.75, 0.8, 1., 3.};
        float xbins_right[] = {150, 1000};
        float ybins_right[] = {0, 0.7, 0.8, 0.85, 1., 3.};

//	float xbins_left[] = {0, 30, 40, 50, 60, 70, 80, 100, 150};
//	float ybins_left[]= {0, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 1., 3.};
//	float xbins_right[] = {150, 1000};
//	float ybins_right[] = {0, 0.6, 0.7, 0.8, 0.85, 1., 3.};


	//	float xbins_left[] = {0, 30, 40, 50, 60, 70, 80, 100, 150};
	//	float ybins_left[]= {0, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 1., 3.};
	//	float xbins_right[] = {150, 1000};
	//	float ybins_right[] = {0, 0.6, 0.7, 0.8, 0.85, 1., 3.};
	const int nxaxis_left = sizeof(xbins_left)/sizeof(float);
	const int nyaxis_left = sizeof(ybins_left)/sizeof(float);

	const int nxaxis_right = sizeof(xbins_right)/sizeof(float);
	const int nyaxis_right = sizeof(ybins_right)/sizeof(float);


	h_Count_Taus = fs->make<TH1D>("h_Count_Taus","h_Count_Taus",10,0,10); 
	h_Count_Taus->Sumw2();
	h_Events_Before_Skim = fs->make<TH1D>("h_Events_Before_Skim","h_Events_Before_Skim",2,0,2);
	h_Events_After_Skim = fs->make<TH1D>("h_Events_After_Skim","h_Events_After_Skim",2,0,2);
	h_Events_After_GenFilter = fs->make<TH1D>("h_Events_After_GenFilter","h_Events_After_GenFilter",2,0,2);
	h_Fill_Mass_Gen_toChk = fs->make<TH1D>("h_Fill_Mass_Gen_toChk","h_Fill_Mass_Gen_toChk",3000,0,3000);
	h_Fill_Mass_Gen_toChk_before = fs->make<TH1D>("h_Fill_Mass_Gen_toChk_before","h_Fill_Mass_Gen_toChk_before",3000,0,3000);

	h_Fill_Mass_Gen_toChk->Sumw2();
	h_Fill_Mass_Gen_toChk_before->Sumw2();
	h_Events_Before_Skim->Sumw2();
	h_Events_After_Skim->Sumw2();
	h_Events_After_GenFilter->Sumw2();

	h_Fill_RawEvents  = fs->make<TH1D>("h_Fill_RawEvents","h_Fill_RawEvents",2,0,2);
	h_Fill_WeightedEvents  = fs->make<TH1D>("h_Fill_WeightedEvents","h_Fill_WeightedEvents",2,0,2);
	h_Fill_RawEvents->Sumw2();
	h_Fill_WeightedEvents->Sumw2();

	h_Events_After_Trigger = fs->make<TH1D>("h_Events_After_Trigger","h_Events_After_Trigger",2,0,2);
	h_Fill_TauPass = fs->make<TH1D>("h_Fill_TauPass","h_Fill_TauPass",2,0,2);
	h_Fill_ElePass = fs->make<TH1D>("h_Fill_ElePass","h_Fill_ElePass",2,0,2);
	h_Events_After_MetFilters  = fs->make<TH1D>("h_Events_After_MetFilters","h_Events_After_MetFilters",2,0,2);
	h_Fill_MuonExtra  = fs->make<TH1D>("h_Fill_MuonExtra","h_Fill_MuonExtra",2,0,2);
	h_Events_After_Trigger->Sumw2();
	h_Fill_TauPass->Sumw2();
	h_Fill_ElePass->Sumw2();
	h_Events_After_MetFilters->Sumw2();
	h_Fill_MuonExtra->Sumw2();

	//	h_Fill_bjet_pt_ETauFake = fs->make<TH1D>("h_Fill_bjet_pt_ETauFake","h_Fill_bjet_pt_ETauFake", 1000,0,1000);
	//	h_Fill_bjet_eta_ETauFake = fs->make<TH1D>("h_Fill_bjet_eta_ETauFake","h_Fill_bjet_eta_ETauFake",100,-5,5);
	//	h_Fill_bjet_disc_ETauFake = fs->make<TH1D>("h_Fill_bjet_disc_ETauFake","h_Fill_bjet_disc_ETauFake",100,0,1);

	//	h_Fill_bjet_pt_ETauFake->Sumw2();
	//	h_Fill_bjet_eta_ETauFake->Sumw2();
	//	h_Fill_bjet_disc_ETauFake->Sumw2();



	string h_names[] = {"taupt_jetpt_pass_OS","taupt_jetpt_fail_OS","taupt_jetpt_pass_SS","taupt_jetpt_fail_SS"}; 
	unsigned int numH_name = sizeof(h_names)/sizeof(string);

	string dms[] = {"DM0", "DM1", "DM10"};
	unsigned int numH_dms = sizeof(dms)/sizeof(string);

	string etaarr[] = {"barrel","endcap"};
	unsigned int num_etaarr = sizeof(etaarr)/sizeof(string);

	for (unsigned int i = 0; i< numH_name; ++i) {
		for (unsigned int k = 0; k< numH_dms ; ++k) {
			for (unsigned int l = 0; l< num_etaarr; ++l) {
				TString nname = h_names[i]+"_"+dms[k]+"_"+etaarr[l];

				TString nname1 = h_names[i]+"_"+dms[k]+"_"+etaarr[l]+"tau_pt";
				TString nname2 = h_names[i]+"_"+dms[k]+"_"+etaarr[l]+"jet_pt";
				TString nname3 = h_names[i]+"_"+dms[k]+"_"+etaarr[l]+"_ratio";
				TString nname4 = h_names[i]+"_"+dms[k]+"_"+etaarr[l]+"_ratio_1";
				TString nname5 = h_names[i]+"_"+dms[k]+"_"+etaarr[l]+"_ratio_2";
				//	TString nname6 = h_names[i]+"_"+dms[k]+"_"+etaarr[l]+"_ratio_3";
				//	TString nname7 = h_names[i]+"_"+dms[k]+"_"+etaarr[l]+"_ratio_4";
				//	TString nname8 = h_names[i]+"_"+dms[k]+"_"+etaarr[l]+"_ratio_5";


				hh[i][k][l] =  fs->make<TH2F>(nname, nname, 1000, 0, 1000, 1000, 0, 1000) ;
				h1[i][k][l] = fs->make<TH1F>(nname1, nname1, 1000, 0, 1000) ;
				h2[i][k][l] = fs->make<TH1F>(nname2, nname2, 1000, 0, 1000) ;
				h3[i][k][l] = fs->make<TH1F>(nname3, nname3, 1000, 0, 10) ;
				//	h4[i][k][l] = fs->make<TH1F>(nname4, nname4, 1000, 0, 10) ;
				//	h5[i][k][l] = fs->make<TH1F>(nname5, nname5, 1000, 0, 10) ;
				//	h6[i][k][l] = fs->make<TH1F>(nname6, nname6, 1000, 0, 10) ;
				//	h7[i][k][l] = fs->make<TH1F>(nname7, nname7, 1000, 0, 10) ;
				//	h8[i][k][l] = fs->make<TH1F>(nname8, nname8, 1000, 0, 10) ;

				h4[i][k][l] = fs->make<TH2F>(nname4, nname4, nxaxis_left-1, xbins_left,nyaxis_left-1, ybins_left) ;
				h5[i][k][l] = fs->make<TH2F>(nname5, nname5, nxaxis_right-1, xbins_right, nyaxis_right-1, ybins_right) ;


				hh[i][k][l]->Sumw2();
				h1[i][k][l]->Sumw2();
				h2[i][k][l]->Sumw2();
				h3[i][k][l]->Sumw2();
				h4[i][k][l]->Sumw2();
				h5[i][k][l]->Sumw2();
				//	h6[i][k][l]->Sumw2();
				//	h7[i][k][l]->Sumw2();
				//	h8[i][k][l]->Sumw2();

			}
		}
	}


}


// ------------ method called once each job just after ending the event loop  ------------
	void 
LowMTRegion_2018::endJob() 
{

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
LowMTRegion_2018::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}








bool LowMTRegion_2018::OverLap05(TLorentzVector l1 , TLorentzVector l2, float conesize) {
	if(dR(l1.Eta(), l1.Phi(), l2.Eta(), l2.Phi()) <= conesize) return true;
	else return false;
}  

float LowMTRegion_2018::deltaPhi( float a, float b) {
	float result = a-b;
	while (result > M_PI) result -= 2* M_PI;
	while (result <= -M_PI) result += 2* M_PI;
	return (fabs(result));

} 

float LowMTRegion_2018::dR(float l1eta, float l1phi, float l2eta, float l2phi ) {
	float deta = l1eta - l2eta;
	float dphi = deltaPhi(l1phi,l2phi);
	return sqrt(deta*deta + dphi*dphi);
}




float LowMTRegion_2018::mTCalculation(float metx, float mety, float mupx, float mupy, float mupt){
	float mt = -1;
	float pX = mupx+metx;
	float pY = mupy+mety;
	float et = mupt + TMath::Sqrt(metx*metx + mety*mety);
	mt = TMath::Sqrt(et*et-(pX*pX + pY*pY));
	return mt;

}      

/*
   float LowMTRegion_2018::PZetaVis( int muindex, int tauindex){
   float pzetavis;
   pzetavis = 999;
   TLorentzVector tau, mu;  
   tau.SetPtEtaPhiE(tree->tau_pt->at(tauindex),tree->tau_eta->at(tauindex),tree->tau_phi->at(tauindex),tree->tau_energy->at(tauindex));
   mu.SetPtEtaPhiM(tree->->at(muindex), tree->mu_gt_eta->at(muindex), tree->mu_gt_phi->at(muindex),mu_mass);
   float zetax = TMath::Cos(mu.Phi()) + TMath::Cos(tau.Phi()) ;
   float zetay = TMath::Sin(mu.Phi()) + TMath::Sin(tau.Phi()) ;
   float zetaR = TMath::Sqrt(pow(zetax,2) + pow(zetay,2));
   zetax = zetax/zetaR;
   zetay = zetay/zetaR;

   float visPx = mu.Px() + tau.Px() ;
   float visPy = mu.Py() + tau.Py() ;

   pzetavis = visPx*zetax + visPy*zetay;
   return pzetavis;

   }

   float LowMTRegion_2018::PZeta(int muindex, int tauindex , float metpx, float metpy){
   float pzeta;
   pzeta = 999;
   TLorentzVector tau, mu;                                         
   tau.SetPxPyPzE(tree->tau_px->at(tauindex),tree->tau_py->at(tauindex),tree->tau_pz->at(tauindex),tree->tau_energy->at(tauindex));
   mu.SetPtEtaPhiM(tree->mu_gt_pt->at(muindex), tree->mu_gt_eta->at(muindex), tree->mu_gt_phi->at(muindex),mu_mass);

   float zetax = TMath::Cos(mu.Phi()) + TMath::Cos(tau.Phi()) ;
   float zetay = TMath::Sin(mu.Phi()) + TMath::Sin(tau.Phi()) ;
   float zetaR = TMath::Sqrt(pow(zetax,2) + pow(zetay,2)); 
   zetax = zetax/zetaR;                                    
   zetay = zetay/zetaR;         

   float vPx = mu.Px() + tau.Px()+metpx ;
   float vPy = mu.Py() + tau.Py()+metpy ;

   pzeta = vPx*zetax + vPy*zetay;
   return pzeta;


   }    
   */

void LowMTRegion_2018::HistoFiller(TH1F *histo, double value, double weight){                                    

	histo->Fill(value, weight);

}  



void LowMTRegion_2018::initializePileupInfo() {

	// Filenames must be c_strings below. Here is the conversion from strings to c_strings
	//         //   // As you can see above cstr1 corresponds to MC and cstr2 corresponds to data.
	//
	char * cstr1;
	char * cstr2;

	//Filenames must be c_strings below. Here is the conversion from strings
	cstr1 = new char [MCHistosForPU.size()+1];
	strcpy (cstr1, MCHistosForPU.c_str());
	cstr2 = new char [DataHistosForPU.size()+1];
	strcpy (cstr2, DataHistosForPU.c_str());

	TFile *file1 = TFile::Open(cstr1);
	TH1* histmc = dynamic_cast<TH1*>(file1->Get("pileup"));
	if(!histmc) {throw std::runtime_error("failed to extract histogram");}
	for(int bin=0; bin<=(histmc->GetXaxis()->GetNbins() + 1); bin++) {
		hPUmc->SetBinContent(bin,histmc->GetBinContent(bin));
	}      
	file1->Close();

	TFile *file2 = TFile::Open(cstr2);
	//file:/uscms_data/d3/aman30/RUN2_50ns/CMSSW_7_4_15/src/Plots/Code/Data50ns.root");
	TH1* histdata = dynamic_cast<TH1*>(file2->Get("pileup"));
	if(!histdata) {throw std::runtime_error("failed to extract histogram");}
	for(int bin=0; bin<=(histdata->GetXaxis()->GetNbins() + 1); bin++) {
		hPUdata->SetBinContent(bin,histdata->GetBinContent(bin));
	}
	file2->Close();
}        

double LowMTRegion_2018::getPileupWeight(float ntruePUInt) {
	int bin;
	double MCintegral;
	double MCvalue;
	double Dataintegral;
	double Datavalue;

	bin = hPUmc->GetBin(ntruePUInt+1);
	MCvalue = hPUmc->GetBinContent(bin);
	MCintegral = hPUmc->Integral();
	Datavalue = hPUdata->GetBinContent(bin);
	Dataintegral = hPUdata->Integral();
	if((MCvalue * Dataintegral) != 0) {pu_weight = (Datavalue * MCintegral) / (MCvalue * Dataintegral);}
	else {pu_weight = 1.0;}
	return pu_weight;
}       


bool LowMTRegion_2018::MatchingToGenMuons(TLorentzVector tau, int &genindex){
	double DRact = 0.5;
	genindex=-1;  
	bool ismatched=false;
	for (unsigned int iGen = 0; iGen < tree->mc_px->size(); iGen++){

		TLorentzVector gen_part;
		gen_part.SetPxPyPzE(tree->mc_px->at(iGen),tree->mc_py->at(iGen),tree->mc_pz->at(iGen),tree->mc_energy->at(iGen));
		bool isMuon = abs(tree->mc_pdgId->at(iGen))==13  ? true : false ;    //change names to Muon
		// muon from Z == not checked for now
		if(isMuon) {
			if(tau.DeltaR(gen_part) < DRact) {

				DRact=tau.DeltaR(gen_part);
				genindex=iGen;
				ismatched=true;

			}

		}     
	}             
	return ismatched;
}


int LowMTRegion_2018::GenTaus(){

	int HadronicTau = 0;
	for (unsigned int iGen = 0; iGen < gen_tau_had_pt.size(); iGen++){
		TLorentzVector gen_part;
		gen_part.SetPtEtaPhiE(gen_tau_had_pt.at(iGen), gen_tau_had_eta.at(iGen), gen_tau_had_phi.at(iGen), gen_tau_had_energy.at(iGen));
		if( gen_part.Pt() > 15.) { HadronicTau++; }


	}
	return HadronicTau;
}

bool LowMTRegion_2018::MatchingToGenTaus(TLorentzVector tau, int &genindex){
	double DRact = 0.2;
	genindex=-1;
	bool ismatched=false;
	for (unsigned int iGen = 0; iGen < tree->mc_tau_had_pt->size(); iGen++){
		TLorentzVector gen_part3;
		gen_part3.SetPtEtaPhiE(tree->mc_tau_had_pt->at(iGen), tree->mc_tau_had_eta->at(iGen), tree->mc_tau_had_phi->at(iGen), tree->mc_tau_had_energy->at(iGen));
		if( gen_part3.Pt() > 15.) {
			if(tau.DeltaR(gen_part3) < DRact) {

				DRact=tau.DeltaR(gen_part3);
				genindex=iGen;
				ismatched=true;

			}
		}
	}
	return ismatched;
}


bool LowMTRegion_2018::FillChain(TChain *chain, const char* inputFileList) {

	ifstream infile(inputFileList);
#if 0
	std::string buffer;
#endif
	if (!infile.is_open()) {
		std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << std::endl;
		return kFALSE;
	}  

	std::cout << "TreeUtilities : FillChain " << std::endl;
	char buffer[255];
	while (infile) {
		infile.getline(buffer, 255);                // delim defaults to '\n'
		if (!infile.good()) break;
		std::cout << "Adding " << buffer << " to chain" << std::endl;
		chain->Add(buffer);
	}  
#if 0
	while (1) {
		infile >> buffer;
		if (!infile.good()) break;
		chain->Add(buffer.c_str());
	}  
#endif
	std::cout << "No. of Entries in this tree : " << chain->GetEntries() << std::endl;
	infile.close();
	return kTRUE;
} 

float LowMTRegion_2018::GetEfficiency(float eta, float pt, TH2F *hist) 
{

	double sf = 1.;
	int nbins = hist->GetXaxis()->GetNbins();

	if (eta > (hist->GetXaxis()->GetBinLowEdge(nbins) +  hist->GetXaxis()->GetBinWidth(nbins)))
		eta =  hist->GetXaxis()->GetBinLowEdge(nbins) + (hist->GetXaxis()->GetBinWidth(nbins)/2.0);

	nbins = hist->GetYaxis()->GetNbins();
	if (pt > (hist->GetYaxis()->GetBinLowEdge(nbins) +  hist->GetYaxis()->GetBinWidth(nbins)))
		pt =  hist->GetYaxis()->GetBinLowEdge(nbins) + (hist->GetYaxis()->GetBinWidth(nbins)/2.0);


	Int_t binX = hist->GetXaxis()->FindFixBin(eta);
	Int_t binY = hist->GetYaxis()->FindFixBin(pt);
	if (pt >= 0 && pt <=1000){
		sf = hist->GetBinContent(binX, binY);}

	return sf;
}

double LowMTRegion_2018::GetCollinearMass(const TLorentzVector &tau, const TLorentzVector &mu,  const TLorentzVector MET) {

	double METproj=fabs((MET.Px()*tau.Px()+MET.Py()*tau.Py())/tau.Pt());
	double xth=1;
	if((tau.Pt()+METproj)!=0) xth=tau.Pt()/(tau.Pt()+METproj);
	//now calculate the visibsle mass
	double mass_vis=(tau+mu).M();
	return mass_vis/sqrt(xth);

}


double LowMTRegion_2018::ElePtBinSF_CR5(double elept, double eleSCeta) {
	double SF=1;

	if(fabs(eleSCeta) < 1.4442) { if (elept >= 35. && elept < 131.6) { SF= 0.140 - (0.0029*elept) + (2.56*0.00001* pow(elept,2)) - (8.48 * 0.00000001 * pow(elept,3));} else if (elept >= 131.6 && elept < 359.3) { SF=0.020 - (0.00013*elept) + (3.50*0.0000001* pow(elept,2)) - (2.90 * 0.0000000001 * pow(elept,3));} else { SF=0.00514 + (4.73*0.0000001*elept);}} 

	if(fabs(eleSCeta) > 1.566 && fabs(eleSCeta) < 2.0) { if (elept >= 35. && elept < 125) {SF= 0.1012 - (0.00094*elept) + (3.37*0.000001* pow(elept,2));} else if (elept >= 125.0 && elept < 226.3) {SF= 0.0488 - (11.37*0.00001*elept);} else { SF =0.0241 + (1.24*0.000001*elept);}} 

	if(fabs(eleSCeta) > 2.0 && fabs(eleSCeta) < 2.5) {  if (elept >= 35. && elept < 152) { SF= 0.0622 - (0.00012*elept);} else {SF= 0.0387;}}

	return SF;
}

double LowMTRegion_2018::MutoTauFR(TLorentzVector tau) {
	int Otype=0;
	double DR=0.2;
	double SF_fake =1.;

	// muon faking tau
	for (unsigned int iGen = 0; iGen < tree->mc_px->size(); iGen++){

		TLorentzVector gen_part;
		gen_part.SetPxPyPzE(tree->mc_px->at(iGen),tree->mc_py->at(iGen),tree->mc_pz->at(iGen),tree->mc_energy->at(iGen));


		bool isMuonPrompt = ( gen_part.Pt() >  8. && (abs(tree->mc_pdgId->at(iGen))==13 ) && (tree->mc_status_flags->at(iGen) ==1 )) ? true : false ;

		bool isElectronPrompt = ( gen_part.Pt() >  8. && ( abs(tree->mc_pdgId->at(iGen))==11 ) && ( tree->mc_status_flags->at(iGen) ==1)) ? true : false ;

		bool isElectronTau = ( gen_part.Pt() >  8. && (abs(tree->mc_pdgId->at(iGen))==11 ) && (tree->mc_status_flags->at(iGen) ==6))   ? true : false ;


		bool isMuonTau = ( gen_part.Pt() >  8. && (abs(tree->mc_pdgId->at(iGen))== 13   ) && ( tree->mc_status_flags->at(iGen) == 6))   ? true : false ;
		if(isMuonPrompt || isElectronTau || isElectronPrompt || isMuonTau ) {
			if( tau.DeltaR(gen_part) < DR ) { if(isMuonPrompt) {Otype=1;}
				else if(isElectronPrompt) {Otype=2;} else if( isMuonTau) {Otype=3;} else if(isElectronTau) {Otype=4;} 
			}
		}
	} 

	if(Otype == 1 || Otype == 3) {

		if( fabs(tau.Eta()) > 0 && fabs(tau.Eta()) < 0.4) { SF_fake = 1.06; }
		if( fabs(tau.Eta()) > 0.4 && fabs(tau.Eta()) < 0.8) { SF_fake = 1.02; }
		if( fabs(tau.Eta()) > 0.8 && fabs(tau.Eta()) < 1.2) { SF_fake = 1.10; }
		if( fabs(tau.Eta()) > 1.2 && fabs(tau.Eta()) < 1.7) { SF_fake = 1.03; }
		if( fabs(tau.Eta()) > 1.7 && fabs(tau.Eta()) < 2.3) { SF_fake = 1.94; }



	}  else if ( Otype == 2 || Otype == 4) {

		if( fabs(tau.Eta())  < 1.460) { SF_fake = 1.96; }
		if( fabs(tau.Eta())  > 1.558) { SF_fake = 1.66; }


	} else { SF_fake = 1.;}

	return SF_fake;
}
bool LowMTRegion_2018::matchedToGenObjetcs(TLorentzVector tau, unsigned int &typey) {
	bool ismatched= false;
	double DR=0.2;
	typey = 0;

	for (unsigned int iGen = 0; iGen < tree->mc_px->size(); iGen++){

		TLorentzVector gen_part;
		gen_part.SetPxPyPzE(tree->mc_px->at(iGen),tree->mc_py->at(iGen),tree->mc_pz->at(iGen),tree->mc_energy->at(iGen));


		bool isMuonPrompt = ( gen_part.Pt() >  8. && (abs(tree->mc_pdgId->at(iGen))==13 ) && (tree->mc_status_flags->at(iGen) ==1 )) ? true : false ;

		bool isElectronPrompt = ( gen_part.Pt() >  8. && ( abs(tree->mc_pdgId->at(iGen))==11 ) && ( tree->mc_status_flags->at(iGen) ==1)) ? true : false ;

		bool isElectronTau = ( gen_part.Pt() >  8. && (abs(tree->mc_pdgId->at(iGen))==11 ) && (tree->mc_status_flags->at(iGen) ==6))   ? true : false ;


		bool isMuonTau = ( gen_part.Pt() >  8. && (abs(tree->mc_pdgId->at(iGen))== 13   ) && ( tree->mc_status_flags->at(iGen) == 6))   ? true : false ;

		if(isMuonPrompt || isElectronTau || isElectronPrompt || isMuonTau ) {
			if( tau.DeltaR(gen_part) < DR ) { if(isMuonPrompt) {typey=1;} 
				else if(isElectronPrompt) {typey=2;} else if( isMuonTau) {typey=3;} else if(isElectronTau) {typey=4;} ismatched = true;  break;}


		}
	}
	return ismatched;

}


double LowMTRegion_2018::FakeRate_SSMtLow(double taupt, double jetpt, TH2F *hist) {
	double SF=0.2;
	double reweight = 0;
	int iBin = hist->FindBin(taupt, jetpt);
	SF = hist->GetBinContent(iBin);
	if (SF != 1) reweight = SF/(1-SF);
	return reweight;

}




bool LowMTRegion_2018::TauMatchedToJet(TLorentzVector tau_p4, unsigned int &matched_jet_indx) {

	bool matched_to_reco_jet=false;
	TLorentzVector jet_p4(0.,0.,0.,0.);
	for (unsigned int ijet = 0; ijet < tree->jet_pt->size(); ijet++){
		if(!(fabs(tree->jet_eta->at(ijet)) < 2.3)) continue;
		if(!(tree->jet_isJetID_2018->at(ijet))) continue;
		TLorentzVector jet_p4_tmp;
		jet_p4_tmp.SetPtEtaPhiE(tree->jet_pt->at(ijet), tree->jet_eta->at(ijet), tree->jet_phi->at(ijet), tree->jet_energy->at(ijet));
		if(!(tau_p4.DeltaR(jet_p4_tmp) < 0.2)) continue;
		matched_to_reco_jet=true;
		jet_p4=jet_p4_tmp;
		matched_jet_indx = ijet;
		break;

	}
	return matched_to_reco_jet;
}

DEFINE_FWK_MODULE(LowMTRegion_2018);
