// change done 7 July 2018 -- add >= 1 bjets requirement
#include "LowMTRegion_2017.h"
using namespace std;

LowMTRegion_2017::LowMTRegion_2017(const edm::ParameterSet& iConfig)

{

	isoTau = iConfig.getParameter<string>("isolationTau");
	InputFile = iConfig.getParameter<string>("InputFile");
	DYSample = iConfig.getParameter<string>("DYSample");                                                                                   

	MuonIDIsoSF = iConfig.getParameter<string>("MuonIDIsoSF");
	MuonTriggerSF = iConfig.getParameter<string>("MuonTriggerSF");
	FakeRateSSfile = iConfig.getParameter<string>("FakeRateSSfile");
	FakeRateDYJetsfile = iConfig.getParameter<string>("FakeRateDYJetsfile");
	LooseIsoWPFile  = iConfig.getParameter<string>("LooseIsoWPFile");
	LowPtTES = iConfig.getParameter<string>("LowPtTES");
	HighPtTES = iConfig.getParameter<string>("HighPtTES");
	FESDM = iConfig.getParameter<string>("FESDM");

	size_t npos = -1;
	if(DYSample.find("SE") != npos || DYSample.find("SingleElectron") != npos || DYSample.find("singlephoton") != npos ) { isdata = true; } else { isdata = false;}
	std::cout<<"file:"<<DYSample<<"\tis data \t"<<isdata<<std::endl;
	muscale=1.00;

	TFile *fMuSF = TFile::Open(Form(MuonIDIsoSF.c_str()));
//	eTrigger_barrel = (TGraphAsymmErrors*)(fMuSF->Get("barrel"));
//	eTrigger_endcap = (TGraphAsymmErrors*) (fMuSF->Get("endcap"));
	fhDMuMediumSF = (TH2F*)(fMuSF->Get("SF_TH2F_Barrel"));
        fhDMuMediumSF_endcap = (TH2F*)(fMuSF->Get("SF_TH2F_EndCap"));
	delete fMuSF;

	isOS = iConfig.getParameter<bool>("isOS");
	isSS = iConfig.getParameter<bool>("isSS");

	sftool_loose = new TauIDSFTool(LooseIsoWPFile, 2017, "DeepTau2017v2p1" , "VVVLoose", false);
	sftool_tight = new TauIDSFTool(LooseIsoWPFile,2017,"DeepTau2017v2p1","Tight",false);

	TFile *f_lowptTES = TFile::Open(Form(LowPtTES.c_str()));
	histo_LowPtTES = (TH1F*)(f_lowptTES->Get("tes"));
	delete f_lowptTES;

	TFile *f_highptTES = TFile::Open(Form(HighPtTES.c_str()));
	histo_HighPtTES = (TH1F*)(f_highptTES->Get("tes"));
	delete f_highptTES;

	TFile *f_FESDM =  TFile::Open(Form(FESDM.c_str()));
	graph_FES = (TGraph*)(f_FESDM->Get("fes"));
	delete f_FESDM;        
}


LowMTRegion_2017::~LowMTRegion_2017()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)
	delete sftool_loose;
	delete sftool_tight;
}


//
// member functions
//

// ------------ method called for each event  ------------
	void
LowMTRegion_2017::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
		bool trigger_ele(false), trigger_photon(false);
		h_Events_Before_Skim->Fill(1.);
		// gen filters
		reject_event=false;
		size_t npos = -1;
		Mass=0.;

		if(!isdata &&  (DYSample.find("DYinc") != npos || DYSample.find("DY") != npos ||  DYSample.find("TT") != npos || DYSample.find("WW")  != npos || DYSample.find("TTinc") != npos ||  DYSample.find("WWinc") != npos  )) {


			vector<TLorentzVector> lep;
			lep.clear();
			reject_event=false;

			for (unsigned int iLHE = 0; iLHE < tree->LHE_Pt->size(); ++iLHE) {
				if ( (fabs(tree->LHE_pdgid->at(iLHE)) == 11 || fabs(tree->LHE_pdgid->at(iLHE)) == 13 || fabs(tree->LHE_pdgid->at(iLHE)) == 15)) {
					TLorentzVector l1_p4 ;
					l1_p4.SetPtEtaPhiE(tree->LHE_Pt->at(iLHE),tree->LHE_Eta->at(iLHE),tree->LHE_Phi->at(iLHE),tree->LHE_E->at(iLHE));
					lep.push_back(l1_p4);
				}
			}


			if( lep.size() ==  2 ) { //continue;
				Mass =  (lep.at(0)+ lep.at(1)).M();
				if( DYSample.find("DYinc") != npos || DYSample.find("DYJets_inc")  != npos ) {  if (  (lep.at(0)+ lep.at(1)).M() > 100.) { reject_event = true; Mass=  (lep.at(0)+ lep.at(1)).M();}}

				if( DYSample.find("WWinc") != npos || DYSample.find("WW_inc") != npos ) { if(  (lep.at(0)+ lep.at(1)).M() > 200.) {reject_event = true; Mass=  (lep.at(0)+ lep.at(1)).M(); }}

				if( DYSample.find("TTinc") != npos || DYSample.find("TTJets_inc") != npos || DYSample.find("TTTo2L2Nu_TuneCUETP8M2_ttHtranche3") != npos) { if(  (lep.at(0)+ lep.at(1)).M() > 500.) {reject_event = true; Mass=  (lep.at(0)+ lep.at(1)).M(); }}
			} else { reject_event = false ; }
		}


		if(!isdata &&  (DYSample.find("WJets") != npos || DYSample.find("WJetsinc") != npos )) {

			reject_event=false;

			TLorentzVector l_p4, nu_p4, lnu_p4;
			l_p4.SetPxPyPzE(0, 0, 0, 0);
			nu_p4.SetPxPyPzE(0, 0, 0, 0);
			lnu_p4.SetPxPyPzE(0, 0, 0, 0);
			int l_pdgid = 0, nu_pdgid = 0;
			Mass = 0;
			for (unsigned int iLHE = 0; iLHE < tree->LHE_Pt->size(); ++iLHE) {

				if (  (fabs(tree->LHE_pdgid->at(iLHE)) == 11 || fabs(tree->LHE_pdgid->at(iLHE)) == 13 || fabs(tree->LHE_pdgid->at(iLHE)) == 15))  {
					l_p4.SetPtEtaPhiE(tree->LHE_Pt->at(iLHE),tree->LHE_Eta->at(iLHE),tree->LHE_Phi->at(iLHE),tree->LHE_E->at(iLHE));
					l_pdgid = tree->LHE_pdgid->at(iLHE);
				}
				else if (abs(tree->LHE_pdgid->at(iLHE)) == 12 || abs(tree->LHE_pdgid->at(iLHE)) == 14 || abs(tree->LHE_pdgid->at(iLHE)) == 16) {
					nu_p4.SetPtEtaPhiE(tree->LHE_Pt->at(iLHE),tree->LHE_Eta->at(iLHE),tree->LHE_Phi->at(iLHE),tree->LHE_E->at(iLHE));
					nu_pdgid = tree->LHE_pdgid->at(iLHE);
				}
			}


			if(fabs(l_pdgid)+1 == fabs(nu_pdgid)) {
				lnu_p4 = l_p4 + nu_p4;
				Mass=  lnu_p4.Pt();

				if(DYSample.find("WJetsinc") != npos  || DYSample.find("WJetsToLNu_TuneCUETP8M1") != npos) {if (lnu_p4.Pt() > 100) { reject_event = true; Mass=  lnu_p4.Pt();}}

			}  else { reject_event = false;}
		}

		if(reject_event) continue;

		if(!isdata) h_Fill_Mass_Gen_toChk->Fill(Mass,tree->mc_w_sign);

		if(DYSample.find("singlephoton") != npos && isdata ) {

			if(!(!tree->trig_HLT_Ele35_WPTight_Gsf_accept && (!tree->trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_accept) &&  (tree->trig_HLT_Photon200_accept))) continue;
		}

		if((DYSample.find("SE") != npos || DYSample.find("SingleElectron") != npos ) && isdata ){
			if( !( tree->trig_HLT_Ele35_WPTight_Gsf_accept || tree->trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_accept )) continue;
		}

		if(!isdata) {
			if( !( tree->trig_HLT_Ele35_WPTight_Gsf_accept || tree->trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_accept || tree->trig_HLT_Photon200_accept )) continue;
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
		if(isdata) if(!tree->trig_Flag_eeBadScFilter_accept) continue;

		if(isdata) h_Events_After_MetFilters->Fill(1.);
		if(!isdata) h_Events_After_MetFilters->Fill(1., tree->mc_w_sign);


		// loop for OS and SS events
		// events filters

		sumweight =0;
		weighthis=1.;
		//std::cout<<"DYSample:"<< DYSample << std::endl;
		if(!isdata) {pu_weight =PU_2017_Rereco::MC_pileup_weight(tree->mc_trueNumInteractions, DYSample, "Rereco_all"); } else {pu_weight = 1;}

		sumweight = sumweight+tree->mc_weight;
		// prefiring weight 
		double pref_wt = 1.;
		if(!isdata) pref_wt =  tree->ev_prefiringweight;


		event_wt_NNPDF = 1.;
		if(!isdata &&  ( DYSample.find("TT") != npos ||  DYSample.find("TTinc") != npos || DYSample.find("TT_had") != npos || DYSample.find("TT_2l2nu") != npos || DYSample.find("TT_semilep") != npos || DYSample.find("DY") != npos )) {

			TLorentzVector ele_gen, tau_gen;
			int use_ele_status = 1;
			for(unsigned int iMC=0 ; iMC<tree->mc_n ; ++iMC){
				if( abs(tree->mc_pdgId->at(iMC)) == 11 && abs(tree->mc_status->at(iMC)) == 23) {use_ele_status=23; break; }

			}
			int ele_tau_counter = 0;
			bool ele_found=false;
			bool tau_found=false;
			for(unsigned int iMC=0 ; iMC<tree->mc_n ; ++iMC){
				TLorentzVector genP;
				genP.SetPtEtaPhiE( tree->mc_pt->at(iMC), tree->mc_eta->at(iMC), tree->mc_phi->at(iMC), tree->mc_energy->at(iMC));
				if( abs(tree->mc_pdgId->at(iMC)) == 11 && abs(tree->mc_status->at(iMC)) == use_ele_status) { ele_gen = genP; ele_found=true; }


				if( abs(tree->mc_pdgId->at(iMC)) == 15 && abs(tree->mc_status->at(iMC)) == 2) { tau_gen = genP; tau_found=true;}
				ele_tau_counter++ ;
			}
			if ( ele_found && tau_found && ele_tau_counter >= 2) {
				double Gen_Z_mass=(ele_gen+tau_gen).M();

				double Gen_Led_Et;
				int region=0;

				if(fabs(ele_gen.Eta())<1.4442 && fabs(tau_gen.Eta())<1.46)region=1;
				else if ((fabs(tau_gen.Eta())>1.46 && fabs(ele_gen.Eta())<1.4442) || (fabs(tau_gen.Eta())<1.46 && fabs(ele_gen.Eta())> 1.566))region=2;
				else if(fabs(ele_gen.Eta())>1.566 && fabs(tau_gen.Eta())>1.46)region=3;

				if(ele_gen.Pt() > tau_gen.Pt() ) { Gen_Led_Et = ele_gen.Pt() ;} else {Gen_Led_Et = tau_gen.Pt();}
				if(DYSample.find("DY") != npos ) {
					if(Gen_Z_mass<120){
						if(region==1)      event_wt_NNPDF=(Gen_Led_Et<150) ? 3.596-0.2076 *Gen_Led_Et+0.005795*pow(Gen_Led_Et,2)-7.421e-05*pow(Gen_Led_Et,3)+4.447e-07*pow(Gen_Led_Et,4)-1.008e-9 *pow(Gen_Led_Et,5) : 0.969125;
						else if(region==2) event_wt_NNPDF=(Gen_Led_Et<150) ? 2.066-0.09495*Gen_Led_Et+0.002664*pow(Gen_Led_Et,2)-3.242e-05*pow(Gen_Led_Et,3)+1.755e-07*pow(Gen_Led_Et,4)-3.424e-10*pow(Gen_Led_Et,5) : 1.191875;
						else if(region==3) event_wt_NNPDF=(Gen_Led_Et<150) ? 2.865-0.1443 *Gen_Led_Et+0.003799*pow(Gen_Led_Et,2)-4.482e-05*pow(Gen_Led_Et,3)+2.429e-07*pow(Gen_Led_Et,4)-4.93e-10 *pow(Gen_Led_Et,5) : 0.9609375;
					}
					else{
						if(region==1)      event_wt_NNPDF=(Gen_Z_mass<5000) ? 0.8934+0.0002193 *Gen_Z_mass-1.961e-7*pow(Gen_Z_mass,2)+8.704e-11*pow(Gen_Z_mass,3)-1.551e-14*pow(Gen_Z_mass,4)+1.112e-18*pow(Gen_Z_mass,5) : 1.74865;
						else if(region==2) event_wt_NNPDF=(Gen_Z_mass<5000) ? 0.8989+0.000182  *Gen_Z_mass-1.839e-7*pow(Gen_Z_mass,2)+1.026e-10*pow(Gen_Z_mass,3)-2.361e-14*pow(Gen_Z_mass,4)+1.927e-18*pow(Gen_Z_mass,5) : 1.302025;
						else if(region==3) event_wt_NNPDF=(Gen_Z_mass<5000) ? 0.9328-7.23e-6   *Gen_Z_mass+3.904e-9*pow(Gen_Z_mass,2)+2.454e-11*pow(Gen_Z_mass,3)-1.038e-14*pow(Gen_Z_mass,4)+1.543e-18*pow(Gen_Z_mass,5) : 2.396125;
					}

				} else {
					event_wt_NNPDF = GetTTWeight(Gen_Z_mass) ;

				}
			}

		}

		if(!isdata) weighthis = tree->mc_w_sign*pu_weight*pref_wt*event_wt_NNPDF; 
		if(isdata) { weighthis = 1.;}

		if(!isdata &&  ( DYSample.find("TT") != npos || DYSample.find("TTinc") != npos )) {

			double weight_top1 = 1.0;
			//double weight_top2 = 1.0;
			int tops=0;
			double tmp_t1 , tmp_t2;
			tmp_t1=0;
			tmp_t2=0;
			for (unsigned int iLHE = 0; iLHE < tree->LHE_Pt->size(); ++iLHE) {


				if ( fabs(tree->LHE_pdgid->at(iLHE)) != 6) continue;

				bool top1_found=false;
				bool top2_found=false;

				if ( tree->LHE_pdgid->at(iLHE) == 6) { tmp_t1 = exp(0.0615-0.0005*tree->LHE_Pt->at(iLHE)); top1_found=true; }

				if ( tree->LHE_pdgid->at(iLHE) == -6) { tmp_t2 = exp(0.0615-0.0005*tree->LHE_Pt->at(iLHE)); top2_found=true; }

				if(top1_found || top2_found) tops++;
			}
			if(tops == 2) weight_top1 = sqrt(tmp_t1 * tmp_t2);


			weighthis = weighthis* weight_top1;
		}


		//// about NNPDF reweighting


		int electrons_no = 0;
		for (unsigned int dau1index = 0; dau1index < tree->gsf_caloEnergy->size(); dau1index++){
			float ET1 = tree->gsf_caloEnergy->at(dau1index)*sin(2.*atan(exp(-1.*tree->gsf_eta->at(dau1index)))) ;


			TLorentzVector DauEle ;
			DauEle.SetPtEtaPhiM(tree->gsf_pt->at(dau1index), tree->gsf_eta->at(dau1index), tree->gsf_phi->at(dau1index),m_el);
			if(DauEle.Pt() > 10. && fabs(DauEle.Eta() ) < 2.5 && tree->gsf_VID_cutBasedElectronID_Fall17_94X_V1_veto->at(dau1index) ) {

				for (unsigned int dau2index = 0; dau2index < tree->gsf_caloEnergy->size(); dau2index++){
					if(dau1index == dau2index) continue;
					TLorentzVector DauEle2 ;
					DauEle2.SetPtEtaPhiM(tree->gsf_pt->at(dau2index), tree->gsf_eta->at(dau2index), tree->gsf_phi->at(dau2index),m_el);
					if(DauEle.DeltaR(DauEle2) < 0.5) continue;


					if(DauEle2.Pt() > 10. && fabs(DauEle2.Eta() ) < 2.5 && tree->gsf_VID_cutBasedElectronID_Fall17_94X_V1_veto->at(dau2index) ) {

						electrons_no++;

					}
				}
			}
		}
		if(electrons_no > 0) continue;




		vector<unsigned int> ele_indexes, tau_indexes;
		ele_indexes.clear();
		tau_indexes.clear();

		for (unsigned int dau1index = 0; dau1index < tree->gsf_caloEnergy->size(); dau1index++){
			float ET1 = tree->gsf_caloEnergy->at(dau1index)*sin(2.*atan(exp(-1.*tree->gsf_eta->at(dau1index)))) ;


			TLorentzVector DauEle ;
			DauEle.SetPtEtaPhiM(tree->gsf_pt->at(dau1index), tree->gsf_eta->at(dau1index), tree->gsf_phi->at(dau1index),m_el);
			if(!( DauEle.Pt() > 50.  && fabs(tree->gsf_sc_eta->at(dau1index)) < 2.5) ) continue;
			if(!(tree->gsf_VID_heepElectronID_HEEPV70->at(dau1index))) continue;

			ele_indexes.push_back(dau1index);

		}

		if(ele_indexes.size() == 0 ) continue;

		if(isdata) h_Fill_ElePass->Fill(1.);
		if(!isdata) h_Fill_ElePass->Fill(1., tree->mc_w_sign);

		for (unsigned int dau2index = 0; dau2index < tree->tau_pt->size(); dau2index++){

			if(!(tree->tau_pt->at(dau2index) > 30.)) continue;
			if(!(fabs(tree->tau_eta->at(dau2index)) < 2.3)) continue;

			if(!(tree->tau_decayModeFindingNewDMs->at(dau2index) > 0.5)) continue;
			if(tree->tau_decayMode->at(dau2index) == 5 || tree->tau_decayMode->at(dau2index) == 6 ) continue;

			if(!(tree->tau_byVVVLooseDeepTau2017v2p1VSjet->at(dau2index) > 0.5)) continue;
			if(!(tree->tau_byLooseDeepTau2017v2p1VSe->at(dau2index) > 0.5)) continue;
			if(!(tree->tau_byTightDeepTau2017v2p1VSmu->at(dau2index) > 0.5)) continue;



			TLorentzVector DauTau ;
			DauTau.SetPxPyPzE(tree->tau_px->at(dau2index), tree->tau_py->at(dau2index), tree->tau_pz->at(dau2index),tree->tau_energy->at(dau2index));

			tau_indexes.push_back(dau2index);
		}

		if(tau_indexes.size() == 0 ) continue;

		if(isdata) h_Fill_TauPass->Fill(1.);
		if(!isdata) h_Fill_TauPass->Fill(1., tree->mc_w_sign);

		int muon_counter = 0;
		for(unsigned int imuon=0; imuon < tree->mu_ibt_pt->size(); imuon++) {
			if(!(tree->mu_ibt_pt->at(imuon) > 35.)) continue;
			if(!(fabs(tree->mu_ibt_eta->at(imuon)) < 2.4)) continue;
			if(!(tree->mu_isHighPtMuon->at(imuon) == 1)) continue;//mu_isHighPtMuon->at(imuon))) continue;
			if(!(tree->mu_isoTrackerBased03->at(imuon) < 0.15)) continue;
			muon_counter++;

		}
		if(muon_counter != 0 ) continue;


		if(isdata) h_Fill_MuonExtra->Fill(1.);
		if(!isdata) h_Fill_MuonExtra->Fill(1., tree->mc_w_sign);

		double masspair = 0;
		TLorentzVector FirstObj(0.,0.,0.,0), SecondObj(0.,0.,0.,0);
		unsigned int first_index = -1;
		unsigned int second_index = -1;

		TLorentzVector FirstObj_FF(0.,0.,0.,0), SecondObj_FF(0.,0.,0.,0);
		first_index_FF= -1;
		second_index_FF = -1;
		double masspair_FF= 0;
		bool found_one_ele = false;

		for (unsigned int dau1index = 0; dau1index < tree->gsf_caloEnergy->size(); dau1index++){
			if(found_one_ele) continue;

			float ET1 = tree->gsf_caloEnergy->at(dau1index)*sin(2.*atan(exp(-1.*tree->gsf_eta->at(dau1index)))) ;

			TLorentzVector DauEle ;
			DauEle.SetPtEtaPhiM(tree->gsf_pt->at(dau1index), tree->gsf_eta->at(dau1index), tree->gsf_phi->at(dau1index),m_el);
			if(!( DauEle.Pt() > 50.  && fabs(tree->gsf_sc_eta->at(dau1index)) < 2.5) ) continue;
			if(!(tree->gsf_VID_heepElectronID_HEEPV70->at(dau1index))) continue;



			for (unsigned int dau2index = 0; dau2index < tree->tau_pt->size(); dau2index++){

				if(!(tree->tau_pt->at(dau2index) > 30.)) continue;
				if(!(fabs(tree->tau_eta->at(dau2index)) < 2.3)) continue;


				if(!(tree->tau_decayModeFindingNewDMs->at(dau2index) > 0.5)) continue;
				if(tree->tau_decayMode->at(dau2index) == 5 || tree->tau_decayMode->at(dau2index) == 6 ) continue;

				if(!(tree->tau_byVVVLooseDeepTau2017v2p1VSjet->at(dau2index) > 0.5)) continue;
				if(!(tree->tau_byLooseDeepTau2017v2p1VSe->at(dau2index) > 0.5)) continue;
				if(!(tree->tau_byTightDeepTau2017v2p1VSmu->at(dau2index) > 0.5)) continue;




				TLorentzVector DauTau ;
				DauTau.SetPxPyPzE(tree->tau_px->at(dau2index), tree->tau_py->at(dau2index), tree->tau_pz->at(dau2index),tree->tau_energy->at(dau2index));

				TLorentzVector DauTau_raw ;
				DauTau_raw = DauTau;

				if(!isdata) {
					if(tree->tau_decayMode->at(dau2index) == 0 || tree->tau_decayMode->at(dau2index) == 1) {

						DauTau =  ElectronTauFakeEnergyScale(DauTau, tree->tau_decayMode->at(dau2index), "nom");
					}
				}

				if(!isdata) DauTau = TauEnergyScale(DauTau, tree->tau_decayMode->at(dau2index), "normal");

				TVector2 METcorr;
				METcorr.SetMagPhi(tree->MET_FinalCollection_Pt,tree->MET_FinalCollection_phi);

				TLorentzVector metV;
				metV.SetPxPyPzE(METcorr.Px(), METcorr.Py(), 0. , tree->MET_FinalCollection_Pt);


				TLorentzVector met_corrected;
				met_corrected = metV -  DauTau + DauTau_raw;

				double mt_cut = mTCalculation(met_corrected.Px(), met_corrected.Py(), DauEle.Px(), DauEle.Py(), DauEle.Pt());
				if(mt_cut > 120.) continue;



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
					if(!(tree->jet_isJetID_2017->at(ijet))) continue;
					TLorentzVector jet_p4_tmp;
					jet_p4_tmp.SetPtEtaPhiE(tree->jet_pt->at(ijet), tree->jet_eta->at(ijet), tree->jet_phi->at(ijet), tree->jet_energy->at(ijet));
					if(!(DauTau.DeltaR(jet_p4_tmp) < 0.2)) continue;
					matched_to_reco_jet=true;
					jet_p4=jet_p4_tmp;
					break;

				}

				if(!(matched_to_reco_jet)) continue;
				int j_dm = -1, k_eta = -1, k_pt = -1;
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
//				else if (fabs(tree->tau_eta->at(dau2index)) > 1.558) {
//					k_eta = 1;
//				}
				else {
//					continue;
					k_eta = 1;
				}

				double weight_sd = 1.;
				double fake_wt = 1.;
				// HEEP SF
//				if( (!isdata) )  {  if( tree->gsf_isEB->at(dau1index)) {weight_sd = 0.967 * eTrigger_barrel->Eval(DauEle.Pt()); } else { weight_sd =  0.973 * eTrigger_endcap->Eval(DauEle.Pt()); }}
//
				if( (!isdata) )  {  if( tree->gsf_isEB->at(dau1index)) {weight_sd = 0.967 * GetEfficiency(tree->gsf_sc_eta->at(dau1index), DauEle.Pt()) ;} else { weight_sd =  0.973 * GetEfficiency(tree->gsf_sc_eta->at(dau1index), DauEle.Pt()) ;}}
 

				//0.89 is tauID SF
				if(!isdata) { fake_wt =  MutoTauFR(DauTau,"normal","");}

				double final_weight = weight_sd*fake_wt*weighthis;
				found_one_ele = true;

				if(tree->gsf_charge->at(dau1index) * tree->tau_charge->at(dau2index) == -1 ) {
					if (tree->tau_byTightDeepTau2017v2p1VSjet->at(dau2index) > 0.5) {
						int hadtau=0;
						if(!isdata){
							int genindex_y;
							genindex_y = -1;
							if(MatchingToGenTaus(DauTau, genindex_y)  && genindex_y!=-1) {hadtau=5;}
						}
						hh[0][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index), jet_p4.Pt(), final_weight*sftool_tight->getSFvsPT(tree->tau_pt->at(dau2index),hadtau));
						h1[0][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index), final_weight*sftool_tight->getSFvsPT(tree->tau_pt->at(dau2index),hadtau));
						h2[0][j_dm][k_eta]->Fill(jet_p4.Pt(), final_weight*sftool_tight->getSFvsPT(tree->tau_pt->at(dau2index)));
						h3[0][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index)/jet_p4.Pt(), final_weight*sftool_tight->getSFvsPT(tree->tau_pt->at(dau2index),hadtau));

						if(tree->tau_pt->at(dau2index) < 150.) { h4[0][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index), tree->tau_pt->at(dau2index)/jet_p4.Pt(), final_weight*sftool_tight->getSFvsPT(tree->tau_pt->at(dau2index),hadtau)); }
						if(tree->tau_pt->at(dau2index)  >= 150. ) { h5[0][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index), tree->tau_pt->at(dau2index)/jet_p4.Pt(), final_weight*sftool_tight->getSFvsPT(tree->tau_pt->at(dau2index),hadtau)); }



					}
					if ((tree->tau_byTightDeepTau2017v2p1VSjet->at(dau2index) < 0.5) && (tree->tau_byVVVLooseDeepTau2017v2p1VSjet->at(dau2index) > 0.5)) {

						int hadtau=0;
						if(!isdata){
							int genindex_y;
							genindex_y = -1;
							if(MatchingToGenTaus(DauTau, genindex_y)  && genindex_y!=-1) {hadtau=5;}
						}

						hh[1][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index), jet_p4.Pt(), final_weight*sftool_loose->getSFvsPT(tree->tau_pt->at(dau2index),hadtau));
						h1[1][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index), final_weight*sftool_loose->getSFvsPT(tree->tau_pt->at(dau2index),hadtau));
						h2[1][j_dm][k_eta]->Fill(jet_p4.Pt(), final_weight*sftool_loose->getSFvsPT(tree->tau_pt->at(dau2index),hadtau));
						h3[1][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index)/jet_p4.Pt(), final_weight*sftool_loose->getSFvsPT(tree->tau_pt->at(dau2index),hadtau));

						if(tree->tau_pt->at(dau2index) < 150. ) { h4[1][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index), tree->tau_pt->at(dau2index)/jet_p4.Pt(), final_weight*sftool_loose->getSFvsPT(tree->tau_pt->at(dau2index),hadtau)); }
						if(tree->tau_pt->at(dau2index)  >= 150. ) { h5[1][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index), tree->tau_pt->at(dau2index)/jet_p4.Pt(), final_weight*sftool_loose->getSFvsPT(tree->tau_pt->at(dau2index),hadtau)); }



					}

				} else if (tree->gsf_charge->at(dau1index) * tree->tau_charge->at(dau2index) == 1) {


					if (tree->tau_byTightDeepTau2017v2p1VSjet->at(dau2index) > 0.5) {

						int hadtau_yx=0;
						if(!isdata){
							int genindex_yx;
							genindex_yx = -1;
							if(MatchingToGenTaus(DauTau, genindex_yx)  && genindex_yx!=-1) {hadtau_yx=5;}
						}
						hh[2][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index), jet_p4.Pt(), final_weight*sftool_tight->getSFvsPT(tree->tau_pt->at(dau2index),hadtau_yx));
						h1[2][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index), final_weight*sftool_tight->getSFvsPT(tree->tau_pt->at(dau2index),hadtau_yx));
						h2[2][j_dm][k_eta]->Fill(jet_p4.Pt(), final_weight*sftool_tight->getSFvsPT(tree->tau_pt->at(dau2index),hadtau_yx));
						h3[2][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index)/jet_p4.Pt(), final_weight*sftool_tight->getSFvsPT(tree->tau_pt->at(dau2index),hadtau_yx));

						if(tree->tau_pt->at(dau2index)   < 150.) { h4[2][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index), tree->tau_pt->at(dau2index)/jet_p4.Pt(), final_weight*sftool_tight->getSFvsPT(tree->tau_pt->at(dau2index),hadtau_yx)); }
						if(tree->tau_pt->at(dau2index)  >= 150. ) { h5[2][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index), tree->tau_pt->at(dau2index)/jet_p4.Pt(), final_weight*sftool_tight->getSFvsPT(tree->tau_pt->at(dau2index),hadtau_yx)); }



					}

					if ((tree->tau_byTightDeepTau2017v2p1VSjet->at(dau2index) < 0.5) && (tree->tau_byVVVLooseDeepTau2017v2p1VSjet->at(dau2index) > 0.5)) {

						int hadtau_yx=0;
						if(!isdata){
							int genindex_yx;
							genindex_yx = -1;
							if(MatchingToGenTaus(DauTau, genindex_yx)  && genindex_yx!=-1) {hadtau_yx=5;}
						}

						hh[3][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index), jet_p4.Pt(), final_weight*sftool_loose->getSFvsPT(tree->tau_pt->at(dau2index),hadtau_yx));
						h1[3][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index), final_weight*sftool_loose->getSFvsPT(tree->tau_pt->at(dau2index),hadtau_yx));
						h2[3][j_dm][k_eta]->Fill(jet_p4.Pt(), final_weight*sftool_loose->getSFvsPT(tree->tau_pt->at(dau2index),hadtau_yx));
						h3[3][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index)/jet_p4.Pt(), final_weight*sftool_loose->getSFvsPT(tree->tau_pt->at(dau2index),hadtau_yx));

						if(tree->tau_pt->at(dau2index)  < 150.) { h4[3][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index), tree->tau_pt->at(dau2index)/jet_p4.Pt(), final_weight*sftool_loose->getSFvsPT(tree->tau_pt->at(dau2index),hadtau_yx)); }
						if(tree->tau_pt->at(dau2index)  >= 150. ) { h5[3][j_dm][k_eta]->Fill(tree->tau_pt->at(dau2index), tree->tau_pt->at(dau2index)/jet_p4.Pt(), final_weight*sftool_loose->getSFvsPT(tree->tau_pt->at(dau2index),hadtau_yx)); }



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
LowMTRegion_2017::beginJob()
{
	//	fIn = TFile::Open(InputFile.c_str());
	std::cout<<"InputFile::"<<InputFile<< std::endl;
	//	file_db1.open(InputFile);                                                                                       
	float xbins_left[] = {0, 30, 40, 50, 80, 150};
	float ybins_left[]= {0, 0.5, 0.6, 0.65, 0.7, 0.75, 1., 5.};
	float xbins_right[] = {150, 1000};
	float ybins_right[] = {0, 0.7, 1., 5.};


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


				hh[i][k][l] =  fs->make<TH2F>(nname, nname, 1000, 0, 1000, 1000, 0, 1000) ;
				h1[i][k][l] = fs->make<TH1F>(nname1, nname1, 1000, 0, 1000) ;
				h2[i][k][l] = fs->make<TH1F>(nname2, nname2, 1000, 0, 1000) ;
				h3[i][k][l] = fs->make<TH1F>(nname3, nname3, 1000, 0, 10) ;

				//				h4[i][k][l] = fs->make<TH2F>(nname4, nname4, 1000, 0, 1000, 100, 0, 5); // nxaxis_left-1, xbins_left,nyaxis_left-1, ybins_left) ;
				//				h5[i][k][l] = fs->make<TH2F>(nname5, nname5, 1000, 0, 1000, 100, 0, 5); //nxaxis_right-1, xbins_right, nyaxis_right-1, ybins_right) ;

				h4[i][k][l] = fs->make<TH2F>(nname4, nname4, nxaxis_left-1, xbins_left,nyaxis_left-1, ybins_left) ;
				h5[i][k][l] = fs->make<TH2F>(nname5, nname5, nxaxis_right-1, xbins_right, nyaxis_right-1, ybins_right) ;

				hh[i][k][l]->Sumw2();
				h1[i][k][l]->Sumw2();
				h2[i][k][l]->Sumw2();
				h3[i][k][l]->Sumw2();
				h4[i][k][l]->Sumw2();
				h5[i][k][l]->Sumw2();

			}
		}
	}


}


// ------------ method called once each job just after ending the event loop  ------------
	void 
LowMTRegion_2017::endJob() 
{

	std::cout<<"Total No of nEventsRaw"<<":\t"<<nEventsRaw<<"\t: Stored events : \t"<<nEventsStored<<"\t weighted events :\t"<<mc_nEventsWeighted<<"\t sum of weights : \t"<<nEventsiihe<<std::endl;
//	delete eTrigger_barrel ;
//	delete eTrigger_endcap;
        delete fhDMuMediumSF;
        delete fhDMuMediumSF_endcap;
	delete histo_LowPtTES; delete histo_HighPtTES;
	delete                 graph_FES ;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
LowMTRegion_2017::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}


TLorentzVector LowMTRegion_2017::ElectronTauFakeEnergyScale(TLorentzVector tau, int taudecaymode, string type_of){

	double DR=0.2;
	TLorentzVector newtau;
	newtau= tau;
	for (unsigned int iGen = 0; iGen < tree->mc_px->size(); iGen++){

		TLorentzVector gen_part;
		gen_part.SetPxPyPzE(tree->mc_px->at(iGen),tree->mc_py->at(iGen),tree->mc_pz->at(iGen),tree->mc_energy->at(iGen));



		bool isElectronPrompt = ( gen_part.Pt() >  8. && ( abs(tree->mc_pdgId->at(iGen))==11 ) && (tree->mc_status_flags->at(iGen) >> 0 & 1  )) ? true : false ;

		bool isElectronTau = ( gen_part.Pt() >  8. && (abs(tree->mc_pdgId->at(iGen))==11 ) && ( tree->mc_status_flags->at(iGen) >> 5 & 1 ))  ? true : false ;


		if( (isElectronTau || isElectronPrompt) && ( tau.DeltaR(gen_part) < DR )) {
			if(taudecaymode == 0 || taudecaymode == 1) {

				newtau = tau*getEFakeTESSF(tau.Eta(), taudecaymode, type_of) ;


			}
		}
	}
	return newtau;

}



double LowMTRegion_2017::getEFakeTESSF(double taueta, int taudecaymode, string type_of){
	double sf_etauFS = 1;
	int i=-1;
	if(taudecaymode == 0 && fabs(taueta) < 1.5) { i =0;}
	else if(taudecaymode == 1 && fabs(taueta) < 1.5) { i =1;}
	else if(taudecaymode == 0 && (!(fabs(taueta) < 1.5))) { i =2;}
	else if(taudecaymode == 1 && (!(fabs(taueta) < 1.5))) { i =3;}

	double value_norm = graph_FES->GetY()[i];
	double value_high = graph_FES->GetErrorYhigh(i);
	double value_low = graph_FES->GetErrorYlow(i);
	if(i !=-1) {
		if(type_of == "nom") sf_etauFS = value_norm;
		else if (type_of == "Up") sf_etauFS = value_norm+value_high;
		else if (type_of == "Down") sf_etauFS =value_norm-value_low;

		return sf_etauFS;
	} else return 1;
}

double LowMTRegion_2017::GetTopQscale(TString process, float Mll, TString var){
	if (Mll <0 || var=="nom") {
		return 1;
	}
	else {
		double val_unc = 0.135 - 5.981*pow(10,-5)*Mll + 1.807*pow(10,-7)*pow(Mll,2) - 1.815*pow(10,-10)*pow(Mll,3) + 7.875*pow(10,-14)*pow(Mll,4) - 1.229*pow(10,-17)*pow(Mll,5);

		double weight = 0;
		if (var=="nom") weight = 1;
		else if (var=="up") weight = 1+fabs(val_unc);
		else if (var=="down") weight = 1-fabs(val_unc);

		return weight;
	}
}


double LowMTRegion_2017::GetTopPDF(TString process, float Mll, TString var){
	if (Mll <0 || var=="nom") {
		return 1;
	}
	else {
		double val_unc = 0.49 - 0.0007795*Mll + 1.59*pow(10,-6)*pow(Mll,2) - 1.166*pow(10,-9)*pow(Mll,3) + 3.93*pow(10,-13)*pow(Mll,4) - 4.72*pow(10,-17)*pow(Mll,5);

		double weight = 0;
		if (var=="nom") weight = 1;
		else if (var=="up") weight = 1+fabs(val_unc);
		else if (var=="down") weight = 1-fabs(val_unc);

		return weight;
	}
}

double LowMTRegion_2017::GetWWNLO(TString process, float Mll, TString var){
	if (Mll <0 || var=="nom") {
		return 1;
	}
	else {

		double val_unc = 1.0 - ( 0.993 - 2.001*pow(10,-4)*Mll + 2.838*pow(10,-8)*pow(Mll,2));

		double weight = 0;
		if (var=="nom") weight = 1;
		else if (var=="up") weight = 1+fabs(val_unc);
		else if (var=="down") weight = 1-fabs(val_unc);

		return weight;


	}

}




bool LowMTRegion_2017::OverLap05(TLorentzVector l1 , TLorentzVector l2, float conesize) {
	if(dR(l1.Eta(), l1.Phi(), l2.Eta(), l2.Phi()) <= conesize) return true;
	else return false;
}  

float LowMTRegion_2017::deltaPhi( float a, float b) {
	float result = a-b;
	while (result > M_PI) result -= 2* M_PI;
	while (result <= -M_PI) result += 2* M_PI;
	return (fabs(result));

} 

float LowMTRegion_2017::dR(float l1eta, float l1phi, float l2eta, float l2phi ) {
	float deta = l1eta - l2eta;
	float dphi = deltaPhi(l1phi,l2phi);
	return sqrt(deta*deta + dphi*dphi);
}




float LowMTRegion_2017::mTCalculation(float metx, float mety, float mupx, float mupy, float mupt){
	float mt = -1;
	float pX = mupx+metx;
	float pY = mupy+mety;
	float et = mupt + TMath::Sqrt(metx*metx + mety*mety);
	mt = TMath::Sqrt(et*et-(pX*pX + pY*pY));
	return mt;

}     

TLorentzVector LowMTRegion_2017::TauEnergyScale(TLorentzVector tau, int taudecaymode, string type_of){

	double DRact = 0.2;
	TLorentzVector newtau;
	newtau= tau;
	bool ismatched=false;
	for (unsigned int iGen = 0; iGen < tree->mc_tau_had_pt->size(); iGen++){
		TLorentzVector gen_part3;
		gen_part3.SetPtEtaPhiE(tree->mc_tau_had_pt->at(iGen), tree->mc_tau_had_eta->at(iGen), tree->mc_tau_had_phi->at(iGen), tree->mc_tau_had_energy->at(iGen));
		if( gen_part3.Pt() > 15.) {
			if(tau.DeltaR(gen_part3) < DRact) {

				DRact=tau.DeltaR(gen_part3);
				ismatched=true;

			}
		}
	}


	if(ismatched) {
		double taupt_nom = tau.Pt();
		if(type_of == "normal") { newtau = tau*getTESSF(taupt_nom, taudecaymode, "normal");}
		else if ( type_of == "upDM0" && taudecaymode == 0) { newtau = tau*getTESSF(taupt_nom, taudecaymode, "Up");}
		else if ( type_of == "upDM1" && taudecaymode == 1) { newtau = tau*getTESSF(taupt_nom, taudecaymode, "Up");}
		else if ( type_of == "upDM10" && taudecaymode == 10) { newtau = tau*getTESSF(taupt_nom, taudecaymode, "Up");}
		else if ( type_of == "upDM11" && taudecaymode == 11) { newtau = tau*getTESSF(taupt_nom, taudecaymode, "Up");}
		else if ( type_of == "downDM0" && taudecaymode == 0) { newtau = tau*getTESSF(taupt_nom, taudecaymode, "Down");}
		else if ( type_of == "downDM1" && taudecaymode == 1) { newtau = tau*getTESSF(taupt_nom, taudecaymode, "Down");}
		else if ( type_of == "downDM10" && taudecaymode == 10) { newtau = tau*getTESSF(taupt_nom, taudecaymode, "Down");}
		else if ( type_of == "downDM11" && taudecaymode == 11) { newtau = tau*getTESSF(taupt_nom, taudecaymode, "Down");}
		else if ( type_of == "FullUp" ) {

			newtau = tau * (1. + 0.03);
		} else if ( type_of == "FullDown" ) {

			newtau = tau * (1. - 0.03);
		}
	}
	return newtau;
}

double LowMTRegion_2017::getTESSF(double taupt, int taudecaymode, string type_of){

	double low_pt_tau_thr = 34.;
	double high_pt_tau_thr = 170.;
	double sf_tes = 1.;
	if(type_of != "None") {
		int bin = histo_LowPtTES->GetXaxis()->FindBin(taudecaymode);
		sf_tes = histo_LowPtTES->GetBinContent(bin);

		double err = 0;
		if ( taupt >= high_pt_tau_thr ) {
			int bin_high = histo_HighPtTES->GetXaxis()->FindBin(taudecaymode);
			err = histo_HighPtTES->GetBinError(bin_high);

		} else if ( taupt > low_pt_tau_thr) {
			int bin_high = histo_HighPtTES->GetXaxis()->FindBin(taudecaymode);
			double err_high = histo_HighPtTES->GetBinError(bin_high);
			double err_low = histo_LowPtTES->GetBinError(bin);
			err = err_low + (err_high-err_low)/(high_pt_tau_thr-low_pt_tau_thr)*(taupt-low_pt_tau_thr) ;
		} else {
			err = histo_LowPtTES->GetBinError(bin);
		}

		if (type_of == "Up") { sf_tes = sf_tes+err;}
		else if(type_of == "Down") { sf_tes = sf_tes-err; }
		else {sf_tes = sf_tes*1.;}
	} // closing of "None"

	return sf_tes;
}



/*
   float LowMTRegion_2017::PZetaVis( int muindex, int tauindex){
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

   float LowMTRegion_2017::PZeta(int muindex, int tauindex , float metpx, float metpy){
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

void LowMTRegion_2017::HistoFiller(TH1F *histo, double value, double weight){                                    

	histo->Fill(value, weight);

}  



void LowMTRegion_2017::initializePileupInfo() {

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

double LowMTRegion_2017::getPileupWeight(float ntruePUInt) {
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


bool LowMTRegion_2017::MatchingToGenMuons(TLorentzVector tau, int &genindex){
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


int LowMTRegion_2017::GenTaus(){

	int HadronicTau = 0;
	for (unsigned int iGen = 0; iGen < gen_tau_had_pt.size(); iGen++){
		TLorentzVector gen_part;
		gen_part.SetPtEtaPhiE(gen_tau_had_pt.at(iGen), gen_tau_had_eta.at(iGen), gen_tau_had_phi.at(iGen), gen_tau_had_energy.at(iGen));
		if( gen_part.Pt() > 15.) { HadronicTau++; }


	}
	return HadronicTau;
}

bool LowMTRegion_2017::MatchingToGenTaus(TLorentzVector tau, int &genindex){
	double DRact = 0.2;
	genindex=-1;
	bool ismatched=false;
	for (unsigned int iGen = 0; iGen < tree->mc_tau_had_pt->size(); iGen++){
		TLorentzVector gen_part3;
		gen_part3.SetPtEtaPhiE(tree->mc_tau_had_pt->at(iGen), tree->mc_tau_had_eta->at(iGen), tree->mc_tau_had_phi->at(iGen), tree->mc_tau_had_energy->at(iGen));

		bool ismatched_to_lepton=false;
		for (unsigned int iMC2 = 0; iMC2 < tree->mc_tau_pt->size(); ++iMC2) {
			TLorentzVector pp4;
			pp4.SetPtEtaPhiE(tree->mc_tau_pt->at(iMC2), tree->mc_tau_eta->at(iMC2),tree->mc_tau_phi->at(iMC2), tree->mc_tau_energy->at(iMC2));
			if ((gen_part3.DeltaR(pp4) < 0.2 ) && (tree->mc_tau_decay->at(iMC2) == 0 || tree->mc_tau_decay->at(iMC2) == 1)) ismatched_to_lepton=true;
		}

		if(ismatched_to_lepton) continue;

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


bool LowMTRegion_2017::FillChain(TChain *chain, const char* inputFileList) {

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

float LowMTRegion_2017::GetEfficiency(float eta, float pt) 
{

double sf = 1.;


        if( fabs(eta) < 1.5) {

                int iBin = fhDMuMediumSF->FindBin(eta, pt);
                sf = fhDMuMediumSF->GetBinContent(iBin);

        } else if( fabs(eta) > 1.5) {

                int iBin = fhDMuMediumSF_endcap->FindBin(eta, pt);
                sf = fhDMuMediumSF_endcap->GetBinContent(iBin);

        }

/*	double sf = 1.;
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
*/
	return sf;
}

double LowMTRegion_2017::GetCollinearMass(const TLorentzVector &tau, const TLorentzVector &mu,  const TLorentzVector MET) {

	double METproj=fabs((MET.Px()*tau.Px()+MET.Py()*tau.Py())/tau.Pt());
	double xth=1;
	if((tau.Pt()+METproj)!=0) xth=tau.Pt()/(tau.Pt()+METproj);
	//now calculate the visibsle mass
	double mass_vis=(tau+mu).M();
	return mass_vis/sqrt(xth);

}


double LowMTRegion_2017::ElePtBinSF_CR5(double elept, double eleSCeta) {
	double SF=1;

	if(fabs(eleSCeta) < 1.4442) { if (elept >= 35. && elept < 131.6) { SF= 0.140 - (0.0029*elept) + (2.56*0.00001* pow(elept,2)) - (8.48 * 0.00000001 * pow(elept,3));} else if (elept >= 131.6 && elept < 359.3) { SF=0.020 - (0.00013*elept) + (3.50*0.0000001* pow(elept,2)) - (2.90 * 0.0000000001 * pow(elept,3));} else { SF=0.00514 + (4.73*0.0000001*elept);}} 

	if(fabs(eleSCeta) > 1.566 && fabs(eleSCeta) < 2.0) { if (elept >= 35. && elept < 125) {SF= 0.1012 - (0.00094*elept) + (3.37*0.000001* pow(elept,2));} else if (elept >= 125.0 && elept < 226.3) {SF= 0.0488 - (11.37*0.00001*elept);} else { SF =0.0241 + (1.24*0.000001*elept);}} 

	if(fabs(eleSCeta) > 2.0 && fabs(eleSCeta) < 2.5) {  if (elept >= 35. && elept < 152) { SF= 0.0622 - (0.00012*elept);} else {SF= 0.0387;}}

	return SF;
}

double LowMTRegion_2017::MutoTauFR(TLorentzVector tau, string type_unc, string apply_to_what) {
	int Otype=0;
	double DR=0.2;
	double SF_fake =1.;
	double SF_fake_1 =1.;
for (unsigned int iGen = 0; iGen < tree->mc_px->size(); iGen++){

                TLorentzVector gen_part;
                gen_part.SetPxPyPzE(tree->mc_px->at(iGen),tree->mc_py->at(iGen),tree->mc_pz->at(iGen),tree->mc_energy->at(iGen));

                bool isMuonPrompt = ( gen_part.Pt() >  8. && (abs(tree->mc_pdgId->at(iGen))==13 ) && (tree->mc_status_flags->at(iGen) >> 0 & 1  )) ? true : false ;

                bool isElectronPrompt = ( gen_part.Pt() >  8. && ( abs(tree->mc_pdgId->at(iGen))==11 ) && (tree->mc_status_flags->at(iGen) >> 0 & 1  )) ? true : false ;

                bool isElectronTau = ( gen_part.Pt() >  8. && (abs(tree->mc_pdgId->at(iGen))==11 ) && ( tree->mc_status_flags->at(iGen) >> 5 & 1 ))  ? true : false ;


                bool isMuonTau = ( gen_part.Pt() >  8. && (abs(tree->mc_pdgId->at(iGen))== 13   ) && ( tree->mc_status_flags->at(iGen) >> 5 & 1 ))   ? true : false ;


                if(isMuonPrompt || isElectronTau || isElectronPrompt || isMuonTau ) {
                        if( tau.DeltaR(gen_part) < DR ) { if(isMuonPrompt) {Otype=1;}
                                else if(isElectronPrompt) {Otype=2;} else if( isMuonTau) {Otype=3;} else if(isElectronTau) {Otype=4;}
                        }
                }
        }

  if(Otype == 1 || Otype == 3) {

                double sf_r_1 = 0.970189; double err_r_1 = 0.152359;
                double sf_r_2 = 0.643275; double err_r_2 = 0.247892;
                double sf_r_3 = 0.663525; double err_r_3 = 0.255819;
                double sf_r_4 = 1.08482 ; double err_r_4 = 0.45344;
                double sf_r_5 = 4.1917 ; double err_r_5 = 0.811688;
                if( fabs(tau.Eta()) > 0 && fabs(tau.Eta()) < 0.4) { if(type_unc == "normal") {SF_fake = sf_r_1;} else if(type_unc == "up" && apply_to_what == "mutau") {SF_fake= sf_r_1+err_r_1;} else if(type_unc == "down" && apply_to_what == "mutau") { SF_fake = sf_r_1-err_r_1;}  } //1.47; }     // 1.22
                if( fabs(tau.Eta()) > 0.4 && fabs(tau.Eta()) < 0.8) { if(type_unc == "normal")  {SF_fake = sf_r_2;} else if(type_unc == "up" && apply_to_what == "mutau" ) {SF_fake= sf_r_2+err_r_2;} else if(type_unc == "down" && apply_to_what == "mutau" ) { SF_fake = sf_r_2-err_r_2;}}
                if( fabs(tau.Eta()) > 0.8 && fabs(tau.Eta()) < 1.2) { if(type_unc == "normal")  {SF_fake = sf_r_3;} else if(type_unc == "up" && apply_to_what == "mutau") {SF_fake= sf_r_3+err_r_3;} else if(type_unc == "down" && apply_to_what == "mutau" ) { SF_fake = sf_r_3-err_r_3;}}
                if( fabs(tau.Eta()) > 1.2 && fabs(tau.Eta()) < 1.7) { if(type_unc == "normal")  {SF_fake = sf_r_4;} else if(type_unc == "up" && apply_to_what == "mutau") {SF_fake= sf_r_4+err_r_4;} else if(type_unc == "down" && apply_to_what == "mutau") { SF_fake = sf_r_4-err_r_4;}}
                if( fabs(tau.Eta()) > 1.7 && fabs(tau.Eta()) < 2.3) { if(type_unc == "normal") {SF_fake = sf_r_5;} else if(type_unc == "up" && apply_to_what == "mutau") {SF_fake= sf_r_5+err_r_5;} else if(type_unc == "down" && apply_to_what == "mutau" ) { SF_fake = sf_r_5-err_r_5;}}



}  else if ( Otype == 2 || Otype == 4) {

        double sf_e_1 = 0.96 ; double err_e_1 = 0.11;
        double sf_e_2 =0.91;  double err_e_2 = 0.2;
        if( fabs(tau.Eta())  < 1.460) {  if(type_unc == "normal")  {SF_fake = sf_e_1;} else if(type_unc == "up" && apply_to_what == "etau") {SF_fake= sf_e_1+err_e_1;} else if(type_unc == "down" && apply_to_what == "etau") { SF_fake = sf_e_1-err_e_1;}}
        if( fabs(tau.Eta())  > 1.558) { if(type_unc == "normal")  {SF_fake = sf_e_2;} else if(type_unc == "up" && apply_to_what == "etau" ) {SF_fake= sf_e_2+err_e_2;} else if(type_unc == "down" && apply_to_what == "etau") { SF_fake = sf_e_2-err_e_2;}}



} else { SF_fake = 1.;}

return SF_fake;

	// muon faking tau
}
bool LowMTRegion_2017::matchedToGenObjetcs(TLorentzVector tau, unsigned int &typey) {
	bool ismatched= false;
	double DR=0.2;
	typey = 0;

	for (unsigned int iGen = 0; iGen < tree->mc_px->size(); iGen++){

		TLorentzVector gen_part;
		gen_part.SetPxPyPzE(tree->mc_px->at(iGen),tree->mc_py->at(iGen),tree->mc_pz->at(iGen),tree->mc_energy->at(iGen));


		bool isMuonPrompt = ( gen_part.Pt() >  8. && (abs(tree->mc_pdgId->at(iGen))==13 )) ? true : false ;

		bool isElectronPrompt = ( gen_part.Pt() >  8. && ( abs(tree->mc_pdgId->at(iGen))==11 )  ) ? true : false ;

		bool isElectronTau = ( gen_part.Pt() >  8. && (abs(tree->mc_pdgId->at(iGen))==11 ) )   ? true : false ;


		bool isMuonTau = ( gen_part.Pt() >  8. && (abs(tree->mc_pdgId->at(iGen))== 13   ) )   ? true : false ;

		if(isMuonPrompt || isElectronTau || isElectronPrompt || isMuonTau ) {
			if( tau.DeltaR(gen_part) < DR ) { if(isMuonPrompt) {typey=1;} 
				else if(isElectronPrompt) {typey=2;} else if( isMuonTau) {typey=3;} else if(isElectronTau) {typey=4;} ismatched = true;  break;}


		}
	}
	return ismatched;

}


double LowMTRegion_2017::FakeRate_SSMtLow(double taupt, double jetpt, TH2F *hist) {
	double SF=0.2;
	double reweight = 0;
	int iBin = hist->FindBin(taupt, jetpt);
	SF = hist->GetBinContent(iBin);
	if (SF != 1) reweight = SF/(1-SF);
	return reweight;

}

double LowMTRegion_2017::GetTTWeight(double ttmass) {
	double SFwt = 1;
	double event_wt=1;
	if (  ttmass >= 1600.0 && ttmass < 1850.0) event_wt = 0.682164780174015;
	if ( ttmass >= 140.0 && ttmass <  160.0) event_wt  =  0.9729072028001798;
	if ( ttmass >= 120.0 && ttmass <  140.0) event_wt  =  0.9797922180747272;
	if ( ttmass >= 760.0 && ttmass <  840.0) event_wt  =  0.8864372437465602 ;
	if ( ttmass >= 50.0 && ttmass <  60.0) event_wt  =  0.9813981902499459 ;
	if ( ttmass >= 2100.0 && ttmass <  2600.0) event_wt  =  0.7427499671322355;
	if ( ttmass >= 440.0 && ttmass <  480.0) event_wt  =  0.9750529374682331 ;
	if ( ttmass >= 1450.0 && ttmass <  1600.0) event_wt  =  0.7714163491726237;
	if ( ttmass >= 360.0 && ttmass <  400.0) event_wt  =  0.977370893836097 ;
	if ( ttmass >= 60.0 && ttmass <  70.0) event_wt  =  0.9611380009719884;
	if ( ttmass >= 520.0 && ttmass <  560.0) event_wt  =  0.8977197712106587;
	if ( ttmass >= 1150.0 && ttmass <  1300.0) event_wt  =  0.7666836900323273;
	if ( ttmass >= 560.0 && ttmass <  600.0) event_wt  =  0.8973398844702225;
	if ( ttmass >= 90.0 && ttmass <  100.0) event_wt  =  0.9652106763561925 ;
	if ( ttmass >= 0.0 && ttmass <  50.0) event_wt  =  1.0 ;
	if ( ttmass >= 480.0 && ttmass <  520.0) event_wt  =  0.9542435445703805 ;
	if ( ttmass >= 160.0 && ttmass <  180.0) event_wt  =  0.9792590199734508 ;
	if ( ttmass >= 100.0 && ttmass <  120.0) event_wt  =  0.9700352340530081 ;
	if ( ttmass >= 680.0 && ttmass <  760.0) event_wt  =  0.8843463646929697 ;
	if ( ttmass >= 1000.0 && ttmass <  1150.0) event_wt  =  0.8318764104972942 ;
	if ( ttmass >= 600.0 && ttmass <  680.0) event_wt  =  0.8989140828429507 ;
	if ( ttmass >= 80.0 && ttmass <  90.0) event_wt  =  0.9590738695444241 ;
	if ( ttmass >= 400.0 && ttmass <  440.0) event_wt  =  0.9694798809384655;
	if ( ttmass >= 280.0 && ttmass <  320.0) event_wt  =  0.9702775331129567 ;
	if ( ttmass >= 200.0 && ttmass <  240.0) event_wt  =  0.9739080265098249 ;
	if ( ttmass >= 320.0 && ttmass <  360.0) event_wt  =  0.9772284832734035 ;
	if ( ttmass >= 1850.0 && ttmass <  2100.0) event_wt  =  0.7797289419934624 ;
	if ( ttmass >= 840.0 && ttmass <  920.0) event_wt  =  0.841430547544908 ;
	if ( ttmass >= 1300.0 && ttmass <  1450.0) event_wt  =  0.7910389715317062;
	if ( ttmass >= 180.0 && ttmass <  200.0) event_wt  =  0.9682079938409833 ;
	if ( ttmass >= 920.0 && ttmass <  1000.0) event_wt  =  0.8423712691558474 ;
	if ( ttmass >= 240.0 && ttmass <  280.0) event_wt  =  0.9729173152048058 ;
	if ( ttmass >= 70.0 && ttmass <  80.0) event_wt  =  0.9586548879419882;

	SFwt = event_wt;
	return SFwt;
}


bool LowMTRegion_2017::TauMatchedToJet(TLorentzVector tau_p4, unsigned int &matched_jet_indx) {

	bool matched_to_reco_jet=false;
	TLorentzVector jet_p4(0.,0.,0.,0.);
	for (unsigned int ijet = 0; ijet < tree->jet_pt->size(); ijet++){
		if(!(fabs(tree->jet_eta->at(ijet)) < 2.3)) continue;
		if(!(tree->jet_isJetID_2017->at(ijet))) continue;
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

DEFINE_FWK_MODULE(LowMTRegion_2017);
