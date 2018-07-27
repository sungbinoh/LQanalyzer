/***************************************************************************
 * @Project: LQFR_syst_cal Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/


/// Local includes
#include "FR_syst_cal.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"

//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (FR_syst_cal);


/**
 *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
 *
 */
FR_syst_cal::FR_syst_cal() :  AnalyzerCore(),  out_electrons(0) {


  // To have the correct name in the log:                                                                                                                            
  SetLogName("FR_syst_cal");

  Message("In FR_syst_cal constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

}


void FR_syst_cal::InitialiseAnalysis() throw( LQError ) {
  
  /// Initialise histograms
  MakeHistograms();  
  //
  // You can out put messages simply with Message function. Message( "comment", output_level)   output_level can be VERBOSE/INFO/DEBUG/WARNING 
  // You can also use m_logger << level << "comment" << int/double  << LQLogger::endmsg;
  //

  MakeCleverHistograms(sighist_m, "SingleMuonSNU_prompt");
  MakeCleverHistograms(sighist_m, "SingleMuonGent_loose_prompt");
  MakeCleverHistograms(sighist_m, "SingleMuonGent_tight_prompt");
  return;
}


void FR_syst_cal::ExecuteEvents()throw( LQError ){
  
  ListTriggersAvailable();

  /*bool SBt_JSf = CheckEventComparison("suoh","JS_fake_mu22119_periodC_SKDoubleMuon_hnfake_cat_v8-0-7_FRCalculator_Mu_dxysig_DILEP",
          "suoh","SB_fake_mu22118_periodC_SKDoubleMuon_hnfake_cat_v8-0-7_FR_syst_cal", false);
  bool SBf_JSt = CheckEventComparison("suoh","JS_fake_mu22119_periodC_SKDoubleMuon_hnfake_cat_v8-0-7_FRCalculator_Mu_dxysig_DILEP",
                                      "suoh","SB_fake_mu22118_periodC_SKDoubleMuon_hnfake_cat_v8-0-7_FR_syst_cal", true);
  
  if(!SBt_JSf && !SBf_JSt) return;
  cout << "====================================" << endl;
  cout << "RunNumber : " << eventbase->GetEvent().RunNumber() << ", EventNumber : " << eventbase->GetEvent().EventNumber() << endl;
  if(SBt_JSf){
  cout << "(SB && !JS)" << endl;
    FillEventComparisonFile("JSf_SBt_diff");
  }
  if(SBf_JSt){
  cout << "(!SB && JS)" << endl;
    FillEventComparisonFile("JSt_SBf_diff");
  }
  */
  
  
  //if(eventbase->GetEvent().RunNumber() != 276282 && eventbase->GetEvent().RunNumber() != 276283) return;
  
  //// Initial event cuts
  /// MET FIleters
  if(!PassMETFilter()) return;
  
  /// Require good promary vertex
  if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex
  numberVertices = eventbase->GetEvent().nVertices();
  
  std::vector<TString> HLTs_single;
  HLTs_single.push_back("HLT_Mu20_v");
  HLTs_single.push_back("HLT_Mu27_v");
  HLTs_single.push_back("HLT_Mu50_v");
  if(!PassTriggerOR(HLTs_single)) return;
  
  std::vector<snu::KElectron> elColl = GetElectrons("ELECTRON_HN_VETO");
  if(elColl.size() != 0) return; // no e 
  
  if(!isData) weight*=MCweight;
  
  //cout << weight << endl;
  std::vector<snu::KMuon> muColl = GetMuons("MUON_SUSY_VETO",true);
  CorrectedMETRochester(muColl);
  snu::KEvent Evt = eventbase->GetEvent();
  METauto = Evt.PFMET();
  METphiauto = Evt.METPhi();
  //cout << "after correcting" << endl;
  //cout << "MET : " << Evt.PFMET() << ", METphi : " << Evt.METPhi() << endl;
  
  std::vector<snu::KMuon> loose_mu_noMiniiso = GetMuons("MUON_SUSY_VETO",true);
  std::vector<snu::KMuon> loose_mu_central;
  for(unsigned int i = 0; i < loose_mu_noMiniiso.size(); i++){
    if(loose_mu_noMiniiso.at(i).RelMiniIso() < 0.6) loose_mu_central.push_back(loose_mu_noMiniiso.at(i));
  }
  
  std::vector<snu::KMuon> loosest_mu = GetMuons("MUON_SUSY_VETO",true);
  double reliso_cut[] = {0.6, 0.3, 0.4};
  double dxy_cut[] = {999., 0.1, 1.0, 10.};
  double dz_cut[] = {999., 1.0, 10};
  double dxy_sip_cut[] = {999., 8., 10., 15.};

  int reliso_i[]  = {1, 2, 0, 0, 0, 0, 0, 0, 0, 0};
  int dxy_i[]     = {0, 0, 1, 2, 3, 0, 0, 0, 0, 0};
  int dz_i[]      = {0, 0, 0, 0, 0, 1, 2, 0, 0, 0};
  int dxy_sip_i[] = {0, 0, 0, 0, 0, 0, 0, 1, 2, 3};

  TString cut_string[] = {"iso_p3", "iso_p4", "dxy_p1", "dxy_3", "dxy_10", "dz_1", "dz_10", "SIP_8", "SIP_10", "SIP_15"};
  
  for(int i = 0; i < 10; i ++){
    std::vector<snu::KMuon> this_mu;
    for(int j = 0; j < loosest_mu.size(); j++){
      if(loosest_mu.at(j).RelMiniIso() < reliso_cut[reliso_i[i]] && fabs(loosest_mu.at(j).dXY()) < dxy_cut[dxy_i[i]] && fabs(loosest_mu.at(j).dZ()) < dz_cut[dz_i[i]]
	 && fabs(loosest_mu.at(j).dXYSig()) < dxy_sip_cut[dxy_sip_i[i]]){
	this_mu.push_back(loosest_mu.at(j));
      }//if pass id cuts
    }//for j
    RunFakes("MUON_SUSY_TIGHT", "MUON_SUSY_VETO_" + cut_string[i], this_mu, false);
  }//for i
  
  
  //RunFakes("MUON_HN_TIGHT", "MUON_HN_LOOSE");
  RunFakes("MUON_SUSY_TIGHT", "MUON_SUSY_VETO_central", loose_mu_central, false);
}



void FR_syst_cal::RunFakes(TString tight_ID, TString loose_ID, std::vector<snu::KMuon> loose_mu, bool debugging){
  
  //bool debugging = std::find(k_flags.begin(), k_flags.end(), "Debug") != k_flags.end();
  
  
  FillHist("cutflow_22", 0.5, 1., 0., 10., 10);
  
  std::map< TString, std::vector<double> > HLT_ptrange; //defining muon pt range for each single muon trigger
  std::map< TString, double > HLT_ptmin;
  HLT_ptrange["HLT_Mu20_v"].push_back(35.);
  HLT_ptrange["HLT_Mu20_v"].push_back(45.);
  HLT_ptmin["HLT_Mu20_v"] = 23.;
  
  HLT_ptrange["HLT_Mu27_v"].push_back(45.);
  HLT_ptrange["HLT_Mu27_v"].push_back(80.);
  HLT_ptmin["HLT_Mu27_v"] = 30.;
  
  HLT_ptrange["HLT_Mu50_v"].push_back(80.);
  HLT_ptrange["HLT_Mu50_v"].push_back(9999.);
  HLT_ptmin["HLT_Mu50_v"] = 55.;
  
  std::vector<TString> AllHLTs;
  for(std::map< TString, std::vector<double> >::iterator it=HLT_ptrange.begin(); it!=HLT_ptrange.end(); it++){
    AllHLTs.push_back(it->first);
  }
  
  bool DijetFake = std::find(k_flags.begin(), k_flags.end(), "DijetFake") != k_flags.end(); //flag for QCD MC fr cal
  bool DijetPrompt = std::find(k_flags.begin(), k_flags.end(), "DijetPrompt") != k_flags.end(); //flag for MC subtraction
  
  std::vector<snu::KJet> jets = GetJets("JET_HN");
  int For_HLT_Mu3_PFJet40_v = 0;
  int n_bjets = 0;
  for(unsigned int i = 0; i < jets.size(); i++){
    if( IsBTagged(jets.at(i), snu::KJet::CSVv2, snu::KJet::Medium) ) n_bjets++;
    if( jets.at(i).Pt() > 50. ) For_HLT_Mu3_PFJet40_v++;
  }

  float pileup_reweight=(1.0);
  if(!k_isdata){
    pileup_reweight = mcdata_correction->CatPileupWeight(eventbase->GetEvent(),0);
  }

  


  FillHist("cutflow_22", 1.5, 1., 0., 10., 10);
  
  std::vector<snu::KMuon> veto_mu = GetMuons("MUON_HN_VETO",true);
  //if(loose_mu.size() != 1) return; //require only 1 muon
  if(debugging) cout << "veto_mu.size() != 1 : " << veto_mu.size() <<endl;
  if(veto_mu.size() != 1) return;
  if(debugging) cout << "jets.size() : " << jets.size() << endl;
  if(jets.size() < 1) return; // Njet >= 1
  
  std::vector<snu::KMuon> hnloose;//keep only prompt muon for MC subtraction, keep only fake muon for QCD fr cal
  hnloose.clear();
 
  if(debugging) cout << "loose_mu.size() : " << loose_mu.size() << endl;
  if(loose_mu.size() != 1) return;
  
  for(unsigned int i=0; i<loose_mu.size(); i++){
    if(k_isdata) hnloose.push_back( loose_mu.at(i) );
    if(!k_isdata){
      if(DijetFake){
	if( !TruthMatched( loose_mu.at(i)) ){
	  hnloose.push_back( loose_mu.at(i) );
	}
      }
      else if(DijetPrompt){
	if( TruthMatched( loose_mu.at(i)) ){
	  hnloose.push_back( loose_mu.at(i) );
	}
      }
      else{
	return;
      }
    }
  }


  if(debugging) cout << "hnloose.size() : " << hnloose.size() << endl; 
  if(hnloose.size() != 1) return;
  
  FillHist("cutflow_22", 2.5, 1., 0., 10., 10);
  
  snu::KMuon muon = hnloose.at(0);
  TLorentzVector metvec;
  metvec.SetPtEtaPhiE( METauto, 0, METphiauto, METauto);
  double MTval = AnalyzerCore::MT( muon, metvec );
  if(debugging) cout << "METauto : " << METauto << ", MTval : " << MTval << endl;
  
  TString current_trigger;
  double weight_by_pt(0.);//calculate weight for this event                                                                                                                        
  for(std::map< TString, std::vector<double> >::iterator it=HLT_ptrange.begin(); it!=HLT_ptrange.end(); it++){
    double tmp = GetTriggerWeightByPtRange(it->first, it->second, HLT_ptmin[it->first], hnloose, For_HLT_Mu3_PFJet40_v);
    if(tmp > 0.1){
      weight_by_pt = tmp;
      current_trigger = it->first;
      break;
    }
  }
  double Nvtx_reweight = 1.;
  snu::KEvent Event = eventbase->GetEvent();
  int nvtx = Event.nVertices();
  if(!k_isdata && weight_by_pt > 0.1){
    Nvtx_reweight = mcdata_correction -> GetVtxReweight(current_trigger, nvtx);
  }
  
  if(DijetPrompt && (METauto > 20 || MTval > 25) ) return;//give event cut for MC subtraction case
  
  double this_weight = weight*pileup_reweight;
  bool IsThisTight = PassID( muon, tight_ID) && (muon.RelMiniIso() < 0.2);

  double conept = muon.Pt()*(1 + max(0.,(muon.RelMiniIso() - 0.2) ) );
  
  int awayjet_syst_n = 4;
  if(!loose_ID.Contains("MUON_SUSY_VETO_central")){
    awayjet_syst_n = 1;
  }
  awayjet_syst_n = 1;
  
  double awayjet_pt[] = {40., 20., 30., 60}; //away jet pt cuts for systematics, central value is 40 GeV 
  double awayjet_dphi[] = {1.5, 2.0, 3.0};
  TString awayjet_dphi_str[] = {"1p5", "2p0", "3p0"};
  
  if(debugging) cout << "start awayjet for loop" << endl;
  for(int i = 0; i < awayjet_syst_n; i ++){
    double AwayjetPt = awayjet_pt[i];
    bool histfilled = false; //fill only once for a event with many away jets
    
    FillHist("Weight_VS_Pt_Awayjet_"+TString::Itoa(AwayjetPt,10)+"_all_" + loose_ID, muon.Pt(), weight_by_pt, 1., 0., 100., 100., 0., 300., 60.);
    FillHist("Weight_VS_Pt_cone_Awayjet_"+TString::Itoa(AwayjetPt,10)+"_all_" + loose_ID, conept, weight_by_pt, 1., 0., 100., 100., 0., 300., 60.);
        
    for(unsigned int j = 0; j < jets.size(); j++){
      if(histfilled) break;
      snu::KJet jet = jets.at(j);
      if(debugging) cout << "AwayjetPt : " << AwayjetPt << ", jet.Pt(): " << jet.Pt() << ". jets.size() : " << jets.size() << endl;
      if( jet.Pt() < AwayjetPt ) continue;
      
      double dPhi = muon.DeltaPhi( jet );
      
      bool UseEvent = false;
      if(debugging) cout << "dPhi : " << dPhi << ". (jet.Pt()/muon.Pt() : " << (jet.Pt()/muon.Pt()) << endl;
      UseEvent = (dPhi > 2.5) && (jet.Pt()/muon.Pt() > 1.);
      
      if(loose_ID.Contains("MUON_SUSY_VETO_central") && i == 0){
	FillHist("dPhi_jet_lepton_loose", dPhi, this_weight * weight_by_pt, 0., 3.5, 350);
	FillHist("Pt_jet_over_Pt_lepton_loose", jet.Pt()/muon.Pt(), this_weight * weight_by_pt, 0., 5., 500);
      }
      if(loose_ID.Contains("MUON_SUSY_VETO_central") && i == 0 && IsThisTight){
        FillHist("dPhi_jet_lepton_tight", dPhi, this_weight * weight_by_pt, 0., 3.5, 350);
        FillHist("Pt_jet_over_Pt_lepton_tight", jet.Pt()/muon.Pt(), this_weight * weight_by_pt, 0., 5., 500);
      }
      
      if(debugging) cout <<  "UseEvent : " << UseEvent << endl;
      
      if( UseEvent ){
	if(loose_ID.Contains("LOOSEv7_SIP3") && i == 0){
	  //FillEventComparisonFile("SB_fake_mu");
	  //FillHist("PeroidB_RunNumber_L", eventbase->GetEvent().RunNumber(), 1., 275650, 276290, 640);
	  //if(IsThisTight) FillHist("PeroidB_RunNumber_T", eventbase->GetEvent().RunNumber(), 1., 275650, 276290, 640);
	    
	  if(debugging){
	    int Run_N = eventbase->GetEvent().RunNumber();
	    int Event_N = eventbase->GetEvent().EventNumber();
	    cout << "##########" << endl;
	    cout << "AwayjetPt : " << AwayjetPt << endl;
	    cout << loose_ID << endl;
	    cout << "Run N : " << Run_N << endl;
	    cout << "Event N : " << Event_N << endl;
	    cout << "-> MET : " << METauto << ", METphi : " << METphiauto << endl;  
	  }
	}
	TString current_hist_name = "SingleMuonTrigger_Dijet_Awayjet_" + TString::Itoa(AwayjetPt,10)+"_all_" + loose_ID;

	FillDenAndNum("SingleMuonTrigger_Dijet_Awayjet_"+TString::Itoa(AwayjetPt,10)+"_all_" + loose_ID, muon, this_weight * weight_by_pt, IsThisTight);
	
	if(IsThisTight){
	  double muon_SF = mcdata_correction->MuonScaleFactor("MUON_HN_TIGHT", hnloose, 0);
	  double trigger_SF = 1.;
	  if(current_trigger.Contains("HLT_Mu17")){
	    trigger_SF = 0.975082;
	  }
	  else if(current_trigger.Contains("HLT_Mu8")){
	    trigger_SF = 0.961195;      
	  }
	  else if(current_trigger.Contains("HLT_Mu8")){
	    trigger_SF = 1.0166;
	  }
	  else trigger_SF = 1.;
	  double SFs = muon_SF * trigger_SF;
	  if(k_isdata) SFs = 1.;
	  FillDenAndNum(current_hist_name + "_Nvtx_reweight", muon, this_weight * weight_by_pt * Nvtx_reweight * SFs, IsThisTight);
	}// if istight = apply SFs for tight lepton
	else FillDenAndNum(current_hist_name + "_Nvtx_reweight", muon, this_weight * weight_by_pt * Nvtx_reweight, IsThisTight);
	
	histfilled = true;
      }
    } // END Tag jet loop
  }//for loop away jet PTs

  if(!loose_ID.Contains("MUON_HN_LOOSEv7_SIP3")) return;
  
  //Away jet dphi syst
  for(int i = 0; i < 3; i ++){
    double Awayjet_dphi = awayjet_dphi[i];
    bool histfilled = false;
    for(unsigned int j = 0; j < jets.size(); j++){
      if(histfilled) break;
      snu::KJet jet = jets.at(j);
      if( jet.Pt() < 40. ) continue;
      double dPhi = muon.DeltaPhi( jet );

      bool UseEvent = false;
      UseEvent = (dPhi > Awayjet_dphi) && (jet.Pt()/muon.Pt() > 1.);
      if( UseEvent ){
        FillDenAndNum("SingleMuonTrigger_Dijet_Awayjet_Dphi"+ awayjet_dphi_str[i] +"_all_" + loose_ID, muon, this_weight * weight_by_pt, IsThisTight);
	histfilled = true;
      }
    }// END Tag jet loop
  }//for lppt away jet dphi
  

  //cout << "=====================" << endl;

}


double FR_syst_cal::GetTriggerWeightByPtRange(TString hltname, vector<double> ptrange, double trigger_safe_pt, std::vector<snu::KMuon> muons, int npfjet50){

  double prescale_trigger = 0.;
  
  if(muons.size()==1){
    snu::KMuon muon = muons.at(0);
    double min_pt_cone = ptrange.at(0);
    double max_pt_cone = ptrange.at(1);
    bool SafePt = (muon.MiniAODPt() > trigger_safe_pt);

    if(SafePt && PassTrigger(hltname)){
      double TightISO = 0.2;
      double conept = muon.Pt()*(1 + max(0.,(muon.RelMiniIso()-TightISO) ) );

      if(conept >= min_pt_cone && conept < max_pt_cone){
        if(hltname=="HLT_Mu3_PFJet40_v"){
	  prescale_trigger = WeightByTrigger(hltname, TargetLumi) ;
	}
	else if(hltname=="HLT_Mu8_TrkIsoVVL_v"){
	  prescale_trigger = WeightByTrigger(hltname, TargetLumi) ;
	}
	else prescale_trigger = WeightByTrigger(hltname, TargetLumi) ;
      }
      
      if(hltname=="HLT_Mu3_PFJet40_v"){
        if(npfjet50==0) prescale_trigger = 0.;
      }
    }
    
  }
  
  if(prescale_trigger == 0.) return 0.;
  if(k_isdata) return 1.;
  return prescale_trigger;

}


void FR_syst_cal::FillDenAndNum(TString prefix, snu::KMuon muon, double thisweight, bool isTight){
  
  //float etaarray [] = {0.0, 0.8, 1.479, 2.0, 2.5};
  //float ptarray [] = {0., 5., 12., 15., 20., 25., 30., 35., 45., 60., 100., 200.};
  
  //new binning for new pt ranges 
  float etaarray [] = {0.0, 0.8, 1.479, 2.5};
  //float ptarray [] = {0., 5., 10., 15., 25., 35., 50., 70.};
  //float ptarray [] = {0., 10., 20., 30., 40., 50., 60., 70.};
  float ptarray [] = {5., 10., 20., 30., 40., 50., 60., 70.};
  //float ptaaray [] = {5., 12., };
  

  float etaarray_2 [] = {0.0, 1.479, 2.5};
  float ptarray_2 [] = {10.,15.,40.,200.};
  
  double TightISO = 0.2;
  double conept = muon.Pt()*(1 + max(0.,(muon.RelMiniIso()-TightISO) ) );

  /*
  FillHist(prefix+"_eta_F0", muon.Eta(), thisweight, -3., 3., 30);
  FillHist(prefix+"_pt_F0", muon.Pt(), thisweight, 0., 200., 200);
  FillHist(prefix+"_pt_cone_F0", conept, thisweight, 0., 200., 200);
  FillHist(prefix+"_RelIso_F0", muon.RelIso04(), thisweight, 0., 1., 100);
  FillHist(prefix+"_Chi2_F0", muon.GlobalChi2(), thisweight, 0., 50., 50);
  FillHist(prefix+"_dXY_F0", fabs(muon.dXY()), thisweight, 0., 1., 1000);
  FillHist(prefix+"_dXYSig_F0", fabs(muon.dXYSig()), thisweight, 0., 15., 150);
  FillHist(prefix+"_dZ_F0", fabs(muon.dZ()), thisweight, 0., 0.5, 50);
  FillHist(prefix+"_onebin_F0", 0.5, thisweight, 0., 1., 1);
  */
  //FillHist(prefix+"_events_pt_vs_eta_F0", muon.Pt(), fabs(muon.Eta()), thisweight, ptarray, 11, etaarray, 4);
  //FillHist(prefix+"_events_pt_cone_vs_eta_F0", conept, fabs(muon.Eta()), thisweight, ptarray, 11, etaarray, 4);
  FillHist(prefix+"_events_pt_vs_eta_F0", muon.Pt(), fabs(muon.Eta()), thisweight, ptarray, 7, etaarray, 3);
  FillHist(prefix+"_events_pt_cone_vs_eta_F0", conept, fabs(muon.Eta()), thisweight, ptarray, 7, etaarray, 3);
  FillHist(prefix+"_PFMET_F0", METauto, thisweight, 0., 1000., 1000);
  
  if( isTight ){
    /*
    FillHist(prefix+"_eta_F", muon.Eta(), thisweight, -3., 3., 30);
    FillHist(prefix+"_pt_F", muon.Pt(), thisweight, 0., 200., 200);
    FillHist(prefix+"_pt_cone_F", conept, thisweight, 0., 200., 200);
    FillHist(prefix+"_RelIso_F", muon.RelIso04(), thisweight, 0., 1., 100);
    FillHist(prefix+"_Chi2_F", muon.GlobalChi2(), thisweight, 0., 50., 50);
    FillHist(prefix+"_dXY_F", fabs(muon.dXY()), thisweight, 0., 1., 1000);
    FillHist(prefix+"_dXYSig_F", fabs(muon.dXYSig()), thisweight, 0., 15., 150);
    FillHist(prefix+"_dZ_F", fabs(muon.dZ()), thisweight, 0., 0.5, 50);
    FillHist(prefix+"_onebin_F", 0.5, thisweight, 0., 1., 1);
    */
    //FillHist(prefix+"_events_pt_vs_eta_F", muon.Pt(), fabs(muon.Eta()), thisweight, ptarray, 11, etaarray, 4);
    //FillHist(prefix+"_events_pt_cone_vs_eta_F", conept, fabs(muon.Eta()), thisweight, ptarray, 11, etaarray, 4);
    FillHist(prefix+"_events_pt_vs_eta_F", muon.Pt(), fabs(muon.Eta()), thisweight, ptarray, 7, etaarray, 3);
    FillHist(prefix+"_events_pt_cone_vs_eta_F", conept, fabs(muon.Eta()), thisweight, ptarray, 7, etaarray, 3);
    FillHist(prefix+"_PFMET_F", METauto, thisweight, 0., 1000., 1000);
    
  }
  

}


void FR_syst_cal::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);
  m_logger<< INFO << "Number of events that pass 1 7GeV trigger = " << n_17_pass  << LQLogger::endmsg;
  m_logger<< INFO << "Number of events that pass 17 GeV + jet trigger = " << n_17_jet_pass  << LQLogger::endmsg;
  m_logger<< INFO << "Number of events that pass 17 GeV || jet trigger = " << n_17_17_jet_pass  << LQLogger::endmsg;

}

void FR_syst_cal::BeginCycle() throw( LQError ){
  
  Message("In begin Cycle", INFO);
  
  //If you wish to output variables to output file use DeclareVariable
  // clear these variables in ::ClearOutputVectors function
  //DeclareVariable(obj, label, treename );
  //DeclareVariable(obj, label ); //-> will use default treename: LQTree
  //DeclareVariable(out_electrons, "Signal_Electrons", "LQTree");
  //DeclareVariable(out_muons, "Signal_Muons");

  n_17_jet_pass=0;
  n_17_17_jet_pass=0;
  n_17_pass=0;

  
  return;
  
}

FR_syst_cal::~FR_syst_cal() {
  
  Message("In FR_syst_cal Destructor" , INFO);
  
}



void FR_syst_cal::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


///############### THESE ARE FUNCTIONS SPECIFIC TO THIS CYCLE

void FR_syst_cal::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this FR_syst_calCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void FR_syst_cal::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}




