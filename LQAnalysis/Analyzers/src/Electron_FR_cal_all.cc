/***************************************************************************
 * @Project: LQElectron_FR_cal_all Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/


/// Local includes
#include "Electron_FR_cal_all.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"

//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (Electron_FR_cal_all);


/**
 *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
 *
 */
Electron_FR_cal_all::Electron_FR_cal_all() :  AnalyzerCore(),  out_electrons(0) {


  // To have the correct name in the log:                                                                                                                            
  SetLogName("Electron_FR_cal_all");

  Message("In Electron_FR_cal_all constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

}


void Electron_FR_cal_all::InitialiseAnalysis() throw( LQError ) {
  
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


void Electron_FR_cal_all::ExecuteEvents()throw( LQError ){
  
  cout << "ahah" << endl;
  //cout << "--------------------------------" << endl;
  
  /*bool SBt_JSf = CheckEventComparison("suoh","JS_fake_mu22119_periodC_SKDoubleMuon_hnfake_cat_v8-0-7_FRCalculator_Mu_dxysig_DILEP",
				      "suoh","SB_fake_mu22118_periodC_SKDoubleMuon_hnfake_cat_v8-0-7_Electron_FR_cal_all", false);
  bool SBf_JSt = CheckEventComparison("suoh","JS_fake_mu22119_periodC_SKDoubleMuon_hnfake_cat_v8-0-7_FRCalculator_Mu_dxysig_DILEP",
                                      "suoh","SB_fake_mu22118_periodC_SKDoubleMuon_hnfake_cat_v8-0-7_Electron_FR_cal_all", true);
  
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
  
  cout << "RunNumber : " << eventbase->GetEvent().RunNumber() << ", EventNumber : " << eventbase->GetEvent().EventNumber() << endl;
  
  //if(eventbase->GetEvent().RunNumber() != 276282 && eventbase->GetEvent().RunNumber() != 276283) return;
  
  //// Initial event cuts
  /// MET FIleters
  if(!PassMETFilter()) return;
  
  /// Require good promary vertex
  if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex
  numberVertices = eventbase->GetEvent().nVertices();
  
  std::vector<TString> HLTs_single;
  HLTs_single.push_back("HLT_Photon22_v");
  HLTs_single.push_back("HLT_Photon30_v");
  HLTs_single.push_back("HLT_Photon36_v");
  HLTs_single.push_back("HLT_Photon50_v");
  HLTs_single.push_back("HLT_Photon75_v");
  HLTs_single.push_back("HLT_Photon90_v");
  HLTs_single.push_back("HLT_Photon120_v");
  HLTs_single.push_back("HLT_Photon175_v");

  FillHist("HLT_Pass", 0.5, 1., 0., 20., 20);
  std::vector<snu::KElectron> electron_nocut = GetElectrons(true, true, "ELECTRON_NOCUT");
  unsigned int N_electron_nocut = 0;
  N_electron_nocut = electron_nocut.size();
  double el_pt_max = 0.;
  if(N_electron_nocut > 0){
    for(int i_el_nocut = 0; i_el_nocut < N_electron_nocut; i_el_nocut++){
      double current_el_pt = electron_nocut.at(i_el_nocut).Pt();
      if(current_el_pt > el_pt_max) el_pt_max = current_el_pt;
    }
  }
  
  if(PassTrigger("HLT_Photon22_v")){
    if(el_pt_max > 190) FillHist("HLT_Pass", 1.5, 1., 0., 20., 20); 
    if(N_electron_nocut > 0) FillHist("El_Pt_MAX_HLT_Photon22_v", el_pt_max, 1., 0., 500., 500);
  }  
  if(PassTrigger("HLT_Photon30_v")){
    if(el_pt_max > 190) FillHist("HLT_Pass", 2.5, 1., 0., 20., 20);
    if(N_electron_nocut > 0) FillHist("El_Pt_MAX_HLT_Photon30_v", el_pt_max, 1., 0., 500., 500);
  }  
  if(PassTrigger("HLT_Photon36_v")){ 
    if(el_pt_max > 190) FillHist("HLT_Pass", 3.5, 1., 0., 20., 20);
    if(N_electron_nocut > 0) FillHist("El_Pt_MAX_HLT_Photon36_v", el_pt_max, 1., 0., 500., 500);
  }  
  if(PassTrigger("HLT_Photon50_v")){
    if(el_pt_max > 190) FillHist("HLT_Pass", 4.5, 1., 0., 20., 20);
    if(N_electron_nocut > 0) FillHist("El_Pt_MAX_HLT_Photon50_v", el_pt_max, 1., 0., 500., 500);
  }  
  if(PassTrigger("HLT_Photon75_v")){
    if(el_pt_max > 190) FillHist("HLT_Pass", 5.5, 1., 0., 20., 20);
    if(N_electron_nocut > 0) FillHist("El_Pt_MAX_HLT_Photon75_v", el_pt_max, 1., 0., 500., 500);
  }
  if(PassTrigger("HLT_Photon90_v")){
    if(el_pt_max > 190) FillHist("HLT_Pass", 6.5, 1., 0., 20., 20);
    if(N_electron_nocut > 0) FillHist("El_Pt_MAX_HLT_Photon90_v", el_pt_max, 1., 0., 500., 500);
  }  
  if(PassTrigger("HLT_Photon120_v")){
    if(el_pt_max > 190) FillHist("HLT_Pass", 7.5, 1., 0., 20., 20);
    if(N_electron_nocut > 0) FillHist("El_Pt_MAX_HLT_Photon120_v", el_pt_max, 1., 0., 500., 500);
  }
  if(PassTrigger("HLT_Photon175_v")){
    if(el_pt_max > 190) FillHist("HLT_Pass", 8.5, 1., 0., 20., 20);
    if(N_electron_nocut > 0) FillHist("El_Pt_MAX_HLT_Photon175_v", el_pt_max, 1., 0., 500., 500);
  }
  
  
  
  if(!PassTriggerOR(HLTs_single)) return;
  
  std::vector<snu::KMuon> muColl = GetMuons("MUON_HN_VETO");
  if(muColl.size() != 0) return;//no mu
  
  if(!isData) weight*=MCweight;
  
  //cout << weight << endl;
  std::vector<snu::KElectron> elColl = GetElectrons("ELECTRON_SUSY_HNPAIR_VETO",true);
  CorrectedMETRochester(muColl);
  snu::KEvent Evt = eventbase->GetEvent();
  METauto = Evt.PFMET();
  METphiauto = Evt.METPhi();
  //cout << "after correcting" << endl;
  //cout << "MET : " << Evt.PFMET() << ", METphi : " << Evt.METPhi() << endl;
  
  std::vector<snu::KElectron> loose_electron = GetElectrons(false, true, "ELECTRON_SUSY_HNPAIR_LOOSE");
  std::vector<snu::KElectron> loosest_electron = GetElectrons(false, true, "ELECTRON_SUSY_HNPAIR_VETO");
    
  double reliso_cut[] = {0.6, 0.3, 0.4};
  double dxy_b_cut[] = {0.05, 0.03, 0.07, 0.09};
  double dxy_e_cut[] = {0.1, 0.08, 0.12, 0.14};
  double dz_b_cut[] = {0.1, 0.05, 0.15};
  double dz_e_cut[] = {0.2, 0.15, 0.25};
  //double dxy_sip_cut[] = {999., 8., 10., 15.};
  
  int reliso_i[]  = {1, 2, 0, 0, 0, 0, 0, 0, 0, 0};
  int dxy_i[]   = {0, 0, 1, 2, 3, 0, 0, 0, 0, 0};
  int dz_i[]    = {0, 0, 0, 0, 0, 1, 2, 0, 0, 0};
  //int dxy_sip_i[] = {0, 0, 0, 0, 0, 0, 0, 1, 2, 3};

  TString cut_string[] = {"iso_p3", "iso_p4", "dxy_minus_p02", "dxy_plus_p02", "dxy_plus_p04", "dz_minus_p05", "dz_plus_p05"};
  
  int N_cut_syst = 7;
  for(int i_cut_syst = 0; i_cut_syst < N_cut_syst; i_cut_syst ++){
    std::vector<snu::KElectron> this_el;
    for(int j_loose_el = 0; j_loose_el < loosest_electron.size(); j_loose_el ++){
      if(loosest_electron.at(j_loose_el).PFRelMiniIso(false) < reliso_cut[reliso_i[i_cut_syst]] ){
	if( fabs( loosest_electron.at(j_loose_el).SCEta() ) < 1.479 && fabs(loosest_electron.at(j_loose_el).dxy()) < dxy_b_cut[dxy_i[i_cut_syst]] && fabs(loosest_electron.at(j_loose_el).dz()) < dz_b_cut[dz_i[i_cut_syst]]){//if berrel & passing IP cuts
	  this_el.push_back(loosest_electron.at(j_loose_el));
	}
	if( fabs( loosest_electron.at(j_loose_el).SCEta() ) > 1.479 && fabs(loosest_electron.at(j_loose_el).dxy()) < dxy_e_cut[dxy_i[i_cut_syst]] && fabs(loosest_electron.at(j_loose_el).dz()) < dz_e_cut[dz_i[i_cut_syst]]){//if endcap & passing IP cuts
	  this_el.push_back(loosest_electron.at(j_loose_el));
	}
      }//if pass iso cut
    }//for j
    //cout << "RunFakes_ID" << endl;
    RunFakes_ID("ELECTRON_SUSY_HNPAIR_TIGHT", "ELECTRON_SUSY_HNPAIR_LOOSE_" + cut_string[i_cut_syst], this_el, false);
  }//for i
  
  //cout << "RunFakes start" << endl;
  RunFakes("ELECTRON_SUSY_HNPAIR_TIGHT", "ELECTRON_SUSY_HNPAIR_LOOSE_central", loose_electron, false);
  
  
}

//function to run for loose ID systematics : FR only for away jet pt 40 GeV
void Electron_FR_cal_all::RunFakes_ID(TString tight_ID, TString loose_ID, std::vector<snu::KElectron> loose_el, bool debugging){
  std::map< TString, std::vector<double> > HLT_ptrange; //defining muon pt range for each single muon trigger
  std::map< TString, double > HLT_ptmin;

  HLT_ptrange["HLT_Photon22_v"].push_back(40.);
  HLT_ptrange["HLT_Photon22_v"].push_back(55.);
  HLT_ptmin["HLT_Photon22_v"] = 25.;
  
  HLT_ptrange["HLT_Photon30_v"].push_back(55.);
  HLT_ptrange["HLT_Photon30_v"].push_back(65.);
  HLT_ptmin["HLT_Photon30_v"] = 35.;

  HLT_ptrange["HLT_Photon36_v"].push_back(65.);
  HLT_ptrange["HLT_Photon36_v"].push_back(85.);
  HLT_ptmin["HLT_Photon36_v"] = 40.;

  HLT_ptrange["HLT_Photon50_v"].push_back(85.);
  HLT_ptrange["HLT_Photon50_v"].push_back(120.);
  HLT_ptmin["HLT_Photon50_v"] = 55.;

  HLT_ptrange["HLT_Photon75_v"].push_back(120.);
  HLT_ptrange["HLT_Photon75_v"].push_back(150.);
  HLT_ptmin["HLT_Photon75_v"] = 80.;

  HLT_ptrange["HLT_Photon90_v"].push_back(150.);
  HLT_ptrange["HLT_Photon90_v"].push_back(195.);
  HLT_ptmin["HLT_Photon90_v"] = 100.;

  HLT_ptrange["HLT_Photon120_v"].push_back(195.);
  HLT_ptrange["HLT_Photon120_v"].push_back(280.);
  HLT_ptmin["HLT_Photon120_v"] = 130.;

  HLT_ptrange["HLT_Photon175_v"].push_back(280.);
  HLT_ptrange["HLT_Photon175_v"].push_back(9999.);
  HLT_ptmin["HLT_Photon175_v"] = 190.;
  
  
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
  
  if(debugging) cout << "jets.size() : " << jets.size() << endl;
  if(jets.size() < 1) return; // Njet >= 1
  
  std::vector<snu::KElectron> hnloose;//keep only prompt electron for MC subtraction, keep only fake electron for QCD fr cal
  hnloose.clear();
  
  if(debugging) cout << "loose_el.size() : " << loose_el.size() << endl;
  if(loose_el.size() != 1) return;
  
  for(unsigned int i=0; i<loose_el.size(); i++){
    if(k_isdata) hnloose.push_back( loose_el.at(i) );
    if(!k_isdata){
      if(DijetFake){
	if( !TruthMatched( loose_el.at(i), false) ){
          hnloose.push_back( loose_el.at(i) );
        }
      }
      else if(DijetPrompt){
        if( TruthMatched( loose_el.at(i), false) ){
	  //if(loose_el.at(i).GetType() == 1 || loose_el.at(i).GetType() == 8){  
	  hnloose.push_back( loose_el.at(i) );
	}
      }
      else{
        return;
      }
    }
  }
  
  if(debugging) cout << "hnloose.size() : " << hnloose.size() << endl;
  if(hnloose.size() != 1) return;
  
  snu::KElectron electron = hnloose.at(0);
  TLorentzVector metvec;
  metvec.SetPtEtaPhiE( METauto, 0, METphiauto, METauto);
  double MTval = AnalyzerCore::MT( electron, metvec );
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

  //if(DijetPrompt && (METauto > 20 || MTval > 25) ) return;//give event cut for MC subtraction case
  if(DijetPrompt && (METauto > 80 || MTval > 25) ) return;
  

  double this_weight = weight*pileup_reweight;
  bool IsThisTight = PassID( electron, tight_ID) && (electron.PFRelMiniIso(false) < 0.1);
  double conept = electron.Pt()*(1 + max(0.,(electron.PFRelMiniIso(false) - 0.1) ) );
  
  double AwayjetPt = 40.;
  bool histfilled = false;
  for(unsigned int j = 0; j < jets.size(); j++){//for loop over all jets in the event
    if(histfilled) break;
    snu::KJet jet = jets.at(j);
    if(debugging) cout << "AwayjetPt : " << AwayjetPt << ", jet.Pt(): " << jet.Pt() << ". jets.size() : " << jets.size() << endl;
    if( jet.Pt() < AwayjetPt ) continue;

    double dPhi = electron.DeltaPhi( jet );

    bool UseEvent = false;
    if(debugging) cout << "dPhi : " << dPhi << ". (jet.Pt()/electron.Pt() : " << (jet.Pt()/electron.Pt()) << endl;
    UseEvent = (dPhi > 2.5) && (jet.ChargedEMEnergyFraction() < 0.65) && (jet.Pt()/electron.Pt() > 1.0);
    
    if(debugging) cout <<  "UseEvent : " << UseEvent << endl;
    
    if( UseEvent ){
      TString current_hist_name = "SingleElectronTrigger_Dijet_Awayjet_" + TString::Itoa(AwayjetPt,10)+"_all_" + loose_ID;
      double electron_RecoSF = mcdata_correction->ElectronRecoScaleFactor(hnloose);
            
      if(IsThisTight){
	//double electron_id_sf = mcdata_correction->ElectronScaleFactor("", hnloose, 0);
	double electron_id_sf = 1.;
	double electron_SF = electron_id_sf * electron_RecoSF;
	double SFs = electron_SF;
	if(k_isdata) SFs = 1.;
	//FillDenAndNum(current_hist_name + "_Nvtx_reweight", muon, this_weight * weight_by_pt * Nvtx_reweight * SFs, IsThisTight);
	FillDenAndNum(current_hist_name + "_Nvtx_reweight", electron, this_weight * weight_by_pt * SFs, IsThisTight);

      }// if istight = apply SFs for tight lepton
      else FillDenAndNum(current_hist_name + "_Nvtx_reweight", electron, this_weight * weight_by_pt * electron_RecoSF, IsThisTight);

      histfilled = true;
    }//end of if(UseEvent)
  }//end for jets
}
////////////////////////////////////////////////////////RunFakes_ID ends

void Electron_FR_cal_all::RunFakes(TString tight_ID, TString loose_ID, std::vector<snu::KElectron> loose_el, bool debugging){
  
  std::map< TString, std::vector<double> > HLT_ptrange; //defining muon pt range for each single muon trigger
  std::map< TString, double > HLT_ptmin;
  HLT_ptrange["HLT_Photon22_v"].push_back(40.);
  HLT_ptrange["HLT_Photon22_v"].push_back(55.);
  HLT_ptmin["HLT_Photon22_v"] = 25.;
  
  HLT_ptrange["HLT_Photon30_v"].push_back(55.);
  HLT_ptrange["HLT_Photon30_v"].push_back(65.);
  HLT_ptmin["HLT_Photon30_v"] = 35.;
  
  HLT_ptrange["HLT_Photon36_v"].push_back(65.);
  HLT_ptrange["HLT_Photon36_v"].push_back(85.);
  HLT_ptmin["HLT_Photon36_v"] = 40.;
  
  HLT_ptrange["HLT_Photon50_v"].push_back(85.);
  HLT_ptrange["HLT_Photon50_v"].push_back(120.);
  HLT_ptmin["HLT_Photon50_v"] = 55.;
  
  HLT_ptrange["HLT_Photon75_v"].push_back(120.);
  HLT_ptrange["HLT_Photon75_v"].push_back(150.);
  HLT_ptmin["HLT_Photon75_v"] = 80.;

  HLT_ptrange["HLT_Photon90_v"].push_back(150.);
  HLT_ptrange["HLT_Photon90_v"].push_back(195.);
  HLT_ptmin["HLT_Photon90_v"] = 100.;

  HLT_ptrange["HLT_Photon120_v"].push_back(195.);
  HLT_ptrange["HLT_Photon120_v"].push_back(280.);
  HLT_ptmin["HLT_Photon120_v"] = 130.;

  HLT_ptrange["HLT_Photon175_v"].push_back(280.);
  HLT_ptrange["HLT_Photon175_v"].push_back(9999.);
  HLT_ptmin["HLT_Photon175_v"] = 190.;
  
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
  
  //normalize test
  if(DijetPrompt){
    double m_Z = 91.1876;
    std::vector< snu::KElectron > tightelectrons = GetElectrons(false, true, "ELECTRON_SUSY_HNPAIR_TIGHT");
    
    double electron_id_iso_sf = 1.;
    double electron_id_iso_sf_syst[2] = {1., 1.};
    double electron_RecoSF = 1.;
    
    if(!k_isdata){
      electron_id_iso_sf = 1.;
      electron_id_iso_sf_syst[0] = 1.;
      electron_id_iso_sf_syst[1] = 1.;
      electron_RecoSF = mcdata_correction->ElectronRecoScaleFactor(tightelectrons, 0);
    }
    TString electron_id_iso_sf_string[2] = {"ID_Up", "ID_Down"};
        
    TLorentzVector metvec;
    metvec.SetPtEtaPhiE( METauto, 0, METphiauto, METauto );
    
    for(unsigned int i = 0; i < AllHLTs.size(); i++){
      TString ThisTrigger = AllHLTs.at(i);
      snu::KEvent Event = eventbase->GetEvent();
      int nvtx = Event.nVertices();
      
      if( PassTrigger(ThisTrigger) && (loose_el.size() > 0 ) ){
	double max_loose_el_pt = 0.;
	for(unsigned int i_loose_e = 0; i_loose_e < loose_el.size(); i_loose_e ++){
	  if(loose_el.at(i_loose_e).Pt() > max_loose_el_pt) max_loose_el_pt = loose_el.at(i_loose_e).Pt();
	}
	
	if(HLT_ptmin[ThisTrigger] > max_loose_el_pt) continue;
	double triggerweight = 0.;
	if(ThisTrigger.Contains("_Photon22_") ) triggerweight = 1.599;
	else if(ThisTrigger.Contains("_Photon30_") ) triggerweight = 6.615;
	else if(ThisTrigger.Contains("_Photon36_") ) triggerweight = 13.141;
        else if(ThisTrigger.Contains("_Photon50_") ) triggerweight = 31.008;
        else if(ThisTrigger.Contains("_Photon75_") ) triggerweight = 134.618;
        else if(ThisTrigger.Contains("_Photon90_") ) triggerweight = 264.215;
        else if(ThisTrigger.Contains("_Photon120_") ) triggerweight = 537.352;
        else if(ThisTrigger.Contains("_Photon175_") ) triggerweight = 35860.065;
	else triggerweight = 0.;
	
	double this_weight = weight * electron_id_iso_sf * electron_RecoSF * triggerweight * pileup_reweight;
        double weight_noID = weight * electron_RecoSF * triggerweight * pileup_reweight;
	
	if(k_isdata){
	  this_weight = 1.;
	  weight_noID = 1.;
	}
	
	snu::KEvent Event = eventbase->GetEvent();
	int nvtx = Event.nVertices();
	
	//==== 1) Z-Peak
        if(tightelectrons.size() == 2){
          double mll = (tightelectrons.at(0) + tightelectrons.at(1)).M();
          //if( (tightelectrons.at(0).Pt() > HLT_ptmin[ThisTrigger] ) && (tightelectrons.at(1).Pt() > 20.) && (fabs(mll-m_Z) < 10.) ){
	  if( (tightelectrons.at(0).Pt() > HLT_ptmin[ThisTrigger]) && (tightelectrons.at(1).Pt() > 20.) && (fabs(mll-m_Z) < 10.) ){
	    
	    FillHist(ThisTrigger+"_ZPeak_mll", mll, this_weight, 0., 200., 200);
            FillHist(ThisTrigger+"_ZPeak_leadpt", tightelectrons.at(0).Pt(), this_weight, 0., 500., 500);
            FillHist(ThisTrigger+"_ZPeak_subleadpt", tightelectrons.at(1).Pt(), this_weight, 0., 500., 500);
	    FillHist(ThisTrigger+"_ZPeak_Nvtx", nvtx, this_weight, 0., 70., 70);
	    FillHist(ThisTrigger+"_ZPeak_PFMET", METauto, this_weight, 0., 500., 500);
	    
	    for(int ID_i = 0; ID_i < 2; ID_i++){
	      FillHist(ThisTrigger+"_ZPeak_mll_Nvtx_reweight_" + electron_id_iso_sf_string[ID_i], mll, weight_noID * electron_id_iso_sf_syst[ID_i], 0., 200., 200);
              FillHist(ThisTrigger+"_ZPeak_leadpt_Nvtx_reweight_" + electron_id_iso_sf_string[ID_i], tightelectrons.at(0).Pt(), weight_noID * electron_id_iso_sf_syst[ID_i], 0., 500., 500);
              FillHist(ThisTrigger+"_ZPeak_subleadpt_Nvtx_reweight_" + electron_id_iso_sf_string[ID_i], tightelectrons.at(1).Pt(), weight_noID * electron_id_iso_sf_syst[ID_i], 0., 500., 500);
              FillHist(ThisTrigger+"_ZPeak_Nvtx_Nvtx_reweight_" + electron_id_iso_sf_string[ID_i], nvtx, weight_noID * electron_id_iso_sf_syst[ID_i], 0., 70., 70);
              FillHist(ThisTrigger+"_ZPeak_PFMET_Nvtx_reweight_" + electron_id_iso_sf_string[ID_i], METauto, weight_noID * electron_id_iso_sf_syst[ID_i], 0., 500., 500);
	    }
	    double MTval = AnalyzerCore::MT( loose_el.at(0), metvec );
	    if( (MTval > 60.) && (MTval < 100.)){//MT cut applied Z region
	      FillHist(ThisTrigger+"_ZPeak_MTcut_mll", mll, this_weight, 0., 200., 200);
	      FillHist(ThisTrigger+"_ZPeak_MTcut_leadpt", tightelectrons.at(0).Pt(), this_weight, 0., 500., 500);
	      FillHist(ThisTrigger+"_ZPeak_MTcut_subleadpt", tightelectrons.at(1).Pt(), this_weight, 0., 500., 500);
	      FillHist(ThisTrigger+"_ZPeak_MTcut_Nvtx", nvtx, this_weight, 0., 70., 70);
	      FillHist(ThisTrigger+"_ZPeak_MTcut_PFMET", METauto, this_weight, 0., 500., 500);
	    }//MT cut
	  }
	}
	
	//==== 2) W
	if(loose_el.size() ==1 && loose_el.at(0).Pt() > HLT_ptmin[ThisTrigger]){
	  double MTval = AnalyzerCore::MT( loose_el.at(0), metvec );
	  double current_conept = loose_el.at(0).Pt()*(1 + max(0.,(loose_el.at(0).PFRelMiniIso(false)-0.1) ) );
	  
	  if((MTval > 60.) && (MTval < 100.)){
	    FillHist(ThisTrigger + "_Ptcone_MET_inclusive_loose", current_conept, weight_noID, 0., 200., 200);
	    FillHist(ThisTrigger + "_Eta_MET_inclusive_loose", loose_el.at(0).Eta(),  weight_noID, -3., 3., 100);
	    FillHist(ThisTrigger + "_MET_MET_inclusive_tight_loose", METauto, weight_noID, 0., 500., 500);
	    FillHist(ThisTrigger + "_MT_MET_inclusive_tight_loose", MTval, weight_noID, 0., 500., 500);
	    
	    bool IsThisTight = PassID( loose_el.at(0), tight_ID);
	    if(IsThisTight){
	      FillHist(ThisTrigger + "_Ptcone_MET_inclusive_tight", current_conept, this_weight, 0., 200., 200);
	      FillHist(ThisTrigger + "_Eta_MET_inclusive_tight", loose_el.at(0).Eta(),  this_weight, -3., 3., 100);
	      FillHist(ThisTrigger + "_MET_MET_inclusive_tight", METauto, this_weight, 0., 500., 500);
	      FillHist(ThisTrigger + "_MT_MET_inclusive_tight", MTval, this_weight, 0., 500., 500);
	    }
	  }	
	}
	
	if(tightelectrons.size() == 1 && tightelectrons.at(0).Pt() > HLT_ptmin[ThisTrigger]){
          double MTval = AnalyzerCore::MT( tightelectrons.at(0), metvec );
          FillHist(ThisTrigger + "MT_MET_MT_inclusive_tight", MTval, this_weight, 0., 500., 500);
	  
	  if( (METauto > 40.) && (MTval > 60.) && (MTval < 100.) && (tightelectrons.at(0).Pt() > 20.)){
	    FillHist(ThisTrigger+"_W_John_PFMET", METauto, this_weight, 0., 500., 500);
            FillHist(ThisTrigger+"_W_John_MT", MTval, this_weight, 0., 500., 500);
            FillHist(ThisTrigger+"_W_John_leadpt", tightelectrons.at(0).Pt(), this_weight, 0., 500., 500);
	    FillHist(ThisTrigger+"_W_John_Nvtx", nvtx, this_weight, 0., 70., 70);
	    
	    for(int ID_i = 0; ID_i < 2; ID_i++){
	      FillHist(ThisTrigger+"_W_John_PFMET_Nvtx_reweight_" + electron_id_iso_sf_string[ID_i], METauto, weight_noID * electron_id_iso_sf_syst[ID_i], 0., 500., 500);
              FillHist(ThisTrigger+"_W_John_MT_Nvtx_reweight_" + electron_id_iso_sf_string[ID_i], MTval, weight_noID * electron_id_iso_sf_syst[ID_i], 0., 500., 500);
              FillHist(ThisTrigger+"_W_John_leadpt_Nvtx_reweight_" + electron_id_iso_sf_string[ID_i], tightelectrons.at(0).Pt(), weight_noID * electron_id_iso_sf_syst[ID_i], 0., 500., 500);
              FillHist(ThisTrigger+"_W_John_Nvtx_Nvtx_reweight_" + electron_id_iso_sf_string[ID_i], nvtx, weight_noID * electron_id_iso_sf_syst[ID_i], 0., 70., 70);
	    }
	  }	      
        }
      }
    }
  }
  //normaize test ends
  //if(debugging) cout << "trigger" << endl;
  
  if(debugging) cout << "jets.size() : " << jets.size() << endl;
  if(jets.size() < 1) return; // Njet >= 1
  
  std::vector<snu::KElectron> hnloose;//keep only prompt muon for MC subtraction, keep only fake muon for QCD fr cal
  hnloose.clear();
 
  if(debugging) cout << "loose_el.size() : " << loose_el.size() << endl;
  if(loose_el.size() != 1) return;
  
  for(unsigned int i=0; i<loose_el.size(); i++){
    if(k_isdata) hnloose.push_back( loose_el.at(i) );
    if(!k_isdata){
      if(DijetFake){
	if( !TruthMatched( loose_el.at(i), false) ){
          hnloose.push_back( loose_el.at(i) );
        }
      }
      else if(DijetPrompt){
        if( TruthMatched( loose_el.at(i), false) ){
	  //if(loose_el.at(i).GetType() == 1 || loose_el.at(i).GetType() == 8){  
	  hnloose.push_back( loose_el.at(i) );
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
  
  snu::KElectron electron = hnloose.at(0);
  TLorentzVector metvec;
  metvec.SetPtEtaPhiE( METauto, 0, METphiauto, METauto);
  double MTval = AnalyzerCore::MT( electron, metvec );
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
  
  if(DijetPrompt && (METauto > 80 || MTval > 25) ) return;//give event cut for MC subtraction case
  //test for MET , MT cut on FR
  
  double this_weight = weight*pileup_reweight;
  bool IsThisTight = PassID( electron, tight_ID) && (electron.PFRelMiniIso(false) < 0.1);
  double conept = electron.Pt()*(1 + max(0.,(electron.PFRelMiniIso(false) - 0.1) ) );
  
  int awayjet_syst_n = 11;
  double awayjet_pt[] = {40., 20., 30., 60}; //away jet pt cuts for systematics, central value is 40 GeV 
  double awayjet_dphi[] = {2.5, 1.5, 2.0, 3.0};//away jet dphi cut for systematics, central value is 2.5
  double pj_over_pl[] = {1.0, 0.8, 1.2, 1.5, 2.0};//pt(j)/pt(l) cut for systematics, central value is 1.5
  int awayjet_pt_i[] =   {1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0};
  int awayjet_dphi_i[] = {0, 0, 0, 1, 2, 3, 0, 0, 0, 0, 0};
  int pj_over_pl_i[] =   {0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 0};
  //TString awayjet_syst_str[] = {"awayjet_pt_20", "awayjet_pt_30", "awayjet_pt_60", "awayjet_dphi_1p5", "awayjet_dphi_2", "awayjet_dphi_3", "pj_over_pl_1", "pj_over_pl_2",
  //"central"};
  TString awayjet_syst_str[] = {"awayjet_pt_20", "awayjet_pt_30", "awayjet_pt_60", "awayjet_dphi_1p5", "awayjet_dphi_2", "awayjet_dphi_3",
                                "pj_over_pl_0p8", "pj_over_pl_1p2", "pj_over_pl_1p5", "pj_over_pl_2",
                                "central"};
  
  if(debugging) cout << "start awayjet for loop" << endl;
  //awayjet systematic for loop starts
  for(int i = 0; i < awayjet_syst_n; i ++){
    //i = 8; //run only for central
    //cuts for dijet topology, and string of current syst.
    double Awayjet_Pt_cut = awayjet_pt[ awayjet_pt_i[i] ];
    double Awayjet_dphi_cut = awayjet_dphi[ awayjet_dphi_i[i] ];
    double pj_over_pl_cut = pj_over_pl[ pj_over_pl_i[i] ];
    TString current_syst_str = awayjet_syst_str[i];
    
    bool histfilled = false; 
    
    for(unsigned int j = 0; j < jets.size(); j++){//for loop all jets
      if(histfilled) break;
      snu::KJet jet = jets.at(j);
      if(debugging) cout << "AwayjetPt : " << Awayjet_Pt_cut << ", jet.Pt(): " << jet.Pt() << ". jets.size() : " << jets.size() << endl;
      if( jet.Pt() < Awayjet_Pt_cut ) continue;
      double dPhi = electron.DeltaPhi( jet );
      
      bool UseEvent = false;
      if(debugging) cout << "dPhi : " << dPhi << ". (jet.Pt()/electron.Pt() : " << (jet.Pt()/electron.Pt()) << endl;
      UseEvent = (dPhi > 2.5) && (jet.ChargedEMEnergyFraction() < 0.65) && (jet.Pt()/electron.Pt() > 1.0);
      
      if(i == 10){//if central
	FillHist("dPhi_jet_lepton_loosse", dPhi, this_weight * weight_by_pt, 0., 3.5, 350);
        FillHist("Pt_jet_over_Pt_lepton_loose", jet.Pt()/electron.Pt(), this_weight * weight_by_pt, 0., 5., 500);
        FillHist("ChargedEMEnergyFraction_loose", jet.ChargedEMEnergyFraction(), this_weight * weight_by_pt, 0., 2., 200);
	if(IsThisTight){
	  FillHist("dPhi_jet_lepton_tight", dPhi, this_weight * weight_by_pt, 0., 3.5, 350);
	  FillHist("Pt_jet_over_Pt_lepton_tight", jet.Pt()/electron.Pt(), this_weight * weight_by_pt, 0., 5., 500);
	  FillHist("ChargedEMEnergyFraction_tight", jet.ChargedEMEnergyFraction(), this_weight * weight_by_pt, 0., 2., 200);
	}
      }
      
      if(debugging) cout <<  "UseEvent : " << UseEvent << endl;
      
      if( UseEvent ){//if pass cuts
	TString current_hist_name = "SingleElectronTrigger_Dijet_" + current_syst_str + "_all_" + loose_ID;
	double electron_RecoSF = mcdata_correction->ElectronRecoScaleFactor(hnloose);
	
	if(i == 10){//if central
	  FillHist("MET_VS_MT_UseEvent_loose", METauto, MTval, this_weight * weight_by_pt * electron_RecoSF, 0., 200., 20., 0., 200., 20);
          FillHist("MET_VS_Ptcone_UseEvent_loose", METauto, conept, this_weight * weight_by_pt * electron_RecoSF, 0., 200., 20., 0., 200., 20);
          FillHist("MET_VS_Eta_UseEvent_loose", METauto, electron.Eta(), this_weight * weight_by_pt * electron_RecoSF, 0., 200., 20., -3., 3., 100);
          FillHist("Ptcone_UseEvent_loose", conept, this_weight * weight_by_pt * electron_RecoSF, 0., 500., 500);
	  FillHist("Eta_UseEvent_loose", electron.Eta(), this_weight * weight_by_pt * electron_RecoSF, -3., 3., 100);
	}
	
	if(IsThisTight){
	  double electron_id_sf = 1.;
	  double electron_SF = electron_id_sf * electron_RecoSF;
	  double trigger_SF = 1.;
	  double SFs = electron_SF * trigger_SF;
          if(k_isdata) SFs = 1.;
          FillDenAndNum(current_hist_name + "_Nvtx_reweight", electron, this_weight * weight_by_pt * SFs, IsThisTight);
	  
	  snu::KEvent Evt_current = eventbase->GetEvent();
	  
	  if( i == 11){
	    cout << "RunNumber = " << eventbase->GetEvent().RunNumber() << endl;
	    cout << "EventNumber = " << eventbase->GetEvent().EventNumber() << endl;
	    cout << "trigger = " << current_trigger << endl;
	    cout << "electron.Pt() = " << electron.Pt() << endl;
	    cout << "electron.PFRelIso(0.3) = " << electron.PFRelIso(0.3) << endl;
	    cout << "electron.dXYSig() = " << electron.dxySig() << endl;
	    cout << "electron.SCEta() = " << electron.SCEta() << endl;
	    cout << "jet.Pt() = " << jet.Pt() << endl;
	    cout << "AwayjetPt = " << Awayjet_Pt_cut << endl;
	    cout << "dPhi = " << dPhi << endl;
	    cout << "jet.Pt()/electron.Pt() = " << jet.Pt()/electron.Pt() << endl;
	    cout << "jet.ChargedEMEnergyFraction() = " << jet.ChargedEMEnergyFraction() << endl;
	    cout << "METauto = " << METauto << endl;
	    cout << "MTval = " << MTval << endl;
	    cout << "this_weight = " << this_weight << endl;
	    cout << "muon_scale_sf = " << electron_SF << endl;
	    cout << "electron_id_sf = " << electron_id_sf << endl;
	    cout << "weight_by_pt = " << weight_by_pt << endl;
	    cout << "weight = " << weight << endl;
	    cout << "--> UseEvent = " << UseEvent << endl << endl;
	    cout << "total weight : " << this_weight * weight_by_pt * SFs << endl;
	  }
	  
	  if(i == 10){//if central
	    FillHist("MET_VS_MT_UseEvent_tight", METauto, MTval, this_weight * weight_by_pt * SFs, 0., 200., 20., 0., 200., 20);
	    FillHist("MET_VS_Ptcone_UseEvent_tight", METauto, conept, this_weight * weight_by_pt * SFs, 0., 200., 20., 0., 200., 20);
	    FillHist("MET_VS_Eta_UseEvent_tight", METauto, electron.Eta(), this_weight * weight_by_pt * SFs, 0., 200., 20., -3., 3., 100);
	    FillHist("Ptcone_UseEvent_tight", conept, this_weight * weight_by_pt * SFs, 0., 500., 500);
	    FillHist("Eta_UseEvent_tight", electron.Eta(), this_weight * weight_by_pt * SFs, -3., 3., 100);
	  }
	}// if istight = apply SFs for tight lepton
	else FillDenAndNum(current_hist_name + "_Nvtx_reweight", electron, this_weight * weight_by_pt * electron_RecoSF, IsThisTight);
		
        histfilled = true;
      }
    }//for loop all jets ends
  }//awayjet systematic for loop ends
  
  
}
//////////////////////////////////RunFakes ends


///functionto get JES, JER, lepton energy, lepton ID SF syst for prompt contribution
void Electron_FR_cal_all::RunFakes_subtraction(TString tight_ID, TString loose_ID, std::vector<snu::KElectron> loose_el, bool debugging){
  std::map< TString, std::vector<double> > HLT_ptrange; //defining muon pt range for each single muon trigger
  std::map< TString, double > HLT_ptmin;
  HLT_ptrange["HLT_Photon22_v"].push_back(40.);
  HLT_ptrange["HLT_Photon22_v"].push_back(55.);
  HLT_ptmin["HLT_Photon22_v"] = 25.;

  HLT_ptrange["HLT_Photon30_v"].push_back(55.);
  HLT_ptrange["HLT_Photon30_v"].push_back(65.);
  HLT_ptmin["HLT_Photon30_v"] = 35.;

  HLT_ptrange["HLT_Photon36_v"].push_back(65.);
  HLT_ptrange["HLT_Photon36_v"].push_back(85.);
  HLT_ptmin["HLT_Photon36_v"] = 40.;

  HLT_ptrange["HLT_Photon50_v"].push_back(85.);
  HLT_ptrange["HLT_Photon50_v"].push_back(120.);
  HLT_ptmin["HLT_Photon50_v"] = 55.;

  HLT_ptrange["HLT_Photon75_v"].push_back(120.);
  HLT_ptrange["HLT_Photon75_v"].push_back(150.);
  HLT_ptmin["HLT_Photon75_v"] = 80.;

  HLT_ptrange["HLT_Photon90_v"].push_back(150.);
  HLT_ptrange["HLT_Photon90_v"].push_back(195.);
  HLT_ptmin["HLT_Photon90_v"] = 100.;

  HLT_ptrange["HLT_Photon120_v"].push_back(195.);
  HLT_ptrange["HLT_Photon120_v"].push_back(280.);
  HLT_ptmin["HLT_Photon120_v"] = 130.;

  HLT_ptrange["HLT_Photon175_v"].push_back(280.);
  HLT_ptrange["HLT_Photon175_v"].push_back(9999.);
  HLT_ptmin["HLT_Photon175_v"] = 190.;
  

  std::vector<TString> AllHLTs;
  for(std::map< TString, std::vector<double> >::iterator it=HLT_ptrange.begin(); it!=HLT_ptrange.end(); it++){
    AllHLTs.push_back(it->first);
  }
  
  bool DijetFake = std::find(k_flags.begin(), k_flags.end(), "DijetFake") != k_flags.end(); //flag for QCD MC fr cal
  bool DijetPrompt = std::find(k_flags.begin(), k_flags.end(), "DijetPrompt") != k_flags.end(); //flag for MC subtraction
  
  //call loose particles
  TString MuonLooseID_loosest = "MUON_SUSY_VETO";
  TString MuonVetoID_loosest = "MUON_SUSY_VETO";
  TString ElectronLooseID_loosest = "ELECTRON_SUSY_HNPAIR_LOOSE_loosest";
  TString ElectronVetoID_loosest = "ELECTRON_SUSY_HNPAIR_VETO_loosest";
  //==== Muons
  std::vector<snu::KMuon> muons_loosest = GetMuons(MuonLooseID_loosest, true);
  std::vector<snu::KMuon> muons_veto_loosest = GetMuons(MuonVetoID_loosest, true);
  std::vector<snu::KMuon> muons_veto;
  //==== Electrons
  std::vector<snu::KElectron> electrons_loosest = GetElectrons(false, true, ElectronLooseID_loosest);
  std::vector<snu::KElectron> electrons_veto_loosest = GetElectrons(true, true, ElectronVetoID_loosest);
  //==== Jets
  std::vector<snu::KJet> jets_eta5_nolepveto_loosest = GetJets("JET_HN_eta5_nolepveto_loosest", 10., 5.);
  
  TString syst_str[] = {"_central_subt", "_ElectronEn_up", "_ElectronEn_down", "_JetEn_up", "_JetEn_down", "_JetRes_up", "_JetRes_down", "_Unclustered_up", "_Unclustered_down",
			"_ElectronIDSF_up", "_ElectronIDSF_down"};
  
  int n_syst = 11;
  snu::KEvent Evt = eventbase->GetEvent();
  
  for(int i = 0; i < n_syst; i++){//for loop systemati categories
    TString this_syst = syst_str[i];
    
    if(i == 2) METauto = Evt.PFMETShifted(snu::KEvent::JetEn, snu::KEvent::up);//_JetEn_up
    if(i == 3) METauto = Evt.PFMETShifted(snu::KEvent::JetEn, snu::KEvent::down);//_JetEn_down
    if(i == 4) METauto = Evt.PFMETShifted(snu::KEvent::JetRes, snu::KEvent::up);//_JetRes_up
    if(i == 5) METauto = Evt.PFMETShifted(snu::KEvent::JetRes, snu::KEvent::down);//_JetRes_down
    if(i == 6) METauto = Evt.PFMETShifted(snu::KEvent::Unclustered, snu::KEvent::up);//_Unclustered_up
    if(i == 7) METauto = Evt.PFMETShifted(snu::KEvent::Unclustered, snu::KEvent::down);//_Unclustered_down
    
    //====================
    //==== Make Electron
    //====================
    std::vector<snu::KElectron> electrons, electrons_veto;
    std::vector<snu::KElectron> electrons_notshifted;
    
    double this_ElectronLooseRelIso = 0.6;
    double this_ElectronVetoRelIso = 0.6;

    int ElEnDir = 0;
    if(this_syst == "_ElectronEn_up"){
      ElEnDir = 1;
      //==== Signal Leptons
      for(unsigned int j=0; j<electrons_loosest.size(); j++){
	snu::KElectron this_electron = electrons_loosest.at(j);
	this_electron.SetPtEtaPhiM( this_electron.Pt()*this_electron.PtShiftedUp(), this_electron.Eta(), this_electron.Phi(), this_electron.M() );
	double new_RelIso = this_electron.PFRelMiniIso(false)/this_electron.PtShiftedUp();
	this_electron.SetPFRelMiniIsoRho(new_RelIso);
	if( this_electron.Pt() >= 10. && new_RelIso < this_ElectronLooseRelIso ) electrons.push_back( this_electron );
      }
      //cout << this_syst << endl;
      //cout <<"El E Up" << endl;
      //cout << "electrons.size() : " << electrons.size() << endl;
      //if(electrons.size() > 0) cout << "1st electron's pt : " << electrons.at(0).Pt() << endl; 
      
      //==== Veto Leptons
      for(unsigned int j=0; j<electrons_veto_loosest.size(); j++){
	snu::KElectron this_electron = electrons_veto_loosest.at(j);
	this_electron.SetPtEtaPhiM( this_electron.Pt()*this_electron.PtShiftedUp(), this_electron.Eta(), this_electron.Phi(), this_electron.M() );
	double new_RelIso = this_electron.PFRelMiniIso(false)/this_electron.PtShiftedUp();
	this_electron.SetPFRelMiniIsoRho(new_RelIso);
	if( this_electron.Pt() >= 10. && new_RelIso < this_ElectronVetoRelIso ) electrons_veto.push_back( this_electron );
      }
    }
    else if(this_syst == "_ElectronEn_down"){
      ElEnDir = -1;
      //==== Signal Leptons
      for(unsigned int j=0; j<electrons_loosest.size(); j++){
	snu::KElectron this_electron = electrons_loosest.at(j);
	this_electron.SetPtEtaPhiM( this_electron.Pt()*this_electron.PtShiftedDown(), this_electron.Eta(), this_electron.Phi(), this_electron.M() );
	double new_RelIso = this_electron.PFRelMiniIso(false)/this_electron.PtShiftedDown();
	this_electron.SetPFRelMiniIsoRho(new_RelIso);
	if( this_electron.Pt() >= 10. && new_RelIso < this_ElectronLooseRelIso ) electrons.push_back( this_electron );
      }
      //==== Veto Leptons
      for(unsigned int j=0; j<electrons_veto_loosest.size(); j++){
	snu::KElectron this_electron = electrons_veto_loosest.at(j);
	this_electron.SetPtEtaPhiM( this_electron.Pt()*this_electron.PtShiftedDown(), this_electron.Eta(), this_electron.Phi(), this_electron.M() );
	double new_RelIso = this_electron.PFRelMiniIso(false)/this_electron.PtShiftedDown();
	this_electron.SetPFRelMiniIsoRho(new_RelIso);
	if( this_electron.Pt() >= 10. && new_RelIso < this_ElectronVetoRelIso ) electrons_veto.push_back( this_electron );
      }
    }
    //==== normal electrons
    else{
      //==== Signal Leptons
      for(unsigned int j=0; j<electrons_loosest.size(); j++){
	snu::KElectron this_electron = electrons_loosest.at(j);
	if( this_electron.Pt() >= 10. && this_electron.PFRelMiniIso(false) < this_ElectronLooseRelIso ) electrons.push_back( this_electron );
      }
      //==== Veto Leptons
      for(unsigned int j=0; j<electrons_veto_loosest.size(); j++){
	snu::KElectron this_electron = electrons_veto_loosest.at(j);
	if( this_electron.Pt() >= 10. && this_electron.PFRelMiniIso(false) < this_ElectronVetoRelIso ) electrons_veto.push_back( this_electron );
      }
    }
    
    JSCorrectedMETElectron(ElEnDir, electrons_veto_loosest, METauto, METphiauto);
    
    //===============
    //==== Make Jet
    //===============
    std::vector<snu::KJet> jets_eta5_nolepveto; // eta < 5, NO lepton-veto
    if( (this_syst == "_JetEn_up") || (this_syst == "_JetRes_up") ){
      for(unsigned int j=0; j<jets_eta5_nolepveto_loosest.size(); j++){
	snu::KJet this_jet = jets_eta5_nolepveto_loosest.at(j);

	double this_scaling = 1.;
	if(this_syst == "_JetEn_up") this_scaling = this_jet.ScaledUpEnergy();
	if(this_syst == "_JetRes_up") this_scaling = this_jet.SmearedResUp()/this_jet.SmearedRes();

	this_jet *= this_scaling;

	if(this_jet.Pt() >= 20.) jets_eta5_nolepveto.push_back(this_jet);
      }
    }
    else if( (this_syst == "_JetEn_down") || (this_syst == "_JetRes_down") ){
      for(unsigned int j=0; j<jets_eta5_nolepveto_loosest.size(); j++){
	snu::KJet this_jet = jets_eta5_nolepveto_loosest.at(j);

	double this_scaling = 1.;
	if(this_syst == "_JetEn_down") this_scaling = this_jet.ScaledDownEnergy();
	if(this_syst == "_JetRes_down") this_scaling = this_jet.SmearedResDown()/this_jet.SmearedRes();;

	this_jet *= this_scaling;

	if(this_jet.Pt() >= 20.) jets_eta5_nolepveto.push_back(this_jet);
      }
    }
    else{
      for(unsigned int j=0; j<jets_eta5_nolepveto_loosest.size(); j++){
	snu::KJet this_jet = jets_eta5_nolepveto_loosest.at(j);
	if(this_jet.Pt() >= 20.) jets_eta5_nolepveto.push_back(this_jet);
      }
    }

    std::vector<snu::KJet> jets; // eta < 2.5, lepton-veto, normla jet selection to be used as away jets
    for(unsigned int j=0; j<jets_eta5_nolepveto.size(); j++){
      snu::KJet this_jet = jets_eta5_nolepveto.at(j);
      bool IsNormalJet = fabs( this_jet.Eta() ) < 2.5;
      bool lepinside = HasLeptonInsideJet(this_jet, muons_veto, electrons_veto);
            
      if(IsNormalJet && !lepinside) jets.push_back( this_jet );
    }
    
    
    //now we have jets(normal jet selection with systematic up & down) and muons(loose ID with up & down E) & muons_veto, also METauto & METphiauto with correct shift
    int For_HLT_Mu3_PFJet40_v = 0;
    int n_bjets = 0;
    /*
    for(unsigned int i = 0; i < jets.size(); i++){
      if( IsBTagged(jets.at(i), snu::KJet::CSVv2, snu::KJet::Medium) ) n_bjets++;
      if( jets.at(i).Pt() > 50. ) For_HLT_Mu3_PFJet40_v++;
    }
    */
    if(debugging) cout << "electrons_veto.size() != 1 : " << electrons_veto.size() <<endl;
    if(electrons_veto.size() != 1) continue; //return; // one e 
    if(debugging) cout << "jets.size() : " << jets.size() << endl;
    if(jets.size() < 1) continue; //return; // Njet >= 1
    
    std::vector<snu::KElectron> hnloose;//keep only prompt muon for MC subtraction, keep only fake muon for QCD fr cal
    hnloose.clear();
    if(debugging) cout << "electrons.size() : " << electrons.size() << endl;
    if(electrons.size() != 1) continue; //return;
    for(unsigned int i=0; i<electrons.size(); i++){
      if(k_isdata) hnloose.push_back( electrons.at(i) );
      if(!k_isdata){
	if(DijetFake){
	  if( !TruthMatched( electrons.at(i), false) ){
	    hnloose.push_back( electrons.at(i) );
	  }
	}
	else if(DijetPrompt){
	  if(debugging) cout << "DijetPrompt" << endl;
	  if( TruthMatched( electrons.at(i), false) ){
	    //if(muons.at(i).GetType() == 1 || muons.at(i).GetType() == 8){
	    if(debugging) cout << "truth" << endl;
	    hnloose.push_back( electrons.at(i) );
	  }
	}
	else{
	  continue; //return;
	}
      }
    }
    
    if(debugging) cout << "hnloose.size() : " << hnloose.size() << endl;
    if(hnloose.size() != 1) continue; //return;
    
    snu::KElectron electron = hnloose.at(0);
    TLorentzVector metvec;
    metvec.SetPtEtaPhiE( METauto, 0, METphiauto, METauto);
    double MTval = AnalyzerCore::MT( electron, metvec );
    if(debugging) cout << "METauto : " << METauto << ", MTval : " << MTval << endl;
    
    /////////////////////////////////////weights
    //---------------------------------------------
    //Weight by pt range for each HLT
    //---------------------------------------------
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
    //---------------------------------------------
    //pileup_reweight (69 mb minbias X-section)
    //---------------------------------------------
    float pileup_reweight=(1.0);
    if(!k_isdata){
      pileup_reweight = mcdata_correction->CatPileupWeight(eventbase->GetEvent(),0);
    }
    //---------------------------------------------
    //Additional Nvtx correction
    //---------------------------------------------
    double Nvtx_reweight = 1.;
    //---------------------------------------------
    //Constant Normalization Factor
    //---------------------------------------------
    double constant_SF = 1.;
    //if(k_isdata) constant_SF = 1.;
    //---------------------------------------------
    //Muon ID scale-factor, ID SF is applied only for tight muon not for loose
    //---------------------------------------------
    
    bool IsThisTight = PassID( electron, tight_ID);
    
    double electron_SF = 1.;
    double electron_RecoSF = mcdata_correction->ElectronRecoScaleFactor(hnloose);

    if(IsThisTight && !k_isdata){
      
      //fix me
      //if(this_syst == "_ElectronIDSF_up")  electron_SF = mcdata_correction->SBElectronScaleFactor("ELECTRON_HN_TIGHTv4", hnloose, 1);
      //else if(this_syst == "_ElectronIDSF_down")  electron_SF = mcdata_correction->SBElectronScaleFactor("ELECTRON_HN_TIGHTv4", hnloose, -1);
      //else electron_SF = mcdata_correction->ElectronScaleFactor("ELECTRON_HN_TIGHTv4", hnloose, 0);
    
      if(this_syst == "_ElectronIDSF_up")  electron_SF = 1.;
      else if(this_syst == "_ElectronIDSF_down")  electron_SF = 1.;
      else electron_SF = 1.;
    }
    
    double SFs = electron_RecoSF;
    if(IsThisTight) SFs *= electron_SF;
    double this_weight = weight*pileup_reweight*SFs;
    
    if(DijetPrompt && (METauto > 80 || MTval > 25) ) continue; //return;//give event cut for MC subtraction case
    
    double AwayjetPt = 40.;
    bool histfilled = false; 
    for(unsigned int j = 0; j < jets.size(); j++){//for loop all jets in the event
      if(histfilled) break;
      snu::KJet jet = jets.at(j);
      if( jet.Pt() < AwayjetPt ) continue;
      
      double dPhi = electron.DeltaPhi( jet );
      
      bool UseEvent = false;
      if(debugging) cout << "dPhi : " << dPhi << ". (jet.Pt()/electron.Pt() : " << (jet.Pt()/electron.Pt()) << endl;
      UseEvent = (dPhi > 2.5) && (jet.ChargedEMEnergyFraction() < 0.65) && (jet.Pt()/electron.Pt() > 1.);

      if(debugging) cout <<  "UseEvent : " << UseEvent << endl;
      
      if( UseEvent ){
	TString current_hist_name =  "SingleElectronTrigger_Dijet" + this_syst +"_all_" + loose_ID;
	FillDenAndNum(current_hist_name, electron, this_weight * weight_by_pt, IsThisTight);
	histfilled = true;
      }
    }//for loop all jet ends
  }//for loop systematic categories ends 
}//RunFakes_subtraction ends


double Electron_FR_cal_all::GetTriggerWeightByPtRange(TString hltname, vector<double> ptrange, double trigger_safe_pt, std::vector<snu::KElectron> electrons, int npfjet50){

  double prescale_trigger = 0.;
  
  if(electrons.size()==1){
    snu::KElectron electron = electrons.at(0);
    double min_pt_cone = ptrange.at(0);
    double max_pt_cone = ptrange.at(1);
    bool SafePt = (electron.Pt() > trigger_safe_pt);
    
    if(SafePt && PassTrigger(hltname)){
      double TightISO = 0.1;
      double conept = electron.Pt()*(1 + max(0.,(electron.PFRelMiniIso(false)-TightISO) ) );

      if(conept >= min_pt_cone && conept < max_pt_cone){
	if(hltname=="HLT_Photon22_v"){
	  prescale_trigger = 1.599;
	}
	else if(hltname=="HLT_Photon30_v"){
	  prescale_trigger = 6.615;
	}
	else if(hltname=="HLT_Photon36_v"){
          prescale_trigger = 13.141;
	}
	else if(hltname=="HLT_Photon50_v"){
	  prescale_trigger = 31.008;
	}
	else if(hltname=="HLT_Photon75_v"){
	  prescale_trigger = 134.618;
	}
	else if(hltname=="HLT_Photon90_v"){
	  prescale_trigger = 264.215;
	}
	else if(hltname=="HLT_Photon120_v"){
	  prescale_trigger = 537.352;
	}
	else if(hltname=="HLT_Photon175_v"){
	  prescale_trigger = 35860.065;
	}
	else{
	  prescale_trigger = 0.;
	}
	
      }
      
    }
    
  }
  
  if(prescale_trigger == 0.) return 0.;
  if(k_isdata) return 1.;
  return prescale_trigger;

}

void Electron_FR_cal_all::FillDenAndNum(TString prefix, snu::KElectron electron, double thisweight, bool isTight){

  int n_pt_bin = 12;
  
  //float etaarray [] = {0.0, 0.8, 1.479, 2.0, 2.5};
  //float ptarray [] = {0., 5., 12., 15., 20., 25., 30., 35., 45., 60., 100., 200.};
  
  //new binning for new pt ranges 
  float etaarray [] = {0.0, 0.8, 1.479, 2.4};
  //float ptarray [] = {0., 5., 10., 15., 25., 35., 50., 70.};
  //float ptarray [] = {0., 10., 20., 30., 40., 50., 60., 70.};
  //float ptarray [] = {5., 10., 20., 30., 40., 50., 60., 70.};
  //float ptarray [] = {35, 45, 75, 80, 110, 150, 200, 250, 300, 350, 500, 1000, 1500, 2000};
  float ptarray [] = {40, 55, 65, 85, 120, 150, 195, 280, 350, 500, 1000, 1500, 2000};

  
  //float ptaaray [] = {5., 12., };
  
    
  float etaarray_2 [] = {0.0, 1.479, 2.5};
  float ptarray_2 [] = {10.,15.,40.,200.};
  
  double TightISO = 0.1;
  double conept = electron.Pt()*(1 + max(0.,(electron.PFRelMiniIso(false)-TightISO) ) );

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
  FillHist(prefix+"_events_pt_vs_eta_F0", electron.Pt(), fabs(electron.Eta()), thisweight, ptarray, n_pt_bin, etaarray, 3);
  FillHist(prefix+"_events_pt_cone_vs_eta_F0", conept, fabs(electron.Eta()), thisweight, ptarray, n_pt_bin, etaarray, 3);
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
    FillHist(prefix+"_events_pt_vs_eta_F", electron.Pt(), fabs(electron.Eta()), thisweight, ptarray, n_pt_bin, etaarray, 3);
    FillHist(prefix+"_events_pt_cone_vs_eta_F", conept, fabs(electron.Eta()), thisweight, ptarray, n_pt_bin, etaarray, 3);
    FillHist(prefix+"_PFMET_F", METauto, thisweight, 0., 1000., 1000);
  }
  

}


void Electron_FR_cal_all::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);
  m_logger<< INFO << "Number of events that pass 1 7GeV trigger = " << n_17_pass  << LQLogger::endmsg;
  m_logger<< INFO << "Number of events that pass 17 GeV + jet trigger = " << n_17_jet_pass  << LQLogger::endmsg;
  m_logger<< INFO << "Number of events that pass 17 GeV || jet trigger = " << n_17_17_jet_pass  << LQLogger::endmsg;

}

void Electron_FR_cal_all::BeginCycle() throw( LQError ){
  
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

Electron_FR_cal_all::~Electron_FR_cal_all() {
  
  Message("In Electron_FR_cal_all Destructor" , INFO);
  
}



void Electron_FR_cal_all::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


///############### THESE ARE FUNCTIONS SPECIFIC TO THIS CYCLE

void Electron_FR_cal_all::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this Electron_FR_cal_allCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void Electron_FR_cal_all::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}


