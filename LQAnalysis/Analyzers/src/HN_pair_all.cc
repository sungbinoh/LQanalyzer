// $Id: HN_pair_all.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQHN_pair_all Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "HN_pair_all.h"

//Core includes
#include "Reweight.h"
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (HN_pair_all);

/**
 *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
 *
 */
HN_pair_all::HN_pair_all() :  AnalyzerCore(), out_muons(0)  {
  // To have the correct name in the log:                                                                                                                            
  SetLogName("HN_pair_all");
  
  Message("In HN_pair_all constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
  MakeCleverHistograms(hnpairmm, "Jet_study");
  MakeCleverHistograms(hnpairmm, "SR1_DiMu");
  MakeCleverHistograms(hnpairmm, "SR1_DiEle");
  MakeCleverHistograms(hnpairmm, "SR1_EMu");
  MakeCleverHistograms(hnpairmm, "CR1_DiMu");
  MakeCleverHistograms(hnpairmm, "CR1_DiEle");
  MakeCleverHistograms(hnpairmm, "CR1_EMu");
  MakeCleverHistograms(hnpairmm, "CR2_DiMu");
  MakeCleverHistograms(hnpairmm, "CR2_DiEle");
  MakeCleverHistograms(hnpairmm, "CR2_EMu");
  MakeCleverHistograms(hnpairmm, "CR3_DiMu");
  MakeCleverHistograms(hnpairmm, "CR3_DiEle");
  MakeCleverHistograms(hnpairmm, "CR3_EMu");
  MakeCleverHistograms(hnpairmm, "CR4_DiMu");
  MakeCleverHistograms(hnpairmm, "CR4_DiEle");
  MakeCleverHistograms(hnpairmm, "CR4_EMu");
  MakeCleverHistograms(hnpairmm, "CR5_DiMu");
  MakeCleverHistograms(hnpairmm, "CR5_DiEle");
  MakeCleverHistograms(hnpairmm, "CR5_EMu");
}


void HN_pair_all::InitialiseAnalysis() throw( LQError ) {

  /// Initialise histograms
  MakeHistograms();  
  //
  // You can out put messages simply with Message function. Message( "comment", output_level)   output_level can be VERBOSE/INFO/DEBUG/WARNING 
  // You can also use m_logger << level << "comment" << int/double  << LQLogger::endmsg;
  //
  
  Message("Making clever hists for Z ->ll test code", INFO);
  
  /// only available in v7-6-X branch and newer
  //// default lumimask is silver ////
  //// In v7-6-2-(current) the default is changed to gold (since METNoHF bug)
  ///When METNoHF isfixed the default will be back to silver
  /// set to gold if you want to use gold json in analysis
  /// To set uncomment the line below:

  /*
  TDirectory* origDir = gDirectory;
  
  TString lqdir = getenv("LQANALYZER_DIR");

  TFile *file_data_fake = new TFile( lqdir+"/data/Fake/80X/Muon_Data_fake_Rate_syst.root");
  TIter next(gDirectory->GetListOfKeys());
  TKey *key;
  vector<TString> histnames;
  vector<TString> histnames_2D;
  while ((key = (TKey*)next())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (cl->InheritsFrom("TH2D")) histnames_2D.push_back(key -> GetName());
    if (!cl->InheritsFrom("TH1")) continue;
    histnames.push_back(key -> GetName());
  }
  gROOT->cd();
  
  TDirectory* tempDir = 0;
  tempDir = gROOT->mkdir("temp_dir");
  tempDir->cd();  

  Message("saving hist strings", INFO);
  for(int i = 0; i < histnames_2D.size(); i ++){
    Message(histnames_2D.at(i), INFO);
    hist_Muon_FR_syst[histnames_2D.at(i)] = (TH2D*)file_data_fake -> Get(histnames_2D.at(i))->Clone();
  }
  Message("saving hist strings ends", INFO);
  file_data_fake -> Close();
  delete file_data_fake;
  origDir->cd();
  */
  
  return;
}

void HN_pair_all::ExecuteEvents()throw( LQError ){
  
  /// Apply the gen weight 
  if(!isData) weight*=MCweight;
  
  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;
  
  //pileup reweight
  // -- Get global weights                                                                                                                                                                                                                                                      
  float pileup_reweight=(1.0);
  if(!k_isdata){
    pileup_reweight = mcdata_correction->CatPileupWeight(eventbase->GetEvent(),0);
  }

  /////////////////////////////////////////////////////////////////
  // -- weight = 1.0      to see signal eff
  /////////////////////////////////////////////////////////////////
  if(!isData) weight *= pileup_reweight;
  //if(!isData) weight = pileup_reweight; // only for signal
  
  
  // -- call truth particle
  std::vector<snu::KTruth> truthColl;
  eventbase->GetTruthSel()->Selection(truthColl);
  
  bool has_Nmix = false;
  int another_neutrino_1 = 9999999, another_neutrino_2 = 9999999;
  
  
  bool mumu_signal = std::find(k_flags.begin(), k_flags.end(), "hn_pair_mm") != k_flags.end();
  bool ee_signal = std::find(k_flags.begin(), k_flags.end(), "hn_pair_ee") != k_flags.end();
  bool tau_veto = std::find(k_flags.begin(), k_flags.end(), "tau_veto") != k_flags.end();
  
  if(mumu_signal){
    another_neutrino_1 = 9900016;
    another_neutrino_2 = 9900012;
  }
  if(ee_signal){
    another_neutrino_1 = 9900014;
    another_neutrino_2 = 9900016;
  }
  
  for(int i = 0; i < truthColl.size(); i++){
    if( (abs(truthColl.at(i).PdgId()) == another_neutrino_1) || (abs(truthColl.at(i).PdgId()) == another_neutrino_2) ) has_Nmix = true;
  }
  if(has_Nmix) return;
  
  //FillCutFlow("NoCut", weight * 35863.3);
  FillCutFlow("NoCut", weight);
  
  FillHist("signal_eff", 0.5, weight * 35863.3, 0., 40., 40); //after lepton skim
  

  if(isData) FillHist("Nvtx_nocut_data",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);
  else  FillHist("Nvtx_nocut_mc",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);
  
  ///#### CAT:::PassBasicEventCuts is updated: uses selections as described in https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFilters: If you see this is out of date please comment
  if(!PassMETFilter()) return;     /// Initial event cuts :
  FillCutFlow("EventCut", weight);
  FillHist("signal_eff", 1.5, weight * 35863.3, 0., 40., 40); //after PassMETFilter
  
  
  /// #### CAT::: triggers stored are all HLT_Ele/HLT_DoubleEle/HLT_Mu/HLT_TkMu/HLT_Photon/HLT_DoublePhoton
  if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex
  FillHist("signal_eff", 2.5, weight * 35863.3, 0., 40., 40); //after good primary vtx cut
  
  m_logger << DEBUG << "passedTrigger "<< LQLogger::endmsg;


  // -- Store Trigger Boolean
  TString mu50_trig1 = "HLT_Mu50_v";
  TString mu50_trig2 = "HLT_TkMu50_v";
  vector<TString> mu50_trignames;
  mu50_trignames.push_back(mu50_trig1);
  mu50_trignames.push_back(mu50_trig2);
  
  //TString diele_trig1="HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v";
  TString diele_trig1="HLT_DoublePhoton60_v";

  
  bool mu50_pass = PassTriggerOR(mu50_trignames); // for mumu and emu channel
  bool diele_pass = PassTrigger(diele_trig1); // for ee channel

  if(mu50_pass)  FillHist("signal_eff", 3.5, weight * 35863.3, 0., 40., 40);
  if(diele_pass)  FillHist("signal_eff", 4.5, weight * 35863.3, 0., 40., 40);


  if(!mu50_pass && !diele_pass) return;
  //cout << "dimu_pass : " << dimu_pass << ". mu50_pass ; " << mu50_pass << endl;
  //if(mu50_pass) cout << "mu50_pass" << endl;
  
  // -- Get Veto Electrons, and store number of them
  //std::vector<snu::KElectron> electrons_susy_veto = GetElectrons("ELECTRON_HN_VETO");
  std::vector<snu::KElectron> electrons_susy_veto = GetElectrons("ELECTRON_SUSY_VETO");
  std::vector<snu::KElectron> electrons_veto;
  for(int i = 0; i < electrons_susy_veto.size(); i++){
    bool pass_eta_cut = fabs(electrons_susy_veto.at(i).SCEta()) < 1.4 || fabs(electrons_susy_veto.at(i).SCEta()) > 1.6;
    if( electrons_susy_veto.at(i).PFRelMiniIso() < 0.40 && pass_eta_cut)  electrons_veto.push_back(electrons_susy_veto.at(i));
  }
  int N_veto_ele = electrons_veto.size();
  
  // -- Get Veto Muons, return if there are not exactly two veto muons
  //std::vector<snu::KMuon> muons_susy_veto = GetMuons("MUON_HN_VETO");
  std::vector<snu::KMuon> muons_susy_veto = GetMuons("MUON_SUSY_VETO");
  std::vector<snu::KMuon> muons_veto;
  for(int i = 0; i < muons_susy_veto.size(); i++){
    bool pass_eta_cut = fabs(muons_susy_veto.at(i).Eta()) < 1.4 || fabs(muons_susy_veto.at(i).Eta()) > 1.6;
    if( muons_susy_veto.at(i).RelMiniIso() < 0.40 && pass_eta_cut) muons_veto.push_back(muons_susy_veto.at(i));
  }
  CorrectedMETRochester(muons_veto);
  int N_veto_muon = muons_veto.size();
  
  // -- Get AK8 jets
  std::vector<snu::KFatJet> fatjets_loosest_UnSmeared = GetFatJets("FATJET_HN_nolepveto");
  std::vector<snu::KFatJet> fatjets_loosest = GetCorrectedFatJet(fatjets_loosest_UnSmeared); // Smear both energy and mass
  std::vector<snu::KFatJet> fatjets;
  for(int i_fat = 0; i_fat < fatjets_loosest.size(); i_fat++){
    snu::KFatJet this_jet = fatjets_loosest.at(i_fat);
    if(JSFatJetID(this_jet) && this_jet.Pt() >= 300.) fatjets.push_back(this_jet);
  }
  
  // -- Get AK4 jets
  std::vector<snu::KJet> jets_eta5_nolepveto_loosest = GetJets("JET_HN_eta5_nolepveto", 30., 5.);
  std::vector<snu::KJet> jets_lepveto;
  
  for(int i_jet = 0; i_jet < jets_eta5_nolepveto_loosest.size(); i_jet++){
    double NormalJetMaxEta = 2.7;
    
    snu::KJet this_jet = jets_eta5_nolepveto_loosest.at(i_jet);
    bool PassPUID = this_jet.PassPileUpMVA("Loose");
    bool IsNormalJet = fabs( this_jet.Eta() ) < NormalJetMaxEta;
    bool AwayFromFatJet = IsAwayFromFatJet(this_jet, fatjets);
    bool lepinside = HasLeptonInsideJet(this_jet, muons_veto, electrons_veto);
    
    if(IsNormalJet && PassPUID && AwayFromFatJet && !lepinside) jets_lepveto.push_back( this_jet ); // normal jets without lepton veto & fatjet veto
  }
  
  // -- Call MET
  float MET = eventbase->GetEvent().PFMET();

  // -- Get global weights
  double trigger_weight = WeightByTrigger("HLT_Mu50_v", TargetLumi);//weight with any unprescaled trigger

  //cout << "weight : " << weight << ", trigger_weight : " << trigger_weight << ", trigger_weight_dimu : " << trigger_weight_dimu << endl;
  
  if(!isData) weight *= trigger_weight;

  //cout << "weight : " << weight << ", dimu_pass : " << dimu_pass << endl;
  
  /////////////////////////////////////////////////////////////////
  // -- weight = 1.0      to see signal eff
  /////////////////////////////////////////////////////////////////
 
  // -- Check if run fake
  bool NonPromptRun = std::find(k_flags.begin(), k_flags.end(), "RunFake") != k_flags.end();
  
  //TString muon_tight_id = "MUON_HN_TIGHT";
  TString muon_loose_id = "MUON_HN_LOOSEv7_SIP3";
  TString muon_tight_id = "MUON_SUSY_TIGHT";
  
  std::vector<snu::KMuon> muon_sch_tight = GetMuons("MUON_HN_TIGHT",true);
  std::vector<snu::KMuon> muon_Nocut = GetMuons("MUON_NOCUT",true);
  std::vector<snu::KMuon> muons;
  for(int i = 0; i < muon_Nocut.size(); i++){
    //if(abs(muon_Nocut.at(i).MotherPdgId()) == 15) cout << muon_Nocut.at(i).MotherPdgId() << endl;
    bool tau_pass = !tau_veto || abs(muon_Nocut.at(i).MotherPdgId()) != 15;
    bool pass_eta_cut = fabs(muon_Nocut.at(i).Eta()) < 1.4 || fabs(muon_Nocut.at(i).Eta()) > 1.6;
    if( !NonPromptRun && PassID(muon_Nocut.at(i), muon_tight_id) && (muon_Nocut.at(i).RelMiniIso() < 0.20) && pass_eta_cut && tau_pass) muons.push_back(muon_Nocut.at(i));
    if( NonPromptRun  && PassID(muon_Nocut.at(i), muon_loose_id) && (muon_Nocut.at(i).RelMiniIso() < 0.20) && pass_eta_cut && tau_pass) muons.push_back(muon_Nocut.at(i));//store loose muon for fake bkg
  }

  std::vector<snu::KElectron> electron_sch_tight = GetElectrons("ELECTRON_HN_TIGHT",true);
  TString electron_tight_id = "ELECTRON_SUSY_TIGHT";
  TString electron_loose_id = "ELECTRON_HN_FAKELOOSEv7";
  std::vector<snu::KElectron> electron_Nocut = GetElectrons("ELECTRON_NOCUT", true);
  std::vector<snu::KElectron> electrons;
  for(int i = 0; i < electron_Nocut.size(); i++){
    bool tau_pass = !tau_veto || abs(electron_Nocut.at(i).MotherPdgId()) != 15;
    bool pass_eta_cut = fabs(electron_Nocut.at(i).SCEta()) < 1.4 || fabs(electron_Nocut.at(i).SCEta()) > 1.6;
    if( !NonPromptRun && PassID(electron_Nocut.at(i), electron_tight_id) && (electron_Nocut.at(i).PFRelMiniIso() < 0.10) && pass_eta_cut){ //&& (electron_Nocut.at(i).MissingHits() == 0) ){
      if( fabs(electron_Nocut.at(i).SCEta()) < 0.8 && electron_Nocut.at(i).MVA() > -0.85  ) electrons.push_back(electron_Nocut.at(i));
      else if( fabs(electron_Nocut.at(i).SCEta()) > 0.8 && fabs(electron_Nocut.at(i).SCEta()) < 1.479 && electron_Nocut.at(i).MVA() > -0.91  ) electrons.push_back(electron_Nocut.at(i));
      else if( fabs(electron_Nocut.at(i).SCEta()) > 1.479 && fabs(electron_Nocut.at(i).SCEta()) < 2.5 && electron_Nocut.at(i).MVA() > -0.83  ) electrons.push_back(electron_Nocut.at(i));
      else continue;
    }
    if( NonPromptRun  && PassID(electron_Nocut.at(i), electron_loose_id) && (electron_Nocut.at(i).PFRelMiniIso() < 0.10) && pass_eta_cut){ //&& (electron_Nocut.at(i).MissingHits() == 0)  ){
      if( fabs(electron_Nocut.at(i).SCEta()) < 0.8 && electron_Nocut.at(i).MVA() > -0.85  ) electrons.push_back(electron_Nocut.at(i));
      else if( fabs(electron_Nocut.at(i).SCEta()) > 0.8 && fabs(electron_Nocut.at(i).SCEta()) < 1.479 && electron_Nocut.at(i).MVA() > -0.91  ) electrons.push_back(electron_Nocut.at(i));
      else if( fabs(electron_Nocut.at(i).SCEta()) > 1.479 && fabs(electron_Nocut.at(i).SCEta()) < 2.5 && electron_Nocut.at(i).MVA() > -0.83  ) electrons.push_back(electron_Nocut.at(i));
      else continue;
    }
  }

  int N_muon = muons.size();
  int N_electron = electrons.size();
  
  /*
  if(mu50_pass){
    //FillCLHist(hnpairmm, "Jet_study", eventbase->GetEvent(), muons_veto, electrons_veto, jets_nolepveto, weight, 0);
  }
  */
  //cout << "---------------" << endl;
  //cout << "weight : " << weight << endl;
  
  Signal_region_1("DiMu", mu50_pass, jets_lepveto, fatjets, electrons, muons, N_electron, N_veto_ele, N_muon, N_veto_muon, false);
  Signal_region_1("EMu", mu50_pass, jets_lepveto, fatjets, electrons, muons, N_electron, N_veto_ele, N_muon, N_veto_muon, false);
  Signal_region_1("DiEle", diele_pass, jets_lepveto, fatjets, electrons, muons, N_electron, N_veto_ele, N_muon, N_veto_muon, false);
  
  Control_region_1("DiMu", mu50_pass, jets_lepveto, fatjets, electrons, muons, N_electron, N_veto_ele, N_muon, N_veto_muon, false);
  Control_region_1("EMu", mu50_pass, jets_lepveto, fatjets, electrons, muons, N_electron, N_veto_ele, N_muon, N_veto_muon, false);
  Control_region_1("DiEle", diele_pass, jets_lepveto, fatjets, electrons, muons, N_electron, N_veto_ele, N_muon, N_veto_muon, false);
  
  Control_region_2("DiMu", mu50_pass, jets_lepveto, fatjets, electrons, muons, N_electron, N_veto_ele, N_muon, N_veto_muon, false);
  Control_region_2("EMu", mu50_pass, jets_lepveto, fatjets, electrons, muons, N_electron, N_veto_ele, N_muon, N_veto_muon, false);
  Control_region_2("DiEle", diele_pass, jets_lepveto, fatjets, electrons, muons, N_electron, N_veto_ele, N_muon, N_veto_muon, false);

  Control_region_3("DiMu", mu50_pass, jets_lepveto, fatjets, electrons, muons, N_electron, N_veto_ele, N_muon, N_veto_muon, false);
  Control_region_3("EMu", mu50_pass, jets_lepveto, fatjets, electrons, muons, N_electron, N_veto_ele, N_muon, N_veto_muon, false);
  Control_region_3("DiEle", diele_pass, jets_lepveto, fatjets, electrons, muons, N_electron, N_veto_ele, N_muon, N_veto_muon, false);

  Control_region_4("DiMu", mu50_pass, jets_lepveto, fatjets, electrons, muons, N_electron, N_veto_ele, N_muon, N_veto_muon, false);
  Control_region_4("EMu", mu50_pass, jets_lepveto, fatjets, electrons, muons, N_electron, N_veto_ele, N_muon, N_veto_muon, false);
  Control_region_4("DiEle", diele_pass, jets_lepveto, fatjets, electrons, muons, N_electron, N_veto_ele, N_muon, N_veto_muon, false);
  
  /*
  // -- check with s-ch ID
  Signal_region_1("DiMu", mu50_pass, jets_lepveto, fatjets, electron_sch_tight, muon_sch_tight, electron_sch_tight.size(), N_veto_ele, muon_sch_tight.size(), N_veto_muon, false);
  Signal_region_1("EMu", mu50_pass, jets_lepveto, fatjets, electron_sch_tight, muon_sch_tight, electron_sch_tight.size(), N_veto_ele, muon_sch_tight.size(), N_veto_muon, false);
  Signal_region_1("DiEle", diele_pass, jets_lepveto, fatjets, electron_sch_tight, muon_sch_tight, electron_sch_tight.size(), N_veto_ele, muon_sch_tight.size(), N_veto_muon, false);

  Control_region_1("DiMu", mu50_pass, jets_lepveto, fatjets, electron_sch_tight, muon_sch_tight, electron_sch_tight.size(), N_veto_ele, muon_sch_tight.size(), N_veto_muon, false);
  Control_region_1("EMu", mu50_pass, jets_lepveto, fatjets, electron_sch_tight, muon_sch_tight, electron_sch_tight.size(), N_veto_ele, muon_sch_tight.size(), N_veto_muon, false);
  Control_region_1("DiEle", diele_pass, jets_lepveto, fatjets, electron_sch_tight, muon_sch_tight, electron_sch_tight.size(), N_veto_ele, muon_sch_tight.size(), N_veto_muon, false);

  Control_region_2("DiMu", mu50_pass, jets_lepveto, fatjets, electron_sch_tight, muon_sch_tight, electron_sch_tight.size(), N_veto_ele, muon_sch_tight.size(), N_veto_muon, false);
  Control_region_2("EMu", mu50_pass, jets_lepveto, fatjets, electron_sch_tight, muon_sch_tight, electron_sch_tight.size(), N_veto_ele, muon_sch_tight.size(), N_veto_muon, false);
  Control_region_2("DiEle", diele_pass, jets_lepveto, fatjets, electron_sch_tight, muon_sch_tight, electron_sch_tight.size(), N_veto_ele, muon_sch_tight.size(), N_veto_muon, false);
  

  Control_region_3("DiMu", mu50_pass, jets_lepveto, fatjets, electron_sch_tight, muon_sch_tight, electron_sch_tight.size(), N_veto_ele, muon_sch_tight.size(), N_veto_muon, false);
  Control_region_3("EMu", mu50_pass, jets_lepveto, fatjets, electron_sch_tight, muon_sch_tight, electron_sch_tight.size(), N_veto_ele, muon_sch_tight.size(), N_veto_muon, false);
  Control_region_3("DiEle", diele_pass, jets_lepveto, fatjets, electron_sch_tight, muon_sch_tight, electron_sch_tight.size(), N_veto_ele, muon_sch_tight.size(), N_veto_muon, false);

  Control_region_4("DiMu", mu50_pass, jets_lepveto, fatjets, electron_sch_tight, muon_sch_tight, electron_sch_tight.size(), N_veto_ele, muon_sch_tight.size(), N_veto_muon, false);
  Control_region_4("EMu", mu50_pass, jets_lepveto, fatjets, electron_sch_tight, muon_sch_tight, electron_sch_tight.size(), N_veto_ele, muon_sch_tight.size(), N_veto_muon, false);
  Control_region_4("DiEle", diele_pass, jets_lepveto, fatjets, electron_sch_tight, muon_sch_tight, electron_sch_tight.size(), N_veto_ele, muon_sch_tight.size(), N_veto_muon, false);
  */
  
  return;
}// End of execute event loop

void HN_pair_all::Signal_region_1(TString channel, bool trigger_pass, std::vector<snu::KJet> jets, std::vector<snu::KFatJet> fatjets, std::vector<snu::KElectron> electrons, std::vector<snu::KMuon> muons, int N_electron, int N_veto_ele, int N_muon, int N_veto_muon, bool NonPromptRun){
  // -- channel = DiEle, DiMu, EMu
  bool debug = false;
  
  // -- Trigger pass
  if(!trigger_pass) return;

   
  // -- check N lepton condition.
  bool pass_N_lep = false;
  std::vector<KLepton> Leptons;
  if(channel.Contains("DiEle")){
    if(N_electron == 2 && N_veto_ele == 2 && N_veto_muon == 0){
      FillHist("signal_eff", 5.5, weight, 0., 40., 40); 
      pass_N_lep =true;
      Leptons.push_back(electrons.at(0));
      Leptons.push_back(electrons.at(1));
    }
  }
  else if(channel.Contains("DiMu")){
    if(N_muon == 2 && N_veto_muon == 2 && N_veto_ele == 0){
      FillHist("signal_eff", 6.5, weight, 0., 40., 40);
      pass_N_lep =true;
      Leptons.push_back(muons.at(0));
      Leptons.push_back(muons.at(1));
    }
  }
  else if(channel.Contains("EMu")){
    if(N_muon == 1 && N_veto_muon ==1 && N_electron == 1 && N_veto_ele == 1){
      pass_N_lep =true;
      FillHist("signal_eff", 7.5, weight, 0., 40., 40);
      if(electrons.at(0).Pt() > muons.at(0).Pt()){
	Leptons.push_back(electrons.at(0));
	Leptons.push_back(muons.at(0));
      }
      else{
	Leptons.push_back(muons.at(0));
	Leptons.push_back(electrons.at(0));
      }
    }
  }
  else pass_N_lep = false;
 
  if(!pass_N_lep) return;
  
  //if(debug) cout << "1 , trigger pass : " << trigger_pass << endl;
    
  // -- return ! Pt(1st muon) > 60 GeV && Pt(2nd muon) > 20 GeV
  double Lep_1st_Pt, Lep_2nd_Pt;
  Lep_1st_Pt = Leptons.at(0).Pt();
  Lep_2nd_Pt = Leptons.at(1).Pt();

  if(Lep_1st_Pt > 53 && Lep_2nd_Pt > 20){
    if(channel.Contains("DiMu")) FillHist("signal_eff", 29.5, weight, 0., 40., 40);
  }
  if(Lep_1st_Pt > 53 && Lep_2nd_Pt > 53){
    if(channel.Contains("DiMu")) FillHist("signal_eff", 30.5, weight, 0., 40., 40);
  }
  if(Lep_1st_Pt > 65 && Lep_2nd_Pt > 53){
    if(channel.Contains("DiMu")) FillHist("signal_eff", 31.5, weight, 0., 40., 40);
  }

  if(Lep_1st_Pt < 65 || Lep_2nd_Pt< 65) return;
  if(channel.Contains("DiEle")){
    FillHist("signal_eff", 8.5, weight, 0., 40., 40);
  }
  else if(channel.Contains("DiMu")){
    FillHist("signal_eff", 9.5, weight, 0., 40., 40);
  }
  else if(channel.Contains("EMu")){
    FillHist("signal_eff", 10.5, weight, 0., 40., 40);
  }
  else return;


  if(debug) cout << "2"<< endl;

 
  std::vector<snu::KJet> jets_pt40;
  for(int i = 0; i < jets.size(); i++){
    if(jets.at(i).Pt() > 40) jets_pt40.push_back(jets.at(i));
  }

  int nbjet=0;
  for(unsigned int ij=0; ij < jets_pt40.size(); ij++){
    if( IsBTagged(jets_pt40.at(ij),   snu::KJet::CSVv2, snu::KJet::Medium, -1, 0))  nbjet++;
  }
  if(nbjet != 0) return;


  snu::KParticle ll = Leptons.at(0) + Leptons.at(1);
  if(ll.M() < 150) return;
  if(channel.Contains("DiEle")){
    FillHist("signal_eff",11.5, weight, 0., 40., 40);
  }
  else if(channel.Contains("DiMu")){
    FillHist("signal_eff", 12.5, weight, 0., 40., 40);
  }
  else if(channel.Contains("EMu")){
    FillHist("signal_eff", 13.5, weight, 0., 40., 40);
  }
  else return;

 
  snu::KParticle lljjjj;

  // -- Check jet condition
  if(fatjets.size() == 0 && jets_pt40.size() > 3){
    lljjjj = Leptons.at(0) + Leptons.at(1) + jets_pt40.at(0) + jets_pt40.at(1) + jets_pt40.at(2) + jets_pt40.at(3);
  }
  else if(fatjets.size() == 1 && jets_pt40.size() > 1){
    bool lep_in_ak8 = false;

    double dR_fat_lep1 = fatjets.at(0).DeltaR(Leptons.at(0));
    double dR_fat_lep2 = fatjets.at(0).DeltaR(Leptons.at(1));

    KLepton closest_lep;
    KLepton further_lep;
    
    if(dR_fat_lep1 < dR_fat_lep2){
      closest_lep = Leptons.at(0);
      further_lep = Leptons.at(1);
    }
    else{
      closest_lep = Leptons.at(1);
      further_lep = Leptons.at(0);
    }

    if(fatjets.at(0).DeltaR(closest_lep) < 0.8) lep_in_ak8 = true;
    
    if(lep_in_ak8) lljjjj = fatjets.at(0) + further_lep + jets_pt40.at(0) + jets_pt40.at(1);
    else lljjjj = closest_lep + further_lep + fatjets.at(0) + jets_pt40.at(0) + jets_pt40.at(1);
  }
  else if(fatjets.size() > 1){
    bool lep1_in_ak8 = false;
    bool lep2_in_ak8 = false;
    
    double dR_fat_lep1 = fatjets.at(0).DeltaR(Leptons.at(0));
    double dR_fat_lep2 = fatjets.at(0).DeltaR(Leptons.at(0));

    KLepton closest_lep;
    KLepton further_lep;
    if(dR_fat_lep1 < dR_fat_lep2){
      closest_lep = Leptons.at(0);
      further_lep = Leptons.at(0);
    }
    else{
      closest_lep = Leptons.at(0);
      further_lep = Leptons.at(0);
    }

    if(fatjets.at(0).DeltaR(closest_lep) < 0.8) lep1_in_ak8 = true;
    if(fatjets.at(1).DeltaR(further_lep) < 0.8) lep2_in_ak8 = true;
    
    if(lep1_in_ak8 && lep2_in_ak8) lljjjj = fatjets.at(0) + fatjets.at(1);
    else if(lep1_in_ak8 && !lep2_in_ak8) lljjjj = fatjets.at(0) + fatjets.at(1) + further_lep;
    else if(!lep1_in_ak8 && lep2_in_ak8) lljjjj = fatjets.at(0) + fatjets.at(1) + closest_lep;
    else lljjjj = fatjets.at(0) + fatjets.at(1) + further_lep + closest_lep;
  }
  else return;
  if(channel.Contains("DiEle")){
    FillHist("signal_eff", 14.5, weight, 0., 40., 40);
  }
  else if(channel.Contains("DiMu")){
    FillHist("signal_eff", 15.5, weight, 0., 40., 40);
  }
  else if(channel.Contains("EMu")){
    FillHist("signal_eff", 16.5, weight, 0., 40., 40);
  }
  else return;


  if(lljjjj.M() <300) return;
  if(channel.Contains("DiEle")){
    FillHist("signal_eff", 17.5, weight, 0., 40., 40);
  }
  else if(channel.Contains("DiMu")){
    FillHist("signal_eff", 18.5, weight, 0., 40., 40);
  }
  else if(channel.Contains("EMu")){
    FillHist("signal_eff", 19.5, weight, 0., 40., 40);
  }
  else return;
  
  
  // -- So far, we have lepton, jets, and muons as jets_pt30, muons 

  // -- Call additional Weights & SFs
  //double trigger_eff_Data = mcdata_correction->TriggerEfficiencyLegByLeg(electrons_veto, "", muons, "MUON_HN_TIGHT", 0, 0, 0);
  //double trigger_eff_MC = mcdata_correction->TriggerEfficiencyLegByLeg(electrons_veto, "", muons, "MUON_HN_TIGHT", 0, 1, 0);
  //double trigger_sf = trigger_eff_Data/trigger_eff_MC;
  double trigger_sf = 1.; // 1st muon Pt > 60 GeV with HLT_Mu50 would lead to trigger SF ~ 1. Approx. for now
  
  double muon_id_iso_sf = mcdata_correction->MuonScaleFactor("MUON_SUSY_TIGHT", muons, 0);
  double MuTrkEffSF = mcdata_correction->MuonTrackingEffScaleFactor(muons);
  //cout << ", muon_1_pt : " << muons.at(0).Pt() << ", muon_2_pt : " << muons.at(1).Pt() << ", id_sf : " << muon_id_iso_sf << endl;
  double electron_sf = mcdata_correction->ElectronScaleFactor("ELECTRON_SUSY_TIGHT", electrons, 0);
  double electron_RecoSF = mcdata_correction->ElectronRecoScaleFactor(electrons);
  

  cout << "N_electron : " << N_electron << ", electron_sf : " << electron_sf << ", N_muon : " << N_muon << ", muon_id_iso_sf : " << muon_id_iso_sf << endl;


  double current_weight = weight;
  if(!isData){
    current_weight = current_weight * trigger_sf * muon_id_iso_sf * MuTrkEffSF * electron_sf * electron_RecoSF;
  }
  
  if(NonPromptRun){
    double FR_weight = GetFRWeight_SB(muons, "MUON_SUSY_TIGHT");
    current_weight = FR_weight;
  }
  
  FillCLHist(hnpairmm, "SR1_" + channel, eventbase->GetEvent(), Leptons, N_electron, N_muon, jets_pt40, fatjets, current_weight, 0); 
  
}

void HN_pair_all::Control_region_1(TString channel, bool trigger_pass, std::vector<snu::KJet> jets, std::vector<snu::KFatJet> fatjets, std::vector<snu::KElectron> electrons, std::vector<snu::KMuon> muons, int N_electron, int N_veto_ele, int N_muon, int N_veto_muon, bool NonPromptRun){
  // -- This functio covers two control regions
  // -- 1. M(ll) < 130 GeV without M(lljjjj) cut
  // -- 2. M(lljjjj) < 300 GeV

  bool debug = false;
  
  // -- Trigger pass
  if(!trigger_pass) return;

    
  //check N lepton condition
  bool pass_N_lep = false;
  std::vector<KLepton> Leptons;
  if(channel.Contains("DiEle")){
    if(N_electron == 2 && N_veto_ele == 2 && N_veto_muon == 0){
      pass_N_lep =true;
      Leptons.push_back(electrons.at(0));
      Leptons.push_back(electrons.at(1));
    }
  }
  else if(channel.Contains("DiMu")){
    if(N_muon == 2 && N_veto_muon == 2 && N_veto_ele == 0){
      pass_N_lep =true;
      Leptons.push_back(muons.at(0));
      Leptons.push_back(muons.at(1));
    }
  }
  else if(channel.Contains("EMu")){
    if(N_muon == 1 && N_veto_muon ==1 && N_electron == 1 && N_veto_ele == 1){
      pass_N_lep =true;
      if(electrons.at(0).Pt() > muons.at(0).Pt()){
	Leptons.push_back(electrons.at(0));
        Leptons.push_back(muons.at(0));
      }
      else{
	Leptons.push_back(muons.at(0));
	Leptons.push_back(electrons.at(0));
      }
    }
  }
  else pass_N_lep = false;

  if(!pass_N_lep) return;

  
  //if(debug) cout << "1 , trigger pass : " << trigger_pass << endl;
  // -- return ! Pt(1st muon) > 60 GeV && Pt(2nd muon) > 20 GeV
  double Lep_1st_Pt, Lep_2nd_Pt;
  Lep_1st_Pt = Leptons.at(0).Pt();
  Lep_2nd_Pt = Leptons.at(1).Pt();
  if(Lep_1st_Pt < 65 || Lep_2nd_Pt< 65) return;

  if(debug) cout << "2"<< endl;


  std::vector<snu::KJet> jets_pt40;
  for(int i = 0; i < jets.size(); i++){
    if(jets.at(i).Pt() > 40) jets_pt40.push_back(jets.at(i));
  }

  int nbjet=0;
  for(unsigned int ij=0; ij <jets_pt40.size(); ij++){
    if( IsBTagged(jets_pt40.at(ij),   snu::KJet::CSVv2, snu::KJet::Medium, -1, 0))  nbjet++;
  }
  if(nbjet == 0) return;

  snu::KParticle ll = Leptons.at(0) + Leptons.at(1);
  snu::KParticle lljjjj;

  if(fatjets.size() == 0 && jets_pt40.size() > 3){
    lljjjj = Leptons.at(0) + Leptons.at(1) + jets_pt40.at(0) + jets_pt40.at(1) + jets_pt40.at(2) + jets_pt40.at(3);
  }
  else if(fatjets.size() == 1 && jets_pt40.size() > 1){
    bool lep_in_ak8 = false;

    double dR_fat_lep1 = fatjets.at(0).DeltaR(Leptons.at(0));
    double dR_fat_lep2 = fatjets.at(0).DeltaR(Leptons.at(1));

    KLepton closest_lep;
    KLepton further_lep;

    if(dR_fat_lep1 < dR_fat_lep2){
      closest_lep = Leptons.at(0);
      further_lep = Leptons.at(1);
    }
    else{
      closest_lep = Leptons.at(1);
      further_lep = Leptons.at(0);
    }

    if(fatjets.at(0).DeltaR(closest_lep) < 0.8) lep_in_ak8 = true;

    if(lep_in_ak8) lljjjj = fatjets.at(0) + further_lep + jets_pt40.at(0) + jets_pt40.at(1);
    else lljjjj = closest_lep + further_lep + fatjets.at(0) + jets_pt40.at(0) + jets_pt40.at(1);
  }
  else if(fatjets.size() > 1){
    
    bool lep1_in_ak8 = false;
    bool lep2_in_ak8 = false;
    
    double dR_fat_lep1 = fatjets.at(0).DeltaR(Leptons.at(0));
    double dR_fat_lep2 = fatjets.at(0).DeltaR(Leptons.at(0));

    KLepton closest_lep;
    KLepton further_lep;
    if(dR_fat_lep1 < dR_fat_lep2){
      closest_lep = Leptons.at(0);
      further_lep = Leptons.at(0);
    }
    else{
      closest_lep = Leptons.at(0);
    further_lep = Leptons.at(0);
    }
    
    if(fatjets.at(0).DeltaR(closest_lep) < 0.8) lep1_in_ak8 = true;
    if(fatjets.at(1).DeltaR(further_lep) < 0.8) lep2_in_ak8 = true;
    
    if(lep1_in_ak8 && lep2_in_ak8) lljjjj = fatjets.at(0) + fatjets.at(1);
    else if(lep1_in_ak8 && !lep2_in_ak8) lljjjj = fatjets.at(0) + fatjets.at(1) + further_lep;
    else if(!lep1_in_ak8 && lep2_in_ak8) lljjjj = fatjets.at(0) + fatjets.at(1) + closest_lep;
    else lljjjj = fatjets.at(0) + fatjets.at(1) + further_lep + closest_lep;
  }
  else return;
  

  bool CR1 = false, CR2 = false;
  TString which_CR;
  if(ll.M() < 150){
    CR1 = true;
    which_CR = "CR1_";
    if(channel.Contains("DiEle")){
      FillHist("signal_eff", 20.5, weight, 0., 40., 40);
    }
    else if(channel.Contains("DiMu")){
      FillHist("signal_eff", 21.5, weight, 0., 40., 40);
    }
    else if(channel.Contains("EMu")){
      FillHist("signal_eff", 22.5, weight, 0., 40., 40);
    }
    else return;
  }
  else if(ll.M() > 150 && lljjjj.M() < 300 ){
    CR2 = true;
    which_CR = "CR2_";
   
  }
  else return;
  // -- So far, we have jets, and muons as jets_pt30, muons
  // -- Call additional Weights & SFs 
  //double trigger_eff_Data = mcdata_correction->TriggerEfficiencyLegByLeg(electrons_veto, "", muons, "MUON_HN_TIGHT", 0, 0, 0);
  //double trigger_eff_MC = mcdata_correction->TriggerEfficiencyLegByLeg(electrons_veto, "", muons, "MUON_HN_TIGHT", 0, 1, 0);
  //double trigger_sf = trigger_eff_Data/trigger_eff_MC;
  double trigger_sf = 1.; // 1st muon Pt > 60 GeV with HLT_Mu50 would lead to trigger SF ~ 1. Approx. for now

  double muon_id_iso_sf = mcdata_correction->MuonScaleFactor("MUON_SUSY_TIGHT", muons, 0);
  double MuTrkEffSF = mcdata_correction->MuonTrackingEffScaleFactor(muons);
  //cout << muon_tight_id << ", muon_1_pt : " << muons.at(0).Pt() << ", muon_2_pt : " << muons.at(1).Pt() << ", id_sf : " << muon_id_iso_sf << endl;
  double electron_sf = mcdata_correction->ElectronScaleFactor("ELECTRON_SUSY_TIGHT", electrons, 0);
  double electron_RecoSF = mcdata_correction->ElectronRecoScaleFactor(electrons);

  double current_weight = weight;
  if(!isData){
    current_weight = current_weight * trigger_sf * muon_id_iso_sf * MuTrkEffSF * electron_sf * electron_RecoSF;
  }

  if(NonPromptRun){
    double FR_weight = GetFRWeight_SB(muons, "MUON_SUSY_TIGHT");
    current_weight = FR_weight;
  }

  FillCLHist(hnpairmm, which_CR + channel, eventbase->GetEvent(), Leptons, N_electron, N_muon, jets_pt40, fatjets, current_weight, 0);
}

void HN_pair_all::Control_region_2(TString channel, bool trigger_pass, std::vector<snu::KJet> jets, std::vector<snu::KFatJet> fatjets, std::vector<snu::KElectron> electrons, std::vector<snu::KMuon> muons, int N_electron, int N_veto_ele, int N_muon, int N_veto_muon, bool NonPromptRun){
  // -- 1. No mass cut for Njet = 0 or 1
  // -- To check overall normalization
  bool debug = false;

  if(!trigger_pass) return;

  //if(jets.size() > 1) return;

  //check N lepton condition
  bool pass_N_lep = false;
  std::vector<KLepton> Leptons;
  if(channel.Contains("DiEle")){
    if(N_electron == 2 && N_veto_ele == 2 && N_veto_muon == 0){
      pass_N_lep =true;
      Leptons.push_back(electrons.at(0));
      Leptons.push_back(electrons.at(1));
    }
  }
  else if(channel.Contains("DiMu")){
    if(N_muon == 2 && N_veto_muon == 2 && N_veto_ele == 0){
      pass_N_lep =true;
      Leptons.push_back(muons.at(0));
      Leptons.push_back(muons.at(1));
    }
  }
  else if(channel.Contains("EMu")){
    if(N_muon == 1 && N_veto_muon == 1 && N_electron == 1 && N_veto_ele == 1){
      pass_N_lep =true;
      if(electrons.at(0).Pt() > muons.at(0).Pt()){
        Leptons.push_back(electrons.at(0));
        Leptons.push_back(muons.at(0));
      }
      else{
	Leptons.push_back(muons.at(0));
        Leptons.push_back(electrons.at(0));
      }
    }
  }
  else pass_N_lep = false;
  if(!pass_N_lep) return;
  
  double Lep_1st_Pt, Lep_2nd_Pt;
  Lep_1st_Pt = Leptons.at(0).Pt();
  Lep_2nd_Pt = Leptons.at(1).Pt();
  if(Lep_1st_Pt < 65 || Lep_2nd_Pt< 65) return;


  snu::KParticle ll = Leptons.at(0) + Leptons.at(1);
  if(ll.M() < 10) return;
  
  if(debug) cout << "2"<< endl;


  std::vector<snu::KJet> jets_pt40;
  for(int i = 0; i < jets.size(); i++){
    if(jets.at(i).Pt() > 40) jets_pt40.push_back(jets.at(i));
  }

  if(jets_pt40.size() > 1) return;
  
  double trigger_sf = 1.; // 1st muon Pt > 60 GeV with HLT_Mu50 would lead to trigger SF ~ 1. Approx. for now
  double muon_id_iso_sf = mcdata_correction->MuonScaleFactor("MUON_SUSY_TIGHT", muons, 0);
  double MuTrkEffSF = mcdata_correction->MuonTrackingEffScaleFactor(muons);
  //cout << muon_tight_id << ", muon_1_pt : " << muons.at(0).Pt() << ", muon_2_pt : " << muons.at(1).Pt() << ", id_sf : " << muon_id_iso_sf << endl;
  double electron_sf = mcdata_correction->ElectronScaleFactor("ELECTRON_SUSY_TIGHT", electrons, 0);
  double electron_RecoSF = mcdata_correction->ElectronRecoScaleFactor(electrons);

  double current_weight = weight;
  if(!isData){
    current_weight = current_weight * trigger_sf * muon_id_iso_sf * MuTrkEffSF * electron_sf * electron_RecoSF;
  }

  if(NonPromptRun){
    double FR_weight = GetFRWeight_SB(muons, "MUON_SUSY_TIGHT");
    current_weight = FR_weight;
  }

  FillCLHist(hnpairmm, "CR3_" + channel, eventbase->GetEvent(), Leptons, N_electron, N_muon, jets_pt40, fatjets, current_weight, 0);
  

}

void HN_pair_all::Control_region_3(TString channel, bool trigger_pass, std::vector<snu::KJet> jets, std::vector<snu::KFatJet> fatjets, std::vector<snu::KElectron> electrons, std::vector<snu::KMuon> muons, int N_electron, int N_veto_ele, int N_muon, int N_veto_muon, bool NonPromptRun){
  //CR region for Z bkg, |m(Z) - m(ll)| < 10 GeV without additional cuts

  if(!trigger_pass) return;
  
  bool pass_N_lep = false;
  std::vector<KLepton> Leptons;
  if(channel.Contains("DiEle")){
    if(N_electron == 2 && N_veto_ele == 2 && N_veto_muon == 0){
      pass_N_lep =true;
      Leptons.push_back(electrons.at(0));
      Leptons.push_back(electrons.at(1));
    }
  }
  else if(channel.Contains("DiMu")){
    if(N_muon == 2 && N_veto_muon == 2 && N_veto_ele == 0){
      pass_N_lep =true;
      Leptons.push_back(muons.at(0));
      Leptons.push_back(muons.at(1));
    }
  }
  else if(channel.Contains("EMu")){
    if(N_muon == 1 && N_veto_muon == 1 && N_electron == 1 && N_veto_ele == 1){
      pass_N_lep =true;
      if(electrons.at(0).Pt() > muons.at(0).Pt()){
        Leptons.push_back(electrons.at(0));
        Leptons.push_back(muons.at(0));
      }
      else{
        Leptons.push_back(muons.at(0));
        Leptons.push_back(electrons.at(0));
      }
    }
  }
  else pass_N_lep = false;
  if(!pass_N_lep) return;
  
  double Lep_1st_Pt, Lep_2nd_Pt;
  Lep_1st_Pt = Leptons.at(0).Pt();
  Lep_2nd_Pt = Leptons.at(1).Pt();
  if(Lep_1st_Pt < 65 || Lep_2nd_Pt< 65) return;

  snu::KParticle ll = Leptons.at(0) + Leptons.at(1);
  double M_Z = 91.2;
  if(fabs( ll.M() - M_Z ) > 10) return;

  if(channel.Contains("DiEle")){
    FillHist("signal_eff", 23.5, weight, 0., 40., 40);
  }
  else if(channel.Contains("DiMu")){
    FillHist("signal_eff", 24.5, weight, 0., 40., 40);
  }
  else if(channel.Contains("EMu")){
    FillHist("signal_eff", 25.5, weight, 0., 40., 40);
  }
  else return;


  std::vector<snu::KJet> jets_pt40;
  for(int i = 0; i < jets.size(); i++){
    if(jets.at(i).Pt() > 40) jets_pt40.push_back(jets.at(i));
  }

  int nbjet=0;
  for(unsigned int ij=0; ij <jets_pt40.size(); ij++){
    if( IsBTagged(jets_pt40.at(ij),   snu::KJet::CSVv2, snu::KJet::Medium, -1, 0))  nbjet++;
  }
  
  
  double trigger_sf = 1.; // 1st muon Pt > 60 GeV with HLT_Mu50 would lead to trigger SF ~ 1. Approx. for now                                                                                                                                                                   
  double muon_id_iso_sf = mcdata_correction->MuonScaleFactor("MUON_SUSY_TIGHT", muons, 0);
  double MuTrkEffSF = mcdata_correction->MuonTrackingEffScaleFactor(muons);
  double electron_sf = mcdata_correction->ElectronScaleFactor("ELECTRON_SUSY_TIGHT", electrons, 0);
  double electron_RecoSF = mcdata_correction->ElectronRecoScaleFactor(electrons);

  double current_weight = weight;
  if(!isData){
    current_weight = current_weight * trigger_sf * muon_id_iso_sf * MuTrkEffSF * electron_sf * electron_RecoSF;
  }

  if(NonPromptRun){
    double FR_weight = GetFRWeight_SB(muons, "MUON_SUSY_TIGHT");
    current_weight = FR_weight;
  }

  FillCLHist(hnpairmm, "CR4_" + channel, eventbase->GetEvent(), Leptons, N_electron, N_muon, jets_pt40, fatjets, current_weight, nbjet);

}


void HN_pair_all::Control_region_4(TString channel, bool trigger_pass, std::vector<snu::KJet> jets, std::vector<snu::KFatJet> fatjets, std::vector<snu::KElectron> electrons, std::vector<snu::KMuon> muons, int N_electron, int N_veto_ele, int N_muon, int N_veto_muon, bool NonPromptRun){
  //CR of ttbar dominant region. Njet > 1, Nbejt > 0, MEt > 40 GeV
  
  if(!trigger_pass) return;

  bool pass_N_lep = false;
  std::vector<KLepton> Leptons;
  if(channel.Contains("DiEle")){
    if(N_electron == 2 && N_veto_ele == 2 && N_veto_muon == 0){
      pass_N_lep =true;
      Leptons.push_back(electrons.at(0));
      Leptons.push_back(electrons.at(1));
    }
  }
  else if(channel.Contains("DiMu")){
    if(N_muon == 2 && N_veto_muon == 2 && N_veto_ele == 0){
      pass_N_lep =true;
      Leptons.push_back(muons.at(0));
      Leptons.push_back(muons.at(1));
    }
  }
  else if(channel.Contains("EMu")){
    if(N_muon == 1 && N_veto_muon == 1 && N_electron == 1 && N_veto_ele == 1){
      pass_N_lep =true;
      if(electrons.at(0).Pt() > muons.at(0).Pt()){
        Leptons.push_back(electrons.at(0));
        Leptons.push_back(muons.at(0));
      }
      else{
        Leptons.push_back(muons.at(0));
        Leptons.push_back(electrons.at(0));
      }
    }
  }
  else pass_N_lep = false;
  if(!pass_N_lep) return;

  double Lep_1st_Pt, Lep_2nd_Pt;
  Lep_1st_Pt = Leptons.at(0).Pt();
  Lep_2nd_Pt = Leptons.at(1).Pt();
  if(Lep_1st_Pt < 65 || Lep_2nd_Pt< 65) return;
  
  std::vector<snu::KJet> jets_pt40;
  for(int i = 0; i < jets.size(); i++){
    if(jets.at(i).Pt() > 40) jets_pt40.push_back(jets.at(i));
  }
    
  if(jets_pt40.size() < 2) return;
  
  int nbjet=0;
  for(unsigned int ij=0; ij < jets_pt40.size(); ij++){
    if( IsBTagged(jets_pt40.at(ij),   snu::KJet::CSVv2, snu::KJet::Medium, -1, 0))  nbjet++;
  }
  if(nbjet < 1) return;

  float MET = eventbase->GetEvent().PFMET();
  if(MET < 40) return;

  // -- fill cutflow
  if(channel.Contains("DiEle")){
    FillHist("signal_eff", 26.5, weight, 0., 40., 40);
  }
  else if(channel.Contains("DiMu")){
    FillHist("signal_eff", 27.5, weight, 0., 40., 40);
  }
  else if(channel.Contains("EMu")){
    FillHist("signal_eff", 28.5, weight, 0., 40., 40);
  }
  else return;



  double trigger_sf = 1.; // 1st muon Pt > 60 GeV with HLT_Mu50 would lead to trigger SF ~ 1. Approx. for now
  double muon_id_iso_sf = mcdata_correction->MuonScaleFactor("MUON_SUSY_TIGHT", muons, 0);
  double MuTrkEffSF = mcdata_correction->MuonTrackingEffScaleFactor(muons);
  double electron_sf = mcdata_correction->ElectronScaleFactor("ELECTRON_SUSY_TIGHT", electrons, 0);
  double electron_RecoSF = mcdata_correction->ElectronRecoScaleFactor(electrons);
  

  double current_weight = weight;
  if(!isData){
    current_weight = current_weight * trigger_sf * muon_id_iso_sf * MuTrkEffSF * electron_sf * electron_RecoSF;
  }

  if(NonPromptRun){
    double FR_weight = GetFRWeight_SB(muons, "MUON_SUSY_TIGHT");
    current_weight = FR_weight;
  }

  FillCLHist(hnpairmm, "CR5_" + channel, eventbase->GetEvent(), Leptons, N_electron, N_muon, jets_pt40, fatjets, current_weight, nbjet);

}



void HN_pair_all::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void HN_pair_all::BeginCycle() throw( LQError ){
  
  Message("In begin Cycle", INFO);
  
  string analysisdir = getenv("FILEDIR");  
  if(!k_isdata) //reweightPU = new Reweight((analysisdir + "SNUCAT_Pileup.root").c_str());

    //
    //If you wish to output variables to output file use DeclareVariable
    // clear these variables in ::ClearOutputVectors function
    //DeclareVariable(obj, label, treename );
    //DeclareVariable(obj, label ); //-> will use default treename: LQTree
    //  DeclareVariable(out_electrons, "Signal_Electrons", "LQTree");
    //  DeclareVariable(out_muons, "Signal_Muons");

  
    return;
  
}


HN_pair_all::~HN_pair_all() {
  
  Message("In HN_pair_all Destructor" , INFO);
  //if(!k_isdata)delete reweightPU;
  
}



void HN_pair_all::FillCutFlow(TString cut, float weight){

  
  if(GetHist("cutflow")) {
    GetHist("cutflow")->Fill(cut,weight);
   
  }
  else{
    AnalyzerCore::MakeHistograms("cutflow", 5,0.,5.);

    GetHist("cutflow")->GetXaxis()->SetBinLabel(1,"NoCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(2,"EventCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(3,"TriggerCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(4,"VertexCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(5,"DiMu_tight");
   
    
  }
}


void HN_pair_all::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


void HN_pair_all::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this HN_pair_allCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void HN_pair_all::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}


bool HN_pair_all::JSFatJetID(snu::KFatJet fatjet){
  
  if( !(fatjet.SoftDropMass() > 60) ) return false;
  //if( !( fatjet.Tau2()/fatjet.Tau1() < 0.60 ) ) return false;

  return true;

}

bool HN_pair_all::IsAwayFromFatJet(snu::KJet jet, vector<snu::KFatJet> fatjets){

  for(unsigned int i=0; i<fatjets.size(); i++){
    if( jet.DeltaR( fatjets.at(i) ) < 1.0 ) return false;
  }

  return true;

}

double HN_pair_all::GetFRWeight_SB(std::vector<snu::KMuon> muon_looseColl, TString tight_ID){
  
  double mu_1_cone_pt = muon_looseColl.at(0).Pt()*(1 + max(0.,(muon_looseColl.at(0).RelIso04()-0.07) ) );
  double mu_2_cone_pt = muon_looseColl.at(1).Pt()*(1 + max(0.,(muon_looseColl.at(1).RelIso04()-0.07) ) );
  
  double fr1 = 0.;
  double fr2 = 0.;
  
  if(mu_1_cone_pt < 60){
    int binx = hist_Muon_FR_syst["FR_ptcone_central"]->FindBin(mu_1_cone_pt, fabs(muon_looseColl.at(0).Eta()));
    fr1 = float(hist_Muon_FR_syst["FR_ptcone_central"]->GetBinContent(binx));
  }
  else{
    int binx = hist_Muon_FR_syst["FR_ptcone_central"]->FindBin(55., fabs(muon_looseColl.at(0).Eta()));
    fr1 = float(hist_Muon_FR_syst["FR_ptcone_central"]->GetBinContent(binx));
  }
  if(mu_2_cone_pt < 60){
    int binx = hist_Muon_FR_syst["FR_ptcone_central"]->FindBin(mu_2_cone_pt, fabs(muon_looseColl.at(1).Eta()));
    fr2 = float(hist_Muon_FR_syst["FR_ptcone_central"]->GetBinContent(binx));
  }
  else{
    int binx = hist_Muon_FR_syst["FR_ptcone_central"]->FindBin(55., fabs(muon_looseColl.at(1).Eta()));
    fr2 = float(hist_Muon_FR_syst["FR_ptcone_central"]->GetBinContent(binx));
  }

  bool mu_1_tight = PassID( muon_looseColl.at(0), tight_ID);
  bool mu_2_tight = PassID( muon_looseColl.at(1), tight_ID);
  
  double alpha = 0.;
  double termTT = 0.;
  double termTL = 0.;
  double termLT = 0.;
  double termLL = 0.;
  alpha = 1./((1.- fr1)*(1.- fr2));
  termTT = 0;
  termTL = alpha*(fr2*(1.-fr1));
  termLT = alpha*(fr1*(1.-fr2));
  termLL = -1.*alpha*(fr2*(fr1));
  
  if(mu_1_tight && !mu_2_tight){//TL
    return termTL;
  }
  else if(!mu_1_tight&& mu_2_tight){//LT
    return termLT;
  }
  else if(!mu_1_tight&& !mu_2_tight){//LL
    return termLL;
  }
  else return 0.;//TT
}



