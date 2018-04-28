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
  MakeCleverHistograms(hnpairmm, "SR1_Emu");
  MakeCleverHistograms(hnpairmm, "CR1_DiMu");
  MakeCleverHistograms(hnpairmm, "CR1_DiEle");
  MakeCleverHistograms(hnpairmm, "CR1_Emu");
  MakeCleverHistograms(hnpairmm, "CR2_DiMu");
  MakeCleverHistograms(hnpairmm, "CR2_DiEle");
  MakeCleverHistograms(hnpairmm, "CR2_Emu");
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
  //if(!isData) weight*=MCweight;
  
  FillHist("signal_eff", 1.5, 1., 0., 10., 10); //after lepton skim
  
  //return;
  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;
  
  FillCutFlow("NoCut", weight);
  
  if(isData) FillHist("Nvtx_nocut_data",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);
  else  FillHist("Nvtx_nocut_mc",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);
  
  ///#### CAT:::PassBasicEventCuts is updated: uses selections as described in https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFilters: If you see this is out of date please comment
  if(!PassMETFilter()) return;     /// Initial event cuts :
  FillCutFlow("EventCut", weight);
  FillHist("signal_eff", 2.5, 1., 0., 10., 10); //after PassMETFilter
  
  
  /// #### CAT::: triggers stored are all HLT_Ele/HLT_DoubleEle/HLT_Mu/HLT_TkMu/HLT_Photon/HLT_DoublePhoton
  if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex
  FillHist("signal_eff", 3.5, 1., 0., 10., 10); //after good primary vtx cut
  
  m_logger << DEBUG << "passedTrigger "<< LQLogger::endmsg;


  // -- Store Trigger Boolean
  TString mu50_trig1 = "HLT_Mu50_v";
  TString mu50_trig2 = "HLT_TkMu50_v";
  vector<TString> mu50_trignames;
  mu50_trignames.push_back(mu50_trig1);
  mu50_trignames.push_back(mu50_trig2);
  
  TString diele_trig1="HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v";
  
  bool mu50_pass = PassTriggerOR(mu50_trignames); // for mumu and emu channel
  bool diele_pass = PassTrigger(diele_trig1); // for ee channel

  if(!mu50_pass && !diele_pass) return;
  //cout << "dimu_pass : " << dimu_pass << ". mu50_pass ; " << mu50_pass << endl;
  //if(mu50_pass) cout << "mu50_pass" << endl;
  
  // -- Get Veto Electrons, and store number of them
  std::vector<snu::KElectron> electrons_veto = GetElectrons("ELECTRON_HN_VETO");
  int N_veto_ele = electrons_veto.size();
  
  // -- Get Veto Muons, return if there are not exactly two veto muons
  std::vector<snu::KMuon> muons_veto = GetMuons("MUON_HN_VETO");
  CorrectedMETRochester(muons_veto);
  int N_veto_muon = muons_veto.size();
  
  
  // -- Get AK4 jets
  std::vector<snu::KJet> jets_eta5_nolepveto_loosest = GetJets("JET_HN_eta5_nolepveto", 30., 5.);
  std::vector<snu::KJet> jets_nolepveto;
  
  for(int i_jet = 0; i_jet < jets_eta5_nolepveto_loosest.size(); i_jet++){
    double NormalJetMaxEta = 2.7;
    
    snu::KJet this_jet = jets_eta5_nolepveto_loosest.at(i_jet);
    bool PassPUID = this_jet.PassPileUpMVA("Loose");
    bool IsNormalJet = fabs( this_jet.Eta() ) < NormalJetMaxEta;
    //bool lepinside = HasLeptonInsideJet(this_jet, muons_veto, electrons_veto);
    
    if(IsNormalJet && PassPUID) jets_nolepveto.push_back( this_jet ); // normal jets without lepton veto & fatjet veto
  }
  

  // -- Get size of jets
  int N_jet = jets_nolepveto.size();
  
  
  // -- call truth particles
  std::vector<snu::KTruth> truthColl;
  eventbase->GetTruthSel()->Selection(truthColl);
  
  
  FillHist("signal_eff", 4.5, 1., 0., 10., 10); //after tirgger
  
  // -- Call MET
  float MET = eventbase->GetEvent().PFMET();

  // -- Get global weights
  if(!isData) weight*=MCweight;
  float pileup_reweight=(1.0);
  if(!k_isdata){
    pileup_reweight = mcdata_correction->CatPileupWeight(eventbase->GetEvent(),0);
  }
  double trigger_weight = WeightByTrigger("HLT_Mu50_v", TargetLumi);//weight with any unprescaled trigger

  //cout << "weight : " << weight << ", trigger_weight : " << trigger_weight << ", trigger_weight_dimu : " << trigger_weight_dimu << endl;
  
  if(!isData) weight *= trigger_weight;

  //cout << "weight : " << weight << ", dimu_pass : " << dimu_pass << endl;
  
  /////////////////////////////////////////////////////////////////
  // -- weight = 1.0      to see signal eff
  /////////////////////////////////////////////////////////////////
  weight *= pileup_reweight; //for Data VS MC comparison
  //weight = pileup_reweight * MCweight; //to see signal eff
  // -- Check if run fake
  bool NonPromptRun = std::find(k_flags.begin(), k_flags.end(), "RunFake") != k_flags.end();
  
  //TString muon_tight_id = "MUON_HN_TIGHT";
  TString muon_loose_id = "MUON_HN_LOOSEv7_SIP3";
  TString muon_tight_id = "MUON_SUSY_TIGHT";
  std::vector<snu::KMuon> muon_Nocut = GetMuons("MUON_NOCUT",true);
  CorrectedMETRochester(muon_Nocut);
  std::vector<snu::KMuon> muons;
  for(int i = 0; i < muon_Nocut.size(); i++){
    if( !NonPromptRun && PassID(muon_Nocut.at(i), muon_tight_id) && (muon_Nocut.at(i).RelMiniIso() < 0.20) ) muons.push_back(muon_Nocut.at(i));
    if( NonPromptRun  && PassID(muon_Nocut.at(i), muon_loose_id) && (muon_Nocut.at(i).RelMiniIso() < 0.20) ) muons.push_back(muon_Nocut.at(i));//store loose muon for fake bkg
  }


  TString electron_tight_id = "ELECTRON_SUSY_TIGHT";
  TString electron_loose_id = "ELECTRON_HN_FAKELOOSEv7";
  std::vector<snu::KElectron> electron_Nocut = GetElectrons("ELECTRON_NOCUT", true);
  std::vector<snu::KElectron> electrons;
  for(int i = 0; i < electron_Nocut.size(); i++){
    if( !NonPromptRun && PassID(electron_Nocut.at(i), electron_tight_id) && (electron_Nocut.at(i).PFRelMiniIso() < 0.10) ){
      if( fabs(electron_Nocut.at(i).Eta()) < 0.8 && electron_Nocut.at(i).MVA() > -0.85  ) electrons.push_back(electron_Nocut.at(i));
      else if( fabs(electron_Nocut.at(i).Eta()) > 0.8 && fabs(electron_Nocut.at(i).Eta()) < 1.479 && electron_Nocut.at(i).MVA() > -0.91  ) electrons.push_back(electron_Nocut.at(i));
      else if( fabs(electron_Nocut.at(i).Eta()) > 1.479 && fabs(electron_Nocut.at(i).Eta()) < 2.5 && electron_Nocut.at(i).MVA() > -0.83  ) electrons.push_back(electron_Nocut.at(i));
      else continue;
    }
    if( NonPromptRun  && PassID(electron_Nocut.at(i), electron_loose_id) && (electron_Nocut.at(i).PFRelMiniIso() < 0.10) ){
      if( fabs(electron_Nocut.at(i).Eta()) < 0.8 && electron_Nocut.at(i).MVA() > -0.85  ) electrons.push_back(electron_Nocut.at(i));
      else if( fabs(electron_Nocut.at(i).Eta()) > 0.8 && fabs(electron_Nocut.at(i).Eta()) < 1.479 && electron_Nocut.at(i).MVA() > -0.91  ) electrons.push_back(electron_Nocut.at(i));
      else if( fabs(electron_Nocut.at(i).Eta()) > 1.479 && fabs(electron_Nocut.at(i).Eta()) < 2.5 && electron_Nocut.at(i).MVA() > -0.83  ) electrons.push_back(electron_Nocut.at(i));
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
  
  Signal_region_1("DiMu", mu50_pass, jets_nolepveto, electrons, muons, N_electron, N_veto_ele, N_muon, N_veto_muon, false);
  Signal_region_1("EMu", mu50_pass, jets_nolepveto, electrons, muons, N_electron, N_veto_ele, N_muon, N_veto_muon, false);
  Signal_region_1("DiEle", diele_pass, jets_nolepveto, electrons, muons, N_electron, N_veto_ele, N_muon, N_veto_muon, false);

  
  Control_region_1("DiMu", mu50_pass, jets_nolepveto, electrons, muons, N_electron, N_veto_ele, N_muon, N_veto_muon, false);
  Control_region_1("EMu", mu50_pass, jets_nolepveto, electrons, muons, N_electron, N_veto_ele, N_muon, N_veto_muon, false);
  Control_region_1("DiEle", diele_pass, jets_nolepveto, electrons, muons, N_electron, N_veto_ele, N_muon, N_veto_muon, false);
  
  return;
}// End of execute event loop

void HN_pair_all::Signal_region_1(TString channel, bool trigger_pass, std::vector<snu::KJet> jets, std::vector<snu::KElectron> electrons, std::vector<snu::KMuon> muons, int N_electron, int N_veto_ele, int N_muon, int N_veto_muon, bool NonPromptRun){
  //channel = DiEle, DiMu, EMu
  bool debug = false;
  
  if(!trigger_pass) return;
  if(jets.size() < 4) return;
  
   //check N lepton condition.
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
  if(Lep_1st_Pt < 60 || Lep_2nd_Pt< 53) return;

  if(debug) cout << "2"<< endl;

 
  // -- get lepton cleaned jet collection and apply pt > 30 cut again
  std::vector<snu::KJet> cleaned_jets;
  cleaned_jets = Remove_lepton_from_jet(jets, Leptons);

  if(debug) cout << "3" << endl;


  std::vector<snu::KJet> jets_pt40;
  for(int i = 0; i < cleaned_jets.size(); i++){
    if(cleaned_jets.at(i).Pt() > 40) jets_pt40.push_back(cleaned_jets.at(i));
  }

  if(jets_pt40.size() <4) return;
  
  snu::KParticle ll = Leptons.at(0) + Leptons.at(1);
  snu::KParticle lljjjj = Leptons.at(0) + Leptons.at(1) + jets_pt40.at(0) + jets_pt40.at(1) + jets_pt40.at(2) + jets_pt40.at(3);
  
  if(ll.M() < 150 || lljjjj.M() <300) return;
  
  // -- So far, we have lepton cleaned jets, and muons as jets_pt30, muons 
  

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
  

  //cout << "N_electron : " << N_electron << ", electron_sf : " << electron_sf << ", N_muon : " << N_muon << ", muon_id_iso_sf : " << muon_id_iso_sf << endl;


  double current_weight = weight;
  if(!isData){
    current_weight = current_weight * trigger_sf * muon_id_iso_sf * MuTrkEffSF * electron_sf * electron_RecoSF;
  }
  
  if(NonPromptRun){
    double FR_weight = GetFRWeight_SB(muons, "MUON_SUSY_TIGHT");
    current_weight = FR_weight;
  }
  
  FillCLHist(hnpairmm, "SR1_" + channel, eventbase->GetEvent(), Leptons, N_electron, N_muon, jets_pt40, current_weight, 0); 
  
}

void HN_pair_all::Control_region_1(TString channel, bool trigger_pass, std::vector<snu::KJet> jets, std::vector<snu::KElectron> electrons, std::vector<snu::KMuon> muons, int N_electron, int N_veto_ele, int N_muon, int N_veto_muon, bool NonPromptRun){
  // -- This functio covers two control regions
  // -- 1. M(ll) < 130 GeV without M(lljjjj) cut
  // -- 2. M(lljjjj) < 300 GeV

  bool debug = false;
  
  if(!trigger_pass) return;

  if(jets.size() < 4) return;

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
  if(Lep_1st_Pt < 60 || Lep_2nd_Pt< 53) return;

  if(debug) cout << "2"<< endl;


  // -- get lepton cleaned jet collection and apply pt > 30 cut again
  std::vector<snu::KJet> cleaned_jets;
  cleaned_jets = Remove_lepton_from_jet(jets, Leptons);

  if(debug) cout << "3" << endl;


  std::vector<snu::KJet> jets_pt40;
  for(int i = 0; i < cleaned_jets.size(); i++){
    if(cleaned_jets.at(i).Pt() > 40) jets_pt40.push_back(cleaned_jets.at(i));
  }

  if(jets_pt40.size() <4) return;

  snu::KParticle ll = Leptons.at(0) + Leptons.at(1);
  snu::KParticle lljjjj = Leptons.at(0) + Leptons.at(1) + jets_pt40.at(0) + jets_pt40.at(1) + jets_pt40.at(2) + jets_pt40.at(3);

  bool CR1 = false, CR2 = false;
  TString which_CR;
  if(ll.M() < 150){
    CR1 = true;
    which_CR = "CR1_";
  }
  else if(ll.M() > 150 && lljjjj.M() < 300 ){
    CR2 = true;
    which_CR = "CR2_";
  }
  else return;
  // -- So far, we have lepton cleaned jets, and muons as jets_pt30, muons
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

  FillCLHist(hnpairmm, which_CR + channel, eventbase->GetEvent(), Leptons, N_electron, N_muon, jets_pt40, current_weight, 0);
}




std::vector<snu::KJet> HN_pair_all::Remove_lepton_from_jet(std::vector<snu::KJet> jets, std::vector<KLepton> Leptons){
  //this function removes muon's 4 momentum from Non-lepton veto jet's 4 momentum if muon is inside of it (dR < 0.4)
  
  if(Leptons.size() != 2) return jets;
  
  std::vector<snu::KJet> cleaned_jets;  
  int index_closest_jet_mu_1, index_closest_jet_mu_2;
  double dR_min_mu1 = 999., dR_min_mu2 = 999.;
  bool mu_1_inside_jet = false;
  bool mu_2_inside_jet = false;
  for(int i = 0; i < jets.size(); i++){
    double current_dR_mu1 = jets.at(i).DeltaR(Leptons.at(0));
    double current_dR_mu2 = jets.at(i).DeltaR(Leptons.at(1));
    if(current_dR_mu1 < dR_min_mu1){
      index_closest_jet_mu_1 = i;
      dR_min_mu1 = current_dR_mu1;
    }
    if(current_dR_mu2 < dR_min_mu2){
      index_closest_jet_mu_2 = i;
      dR_min_mu2 = current_dR_mu2;
    }
  }
  
  if(dR_min_mu1 < 0.4) mu_1_inside_jet = true;
  if(dR_min_mu2 < 0.4) mu_2_inside_jet = true;
  
  bool from_same_jet = (index_closest_jet_mu_1 == index_closest_jet_mu_2);//two leptons are from same AK4jet
  
  if(from_same_jet){//when two muons are from same jet
    for(int i = 0; i < jets.size(); i++){
      if(i == index_closest_jet_mu_1){
	snu::KJet jet_mu1;
	if(mu_1_inside_jet) jet_mu1 = Clean_jet_lepton(jets.at(i), Leptons.at(0));
	else jet_mu1 = jets.at(i);
	
	snu::KJet jet_mu2;
	if(mu_2_inside_jet) jet_mu2 = Clean_jet_lepton(jet_mu1, Leptons.at(1));
	else jet_mu2 = jet_mu1;
	
	cleaned_jets.push_back(jet_mu2);
      }
      else cleaned_jets.push_back(jets.at(i));
    }

    return cleaned_jets;
  }
  else if(mu_1_inside_jet || mu_2_inside_jet){
    for(int i = 0; i < jets.size(); i++){
      if(i == index_closest_jet_mu_1 && mu_1_inside_jet){
	snu::KJet jet_mu1;
	jet_mu1 = Clean_jet_lepton(jets.at(i), Leptons.at(0));
	cleaned_jets.push_back(jet_mu1);
      }
      else if(i == index_closest_jet_mu_2 && mu_2_inside_jet){
	snu::KJet jet_mu2;
	jet_mu2 = Clean_jet_lepton(jets.at(i), Leptons.at(1));
	cleaned_jets.push_back(jet_mu2);
      }
      else cleaned_jets.push_back(jets.at(i));
    }
    
    return cleaned_jets;
  }
  else return jets;

}

snu::KJet HN_pair_all::Clean_jet_lepton(snu::KJet jet, KLepton Lepton){
  double energy = jet.E() - Lepton.E();
  double phi = ( jet.E() * jet.Phi() - Lepton.E() * Lepton.Phi() ) / energy;
  double eta = ( jet.E() * jet.Eta() - Lepton.E() * Lepton.Eta() ) / energy;
  double cos_theta = TMath::CosH(eta);
  double momentum = TMath::Sqrt( energy * energy - jet.M() * jet.M() );//keep mass of jet as its original's one
  double Pt = momentum * cos_theta;
  jet.SetPtEtaPhiE(Pt, eta, phi, energy);
  
  return jet;
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
  
  if( !(fatjet.PrunedMass() > 40) ) return false;
  if( !(fatjet.PrunedMass() < 130) ) return false;
  if( !( fatjet.Tau2()/fatjet.Tau1() < 0.60 ) ) return false;

  return true;

}

bool HN_pair_all::IsAwayFromFatJet(snu::KJet jet, vector<snu::KFatJet> fatjets){

  for(unsigned int i=0; i<fatjets.size(); i++){
    if( jet.DeltaR( fatjets.at(i) ) < 0.8 ) return false;
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



