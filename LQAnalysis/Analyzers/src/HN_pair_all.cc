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
  MakeCleverHistograms(hnpairmm, "SR_OS_DiMu");
  MakeCleverHistograms(hnpairmm, "SR_OS_DiEle");
  MakeCleverHistograms(hnpairmm, "SR_OS_EMu");
  MakeCleverHistograms(hnpairmm, "CR1_DiMu");
  MakeCleverHistograms(hnpairmm, "CR1_DiEle");
  MakeCleverHistograms(hnpairmm, "CR1_EMu");
  MakeCleverHistograms(hnpairmm, "CR_Z_mass_DiMu");
  MakeCleverHistograms(hnpairmm, "CR_Z_mass_DiEle");
  MakeCleverHistograms(hnpairmm, "CR_Z_mass_EMu");
  MakeCleverHistograms(hnpairmm, "CR_ttbar_dom_DiMu");
  MakeCleverHistograms(hnpairmm, "CR_ttbar_dom_DiEle");
  MakeCleverHistograms(hnpairmm, "CR_ttbar_dom_EMu");
        
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


  return;
}

void HN_pair_all::ExecuteEvents()throw( LQError ){
  //ListTriggersAvailable();

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

  //remove flavour mixing for signal samples
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
  
  
  // ==== Collect systematically shifted objects
  
  // == Get Electron selection with syst flags
  std::vector<snu::KElectron> electron_all = GetElectrons("ELECTRON_PTETA");
  std::vector<snu::KElectron> electron_all_Scale_Up   = ElectronUsePtScale(electron_all,  1);
  std::vector<snu::KElectron> electron_all_Scale_Down = ElectronUsePtScale(electron_all, -1);

  std::vector<std::vector<snu::KElectron>> electrons_all_syst;
  electrons_all_syst.clear();
  electrons_all_syst.clear.push_back(electron_all);
  electrons_all_syst.clear.push_back(electron_all_Scale_Up);
  electrons_all_syst.clear.push_back(electron_all_Scale_Down);
  
  // == Get Muon selection with syst flags
  std::vector<snu::KMuon> muon_all = GetMuons("MUON_PTETA");
  std::vector<snu::KMuon> muon_all_Scale_Up   = MuonUsePtScale(muon_all,  1);
  std::vector<snu::KMuon> muon_all_Scale_Down = MuonUsePtScale(muon_all, -1);
  
  std::vector<std::vector<snu::KMuon>> muns_all_syst;
  muns_all_syst.clear();
  muns_all_syst.push_back(muon_all);
  muns_all_syst.push_back(muon_all_Scale_Up);
  muns_all_syst.push_back(muon_all_Scale_Down);
  
  // == Get Jet selection with syst flags
  std::vector<snu::KJet> Jets_all = GetJets("JET_NOCUT", 20., 5.);
  std::vector<snu::KJet> Jets_all_Scale_Up   = JetUsePtScale(Jets_all,  1);
  std::vector<snu::KJet> Jets_all_Scale_Down = JetUsePtScale(Jets_all, -1);
  std::vector<snu::KJet> Jets_all_Res_Up;
  if(!isData) Jets_all_Res_Up= JetUsePtRes(  Jets_all,  1);
  else Jets_all_Res_Up = GetJets("JET_NOCUT", 20., 5.);
  std::vector<snu::KJet> Jets_all_Res_Down;
  if(!isData) Jets_all_Res_Down = JetUsePtRes(  Jets_all, -1);
  else Jets_all_Res_Down = GetJets("JET_NOCUT", 20., 5.);
  
  std::vector<std::vector<snu::KJet>> Jets_all_syst;
  Jets_all_syst.clear();
  Jets_all_syst.push_back(Jets_all);
  Jets_all_syst.push_back(Jets_all_Scale_Up);
  Jets_all_syst.push_back(Jets_all_Scale_Down);
  Jets_all_syst.push_back(Jets_all_Res_Up);
  Jets_all_syst.push_back(Jets_all_Res_Down);
  
  // == Get FatJet selection with syst flags
  std::vector<snu::KFatJet> FatJets_all = GetFatJets("FATJET_NOCUT");
  std::vector<snu::KFatJet> FatJets_all_Scale_Up   = FatJetUsePtScale(FatJets_all,  1);
  std::vector<snu::KFatJet> FatJets_all_Scale_Down = FatJetUsePtScale(FatJets_all, -1);
  std::vector<snu::KFatJet> FatJets_all_Res_Up;
  if(!isData) FatJets_all_Res_Up = FatJetUsePtRes(FatJets_all,  1);
  else FatJets_all_Res_Up = GetFatJets("FATJET_NOCUT");
  std::vector<snu::KFatJet> FatJets_all_Res_Down;
  if(!isData) FatJets_all_Res_Down = FatJetUsePtRes(FatJets_all,  -1);
  else FatJets_all_Res_Down = GetFatJets("FATJET_NOCUT");

  std::vector<std::vector<snu::KFatJet>> FatJets_all_syst;
  FatJets_all_syst.clear();
  FatJets_all_syst.push_back(FatJets_all);
  FatJets_all_syst.push_back(FatJets_all_Scale_Up);
  FatJets_all_syst.push_back(FatJets_all_Scale_Down);
  FatJets_all_syst.push_back(FatJets_all_Res_Up);
  FatJets_all_syst.push_back(FatJets_all_Res_Down);
  
  
  // ==== Define Systematic flag strings and order  
  const int N_systs = 9;
  TString syst_flags[N_systs] = {"central", "ElectronScaleUp", "ElectronScaleDown", "MuonScaleUp", "MuonScaleDown", "JetsScaleUp", "JetsScaleDown", "JetsResUp", "JetsResDown"};
  int electron_index[N_systs] = {0        , 1                , 2                  , 0            , 0              , 0            , 0              , 0          , 0            };
  int muon_index[N_systs]     = {0        , 0                , 0                  , 1            , 2              , 0            , 0              , 0          , 0            };
  int jet_index[N_systs]      = {0        , 0                , 0                  , 0            , 0              , 1            , 2              , 3          , 4            };
  
  for(int i_syst = 0; i_syst < N_systs; i_syst++){
    ExecuteEventFromSyst( electrons_all_syst.at(electron_index[i_syst]), muons_all_syst(muon_index[i_syst]), Jets_all_syst.at(jet_index[i_syst]), FatJets_all_syst.at(jet_index[i_syst]), syst_flags[i_syst] );
  }
  
  
}    


void HN_pair_all::ExecuteEventFromSyst(std::vector<snu::KElectron> electron_all, std::vector<snu::KMuon> muons_all, std::vector<snu::KJet> jets_all, std::vector<snu::KFatJet> fatjets_all, TString syst_flag){
    
    
  
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
  std::vector<snu::KElectron> electrons_susy_veto = GetElectrons("ELECTRON_SUSY_HNPAIR_VETO");
  std::vector<snu::KElectron> electrons_veto;
  for(int i = 0; i < electrons_susy_veto.size(); i++){
    if( electrons_susy_veto.at(i).PFRelMiniIso(false) < 0.60) {
      electrons_veto.push_back(electrons_susy_veto.at(i));
    }
  }
  int N_veto_ele = electrons_veto.size();
    
  // -- Get Veto Muons, return if there are not exactly two veto muons
  //std::vector<snu::KMuon> muons_susy_veto = GetMuons("MUON_HN_VETO");
  std::vector<snu::KMuon> muons_susy_veto = GetMuons("MUON_SUSY_VETO");
  std::vector<snu::KMuon> muons_veto;
  for(int i = 0; i < muons_susy_veto.size(); i++){
    if( muons_susy_veto.at(i).PFRelMiniIsoRho() < 0.60) muons_veto.push_back(muons_susy_veto.at(i));
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
  TString muon_loose_id = "MUON_SUSY_VETO";
  TString muon_tight_id = "MUON_SUSY_TIGHT";
  
  std::vector<snu::KMuon> muon_Nocut = GetMuons("MUON_NOCUT",true);
  std::vector<snu::KMuon> muons;
  std::vector<snu::KMuon> muons_SS;
  for(int i = 0; i < muon_Nocut.size(); i++){
    //if(abs(muon_Nocut.at(i).MotherPdgId()) == 15) cout << muon_Nocut.at(i).MotherPdgId() << endl;
    bool tau_pass = !tau_veto || abs(muon_Nocut.at(i).MotherPdgId()) != 15;
    bool pass_eta_cut = fabs(muon_Nocut.at(i).Eta()) < 1.4 || fabs(muon_Nocut.at(i).Eta()) > 1.6;
    if( !NonPromptRun && PassID(muon_Nocut.at(i), muon_tight_id) && (muon_Nocut.at(i).PFRelMiniIsoRho() < 0.20) && tau_pass){
      if(pass_eta_cut) muons.push_back(muon_Nocut.at(i));
      muons_SS.push_back(muon_Nocut.at(i));
    }
    if( NonPromptRun  && PassID(muon_Nocut.at(i), muon_loose_id) && (muon_Nocut.at(i).PFRelMiniIsoRho() < 0.60) && tau_pass){
      if(pass_eta_cut) muons.push_back(muon_Nocut.at(i));//store loose muon for fake bkg
      muons_SS.push_back(muon_Nocut.at(i));
    }
  }
  
  std::vector<snu::KElectron> electron_sch_tight = GetElectrons("ELECTRON_HN_TIGHT",true);
  TString electron_tight_id = "ELECTRON_SUSY_HNPAIR_TIGHT";
  //TString electron_tight_id = "ELECTRON_HN_TIGHTv4_miniiso";
  TString electron_loose_id = "ELECTRON_SUSY_HNPAIR_VETO";
  std::vector<snu::KElectron> electron_Nocut = GetElectrons("ELECTRON_NOCUT", true);
  std::vector<snu::KElectron> electrons;
  std::vector<snu::KElectron> electrons_SS;
  for(int i = 0; i < electron_Nocut.size(); i++){
    bool tau_pass = !tau_veto || abs(electron_Nocut.at(i).MotherPdgId()) != 15;
    bool pass_eta_cut = fabs(electron_Nocut.at(i).SCEta()) < 1.4 || fabs(electron_Nocut.at(i).SCEta()) > 1.6;
    if( !NonPromptRun && PassID(electron_Nocut.at(i), electron_tight_id) && tau_pass){ //&& (electron_Nocut.at(i).MissingHits() == 0) ){
      if(pass_eta_cut) electrons.push_back(electron_Nocut.at(i));
      electrons_SS.push_back(electron_Nocut.at(i));
    }
    if( NonPromptRun  && PassID(electron_Nocut.at(i), electron_loose_id) && tau_pass){ //&& (electron_Nocut.at(i).MissingHits() == 0)  ){
      if(pass_eta_cut) electrons.push_back(electron_Nocut.at(i));
      electrons_SS.push_back(electron_Nocut.at(i));
    }
  }

  int N_muon = muons.size();
  int N_electron = electrons.size();
  int N_muon_SS = muons_SS.size();
  int N_electron_SS = electrons_SS.size();
  
  
  // -- run Jobs with OS lepton definition, CR3 : Z maas, CR4 : ttbar
  
  SR("OS_DiMu", mu50_pass, jets_lepveto, fatjets, electrons, muons, N_electron, N_veto_ele, N_muon, N_veto_muon, NonPromptRun);
  SR("OS_EMu", mu50_pass, jets_lepveto, fatjets, electrons, muons, N_electron, N_veto_ele, N_muon, N_veto_muon, NonPromptRun);
  SR("OS_DiEle", diele_pass, jets_lepveto, fatjets, electrons, muons, N_electron, N_veto_ele, N_muon, N_veto_muon, NonPromptRun);
  
  Control_region_1("OS_DiMu", mu50_pass, jets_lepveto, fatjets, electrons, muons, N_electron, N_veto_ele, N_muon, N_veto_muon, NonPromptRun);
  Control_region_1("OS_EMu", mu50_pass, jets_lepveto, fatjets, electrons, muons, N_electron, N_veto_ele, N_muon, N_veto_muon, NonPromptRun);
  Control_region_1("OS_DiEle", diele_pass, jets_lepveto, fatjets, electrons, muons, N_electron, N_veto_ele, N_muon, N_veto_muon, NonPromptRun);
  
  CR_Z_mass("OS_DiMu", mu50_pass, jets_lepveto, fatjets, electrons, muons, N_electron, N_veto_ele, N_muon, N_veto_muon, NonPromptRun);
  CR_Z_mass("OS_EMu", mu50_pass, jets_lepveto, fatjets, electrons, muons, N_electron, N_veto_ele, N_muon, N_veto_muon, NonPromptRun);
  CR_Z_mass("OS_DiEle", diele_pass, jets_lepveto, fatjets, electrons, muons, N_electron, N_veto_ele, N_muon, N_veto_muon, NonPromptRun);

  CR_ttbar_dom("OS_DiMu", mu50_pass, jets_lepveto, fatjets, electrons, muons, N_electron, N_veto_ele, N_muon, N_veto_muon, NonPromptRun);
  CR_ttbar_dom("OS_EMu", mu50_pass, jets_lepveto, fatjets, electrons, muons, N_electron, N_veto_ele, N_muon, N_veto_muon, NonPromptRun);
  CR_ttbar_dom("OS_DiEle", diele_pass, jets_lepveto, fatjets, electrons, muons, N_electron, N_veto_ele, N_muon, N_veto_muon, NonPromptRun);

  
  return;
  
  }// End of execute event loop

void HN_pair_all::SR(TString channel, bool trigger_pass, std::vector<snu::KJet> jets, std::vector<snu::KFatJet> fatjets, std::vector<snu::KElectron> electrons, std::vector<snu::KMuon> muons, int N_electron, int N_veto_ele, int N_muon, int N_veto_muon, bool NonPromptRun){
  // -- channel = DiEle, DiMu, EMu
  //cout << "start Signal_region_1 : " << channel << endl;
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
      //cout << "MM Leptons.at(0).LeptonFlavour() : " << Leptons.at(0).LeptonFlavour() << ", Leptons.at(1).LeptonFlavour() : " << Leptons.at(1).LeptonFlavour() << endl;
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
      //cout << "Leptons.at(0).LeptonFlavour() : " << Leptons.at(0).LeptonFlavour() << ", Leptons.at(1).LeptonFlavour() : " << Leptons.at(1).LeptonFlavour() << endl;
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
  //if(Lep_1st_Pt < 75 || Lep_2nd_Pt< 75) return; //for 2017
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
  snu::KParticle HN_1, HN_2;
  snu::KParticle leading_HN, second_HN;
  // -- Check jet condition
  if(fatjets.size() == 0 && jets_pt40.size() > 3){
    float min_del_HNHN = 9999.;
    float min_del_R_mu1j = 9999.;
    float min_del_R_mu2j = 9999.;
    
    for(int j1 = 0; j1 < 4; j1 ++){
      if(Leptons.at(0).DeltaR(jets_pt40.at(j1)) < min_del_R_mu1j) min_del_R_mu1j = Leptons.at(0).DeltaR(jets_pt40.at(j1));
      if(Leptons.at(1).DeltaR(jets_pt40.at(j1)) < min_del_R_mu2j) min_del_R_mu2j = Leptons.at(1).DeltaR(jets_pt40.at(j1));

      for(int j2 = 0; j2 < 4; j2 ++){
	snu::KParticle all_jet = jets_pt40.at(0) + jets_pt40.at(1) + jets_pt40.at(2) + jets_pt40.at(3);
	snu::KParticle current_HN1 = Leptons.at(0) + jets_pt40.at(j1) + jets_pt40.at(j2);
        all_jet = all_jet - jets_pt40.at(j1) - jets_pt40.at(j2);
	snu::KParticle current_HN2 = Leptons.at(1) + all_jet;

        if(fabs(current_HN1.M() - current_HN2.M()) < min_del_HNHN && j1 != j2){
          min_del_HNHN = fabs(current_HN1.M() - current_HN2.M());
          HN_1 = current_HN1;
          HN_2 = current_HN2;
        }
        current_HN1 = Leptons.at(0) + all_jet;
        current_HN2 = Leptons.at(1) + jets_pt40.at(j1) + jets_pt40.at(j2);
        if(fabs(current_HN1.M() - current_HN2.M()) < min_del_HNHN && j1 != j2){
          min_del_HNHN = fabs(current_HN1.M() - current_HN2.M());
          HN_1 = current_HN1;
          HN_2 = current_HN2;
        }
      }
    }
    if(HN_1.Pt() > HN_2.Pt()){
      leading_HN = HN_1;
      second_HN = HN_2;
    }
    if(HN_1.Pt() < HN_2.Pt()){
      leading_HN = HN_2;
      second_HN = HN_1;
    }

    lljjjj = leading_HN + second_HN;
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
    
    HN_1 = further_lep + jets_pt40.at(0) + jets_pt40.at(1);
    
    if(fatjets[0].DeltaR(closest_lep) < 0.8) lep_in_ak8 = true;
    if(lep_in_ak8) HN_2 = fatjets.at(0);
    else HN_2 = fatjets.at(0) + closest_lep;
    
    if(HN_1.Pt() > HN_2.Pt()){
      leading_HN = HN_1;
      second_HN = HN_2;
    }
    if(HN_1.Pt() < HN_2.Pt()){
      leading_HN = HN_2;
      second_HN = HN_1;
    }
    
    lljjjj = leading_HN + second_HN;
  }
  else if(fatjets.size() > 1){
    bool lep1_in_ak8 = false;
    bool lep2_in_ak8 = false;
    
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
    
    if(fatjets.at(0).DeltaR(closest_lep) < 0.8) lep1_in_ak8 = true;
    if(fatjets.at(1).DeltaR(further_lep) < 0.8) lep2_in_ak8 = true;
    
    if(lep1_in_ak8) HN_1 = fatjets.at(0);
    else HN_1 = fatjets.at(0) + closest_lep;
    
    if(lep2_in_ak8) HN_2 = fatjets.at(1);
    else HN_2 = fatjets.at(1) + further_lep;
    
    if(HN_1.Pt() > HN_2.Pt()){
      leading_HN = HN_1;
      second_HN = HN_2;
    }
    if(HN_1.Pt() < HN_2.Pt()){
      leading_HN = HN_2;
      second_HN = HN_1;
    }
    
    lljjjj = leading_HN + second_HN;
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
  

  //cout << "N_electron : " << N_electron << ", electron_sf : " << electron_sf << ", N_muon : " << N_muon << ", muon_id_iso_sf : " << muon_id_iso_sf << endl;


  double current_weight = weight;
  if(!isData){
    current_weight = current_weight * trigger_sf * muon_id_iso_sf * MuTrkEffSF * electron_sf * electron_RecoSF;
  }
  
  if(NonPromptRun){
    current_weight = weight_fake;
  }
  

  
  

  FillCLHist(hnpairmm, "SR1_" + channel, eventbase->GetEvent(), Leptons, N_electron, N_muon, jets_pt40, fatjets, current_weight, 0); 
  //if(leading_HN.M() > 200) FillCLHist(hnpairmm, "SR1_" + channel + "_HNcut", eventbase->GetEvent(), Leptons, N_electron, N_muon, jets_pt40, fatjets, current_weight, 0);

  if(channel.Contains("SS_EMu") || channel.Contains("SS_DiEle")){
    FillCLHist(hnpairmm, "SR1_" + channel + "_CF", eventbase->GetEvent(), Leptons, N_electron, N_muon, jets_pt40, fatjets, weight_CF, 0);
    //if(leading_HN.M() > 200) FillCLHist(hnpairmm, "SR1_" + channel + "_CF_HNcut", eventbase->GetEvent(), Leptons, N_electron, N_muon, jets_pt40, fatjets, weight_CF, 0);
  }
}

void HN_pair_all::Control_region_1(TString channel, bool trigger_pass, std::vector<snu::KJet> jets, std::vector<snu::KFatJet> fatjets, std::vector<snu::KElectron> electrons, std::vector<snu::KMuon> muons, int N_electron, int N_veto_ele, int N_muon, int N_veto_muon, bool NonPromptRun){
  // -- This functio covers two control regions
  // -- 1. M(ll) < 130 GeV without M(lljjjj) cut

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
    current_weight = weight_fake;
  }

  FillCLHist(hnpairmm, which_CR + channel, eventbase->GetEvent(), Leptons, N_electron, N_muon, jets_pt40, fatjets, current_weight, 0);
  if(channel.Contains("SS_EMu") || channel.Contains("SS_DiEle"))  FillCLHist(hnpairmm, which_CR + channel + "_CF", eventbase->GetEvent(), Leptons, N_electron, N_muon, jets_pt40, fatjets, weight_CF, 0);

}


void HN_pair_all::CR_Z_mass(TString channel, bool trigger_pass, std::vector<snu::KJet> jets, std::vector<snu::KFatJet> fatjets, std::vector<snu::KElectron> electrons, std::vector<snu::KMuon> muons, int N_electron, int N_veto_ele, int N_muon, int N_veto_muon, bool NonPromptRun){
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
    current_weight = weight_fake;
  }

  FillCLHist(hnpairmm, "CR3_" + channel, eventbase->GetEvent(), Leptons, N_electron, N_muon, jets_pt40, fatjets, current_weight, nbjet);
  if(channel.Contains("SS_EMu") || channel.Contains("SS_DiEle"))  FillCLHist(hnpairmm, "CR3_" + channel + "_CF", eventbase->GetEvent(), Leptons, N_electron, N_muon, jets_pt40, fatjets, weight_CF, 0);

}


void HN_pair_all::CR_ttbar_dom(TString channel, bool trigger_pass, std::vector<snu::KJet> jets, std::vector<snu::KFatJet> fatjets, std::vector<snu::KElectron> electrons, std::vector<snu::KMuon> muons, int N_electron, int N_veto_ele, int N_muon, int N_veto_muon, bool NonPromptRun){
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
    current_weight = weight_fake;
  }

  FillCLHist(hnpairmm, "CR4_" + channel, eventbase->GetEvent(), Leptons, N_electron, N_muon, jets_pt40, fatjets, current_weight, nbjet);
  if(channel.Contains("SS_EMu") || channel.Contains("SS_DiEle"))  FillCLHist(hnpairmm, "CR4_" + channel + "_CF", eventbase->GetEvent(), Leptons, N_electron, N_muon, jets_pt40, fatjets, weight_CF, 0);
  
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

vector<snu::KParticle> HN_pair_all::RecoPairN(std::vector<KLepton> lepptrs, vector<snu::KFatJet> fatjets, vector<snu::KJet> jets){

  vector<snu::KParticle> out;
  out.clear();

  
  if(fatjets.size()==0 && jets.size()>3 && lepptrs.size() == 2){

    Particle Dummy_AllJets = jets.at(0)+jets.at(1)+jets.at(2)+jets.at(3);

    double mindM = 999999999;
    Particle temp_N[2];

    for(int i=1; i<=3; i++){

      Particle TwoJet[2];
      TwoJet[0] = jets.at(0)+jets.at(i);
      TwoJet[1] = Dummy_AllJets-TwoJet[0];

      Particle N_00 = lepptrs.at(0)+TwoJet[0];
      Particle N_11 = lepptrs.at(1)+TwoJet[1];
      if( fabs(N_00.M()-N_11.M()) < mindM ){
        mindM = fabs(N_00.M()-N_11.M());
	temp_N[0] = N_00;
        temp_N[1] = N_11;
      }

      Particle N_01 = lepptrs.at(0)+TwoJet[1];
      Particle N_10 = lepptrs.at(1)+TwoJet[0];
      if( fabs(N_01.M()-N_10.M()) < mindM ){
	mindM = fabs(N_01.M()-N_10.M());
        temp_N[0] = N_01;
	temp_N[1] = N_10;
      }

    }

    out.push_back(temp_N[0]);
    out.push_back(temp_N[1]);

  }
  else if(fatjets.size()==1 && jets.size() > 1){

    if(lepptrs.size() == 1){

      Particle temp_N[2];

      FatJet fatjet = fatjets.at(0);
      temp_N[0] = fatjet;
      temp_N[1] = lepptrs.at(0) + jets.at(0) + jets.at(1);

      if(fatjet.DeltaR( lepptrs.at(0)) > 1.1) {
        out.push_back(temp_N[0]);
        out.push_back(temp_N[1]);
      }

    }

    if(lepptrs.size() == 2){
      Particle temp_N[2];

      FatJet fatjet = fatjets.at(0);
      if(fatjet.DeltaR( lepptrs.at(0) ) < fatjet.DeltaR( lepptrs.at(1) )){
        temp_N[0] = AddFatJetAndLepton(fatjet, lepptrs.at(0));
        temp_N[1] = lepptrs.at(1) + jets.at(0) + jets.at(1);
      }
      else{
        temp_N[0] = AddFatJetAndLepton(fatjet, lepptrs.at(1));
        temp_N[1] = lepptrs.at(0) + jets.at(0) + jets.at(1);
      }

      out.push_back(temp_N[0]);
      out.push_back(temp_N[1]);

    }


  }
  else if(fatjets.size()>1){

    if(lepptrs.size() == 0){
      Particle temp_N[2];

      temp_N[0] = fatjets.at(0);
      temp_N[1] = fatjets.at(1);

      out.push_back(temp_N[0]);
      out.push_back(temp_N[1]);
    }
    if(lepptrs.size() == 1){

      FatJet fatjet = fatjets.at(0);
      Particle temp_N[2];
      if(fatjet.DeltaR( lepptrs.at(0) ) < fatjet.DeltaR( lepptrs.at(0) )){
	temp_N[0] = AddFatJetAndLepton(fatjet, lepptrs.at(0));
        temp_N[1] = fatjets.at(1);
      }
      else{
        temp_N[0] = fatjet;
        temp_N[1] = AddFatJetAndLepton(fatjets.at(1), lepptrs.at(0));
      }
    }
    if(lepptrs.size() == 2){

      Particle temp_N[2];

      FatJet fatjet = fatjets.at(0); // Leading FatJet this time                                                                                                                                                                                                                                                                                                            

      if(fatjet.DeltaR( lepptrs.at(0) ) < fatjet.DeltaR( lepptrs.at(1) )){
        temp_N[0] = AddFatJetAndLepton(fatjet,        lepptrs.at(0));
        temp_N[1] = AddFatJetAndLepton(fatjets.at(1), lepptrs.at(1));
      }
      else{
        temp_N[0] = AddFatJetAndLepton(fatjets.at(1), lepptrs.at(0));
        temp_N[1] = AddFatJetAndLepton(fatjet,        lepptrs.at(1));
      }

      out.push_back(temp_N[0]);
      out.push_back(temp_N[1]);

    }


  }
  else{

  }

  return out;



}



// ==== Functions to call object with shift for systematics
std::vector<snu::KElectron> HN_pair_all::ElectronUsePtScale(std::vector<snu::KElectron> electrons, int up_down){
  //up_down = -1 down
  //up_down = 1 up

  std::vector<snu::KElectron> out;
  for(unsigned int i=0; i<electrons.size(); i++){
    snu::KElectron this_electron= electrons.at(i);
    if(up_down == -1) this_electron.SetPtEtaPhiM( this_electron.Pt()*this_electron.PtShiftedDown(), this_electron.Eta(), this_electron.Phi(), this_electron.M() );
    else this_electron.SetPtEtaPhiM( this_electron.Pt()*this_electron.PtShiftedUp(), this_electron.Eta(), this_electron.Phi(), this_electron.M() );

    out.push_back(this_electron);
  }

  return out;
}

std::vector<snu::KMuon> HN_pair_all::MuonUsePtScale(std::vector<snu::KMuon> muons, int up_down){
  //up_down = -1 down
  //up_down = 1 up

  std::vector<snu::KMuon> out;
  for(unsigned int i=0; i<muons.size(); i++){
    snu::KMuon this_muon= muons.at(i);
    if(up_down == -1) this_muon.SetPtEtaPhiM( this_muon.Pt()*this_muon.PtShiftedDown(), this_muon.Eta(), this_muon.Phi(), this_muon.M() );
    else this_muon.SetPtEtaPhiM( this_muon.Pt()*this_muon.PtShiftedUp(), this_muon.Eta(), this_muon.Phi(), this_muon.M() );
    
    out.push_back(this_muon);
  }
  
  return out;
}
std::vector<snu::KJet> HN_pair_all::JetUsePtScale(std::vector<snu::KJet> jets, int up_down){
  //up_down = -1 down
  //up_down = 1 up

  std::vector<snu::KJet> out;
  for(unsigned int i=0; i<jets.size(); i++){
    snu::KJet this_jet = jets.at(i);
    
    double this_scale= 1.;
    if(up_down == -1) this_scale = this_jet.ScaledDownEnergy();
    else this_scale = this_jet.ScaledUpEnergy();
    
    this_jet *= this_scaling;
    out.push_back(this_jet);
  }

  return out;
}
std::vector<snu::KJet> HN_pair_all::JetUsePtRes(std::vector<snu::KJet> jets, int up_down){
  //up_down = -1 down
  //up_down = 1 up

  std::vector<snu::KJet> out;
  for(unsigned int i=0; i<jets.size(); i++){
    snu::KJet this_jet = jets.at(i);

    double this_scale= 1.;
    if(up_down == -1) this_scale = this_jet.SmearedResDown()/this_jet.SmearedRes();
    else this_scale = this_jet.SmearedResUp()/this_jet.SmearedRes();

    this_jet *= this_scaling;
    out.push_back(this_jet);
  }

  return out;
}
std::vector<snu::KFatJet> HN_pair_all::FatJetUsePtScale(std::vector<snu::KFatJet> fatjets, int up_down){
  //up_down = -1 down
  //up_down = 1 up

  std::vector<snu::KFatJet> out;
  for(unsigned int i=0; i<fatjets.size(); i++){
    snu::KFatJet this_jet = fatjets.at(i);

    double this_scale= 1.;
    if(up_down == -1) this_scale = this_jet.ScaledDownEnergy();
    else this_scale = this_jet.ScaledUpEnergy();

    this_jet *= this_scaling;
    out.push_back(this_jet);
  }

  return out;
}
std::vector<snu::KFatJet> HN_pair_all::FatJetUsePtRes(std::vector<snu::KFatJet> fatjets, int up_down){
  //up_down = -1 down
  //up_down = 1 up

  std::vector<snu::KFatJet> out;
  for(unsigned int i=0; i<fatjets.size(); i++){
    snu::KFatJet this_jet = fatjets.at(i);

    double this_scale= 1.;
    if(up_down == -1) this_scale = this_jet.SmearedResDown()/this_jet.SmearedRes();
    else this_scale = this_jet.SmearedResUp()/this_jet.SmearedRes();

    this_jet *= this_scaling;
    out.push_back(this_jet);
  }

  return out;
}






















