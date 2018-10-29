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
  electrons_all_syst.push_back(electron_all);
  electrons_all_syst.push_back(electron_all_Scale_Up);
  electrons_all_syst.push_back(electron_all_Scale_Down);
  
  // == Get Muon selection with syst flags
  std::vector<snu::KMuon> muon_all = GetMuons("MUON_PTETA");
  std::vector<snu::KMuon> muon_all_Scale_Up   = MuonUsePtScale(muon_all,  1);
  std::vector<snu::KMuon> muon_all_Scale_Down = MuonUsePtScale(muon_all, -1);
  
  std::vector<std::vector<snu::KMuon>> muons_all_syst;
  muons_all_syst.clear();
  muons_all_syst.push_back(muon_all);
  muons_all_syst.push_back(muon_all_Scale_Up);
  muons_all_syst.push_back(muon_all_Scale_Down);
  
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
    ExecuteEventFromSyst( electrons_all_syst.at(electron_index[i_syst]), muons_all_syst.at(muon_index[i_syst]), Jets_all_syst.at(jet_index[i_syst]), FatJets_all_syst.at(jet_index[i_syst]), syst_flags[i_syst] );
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
  std::vector<snu::KElectron> electrons_veto = SelectElectronsWithID(electron_all, "Loose", 35., 2.5);
  std::vector<snu::KElectron> electrons = SelectElectronsWithID(electron_all, "Tight", 65., 2.5);
  
  int N_veto_ele = electrons_veto.size();
  int N_electron = electrons.size();

  
  // -- Get Veto Muons, return if there are not exactly two veto muons
  //std::vector<snu::KMuon> muons_susy_veto = GetMuons("MUON_HN_VETO");
  std::vector<snu::KMuon> muons_veto = SelectMuonsWithID(muons_all, "HighPt", 10., 2.4);
  CorrectedMETRochester(muons_veto);
  std::vector<snu::KMuon> muons = SelectMuonsWithID(muons_all, "HighPt_TrkIso", 65., 2.4);
  CorrectedMETRochester(muons);
  
  int N_veto_muon = muons_veto.size();
  int N_muon = muons.size();

  // -- Get AK8 jets
  std::vector<snu::KFatJet> fatjets = SelectFatJetsWithID(fatjets_all, 300., 2.7);
  
  // -- Get AK4 jets
  std::vector<snu::KJet> jets_veto = SelectJetsWithID(jets_all, 40., 5.);
  std::vector<snu::KJet> jets_lepveto;
  
  for(int i_jet = 0; i_jet < jets_veto.size(); i_jet++){
    double NormalJetMaxEta = 2.7;
    
    snu::KJet this_jet = jets_veto.at(i_jet);
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
  double current_weight = weight;
  if(!isData) current_weight *= trigger_weight;

  //cout << "weight : " << weight << ", dimu_pass : " << dimu_pass << endl;
  
  /////////////////////////////////////////////////////////////////
  // -- weight = 1.0      to see signal eff
  /////////////////////////////////////////////////////////////////
 
    
  // -- run Jobs with OS lepton definition, CR3 : Z maas, CR4 : ttbar

  bool trig_pass_for_channel = false;
  TString current_channel;
  if(N_veto_ele == 2 && N_veto_muon == 0){
    current_channel = "DiEle_";
    trig_pass_for_channel = diele_pass;
  }
  else if(N_veto_muon == 2 && N_veto_ele == 0){
    current_channel = "DiMu_";
    trig_pass_for_channel = mu50_pass;
  }
  else if(N_veto_muon == 1 && N_veto_ele == 1){
    current_channel = "EMu_";
    trig_pass_for_channel = mu50_pass;
  }
  else return;

  
  SR(current_channel + syst_flag, trig_pass_for_channel, current_weight, jets_lepveto, fatjets, electrons, electrons_veto, muons, muons_veto, N_electron, N_veto_ele, N_muon, N_veto_muon);
  CR_Z_mass(current_channel + syst_flag, trig_pass_for_channel, current_weight, jets_lepveto, fatjets, electrons, electrons_veto, muons, muons_veto, N_electron, N_veto_ele, N_muon, N_veto_muon);
  CR_ttbar_dom(current_channel + syst_flag, trig_pass_for_channel, current_weight, jets_lepveto, fatjets, electrons, electrons_veto, muons, muons_veto, N_electron, N_veto_ele, N_muon, N_veto_muon);
  
  return;
  
}// End of execute event loop

void HN_pair_all::SR(TString channel, bool trigger_pass, double current_weight,  std::vector<snu::KJet> jets, std::vector<snu::KFatJet> fatjets, std::vector<snu::KElectron> electrons, std::vector<snu::KElectron> electrons_veto, std::vector<snu::KMuon> muons, std::vector<snu::KMuon> muons_veto, 
		     int N_electron, int N_veto_ele, int N_muon, int N_veto_muon){
  // -- channel = DiEle, DiMu, EMu
  //cout << "start Signal_region_1 : " << channel << endl;
  bool debug = false;
  
  // -- Trigger pass
  if(!trigger_pass) return;
  
  
  // -- check N lepton condition.
  bool pass_N_lep = false;
  std::vector<KLepton> Leptons_veto;
  std::vector<KLepton> Leptons;
  Leptons_veto.clear();
  Leptons.clear();
  for(unsigned int i = 0; i < electrons_veto.size(); i++)    Leptons_veto.push_back(electrons_veto.at(i));
  for(unsigned int i = 0; i < muons_veto.size(); i++)    Leptons_veto.push_back(muons_veto.at(i));
  for(unsigned int i = 0; i < electrons.size(); i++)    Leptons.push_back(electrons.at(i));
  for(unsigned int i = 0; i < muons.size(); i++)    Leptons.push_back(muons.at(i));
  
  if(Leptons.size() > 2) return;
  
  
  // -- return ! Pt(1st muon) > 60 GeV && Pt(2nd muon) > 20 GeV
  double Lep_1st_Pt, Lep_2nd_Pt;
  Lep_1st_Pt = Leptons_veto.at(0).Pt();
  Lep_2nd_Pt = Leptons_veto.at(1).Pt();
  
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

 
  snu::KParticle ll = Leptons_veto.at(0) + Leptons_veto.at(1);
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
  
  vector<snu::KParticle> Ns = RecoPairN(Leptons, fatjets, jets);
  if(Ns.size() != 2) return;
    
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
  
  snu::KParticle Zp = Ns.at(0) + Ns.at(1);
  if(Ns.at(0).M() < 80 || Ns.at(1).M() < 80) return;
  if(Zp.M() < 300) return;
    
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
  TString Region_str = "SR_";

  FillLeptonPlots(Leptons, Region_str + channel, weight);
  JSFillHist(Region_str + channel, "mZp_" + Region_str + channel, Zp.M(), weight, 6000, 0., 6000.);
  JSFillHist(Region_str + channel, "mN_" + Region_str + channel, Ns.at(0).M(), weight, 5000, 0., 5000.);
  JSFillHist(Region_str + channel, "mN_" + Region_str + channel, Ns.at(1).M(), weight, 5000, 0., 5000.);
  
}


void HN_pair_all::CR_Z_mass(TString channel, bool trigger_pass, double current_weight, std::vector<snu::KJet> jets, std::vector<snu::KFatJet> fatjets, std::vector<snu::KElectron> electrons, std::vector<snu::KElectron> electrons_veto, std::vector<snu::KMuon> muons, std::vector<snu::KMuon> muons_veto,
			    int N_electron, int N_veto_ele, int N_muon, int N_veto_muon){
  //CR region for Z bkg, |m(Z) - m(ll)| < 10 GeV without additional cuts

  if(!trigger_pass) return;
  
  // -- check N lepton condition.
  bool pass_N_lep = false;
  std::vector<KLepton> Leptons_veto;
  std::vector<KLepton> Leptons;
  Leptons_veto.clear();
  Leptons.clear();
  for(unsigned int i = 0; i < electrons_veto.size(); i++)    Leptons_veto.push_back(electrons_veto.at(i));
  for(unsigned int i = 0; i < muons_veto.size(); i++)    Leptons_veto.push_back(muons_veto.at(i));
  for(unsigned int i = 0; i < electrons.size(); i++)    Leptons.push_back(electrons.at(i));
  for(unsigned int i = 0; i < muons.size(); i++)    Leptons.push_back(muons.at(i));
  
  if(Leptons.size() > 2) return;
  
  // -- return ! Pt(1st muon) > 60 GeV && Pt(2nd muon) > 20 GeV
  double Lep_1st_Pt, Lep_2nd_Pt;
  Lep_1st_Pt = Leptons_veto.at(0).Pt();
  Lep_2nd_Pt = Leptons_veto.at(1).Pt();
  
  if(Lep_1st_Pt < 65 || Lep_2nd_Pt< 65) return;
  
  snu::KParticle ll = Leptons_veto.at(0) + Leptons_veto.at(1);
  double M_ll = ll.M();
  if(fabs(ll.M() - M_Z) > 10.) return;
  
  TString Region_str = "CR_Zmass_";
  FillLeptonPlots(Leptons_veto, Region_str + channel, weight);
  JSFillHist(Region_str + channel, "mll_" + Region_str + channel, M_ll, weight, 1000, 0., 1000.);

  if(Leptons.size() == 2){
    snu::KParticle ll_tight = Leptons.at(0) + Leptons.at(1);
    double M_ll_tight = ll_tight.M();
    JSFillHist(Region_str + channel, "mll_tight_" + Region_str + channel, M_ll_tight, weight, 1000, 0., 1000.);
  }

  vector<snu::KParticle> Ns = RecoPairN(Leptons, fatjets, jets);
  if(Ns.size() != 2) return;
    
  snu::KParticle Zp = Ns.at(0) + Ns.at(1);
  JSFillHist(Region_str + channel, "mZp_" + Region_str + channel, Zp.M(), weight, 6000, 0., 6000.);

}


void HN_pair_all::CR_ttbar_dom(TString channel, bool trigger_pass, double current_weight, std::vector<snu::KJet> jets, std::vector<snu::KFatJet> fatjets, std::vector<snu::KElectron> electrons, std::vector<snu::KElectron> electrons_veto, std::vector<snu::KMuon> muons, std::vector<snu::KMuon> muons_veto,
			       int N_electron, int N_veto_ele, int N_muon, int N_veto_muon){
  //CR of ttbar dominant region. Njet > 1, Nbejt > 0, MEt > 40 GeV
  
  if(!trigger_pass) return;

  // -- check N lepton condition.
  bool pass_N_lep = false;
  std::vector<KLepton> Leptons_veto;
  std::vector<KLepton> Leptons;
  Leptons_veto.clear();
  Leptons.clear();
  for(unsigned int i = 0; i < electrons_veto.size(); i++)    Leptons_veto.push_back(electrons_veto.at(i));
  for(unsigned int i = 0; i < muons_veto.size(); i++)    Leptons_veto.push_back(muons_veto.at(i));
  for(unsigned int i = 0; i < electrons.size(); i++)    Leptons.push_back(electrons.at(i));
  for(unsigned int i = 0; i < muons.size(); i++)    Leptons.push_back(muons.at(i));

  if(Leptons.size() > 2) return;
  
  // -- return ! Pt(1st muon) > 65 GeV && Pt(2nd muon) > 65 GeV
  double Lep_1st_Pt, Lep_2nd_Pt;
  Lep_1st_Pt = Leptons_veto.at(0).Pt();
  Lep_2nd_Pt = Leptons_veto.at(1).Pt();
  if(Lep_1st_Pt < 65 || Lep_2nd_Pt< 65) return;

  snu::KParticle ll = Leptons_veto.at(0) + Leptons_veto.at(1);
  double M_ll = ll.M();
  if(ll.M() < 55.) return;
  
  if(jets.size() < 2) return;

  
  int nbjet=0;
  for(unsigned int ij=0; ij < jets.size(); ij++){
    if( IsBTagged(jets.at(ij),   snu::KJet::CSVv2, snu::KJet::Medium, -1, 0))  nbjet++;
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

  TString Region_str = "CR_ttbar_";
  FillLeptonPlots(Leptons_veto, Region_str + channel, weight);
  JSFillHist(Region_str + channel, "mll_" + Region_str + channel, M_ll, weight, 1000, 0., 1000.);
  
  if(Leptons.size() == 2){
    snu::KParticle ll_tight = Leptons.at(0) + Leptons.at(1);
    double M_ll_tight = ll_tight.M();
    JSFillHist(Region_str + channel, "mll_tight_" + Region_str + channel, M_ll_tight, weight, 1000, 0., 1000.);
  }
  
  vector<snu::KParticle> Ns = RecoPairN(Leptons, fatjets, jets);
  if(Ns.size() != 2) return;
  
  snu::KParticle Zp = Ns.at(0) + Ns.at(1);
  JSFillHist(Region_str + channel, "mZp_" + Region_str + channel, Zp.M(), weight, 6000, 0., 6000.);
  
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

    snu::KParticle Dummy_AllJets = jets.at(0)+jets.at(1)+jets.at(2)+jets.at(3);

    double mindM = 999999999;
    snu::KParticle temp_N[2];

    for(int i=1; i<=3; i++){

      snu::KParticle TwoJet[2];
      TwoJet[0] = jets.at(0)+jets.at(i);
      TwoJet[1] = Dummy_AllJets-TwoJet[0];

      snu::KParticle N_00 = lepptrs.at(0)+TwoJet[0];
      snu::KParticle N_11 = lepptrs.at(1)+TwoJet[1];
      if( fabs(N_00.M()-N_11.M()) < mindM ){
        mindM = fabs(N_00.M()-N_11.M());
	temp_N[0] = N_00;
        temp_N[1] = N_11;
      }

      snu::KParticle N_01 = lepptrs.at(0)+TwoJet[1];
      snu::KParticle N_10 = lepptrs.at(1)+TwoJet[0];
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

      snu::KParticle temp_N[2];

      snu::KFatJet fatjet = fatjets.at(0);
      temp_N[0] = fatjet;
      temp_N[1] = lepptrs.at(0) + jets.at(0) + jets.at(1);

      if(fatjet.DeltaR( lepptrs.at(0)) > 1.1) {
        out.push_back(temp_N[0]);
        out.push_back(temp_N[1]);
      }

    }

    if(lepptrs.size() == 2){
      snu::KParticle temp_N[2];

      snu::KFatJet fatjet = fatjets.at(0);
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
      snu::KParticle temp_N[2];

      temp_N[0] = fatjets.at(0);
      temp_N[1] = fatjets.at(1);

      out.push_back(temp_N[0]);
      out.push_back(temp_N[1]);
    }
    if(lepptrs.size() == 1){

      snu::KFatJet fatjet = fatjets.at(0);
      snu::KParticle temp_N[2];
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

      snu::KParticle temp_N[2];

      snu::KFatJet fatjet = fatjets.at(0); // Leading FatJet this time                                                                                                                                                                                                                                                                                                            

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
    
    this_jet *= this_scale;
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

    this_jet *= this_scale;
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
    
    this_jet *= this_scale;
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
    
    this_jet *= this_scale;
    out.push_back(this_jet);
  }
  
  return out;
}

// == Functions to call object passing certain ID , pt cut , and eta cut
std::vector<snu::KElectron> HN_pair_all::SelectElectronsWithID(std::vector<snu::KElectron> electrons, TString id, double ptmin, double fetamax){
  bool use_loose = false;
  if(id.Contains("Loose")) use_loose = true;
  
  std::vector<snu::KElectron> out;
  for(unsigned int i=0; i<electrons.size(); i++){
    snu::KElectron this_electron = electrons.at(i);
    if(!( this_electron.Pt()>ptmin )){
      continue;
    }
    if(!( fabs(this_electron.SCEta())<fetamax )){
      continue;
    }
    if(!use_loose){
      if(!this_electron.PassLoose() && !this_electron.PassHEEP() ){
	continue;
      }
    }
    else{
      if(!this_electron.PassHEEP()){
	continue;
      }
    }
   
    out.push_back(this_electron);
  }
  return out;
  
}

std::vector<snu::KMuon> HN_pair_all::SelectMuonsWithID(std::vector<snu::KMuon> muons, TString id, double ptmin, double fetamax){
  
  std::vector<snu::KMuon> out;
  bool use_trkiso = false;
  if(id.Contains("TrkIso")) use_trkiso = true;
  
  for(unsigned int i=0; i<muons.size(); i++){
    snu::KMuon this_muon= muons.at(i);
    if(!( this_muon.Pt()>ptmin )){
      continue;
    }
    if(!( fabs(this_muon.Eta())<fetamax )){
      continue;
    }
    if(!( this_muon.IsHighPt() )){
      continue;
    }
    if(use_trkiso && this_muon.TrkIso() / this_muon.Pt() > 0.1){;
      continue;
    }
    out.push_back(this_muon);
  }
  return out;
  
}


std::vector<snu::KJet> HN_pair_all::SelectJetsWithID(std::vector<snu::KJet> jets, double ptmin, double fetamax){
  
  std::vector<snu::KJet> out;

  for(unsigned int i=0; i<jets.size(); i++){
    snu::KJet this_jet= jets.at(i);
    if(!( this_jet.Pt()>ptmin )){
      continue;
    }
    if(!( fabs(this_jet.Eta())<fetamax )){
      continue;
    }
    if(!( this_jet.PassTightID() )){
      continue;
    }
    out.push_back(this_jet);
  }
  return out;
  
}
 

std::vector<snu::KFatJet> HN_pair_all::SelectFatJetsWithID(std::vector<snu::KFatJet> fatjets, double ptmin, double fetamax){

  std::vector<snu::KFatJet> out;

  for(unsigned int i=0; i<fatjets.size(); i++){
    snu::KFatJet this_fatjet= fatjets.at(i);
    if(!( this_fatjet.Pt()>ptmin )){
      continue;
    }
    if(!( fabs(this_fatjet.Eta())<fetamax )){
      continue;
    }
    if(!( this_fatjet.PileupJetIDTight() )){
      continue;
    }
    if(this_fatjet.SoftDropMass() < 60.){
      continue;
    }
    out.push_back(this_fatjet);
  }
  return out;

}
















