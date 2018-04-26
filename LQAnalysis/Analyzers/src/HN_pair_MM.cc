// $Id: HN_pair_MM.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQHN_pair_MM Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "HN_pair_MM.h"

//Core includes
#include "Reweight.h"
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (HN_pair_MM);

/**
 *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
 *
 */
HN_pair_MM::HN_pair_MM() :  AnalyzerCore(), out_muons(0)  {
  // To have the correct name in the log:                                                                                                                            
  SetLogName("HN_pair_MM");
  
  Message("In HN_pair_MM constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
  MakeCleverHistograms(hnpairmm, "Jet_study");
  MakeCleverHistograms(hnpairmm, "SR1_MUON_HN_TIGHT");
  MakeCleverHistograms(hnpairmm, "SR1_MUON_SUSY_TIGHT");

}


void HN_pair_MM::InitialiseAnalysis() throw( LQError ) {

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

void HN_pair_MM::ExecuteEvents()throw( LQError ){
  
  /// Apply the gen weight 
  if(!isData) weight*=MCweight;
  
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
  TString dimuon_trigmuon_trig1="HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v";
  TString dimuon_trigmuon_trig2="HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v";
  TString dimuon_trigmuon_trig3="HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v";
  TString dimuon_trigmuon_trig4="HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v";
  vector<TString> trignames;
  trignames.push_back(dimuon_trigmuon_trig1);
  trignames.push_back(dimuon_trigmuon_trig2);
  trignames.push_back(dimuon_trigmuon_trig3);
  trignames.push_back(dimuon_trigmuon_trig4);
  bool dimu_pass = PassTriggerOR(trignames);
  bool mu50_pass = PassTrigger("HLT_Mu50_v");
  //cout << "dimu_pass : " << dimu_pass << ". mu50_pass ; " << mu50_pass << endl;
  //if(mu50_pass) cout << "mu50_pass" << endl;
  
  // -- Get Veto Electrons, return if there is veto electron
  std::vector<snu::KElectron> electrons_veto = GetElectrons("ELECTRON_HN_VETO");
  if(electrons_veto.size() > 0) return;

  // -- Get Veto Muons, return if there are not exactly two veto muons
  std::vector<snu::KMuon> muons_veto = GetMuons("MUON_HN_VETO");
  if(muons_veto.size() != 2) return;
  CorrectedMETRochester(muons_veto);
  
  // -- Get AK8 jets
  std::vector<snu::KFatJet> fatjets_loosest_UnSmeared = GetFatJets("FATJET_HN_nolepveto");
  std::vector<snu::KFatJet> fatjets_loosest = GetCorrectedFatJet(fatjets_loosest_UnSmeared); // Smear both energy and mass
  std::vector<snu::KFatJet> fatjets;
  for(int i_fat = 0; i_fat < fatjets_loosest.size(); i_fat++){
    snu::KFatJet this_jet = fatjets_loosest.at(i_fat);
    if(JSFatJetID(this_jet) && this_jet.Pt() >= 200.) fatjets.push_back(this_jet);
  }

  // -- Get AK4 jets
  std::vector<snu::KJet> jets_eta5_nolepveto_loosest = GetJets("JET_HN_eta5_nolepveto", 30., 5.);
  std::vector<snu::KJet> jets;
  std::vector<snu::KJet> jets_nolepveto;
  std::vector<snu::KJet> jets_InSideFatJet;
  
  for(int i_jet = 0; i_jet < jets_eta5_nolepveto_loosest.size(); i_jet++){
    double NormalJetMaxEta = 2.7;
    
    snu::KJet this_jet = jets_eta5_nolepveto_loosest.at(i_jet);
    bool PassPUID = this_jet.PassPileUpMVA("Loose");
    bool awayfromfatjet = IsAwayFromFatJet(this_jet, fatjets);
    bool IsNormalJet = fabs( this_jet.Eta() ) < NormalJetMaxEta;
    bool lepinside = HasLeptonInsideJet(this_jet, muons_veto, electrons_veto);
    
    if(IsNormalJet && !lepinside && !awayfromfatjet) jets_InSideFatJet.push_back( this_jet ); // jets insdie fatjet
    if(IsNormalJet && !lepinside && PassPUID && awayfromfatjet) jets.push_back( this_jet ); // normal jets
    if(IsNormalJet && PassPUID) jets_nolepveto.push_back( this_jet ); // normal jets without lepton veto & fatjet veto
  }


  // -- Get size of jets
  int N_fatjet = fatjets.size();
  int N_jet = jets.size();
  
  
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
  double trigger_weight = WeightByTrigger("HLT_Mu50_v", TargetLumi);
  if(!isData) weight *= trigger_weight;
  
  /////////////////////////////////////////////////////////////////
  // -- weight = 1.0      to see signal eff
  /////////////////////////////////////////////////////////////////
  weight *= pileup_reweight; //for Data VS MC comparison
  //weight = pileup_reweight * MCweight; //to see signal eff
  // -- Check if run fake
  bool NonPromptRun = std::find(k_flags.begin(), k_flags.end(), "RunFake") != k_flags.end();
  
  TString muon_tight_id = "MUON_HN_TIGHT";
  TString muon_loose_id = "MUON_HN_LOOSEv7_SIP3";

  TString muon_susy_id = "MUON_SUSY_TIGHT";

  
  if(mu50_pass){
    FillCLHist(hnpairmm, "Jet_study", eventbase->GetEvent(), muons_veto, electrons_veto, jets_nolepveto, weight, 0);
  }

  //cout << "---------------" << endl;
  //cout << "weight : " << weight << endl;
  
  Signal_region_1(muon_tight_id, muon_loose_id, mu50_pass, jets_nolepveto, electrons_veto, false);
  Signal_region_1(muon_susy_id, muon_loose_id, mu50_pass, jets_nolepveto, electrons_veto, false);

  Control_region_1(muon_susy_id, muon_loose_id, mu50_pass, jets_nolepveto, electrons_veto, false);
  Control_region_1(muon_tight_id, muon_loose_id, mu50_pass, jets_nolepveto, electrons_veto, false);

  return;
}// End of execute event loop

void HN_pair_MM::Signal_region_1(TString muon_tight_id, TString muon_loose_id, bool trigger_pass, std::vector<snu::KJet> jets, std::vector<snu::KElectron> electrons_veto, bool NonPromptRun){
  
  bool debug = false;
  
  //if(debug) cout << "1 , trigger pass : " << trigger_pass << endl;
  
  if(!trigger_pass) return;
  if(jets.size() < 4) return;
  std::vector<snu::KMuon> muon_Nocut = GetMuons("MUON_NOCUT",true);
  std::vector<snu::KMuon> muons;
  for(int i = 0; i < muon_Nocut.size(); i++){
    if( !NonPromptRun && PassID(muon_Nocut.at(i), muon_tight_id) && (muon_Nocut.at(i).RelMiniIso() < 0.20) ) muons.push_back(muon_Nocut.at(i));
    if( NonPromptRun  && PassID(muon_Nocut.at(i), muon_loose_id) && (muon_Nocut.at(i).RelMiniIso() < 0.20) ) muons.push_back(muon_Nocut.at(i));//store loose muon for fake bkg
  }
  
  // -- return if Nmuon != 2
  int Nmuons;
  Nmuons = muons.size();
  if(Nmuons != 2) return;
    
  // -- return ! Pt(1st muon) > 60 GeV && Pt(2nd muon) > 20 GeV
  double mu_1st_pt, mu_2nd_pt;
  mu_1st_pt = muons.at(0).Pt();
  mu_2nd_pt = muons.at(1).Pt();
  if(mu_1st_pt < 60 || mu_2nd_pt < 20) return;

  if(debug) cout << "2"<< endl;

 
  // -- get lepton cleaned jet collection and apply pt > 30 cut again
  std::vector<snu::KJet> cleaned_jets;
  cleaned_jets = Remove_muon_from_jet(jets, muons);

  if(debug) cout << "3" << endl;


  std::vector<snu::KJet> jets_pt30;
  for(int i = 0; i < cleaned_jets.size(); i++){
    if(cleaned_jets.at(i).Pt() > 30) jets_pt30.push_back(cleaned_jets.at(i));
  }

  if(jets_pt30.size() <4) return;
  
  
  // -- So far, we have lepton cleaned jets, and muons as jets_pt30, muons 
  

  // -- Call additional Weights & SFs
  //double trigger_eff_Data = mcdata_correction->TriggerEfficiencyLegByLeg(electrons_veto, "", muons, "MUON_HN_TIGHT", 0, 0, 0);
  //double trigger_eff_MC = mcdata_correction->TriggerEfficiencyLegByLeg(electrons_veto, "", muons, "MUON_HN_TIGHT", 0, 1, 0);
  //double trigger_sf = trigger_eff_Data/trigger_eff_MC;
  double trigger_sf = 1.; // 1st muon Pt > 60 GeV with HLT_Mu50 would lead to trigger SF ~ 1. Approx. for now
  
  double muon_id_iso_sf = mcdata_correction->MuonScaleFactor(muon_tight_id, muons, 0);
  double MuTrkEffSF = mcdata_correction->MuonTrackingEffScaleFactor(muons);
  //cout << muon_tight_id << ", muon_1_pt : " << muons.at(0).Pt() << ", muon_2_pt : " << muons.at(1).Pt() << ", id_sf : " << muon_id_iso_sf << endl;
  
  double current_weight = weight;
  if(!isData){
    current_weight = current_weight * trigger_sf * muon_id_iso_sf * MuTrkEffSF; // No SFs so far
  }
  
  if(NonPromptRun){
    double FR_weight = GetFRWeight_SB(muons, muon_tight_id);
    current_weight = FR_weight;
  }
  
  FillCLHist(hnpairmm, "SR1_" + muon_tight_id, eventbase->GetEvent(), muons, electrons_veto, jets_pt30, current_weight, 0);
  
}

void HN_pair_MM::Control_region_1(TString muon_tight_id, TString muon_loose_id, bool trigger_pass, std::vector<snu::KJet> jets, std::vector<snu::KElectron> electrons_veto, bool NonPromptRun){
  //control region for EW bkg. OS mumu + 0 or 1 jet
  bool debug = false;
  
  if(!trigger_pass) return;
  
  std::vector<snu::KMuon> muon_Nocut = GetMuons("MUON_NOCUT",true);
  std::vector<snu::KMuon> muons;
  for(int i = 0; i < muon_Nocut.size(); i++){
    if( !NonPromptRun && PassID(muon_Nocut.at(i), muon_tight_id) && (muon_Nocut.at(i).RelMiniIso() < 0.20) ) muons.push_back(muon_Nocut.at(i));
    if( NonPromptRun  && PassID(muon_Nocut.at(i), muon_loose_id) && (muon_Nocut.at(i).RelMiniIso() < 0.20) ) muons.push_back(muon_Nocut.at(i));//store loose muon for fake bkg
  }

  // -- return if Nmuon != 2 && take OS sign
  int Nmuons;
  Nmuons = muons.size();
  if(Nmuons != 2) return;
  if(muons.at(0).Charge() == muons.at(1).Charge()) return;


  // -- return ! Pt(1st muon) > 60 GeV && Pt(2nd muon) > 20 GeV
  double mu_1st_pt, mu_2nd_pt;
  mu_1st_pt = muons.at(0).Pt();
  mu_2nd_pt = muons.at(1).Pt();
  if(mu_1st_pt < 60 || mu_2nd_pt < 20) return;


  std::vector<snu::KJet> cleaned_jets;
  cleaned_jets = Remove_muon_from_jet(jets, muons);

  if(debug) cout << "3" << endl;

  // -- get lepton cleaned jet collection and apply pt > 30 cut again
  std::vector<snu::KJet> jets_pt30;
  for(int i = 0; i < cleaned_jets.size(); i++){
    if(cleaned_jets.at(i).Pt() > 30) jets_pt30.push_back(cleaned_jets.at(i));
  }

  if(jets_pt30.size() > 2) return;
  
  double trigger_sf = 1.; // 1st muon Pt > 60 GeV with HLT_Mu50 would lead to trigger SF ~ 1. Approx. for now                                                                                                                                                                   

  double muon_id_iso_sf = mcdata_correction->MuonScaleFactor(muon_tight_id, muons, 0);
  double MuTrkEffSF = mcdata_correction->MuonTrackingEffScaleFactor(muons);

  double current_weight = weight;
  if(!isData){
    current_weight = current_weight * trigger_sf * muon_id_iso_sf * MuTrkEffSF;
  }

  // -- Define dimu and cut on dimu.M > 10
  snu::KParticle dimu = muons.at(0) + muons.at(1);
  if(dimu.M() < 10) return;
  
  FillHist("CR1_" + muon_tight_id + "_M_mumu", dimu.M(), current_weight, 10., 1010., 200);
  FillHist("CR1_" + muon_tight_id + "_mu1_Pt", muons.at(0).Pt(), current_weight, 50., 1050., 200);
  FillHist("CR1_" + muon_tight_id + "_mu2_Pt", muons.at(1).Pt(), current_weight, 10., 1010., 200);
  FillHist("CR1_" + muon_tight_id + "_mu1_Eta", muons.at(0).Eta(), current_weight, -3., 3., 600);
  FillHist("CR1_" + muon_tight_id + "_mu2_Eta", muons.at(1).Eta(), current_weight, -3., 3., 600);
  FillHist("CR1_" + muon_tight_id + "_passed_events", 1.5, current_weight, 0., 5., 5);
}




std::vector<snu::KJet> HN_pair_MM::Remove_muon_from_jet(std::vector<snu::KJet> jets, std::vector<snu::KMuon> muons){
  //this function removes muon's 4 momentum from Non-lepton veto jet's 4 momentum if muon is inside of it (dR < 0.4)
  
  if(muons.size() != 2) return jets;
  
  std::vector<snu::KJet> cleaned_jets;  
  int index_closest_jet_mu_1, index_closest_jet_mu_2;
  double dR_min_mu1 = 999., dR_min_mu2 = 999.;
  bool mu_1_inside_jet = false;
  bool mu_2_inside_jet = false;
  for(int i = 0; i < jets.size(); i++){
    double current_dR_mu1 = jets.at(i).DeltaR(muons.at(0));
    double current_dR_mu2 = jets.at(i).DeltaR(muons.at(1));
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
	if(mu_1_inside_jet) jet_mu1 = Clean_jet_lepton(jets.at(i), muons.at(0));
	else jet_mu1 = jets.at(i);
	
	snu::KJet jet_mu2;
	if(mu_2_inside_jet) jet_mu2 = Clean_jet_lepton(jet_mu1, muons.at(1));
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
	jet_mu1 = Clean_jet_lepton(jets.at(i), muons.at(0));
	cleaned_jets.push_back(jet_mu1);
      }
      else if(i == index_closest_jet_mu_2 && mu_2_inside_jet){
	snu::KJet jet_mu2;
	jet_mu2 = Clean_jet_lepton(jets.at(i), muons.at(1));
	cleaned_jets.push_back(jet_mu2);
      }
      else cleaned_jets.push_back(jets.at(i));
    }
    
    return cleaned_jets;
  }
  else return jets;

}

snu::KJet HN_pair_MM::Clean_jet_lepton(snu::KJet jet, snu::KMuon muon){
  double energy = jet.E() - muon.E();
  double phi = ( jet.E() * jet.Phi() - muon.E() * muon.Phi() ) / energy;
  double eta = ( jet.E() * jet.Eta() - muon.E() * muon.Eta() ) / energy;
  double cos_theta = TMath::CosH(eta);
  double momentum = TMath::Sqrt( energy * energy - jet.M() * jet.M() );//keep mass of jet as its original's one
  double Pt = momentum * cos_theta;
  jet.SetPtEtaPhiE(Pt, eta, phi, energy);
  
  return jet;
}

void HN_pair_MM::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void HN_pair_MM::BeginCycle() throw( LQError ){
  
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


HN_pair_MM::~HN_pair_MM() {
  
  Message("In HN_pair_MM Destructor" , INFO);
  //if(!k_isdata)delete reweightPU;
  
}



void HN_pair_MM::FillCutFlow(TString cut, float weight){

  
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


void HN_pair_MM::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


void HN_pair_MM::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this HN_pair_MMCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void HN_pair_MM::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}


bool HN_pair_MM::JSFatJetID(snu::KFatJet fatjet){
  
  if( !(fatjet.PrunedMass() > 40) ) return false;
  if( !(fatjet.PrunedMass() < 130) ) return false;
  if( !( fatjet.Tau2()/fatjet.Tau1() < 0.60 ) ) return false;

  return true;

}

bool HN_pair_MM::IsAwayFromFatJet(snu::KJet jet, vector<snu::KFatJet> fatjets){

  for(unsigned int i=0; i<fatjets.size(); i++){
    if( jet.DeltaR( fatjets.at(i) ) < 0.8 ) return false;
  }

  return true;

}

double HN_pair_MM::GetFRWeight_SB(std::vector<snu::KMuon> muon_looseColl, TString tight_ID){
  
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



