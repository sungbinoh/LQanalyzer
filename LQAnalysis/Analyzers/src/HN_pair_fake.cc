// $Id: HN_pair_fake.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQHN_pair_fake Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "HN_pair_fake.h"

//Core includes
#include "Reweight.h"
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (HN_pair_fake);

/**
 *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
 *
 */
HN_pair_fake::HN_pair_fake() :  AnalyzerCore(), out_muons(0)  {
  // To have the correct name in the log:                                                                                                                            
  SetLogName("HN_pair_fake");
  
  Message("In HN_pair_fake constructor", INFO);
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


void HN_pair_fake::InitialiseAnalysis() throw( LQError ) {

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

void HN_pair_fake::ExecuteEvents()throw( LQError ){
  
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
  if(!isData) weight = pileup_reweight;
  //if(!isData) weight = pileup_reweight; // only for signal
  
  
  // -- call truth particle
  std::vector<snu::KTruth> truthColl;
  eventbase->GetTruthSel()->Selection(truthColl);
  
  bool has_Nmix = false;
  int another_neutrino_1 = 9999999, another_neutrino_2 = 9999999;
  
  
  bool mumu_signal = std::find(k_flags.begin(), k_flags.end(), "hn_pair_mm") != k_flags.end();
  bool ee_signal = std::find(k_flags.begin(), k_flags.end(), "hn_pair_ee") != k_flags.end();
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


  //if(!mu50_pass && !diele_pass) return;
  //cout << "dimu_pass : " << dimu_pass << ". mu50_pass ; " << mu50_pass << endl;
  //if(mu50_pass) cout << "mu50_pass" << endl;
  
  // -- Get Veto Electrons, and store number of them
  //std::vector<snu::KElectron> electrons_susy_veto = GetElectrons("ELECTRON_HN_VETO");
  std::vector<snu::KElectron> electrons_susy_veto = GetElectrons("ELECTRON_SUSY_HNPAIR_VETO");
  std::vector<snu::KElectron> electrons_veto;
  for(int i = 0; i < electrons_susy_veto.size(); i++){
    if( electrons_susy_veto.at(i).PFRelMiniIso(false) < 0.60 )  electrons_veto.push_back(electrons_susy_veto.at(i));

  }
  int N_veto_ele = electrons_veto.size();
  
  // -- Get Veto Muons, return if there are not exactly two veto muons
  //std::vector<snu::KMuon> muons_susy_veto = GetMuons("MUON_HN_VETO");
  std::vector<snu::KMuon> muons_susy_veto = GetMuons("MUON_SUSY_VETO");
  std::vector<snu::KMuon> muons_veto;
  for(int i = 0; i < muons_susy_veto.size(); i++){
    if( muons_susy_veto.at(i).PFRelMiniIsoRho() < 0.60 ) muons_veto.push_back(muons_susy_veto.at(i));
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
  
  std::vector<snu::KElectron> electron_sch_tight = GetElectrons("ELECTRON_HN_TIGHT",true);
  //TString electron_tight_id = "ELECTRON_SUSY_TIGHT";
  TString electron_tight_id = "ELECTRON_SUSY_HNPAIR_TIGHT";
  //TString electron_tight_id = "ELECTRON_HN_TIGHTv4_miniiso";
  TString electron_loose_id = "ELECTRON_SUSY_HNPAIR_VETO";
    

  // -- check with s-ch ID
  Measure_FR_muon(muon_tight_id, "MUON_SUSY_VETO", muons_veto, false);
  Measure_FR_electron(electron_tight_id, "ELECTRON_SUSY_HNPAIR_VETO", electrons_veto, false);
  
  
  return;
}// End of execute event loop


void HN_pair_fake::Measure_FR_muon(TString tight_ID, TString loose_ID, std::vector<snu::KMuon> loose_mu, bool debugging){
  
  bool DijetFake = std::find(k_flags.begin(), k_flags.end(), "DijetFake") != k_flags.end(); //flag for QCD MC fr cal
  bool DijetPrompt = std::find(k_flags.begin(), k_flags.end(), "DijetPrompt") != k_flags.end(); //flag for MC subtraction
  
  if(loose_mu.size() != 1) return;

  std::vector<snu::KMuon> hnloose;//keep only prompt muon for MC subtraction, keep only fake muon for QCD fr cal
  hnloose.clear();
  
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
  
  if(hnloose.size() != 1) return;
  snu::KMuon muon = hnloose.at(0);
  if(muon.Pt() < 30) return;
  
  bool IsThisTight = PassID(muon, tight_ID) && (muon.PFRelMiniIsoRho() < 0.20);
  TString current_name = tight_ID + "FR";
  FillDenAndNum(current_name, muon, weight, IsThisTight);
  

}

void HN_pair_fake::Measure_FR_electron(TString tight_ID, TString loose_ID, std::vector<snu::KElectron> loose_el, bool debugging){
  

  bool DijetFake = std::find(k_flags.begin(), k_flags.end(), "DijetFake") != k_flags.end(); //flag for QCD MC fr cal
  bool DijetPrompt = std::find(k_flags.begin(), k_flags.end(), "DijetPrompt") != k_flags.end(); //flag for MC subtraction

  if(loose_el.size() != 1) return;

  std::vector<snu::KElectron> hnloose;//keep only prompt muon for MC subtraction, keep only fake muon for QCD fr cal
  hnloose.clear();

  for(unsigned int i=0; i<loose_el.size(); i++){
    if(k_isdata) hnloose.push_back( loose_el.at(i) );
    if(!k_isdata){
      if(DijetFake){
        if( !TruthMatched( loose_el.at(i), true) ){
          hnloose.push_back( loose_el.at(i) );
        }
      }
      else if(DijetPrompt){
        if( TruthMatched( loose_el.at(i), true) ){
          hnloose.push_back( loose_el.at(i) );
        }
      }
      else{
        return;
      }
    }
  }

  if(hnloose.size() != 1) return;
  snu::KElectron electron = hnloose.at(0);
  if(electron.Pt() < 30) return;
  
  bool IsThisTight = PassID(electron, tight_ID) && (electron.PFRelMiniIso(false) < 0.10);
  TString current_name = tight_ID + "FR";
  FillDenAndNum(current_name, electron, weight, IsThisTight);

}

void HN_pair_fake::FillDenAndNum(TString prefix, snu::KMuon muon, double thisweight, bool isTight){

  float etaarray [] = {0.0, 0.8, 1.479, 2.5};
  float ptarray [] = {30., 40., 50., 100., 150., 200., 300., 400., 500., 600., 700., 800., 900., 1000.};
  int n_eta = 3, n_pt = 13;
  double TightISO = 0.2;
  double conept = muon.Pt()*(1 + max(0.,(muon.PFRelMiniIsoRho()-TightISO) ) );
  
  FillHist(prefix+"_events_pt_vs_eta_L", conept, fabs(muon.Eta()), thisweight, ptarray, n_pt, etaarray, n_eta);
  
  if( isTight ){
    FillHist(prefix+"_events_pt_vs_eta_T", conept, fabs(muon.Eta()), thisweight, ptarray, n_pt, etaarray, n_eta);
  }

}

void HN_pair_fake::FillDenAndNum(TString prefix, snu::KElectron electron, double thisweight, bool isTight){

  float etaarray [] = {0.0, 0.8, 1.479, 2.5};
  float ptarray [] = {30., 40., 50., 100., 150., 200., 300., 400., 500., 600., 700., 800., 900., 1000.};
  int n_eta = 3, n_pt = 13;
  double TightISO = 0.1;
  double conept = electron.Pt()*(1 + max(0.,(electron.PFRelMiniIso(false)-TightISO) ) );
  
  FillHist(prefix+"_events_pt_vs_eta_L", conept, fabs(electron.SCEta()), thisweight, ptarray, n_pt, etaarray, n_eta);

  if( isTight ){
    FillHist(prefix+"_events_pt_vs_eta_T", conept, fabs(electron.SCEta()), thisweight, ptarray, n_pt, etaarray, n_eta);
  }

}




void HN_pair_fake::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void HN_pair_fake::BeginCycle() throw( LQError ){
  
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


HN_pair_fake::~HN_pair_fake() {
  
  Message("In HN_pair_fake Destructor" , INFO);
  //if(!k_isdata)delete reweightPU;
  
}



void HN_pair_fake::FillCutFlow(TString cut, float weight){

  
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


void HN_pair_fake::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


void HN_pair_fake::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this HN_pair_fakeCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void HN_pair_fake::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}


bool HN_pair_fake::JSFatJetID(snu::KFatJet fatjet){
  
  if( !(fatjet.SoftDropMass() > 60) ) return false;
  //if( !( fatjet.Tau2()/fatjet.Tau1() < 0.60 ) ) return false;

  return true;

}

bool HN_pair_fake::IsAwayFromFatJet(snu::KJet jet, vector<snu::KFatJet> fatjets){

  for(unsigned int i=0; i<fatjets.size(); i++){
    if( jet.DeltaR( fatjets.at(i) ) < 1.0 ) return false;
  }

  return true;

}

double HN_pair_fake::GetFRWeight_SB(std::vector<snu::KMuon> muon_looseColl, TString tight_ID){
  
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



