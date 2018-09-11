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
  MakeCleverHistograms(hnpairmm, "SR1_OS_DiMu");
  MakeCleverHistograms(hnpairmm, "SR1_OS_DiEle");
  MakeCleverHistograms(hnpairmm, "SR1_OS_EMu");
  MakeCleverHistograms(hnpairmm, "CR1_OS_DiMu");
  MakeCleverHistograms(hnpairmm, "CR1_OS_DiEle");
  MakeCleverHistograms(hnpairmm, "CR1_OS_EMu");
  MakeCleverHistograms(hnpairmm, "CR2_OS_DiMu");
  MakeCleverHistograms(hnpairmm, "CR2_OS_DiEle");
  MakeCleverHistograms(hnpairmm, "CR2_OS_EMu");
  MakeCleverHistograms(hnpairmm, "CR3_OS_DiMu");
  MakeCleverHistograms(hnpairmm, "CR3_OS_DiEle");
  MakeCleverHistograms(hnpairmm, "CR3_OS_EMu");
  MakeCleverHistograms(hnpairmm, "CR4_OS_DiMu");
  MakeCleverHistograms(hnpairmm, "CR4_OS_DiEle");
  MakeCleverHistograms(hnpairmm, "CR4_OS_EMu");
  MakeCleverHistograms(hnpairmm, "CR5_OS_DiMu");
  MakeCleverHistograms(hnpairmm, "CR5_OS_DiEle");
  MakeCleverHistograms(hnpairmm, "CR5_OS_EMu");
  
  MakeCleverHistograms(hnpairmm, "SR1_SS_DiMu");
  MakeCleverHistograms(hnpairmm, "SR1_SS_DiEle");
  MakeCleverHistograms(hnpairmm, "SR1_SS_EMu");
  MakeCleverHistograms(hnpairmm, "CR1_SS_DiMu");
  MakeCleverHistograms(hnpairmm, "CR1_SS_DiEle");
  MakeCleverHistograms(hnpairmm, "CR1_SS_EMu");
  MakeCleverHistograms(hnpairmm, "CR2_SS_DiMu");
  MakeCleverHistograms(hnpairmm, "CR2_SS_DiEle");
  MakeCleverHistograms(hnpairmm, "CR2_SS_EMu");
  MakeCleverHistograms(hnpairmm, "CR3_SS_DiMu");
  MakeCleverHistograms(hnpairmm, "CR3_SS_DiEle");
  MakeCleverHistograms(hnpairmm, "CR3_SS_EMu");
  MakeCleverHistograms(hnpairmm, "CR4_SS_DiMu");
  MakeCleverHistograms(hnpairmm, "CR4_SS_DiEle");
  MakeCleverHistograms(hnpairmm, "CR4_SS_EMu");
  MakeCleverHistograms(hnpairmm, "CR5_SS_DiMu");
  MakeCleverHistograms(hnpairmm, "CR5_SS_DiEle");
  MakeCleverHistograms(hnpairmm, "CR5_SS_EMu");
  
  // -- CLH for HN mass cut
  MakeCleverHistograms(hnpairmm, "SR1_OS_DiMu_HNcut");
  MakeCleverHistograms(hnpairmm, "SR1_OS_DiEle_HNcut");
  MakeCleverHistograms(hnpairmm, "SR1_OS_EMu_HNcut");
  
  MakeCleverHistograms(hnpairmm, "SR1_SS_DiMu_HNcut");
  MakeCleverHistograms(hnpairmm, "SR1_SS_DiEle_HNcut");
  MakeCleverHistograms(hnpairmm, "SR1_SS_EMu_HNcut");
    
  MakeCleverHistograms(hnpairmm, "SR1_SS_DiEle_CF_HNcut");
  MakeCleverHistograms(hnpairmm, "SR1_SS_EMu_CF_HNcut");
    
  //CLH for CF bkg
  MakeCleverHistograms(hnpairmm, "SR1_SS_DiEle_CF");
  MakeCleverHistograms(hnpairmm, "SR1_SS_EMu_CF");
  MakeCleverHistograms(hnpairmm, "CR1_SS_DiEle_CF");
  MakeCleverHistograms(hnpairmm, "CR1_SS_EMu_CF");
  MakeCleverHistograms(hnpairmm, "CR2_SS_DiEle_CF");
  MakeCleverHistograms(hnpairmm, "CR2_SS_EMu_CF");
  MakeCleverHistograms(hnpairmm, "CR3_SS_DiEle_CF");
  MakeCleverHistograms(hnpairmm, "CR3_SS_EMu_CF");
  MakeCleverHistograms(hnpairmm, "CR4_SS_DiEle_CF");
  MakeCleverHistograms(hnpairmm, "CR4_SS_EMu_CF");
  MakeCleverHistograms(hnpairmm, "CR5_SS_DiEle_CF");
  MakeCleverHistograms(hnpairmm, "CR5_SS_EMu_CF");
    
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


  TDirectory* origDir = gDirectory;
  
  TString lqdir = getenv("LQANALYZER_DIR");

  //TFile *file_data_fake = new TFile( lqdir+"/data/Fake/80X/FakeRate_QCD.root");
  TFile *file_data_fake = new TFile( lqdir+"/data/Fake/80X/Data_driven_FR_syst_DoubleMuon_Nvtx_Reweight.root");
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
    hist_data_driven[histnames_2D.at(i)] = (TH2D*)file_data_fake -> Get(histnames_2D.at(i))->Clone();
  }

  //Message("saving hist strings ends", INFO);
  file_data_fake -> Close();
  delete file_data_fake;


  // -- electron fake rate
  TFile *file_data_fake_e = new TFile( lqdir+"/data/Fake/80X/FakeRate_QCD.root");
  TIter next_e(gDirectory->GetListOfKeys());
  TKey *key_e;
  vector<TString> histnames_e;
  vector<TString> histnames_2D_e;
  while ((key_e = (TKey*)next_e())) {
    TClass *cl = gROOT->GetClass(key_e->GetClassName());
    if (cl->InheritsFrom("TH2D")) histnames_2D_e.push_back(key_e -> GetName());
    if (!cl->InheritsFrom("TH1")) continue;
    histnames_e.push_back(key_e -> GetName());
  }
  gROOT->cd();
  tempDir->cd();
  
  Message("saving hist strings", INFO);
  for(int i = 0; i < histnames_2D_e.size(); i ++){
    Message(histnames_2D_e.at(i), INFO);
    hist_data_driven[histnames_2D_e.at(i)] = (TH2D*)file_data_fake_e -> Get(histnames_2D_e.at(i))->Clone();
  }
  Message("saving hist strings ends", INFO);

  file_data_fake_e -> Close();
  delete file_data_fake_e;




  // -- Charge flip probability
  TFile *file_data_CF = new TFile( lqdir+"/data/Fake/80X/CF_result_MC.root");
  TIter next_CF(gDirectory->GetListOfKeys());
  TKey *key_CF;
  vector<TString> histnames_CF;
  vector<TString> histnames_2D_CF;
  while ((key_CF = (TKey*)next_CF())) {
    TClass *cl = gROOT->GetClass(key_CF->GetClassName());
    if (cl->InheritsFrom("TH2D")) histnames_2D_CF.push_back(key_CF -> GetName());
    if (!cl->InheritsFrom("TH1")) continue;
    histnames_CF.push_back(key_CF -> GetName());
  }
  gROOT->cd();

  tempDir->cd();
  Message("saving hist strings", INFO);
  for(int i = 0; i < histnames_CF.size(); i ++){
    Message(histnames_CF.at(i), INFO);
    hist_CF[histnames_CF.at(i)] = (TH1D*)file_data_CF -> Get(histnames_CF.at(i))->Clone();
  }
  Message("saving hist strings ends", INFO);
  file_data_CF -> Close();
  delete file_data_CF;
  
  
  origDir->cd();
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

  //cout << "truthcoll.size() : " << truthColl.size() << endl;
  cout << "-------------------------------" << endl;
  
  std::vector<unsigned int> HN_index;
  std::vector<unsigned int> Zp_index;
  HN_index.clear();
  Zp_index.clear();
  
  for(unsigned int i_truth = 0; i_truth < truthColl.size(); i_truth++){
    
    //cout << i_truth << " : pdgid = " << truthColl.at(i_truth).PdgId() << endl;
    // -- save HN index
    if(truthColl.at(i_truth).PdgId() == 9900014) HN_index.push_back(i_truth);
    
    // -- save Z' index
    if(truthColl.at(i_truth).PdgId() ==9900023 || truthColl.at(i_truth).PdgId() == 32) Zp_index.push_back(i_truth);
    
  }
  
  std::vector<unsigned int> HN1_daughter_index;
  std::vector<unsigned int> HN2_daughter_index;
  std::vector<TLorentzVector> HN1_quarks;
  std::vector<TLorentzVector> HN2_quarks;
  TLorentzVector HN1_lepton, HN2_lepton;
  bool HN1_lepton_bool = false, HN2_lepton_bool = false;
  HN1_daughter_index.clear();
  HN2_daughter_index.clear();
  for(unsigned int i_truth = 0; i_truth < truthColl.size(); i_truth++){
    if( truthColl.at(i_truth).IndexMother() == HN_index.at(0) ){
      HN1_daughter_index.push_back(i_truth);
      snu::KTruth current_truth = truthColl.at(i_truth);
      TLorentzVector current_vec;
      current_vec.SetPtEtaPhiM(current_truth.Pt(), current_truth.Eta(), current_truth.Phi(), current_truth.M() );
      if(abs(current_truth.PdgId()) < 10) HN1_quarks.push_back(current_vec);
      else{ 
	HN1_lepton = current_vec;
	HN1_lepton_bool = true;
      }
    }
    if( truthColl.at(i_truth).IndexMother() == HN_index.at(1) ){
      HN2_daughter_index.push_back(i_truth);
      snu::KTruth current_truth= truthColl.at(i_truth);
      TLorentzVector current_vec;
      current_vec.SetPtEtaPhiM(current_truth.Pt(), current_truth.Eta(), current_truth.Phi(), current_truth.M() );
      if(abs(current_truth.PdgId()) < 10) HN2_quarks.push_back(current_vec);
      else{
	HN2_lepton = current_vec;
	HN2_lepton_bool = true;
      }
    }
  }
  
  if(HN1_daughter_index.size() != 3 || HN2_daughter_index.size() != 3) return;
  if(HN1_quarks.size() != 2 || HN2_quarks.size() != 2 || !HN1_lepton_bool || !HN2_lepton_bool) return;


  snu::KTruth HN1 = truthColl.at(HN_index.at(0));
  snu::KTruth HN2 = truthColl.at(HN_index.at(1));

  TLorentzVector HN1_vec, HN2_vec;
  HN1_vec.SetPtEtaPhiM(HN1.Pt(), HN1.Eta(), HN1.Phi(), HN1.M() );
  HN2_vec.SetPtEtaPhiM(HN2.Pt(), HN2.Eta(), HN2.Phi(), HN2.M() );
  
  cout << "before boost" << endl;
  cout << "HN1_vec.Px() : " << HN1_vec.Px() << ", HN1_vec.Py() : " << HN1_vec.Py() << " , HN1_vec.Pz() : " << HN1_vec.Pz() << endl;
  cout << "HN2_vec.Px() : " << HN2_vec.Px() << ", HN2_vec.Py() : " << HN2_vec.Py() << " , HN2_vec.Pz() : " << HN2_vec.Pz() << endl;

  TVector3 HN1_boost, HN2_boost;
  HN1_boost = HN1_vec.BoostVector();
  HN2_boost = HN2_vec.BoostVector();
  
  cout << "after boost" << endl;
  
  HN1_vec.Boost(-1 * HN1_boost);
  HN2_vec.Boost(-1 * HN2_boost);

  cout << "HN1_vec.Px() : " << HN1_vec.Px() << ", HN1_vec.Py() : " << HN1_vec.Py() << " , HN1_vec.Pz() : " << HN1_vec.Pz() << endl;
  cout << "HN2_vec.Px() : " << HN2_vec.Px() << ", HN2_vec.Py() : " << HN2_vec.Py() << " , HN2_vec.Pz() : " << HN2_vec.Pz() << endl;
  



  // -- boost daughters of HNs to rest frame of HNs
  HN1_quarks.at(0).Boost(-1 * HN1_boost);
  HN1_quarks.at(1).Boost(-1 * HN1_boost);
  HN2_quarks.at(0).Boost(-1 * HN2_boost);
  HN2_quarks.at(1).Boost(-1 * HN2_boost);
  HN1_lepton.Boost(-1 * HN1_boost);
  HN2_lepton.Boost(-1 * HN2_boost);
  
  


  double HN1_quarks_P_diff = fabs( HN1_quarks.at(0).P() - HN1_quarks.at(1).P() );
  double HN2_quarks_P_diff = fabs( HN2_quarks.at(0).P() - HN2_quarks.at(1).P() );
  double HN1_WR_M = (HN1_quarks.at(0) + HN1_quarks.at(1)).M();
  double HN2_WR_M = (HN2_quarks.at(0) + HN2_quarks.at(1)).M();
  
  FillHist("HN1_quarks_P_diff", HN1_quarks_P_diff, 1., 0., 1000., 1000);
  FillHist("HN2_quarks_P_diff", HN2_quarks_P_diff, 1., 0., 1000., 1000);
  FillHist("HN_quarks_P_diff", HN1_quarks_P_diff, 1., 0., 1000., 1000);
  FillHist("HN_quarks_P_diff", HN2_quarks_P_diff, 1., 0., 1000., 1000);

  FillHist("HN1_WR_M", HN1_WR_M, 1., 0., 1000., 1000);
  FillHist("HN2_WR_M", HN2_WR_M, 1., 0., 1000., 1000);
  FillHist("HN_WR_M", HN1_WR_M, 1., 0., 1000., 1000);
  FillHist("HN_WR_M", HN2_WR_M, 1., 0., 1000., 1000);

  FillHist("HN1_lepton_pt", HN1_lepton.Pt(), 1., 0., 1000., 1000);
  FillHist("HN2_lepton_pt", HN2_lepton.Pt(), 1., 0., 1000., 1000);
  FillHist("HN_lepton_pt", HN1_lepton.Pt(), 1., 0., 1000., 1000);
  FillHist("HN_lepton_pt", HN2_lepton.Pt(), 1., 0., 1000., 1000);
  
  FillHist("HN1_pt", truthColl.at(HN_index.at(0)).Pt(), 1., 0., 1000., 1000);
  FillHist("HN2_pt", truthColl.at(HN_index.at(1)).Pt(), 1., 0., 1000., 1000);
  FillHist("HN_pt", truthColl.at(HN_index.at(0)).Pt(), 1., 0., 1000., 1000);
  FillHist("HN_pt", truthColl.at(HN_index.at(1)).Pt(), 1., 0., 1000., 1000);
  
  cout << "[[before rotation]]" << endl;
  cout << "HN1_quarks.at(0) P = (" << HN1_quarks.at(0).Px() << ", " << HN1_quarks.at(0).Py()  << ", " << HN1_quarks.at(0).Pz()  << ") = " << HN1_quarks.at(0).P() << endl;
  cout << "HN1_quarks.at(1) P = (" << HN1_quarks.at(1).Px() << ", " << HN1_quarks.at(1).Py()  << ", " << HN1_quarks.at(1).Pz()  << ") = " << HN1_quarks.at(1).P() << endl;
  cout << "HN1_lepton       P = (" << HN1_lepton.Px() << ", " << HN1_lepton.Py()  << ", " << HN1_lepton.Pz()  << ") = " << HN1_lepton.P() << endl;
  

  // -- rotate to X-Y plane for two quarts + lepton from same HN
  double HN1_lepton_phi = HN1_lepton.Phi();
  double HN1_lepton_theta = HN1_lepton.Theta();
  
  HN1_quarks.at(0).RotateZ(-1 * HN1_lepton_phi);
  HN1_quarks.at(1).RotateZ(-1 * HN1_lepton_phi);
  HN1_lepton.RotateZ(-1 * HN1_lepton_phi);

  
  HN1_quarks.at(0).RotateY(-1 * HN1_lepton_theta);
  HN1_quarks.at(1).RotateY(-1 * HN1_lepton_theta);
  HN1_lepton.RotateY(-1 * HN1_lepton_theta); // -- now lepton is on the Z-axis
  
  cout << "[[lepton is on Z axis]]" << endl;
  cout << "HN1_quarks.at(0) P = (" << HN1_quarks.at(0).Px() << ", " << HN1_quarks.at(0).Py()  << ", " << HN1_quarks.at(0).Pz()  << ") = " << HN1_quarks.at(0).P() << endl;
  cout << "HN1_quarks.at(1) P = (" << HN1_quarks.at(1).Px() << ", " << HN1_quarks.at(1).Py()  << ", " << HN1_quarks.at(1).Pz()  << ") = " << HN1_quarks.at(1).P() << endl;
  cout << "HN1_lepton       P = (" << HN1_lepton.Px() << ", " << HN1_lepton.Py()  << ", " << HN1_lepton.Pz()  << ") = " << HN1_lepton.P() << endl;


  double HN1_quark1_phi = HN1_quarks.at(0).Phi();
  HN1_quarks.at(0).RotateZ(-1 * HN1_quark1_phi);
  HN1_quarks.at(1).RotateZ(-1 * HN1_quark1_phi);
  HN1_lepton.RotateZ(-1 * HN1_quark1_phi); // -- now three particles are on Z-X plane
  
  HN1_quarks.at(0).RotateX(-0.5 * TMath::Pi());
  HN1_quarks.at(1).RotateX(-0.5 * TMath::Pi());
  HN1_lepton.RotateX(-0.5 * TMath::Pi());
  
  cout << "[[after all rotations]]" << endl;
  cout << "HN1_quarks.at(0) P = (" << HN1_quarks.at(0).Px() << ", " << HN1_quarks.at(0).Py()  << ", " << HN1_quarks.at(0).Pz()  << ") = " << HN1_quarks.at(0).P() << endl;
  cout << "HN1_quarks.at(1) P = (" << HN1_quarks.at(1).Px() << ", " << HN1_quarks.at(1).Py()  << ", " << HN1_quarks.at(1).Pz()  << ") = " << HN1_quarks.at(1).P() << endl;
  cout << "HN1_lepton       P = (" << HN1_lepton.Px() << ", " << HN1_lepton.Py()  << ", " << HN1_lepton.Pz()  << ") = " << HN1_lepton.P() << endl;
  
  

  double HN2_lepton_phi = HN2_lepton.Phi();
  double HN2_lepton_theta = HN2_lepton.Theta();

  HN2_quarks.at(0).RotateZ(-1 * HN2_lepton_phi);
  HN2_quarks.at(1).RotateZ(-1 * HN2_lepton_phi);
  HN2_lepton.RotateZ(-1 * HN2_lepton_phi);


  HN2_quarks.at(0).RotateY(-1 * HN2_lepton_theta);
  HN2_quarks.at(1).RotateY(-1 * HN2_lepton_theta);
  HN2_lepton.RotateY(-1 * HN2_lepton_theta); // -- now lepton is on the Z-axis                                                                                                                                                                                                  

  cout << "[[lepton is on Z axis]]" << endl;
  cout << "HN2_quarks.at(0) P = (" << HN2_quarks.at(0).Px() << ", " << HN2_quarks.at(0).Py()  << ", " << HN2_quarks.at(0).Pz()  << ") = " << HN2_quarks.at(0).P() << endl;
  cout << "HN2_quarks.at(1) P = (" << HN2_quarks.at(1).Px() << ", " << HN2_quarks.at(1).Py()  << ", " << HN2_quarks.at(1).Pz()  << ") = " << HN2_quarks.at(1).P() << endl;
  cout << "HN2_lepton       P = (" << HN2_lepton.Px() << ", " << HN2_lepton.Py()  << ", " << HN2_lepton.Pz()  << ") = " << HN2_lepton.P() << endl;


  double HN2_quark1_phi = HN2_quarks.at(0).Phi();
  HN2_quarks.at(0).RotateZ(-1 * HN2_quark1_phi);
  HN2_quarks.at(1).RotateZ(-1 * HN2_quark1_phi);
  HN2_lepton.RotateZ(-1 * HN2_quark1_phi); // -- now three particles are on Z-X plane                                                                                                                                                                                           

  HN2_quarks.at(0).RotateX(-0.5 * TMath::Pi());
  HN2_quarks.at(1).RotateX(-0.5 * TMath::Pi());
  HN2_lepton.RotateX(-0.5 * TMath::Pi());

  cout << "[[after all rotations]]" << endl;
  cout << "HN2_quarks.at(0) P = (" << HN2_quarks.at(0).Px() << ", " << HN2_quarks.at(0).Py()  << ", " << HN2_quarks.at(0).Pz()  << ") = " << HN2_quarks.at(0).P() << endl;
  cout << "HN2_quarks.at(1) P = (" << HN2_quarks.at(1).Px() << ", " << HN2_quarks.at(1).Py()  << ", " << HN2_quarks.at(1).Pz()  << ") = " << HN2_quarks.at(1).P() << endl;
  cout << "HN2_lepton       P = (" << HN2_lepton.Px() << ", " << HN2_lepton.Py()  << ", " << HN2_lepton.Pz()  << ") = " << HN2_lepton.P() << endl;
  
  FillHist("HN1_lqq_dist", HN1_quarks.at(0).Px(), HN1_quarks.at(0).Py(), 1., -500., 500., 100, -500., 500., 100, "", "");
  FillHist("HN1_lqq_dist", HN1_quarks.at(1).Px(), HN1_quarks.at(1).Py(), 1., -500., 500., 100, -500., 500., 100, "", "");
  FillHist("HN1_lqq_dist", HN1_lepton.Px(), HN1_lepton.Py(), 1., -500., 500., 100, -500., 500., 100, "", "");
 
  FillHist("HN2_lqq_dist", HN2_quarks.at(0).Px(), HN2_quarks.at(0).Py(), 1., -500., 500., 100, -500., 500., 100, "", "");
  FillHist("HN2_lqq_dist", HN2_quarks.at(1).Px(), HN2_quarks.at(1).Py(), 1., -500., 500., 100, -500., 500., 100, "", "");
  FillHist("HN2_lqq_dist", HN2_lepton.Px(), HN2_lepton.Py(), 1., -500., 500., 100, -500., 500., 100, "", "");

  FillHist("HN_lqq_dist", HN1_quarks.at(0).Px(), HN1_quarks.at(0).Py(), 1., -500., 500., 100, -500., 500., 100, "", "");
  FillHist("HN_lqq_dist", HN1_quarks.at(1).Px(), HN1_quarks.at(1).Py(), 1., -500., 500., 100, -500., 500., 100, "", "");
  FillHist("HN_lqq_dist", HN1_lepton.Px(), HN1_lepton.Py(), 1., -500., 500., 100, -500., 500., 100, "", "");
  FillHist("HN_lqq_dist", HN2_quarks.at(0).Px(), HN2_quarks.at(0).Py(), 1., -500., 500., 100, -500., 500., 100, "", "");
  FillHist("HN_lqq_dist", HN2_quarks.at(1).Px(), HN2_quarks.at(1).Py(), 1., -500., 500., 100, -500., 500., 100, "", "");
  FillHist("HN_lqq_dist", HN2_lepton.Px(), HN2_lepton.Py(), 1., -500., 500., 100, -500., 500., 100, "", "");

  FillHist("HN_qq_dist", HN1_quarks.at(0).Px(), HN1_quarks.at(0).Py(), 1., -500., 500., 100, -500., 500., 100, "", "");
  FillHist("HN_qq_dist", HN1_quarks.at(1).Px(), HN1_quarks.at(1).Py(),1., -500., 500., 100, -500., 500., 100, "", "");
  FillHist("HN_qq_dist", HN2_quarks.at(0).Px(), HN2_quarks.at(0).Py(),1., -500., 500., 100, -500., 500., 100, "", "");
  FillHist("HN_qq_dist", HN2_quarks.at(1).Px(), HN2_quarks.at(1).Py(),1., -500., 500., 100, -500., 500., 100, "", "");
  
  
  return;//return to run only truthColl






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
  
  std::vector<snu::KMuon> muon_sch_tight = GetMuons("MUON_HN_TIGHT",true);
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


  // -- Calculate Fake rate weight for SS dilpeton 
  weight_fake = 0;
  if(NonPromptRun){
    std::vector<KLepton> Leptons;
    bool lep_1_tight = false, lep_2_tight = false;
    if(N_electron_SS == 2 && N_veto_ele == 2 && N_veto_muon == 0){
      Leptons.push_back(electrons_SS.at(0));
      Leptons.push_back(electrons_SS.at(1));
      if(PassID(electrons_SS.at(0), electron_tight_id) && (electrons_SS.at(0).PFRelMiniIso(false) < 0.10) ) lep_1_tight = true;
      if(PassID(electrons_SS.at(1), electron_tight_id) && (electrons_SS.at(1).PFRelMiniIso(false) < 0.10) ) lep_2_tight = true;
      
      if(Leptons.at(0).Charge() == Leptons.at(1).Charge()){
	weight_fake = GetFRWeight_SB(Leptons, lep_1_tight, lep_2_tight);
	cout << "weight_fake : " << weight_fake << endl;
      }
      //else return;
    }
    else if(N_muon_SS == 2 && N_veto_muon == 2 && N_veto_ele == 0){
      Leptons.push_back(muons_SS.at(0));
      Leptons.push_back(muons_SS.at(1));
      if(PassID(muons_SS.at(0), muon_tight_id) && (muons_SS.at(0).PFRelMiniIsoRho() < 0.20) ) lep_1_tight = true;
      if(PassID(muons_SS.at(1), muon_tight_id) && (muons_SS.at(1).PFRelMiniIsoRho() < 0.20) ) lep_2_tight = true;
      
      if(Leptons.at(0).Charge() == Leptons.at(1).Charge()) weight_fake = GetFRWeight_SB(Leptons, lep_1_tight, lep_2_tight);
      //else return;
    }
    else if(N_muon_SS == 1 && N_veto_muon ==1 && N_electron_SS == 1 && N_veto_ele == 1){
      Leptons.push_back(electrons_SS.at(0));
      Leptons.push_back(muons_SS.at(0));
      if(PassID(electrons_SS.at(0), electron_tight_id) && (electrons_SS.at(0).PFRelMiniIso(false) < 0.10) ) lep_1_tight = true;
      if(PassID(muons_SS.at(0), muon_tight_id) && (muons_SS.at(0).PFRelMiniIsoRho() < 0.20) ) lep_2_tight = true;
      
      if(Leptons.at(0).Charge() == Leptons.at(1).Charge()) weight_fake = GetFRWeight_SB(Leptons, lep_1_tight, lep_2_tight);
      //else return;
    }
    else weight_fake = 0;
  }
  

  // -- calculate CF weight for EMU or EE channel
  weight_CF = 0;
  std::vector<KLepton> Leptons_CF;
  if(N_electron_SS == 2 && N_veto_ele == 2 && N_veto_muon == 0){
    Leptons_CF.push_back(electrons_SS.at(0));
    Leptons_CF.push_back(electrons_SS.at(1));
    if(Leptons_CF.at(0).Charge() != Leptons_CF.at(1).Charge()){
      weight_CF = GetCFWeight_SB(Leptons_CF);
      //cout << "ee weight_CF : " << weight_CF << endl;
    }
  }
  else if(N_muon_SS == 1 && N_veto_muon ==1 && N_electron_SS == 1 && N_veto_ele == 1){
    Leptons_CF.push_back(electrons_SS.at(0));
    Leptons_CF.push_back(muons_SS.at(0));
    if(Leptons_CF.at(0).Charge() != Leptons_CF.at(1).Charge()) {
      weight_CF = GetCFWeight_SB(Leptons_CF);
      //cout << "emu weight_CF : " << weight_CF << endl; 
    }
  }
  else weight_CF = 0;
  
  // -- run Jobs with OS lepton definition
  
  Signal_region_1("OS_DiMu", mu50_pass, jets_lepveto, fatjets, electrons, muons, N_electron, N_veto_ele, N_muon, N_veto_muon, NonPromptRun);
  Signal_region_1("OS_EMu", mu50_pass, jets_lepveto, fatjets, electrons, muons, N_electron, N_veto_ele, N_muon, N_veto_muon, NonPromptRun);
  Signal_region_1("OS_DiEle", diele_pass, jets_lepveto, fatjets, electrons, muons, N_electron, N_veto_ele, N_muon, N_veto_muon, NonPromptRun);
  
  Control_region_1("OS_DiMu", mu50_pass, jets_lepveto, fatjets, electrons, muons, N_electron, N_veto_ele, N_muon, N_veto_muon, NonPromptRun);
  Control_region_1("OS_EMu", mu50_pass, jets_lepveto, fatjets, electrons, muons, N_electron, N_veto_ele, N_muon, N_veto_muon, NonPromptRun);
  Control_region_1("OS_DiEle", diele_pass, jets_lepveto, fatjets, electrons, muons, N_electron, N_veto_ele, N_muon, N_veto_muon, NonPromptRun);
  
  Control_region_2("OS_DiMu", mu50_pass, jets_lepveto, fatjets, electrons, muons, N_electron, N_veto_ele, N_muon, N_veto_muon, NonPromptRun);
  Control_region_2("OS_EMu", mu50_pass, jets_lepveto, fatjets, electrons, muons, N_electron, N_veto_ele, N_muon, N_veto_muon, NonPromptRun);
  Control_region_2("OS_DiEle", diele_pass, jets_lepveto, fatjets, electrons, muons, N_electron, N_veto_ele, N_muon, N_veto_muon, NonPromptRun);

  Control_region_3("OS_DiMu", mu50_pass, jets_lepveto, fatjets, electrons, muons, N_electron, N_veto_ele, N_muon, N_veto_muon, NonPromptRun);
  Control_region_3("OS_EMu", mu50_pass, jets_lepveto, fatjets, electrons, muons, N_electron, N_veto_ele, N_muon, N_veto_muon, NonPromptRun);
  Control_region_3("OS_DiEle", diele_pass, jets_lepveto, fatjets, electrons, muons, N_electron, N_veto_ele, N_muon, N_veto_muon, NonPromptRun);

  Control_region_4("OS_DiMu", mu50_pass, jets_lepveto, fatjets, electrons, muons, N_electron, N_veto_ele, N_muon, N_veto_muon, NonPromptRun);
  Control_region_4("OS_EMu", mu50_pass, jets_lepveto, fatjets, electrons, muons, N_electron, N_veto_ele, N_muon, N_veto_muon, NonPromptRun);
  Control_region_4("OS_DiEle", diele_pass, jets_lepveto, fatjets, electrons, muons, N_electron, N_veto_ele, N_muon, N_veto_muon, NonPromptRun);

  // -- run Jobs for SS lepton definition
  Signal_region_1("SS_DiMu", mu50_pass, jets_lepveto, fatjets, electrons_SS, muons_SS, N_electron_SS, N_veto_ele, N_muon_SS, N_veto_muon, NonPromptRun);
  Signal_region_1("SS_EMu", mu50_pass, jets_lepveto, fatjets, electrons_SS, muons_SS, N_electron_SS, N_veto_ele, N_muon_SS, N_veto_muon, NonPromptRun);
  Signal_region_1("SS_DiEle", diele_pass, jets_lepveto, fatjets, electrons_SS, muons_SS, N_electron_SS, N_veto_ele, N_muon_SS, N_veto_muon, NonPromptRun);

  Control_region_1("SS_DiMu", mu50_pass, jets_lepveto, fatjets, electrons_SS, muons_SS, N_electron_SS, N_veto_ele, N_muon_SS, N_veto_muon, NonPromptRun);
  Control_region_1("SS_EMu", mu50_pass, jets_lepveto, fatjets, electrons_SS, muons_SS, N_electron_SS, N_veto_ele, N_muon_SS, N_veto_muon, NonPromptRun);
  Control_region_1("SS_DiEle", diele_pass, jets_lepveto, fatjets, electrons_SS, muons_SS, N_electron_SS, N_veto_ele, N_muon_SS, N_veto_muon, NonPromptRun);

  Control_region_2("SS_DiMu", mu50_pass, jets_lepveto, fatjets, electrons_SS, muons_SS, N_electron_SS, N_veto_ele, N_muon_SS, N_veto_muon, NonPromptRun);
  Control_region_2("SS_EMu", mu50_pass, jets_lepveto, fatjets, electrons_SS, muons_SS, N_electron_SS, N_veto_ele, N_muon_SS, N_veto_muon, NonPromptRun);
  Control_region_2("SS_DiEle", diele_pass, jets_lepveto, fatjets, electrons_SS, muons_SS, N_electron_SS, N_veto_ele, N_muon_SS, N_veto_muon, NonPromptRun);

  Control_region_3("SS_DiMu", mu50_pass, jets_lepveto, fatjets, electrons_SS, muons_SS, N_electron_SS, N_veto_ele, N_muon_SS, N_veto_muon, NonPromptRun);
  Control_region_3("SS_EMu", mu50_pass, jets_lepveto, fatjets, electrons_SS, muons_SS, N_electron_SS, N_veto_ele, N_muon_SS, N_veto_muon, NonPromptRun);
  Control_region_3("SS_DiEle", diele_pass, jets_lepveto, fatjets, electrons_SS, muons_SS, N_electron_SS, N_veto_ele, N_muon_SS, N_veto_muon, NonPromptRun);

  Control_region_4("SS_DiMu", mu50_pass, jets_lepveto, fatjets, electrons_SS, muons_SS, N_electron_SS, N_veto_ele, N_muon_SS, N_veto_muon, NonPromptRun);
  Control_region_4("SS_EMu", mu50_pass, jets_lepveto, fatjets, electrons_SS, muons_SS, N_electron_SS, N_veto_ele, N_muon_SS, N_veto_muon, NonPromptRun);
  Control_region_4("SS_DiEle", diele_pass, jets_lepveto, fatjets, electrons_SS, muons_SS, N_electron_SS, N_veto_ele, N_muon_SS, N_veto_muon, NonPromptRun);
 
  //SS emu no Njet cut
  Control_region_5("SS_EMu", mu50_pass, jets_lepveto, fatjets, electrons_SS, muons_SS, N_electron_SS, N_veto_ele, N_muon_SS, N_veto_muon, electron_tight_id, muon_tight_id, NonPromptRun);
  
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
    current_weight = weight_fake;
  }

  FillCLHist(hnpairmm, "CR2_" + channel, eventbase->GetEvent(), Leptons, N_electron, N_muon, jets_pt40, fatjets, current_weight, 0);
  if(channel.Contains("SS_EMu") || channel.Contains("SS_DiEle"))  FillCLHist(hnpairmm, "CR2_" + channel + "_CF", eventbase->GetEvent(), Leptons, N_electron, N_muon, jets_pt40, fatjets, weight_CF, 0);

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
    current_weight = weight_fake;
  }

  FillCLHist(hnpairmm, "CR3_" + channel, eventbase->GetEvent(), Leptons, N_electron, N_muon, jets_pt40, fatjets, current_weight, nbjet);
  if(channel.Contains("SS_EMu") || channel.Contains("SS_DiEle"))  FillCLHist(hnpairmm, "CR3_" + channel + "_CF", eventbase->GetEvent(), Leptons, N_electron, N_muon, jets_pt40, fatjets, weight_CF, 0);

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
    current_weight = weight_fake;
  }

  FillCLHist(hnpairmm, "CR4_" + channel, eventbase->GetEvent(), Leptons, N_electron, N_muon, jets_pt40, fatjets, current_weight, nbjet);
  if(channel.Contains("SS_EMu") || channel.Contains("SS_DiEle"))  FillCLHist(hnpairmm, "CR4_" + channel + "_CF", eventbase->GetEvent(), Leptons, N_electron, N_muon, jets_pt40, fatjets, weight_CF, 0);
  
}

void HN_pair_all::Control_region_5(TString channel, bool trigger_pass, std::vector<snu::KJet> jets, std::vector<snu::KFatJet> fatjets, std::vector<snu::KElectron> electrons, std::vector<snu::KMuon> muons, int N_electron, int N_veto_ele, int N_muon, int N_veto_muon, 
				   TString electron_tight_id, TString muon_tight_id, bool NonPromptRun){
  //SS emu CR for fake control
  //cout << "start Control_region_5" << endl;
  
  if(!trigger_pass) return;
  
  //cout << "0" << endl;

  
  bool pass_N_lep = false;
  bool lep_1_tight, lep_2_tight, lep_3_tight;
  int N_electron_CLH, N_muon_CLH;
  std::vector<KLepton> Leptons;
  if(channel.Contains("EMu")){
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

  //cout << "1" << endl;

  double Lep_1st_Pt, Lep_2nd_Pt;
  Lep_1st_Pt = Leptons.at(0).Pt();
  Lep_2nd_Pt = Leptons.at(1).Pt();
  
  if(Lep_1st_Pt < 65 || Lep_2nd_Pt < 65) return;
  
  std::vector<snu::KJet> jets_pt40;
  for(int i = 0; i < jets.size(); i++){
    if(jets.at(i).Pt() > 40) jets_pt40.push_back(jets.at(i));
  }
  
  //cout << "3" << endl;

  if(jets_pt40.size() < 2) return;
  
  int nbjet=0;
  for(unsigned int ij=0; ij < jets_pt40.size(); ij++){
    if( IsBTagged(jets_pt40.at(ij),   snu::KJet::CSVv2, snu::KJet::Medium, -1, 0))  nbjet++;
  }
  //if(nbjet != 0) return;
  // -- fill cutflow
 

  if(channel.Contains("EMu")){
    FillHist("signal_eff", 32.5, weight, 0., 40., 40);
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
    cout << "CR5 FR weight : " << current_weight << endl;
  }

  //cout << "CL" << endl;
  
  FillCLHist(hnpairmm, "CR5_" + channel, eventbase->GetEvent(), Leptons, N_electron, N_muon, jets_pt40, fatjets, current_weight, nbjet);
  if(channel.Contains("SS_EMu") || channel.Contains("SS_DiEle"))  FillCLHist(hnpairmm, "CR5_" + channel + "_CF", eventbase->GetEvent(), Leptons, N_electron, N_muon, jets_pt40, fatjets, weight_CF, 0);
  
  
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

double HN_pair_all::GetFRWeight_SB(std::vector<KLepton> loose_lepColl, bool lep_1_tight, bool lep_2_tight){
  
  if(loose_lepColl.size() != 2) return 0.;
  //cout << "loose_lepColl.at(0).LeptonFlavour() : " << loose_lepColl.at(0).LeptonFlavour() << endl;
  //muon : LeptonFlavour() = 1, electron LeptonFlavour() = 2
  //double mu_1_cone_pt = muon_looseColl.at(0).Pt()*(1 + max(0.,(muon_looseColl.at(0).RelIso04()-0.07) ) );
  //double mu_2_cone_pt = muon_looseColl.at(1).Pt()*(1 + max(0.,(muon_looseColl.at(1).RelIso04()-0.07) ) );
  double lep_1_pt = loose_lepColl.at(0).Pt();
  double lep_2_pt = loose_lepColl.at(1).Pt();
  double el_miniiso_tight = 0.1;
  double mu_miniiso_tight = 0.2;
  double lep_1_ptcone = 0;
  double lep_2_ptcone = 0;
    
  TString hist_name_1, hist_name_2;
  if(loose_lepColl.at(0).LeptonFlavour() == 1){//1st lepton is muon
    hist_name_1 = "SingleMuonTrigger_Dijet_central_all_MUON_SUSY_VETO_central_Nvtx_reweight_events_pt_cone_vs_eta_F";
    lep_1_ptcone = loose_lepColl.at(0).Pt()*(1 + max(0.,(loose_lepColl.at(0).miniRelIso() - mu_miniiso_tight) ) );
  }
  else if (loose_lepColl.at(0).LeptonFlavour() == 2){//1st lepton is electron
    hist_name_1 = "ELECTRON_SUSY_HNPAIR_TIGHTFR_events_pt_vs_eta_TElectron";
    lep_1_ptcone = loose_lepColl.at(0).Pt()*(1 + max(0.,(loose_lepColl.at(0).miniRelIso() - el_miniiso_tight) ) );
  }
  else return 0.;

  if(loose_lepColl.at(1).LeptonFlavour() == 1){//2nd lepton is muon
    hist_name_2 = "SingleMuonTrigger_Dijet_central_all_MUON_SUSY_VETO_central_Nvtx_reweight_events_pt_cone_vs_eta_F";
    lep_2_ptcone = loose_lepColl.at(1).Pt()*(1 + max(0.,(loose_lepColl.at(1).miniRelIso() - mu_miniiso_tight) ) );
  }
  else if (loose_lepColl.at(1).LeptonFlavour() == 2){//2nd lepton is electron
    hist_name_2 = "ELECTRON_SUSY_HNPAIR_TIGHTFR_events_pt_vs_eta_TElectron";
    lep_2_ptcone = loose_lepColl.at(1).Pt()*(1 + max(0.,(loose_lepColl.at(1).miniRelIso() - el_miniiso_tight) ) );
  }
  else return 0.;
  
  double fr1 = 0.;
  double fr2 = 0.;
  
  if(lep_1_pt < 110.){
    int binx = hist_data_driven[hist_name_1]->FindBin(lep_1_ptcone, fabs(loose_lepColl.at(0).Eta()));
    fr1 = float(hist_data_driven[hist_name_1]->GetBinContent(binx));
  }
  else{
    int binx = hist_data_driven[hist_name_1]->FindBin(100., fabs(loose_lepColl.at(0).Eta()));
    fr1 = float(hist_data_driven[hist_name_1]->GetBinContent(binx));
  }
  if(lep_2_pt < 110.){
    int binx = hist_data_driven[hist_name_2]->FindBin(lep_2_ptcone, fabs(loose_lepColl.at(1).Eta()));
    fr2 = float(hist_data_driven[hist_name_2]->GetBinContent(binx));
  }
  else{
    int binx = hist_data_driven[hist_name_2]->FindBin(100., fabs(loose_lepColl.at(1).Eta()));
    fr2 = float(hist_data_driven[hist_name_2]->GetBinContent(binx));
  }

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
  
  if(lep_1_tight && !lep_2_tight){//TL
    return termTL;
  }
  else if(!lep_1_tight && lep_2_tight){//LT
    return termLT;
  }
  else if(!lep_1_tight && !lep_2_tight){//LL
    return termLL;
  }
  else return 0.;//TT
}


double HN_pair_all::GetFRWeight_tri_SB(std::vector<KLepton> loose_lepColl, bool lep_1_tight, bool lep_2_tight, bool lep_3_tight, bool geterr){
  if(loose_lepColl.size() != 3) return 0.;
  
  double lep_1_pt = loose_lepColl.at(0).Pt();
  double lep_2_pt = loose_lepColl.at(1).Pt();
  double lep_3_pt = loose_lepColl.at(2).Pt();

  double el_miniiso_tight = 0.1;
  double mu_miniiso_tight = 0.2;
  
  double lep_1_ptcone = 0;
  double lep_2_ptcone = 0;
  double lep_3_ptcone = 0;

  
  TString hist_name_1, hist_name_2, hist_name_3;
  if(loose_lepColl.at(0).LeptonFlavour() == 1){//1st lepton is muon
    hist_name_1 = "SingleMuonTrigger_Dijet_central_all_MUON_SUSY_VETO_central_Nvtx_reweight_events_pt_cone_vs_eta_F";
    lep_1_ptcone = loose_lepColl.at(0).Pt()*(1 + max(0.,(loose_lepColl.at(0).miniRelIso() - mu_miniiso_tight) ) );
  }
  else if (loose_lepColl.at(0).LeptonFlavour() == 2){//1st lepton is electron
    hist_name_1 = "ELECTRON_SUSY_HNPAIR_TIGHTFR_events_pt_vs_eta_TElectron";
    lep_1_ptcone = loose_lepColl.at(0).Pt()*(1 + max(0.,(loose_lepColl.at(0).miniRelIso() - el_miniiso_tight) ) );
  }
  else return 0.;

  if(loose_lepColl.at(1).LeptonFlavour() == 1){//2nd lepton is muon
    hist_name_2 = "SingleMuonTrigger_Dijet_central_all_MUON_SUSY_VETO_central_Nvtx_reweight_events_pt_cone_vs_eta_F";
    lep_2_ptcone = loose_lepColl.at(1).Pt()*(1 + max(0.,(loose_lepColl.at(1).miniRelIso() - mu_miniiso_tight) ) );
  }
  else if (loose_lepColl.at(1).LeptonFlavour() == 2){//2nd lepton is electron
    hist_name_2 = "ELECTRON_SUSY_HNPAIR_TIGHTFR_events_pt_vs_eta_TElectron";
    lep_2_ptcone = loose_lepColl.at(1).Pt()*(1 + max(0.,(loose_lepColl.at(1).miniRelIso() - el_miniiso_tight) ) );
  }
  else return 0.;

  if(loose_lepColl.at(2).LeptonFlavour() == 1){//3rd lepton is muon
    hist_name_3 = "SingleMuonTrigger_Dijet_central_all_MUON_SUSY_VETO_central_Nvtx_reweight_events_pt_cone_vs_eta_F";
    lep_3_ptcone = loose_lepColl.at(2).Pt()*(1 + max(0.,(loose_lepColl.at(2).miniRelIso() - mu_miniiso_tight) ) );
  }
  else if (loose_lepColl.at(1).LeptonFlavour() == 2){//3rd lepton is electron
    hist_name_3 = "ELECTRON_SUSY_HNPAIR_TIGHTFR_events_pt_vs_eta_TElectron";
    lep_3_ptcone = loose_lepColl.at(2).Pt()*(1 + max(0.,(loose_lepColl.at(2).miniRelIso() - el_miniiso_tight) ) );
  }
  else return 0.;

  double fr1 = 0.;
  double fr2 = 0.;
  double fr3 = 0.;
  
  if(lep_1_pt < 1000){
    int binx = hist_data_driven[hist_name_1]->FindBin(lep_1_ptcone, fabs(loose_lepColl.at(0).Eta()));
    fr1 = float(hist_data_driven[hist_name_1]->GetBinContent(binx));
  }
  else{
    int binx = hist_data_driven[hist_name_1]->FindBin(990., fabs(loose_lepColl.at(0).Eta()));
    fr1 = float(hist_data_driven[hist_name_1]->GetBinContent(binx));
  }
  if(lep_2_pt < 1000){
    int binx = hist_data_driven[hist_name_2]->FindBin(lep_2_ptcone, fabs(loose_lepColl.at(1).Eta()));
    fr2 = float(hist_data_driven[hist_name_2]->GetBinContent(binx));
  }
  else{
    int binx = hist_data_driven[hist_name_2]->FindBin(990., fabs(loose_lepColl.at(1).Eta()));
    fr2 = float(hist_data_driven[hist_name_2]->GetBinContent(binx));
  }
  if(lep_3_pt < 1000){
    int binx = hist_data_driven[hist_name_3]->FindBin(lep_3_ptcone, fabs(loose_lepColl.at(2).Eta()));
    fr3 = float(hist_data_driven[hist_name_3]->GetBinContent(binx));
  }
  else{
    int binx = hist_data_driven[hist_name_3]->FindBin(990., fabs(loose_lepColl.at(2).Eta()));
    fr3 = float(hist_data_driven[hist_name_3]->GetBinContent(binx));
  }
  
  vector<float> a, fr_onlyLoose;
  
  a.push_back( fr1/(1.-fr1) );
  a.push_back( fr2/(1.-fr2) );
  a.push_back( fr3/(1.-fr3) );

  float this_weight=-1.;
  
  if(!lep_1_tight) fr_onlyLoose.push_back(a.at(0));
  if(!lep_2_tight) fr_onlyLoose.push_back(a.at(1));
  if(!lep_3_tight) fr_onlyLoose.push_back(a.at(2));
  
  if(fr_onlyLoose.size() < 1) return 0;
  
  for(unsigned int i=0; i<fr_onlyLoose.size(); i++){
    this_weight *= -fr_onlyLoose.at(i);
  }

  if(!geterr) return this_weight;

}



double HN_pair_all::GetCFWeight_SB(std::vector<KLepton> Leptons){
  double cf1 = 0., cf2 = 0.;
  
  double lep_1_pt = Leptons.at(0).Pt();
  double lep_2_pt = Leptons.at(1).Pt();
  double lep_1_eta = Leptons.at(0).Eta();
  double lep_2_eta = Leptons.at(1).Eta();
  TString hist_name_1, hist_name_2;
   
  if(fabs(lep_1_eta) < 0.8) hist_name_1 = "CF_TRUEpT_IB_num";
  else if(fabs(lep_1_eta) < 1.5) hist_name_1 = "CF_TRUEpT_OB_num";
  else hist_name_1 = "CF_TRUEpT_EC_num";
  
  if(fabs(lep_2_eta) < 0.8) hist_name_2 = "CF_TRUEpT_IB_num";
  else if(fabs(lep_2_eta) < 1.5) hist_name_2 = "CF_TRUEpT_OB_num";
  else hist_name_2 = "CF_TRUEpT_EC_num";
  

  if(lep_1_pt > 50){
    int binx = hist_CF[hist_name_1] -> FindBin(1./lep_1_pt);
    cf1 = float(hist_CF[hist_name_1] -> GetBinContent(binx));
  }
  else{
    int binx = hist_CF[hist_name_1] -> FindBin(0.019);
    cf1 = float(hist_CF[hist_name_1] -> GetBinContent(binx));
  }
  if(lep_2_pt > 50){
    int binx = hist_CF[hist_name_2] -> FindBin(1./lep_2_pt);
    cf2 = float(hist_CF[hist_name_2] -> GetBinContent(binx));
  }
  else{
    int binx = hist_CF[hist_name_2] -> FindBin(0.019);
    cf2 = float(hist_CF[hist_name_2] -> GetBinContent(binx));
  }

  return cf1 + cf2;
}

