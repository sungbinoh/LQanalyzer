#ifndef HN_pair_fake_h
#define HN_pair_fake_h

#include "AnalyzerCore.h"
#include <TKey.h>

class HN_pair_fake : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  HN_pair_fake();
  ~HN_pair_fake();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();
  void FillCutFlow(TString cut, float w);
  void Signal_region_1(TString channel, bool trigger_pass, std::vector<snu::KJet> jets, std::vector<snu::KFatJet> fatjets, std::vector<snu::KElectron> electrons, std::vector<snu::KMuon> muons, int N_electron, int N_veto_ele, int N_muon, int N_veto_muon, bool NonPromptRun);
  void Measure_FR_muon(TString tag, TString ID, std::vector<snu::KMuon> loose_mu, bool debugging);
  void Measure_FR_electron(TString tag, TString ID, std::vector<snu::KElectron> loose_el, bool debugging);
  void FillDenAndNum(TString prefix, snu::KMuon muon, double thisweight, bool isTight);
  void FillDenAndNum(TString prefix, snu::KElectron electron, double thisweight, bool isTight);
  bool JSFatJetID(snu::KFatJet fatjet);
  bool IsAwayFromFatJet(snu::KJet jet, vector<snu::KFatJet> fatjets);
  double GetFRWeight_SB(std::vector<snu::KMuon> muon_looseColl, TString tight_ID);
  std::map< TString, TH2D* > hist_Muon_FR_syst;
  std::vector<TString> hist_Muon_FR_syst_name;
 private:

  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( HN_pair_fake, 1);
};
#endif
