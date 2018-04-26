#ifndef HN_pair_MM_h
#define HN_pair_MM_h

#include "AnalyzerCore.h"
#include <TKey.h>

class HN_pair_MM : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  HN_pair_MM();
  ~HN_pair_MM();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();
  void FillCutFlow(TString cut, float w);
  void Signal_region_1(TString muon_tight_id, TString muon_loose_id, bool trigger_pass, std::vector<snu::KJet> jets, std::vector<snu::KElectron> electrons_veto, bool NonPromptRun);
  void Control_region_1(TString muon_tight_id, TString muon_loose_id, bool trigger_pass, std::vector<snu::KJet> jets, std::vector<snu::KElectron> electrons_veto, bool NonPromptRun);
  std::vector<snu::KJet> Remove_muon_from_jet(std::vector<snu::KJet> jets, std::vector<snu::KMuon> muons);
  snu::KJet Clean_jet_lepton(snu::KJet jet, snu::KMuon muon);
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


  ClassDef ( HN_pair_MM, 1);
};
#endif
