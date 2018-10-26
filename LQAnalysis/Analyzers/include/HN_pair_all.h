#ifndef HN_pair_all_h
#define HN_pair_all_h

#include "AnalyzerCore.h"
#include <TKey.h>

class HN_pair_all : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  HN_pair_all();
  ~HN_pair_all();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();
  void FillCutFlow(TString cut, float w);
  void ExecuteEventFromSyst(std::vector<snu::KElectron> electron_all, std::vector<snu::KMuon> muons_all, std::vector<snu::KJet> jets_all, std::vector<snu::KFatJet> fatjets_all, TString syst_flag);
  void SR(TString channel, bool trigger_pass, double current_weight, std::vector<snu::KJet> jets, std::vector<snu::KFatJet> fatjets, std::vector<snu::KElectron> electrons, std::vector<snu::KElectron> electrons_veto, std::vector<snu::KMuon> muons, std::vector<snu::KMuon> muons_veto, int N_electron, int N_veto_ele, int N_muon, int N_veto_muon);
  void CR_Z_mass(TString channel, bool trigger_pass, double current_weight, std::vector<snu::KJet> jets, std::vector<snu::KFatJet> fatjets, std::vector<snu::KElectron> electrons, std::vector<snu::KElectron> electrons_veto, std::vector<snu::KMuon> muons, std::vector<snu::KMuon> muons_veto, int N_electron, int N_veto_ele, int N_muon, int N_veto_muon);
  void CR_ttbar_dom(TString channel, bool trigger_pass, double current_weight, std::vector<snu::KJet> jets, std::vector<snu::KFatJet> fatjets, std::vector<snu::KElectron> electrons, std::vector<snu::KElectron> electrons_veto, std::vector<snu::KMuon> muons, std::vector<snu::KMuon> muons_veto, int N_electron, int N_veto_ele, int N_muon, int N_veto_muon);
  bool JSFatJetID(snu::KFatJet fatjet);
  bool IsAwayFromFatJet(snu::KJet jet, vector<snu::KFatJet> fatjets);
  vector<snu::KParticle> RecoPairN(std::vector<KLepton> lepptrs, vector<snu::KFatJet> fatjets, vector<snu::KJet> jets);
  std::vector<snu::KElectron> ElectronUsePtScale(std::vector<snu::KElectron> electrons, int up_down);
  std::vector<snu::KMuon> MuonUsePtScale(std::vector<snu::KMuon> muons, int up_down);
  std::vector<snu::KJet> JetUsePtScale(std::vector<snu::KJet> jets, int up_down);
  std::vector<snu::KJet> JetUsePtRes(std::vector<snu::KJet> jets, int up_down);
  std::vector<snu::KFatJet> FatJetUsePtScale(std::vector<snu::KFatJet> fatjets, int up_down);
  std::vector<snu::KFatJet> FatJetUsePtRes(std::vector<snu::KFatJet> fatjets, int up_down);
  std::vector<snu::KElectron> SelectElectronsWithID(std::vector<snu::KElectron> electrons, TString id, double ptmin, double fetamax);
  std::vector<snu::KMuon> SelectMuonsWithID(std::vector<snu::KMuon> muons, TString id, double ptmin, double fetamax);
  std::vector<snu::KJet> SelectJetsWithID(std::vector<snu::KJet> jets, double ptmin, double fetamax);
  std::vector<snu::KFatJet> SelectFatJetsWithID(std::vector<snu::KFatJet> fatjets, double ptmin, double fetamax);

  
  
 private:

  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( HN_pair_all, 1);
};
#endif
