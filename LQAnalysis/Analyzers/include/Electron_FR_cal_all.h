#ifndef Electron_FR_cal_all_h
#define Electron_FR_cal_all_h

#include "AnalyzerCore.h"


class Electron_FR_cal_all : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  Electron_FR_cal_all();
  ~Electron_FR_cal_all();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void FillDenAndNum(TString prefix, snu::KMuon muon, double thisweight, bool isTight);
  double GetTriggerWeightByPtRange(TString hltname, vector<double> ptrange, double trigger_safe_pt, std::vector<snu::KMuon> muons, int npfjet50);
      
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();

  
  void RunFakes(TString tag, TString ID, std::vector<snu::KMuon> loosemuons, bool debugging);
  void RunFakes_ID(TString tag, TString ID, std::vector<snu::KMuon> loosemuons, bool debugging);
  void RunFakes_subtraction(TString tag, TString ID, std::vector<snu::KMuon> loosemuons, bool debugging);

  
private:
  
  
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_electrons;
  std::vector<snu::KMuon> out_muons;

  int n_17_pass;
  int n_17_17_jet_pass;
  int n_17_jet_pass;
  TH2D *FR_sampleA;
  double METauto, METphiauto;
  ClassDef ( Electron_FR_cal_all, 1);
};
#endif
