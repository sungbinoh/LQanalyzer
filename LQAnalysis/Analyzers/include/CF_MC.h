#ifndef CF_MC_h
#define CF_MC_h

#include "AnalyzerCore.h"
class CF_MC : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  CF_MC();
  ~CF_MC();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();
 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( CF_MC, 1);
};
#endif
