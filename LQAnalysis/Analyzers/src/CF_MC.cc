// $Id: CF_MC.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQCF_MC Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "CF_MC.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (CF_MC);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
CF_MC::CF_MC() :  AnalyzerCore(), out_muons(0)  {
  
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("CF_MC");
  
  Message("In CF_MC constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
  MakeCleverHistograms(sighist_mm,"DiMuon");


}


void CF_MC::InitialiseAnalysis() throw( LQError ) {
  
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


void CF_MC::ExecuteEvents()throw( LQError ){

  /// Apply the gen weight 
  if(!isData) weight*=MCweight;
    
  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;
   
  FillCutFlow("NoCut", weight);
  
  
  if(!PassMETFilter()) return;     /// Initial event cuts : 
  FillCutFlow("EventCut", weight);
  
  /// #### CAT::: triggers stored are all HLT_Ele/HLT_DoubleEle/HLT_Mu/HLT_TkMu/HLT_Photon/HLT_DoublePhoton
  
  if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex                                                                              
  
  float pileup_reweight=(1.0);
  if (!k_isdata) {   pileup_reweight = mcdata_correction->PileupWeightByPeriod(eventbase->GetEvent());}
  
  
  std::vector<snu::KElectron> CFelectrons = GetElectrons(true,false, "ELECTRON_HN_TIGHTv4");
  std::vector<snu::KTruth> truthColl= eventbase->GetTruth();
  
  if(CFelectrons.size() > 0){//if we have 1 or more truthmatched electrons
    for(unsigned int i = 0; i < CFelectrons.size(); i++){
      if(fabs(CFelectrons[i].Eta())< 0.8){
	FillHist("CF_RECOpT_IB_den", 1. / CFelectrons[i].Pt(), pileup_reweight, 0., 0.05, 100); 
	FillHist("CF_TRUEpT_IB_den", 1. / truthColl.at(CFelectrons[i].MCTruthIndex()).Pt(), pileup_reweight, 0., 0.05, 100);
	if(MCIsCF(CFelectrons[i])){
	  FillHist("CF_RECOpT_IB_num", 1. / CFelectrons[i].Pt(), pileup_reweight, 0., 0.05, 100);
	  FillHist("CF_TRUEpT_IB_num", 1. / truthColl.at(CFelectrons[i].MCTruthIndex()).Pt(), pileup_reweight, 0., 0.05, 100);
	}//if CF
      }//if IB
      else if(fabs(CFelectrons[i].Eta())< 1.5){
	FillHist("CF_RECOpT_OB_den", 1. / CFelectrons[i].Pt(), pileup_reweight, 0., 0.05, 100);
	FillHist("CF_TRUEpT_OB_den", 1. / truthColl.at(CFelectrons[i].MCTruthIndex()).Pt(), pileup_reweight, 0., 0.05, 100);
        if(MCIsCF(CFelectrons[i])){
	  FillHist("CF_RECOpT_OB_num", 1. / CFelectrons[i].Pt(), pileup_reweight, 0., 0.05, 100);
          FillHist("CF_TRUEpT_OB_num", 1. / truthColl.at(CFelectrons[i].MCTruthIndex()).Pt(), pileup_reweight, 0., 0.05, 100);
        }//if CF 
      }//if OB
      else{
	FillHist("CF_RECOpT_EC_den", 1. / CFelectrons[i].Pt(), pileup_reweight, 0., 0.05, 100);
	FillHist("CF_TRUEpT_EC_den", 1. / truthColl.at(CFelectrons[i].MCTruthIndex()).Pt(), pileup_reweight, 0., 0.05, 100);
	if(MCIsCF(CFelectrons[i])){
	  FillHist("CF_RECOpT_EC_num", 1. / CFelectrons[i].Pt(), pileup_reweight, 0., 0.05, 100);
          FillHist("CF_TRUEpT_EC_num", 1. / truthColl.at(CFelectrons[i].MCTruthIndex()).Pt(), pileup_reweight, 0., 0.05, 100);
        }//if CF
      }//if EC
    }//for electron size, i
  }//if electron size > 0
  
  
   
   
   
   return;
}// End of execute event loop
  


void CF_MC::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void CF_MC::BeginCycle() throw( LQError ){
  
  Message("In begin Cycle", INFO);
  
  //
  //If you wish to output variables to output file use DeclareVariable
  // clear these variables in ::ClearOutputVectors function
  //DeclareVariable(obj, label, treename );
  //DeclareVariable(obj, label ); //-> will use default treename: LQTree
  //  DeclareVariable(out_electrons, "Signal_Electrons", "LQTree");
  //  DeclareVariable(out_muons, "Signal_Muons");

  
  return;
  
}

CF_MC::~CF_MC() {
  
  Message("In CF_MC Destructor" , INFO);
  
}


void CF_MC::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void CF_MC::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this CF_MCCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void CF_MC::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}



