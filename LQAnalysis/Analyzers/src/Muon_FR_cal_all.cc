/***************************************************************************
 * @Project: LQMuon_FR_cal_all Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/


/// Local includes
#include "Muon_FR_cal_all.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"

//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (Muon_FR_cal_all);


/**
 *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
 *
 */
Muon_FR_cal_all::Muon_FR_cal_all() :  AnalyzerCore(),  out_electrons(0) {


  // To have the correct name in the log:                                                                                                                            
  SetLogName("Muon_FR_cal_all");

  Message("In Muon_FR_cal_all constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

}


void Muon_FR_cal_all::InitialiseAnalysis() throw( LQError ) {
  
  /// Initialise histograms
  MakeHistograms();  
  //
  // You can out put messages simply with Message function. Message( "comment", output_level)   output_level can be VERBOSE/INFO/DEBUG/WARNING 
  // You can also use m_logger << level << "comment" << int/double  << LQLogger::endmsg;
  //

  MakeCleverHistograms(sighist_m, "SingleMuonSNU_prompt");
  MakeCleverHistograms(sighist_m, "SingleMuonGent_loose_prompt");
  MakeCleverHistograms(sighist_m, "SingleMuonGent_tight_prompt");
  return;
}


void Muon_FR_cal_all::ExecuteEvents()throw( LQError ){
  
  //cout << "--------------------------------" << endl;
  //cout << "1 : " << eventbase->GetEvent().EventNumber() << endl;
  /*bool SBt_JSf = CheckEventComparison("suoh","JS_fake_mu22119_periodC_SKDoubleMuon_hnfake_cat_v8-0-7_FRCalculator_Mu_dxysig_DILEP",
  "suoh","SB_fake_mu22118_periodC_SKDoubleMuon_hnfake_cat_v8-0-7_Muon_FR_cal_all", false);
  bool SBf_JSt = CheckEventComparison("suoh","JS_fake_mu22119_periodC_SKDoubleMuon_hnfake_cat_v8-0-7_FRCalculator_Mu_dxysig_DILEP",
  "suoh","SB_fake_mu22118_periodC_SKDoubleMuon_hnfake_cat_v8-0-7_Muon_FR_cal_all", true);
  
  if(!SBt_JSf && !SBf_JSt) return;
  cout << "====================================" << endl;
  cout << "RunNumber : " << eventbase->GetEvent().RunNumber() << ", EventNumber : " << eventbase->GetEvent().EventNumber() << endl;
  if(SBt_JSf){
    cout << "(SB && !JS)" << endl;
    FillEventComparisonFile("JSf_SBt_diff");
  }
  if(SBf_JSt){
    cout << "(!SB && JS)" << endl;
    FillEventComparisonFile("JSt_SBf_diff");
  }
  */
  
  bool SBt_JSf = CheckEventComparison("suoh","DY_808_31378_SKDY50plus_cat_v8-0-8_Muon_FR_cal_all",
				      "suoh","DY_807_31377_SKLowStat_DYJets_cat_v8-0-7_Muon_FR_cal_all", false);
  bool SBf_JSt = CheckEventComparison("suoh","DY_808_31378_SKDY50plus_cat_v8-0-8_Muon_FR_cal_all",
				      "suoh","DY_807_31377_SKLowStat_DYJets_cat_v8-0-7_Muon_FR_cal_all", true);
  bool both_tag = !SBt_JSf && !SBf_JSt;
  
  
  //if(eventbase->GetEvent().RunNumber() != 276282 && eventbase->GetEvent().RunNumber() != 276283) return;
  
  //// Initial event cuts
  /// MET FIleters
  if(!PassMETFilter()) return;
  
  /// Require good promary vertex
  if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex
  numberVertices = eventbase->GetEvent().nVertices();
  
  
  bool pass_50 = PassTrigger("HLT_Mu50_v");
  bool pass_27 = PassTrigger("HLT_Mu27_v");
  bool pass_20 = PassTrigger("HLT_Mu20_v");

  FillHist("HLT_Mu50_v_Pass", 0.5, 1. * MCweight, 0., 5., 5);
  if(pass_20) FillHist("HLT_Mu50_v_Pass", 1.5, 1. * MCweight, 0., 5., 5);
  if(pass_27) FillHist("HLT_Mu50_v_Pass", 2.5, 1. * MCweight, 0., 5., 5);
  if(pass_50) FillHist("HLT_Mu50_v_Pass", 3.5, 1. * MCweight, 0., 5., 5);

  
  std::vector<TString> HLTs_single;
  HLTs_single.push_back("HLT_Mu20_v");
  HLTs_single.push_back("HLT_Mu27_v");
  HLTs_single.push_back("HLT_Mu50_v");
  if(!PassTriggerOR(HLTs_single)) return;
  
  std::vector<snu::KMuon> muColl_test = GetMuons("MUON_POG_TIGHT",true);
  std::vector<snu::KElectron> electron_test = GetElectrons(false, true, "ELECTRON_POG_TIGHT");

  
  //cout << "weight : " << weight << ", MCweight : " << MCweight << endl;

  if(!isData) weight*=MCweight;

  std::vector<TString> AllHLTs_el;
  AllHLTs_el.push_back("HLT_Photon22_v");
  AllHLTs_el.push_back("HLT_Photon30_v");
  AllHLTs_el.push_back("HLT_Photon36_v");
  AllHLTs_el.push_back("HLT_Photon50_v");
  AllHLTs_el.push_back("HLT_Photon75_v");
  AllHLTs_el.push_back("HLT_Photon90_v");
  AllHLTs_el.push_back("HLT_Photon120_v");
  AllHLTs_el.push_back("HLT_Photon175_v");

  std::vector<TString> AllHLTs;
  AllHLTs.push_back("HLT_Mu20_v");
  AllHLTs.push_back("HLT_Mu27_v");
  AllHLTs.push_back("HLT_Mu50_v");
  
  std::map< TString, double > HLT_ptmin;
  HLT_ptmin["HLT_Photon22_v"] = 25.;
  HLT_ptmin["HLT_Photon30_v"] = 35.;
  HLT_ptmin["HLT_Photon36_v"] = 40.;
  HLT_ptmin["HLT_Photon50_v"] = 55.;
  HLT_ptmin["HLT_Photon75_v"] = 80.;
  HLT_ptmin["HLT_Photon90_v"] = 100.;
  HLT_ptmin["HLT_Photon120_v"] = 130.;
  HLT_ptmin["HLT_Photon175_v"] = 190.;
  
  
  for(unsigned int i=0; i<AllHLTs.size(); i++){
    TString ThisTrigger = AllHLTs.at(i);
  
    if(PassTrigger(ThisTrigger)){// to check difference between 807 and 808
      if(muColl_test.size() == 2){
	if(muColl_test.at(0).Pt() > 55 && muColl_test.at(1).Pt() > 25.){
	//cout << "1 : " << eventbase->GetEvent().EventNumber() << endl;
	  
	  double m_Z = 91.1876;
	  double m_mumu = (muColl_test.at(0) + muColl_test.at(1) ).M();
	  
	  if(fabs(m_mumu - m_Z) < 10.){
	    //if(SBf_JSt) cout << "SBf_JSt 1 : " << eventbase->GetEvent().EventNumber() << endl;
	    //if(SBt_JSf) cout << "SBt_JSf 1 : " << eventbase->GetEvent().EventNumber() << endl;
	    /*
	      if(both_tag){
	      cout << "1 : " << eventbase->GetEvent().EventNumber() << endl;
	      cout << "1st lep Pt, eta, phi : " << muColl_test.at(0).Pt() << ", " << muColl_test.at(0).Eta() << ", " << muColl_test.at(0).Phi() << endl;
	      cout << "2nd lep Pt, eta, phi : " << muColl_test.at(1).Pt() << ", " << muColl_test.at(1).Eta() << ", " << muColl_test.at(1).Phi() << endl;
	      }
	    */
	    //FillEventComparisonFile("DY_808_");
	    FillHist("m_mumu_" + ThisTrigger, m_mumu, 1. * MCweight, 0., 200., 200);
	    FillHist("1st_lep_PT_" + ThisTrigger, muColl_test.at(0).Pt(), 1. * MCweight, 0., 500., 500);
	    FillHist("2nd_lep_PT_" + ThisTrigger, muColl_test.at(1).Pt(), 1. * MCweight, 0., 500., 500);
	    FillHist("1st_lep_PT_" + ThisTrigger, muColl_test.at(0).Pt(), 1. * MCweight, 0., 500., 500);
	    FillHist("2nd_lep_PT_" + ThisTrigger, muColl_test.at(1).Pt(), 1. * MCweight, 0., 500., 500);
	    FillHist("1st_lep_PT_" + ThisTrigger, muColl_test.at(0).Pt(), 1. * MCweight, 0., 500., 500);
	    FillHist("2nd_lep_PT_" + ThisTrigger, muColl_test.at(1).Pt(), 1. * MCweight, 0., 500., 500);
	    //cout << "weight : " << weight << endl;
	    FillHist("m_mumu_" + ThisTrigger + "_weight", m_mumu, weight, 0., 200., 200);
	    FillHist("1st_lep_PT_" + ThisTrigger + "_weight", muColl_test.at(0).Pt(), weight, 0., 500., 500);
	    FillHist("2nd_lep_PT_" + ThisTrigger + "_weight", muColl_test.at(1).Pt(), weight, 0., 500., 500);
	    FillHist("1st_lep_PT_" + ThisTrigger + "_weight", muColl_test.at(0).Pt(), weight, 0., 500., 500);
	    FillHist("2nd_lep_PT_" + ThisTrigger + "_weight", muColl_test.at(1).Pt(), weight, 0., 500., 500);
	    FillHist("1st_lep_PT_" + ThisTrigger + "_weight", muColl_test.at(0).Pt(), weight, 0., 500., 500);
	    FillHist("2nd_lep_PT_" + ThisTrigger + "_weight", muColl_test.at(1).Pt(), weight, 0., 500., 500);
	    
	  }
	}
      }
    }
  }
  
  for(unsigned int i=0; i<AllHLTs_el.size(); i++){
    TString ThisTrigger = AllHLTs_el.at(i);

    if(PassTrigger(ThisTrigger)){// to check difference between 807 and 808                                                                                                                                                                                                     
      if(electron_test.size() == 2){
        if(electron_test.at(0).Pt() > HLT_ptmin[ThisTrigger] && electron_test.at(1).Pt() > 25.){
	  //cout << "1 : " << eventbase->GetEvent().EventNumber() << endl;
	  
          double m_Z = 91.1876;
          double m_mumu = (electron_test.at(0) + electron_test.at(1) ).M();

          if(fabs(m_mumu - m_Z) < 10.){
	    cout << "fill electrons" << endl;

	    
            FillHist("m_elel_" + ThisTrigger, m_mumu, 1. * MCweight, 0., 200., 200);
            FillHist("1st_el_PT_" + ThisTrigger, electron_test.at(0).Pt(), 1. * MCweight, 0., 500., 500);
            FillHist("2nd_el_PT_" + ThisTrigger, electron_test.at(1).Pt(), 1. * MCweight, 0., 500., 500);
            FillHist("1st_el_PT_" + ThisTrigger, electron_test.at(0).Pt(), 1. * MCweight, 0., 500., 500);
            FillHist("2nd_el_PT_" + ThisTrigger, electron_test.at(1).Pt(), 1. * MCweight, 0., 500., 500);
            FillHist("1st_el_PT_" + ThisTrigger, electron_test.at(0).Pt(), 1. * MCweight, 0., 500., 500);
            FillHist("2nd_el_PT_" + ThisTrigger, electron_test.at(1).Pt(), 1. * MCweight, 0., 500., 500);
            //cout << "weight : " << weight << endl;
            FillHist("m_elel_" + ThisTrigger + "_weight", m_mumu, weight, 0., 200., 200);
            FillHist("1st_el_PT_" + ThisTrigger + "_weight", electron_test.at(0).Pt(), weight, 0., 500., 500);
            FillHist("2nd_el_PT_" + ThisTrigger + "_weight", electron_test.at(1).Pt(), weight, 0., 500., 500);
            FillHist("1st_el_PT_" + ThisTrigger + "_weight", electron_test.at(0).Pt(), weight, 0., 500., 500);
            FillHist("2nd_el_PT_" + ThisTrigger + "_weight", electron_test.at(1).Pt(), weight, 0., 500., 500);
            FillHist("1st_el_PT_" + ThisTrigger + "_weight", electron_test.at(0).Pt(), weight, 0., 500., 500);
            FillHist("2nd_el_PT_" + ThisTrigger + "_weight", electron_test.at(1).Pt(), weight, 0., 500., 500);

          }
        }
      }
    }
  }



  return;
    

  std::vector<snu::KElectron> elColl = GetElectrons("ELECTRON_HN_VETO");
  if(elColl.size() != 0) return; // no e 

  
  
  
    
  //cout << weight << endl;
  std::vector<snu::KMuon> muColl = GetMuons("MUON_SUSY_VETO",true);
  CorrectedMETRochester(muColl);
  snu::KEvent Evt = eventbase->GetEvent();
  METauto = Evt.PFMET();
  METphiauto = Evt.METPhi();
  //cout << "after correcting" << endl;
  //cout << "MET : " << Evt.PFMET() << ", METphi : " << Evt.METPhi() << endl;
  
  std::vector<snu::KMuon> loose_mu_noMiniiso = GetMuons("MUON_SUSY_VETO",true);
  std::vector<snu::KMuon> loose_mu_central;
  for(unsigned int i = 0; i < loose_mu_noMiniiso.size(); i++){
    if(loose_mu_noMiniiso.at(i).PFRelMiniIsoRho() < 0.6) loose_mu_central.push_back(loose_mu_noMiniiso.at(i));
  }

  std::vector<snu::KMuon> loosest_mu = GetMuons("MUON_SUSY_VETO",true);
  
  double reliso_cut[] = {0.6, 0.3, 0.4};
  double dxy_cut[] = {999., 0.1, 1.0, 10.};
  double dz_cut[] = {999., 1.0, 10};
  double dxy_sip_cut[] = {999., 8., 10., 15.};

  int reliso_i[]  = {1, 2, 0, 0, 0, 0, 0, 0, 0, 0};
  int dxy_i[]     = {0, 0, 1, 2, 3, 0, 0, 0, 0, 0};
  int dz_i[]      = {0, 0, 0, 0, 0, 1, 2, 0, 0, 0};
  int dxy_sip_i[] = {0, 0, 0, 0, 0, 0, 0, 1, 2, 3};

  TString cut_string[] = {"iso_p3", "iso_p4", "dxy_p1", "dxy_3", "dxy_10", "dz_1", "dz_10", "SIP_8", "SIP_10", "SIP_15"};
  
  for(int i = 0; i < 10; i ++){
    std::vector<snu::KMuon> this_mu;
    for(int j = 0; j < loosest_mu.size(); j++){
      if(loosest_mu.at(j).PFRelMiniIsoRho() < reliso_cut[reliso_i[i]] && fabs(loosest_mu.at(j).dXY()) < dxy_cut[dxy_i[i]] && fabs(loosest_mu.at(j).dZ()) < dz_cut[dz_i[i]]
         && fabs(loosest_mu.at(j).dXYSig()) < dxy_sip_cut[dxy_sip_i[i]]){
	this_mu.push_back(loosest_mu.at(j));
      }//if pass id cuts
    }//for j
    //cout << "RunFakes_ID" << endl;
    RunFakes_ID("MUON_SUSY_TIGHT", "MUON_SUSY_VETO_" + cut_string[i], this_mu, false);
  }//for i
  
  //cout << "RunFakes start" << endl;
  RunFakes("MUON_SUSY_TIGHT", "MUON_SUSY_VETO_central", loose_mu_central, false);
  
  
}

//function to run for loose ID systematics : FR only for away jet pt 40 GeV
void Muon_FR_cal_all::RunFakes_ID(TString tight_ID, TString loose_ID, std::vector<snu::KMuon> loose_mu, bool debugging){
  std::map< TString, std::vector<double> > HLT_ptrange; //defining muon pt range for each single muon trigger
  std::map< TString, double > HLT_ptmin;
  HLT_ptrange["HLT_Mu20_v"].push_back(35.);
  HLT_ptrange["HLT_Mu20_v"].push_back(45.);
  HLT_ptmin["HLT_Mu20_v"] = 23.;

  HLT_ptrange["HLT_Mu27_v"].push_back(45.);
  HLT_ptrange["HLT_Mu27_v"].push_back(80.);
  HLT_ptmin["HLT_Mu27_v"] = 30.;

  HLT_ptrange["HLT_Mu50_v"].push_back(80.);
  HLT_ptrange["HLT_Mu50_v"].push_back(9999.);
  HLT_ptmin["HLT_Mu50_v"] = 55.;
  
  std::vector<TString> AllHLTs;
  for(std::map< TString, std::vector<double> >::iterator it=HLT_ptrange.begin(); it!=HLT_ptrange.end(); it++){
    AllHLTs.push_back(it->first);
  }
  
  bool DijetFake = std::find(k_flags.begin(), k_flags.end(), "DijetFake") != k_flags.end(); //flag for QCD MC fr cal
  bool DijetPrompt = std::find(k_flags.begin(), k_flags.end(), "DijetPrompt") != k_flags.end(); //flag for MC subtraction

  std::vector<snu::KJet> jets = GetJets("JET_HN");
  int For_HLT_Mu3_PFJet40_v = 0;
  int n_bjets = 0;
  for(unsigned int i = 0; i < jets.size(); i++){
    if( IsBTagged(jets.at(i), snu::KJet::CSVv2, snu::KJet::Medium) ) n_bjets++;
    if( jets.at(i).Pt() > 50. ) For_HLT_Mu3_PFJet40_v++;
  }

  float pileup_reweight=(1.0);
  if(!k_isdata){
    pileup_reweight = mcdata_correction->CatPileupWeight(eventbase->GetEvent(),0);
  }
  
  if(debugging) cout << "jets.size() : " << jets.size() << endl;
  if(jets.size() < 1) return; // Njet >= 1
  
  std::vector<snu::KMuon> hnloose;//keep only prompt muon for MC subtraction, keep only fake muon for QCD fr cal
  hnloose.clear();

  if(debugging) cout << "loose_mu.size() : " << loose_mu.size() << endl;
  if(loose_mu.size() != 1) return;

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
	  //if(loose_mu.at(i).GetType() == 1 || loose_mu.at(i).GetType() == 8){  
	  hnloose.push_back( loose_mu.at(i) );
	}
      }
      else{
        return;
      }
    }
  }
  
  if(debugging) cout << "hnloose.size() : " << hnloose.size() << endl;
  if(hnloose.size() != 1) return;
  
  snu::KMuon muon = hnloose.at(0);
  TLorentzVector metvec;
  metvec.SetPtEtaPhiE( METauto, 0, METphiauto, METauto);
  double MTval = AnalyzerCore::MT( muon, metvec );
  if(debugging) cout << "METauto : " << METauto << ", MTval : " << MTval << endl;
  
  TString current_trigger;
  double weight_by_pt(0.);//calculate weight for this event
  
  for(std::map< TString, std::vector<double> >::iterator it=HLT_ptrange.begin(); it!=HLT_ptrange.end(); it++){
    double tmp = GetTriggerWeightByPtRange(it->first, it->second, HLT_ptmin[it->first], hnloose, For_HLT_Mu3_PFJet40_v);
    if(tmp > 0.1){
      weight_by_pt = tmp;
      current_trigger = it->first;
      break;
    }
  }
  double Nvtx_reweight = 1.;
  snu::KEvent Event = eventbase->GetEvent();
  int nvtx = Event.nVertices();
  if(!k_isdata && weight_by_pt > 0.1){
    Nvtx_reweight = mcdata_correction -> GetVtxReweight(current_trigger, nvtx);
  }

  //if(DijetPrompt && (METauto > 20 || MTval > 25) ) return;//give event cut for MC subtraction case
  if(DijetPrompt && (METauto > 80 || MTval > 25) ) return;
  

  double this_weight = weight*pileup_reweight;
  bool IsThisTight = PassID( muon, tight_ID) && (muon.PFRelMiniIsoRho() < 0.2);
  double conept = muon.Pt()*(1 + max(0.,(muon.PFRelMiniIsoRho() - 0.2) ) );

  
  double AwayjetPt = 40.;
  bool histfilled = false;
  for(unsigned int j = 0; j < jets.size(); j++){//for loop over all jets in the event
    if(histfilled) break;
    snu::KJet jet = jets.at(j);
    if(debugging) cout << "AwayjetPt : " << AwayjetPt << ", jet.Pt(): " << jet.Pt() << ". jets.size() : " << jets.size() << endl;
    if( jet.Pt() < AwayjetPt ) continue;

    double dPhi = muon.DeltaPhi( jet );

    bool UseEvent = false;
    if(debugging) cout << "dPhi : " << dPhi << ". (jet.Pt()/muon.Pt() : " << (jet.Pt()/muon.Pt()) << endl;
    UseEvent = (dPhi > 2.5) && (jet.Pt()/muon.Pt() > 1.);
    
    if(debugging) cout <<  "UseEvent : " << UseEvent << endl;
    
    if( UseEvent ){
      TString current_hist_name = "SingleMuonTrigger_Dijet_Awayjet_" + TString::Itoa(AwayjetPt,10)+"_all_" + loose_ID;
      //FillDenAndNum("SingleMuonTrigger_Dijet_Awayjet_"+TString::Itoa(AwayjetPt,10)+"_all_" + loose_ID, muon, this_weight * weight_by_pt, IsThisTight);
      double MuTrkEffSF =  mcdata_correction->MuonTrackingEffScaleFactor(hnloose);
            
      if(IsThisTight){
	double muon_SF = mcdata_correction->MuonScaleFactor("MUON_SUSY_TIGHT", hnloose, 0);
	double trigger_SF = 1.;
	double SFs = muon_SF * trigger_SF * MuTrkEffSF;
	if(k_isdata) SFs = 1.;
	//FillDenAndNum(current_hist_name + "_Nvtx_reweight", muon, this_weight * weight_by_pt * Nvtx_reweight * SFs, IsThisTight);
	FillDenAndNum(current_hist_name + "_Nvtx_reweight", muon, this_weight * weight_by_pt * SFs, IsThisTight);

      }// if istight = apply SFs for tight lepton
      else FillDenAndNum(current_hist_name + "_Nvtx_reweight", muon, this_weight * weight_by_pt * MuTrkEffSF, IsThisTight);

      histfilled = true;
    }//end of if(UseEvent)
  }//end for jets
}
////////////////////////////////////////////////////////RunFakes_ID ends

void Muon_FR_cal_all::RunFakes(TString tight_ID, TString loose_ID, std::vector<snu::KMuon> loose_mu, bool debugging){
  
  std::map< TString, std::vector<double> > HLT_ptrange; //defining muon pt range for each single muon trigger
  std::map< TString, double > HLT_ptmin;
  HLT_ptrange["HLT_Mu20_v"].push_back(35.);
  HLT_ptrange["HLT_Mu20_v"].push_back(45.);
  HLT_ptmin["HLT_Mu20_v"] = 23.;

  HLT_ptrange["HLT_Mu27_v"].push_back(45.);
  HLT_ptrange["HLT_Mu27_v"].push_back(80.);
  HLT_ptmin["HLT_Mu27_v"] = 30.;

  HLT_ptrange["HLT_Mu50_v"].push_back(80.);
  HLT_ptrange["HLT_Mu50_v"].push_back(9999.);
  HLT_ptmin["HLT_Mu50_v"] = 55.;
  
  std::vector<TString> AllHLTs;
  for(std::map< TString, std::vector<double> >::iterator it=HLT_ptrange.begin(); it!=HLT_ptrange.end(); it++){
    AllHLTs.push_back(it->first);
  }
  
  bool DijetFake = std::find(k_flags.begin(), k_flags.end(), "DijetFake") != k_flags.end(); //flag for QCD MC fr cal
  bool DijetPrompt = std::find(k_flags.begin(), k_flags.end(), "DijetPrompt") != k_flags.end(); //flag for MC subtraction
  
  std::vector<snu::KJet> jets = GetJets("JET_HN");
  int For_HLT_Mu3_PFJet40_v = 0;
  int n_bjets = 0;
  for(unsigned int i = 0; i < jets.size(); i++){
    if( IsBTagged(jets.at(i), snu::KJet::CSVv2, snu::KJet::Medium) ) n_bjets++;
    if( jets.at(i).Pt() > 50. ) For_HLT_Mu3_PFJet40_v++;
  }

  float pileup_reweight=(1.0);
  if(!k_isdata){
    pileup_reweight = mcdata_correction->CatPileupWeight(eventbase->GetEvent(),0);
  }
  
  //normalize test
  if(DijetPrompt){
    double m_Z = 91.1876;
    
    std::vector< snu::KMuon > tightmuons_noiso = GetMuons("MUON_SUSY_TIGHT", true, 55., 2.4);
    std::vector< snu::KMuon > tightmuons;
    if(tightmuons_noiso.size() > 0){
      for(unsigned int i_mu_noiso = 0; i_mu_noiso < tightmuons_noiso.size(); i_mu_noiso++){
	if(tightmuons_noiso.at(i_mu_noiso).PFRelMiniIsoRho() < 0.2) tightmuons.push_back(tightmuons_noiso.at(i_mu_noiso)) ;
      }
    }
    
    //std::vector< snu::KMuon > tightmuons = GetMuons("MUON_HN_TIGHT", true, 55., 2.4);

    
    double muon_id_iso_sf = 1.;
    double muon_id_iso_sf_syst[2] = {1.};
    double MuTrkEffSF =  1.;
    if(!k_isdata){
      muon_id_iso_sf = mcdata_correction->MuonScaleFactor("MUON_SUSY_TIGHT", tightmuons, 0);
      muon_id_iso_sf_syst[0] = mcdata_correction->MuonScaleFactor("MUON_SUSY_TIGHT", tightmuons, 1);
      muon_id_iso_sf_syst[1] = mcdata_correction->MuonScaleFactor("MUON_SUSY_TIGHT", tightmuons, -1);
      MuTrkEffSF = mcdata_correction->MuonTrackingEffScaleFactor(tightmuons);
    }
    TString muon_id_iso_sf_string[2] = {"ID_Up", "ID_Down"};
    
    
    TLorentzVector metvec;
    metvec.SetPtEtaPhiE( METauto, 0, METphiauto, METauto );
    
    for(unsigned int i=0; i<AllHLTs.size(); i++){
      TString ThisTrigger = AllHLTs.at(i);
      snu::KEvent Event = eventbase->GetEvent();
      int nvtx = Event.nVertices();
      double Nvtx_reweight[2] = {1.};
      TString Nvtx_reweight_string[2] = {"central", "Liverpool"};
      if(!k_isdata){
	//Nvtx_reweight[0] = mcdata_correction -> GetVtxReweight(ThisTrigger, nvtx);
      }//!is_Data, setting Nvtx reweight values
      
      if( PassTrigger(ThisTrigger) && (loose_mu.size()!=0) ){
	if(HLT_ptmin[ThisTrigger] > loose_mu.at(0).Pt()) continue;
	double triggerweight = 0.;
	if(ThisTrigger.Contains("Mu20") ) triggerweight = 140.315;
	else if(ThisTrigger.Contains("Mu27") ) triggerweight = 250.252;
        else if(ThisTrigger.Contains("Mu50") ) triggerweight = 35863.308;
	else triggerweight = 0.;
	
	double this_weight = weight * muon_id_iso_sf * MuTrkEffSF * triggerweight * pileup_reweight;
        double weight_noID = weight * MuTrkEffSF * triggerweight * pileup_reweight;
	
	if(k_isdata){
	  this_weight = 1.;
	  weight_noID = 1.;
	}
	if(ThisTrigger.Contains("PFJet40")){
          if(For_HLT_Mu3_PFJet40_v==0) this_weight = 0.;
        }
	
	snu::KEvent Event = eventbase->GetEvent();
	int nvtx = Event.nVertices();
	
	//==== 1) Z-Peak
        if(tightmuons.size()==2){
          double mll = (tightmuons.at(0)+tightmuons.at(1)).M();
          if( (tightmuons.at(0).Pt() > 20.) && (tightmuons.at(1).Pt() > 10.) && (fabs(mll-m_Z) < 10.) ){
            FillHist(ThisTrigger+"_ZPeak_mll", mll, this_weight, 0., 200., 200);
            FillHist(ThisTrigger+"_ZPeak_leadpt", tightmuons.at(0).Pt(), this_weight, 0., 500., 500);
            FillHist(ThisTrigger+"_ZPeak_subleadpt", tightmuons.at(1).Pt(), this_weight, 0., 500., 500);
	    FillHist(ThisTrigger+"_ZPeak_Nvtx", nvtx, this_weight, 0., 70., 70);
	    FillHist(ThisTrigger+"_ZPeak_PFMET", METauto, this_weight, 0., 500., 500);
	    
	    for(int Nvtx_i = 0; Nvtx_i < 1; Nvtx_i++){
	      FillHist(ThisTrigger+"_ZPeak_mll_Nvtx_reweight_" + Nvtx_reweight_string[Nvtx_i], mll, this_weight * Nvtx_reweight[Nvtx_i], 0., 200., 200);
	      FillHist(ThisTrigger+"_ZPeak_leadpt_Nvtx_reweight_" + Nvtx_reweight_string[Nvtx_i], tightmuons.at(0).Pt(), this_weight * Nvtx_reweight[Nvtx_i], 0., 500., 500);
	      FillHist(ThisTrigger+"_ZPeak_subleadpt_Nvtx_reweight_" + Nvtx_reweight_string[Nvtx_i], tightmuons.at(1).Pt(), this_weight * Nvtx_reweight[Nvtx_i], 0., 500., 500);
	      FillHist(ThisTrigger+"_ZPeak_Nvtx_Nvtx_reweight_" + Nvtx_reweight_string[Nvtx_i], nvtx, this_weight * Nvtx_reweight[Nvtx_i], 0., 70., 70);
	      FillHist(ThisTrigger+"_ZPeak_PFMET_Nvtx_reweight_" + Nvtx_reweight_string[Nvtx_i], METauto, this_weight * Nvtx_reweight[Nvtx_i], 0., 500., 500);
	    }
	    for(int ID_i = 0; ID_i < 2; ID_i++){
	      FillHist(ThisTrigger+"_ZPeak_mll_Nvtx_reweight_" + muon_id_iso_sf_string[ID_i], mll, weight_noID * muon_id_iso_sf_syst[ID_i], 0., 200., 200);
              FillHist(ThisTrigger+"_ZPeak_leadpt_Nvtx_reweight_" + muon_id_iso_sf_string[ID_i], tightmuons.at(0).Pt(), weight_noID * muon_id_iso_sf_syst[ID_i], 0., 500., 500);
              FillHist(ThisTrigger+"_ZPeak_subleadpt_Nvtx_reweight_" + muon_id_iso_sf_string[ID_i], tightmuons.at(1).Pt(), weight_noID * muon_id_iso_sf_syst[ID_i], 0., 500., 500);
              FillHist(ThisTrigger+"_ZPeak_Nvtx_Nvtx_reweight_" + muon_id_iso_sf_string[ID_i], nvtx, weight_noID * muon_id_iso_sf_syst[ID_i], 0., 70., 70);
              FillHist(ThisTrigger+"_ZPeak_PFMET_Nvtx_reweight_" + muon_id_iso_sf_string[ID_i], METauto, weight_noID * muon_id_iso_sf_syst[ID_i], 0., 500., 500);
	    }
	    double MTval = AnalyzerCore::MT( loose_mu.at(0), metvec );
	    if( (MTval > 60.) && (MTval < 100.)){//MT cut applied Z region
	      FillHist(ThisTrigger+"_ZPeak_MTcut_mll", mll, this_weight, 0., 200., 200);
	      FillHist(ThisTrigger+"_ZPeak_MTcut_leadpt", tightmuons.at(0).Pt(), this_weight, 0., 500., 500);
	      FillHist(ThisTrigger+"_ZPeak_MTcut_subleadpt", tightmuons.at(1).Pt(), this_weight, 0., 500., 500);
	      FillHist(ThisTrigger+"_ZPeak_MTcut_Nvtx", nvtx, this_weight, 0., 70., 70);
	      FillHist(ThisTrigger+"_ZPeak_MTcut_PFMET", METauto, this_weight, 0., 500., 500);
	    }//MT cut
	  }
	}
	
	//==== 2) W
	if(loose_mu.size() ==1){
	  double MTval = AnalyzerCore::MT( loose_mu.at(0), metvec );
	  double current_conept = loose_mu.at(0).Pt()*(1 + max(0.,(loose_mu.at(0).PFRelMiniIsoRho()-0.2) ) );
	  
	  if((MTval > 60.) && (MTval < 100.)){
	    FillHist(ThisTrigger + "_Ptcone_MET_inclusive_loose", current_conept, weight_noID, 0., 200., 200);
	    FillHist(ThisTrigger + "_Eta_MET_inclusive_loose", loose_mu.at(0).Eta(),  weight_noID, -3., 3., 100);
	    FillHist(ThisTrigger + "_MET_MET_inclusive_tight_loose", METauto, weight_noID, 0., 500., 500);
	    FillHist(ThisTrigger + "_MT_MET_inclusive_tight_loose", MTval, weight_noID, 0., 500., 500);
	    
	    bool IsThisTight = PassID( loose_mu.at(0), tight_ID);
	    if(IsThisTight){
	      FillHist(ThisTrigger + "_Ptcone_MET_inclusive_tight", current_conept, this_weight, 0., 200., 200);
	      FillHist(ThisTrigger + "_Eta_MET_inclusive_tight", loose_mu.at(0).Eta(),  this_weight, -3., 3., 100);
	      FillHist(ThisTrigger + "_MET_MET_inclusive_tight", METauto, this_weight, 0., 500., 500);
	      FillHist(ThisTrigger + "_MT_MET_inclusive_tight", MTval, this_weight, 0., 500., 500);
	    }
	  }	
	}
	
	if(tightmuons.size()==1){
          double MTval = AnalyzerCore::MT( tightmuons.at(0), metvec );
          FillHist(ThisTrigger + "MT_MET_MT_inclusive_tight", MTval, this_weight, 0., 500., 500);
	  
	  if( (METauto > 40.) && (MTval > 60.) && (MTval < 100.) && (tightmuons.at(0).Pt() > 20.)){
	    FillHist(ThisTrigger+"_W_John_PFMET", METauto, this_weight, 0., 500., 500);
            FillHist(ThisTrigger+"_W_John_MT", MTval, this_weight, 0., 500., 500);
            FillHist(ThisTrigger+"_W_John_leadpt", tightmuons.at(0).Pt(), this_weight, 0., 500., 500);
	    FillHist(ThisTrigger+"_W_John_Nvtx", nvtx, this_weight, 0., 70., 70);

	    for(int Nvtx_i = 0; Nvtx_i < 1; Nvtx_i++){
	      FillHist(ThisTrigger+"_W_John_PFMET_Nvtx_reweight_" + Nvtx_reweight_string[Nvtx_i], METauto, this_weight * Nvtx_reweight[Nvtx_i], 0., 500., 500);
	      FillHist(ThisTrigger+"_W_John_MT_Nvtx_reweight_" + Nvtx_reweight_string[Nvtx_i], MTval, this_weight * Nvtx_reweight[Nvtx_i], 0., 500., 500);
	      FillHist(ThisTrigger+"_W_John_leadpt_Nvtx_reweight_" + Nvtx_reweight_string[Nvtx_i], tightmuons.at(0).Pt(), this_weight * Nvtx_reweight[Nvtx_i], 0., 500., 500);
	      FillHist(ThisTrigger+"_W_John_Nvtx_Nvtx_reweight_" + Nvtx_reweight_string[Nvtx_i], nvtx, this_weight * Nvtx_reweight[Nvtx_i], 0., 70., 70);
	    }
	    for(int ID_i = 0; ID_i < 2; ID_i++){
	      FillHist(ThisTrigger+"_W_John_PFMET_Nvtx_reweight_" + muon_id_iso_sf_string[ID_i], METauto, weight_noID * muon_id_iso_sf_syst[ID_i], 0., 500., 500);
              FillHist(ThisTrigger+"_W_John_MT_Nvtx_reweight_" + muon_id_iso_sf_string[ID_i], MTval, weight_noID * muon_id_iso_sf_syst[ID_i], 0., 500., 500);
              FillHist(ThisTrigger+"_W_John_leadpt_Nvtx_reweight_" + muon_id_iso_sf_string[ID_i], tightmuons.at(0).Pt(), weight_noID * muon_id_iso_sf_syst[ID_i], 0., 500., 500);
              FillHist(ThisTrigger+"_W_John_Nvtx_Nvtx_reweight_" + muon_id_iso_sf_string[ID_i], nvtx, weight_noID * muon_id_iso_sf_syst[ID_i], 0., 70., 70);
	    }
	  }	      
        }
      }
    }
  }
  //normaize test ends
  //if(debugging) cout << "trigger" << endl;
    
  if(debugging) cout << "jets.size() : " << jets.size() << endl;
  if(jets.size() < 1) return; // Njet >= 1
  
  std::vector<snu::KMuon> hnloose;//keep only prompt muon for MC subtraction, keep only fake muon for QCD fr cal
  hnloose.clear();
 
  if(debugging) cout << "loose_mu.size() : " << loose_mu.size() << endl;
  if(loose_mu.size() != 1) return;
  
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
	  //if(loose_mu.at(i).GetType() ==1 || loose_mu.at(i).GetType() == 8){
	  hnloose.push_back( loose_mu.at(i) );
	}
      }
      else{
	return;
      }
    }
  }
    
  if(debugging) cout << "hnloose.size() : " << hnloose.size() << endl; 
  if(hnloose.size() != 1) return;
  
  FillHist("cutflow_22", 2.5, 1., 0., 10., 10);
  
  snu::KMuon muon = hnloose.at(0);
  TLorentzVector metvec;
  metvec.SetPtEtaPhiE( METauto, 0, METphiauto, METauto);
  double MTval = AnalyzerCore::MT( muon, metvec );
  if(debugging) cout << "METauto : " << METauto << ", MTval : " << MTval << endl;
  
  TString current_trigger;
  double weight_by_pt(0.);//calculate weight for this event                                                                                                                        
  for(std::map< TString, std::vector<double> >::iterator it=HLT_ptrange.begin(); it!=HLT_ptrange.end(); it++){
    double tmp = GetTriggerWeightByPtRange(it->first, it->second, HLT_ptmin[it->first], hnloose, For_HLT_Mu3_PFJet40_v);
    if(tmp > 0.1){
      weight_by_pt = tmp;
      current_trigger = it->first;
      break;
    }
  }
  double Nvtx_reweight = 1.;
  snu::KEvent Event = eventbase->GetEvent();
  int nvtx = Event.nVertices();
  if(!k_isdata && weight_by_pt > 0.1){
    Nvtx_reweight = mcdata_correction -> GetVtxReweight(current_trigger, nvtx);
  }
  
  if(DijetPrompt && (METauto > 80 || MTval > 25) ) return;//give event cut for MC subtraction case
  //test for MET , MT cut on FR
    
  double this_weight = weight*pileup_reweight;
  
  bool IsThisTight = PassID( muon, tight_ID) && (muon.PFRelMiniIsoRho() < 0.2);
  double conept = muon.Pt()*(1 + max(0.,(muon.PFRelMiniIsoRho() - 0.2) ) );
  
  int awayjet_syst_n = 11;
  double awayjet_pt[] = {40., 20., 30., 60}; //away jet pt cuts for systematics, central value is 40 GeV 
  double awayjet_dphi[] = {2.5, 1.5, 2.0, 3.0};//away jet dphi cut for systematics, central value is 2.5
  double pj_over_pl[] = {1.0, 0.8, 1.2, 1.5, 2.0};//pt(j)/pt(l) cut for systematics, central value is 1.5
  int awayjet_pt_i[] =   {1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0};
  int awayjet_dphi_i[] = {0, 0, 0, 1, 2, 3, 0, 0, 0, 0, 0};
  int pj_over_pl_i[] =   {0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 0};
  //TString awayjet_syst_str[] = {"awayjet_pt_20", "awayjet_pt_30", "awayjet_pt_60", "awayjet_dphi_1p5", "awayjet_dphi_2", "awayjet_dphi_3", "pj_over_pl_1", "pj_over_pl_2",
  //"central"};
  TString awayjet_syst_str[] = {"awayjet_pt_20", "awayjet_pt_30", "awayjet_pt_60", "awayjet_dphi_1p5", "awayjet_dphi_2", "awayjet_dphi_3",
                                "pj_over_pl_0p8", "pj_over_pl_1p2", "pj_over_pl_1p5", "pj_over_pl_2",
                                "central"};
  
  if(debugging) cout << "start awayjet for loop" << endl;
  //awayjet systematic for loop starts
  for(int i = 0; i < awayjet_syst_n; i ++){
    //i = 8; //run only for central
    //cuts for dijet topology, and string of current syst.
    double Awayjet_Pt_cut = awayjet_pt[ awayjet_pt_i[i] ];
    double Awayjet_dphi_cut = awayjet_dphi[ awayjet_dphi_i[i] ];
    double pj_over_pl_cut = pj_over_pl[ pj_over_pl_i[i] ];
    TString current_syst_str = awayjet_syst_str[i];
    
    bool histfilled = false; 
    
    for(unsigned int j = 0; j < jets.size(); j++){//for loop all jets
      if(histfilled) break;
      snu::KJet jet = jets.at(j);
      if(debugging) cout << "AwayjetPt : " << Awayjet_Pt_cut << ", jet.Pt(): " << jet.Pt() << ". jets.size() : " << jets.size() << endl;
      if( jet.Pt() < Awayjet_Pt_cut ) continue;
      double dPhi = muon.DeltaPhi( jet );
      
      bool UseEvent = false;
      if(debugging) cout << "dPhi : " << dPhi << ". (jet.Pt()/muon.Pt() : " << (jet.Pt()/muon.Pt()) << endl;
      UseEvent = (dPhi > Awayjet_dphi_cut) && (jet.Pt()/muon.Pt() > pj_over_pl_cut);
      
      if(i == 10){//if central
	FillHist("dPhi_jet_lepton_loose", dPhi, this_weight * weight_by_pt, 0., 3.5, 350);
	FillHist("Pt_jet_over_Pt_lepton_loose", jet.Pt()/muon.Pt(), this_weight * weight_by_pt, 0., 5., 500);
	if(IsThisTight){
	  FillHist("dPhi_jet_lepton_tight", dPhi, this_weight * weight_by_pt, 0., 3.5, 350);
	  FillHist("Pt_jet_over_Pt_lepton_tight", jet.Pt()/muon.Pt(), this_weight * weight_by_pt, 0., 5., 500);
	}
      }
      
      if(debugging) cout <<  "UseEvent : " << UseEvent << endl;
      
      if( UseEvent ){//if pass cuts
	TString current_hist_name = "SingleMuonTrigger_Dijet_" + current_syst_str + "_all_" + loose_ID;
	double MuTrkEffSF =  mcdata_correction->MuonTrackingEffScaleFactor(hnloose);
	
	if(i == 10){//if central
	  FillHist("MET_VS_MT_UseEvent_loose", METauto, MTval, this_weight * weight_by_pt * MuTrkEffSF, 0., 200., 20., 0., 200., 20);
	  FillHist("MET_VS_Ptcone_UseEvent_loose", METauto, conept, this_weight * weight_by_pt * MuTrkEffSF, 0., 200., 20., 0., 200., 20);
	  FillHist("MET_VS_Eta_UseEvent_loose", METauto, muon.Eta(), this_weight * weight_by_pt * MuTrkEffSF, 0., 200., 20., -3., 3., 100);
	  FillHist("Ptcone_UseEvent_loose", conept, this_weight * weight_by_pt * MuTrkEffSF, 0., 200., 200);
	  FillHist("Eta_UseEvent_loose", muon.Eta(), this_weight * weight_by_pt * MuTrkEffSF, -3., 3., 100);
	}
	
	if(IsThisTight){
          double muon_SF = mcdata_correction->MuonScaleFactor("MUON_HN_TIGHT", hnloose, 0);
	  double trigger_SF = 1.;
	  double SFs = muon_SF * trigger_SF * MuTrkEffSF;
          if(k_isdata) SFs = 1.;
          //FillDenAndNum(current_hist_name + "_Nvtx_reweight", muon, this_weight * weight_by_pt * Nvtx_reweight * SFs, IsThisTight);
	  FillDenAndNum(current_hist_name + "_Nvtx_reweight", muon, this_weight * weight_by_pt * SFs, IsThisTight);
	  
	  snu::KEvent Evt_current = eventbase->GetEvent();
	  
	  if( i == 11){
	    cout << "RunNumber = " << Evt_current.RunNumber() << endl;
	    cout << "EventNumber = " << Evt_current.EventNumber() << endl;
	    cout << "trigger = " << current_trigger << endl;
	    cout << "muon.Pt() = " << muon.Pt() << endl;
	    cout << "muon.RelIso04() = " << muon.PFRelMiniIsoRho() << endl;
	    cout << "muon.dXYSig() = " << muon.dXYSig() << endl;
	    cout << "jet.Pt() = " << jet.Pt() << endl;
	    cout << "AwayjetPt = " << Awayjet_Pt_cut << endl;
	    cout << "dPhi = " << dPhi << endl;
	    cout << "jet.Pt()/muon.Pt() = " << jet.Pt()/muon.Pt() << endl;
	    cout << "METauto = " << METauto << endl;
	    cout << "MTval = " << MTval << endl;
	    cout << "this_weight = " << this_weight << endl;
	    //cout << "trigger_additional_sf = " << trigger_additional_sf << endl;
	    cout << "muon_scale_sf = " << SFs<< endl;
	    cout << "weight_by_pt_cone = " << weight_by_pt << endl;
	    cout << "weight = " << weight << endl;
	    cout << "total weight : " << this_weight * weight_by_pt * SFs << endl;
	  }
	  
	  if(i == 10){//if central
	    FillHist("MET_VS_MT_UseEvnet_tight", METauto, MTval, this_weight * weight_by_pt * SFs, 0., 200., 20., 0., 200., 20);
	    FillHist("MET_VS_Ptcone_UseEvnet_tight", METauto, conept, this_weight * weight_by_pt * SFs, 0., 200., 20., 0., 200., 20);
	    FillHist("MET_VS_Eta_UseEvnet_tight", METauto, muon.Eta(), this_weight * weight_by_pt * SFs, 0., 200., 20., -3., 3., 100);
	    FillHist("Ptcone_UseEvent_tight", conept, this_weight * weight_by_pt * SFs, 0., 200., 200);
	    FillHist("Eta_UseEvent_tight", muon.Eta(), this_weight * weight_by_pt * SFs, -3., 3., 100);
	  }
	}// if istight = apply SFs for tight lepton
	//else FillDenAndNum(current_hist_name + "_Nvtx_reweight", muon, this_weight * weight_by_pt * Nvtx_reweight, IsThisTight);
	else FillDenAndNum(current_hist_name + "_Nvtx_reweight", muon, this_weight * weight_by_pt * MuTrkEffSF, IsThisTight);
	
	
        histfilled = true;
      }
    }//for loop all jets ends
  }//awayjet systematic for loop ends
  
  
}
//////////////////////////////////RunFakes ends


///functionto get JES, JER, lepton energy, lepton ID SF syst for prompt contribution
void Muon_FR_cal_all::RunFakes_subtraction(TString tight_ID, TString loose_ID, std::vector<snu::KMuon> loose_mu, bool debugging){
  std::map< TString, std::vector<double> > HLT_ptrange; //defining muon pt range for each single muon trigger
  std::map< TString, double > HLT_ptmin;
  HLT_ptrange["HLT_Mu20_v"].push_back(35.);
  HLT_ptrange["HLT_Mu20_v"].push_back(45.);
  HLT_ptmin["HLT_Mu20_v"] = 23.;

  HLT_ptrange["HLT_Mu27_v"].push_back(45.);
  HLT_ptrange["HLT_Mu27_v"].push_back(80.);
  HLT_ptmin["HLT_Mu27_v"] = 30.;

  HLT_ptrange["HLT_Mu50_v"].push_back(80.);
  HLT_ptrange["HLT_Mu50_v"].push_back(9999.);
  HLT_ptmin["HLT_Mu50_v"] = 55.;
  
  std::vector<TString> AllHLTs;
  for(std::map< TString, std::vector<double> >::iterator it=HLT_ptrange.begin(); it!=HLT_ptrange.end(); it++){
    AllHLTs.push_back(it->first);
  }

  bool DijetFake = std::find(k_flags.begin(), k_flags.end(), "DijetFake") != k_flags.end(); //flag for QCD MC fr cal
  bool DijetPrompt = std::find(k_flags.begin(), k_flags.end(), "DijetPrompt") != k_flags.end(); //flag for MC subtraction
  
  //call loose particles
  TString MuonLooseID_loosest = "MUON_SUSY_VETO";
  TString MuonVetoID_loosest = "MUON_SUSY_VETO";
  TString ElectronLooseID_loosest = "ELECTRON_HN_FAKELOOSEv7_loosest";
  TString ElectronVetoID_loosest = "ELECTRON_HN_VETO_loosest";
  //==== Muons
  std::vector<snu::KMuon> muons_loosest = GetMuons(MuonLooseID_loosest, true);
  std::vector<snu::KMuon> muons_veto_loosest = GetMuons(MuonVetoID_loosest, true);
  //==== Electrons
  std::vector<snu::KElectron> electrons_loosest = GetElectrons(false, true, ElectronLooseID_loosest);
  std::vector<snu::KElectron> electrons_veto_loosest = GetElectrons(true, true, ElectronVetoID_loosest);
  std::vector<snu::KElectron> electrons_veto;
  //==== Jets
  std::vector<snu::KJet> jets_eta5_nolepveto_loosest = GetJets("JET_HN_eta5_nolepveto_loosest", 10., 5.);
  
  TString syst_str[] = {"_MuonEn_up", "_MuonEn_down", "_JetEn_up", "_JetEn_down", "_JetRes_up", "_JetRes_down", "_Unclustered_up", "_Unclustered_down",
                        "_MuonIDSF_up", "_MuonIDSF_down"};
  int n_syst = 10;
  snu::KEvent Evt = eventbase->GetEvent();
  
  for(int i = 0; i < n_syst; i++){//for loop systemati categories
    TString this_syst = syst_str[i];
    
    if(i == 2) METauto = Evt.PFMETShifted(snu::KEvent::JetEn, snu::KEvent::up);//_JetEn_up
    if(i == 3) METauto = Evt.PFMETShifted(snu::KEvent::JetEn, snu::KEvent::down);//_JetEn_down
    if(i == 4) METauto = Evt.PFMETShifted(snu::KEvent::JetRes, snu::KEvent::up);//_JetRes_up
    if(i == 5) METauto = Evt.PFMETShifted(snu::KEvent::JetRes, snu::KEvent::down);//_JetRes_down
    if(i == 6) METauto = Evt.PFMETShifted(snu::KEvent::Unclustered, snu::KEvent::up);//_Unclustered_up
    if(i == 7) METauto = Evt.PFMETShifted(snu::KEvent::Unclustered, snu::KEvent::down);//_Unclustered_down
        
    //================
    //==== Make Muon
    //================
    double this_MuonLooseRelIso = 0.6;
    double this_MuonVetoRelIso = 0.6;
    std::vector<snu::KMuon> muons, muons_veto;
    muons.clear();
    muons_veto.clear();
    if(this_syst == "_MuonEn_up"){
      //==== Signal Leptons
      for(unsigned int j=0; j<muons_loosest.size(); j++){
	snu::KMuon this_muon = muons_loosest.at(j);
	this_muon.SetPtEtaPhiM( this_muon.Pt()*this_muon.PtShiftedUp(), this_muon.Eta(), this_muon.Phi(), this_muon.M() );
	double new_RelIso = this_muon.PFRelMiniIsoRho()/this_muon.PtShiftedUp();
	this_muon.SetRelMiniIsoRho(new_RelIso);
	if( this_muon.Pt() >= 10. && new_RelIso < this_MuonLooseRelIso ) muons.push_back( this_muon );
      }
      //==== Veto Leptons
      for(unsigned int j=0; j<muons_veto_loosest.size(); j++){
	snu::KMuon this_muon = muons_veto_loosest.at(j);
	this_muon.SetPtEtaPhiM( this_muon.Pt()*this_muon.PtShiftedUp(), this_muon.Eta(), this_muon.Phi(), this_muon.M() );
	double new_RelIso = this_muon.PFRelMiniIsoRho()/this_muon.PtShiftedUp();
	this_muon.SetRelMiniIsoRho(new_RelIso);
	if( this_muon.Pt() >= 5. && new_RelIso < this_MuonVetoRelIso ) muons_veto.push_back( this_muon );
      }
    }
    else if(this_syst == "_MuonEn_down"){
      //==== Signal Leptons
      for(unsigned int j=0; j<muons_loosest.size(); j++){
	snu::KMuon this_muon = muons_loosest.at(j);
	this_muon.SetPtEtaPhiM( this_muon.Pt()*this_muon.PtShiftedDown(), this_muon.Eta(), this_muon.Phi(), this_muon.M() );
	double new_RelIso = this_muon.PFRelMiniIsoRho()/this_muon.PtShiftedDown();
	this_muon.SetRelMiniIsoRho(new_RelIso);
	if( this_muon.Pt() >= 10. && new_RelIso < this_MuonLooseRelIso ) muons.push_back( this_muon );
      }
      //==== Veto Leptons
      for(unsigned int j=0; j<muons_veto_loosest.size(); j++){
	snu::KMuon this_muon = muons_veto_loosest.at(j);
	this_muon.SetPtEtaPhiM( this_muon.Pt()*this_muon.PtShiftedDown(), this_muon.Eta(), this_muon.Phi(), this_muon.M() );
	double new_RelIso = this_muon.PFRelMiniIsoRho()/this_muon.PtShiftedDown();
	this_muon.SetRelMiniIsoRho(new_RelIso);
	if( this_muon.Pt() >= 5. && new_RelIso < this_MuonVetoRelIso ) muons_veto.push_back( this_muon );
      }
    }
    //==== normal muons
    else{
      //==== Signal Leptons
      for(unsigned int j=0; j<muons_loosest.size(); j++){
	snu::KMuon this_muon = muons_loosest.at(j);
	if( this_muon.Pt() >= 10. && this_muon.PFRelMiniIsoRho() < this_MuonLooseRelIso ) muons.push_back( this_muon );
      }
      //==== Veto Leptons
      for(unsigned int j=0; j<muons_veto_loosest.size(); j++){
	snu::KMuon this_muon = muons_veto_loosest.at(j);
	if( this_muon.Pt() >= 5. && this_muon.PFRelMiniIsoRho() < this_MuonVetoRelIso ) muons_veto.push_back( this_muon );
      }
    }
    //std::sort(muons.begin(), muons.end(), MuonPtComparing);
    
    JSCorrectedMETRochester(muons, METauto, METphiauto);
    
    //===============
    //==== Make Jet
    //===============
    std::vector<snu::KJet> jets_eta5_nolepveto; // eta < 5, NO lepton-veto
    if( (this_syst == "_JetEn_up") || (this_syst == "_JetRes_up") ){
      for(unsigned int j=0; j<jets_eta5_nolepveto_loosest.size(); j++){
	snu::KJet this_jet = jets_eta5_nolepveto_loosest.at(j);

	double this_scaling = 1.;
	if(this_syst == "_JetEn_up") this_scaling = this_jet.ScaledUpEnergy();
	if(this_syst == "_JetRes_up") this_scaling = this_jet.SmearedResUp()/this_jet.SmearedRes();

	this_jet *= this_scaling;

	if(this_jet.Pt() >= 20.) jets_eta5_nolepveto.push_back(this_jet);
      }
    }
    else if( (this_syst == "_JetEn_down") || (this_syst == "_JetRes_down") ){
      for(unsigned int j=0; j<jets_eta5_nolepveto_loosest.size(); j++){
	snu::KJet this_jet = jets_eta5_nolepveto_loosest.at(j);

	double this_scaling = 1.;
	if(this_syst == "_JetEn_down") this_scaling = this_jet.ScaledDownEnergy();
	if(this_syst == "_JetRes_down") this_scaling = this_jet.SmearedResDown()/this_jet.SmearedRes();;

	this_jet *= this_scaling;

	if(this_jet.Pt() >= 20.) jets_eta5_nolepveto.push_back(this_jet);
      }
    }
    else{
      for(unsigned int j=0; j<jets_eta5_nolepveto_loosest.size(); j++){
	snu::KJet this_jet = jets_eta5_nolepveto_loosest.at(j);
	if(this_jet.Pt() >= 20.) jets_eta5_nolepveto.push_back(this_jet);
      }
    }

    std::vector<snu::KJet> jets; // eta < 2.5, lepton-veto, normla jet selection to be used as away jets
    for(unsigned int j=0; j<jets_eta5_nolepveto.size(); j++){
      snu::KJet this_jet = jets_eta5_nolepveto.at(j);
      bool IsNormalJet = fabs( this_jet.Eta() ) < 2.5;
      bool lepinside = HasLeptonInsideJet(this_jet, muons_veto, electrons_veto);
            
      if(IsNormalJet && !lepinside) jets.push_back( this_jet );
    }
    
    
    //now we have jets(normal jet selection with systematic up & down) and muons(loose ID with up & down E) & muons_veto, also METauto & METphiauto with correct shift
    int For_HLT_Mu3_PFJet40_v = 0;
    int n_bjets = 0;
    for(unsigned int i = 0; i < jets.size(); i++){
      if( IsBTagged(jets.at(i), snu::KJet::CSVv2, snu::KJet::Medium) ) n_bjets++;
      if( jets.at(i).Pt() > 50. ) For_HLT_Mu3_PFJet40_v++;
    }
    if(debugging) cout << "muons_veto.size() != 1 : " << muons_veto.size() <<endl;
    if(muons_veto.size() != 1) continue; //return;
    if(debugging) cout << "jets.size() : " << jets.size() << endl;
    if(jets.size() < 1) continue; //return; // Njet >= 1
    
    std::vector<snu::KMuon> hnloose;//keep only prompt muon for MC subtraction, keep only fake muon for QCD fr cal
    hnloose.clear();
    if(debugging) cout << "muons.size() : " << muons.size() << endl;
    if(muons.size() != 1) continue; //return;
    for(unsigned int i=0; i<muons.size(); i++){
      if(k_isdata) hnloose.push_back( muons.at(i) );
      if(!k_isdata){
	if(DijetFake){
	  if( !TruthMatched( muons.at(i)) ){
	    hnloose.push_back( muons.at(i) );
	  }
	}
	else if(DijetPrompt){
	  if(debugging) cout << "DijetPrompt" << endl;
	  if( TruthMatched( muons.at(i)) ){
	    //if(muons.at(i).GetType() == 1 || muons.at(i).GetType() == 8){
	    if(debugging) cout << "truth" << endl;
	    hnloose.push_back( muons.at(i) );
	  }
	}
	else{
	  continue; //return;
	}
      }
    }
    
    if(debugging) cout << "hnloose.size() : " << hnloose.size() << endl;
    if(hnloose.size() != 1) continue; //return;
    
    snu::KMuon muon = hnloose.at(0);
    TLorentzVector metvec;
    metvec.SetPtEtaPhiE( METauto, 0, METphiauto, METauto);
    double MTval = AnalyzerCore::MT( muon, metvec );
    if(debugging) cout << "METauto : " << METauto << ", MTval : " << MTval << endl;
    
    /////////////////////////////////////weights
    //---------------------------------------------
    //Weight by pt range for each HLT
    //---------------------------------------------
    TString current_trigger;
    double weight_by_pt(0.);//calculate weight for this event
    for(std::map< TString, std::vector<double> >::iterator it=HLT_ptrange.begin(); it!=HLT_ptrange.end(); it++){
      double tmp = GetTriggerWeightByPtRange(it->first, it->second, HLT_ptmin[it->first], hnloose, For_HLT_Mu3_PFJet40_v);
      if(tmp > 0.1){
        weight_by_pt = tmp;
        current_trigger = it->first;
        break;
      }
    }
    //---------------------------------------------
    //pileup_reweight (69 mb minbias X-section)
    //---------------------------------------------
    float pileup_reweight=(1.0);
    if(!k_isdata){
      pileup_reweight = mcdata_correction->CatPileupWeight(eventbase->GetEvent(),0);
    }
    //---------------------------------------------
    //Additional Nvtx correction
    //---------------------------------------------
    double Nvtx_reweight = 1.;
    int nvtx = Evt.nVertices();
    if(!k_isdata && weight_by_pt > 0.1){
      Nvtx_reweight = mcdata_correction -> GetVtxReweight(current_trigger, nvtx);
    }
    //---------------------------------------------
    //Constant Normalization Factor
    //---------------------------------------------
    double constant_SF = 1.;
    //if(k_isdata) constant_SF = 1.;
    //---------------------------------------------
    //Muon ID scale-factor, ID SF is applied only for tight muon not for loose
    //---------------------------------------------
    
    bool IsThisTight = PassID( muon, tight_ID) && (muon.PFRelMiniIsoRho() < 0.2);
    
    double muon_SF = 1.;
    double MuTrkEffSF =  mcdata_correction->MuonTrackingEffScaleFactor(hnloose);
    if(IsThisTight && !k_isdata){
      if(this_syst == "_MuonIDSF_up") muon_SF = mcdata_correction->MuonScaleFactor("MUON_HN_TIGHT", hnloose, 1);
      else if(this_syst == "_MuonIDSF_down") muon_SF = mcdata_correction->MuonScaleFactor("MUON_HN_TIGHT", hnloose, -1);
      else muon_SF = mcdata_correction->MuonScaleFactor("MUON_HN_TIGHT", hnloose, 0);
      //muon_SF *= MuTrkEffSF;
    }
    double SFs = MuTrkEffSF;
    if(IsThisTight) SFs *= muon_SF;
    double this_weight = weight*pileup_reweight*SFs;
    
    if(DijetPrompt && (METauto > 80 || MTval > 25) ) continue; //return;//give event cut for MC subtraction case
    double conept = muon.Pt()*(1 + max(0.,(muon.PFRelMiniIsoRho()-0.2) ) );
    
    double AwayjetPt = 40.;
    bool histfilled = false; 
    for(unsigned int j = 0; j < jets.size(); j++){//for loop all jets in the event
      if(histfilled) break;
      snu::KJet jet = jets.at(j);
      if( jet.Pt() < AwayjetPt ) continue;
      
      double dPhi = muon.DeltaPhi( jet );
      
      bool UseEvent = false;
      if(debugging) cout << "dPhi : " << dPhi << ". (jet.Pt()/muon.Pt() : " << (jet.Pt()/muon.Pt()) << endl;
      UseEvent = (dPhi > 2.5) && (jet.Pt()/muon.Pt() > 1.);

      if(debugging) cout <<  "UseEvent : " << UseEvent << endl;
      
      if( UseEvent ){
	TString current_hist_name = "SingleMuonTrigger_Dijet"+ this_syst +"_all_" + loose_ID;
	FillDenAndNum(current_hist_name, muon, this_weight * weight_by_pt, IsThisTight);
	histfilled = true;
      }
    }//for loop all jet ends
  }//for loop systematic categories ends 
}//RunFakes_subtraction ends


double Muon_FR_cal_all::GetTriggerWeightByPtRange(TString hltname, vector<double> ptrange, double trigger_safe_pt, std::vector<snu::KMuon> muons, int npfjet50){

  double prescale_trigger = 0.;
  
  if(muons.size()==1){
    snu::KMuon muon = muons.at(0);
    double min_pt_cone = ptrange.at(0);
    double max_pt_cone = ptrange.at(1);
    bool SafePt = (muon.MiniAODPt() > trigger_safe_pt);

    if(SafePt && PassTrigger(hltname)){
      double TightISO = 0.2;
      double conept = muon.Pt()*(1 + max(0.,(muon.PFRelMiniIsoRho()-TightISO) ) );

      if(conept >= min_pt_cone && conept < max_pt_cone){
        if(hltname=="HLT_Mu20_v"){
	  //prescale_trigger = WeightByTrigger(hltname, TargetLumi) ;
	  prescale_trigger = 140.315;
	}
	else if(hltname=="HLT_Mu27_v"){
	  prescale_trigger = 250.252;
	}
	//else prescale_trigger = WeightByTrigger(hltname, TargetLumi) ;
	else prescale_trigger = 35863.308;
      }
      
    }
    
  }
  
  if(prescale_trigger == 0.) return 0.;
  if(k_isdata) return 1.;
  return prescale_trigger;

}

void Muon_FR_cal_all::FillDenAndNum(TString prefix, snu::KMuon muon, double thisweight, bool isTight){
  
  //float etaarray [] = {0.0, 0.8, 1.479, 2.0, 2.5};
  //float ptarray [] = {0., 5., 12., 15., 20., 25., 30., 35., 45., 60., 100., 200.};
  
  //new binning for new pt ranges 
  float etaarray [] = {0.0, 0.8, 1.479, 2.4};
  //float ptarray [] = {0., 5., 10., 15., 25., 35., 50., 70.};
  //float ptarray [] = {0., 10., 20., 30., 40., 50., 60., 70.};
  //float ptarray [] = {5., 10., 20., 30., 40., 50., 60., 70.};
  //float ptarray [] = {35, 45, 75, 80, 110, 150, 200, 250, 300, 350, 500, 1000, 1500, 2000};
    float ptarray [] = {35, 45, 55, 65, 80,  95,  110, 125, 140, 200, 250, 300,  350,  500, 1000, 1500, 2000};
  
  
  //float ptaaray [] = {5., 12., };
  
    
  float etaarray_2 [] = {0.0, 1.479, 2.5};
  float ptarray_2 [] = {10.,15.,40.,200.};
  
  double TightISO = 0.2;
  double conept = muon.Pt()*(1 + max(0.,(muon.PFRelMiniIsoRho()-TightISO) ) );

  /*
  FillHist(prefix+"_eta_F0", muon.Eta(), thisweight, -3., 3., 30);
  FillHist(prefix+"_pt_F0", muon.Pt(), thisweight, 0., 200., 200);
  FillHist(prefix+"_pt_cone_F0", conept, thisweight, 0., 200., 200);
  FillHist(prefix+"_RelIso_F0", muon.RelIso04(), thisweight, 0., 1., 100);
  FillHist(prefix+"_Chi2_F0", muon.GlobalChi2(), thisweight, 0., 50., 50);
  FillHist(prefix+"_dXY_F0", fabs(muon.dXY()), thisweight, 0., 1., 1000);
  FillHist(prefix+"_dXYSig_F0", fabs(muon.dXYSig()), thisweight, 0., 15., 150);
  FillHist(prefix+"_dZ_F0", fabs(muon.dZ()), thisweight, 0., 0.5, 50);
  FillHist(prefix+"_onebin_F0", 0.5, thisweight, 0., 1., 1);
  */
  //FillHist(prefix+"_events_pt_vs_eta_F0", muon.Pt(), fabs(muon.Eta()), thisweight, ptarray, 11, etaarray, 4);
  //FillHist(prefix+"_events_pt_cone_vs_eta_F0", conept, fabs(muon.Eta()), thisweight, ptarray, 11, etaarray, 4);
  FillHist(prefix+"_events_pt_vs_eta_F0", muon.Pt(), fabs(muon.Eta()), thisweight, ptarray, 16, etaarray, 3);
  FillHist(prefix+"_events_pt_cone_vs_eta_F0", conept, fabs(muon.Eta()), thisweight, ptarray, 16, etaarray, 3);
  FillHist(prefix+"_PFMET_F0", METauto, thisweight, 0., 1000., 1000);
  
  if( isTight ){
    /*
    FillHist(prefix+"_eta_F", muon.Eta(), thisweight, -3., 3., 30);
    FillHist(prefix+"_pt_F", muon.Pt(), thisweight, 0., 200., 200);
    FillHist(prefix+"_pt_cone_F", conept, thisweight, 0., 200., 200);
    FillHist(prefix+"_RelIso_F", muon.RelIso04(), thisweight, 0., 1., 100);
    FillHist(prefix+"_Chi2_F", muon.GlobalChi2(), thisweight, 0., 50., 50);
    FillHist(prefix+"_dXY_F", fabs(muon.dXY()), thisweight, 0., 1., 1000);
    FillHist(prefix+"_dXYSig_F", fabs(muon.dXYSig()), thisweight, 0., 15., 150);
    FillHist(prefix+"_dZ_F", fabs(muon.dZ()), thisweight, 0., 0.5, 50);
    FillHist(prefix+"_onebin_F", 0.5, thisweight, 0., 1., 1);
    */
    //FillHist(prefix+"_events_pt_vs_eta_F", muon.Pt(), fabs(muon.Eta()), thisweight, ptarray, 11, etaarray, 4);
    //FillHist(prefix+"_events_pt_cone_vs_eta_F", conept, fabs(muon.Eta()), thisweight, ptarray, 11, etaarray, 4);
    FillHist(prefix+"_events_pt_vs_eta_F", muon.Pt(), fabs(muon.Eta()), thisweight, ptarray, 16, etaarray, 3);
    FillHist(prefix+"_events_pt_cone_vs_eta_F", conept, fabs(muon.Eta()), thisweight, ptarray, 16, etaarray, 3);
    FillHist(prefix+"_PFMET_F", METauto, thisweight, 0., 1000., 1000);
  }
  

}


void Muon_FR_cal_all::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);
  m_logger<< INFO << "Number of events that pass 1 7GeV trigger = " << n_17_pass  << LQLogger::endmsg;
  m_logger<< INFO << "Number of events that pass 17 GeV + jet trigger = " << n_17_jet_pass  << LQLogger::endmsg;
  m_logger<< INFO << "Number of events that pass 17 GeV || jet trigger = " << n_17_17_jet_pass  << LQLogger::endmsg;

}

void Muon_FR_cal_all::BeginCycle() throw( LQError ){
  
  Message("In begin Cycle", INFO);
  
  //If you wish to output variables to output file use DeclareVariable
  // clear these variables in ::ClearOutputVectors function
  //DeclareVariable(obj, label, treename );
  //DeclareVariable(obj, label ); //-> will use default treename: LQTree
  //DeclareVariable(out_electrons, "Signal_Electrons", "LQTree");
  //DeclareVariable(out_muons, "Signal_Muons");

  n_17_jet_pass=0;
  n_17_17_jet_pass=0;
  n_17_pass=0;

  
  return;
  
}

Muon_FR_cal_all::~Muon_FR_cal_all() {
  
  Message("In Muon_FR_cal_all Destructor" , INFO);
  
}



void Muon_FR_cal_all::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


///############### THESE ARE FUNCTIONS SPECIFIC TO THIS CYCLE

void Muon_FR_cal_all::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this Muon_FR_cal_allCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void Muon_FR_cal_all::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}


