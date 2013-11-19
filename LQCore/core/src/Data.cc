#ifndef DATA_cc
#define DATA_cc

#include "Data.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TList.h>
#include <TFriendElement.h>
#include <TVirtualIndex.h>
#include <iostream>

Data::Data()
{

   //Init(tree);

}

Data::~Data()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t Data::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Long64_t Data::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Data::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).


    // Set object pointer
   HLTKey = 0;
   HLTInsideDatasetTriggerNames = 0;
   HLTOutsideDatasetTriggerNames = 0;
   HLTFilterName = 0;
   ElectronGsfCtfCharge = 0;
   ElectronGsfCtfScPixCharge = 0;
   ElectronGsfScPixCharge = 0;
   ElectronHLTDoubleEleMatched = 0;
   ElectronHLTSingleEleMatched = 0;
   ElectronHLTSingleEleWP80Matched = 0;
   ElectronHasEcalDrivenSeed = 0;
   ElectronHasMatchedConvPhot = 0;
   ElectronHasTrackerDrivenSeed = 0;
   ElectronIsEB = 0;
   ElectronIsEE = 0;
   MuonHLTSingleIsoMuonMatched = 0;
   MuonHLTSingleMuonMatched = 0;
   PhotonHasMatchedConvPhot = 0;
   PhotonHasMatchedPromptEle = 0;
   PhotonHasPixelSeed = 0;
   PhotonIsEBEEGap = 0;
   PhotonIsEBGap = 0;
   PhotonIsEEGap = 0;
   HLTInsideDatasetTriggerDecisions = 0;
   HLTOutsideDatasetTriggerDecisions = 0;
   VertexIsFake = 0;
   CaloMET = 0;
   CaloMETPhi = 0;
   CaloSumET = 0;
   CaloMETPhiType1Cor = 0;
   CaloMETPhiUncorrType1Cor = 0;
   CaloMETType1Cor = 0;
   CaloMETUncorrType1Cor = 0;
   CaloSumETType1Cor = 0;
   CaloSumETUncorrType1Cor = 0;
   ElectronBeamSpotDXY = 0;
   ElectronBeamSpotDXYError = 0;
   ElectronCaloEnergy = 0;
   ElectronConvFitProb = 0;
   ElectronDCotTheta = 0;
   ElectronDeltaEtaTrkSC = 0;
   ElectronDeltaPhiTrkSC = 0;
   ElectronDist = 0;
   ElectronE1x5OverE5x5 = 0;
   ElectronE2x5OverE5x5 = 0;
   ElectronESuperClusterOverP = 0;
   ElectronEcalIsoDR03 = 0;
   ElectronEcalIsoPAT = 0;
   ElectronEnergy = 0;
   ElectronEta = 0;
   ElectronFbrem = 0;
   ElectronHLTDoubleEleMatchEta = 0;
   ElectronHLTDoubleEleMatchPhi = 0;
   ElectronHLTDoubleEleMatchPt = 0;
   ElectronHLTSingleEleMatchEta = 0;
   ElectronHLTSingleEleMatchPhi = 0;
   ElectronHLTSingleEleMatchPt = 0;
   ElectronHLTSingleEleWP80MatchEta = 0;
   ElectronHLTSingleEleWP80MatchPhi = 0;
   ElectronHLTSingleEleWP80MatchPt = 0;
   ElectronHcalIsoD1DR03 = 0;
   ElectronHcalIsoD2DR03 = 0;
   ElectronHcalIsoDR03 = 0;
   ElectronHcalIsoDR03FullCone = 0;
   ElectronHcalIsoPAT = 0;
   ElectronHoE = 0;
   ElectronLeadVtxDistXY = 0;
   ElectronLeadVtxDistZ = 0;
   ElectronMatchedGenParticleEta = 0;
   ElectronMatchedGenParticlePhi = 0;
   ElectronMatchedGenParticlePt = 0;
   ElectronPFChargedHadronIso03 = 0;
   ElectronPFChargedHadronIso04 = 0;
   ElectronPFNeutralHadronIso03 = 0;
   ElectronPFNeutralHadronIso04 = 0;
   ElectronPFPhotonIso03 = 0;
   ElectronPFPhotonIso04 = 0;
   ElectronPhi = 0;
   ElectronPrimaryVertexDXY = 0;
   ElectronPrimaryVertexDXYError = 0;
   ElectronPt = 0;
   ElectronPtHeep = 0;
   ElectronRelIsoPAT = 0;
   ElectronSCEta = 0;
   ElectronSCPhi = 0;
   ElectronSCPt = 0;
   ElectronSCRawEnergy = 0;
   ElectronSigmaEtaEta = 0;
   ElectronSigmaIEtaIEta = 0;
   ElectronTrackPt = 0;
   ElectronTrackValidFractionOfHits = 0;
   ElectronTrackVx = 0;
   ElectronTrackVy = 0;
   ElectronTrackVz = 0;
   ElectronTrkIsoDR03 = 0;
   ElectronTrkIsoPAT = 0;
   ElectronVtxDistXY = 0;
   ElectronVtxDistZ = 0;
   GenWElectronEnergy = 0;
   GenWElectronEta = 0;
   GenWElectronP = 0;
   GenWElectronPhi = 0;
   GenWElectronPt = 0;
   GenWElectronPx = 0;
   GenWElectronPy = 0;
   GenWElectronPz = 0;
   GenWElectronTauVisibleEta = 0;
   GenWElectronTauVisiblePhi = 0;
   GenWElectronTauVisiblePt = 0;
   GenWElectronVX = 0;
   GenWElectronVY = 0;
   GenWElectronVZ = 0;
   GenZElectronEnergy = 0;
   GenZElectronEta = 0;
   GenZElectronP = 0;
   GenZElectronPhi = 0;
   GenZElectronPt = 0;
   GenZElectronPx = 0;
   GenZElectronPy = 0;
   GenZElectronPz = 0;
   GenZElectronTauVisibleEta = 0;
   GenZElectronTauVisiblePhi = 0;
   GenZElectronTauVisiblePt = 0;
   GenZElectronVX = 0;
   GenZElectronVY = 0;
   GenZElectronVZ = 0;
   PDFCTEQWeights = 0;
   PDFMSTWWeights = 0;
   PDFNNPDFWeights = 0;
   GenJetEMF = 0;
   GenJetEnergy = 0;
   GenJetEta = 0;
   GenJetHADF = 0;
   GenJetP = 0;
   GenJetPhi = 0;
   GenJetPt = 0;
   GenMETCalo = 0;
   GenMETPhiCalo = 0;
   GenSumETCalo = 0;
   GenMETPhiTrue = 0;
   GenMETTrue = 0;
   GenSumETTrue = 0;
   GenWMuEnergy = 0;
   GenWMuEta = 0;
   GenWMuP = 0;
   GenWMuPhi = 0;
   GenWMuPt = 0;
   GenWMuPx = 0;
   GenWMuPy = 0;
   GenWMuPz = 0;
   GenWMuTauVisibleEta = 0;
   GenWMuTauVisiblePhi = 0;
   GenWMuTauVisiblePt = 0;
   GenWMuVX = 0;
   GenWMuVY = 0;
   GenWMuVZ = 0;
   GenZMuEnergy = 0;
   GenZMuEta = 0;
   GenZMuP = 0;
   GenZMuPhi = 0;
   GenZMuPt = 0;
   GenZMuPx = 0;
   GenZMuPy = 0;
   GenZMuPz = 0;
   GenZMuTauVisibleEta = 0;
   GenZMuTauVisiblePhi = 0;
   GenZMuTauVisiblePt = 0;
   GenZMuVX = 0;
   GenZMuVY = 0;
   GenZMuVZ = 0;
   GenParticleEnergy = 0;
   GenParticleEta = 0;
   GenParticleP = 0;
   GenParticlePhi = 0;
   GenParticlePt = 0;
   GenParticlePx = 0;
   GenParticlePy = 0;
   GenParticlePz = 0;
   GenParticleTauVisibleEta = 0;
   GenParticleTauVisiblePhi = 0;
   GenParticleTauVisiblePt = 0;
   GenParticleVX = 0;
   GenParticleVY = 0;
   GenParticleVZ = 0;
   GenWTauEnergy = 0;
   GenWTauEta = 0;
   GenWTauP = 0;
   GenWTauPhi = 0;
   GenWTauPt = 0;
   GenWTauPx = 0;
   GenWTauPy = 0;
   GenWTauPz = 0;
   GenWTauTauVisibleEta = 0;
   GenWTauTauVisiblePhi = 0;
   GenWTauTauVisiblePt = 0;
   GenWTauVX = 0;
   GenWTauVY = 0;
   GenWTauVZ = 0;
   GenZTauEnergy = 0;
   GenZTauEta = 0;
   GenZTauP = 0;
   GenZTauPhi = 0;
   GenZTauPt = 0;
   GenZTauPx = 0;
   GenZTauPy = 0;
   GenZTauPz = 0;
   GenZTauTauVisibleEta = 0;
   GenZTauTauVisiblePhi = 0;
   GenZTauTauVisiblePt = 0;
   GenZTauVX = 0;
   GenZTauVY = 0;
   GenZTauVZ = 0;
   HPSTauAgainstElectronDeadECALDiscr = 0;
   HPSTauAgainstElectronLooseDiscr = 0;
   HPSTauAgainstElectronLooseMVA2Discr = 0;
   HPSTauAgainstElectronLooseMVA3Discr = 0;
   HPSTauAgainstElectronMVA2categoryDiscr = 0;
   HPSTauAgainstElectronMVA2rawDiscr = 0;
   HPSTauAgainstElectronMVA3categoryDiscr = 0;
   HPSTauAgainstElectronMVA3rawDiscr = 0;
   HPSTauAgainstElectronMVADiscr = 0;
   HPSTauAgainstElectronMediumDiscr = 0;
   HPSTauAgainstElectronMediumMVA2Discr = 0;
   HPSTauAgainstElectronMediumMVA3Discr = 0;
   HPSTauAgainstElectronTightDiscr = 0;
   HPSTauAgainstElectronTightMVA2Discr = 0;
   HPSTauAgainstElectronTightMVA3Discr = 0;
   HPSTauAgainstElectronVLooseMVA2Discr = 0;
   HPSTauAgainstElectronVTightMVA3Discr = 0;
   HPSTauAgainstMuonLoose2Discr = 0;
   HPSTauAgainstMuonLooseDiscr = 0;
   HPSTauAgainstMuonMedium2Discr = 0;
   HPSTauAgainstMuonMediumDiscr = 0;
   HPSTauAgainstMuonTight2Discr = 0;
   HPSTauAgainstMuonTightDiscr = 0;
   HPSTauBremsRecoveryEOverPLead = 0;
   HPSTauCombinedIsolationDeltaBetaCorr3HitsDiscr = 0;
   HPSTauDecayModeFindingDiscr = 0;
   HPSTauEcalStripSumEOverPLead = 0;
   HPSTauEmFraction = 0;
   HPSTauEt = 0;
   HPSTauEta = 0;
   HPSTauEtaLeadCharged = 0;
   HPSTauEtaetaMoment = 0;
   HPSTauEtaphiMoment = 0;
   HPSTauHcal3x3OverPLead = 0;
   HPSTauHcalMaxOverPLead = 0;
   HPSTauHcalTotOverPLead = 0;
   HPSTauIsolationMVArawDiscr = 0;
   HPSTauIsolationPFChargedHadrCandsPtSum = 0;
   HPSTauIsolationPFGammaCandsEtSum = 0;
   HPSTauLeadPFChargedHadrCandsignedSipt = 0;
   HPSTauLeadVtxDistXY = 0;
   HPSTauLeadVtxDistZ = 0;
   HPSTauLooseCombinedIsolationDeltaBetaCorr3HitsDiscr = 0;
   HPSTauLooseCombinedIsolationDeltaBetaCorrDiscr = 0;
   HPSTauLooseIsolationDeltaBetaCorrDiscr = 0;
   HPSTauLooseIsolationDiscr = 0;
   HPSTauLooseIsolationMVA2Discr = 0;
   HPSTauLooseIsolationMVADiscr = 0;
   HPSTauMatchedGenJetEta = 0;
   HPSTauMatchedGenJetPhi = 0;
   HPSTauMatchedGenJetPt = 0;
   HPSTauMatchedGenParticleEta = 0;
   HPSTauMatchedGenParticlePhi = 0;
   HPSTauMatchedGenParticlePt = 0;
   HPSTauMaximumHCALPFClusterEt = 0;
   HPSTauMediumCombinedIsolationDeltaBetaCorr3HitsDiscr = 0;
   HPSTauMediumCombinedIsolationDeltaBetaCorrDiscr = 0;
   HPSTauMediumIsolationDeltaBetaCorrDiscr = 0;
   HPSTauMediumIsolationDiscr = 0;
   HPSTauMediumIsolationMVA2Discr = 0;
   HPSTauMediumIsolationMVADiscr = 0;
   HPSTauPhi = 0;
   HPSTauPhiLeadCharged = 0;
   HPSTauPhiphiMoment = 0;
   HPSTauPt = 0;
   HPSTauPtLeadCharged = 0;
   HPSTauSignalPFChargedHadrCandsCount = 0;
   HPSTauSignalPFChargedHadrCandsEta = 0;
   HPSTauSignalPFChargedHadrCandsPhi = 0;
   HPSTauSignalPFChargedHadrCandsPt = 0;
   HPSTauSignalPFGammaCandsCount = 0;
   HPSTauSignalPFGammaCandsEta = 0;
   HPSTauSignalPFGammaCandsPhi = 0;
   HPSTauSignalPFGammaCandsPt = 0;
   HPSTauSignalPFNeutrHadrCandsCount = 0;
   HPSTauSignalPFNeutrHadrCandsEta = 0;
   HPSTauSignalPFNeutrHadrCandsPhi = 0;
   HPSTauSignalPFNeutrHadrCandsPt = 0;
   HPSTauTightCombinedIsolationDeltaBetaCorr3HitsDiscr = 0;
   HPSTauTightCombinedIsolationDeltaBetaCorrDiscr = 0;
   HPSTauTightIsolationDeltaBetaCorrDiscr = 0;
   HPSTauTightIsolationDiscr = 0;
   HPSTauTightIsolationMVA2Discr = 0;
   HPSTauTightIsolationMVADiscr = 0;
   HPSTauVLooseCombinedIsolationDeltaBetaCorrDiscr = 0;
   HPSTauVLooseIsolationDeltaBetaCorrDiscr = 0;
   HPSTauVLooseIsolationDiscr = 0;
   HPSTauVtxDistXY = 0;
   HPSTauVtxDistZ = 0;
   MuonBackToBackCompatibility = 0;
   MuonBeamSpotDXY = 0;
   MuonBeamSpotDXYError = 0;
   MuonBestTrackVtxDistXY = 0;
   MuonBestTrackVtxDistZ = 0;
   MuonCocktailEta = 0;
   MuonCocktailEtaError = 0;
   MuonCocktailGlobalChi2 = 0;
   MuonCocktailP = 0;
   MuonCocktailPhi = 0;
   MuonCocktailPhiError = 0;
   MuonCocktailPt = 0;
   MuonCocktailPtError = 0;
   MuonCocktailQOverPError = 0;
   MuonCocktailTrkD0 = 0;
   MuonCocktailTrkD0Error = 0;
   MuonCocktailTrkDz = 0;
   MuonCocktailTrkDzError = 0;
   MuonCocktailTrkValidFractionOfHits = 0;
   MuonCosmicCompatibility = 0;
   MuonEcalIso = 0;
   MuonEcalVetoIso = 0;
   MuonEnergy = 0;
   MuonEta = 0;
   MuonEtaError = 0;
   MuonGlobalChi2 = 0;
   MuonHLTSingleIsoMuonMatchEta = 0;
   MuonHLTSingleIsoMuonMatchPhi = 0;
   MuonHLTSingleIsoMuonMatchPt = 0;
   MuonHLTSingleMuonMatchEta = 0;
   MuonHLTSingleMuonMatchPhi = 0;
   MuonHLTSingleMuonMatchPt = 0;
   MuonHOIso = 0;
   MuonHcalIso = 0;
   MuonHcalVetoIso = 0;
   MuonMatchedGenParticleEta = 0;
   MuonMatchedGenParticlePhi = 0;
   MuonMatchedGenParticlePt = 0;
   MuonOverlapCompatibility = 0;
   MuonP = 0;
   MuonPFIsoR03ChargedHadron = 0;
   MuonPFIsoR03ChargedParticle = 0;
   MuonPFIsoR03NeutralHadron = 0;
   MuonPFIsoR03NeutralHadronHT = 0;
   MuonPFIsoR03PU = 0;
   MuonPFIsoR03Photon = 0;
   MuonPFIsoR03PhotonHT = 0;
   MuonPFIsoR04ChargedHadron = 0;
   MuonPFIsoR04ChargedParticle = 0;
   MuonPFIsoR04NeutralHadron = 0;
   MuonPFIsoR04NeutralHadronHT = 0;
   MuonPFIsoR04PU = 0;
   MuonPFIsoR04Photon = 0;
   MuonPFIsoR04PhotonHT = 0;
   MuonPhi = 0;
   MuonPhiError = 0;
   MuonPrimaryVertexDXY = 0;
   MuonPrimaryVertexDXYError = 0;
   MuonPt = 0;
   MuonPtError = 0;
   MuonQOverPError = 0;
   MuonTimeCompatibility = 0;
   MuonTrackChi2 = 0;
   MuonTrackerIsoSumPT = 0;
   MuonTrkD0 = 0;
   MuonTrkD0Error = 0;
   MuonTrkDz = 0;
   MuonTrkDzError = 0;
   MuonTrkEta = 0;
   MuonTrkEtaError = 0;
   MuonTrkIso = 0;
   MuonTrkPhi = 0;
   MuonTrkPhiError = 0;
   MuonTrkPt = 0;
   MuonTrkPtError = 0;
   MuonTrkValidFractionOfHits = 0;
   MuonTrkVx = 0;
   MuonTrkVy = 0;
   MuonTrkVz = 0;
   MuonVtxDistXY = 0;
   MuonVtxDistZ = 0;
   PFCandEnergyLeptLink = 0;
   PFCandEtaLeptLink = 0;
   PFCandPhiLeptLink = 0;
   PFCandPtLeptLink = 0;
   PFJetBestVertexTrackAssociationFactor = 0;
   PFJetBeta = 0;
   PFJetBetaClassic = 0;
   PFJetBetaStar = 0;
   PFJetBetaStarClassic = 0;
   PFJetChargedEmEnergyFraction = 0;
   PFJetChargedHadronEnergyFraction = 0;
   PFJetChargedMuEnergyFraction = 0;
   PFJetClosestVertexWeighted3DSeparation = 0;
   PFJetClosestVertexWeightedXYSeparation = 0;
   PFJetClosestVertexWeightedZSeparation = 0;
   PFJetCombinedInclusiveSecondaryVertexBTag = 0;
   PFJetCombinedMVABTag = 0;
   PFJetCombinedSecondaryVertexBTag = 0;
   PFJetCombinedSecondaryVertexMVABTag = 0;
   PFJetElectronEnergyFraction = 0;
   PFJetEnergy = 0;
   PFJetEnergyRaw = 0;
   PFJetEta = 0;
   PFJetHFEMEnergyFraction = 0;
   PFJetHFHadronEnergyFraction = 0;
   PFJetJECUnc = 0;
   PFJetJetBProbabilityBTag = 0;
   PFJetJetProbabilityBTag = 0;
   PFJetL1FastJetJEC = 0;
   PFJetL2L3ResJEC = 0;
   PFJetL2RelJEC = 0;
   PFJetL3AbsJEC = 0;
   PFJetMuonEnergyFraction = 0;
   PFJetNeutralEmEnergyFraction = 0;
   PFJetNeutralHadronEnergyFraction = 0;
   PFJetPhi = 0;
   PFJetPhotonEnergyFraction = 0;
   PFJetPt = 0;
   PFJetPtRaw = 0;
   PFJetSimpleSecondaryVertexHighEffBTag = 0;
   PFJetSimpleSecondaryVertexHighPurBTag = 0;
   PFJetSoftElectronByIP3dBTag = 0;
   PFJetSoftElectronByPtBTag = 0;
   PFJetSoftMuonBTag = 0;
   PFJetSoftMuonByIP3dBTag = 0;
   PFJetSoftMuonByPtBTag = 0;
   PFJetTrackCountingHighEffBTag = 0;
   PFJetTrackCountingHighPurBTag = 0;
   PFMET = 0;
   PFMETPhi = 0;
   PFSumET = 0;
   PFMETPhiType01Cor = 0;
   PFMETType01Cor = 0;
   PFSumETType01Cor = 0;
   PFMETPhiType01XYCor = 0;
   PFMETType01XYCor = 0;
   PFSumETType01XYCor = 0;
   PFMETPhiType1Cor = 0;
   PFMETType1Cor = 0;
   PFSumETType1Cor = 0;
   PhotonAlpha = 0;
   PhotonChi2ConvPhot = 0;
   PhotonDPhiTracksAtVtxConvPhot = 0;
   PhotonDistOfMinApproachConvPhot = 0;
   PhotonE2OverE9 = 0;
   PhotonE3x3 = 0;
   PhotonE4SwissCross = 0;
   PhotonE5x5 = 0;
   PhotonEOverPConvPhot = 0;
   PhotonEcalIsoDR03 = 0;
   PhotonEcalIsoDR04 = 0;
   PhotonEnergy = 0;
   PhotonEta = 0;
   PhotonHcalIsoDR03 = 0;
   PhotonHcalIsoDR03FullCone = 0;
   PhotonHcalIsoDR04 = 0;
   PhotonHcalIsoDR04FullCone = 0;
   PhotonHoE = 0;
   PhotonNDofConvPhot = 0;
   PhotonPairCotThetaSeparationConvPhot = 0;
   PhotonPairInvariantMassConvPhot = 0;
   PhotonPairMomentumxConvPhot = 0;
   PhotonPairMomentumyConvPhot = 0;
   PhotonPairMomentumzConvPhot = 0;
   PhotonPhi = 0;
   PhotonPt = 0;
   PhotonSCenergy = 0;
   PhotonSCeta = 0;
   PhotonSCphi = 0;
   PhotonSCseedEnergy = 0;
   PhotonSEtaEta = 0;
   PhotonSEtaPhi = 0;
   PhotonSMajMaj = 0;
   PhotonSMinMin = 0;
   PhotonSPhiPhi = 0;
   PhotonSigmaIEtaIEta = 0;
   PhotonTimeSeed = 0;
   PhotonTrkIsoHollowDR03 = 0;
   PhotonTrkIsoHollowDR04 = 0;
   PhotonTrkIsoSolidDR03 = 0;
   PhotonTrkIsoSolidDR04 = 0;
   PhotonXVtxConvPhot = 0;
   PhotonYVtxConvPhot = 0;
   PhotonZVtxConvPhot = 0;
   TCMET = 0;
   TCMETPhi = 0;
   TCSumET = 0;
   VertexChi2 = 0;
   VertexNDF = 0;
   VertexRho = 0;
   VertexX = 0;
   VertexXErr = 0;
   VertexY = 0;
   VertexYErr = 0;
   VertexZ = 0;
   VertexZErr = 0;
   PileUpInteractionsTrue = 0;
   HLTFilterObjEta = 0;
   HLTFilterObjPhi = 0;
   HLTFilterObjPt = 0;
   ElectronCharge = 0;
   ElectronClassif = 0;
   ElectronMissingHits = 0;
   ElectronMissingHitsEG = 0;
   ElectronNumberOfBrems = 0;
   ElectronOverlaps = 0;
   ElectronPassEGammaIDEoP = 0;
   ElectronPassEGammaIDLoose = 0;
   ElectronPassEGammaIDMedium = 0;
   ElectronPassEGammaIDTight = 0;
   ElectronPassEGammaIDTrigTight = 0;
   ElectronPassEGammaIDTrigWP70 = 0;
   ElectronPassEGammaIDVeto = 0;
   ElectronPassId = 0;
   ElectronPassIsoPAT = 0;
   ElectronVtxIndex = 0;
   GenWElectronMotherIndex = 0;
   GenWElectronNumDaught = 0;
   GenWElectronPdgId = 0;
   GenWElectronStatus = 0;
   GenWElectronTauDecayMode = 0;
   GenZElectronMotherIndex = 0;
   GenZElectronNumDaught = 0;
   GenZElectronPdgId = 0;
   GenZElectronStatus = 0;
   GenZElectronTauDecayMode = 0;
   PileUpInteractions = 0;
   PileUpOriginBX = 0;
   GenWMuMotherIndex = 0;
   GenWMuNumDaught = 0;
   GenWMuPdgId = 0;
   GenWMuStatus = 0;
   GenWMuTauDecayMode = 0;
   GenZMuMotherIndex = 0;
   GenZMuNumDaught = 0;
   GenZMuPdgId = 0;
   GenZMuStatus = 0;
   GenZMuTauDecayMode = 0;
   GenParticleMotherIndex = 0;
   GenParticleNumDaught = 0;
   GenParticlePdgId = 0;
   GenParticleStatus = 0;
   GenParticleTauDecayMode = 0;
   GenWTauMotherIndex = 0;
   GenWTauNumDaught = 0;
   GenWTauPdgId = 0;
   GenWTauStatus = 0;
   GenWTauTauDecayMode = 0;
   GenZTauMotherIndex = 0;
   GenZTauNumDaught = 0;
   GenZTauPdgId = 0;
   GenZTauStatus = 0;
   GenZTauTauDecayMode = 0;
   HPSTauCharge = 0;
   HPSTauDecayMode = 0;
   HPSTauIsCaloTau = 0;
   HPSTauIsPFTau = 0;
   HPSTauVtxIndex = 0;
   MuonBestTrackVtxIndex = 0;
   MuonCharge = 0;
   MuonCocktailCharge = 0;
   MuonCocktailRefitID = 0;
   MuonCocktailTrkHits = 0;
   MuonGlobalTrkValidHits = 0;
   MuonIsGlobal = 0;
   MuonIsPF = 0;
   MuonIsTracker = 0;
   MuonPassID = 0;
   MuonPixelHits = 0;
   MuonSegmentMatches = 0;
   MuonStationMatches = 0;
   MuonTrackLayersWithMeasurement = 0;
   MuonTrkHits = 0;
   MuonTrkHitsTrackerOnly = 0;
   MuonTrkPixelHits = 0;
   MuonVtxIndex = 0;
   PFCandChargeLeptLink = 0;
   PFJetBestVertexTrackAssociationIndex = 0;
   PFJetChargedHadronMultiplicity = 0;
   PFJetChargedMultiplicity = 0;
   PFJetClosestVertex3DIndex = 0;
   PFJetClosestVertexXYIndex = 0;
   PFJetClosestVertexZIndex = 0;
   PFJetElectronMultiplicity = 0;
   PFJetHFEMMultiplicity = 0;
   PFJetHFHadronMultiplicity = 0;
   PFJetMuonMultiplicity = 0;
   PFJetNConstituents = 0;
   PFJetNeutralHadronMultiplicity = 0;
   PFJetNeutralMultiplicity = 0;
   PFJetPartonFlavour = 0;
   PFJetPassLooseID = 0;
   PFJetPassTightID = 0;
   PFJetPhotonMultiplicity = 0;
   PhotonNTracksConvPhot = 0;
   HLTInsideDatasetTriggerPrescales = 0;
   HLTOutsideDatasetTriggerPrescales = 0;
   L1PhysBits = 0;
   L1TechBits = 0;
   VertexNTracks = 0;
   VertexNTracksW05 = 0;
   HLTFilterObjId = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;

   // do we need this at all now that we loop over the branches                                                                                                          
   // one-by-one?                                                                                                                                                        

   // Remove friends if any, for better performance                                                                                                                      
   bool skipFriends = true; // can be made configurable                                                                                                                  
   if( skipFriends ) {
     TList* flist = tree->GetListOfFriends();
     TIter nextf( flist );
     TFriendElement* fe = 0;
     while( ( fe = ( TFriendElement* ) nextf() ) ) {
       flist->Remove( fe );
       delete fe;
       fe = 0;
     }
   }


   bool deleteIndex = true; // can be made configurable                                                                                                                  
   if( deleteIndex ) {
     if( tree->GetTreeIndex() ) {
    
       tree->SetTreeIndex( 0 );
       delete tree->GetTreeIndex();
     }
   }

   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);
   fChain->SetMaxVirtualSize(0); 

   bool setall(true);

   if(setall){
     //Chain->SetBranchAddress("HLTKey", &HLTKey, &b_HLTKey);
     fChain->SetBranchAddress("HLTInsideDatasetTriggerNames", &HLTInsideDatasetTriggerNames, &b_HLTInsideDatasetTriggerNames);     
     //fChain->SetBranchAddress("HLTOutsideDatasetTriggerNames", &HLTOutsideDatasetTriggerNames, &b_HLTOutsideDatasetTriggerNames);
     //fChain->SetBranchAddress("HLTFilterName", &HLTFilterName, &b_HLTFilterName);
   }

   fChain->SetBranchAddress("isData", &isData, &b_isData);   
   //fChain->SetAutoDelete(kTRUE);
   if(setall){
     //fChain->SetBranchAddress("isBPTX0", &isBPTX0, &b_isBPTX0);
     //fChain->SetBranchAddress("isBSCBeamHalo", &isBSCBeamHalo, &b_isBSCBeamHalo);
     //fChain->SetBranchAddress("isBSCMinBias", &isBSCMinBias, &b_isBSCMinBias);
     //fChain->SetBranchAddress("isBeamScraping", &isBeamScraping, &b_isBeamScraping);
   }

   fChain->SetBranchAddress("isPhysDeclared", &isPhysDeclared, &b_isPhysDeclared);
   fChain->SetBranchAddress("isPrimaryVertex", &isPrimaryVertex, &b_isPrimaryVertex);
   fChain->SetBranchAddress("isTrackingFailure", &isTrackingFailure, &b_isTrackingFailure);
   fChain->SetBranchAddress("passBadEESupercrystalFilter", &passBadEESupercrystalFilter, &b_passBadEESupercrystalFilter);
   fChain->SetBranchAddress("passBeamHaloFilterLoose", &passBeamHaloFilterLoose, &b_passBeamHaloFilterLoose);

   if(setall){
     //fChain->SetBranchAddress("passBeamHaloFilterTight", &passBeamHaloFilterTight, &b_passBeamHaloFilterTight);
     //   fChain->SetBranchAddress("passCaloBoundaryDRFilter", &passCaloBoundaryDRFilter, &b_passCaloBoundaryDRFilter);
   }

   fChain->SetBranchAddress("passEcalDeadCellBoundaryEnergyFilter", &passEcalDeadCellBoundaryEnergyFilter, &b_passEcalDeadCellBoundaryEnergyFilter);
   fChain->SetBranchAddress("passEcalDeadCellTriggerPrimitiveFilter", &passEcalDeadCellTriggerPrimitiveFilter, &b_passEcalDeadCellTriggerPrimitiveFilter);
   fChain->SetBranchAddress("passEcalLaserCorrFilter", &passEcalLaserCorrFilter, &b_passEcalLaserCorrFilter);

   if(setall){
     //fChain->SetBranchAddress("passEcalMaskedCellDRFilter", &passEcalMaskedCellDRFilter, &b_passEcalMaskedCellDRFilter);
     //fChain->SetBranchAddress("passHBHENoiseFilter", &passHBHENoiseFilter, &b_passHBHENoiseFilter);
     //fChain->SetBranchAddress("passLogErrorTooManyClusters", &passLogErrorTooManyClusters, &b_passLogErrorTooManyClusters);
     //fChain->SetBranchAddress("passManyStripClus53X", &passManyStripClus53X, &b_passManyStripClus53X);
     //fChain->SetBranchAddress("passTooManyStripClus53X", &passTooManyStripClus53X, &b_passTooManyStripClus53X);
     //fChain->SetBranchAddress("passTrackingFailureFilter", &passTrackingFailureFilter, &b_passTrackingFailureFilter);
     //fChain->SetBranchAddress("hasVeryForwardPFMuon", &hasVeryForwardPFMuon, &b_hasVeryForwardPFMuon);
     //fChain->SetBranchAddress("hasJetWithBadUnc", &hasJetWithBadUnc, &b_hasJetWithBadUnc);
   }
   
   fChain->SetBranchAddress("ElectronGsfCtfCharge", &ElectronGsfCtfCharge, &b_ElectronGsfCtfCharge);
   fChain->SetBranchAddress("ElectronGsfCtfScPixCharge", &ElectronGsfCtfScPixCharge, &b_ElectronGsfCtfScPixCharge);
   if(setall){
     //fChain->SetBranchAddress("ElectronGsfScPixCharge", &ElectronGsfScPixCharge, &b_ElectronGsfScPixCharge);
     //   fChain->SetBranchAddress("ElectronHLTDoubleEleMatched", &ElectronHLTDoubleEleMatched, &b_ElectronHLTDoubleEleMatched);
     //fChain->SetBranchAddress("ElectronHLTSingleEleMatched", &ElectronHLTSingleEleMatched, &b_ElectronHLTSingleEleMatched);
     //fChain->SetBranchAddress("ElectronHLTSingleEleWP80Matched", &ElectronHLTSingleEleWP80Matched, &b_ElectronHLTSingleEleWP80Matched);
   }
   
   fChain->SetBranchAddress("ElectronHasEcalDrivenSeed", &ElectronHasEcalDrivenSeed, &b_ElectronHasEcalDrivenSeed);
   fChain->SetBranchAddress("ElectronHasMatchedConvPhot", &ElectronHasMatchedConvPhot, &b_ElectronHasMatchedConvPhot);
   fChain->SetBranchAddress("ElectronHasTrackerDrivenSeed", &ElectronHasTrackerDrivenSeed, &b_ElectronHasTrackerDrivenSeed);
   fChain->SetBranchAddress("ElectronIsEB", &ElectronIsEB, &b_ElectronIsEB);
   fChain->SetBranchAddress("ElectronIsEE", &ElectronIsEE, &b_ElectronIsEE);

   if(setall){
     //fChain->SetBranchAddress("MuonHLTSingleIsoMuonMatched", &MuonHLTSingleIsoMuonMatched, &b_MuonHLTSingleIsoMuonMatched);
     //fChain->SetBranchAddress("MuonHLTSingleMuonMatched", &MuonHLTSingleMuonMatched, &b_MuonHLTSingleMuonMatched);
     //fChain->SetBranchAddress("PhotonHasMatchedConvPhot", &PhotonHasMatchedConvPhot, &b_PhotonHasMatchedConvPhot);
     //fChain->SetBranchAddress("PhotonHasMatchedPromptEle", &PhotonHasMatchedPromptEle, &b_PhotonHasMatchedPromptEle);
     //fChain->SetBranchAddress("PhotonHasPixelSeed", &PhotonHasPixelSeed, &b_PhotonHasPixelSeed);
     //fChain->SetBranchAddress("PhotonIsEBEEGap", &PhotonIsEBEEGap, &b_PhotonIsEBEEGap);
     //fChain->SetBranchAddress("PhotonIsEBGap", &PhotonIsEBGap, &b_PhotonIsEBGap);
     //fChain->SetBranchAddress("PhotonIsEEGap", &PhotonIsEEGap, &b_PhotonIsEEGap);
     fChain->SetBranchAddress("HLTInsideDatasetTriggerDecisions", &HLTInsideDatasetTriggerDecisions, &b_HLTInsideDatasetTriggerDecisions);
     //fChain->SetBranchAddress("HLTOutsideDatasetTriggerDecisions", &HLTOutsideDatasetTriggerDecisions, &b_HLTOutsideDatasetTriggerDecisions);
     
     //fChain->SetBranchAddress("rhoForHEEP", &rhoForHEEP, &b_rhoForHEEP);
     //fChain->SetBranchAddress("rhoJetsCCPU", &rhoJetsCCPU, &b_rhoJetsCCPU);
     //fChain->SetBranchAddress("rhoJetsCN", &rhoJetsCN, &b_rhoJetsCN);
     //fChain->SetBranchAddress("rhoJetsCNT", &rhoJetsCNT, &b_rhoJetsCNT);
     //fChain->SetBranchAddress("time", &time, &b_time);
     //fChain->SetBranchAddress("PtHat", &PtHat, &b_PtHat);
   }
   fChain->SetBranchAddress("rhoJets", &rhoJets, &b_rhoJets);
   fChain->SetBranchAddress("Weight", &Weight, &b_Weight);      
   fChain->SetBranchAddress("CaloMET", &CaloMET, &b_CaloMET);
   fChain->SetBranchAddress("CaloSumET", &CaloSumET, &b_CaloSumET);
   fChain->SetBranchAddress("CaloMETPhiType1Cor", &CaloMETPhiType1Cor, &b_CaloMETPhiType1Cor);
   fChain->SetBranchAddress("CaloMETPhiUncorrType1Cor", &CaloMETPhiUncorrType1Cor, &b_CaloMETPhiUncorrType1Cor);
   fChain->SetBranchAddress("CaloMETType1Cor", &CaloMETType1Cor, &b_CaloMETType1Cor);
   fChain->SetBranchAddress("CaloMETUncorrType1Cor", &CaloMETUncorrType1Cor, &b_CaloMETUncorrType1Cor);
   fChain->SetBranchAddress("CaloSumETType1Cor", &CaloSumETType1Cor, &b_CaloSumETType1Cor);
   fChain->SetBranchAddress("CaloSumETUncorrType1Cor", &CaloSumETUncorrType1Cor, &b_CaloSumETUncorrType1Cor);
   fChain->SetBranchAddress("ElectronBeamSpotDXY", &ElectronBeamSpotDXY, &b_ElectronBeamSpotDXY);
   fChain->SetBranchAddress("ElectronBeamSpotDXYError", &ElectronBeamSpotDXYError, &b_ElectronBeamSpotDXYError);
   fChain->SetBranchAddress("ElectronCaloEnergy", &ElectronCaloEnergy, &b_ElectronCaloEnergy);
   fChain->SetBranchAddress("ElectronConvFitProb", &ElectronConvFitProb, &b_ElectronConvFitProb);
   fChain->SetBranchAddress("ElectronDCotTheta", &ElectronDCotTheta, &b_ElectronDCotTheta);
   fChain->SetBranchAddress("ElectronDeltaEtaTrkSC", &ElectronDeltaEtaTrkSC, &b_ElectronDeltaEtaTrkSC);
   fChain->SetBranchAddress("ElectronDeltaPhiTrkSC", &ElectronDeltaPhiTrkSC, &b_ElectronDeltaPhiTrkSC);
   fChain->SetBranchAddress("ElectronDist", &ElectronDist, &b_ElectronDist);
   fChain->SetBranchAddress("ElectronE1x5OverE5x5", &ElectronE1x5OverE5x5, &b_ElectronE1x5OverE5x5);
   fChain->SetBranchAddress("ElectronE2x5OverE5x5", &ElectronE2x5OverE5x5, &b_ElectronE2x5OverE5x5);
   fChain->SetBranchAddress("ElectronESuperClusterOverP", &ElectronESuperClusterOverP, &b_ElectronESuperClusterOverP);
   fChain->SetBranchAddress("ElectronEcalIsoDR03", &ElectronEcalIsoDR03, &b_ElectronEcalIsoDR03);
   fChain->SetBranchAddress("ElectronEcalIsoPAT", &ElectronEcalIsoPAT, &b_ElectronEcalIsoPAT);     
   fChain->SetBranchAddress("ElectronEnergy", &ElectronEnergy, &b_ElectronEnergy);
   fChain->SetBranchAddress("ElectronEta", &ElectronEta, &b_ElectronEta);
   if(setall){
     //fChain->SetBranchAddress("ElectronFbrem", &ElectronFbrem, &b_ElectronFbrem);
     //fChain->SetBranchAddress("ElectronHLTDoubleEleMatchEta", &ElectronHLTDoubleEleMatchEta, &b_ElectronHLTDoubleEleMatchEta);
     //fChain->SetBranchAddress("ElectronHLTDoubleEleMatchPhi", &ElectronHLTDoubleEleMatchPhi, &b_ElectronHLTDoubleEleMatchPhi);
     //fChain->SetBranchAddress("ElectronHLTDoubleEleMatchPt", &ElectronHLTDoubleEleMatchPt, &b_ElectronHLTDoubleEleMatchPt);
     //   fChain->SetBranchAddress("ElectronHLTSingleEleMatchEta", &ElectronHLTSingleEleMatchEta, &b_ElectronHLTSingleEleMatchEta);
     //fChain->SetBranchAddress("ElectronHLTSingleEleMatchPhi", &ElectronHLTSingleEleMatchPhi, &b_ElectronHLTSingleEleMatchPhi);
     //fChain->SetBranchAddress("ElectronHLTSingleEleMatchPt", &ElectronHLTSingleEleMatchPt, &b_ElectronHLTSingleEleMatchPt);
     //fChain->SetBranchAddress("ElectronHLTSingleEleWP80MatchEta", &ElectronHLTSingleEleWP80MatchEta, &b_ElectronHLTSingleEleWP80MatchEta);
     //   fChain->SetBranchAddress("ElectronHLTSingleEleWP80MatchPhi", &ElectronHLTSingleEleWP80MatchPhi, &b_ElectronHLTSingleEleWP80MatchPhi);
     //   fChain->SetBranchAddress("ElectronHLTSingleEleWP80MatchPt", &ElectronHLTSingleEleWP80MatchPt, &b_ElectronHLTSingleEleWP80MatchPt);
     //   fChain->SetBranchAddress("ElectronHcalIsoD1DR03", &ElectronHcalIsoD1DR03, &b_ElectronHcalIsoD1DR03);
     //   fChain->SetBranchAddress("ElectronHcalIsoD2DR03", &ElectronHcalIsoD2DR03, &b_ElectronHcalIsoD2DR03);
     //   fChain->SetBranchAddress("ElectronHcalIsoDR03", &ElectronHcalIsoDR03, &b_ElectronHcalIsoDR03);
     //fChain->SetBranchAddress("ElectronHcalIsoDR03FullCone", &ElectronHcalIsoDR03FullCone, &b_ElectronHcalIsoDR03FullCone);
     //   fChain->SetBranchAddress("ElectronHcalIsoPAT", &ElectronHcalIsoPAT, &b_ElectronHcalIsoPAT);
     //fChain->SetBranchAddress("ElectronLeadVtxDistXY", &ElectronLeadVtxDistXY, &b_ElectronLeadVtxDistXY);
     //fChain->SetBranchAddress("ElectronLeadVtxDistZ", &ElectronLeadVtxDistZ, &b_ElectronLeadVtxDistZ);
     //fChain->SetBranchAddress("ElectronMatchedGenParticleEta", &ElectronMatchedGenParticleEta, &b_ElectronMatchedGenParticleEta);
     //fChain->SetBranchAddress("ElectronMatchedGenParticlePhi", &ElectronMatchedGenParticlePhi, &b_ElectronMatchedGenParticlePhi);
     //fChain->SetBranchAddress("ElectronMatchedGenParticlePt", &ElectronMatchedGenParticlePt, &b_ElectronMatchedGenParticlePt);
   }
  
   fChain->SetBranchAddress("ElectronHoE", &ElectronHoE, &b_ElectronHoE);
   fChain->SetBranchAddress("ElectronPFChargedHadronIso03", &ElectronPFChargedHadronIso03, &b_ElectronPFChargedHadronIso03);
   fChain->SetBranchAddress("ElectronPFChargedHadronIso04", &ElectronPFChargedHadronIso04, &b_ElectronPFChargedHadronIso04);
   fChain->SetBranchAddress("ElectronPFNeutralHadronIso03", &ElectronPFNeutralHadronIso03, &b_ElectronPFNeutralHadronIso03);
   fChain->SetBranchAddress("ElectronPFNeutralHadronIso04", &ElectronPFNeutralHadronIso04, &b_ElectronPFNeutralHadronIso04);
   fChain->SetBranchAddress("ElectronPFPhotonIso03", &ElectronPFPhotonIso03, &b_ElectronPFPhotonIso03);
   fChain->SetBranchAddress("ElectronPFPhotonIso04", &ElectronPFPhotonIso04, &b_ElectronPFPhotonIso04);
   fChain->SetBranchAddress("ElectronPhi", &ElectronPhi, &b_ElectronPhi);     
   fChain->SetBranchAddress("ElectronPrimaryVertexDXY", &ElectronPrimaryVertexDXY, &b_ElectronPrimaryVertexDXY);
   fChain->SetBranchAddress("ElectronPrimaryVertexDXYError", &ElectronPrimaryVertexDXYError, &b_ElectronPrimaryVertexDXYError);
   fChain->SetBranchAddress("ElectronPt", &ElectronPt, &b_ElectronPt);
   
   if(setall){
     //fChain->SetBranchAddress("ElectronPtHeep", &ElectronPtHeep, &b_ElectronPtHeep);
     //fChain->SetBranchAddress("ElectronRelIsoPAT", &ElectronRelIsoPAT, &b_ElectronRelIsoPAT);
   }
   fChain->SetBranchAddress("ElectronSCEta", &ElectronSCEta, &b_ElectronSCEta);
   fChain->SetBranchAddress("ElectronSCPhi", &ElectronSCPhi, &b_ElectronSCPhi);
   fChain->SetBranchAddress("ElectronSCPt", &ElectronSCPt, &b_ElectronSCPt);
   fChain->SetBranchAddress("ElectronSCRawEnergy", &ElectronSCRawEnergy, &b_ElectronSCRawEnergy);
   fChain->SetBranchAddress("ElectronSigmaEtaEta", &ElectronSigmaEtaEta, &b_ElectronSigmaEtaEta);     
   fChain->SetBranchAddress("ElectronSigmaIEtaIEta", &ElectronSigmaIEtaIEta, &b_ElectronSigmaIEtaIEta);
   fChain->SetBranchAddress("ElectronTrackPt", &ElectronTrackPt, &b_ElectronTrackPt);
   fChain->SetBranchAddress("ElectronTrackValidFractionOfHits", &ElectronTrackValidFractionOfHits, &b_ElectronTrackValidFractionOfHits);
   fChain->SetBranchAddress("ElectronTrackVx", &ElectronTrackVx, &b_ElectronTrackVx);
   fChain->SetBranchAddress("ElectronTrackVy", &ElectronTrackVy, &b_ElectronTrackVy);
   fChain->SetBranchAddress("ElectronTrackVz", &ElectronTrackVz, &b_ElectronTrackVz);
   if(setall){
     //fChain->SetBranchAddress("ElectronTrkIsoDR03", &ElectronTrkIsoDR03, &b_ElectronTrkIsoDR03);
     //fChain->SetBranchAddress("ElectronTrkIsoPAT", &ElectronTrkIsoPAT, &b_ElectronTrkIsoPAT);
   }
   fChain->SetBranchAddress("ElectronVtxDistXY", &ElectronVtxDistXY, &b_ElectronVtxDistXY);
   fChain->SetBranchAddress("ElectronVtxDistZ", &ElectronVtxDistZ, &b_ElectronVtxDistZ);
   if(setall){
     ////fChain->SetBranchAddress("GenWElectronEnergy", &GenWElectronEnergy, &b_GenWElectronEnergy);
     //fChain->SetBranchAddress("GenWElectronEta", &GenWElectronEta, &b_GenWElectronEta);
     //fChain->SetBranchAddress("GenWElectronP", &GenWElectronP, &b_GenWElectronP);
     //fChain->SetBranchAddress("GenWElectronPhi", &GenWElectronPhi, &b_GenWElectronPhi);
     //fChain->SetBranchAddress("GenWElectronPt", &GenWElectronPt, &b_GenWElectronPt);
     //fChain->SetBranchAddress("GenWElectronPx", &GenWElectronPx, &b_GenWElectronPx);
     //fChain->SetBranchAddress("GenWElectronPy", &GenWElectronPy, &b_GenWElectronPy);
     //fChain->SetBranchAddress("GenWElectronPz", &GenWElectronPz, &b_GenWElectronPz);
     //fChain->SetBranchAddress("GenWElectronTauVisibleEta", &GenWElectronTauVisibleEta, &b_GenWElectronTauVisibleEta);
     //fChain->SetBranchAddress("GenWElectronTauVisiblePhi", &GenWElectronTauVisiblePhi, &b_GenWElectronTauVisiblePhi);
     //fChain->SetBranchAddress("GenWElectronTauVisiblePt", &GenWElectronTauVisiblePt, &b_GenWElectronTauVisiblePt);
     //fChain->SetBranchAddress("GenWElectronVX", &GenWElectronVX, &b_GenWElectronVX);
     //fChain->SetBranchAddress("GenWElectronVY", &GenWElectronVY, &b_GenWElectronVY);
     //fChain->SetBranchAddress("GenWElectronVZ", &GenWElectronVZ, &b_GenWElectronVZ);
     //fChain->SetBranchAddress("GenZElectronEnergy", &GenZElectronEnergy, &b_GenZElectronEnergy);
     //fChain->SetBranchAddress("GenZElectronEta", &GenZElectronEta, &b_GenZElectronEta);
     //fChain->SetBranchAddress("GenZElectronP", &GenZElectronP, &b_GenZElectronP);
     //fChain->SetBranchAddress("GenZElectronPhi", &GenZElectronPhi, &b_GenZElectronPhi);
     //fChain->SetBranchAddress("GenZElectronPt", &GenZElectronPt, &b_GenZElectronPt);
     //fChain->SetBranchAddress("GenZElectronPx", &GenZElectronPx, &b_GenZElectronPx);
     //fChain->SetBranchAddress("GenZElectronPy", &GenZElectronPy, &b_GenZElectronPy);
     //fChain->SetBranchAddress("GenZElectronPz", &GenZElectronPz, &b_GenZElectronPz);
     //fChain->SetBranchAddress("GenZElectronTauVisibleEta", &GenZElectronTauVisibleEta, &b_GenZElectronTauVisibleEta);
     //fChain->SetBranchAddress("GenZElectronTauVisiblePhi", &GenZElectronTauVisiblePhi, &b_GenZElectronTauVisiblePhi);
     //fChain->SetBranchAddress("GenZElectronTauVisiblePt", &GenZElectronTauVisiblePt, &b_GenZElectronTauVisiblePt);
     //fChain->SetBranchAddress("GenZElectronVX", &GenZElectronVX, &b_GenZElectronVX);
     //fChain->SetBranchAddress("GenZElectronVY", &GenZElectronVY, &b_GenZElectronVY);
     //fChain->SetBranchAddress("GenZElectronVZ", &GenZElectronVZ, &b_GenZElectronVZ);
     //fChain->SetBranchAddress("PDFCTEQWeights", &PDFCTEQWeights, &b_PDFCTEQWeights);
     //fChain->SetBranchAddress("PDFMSTWWeights", &PDFMSTWWeights, &b_PDFMSTWWeights);
     //fChain->SetBranchAddress("PDFNNPDFWeights", &PDFNNPDFWeights, &b_PDFNNPDFWeights);
     //fChain->SetBranchAddress("GenJetEMF", &GenJetEMF, &b_GenJetEMF);
     //fChain->SetBranchAddress("GenJetEnergy", &GenJetEnergy, &b_GenJetEnergy);
     //fChain->SetBranchAddress("GenJetEta", &GenJetEta, &b_GenJetEta);
     //fChain->SetBranchAddress("GenJetHADF", &GenJetHADF, &b_GenJetHADF);
     //fChain->SetBranchAddress("GenJetP", &GenJetP, &b_GenJetP);
     //fChain->SetBranchAddress("GenJetPhi", &GenJetPhi, &b_GenJetPhi);
     //fChain->SetBranchAddress("GenJetPt", &GenJetPt, &b_GenJetPt);
     //fChain->SetBranchAddress("GenMETCalo", &GenMETCalo, &b_GenMETCalo);
     //fChain->SetBranchAddress("GenMETPhiCalo", &GenMETPhiCalo, &b_GenMETPhiCalo);
     //fChain->SetBranchAddress("GenSumETCalo", &GenSumETCalo, &b_GenSumETCalo);
     //fChain->SetBranchAddress("GenMETPhiTrue", &GenMETPhiTrue, &b_GenMETPhiTrue);
     //fChain->SetBranchAddress("GenMETTrue", &GenMETTrue, &b_GenMETTrue);
     //fChain->SetBranchAddress("GenSumETTrue", &GenSumETTrue, &b_GenSumETTrue);
     //fChain->SetBranchAddress("GenWMuEnergy", &GenWMuEnergy, &b_GenWMuEnergy);
     //fChain->SetBranchAddress("GenWMuEta", &GenWMuEta, &b_GenWMuEta);
     //fChain->SetBranchAddress("GenWMuP", &GenWMuP, &b_GenWMuP);
     //fChain->SetBranchAddress("GenWMuPhi", &GenWMuPhi, &b_GenWMuPhi);
     //fChain->SetBranchAddress("GenWMuPt", &GenWMuPt, &b_GenWMuPt);
     //fChain->SetBranchAddress("GenWMuPx", &GenWMuPx, &b_GenWMuPx);
     //fChain->SetBranchAddress("GenWMuPy", &GenWMuPy, &b_GenWMuPy);
     //fChain->SetBranchAddress("GenWMuPz", &GenWMuPz, &b_GenWMuPz);
     //fChain->SetBranchAddress("GenWMuTauVisibleEta", &GenWMuTauVisibleEta, &b_GenWMuTauVisibleEta);
     //fChain->SetBranchAddress("GenWMuTauVisiblePhi", &GenWMuTauVisiblePhi, &b_GenWMuTauVisiblePhi);
     //fChain->SetBranchAddress("GenWMuTauVisiblePt", &GenWMuTauVisiblePt, &b_GenWMuTauVisiblePt);
     //fChain->SetBranchAddress("GenWMuVX", &GenWMuVX, &b_GenWMuVX);
     //fChain->SetBranchAddress("GenWMuVY", &GenWMuVY, &b_GenWMuVY);
     //fChain->SetBranchAddress("GenWMuVZ", &GenWMuVZ, &b_GenWMuVZ);
     //fChain->SetBranchAddress("GenZMuEnergy", &GenZMuEnergy, &b_GenZMuEnergy);
     //fChain->SetBranchAddress("GenZMuEta", &GenZMuEta, &b_GenZMuEta);
     //fChain->SetBranchAddress("GenZMuP", &GenZMuP, &b_GenZMuP);
     //fChain->SetBranchAddress("GenZMuPhi", &GenZMuPhi, &b_GenZMuPhi);
     //fChain->SetBranchAddress("GenZMuPt", &GenZMuPt, &b_GenZMuPt);
     //fChain->SetBranchAddress("GenZMuPx", &GenZMuPx, &b_GenZMuPx);
     //fChain->SetBranchAddress("GenZMuPy", &GenZMuPy, &b_GenZMuPy);
     //fChain->SetBranchAddress("GenZMuPz", &GenZMuPz, &b_GenZMuPz);
     //fChain->SetBranchAddress("GenZMuTauVisibleEta", &GenZMuTauVisibleEta, &b_GenZMuTauVisibleEta);
     //fChain->SetBranchAddress("GenZMuTauVisiblePhi", &GenZMuTauVisiblePhi, &b_GenZMuTauVisiblePhi);
     //fChain->SetBranchAddress("GenZMuTauVisiblePt", &GenZMuTauVisiblePt, &b_GenZMuTauVisiblePt);
     //fChain->SetBranchAddress("GenZMuVX", &GenZMuVX, &b_GenZMuVX);
     //fChain->SetBranchAddress("GenZMuVY", &GenZMuVY, &b_GenZMuVY);
     //fChain->SetBranchAddress("GenZMuVZ", &GenZMuVZ, &b_GenZMuVZ);
   }
   
   fChain->SetBranchAddress("GenParticleEnergy", &GenParticleEnergy, &b_GenParticleEnergy);
   fChain->SetBranchAddress("GenParticleEta", &GenParticleEta, &b_GenParticleEta);
   fChain->SetBranchAddress("GenParticleP", &GenParticleP, &b_GenParticleP);
   fChain->SetBranchAddress("GenParticlePhi", &GenParticlePhi, &b_GenParticlePhi);
   fChain->SetBranchAddress("GenParticlePt", &GenParticlePt, &b_GenParticlePt);
   fChain->SetBranchAddress("GenParticlePx", &GenParticlePx, &b_GenParticlePx);
   fChain->SetBranchAddress("GenParticlePy", &GenParticlePy, &b_GenParticlePy);
   fChain->SetBranchAddress("GenParticlePz", &GenParticlePz, &b_GenParticlePz);
   if(setall){
     //fChain->SetBranchAddress("GenParticleTauVisibleEta", &GenParticleTauVisibleEta, &b_GenParticleTauVisibleEta);
     //fChain->SetBranchAddress("GenParticleTauVisiblePhi", &GenParticleTauVisiblePhi, &b_GenParticleTauVisiblePhi);
     //fChain->SetBranchAddress("GenParticleTauVisiblePt", &GenParticleTauVisiblePt, &b_GenParticleTauVisiblePt);
     //fChain->SetBranchAddress("GenParticleVX", &GenParticleVX, &b_GenParticleVX);
     //fChain->SetBranchAddress("GenParticleVY", &GenParticleVY, &b_GenParticleVY);
     //fChain->SetBranchAddress("GenParticleVZ", &GenParticleVZ, &b_GenParticleVZ);
     //fChain->SetBranchAddress("GenWTauEnergy", &GenWTauEnergy, &b_GenWTauEnergy);
     //fChain->SetBranchAddress("GenWTauEta", &GenWTauEta, &b_GenWTauEta);
     //fChain->SetBranchAddress("GenWTauP", &GenWTauP, &b_GenWTauP);
     //fChain->SetBranchAddress("GenWTauPhi", &GenWTauPhi, &b_GenWTauPhi);
     //fChain->SetBranchAddress("GenWTauPt", &GenWTauPt, &b_GenWTauPt);
     //fChain->SetBranchAddress("GenWTauPx", &GenWTauPx, &b_GenWTauPx);
     //fChain->SetBranchAddress("GenWTauPy", &GenWTauPy, &b_GenWTauPy);
     //fChain->SetBranchAddress("GenWTauPz", &GenWTauPz, &b_GenWTauPz);
     //fChain->SetBranchAddress("GenWTauTauVisibleEta", &GenWTauTauVisibleEta, &b_GenWTauTauVisibleEta);
     //fChain->SetBranchAddress("GenWTauTauVisiblePhi", &GenWTauTauVisiblePhi, &b_GenWTauTauVisiblePhi);
     //fChain->SetBranchAddress("GenWTauTauVisiblePt", &GenWTauTauVisiblePt, &b_GenWTauTauVisiblePt);
     //fChain->SetBranchAddress("GenWTauVX", &GenWTauVX, &b_GenWTauVX);
     //fChain->SetBranchAddress("GenWTauVY", &GenWTauVY, &b_GenWTauVY);
     //fChain->SetBranchAddress("GenWTauVZ", &GenWTauVZ, &b_GenWTauVZ);
     //fChain->SetBranchAddress("GenZTauEnergy", &GenZTauEnergy, &b_GenZTauEnergy);
     //fChain->SetBranchAddress("GenZTauEta", &GenZTauEta, &b_GenZTauEta);
     //fChain->SetBranchAddress("GenZTauP", &GenZTauP, &b_GenZTauP);
     //fChain->SetBranchAddress("GenZTauPhi", &GenZTauPhi, &b_GenZTauPhi);
     //fChain->SetBranchAddress("GenZTauPt", &GenZTauPt, &b_GenZTauPt);
     //fChain->SetBranchAddress("GenZTauPx", &GenZTauPx, &b_GenZTauPx);
     //fChain->SetBranchAddress("GenZTauPy", &GenZTauPy, &b_GenZTauPy);
     //fChain->SetBranchAddress("GenZTauPz", &GenZTauPz, &b_GenZTauPz);
     //fChain->SetBranchAddress("GenZTauTauVisibleEta", &GenZTauTauVisibleEta, &b_GenZTauTauVisibleEta);
     //fChain->SetBranchAddress("GenZTauTauVisiblePhi", &GenZTauTauVisiblePhi, &b_GenZTauTauVisiblePhi);
     //fChain->SetBranchAddress("GenZTauTauVisiblePt", &GenZTauTauVisiblePt, &b_GenZTauTauVisiblePt);
     //fChain->SetBranchAddress("GenZTauVX", &GenZTauVX, &b_GenZTauVX);
     //fChain->SetBranchAddress("GenZTauVY", &GenZTauVY, &b_GenZTauVY);
     //fChain->SetBranchAddress("GenZTauVZ", &GenZTauVZ, &b_GenZTauVZ);
     //fChain->SetBranchAddress("HPSTauAgainstElectronDeadECALDiscr", &HPSTauAgainstElectronDeadECALDiscr, &b_HPSTauAgainstElectronDeadECALDiscr);
     //fChain->SetBranchAddress("HPSTauAgainstElectronLooseDiscr", &HPSTauAgainstElectronLooseDiscr, &b_HPSTauAgainstElectronLooseDiscr);
     //fChain->SetBranchAddress("HPSTauAgainstElectronLooseMVA2Discr", &HPSTauAgainstElectronLooseMVA2Discr, &b_HPSTauAgainstElectronLooseMVA2Discr);
     //fChain->SetBranchAddress("HPSTauAgainstElectronLooseMVA3Discr", &HPSTauAgainstElectronLooseMVA3Discr, &b_HPSTauAgainstElectronLooseMVA3Discr);
     //fChain->SetBranchAddress("HPSTauAgainstElectronMVA2categoryDiscr", &HPSTauAgainstElectronMVA2categoryDiscr, &b_HPSTauAgainstElectronMVA2categoryDiscr);
     //fChain->SetBranchAddress("HPSTauAgainstElectronMVA2rawDiscr", &HPSTauAgainstElectronMVA2rawDiscr, &b_HPSTauAgainstElectronMVA2rawDiscr);
     //fChain->SetBranchAddress("HPSTauAgainstElectronMVA3categoryDiscr", &HPSTauAgainstElectronMVA3categoryDiscr, &b_HPSTauAgainstElectronMVA3categoryDiscr);
     //fChain->SetBranchAddress("HPSTauAgainstElectronMVA3rawDiscr", &HPSTauAgainstElectronMVA3rawDiscr, &b_HPSTauAgainstElectronMVA3rawDiscr);
     //fChain->SetBranchAddress("HPSTauAgainstElectronMVADiscr", &HPSTauAgainstElectronMVADiscr, &b_HPSTauAgainstElectronMVADiscr);
     //fChain->SetBranchAddress("HPSTauAgainstElectronMediumDiscr", &HPSTauAgainstElectronMediumDiscr, &b_HPSTauAgainstElectronMediumDiscr);
     //fChain->SetBranchAddress("HPSTauAgainstElectronMediumMVA2Discr", &HPSTauAgainstElectronMediumMVA2Discr, &b_HPSTauAgainstElectronMediumMVA2Discr);
     //fChain->SetBranchAddress("HPSTauAgainstElectronMediumMVA3Discr", &HPSTauAgainstElectronMediumMVA3Discr, &b_HPSTauAgainstElectronMediumMVA3Discr);
     //fChain->SetBranchAddress("HPSTauAgainstElectronTightDiscr", &HPSTauAgainstElectronTightDiscr, &b_HPSTauAgainstElectronTightDiscr);
     //fChain->SetBranchAddress("HPSTauAgainstElectronTightMVA2Discr", &HPSTauAgainstElectronTightMVA2Discr, &b_HPSTauAgainstElectronTightMVA2Discr);
     //fChain->SetBranchAddress("HPSTauAgainstElectronTightMVA3Discr", &HPSTauAgainstElectronTightMVA3Discr, &b_HPSTauAgainstElectronTightMVA3Discr);
     //fChain->SetBranchAddress("HPSTauAgainstElectronVLooseMVA2Discr", &HPSTauAgainstElectronVLooseMVA2Discr, &b_HPSTauAgainstElectronVLooseMVA2Discr);
     //fChain->SetBranchAddress("HPSTauAgainstElectronVTightMVA3Discr", &HPSTauAgainstElectronVTightMVA3Discr, &b_HPSTauAgainstElectronVTightMVA3Discr);
     //fChain->SetBranchAddress("HPSTauAgainstMuonLoose2Discr", &HPSTauAgainstMuonLoose2Discr, &b_HPSTauAgainstMuonLoose2Discr);
     //fChain->SetBranchAddress("HPSTauAgainstMuonLooseDiscr", &HPSTauAgainstMuonLooseDiscr, &b_HPSTauAgainstMuonLooseDiscr);
     //fChain->SetBranchAddress("HPSTauAgainstMuonMedium2Discr", &HPSTauAgainstMuonMedium2Discr, &b_HPSTauAgainstMuonMedium2Discr);
     //fChain->SetBranchAddress("HPSTauAgainstMuonMediumDiscr", &HPSTauAgainstMuonMediumDiscr, &b_HPSTauAgainstMuonMediumDiscr);
     //fChain->SetBranchAddress("HPSTauAgainstMuonTight2Discr", &HPSTauAgainstMuonTight2Discr, &b_HPSTauAgainstMuonTight2Discr);
     //fChain->SetBranchAddress("HPSTauAgainstMuonTightDiscr", &HPSTauAgainstMuonTightDiscr, &b_HPSTauAgainstMuonTightDiscr);
     //fChain->SetBranchAddress("HPSTauBremsRecoveryEOverPLead", &HPSTauBremsRecoveryEOverPLead, &b_HPSTauBremsRecoveryEOverPLead);
     //fChain->SetBranchAddress("HPSTauCombinedIsolationDeltaBetaCorr3HitsDiscr", &HPSTauCombinedIsolationDeltaBetaCorr3HitsDiscr, &b_HPSTauCombinedIsolationDeltaBetaCorr3HitsDiscr);
     //fChain->SetBranchAddress("HPSTauDecayModeFindingDiscr", &HPSTauDecayModeFindingDiscr, &b_HPSTauDecayModeFindingDiscr);
     //fChain->SetBranchAddress("HPSTauEcalStripSumEOverPLead", &HPSTauEcalStripSumEOverPLead, &b_HPSTauEcalStripSumEOverPLead);
     //fChain->SetBranchAddress("HPSTauEmFraction", &HPSTauEmFraction, &b_HPSTauEmFraction);
     fChain->SetBranchAddress("HPSTauEt", &HPSTauEt, &b_HPSTauEt);
     fChain->SetBranchAddress("HPSTauEta", &HPSTauEta, &b_HPSTauEta);
     //fChain->SetBranchAddress("HPSTauEtaLeadCharged", &HPSTauEtaLeadCharged, &b_HPSTauEtaLeadCharged);
     //fChain->SetBranchAddress("HPSTauEtaetaMoment", &HPSTauEtaetaMoment, &b_HPSTauEtaetaMoment);
     //fChain->SetBranchAddress("HPSTauEtaphiMoment", &HPSTauEtaphiMoment, &b_HPSTauEtaphiMoment);
     //fChain->SetBranchAddress("HPSTauHcal3x3OverPLead", &HPSTauHcal3x3OverPLead, &b_HPSTauHcal3x3OverPLead);
     //fChain->SetBranchAddress("HPSTauHcalMaxOverPLead", &HPSTauHcalMaxOverPLead, &b_HPSTauHcalMaxOverPLead);
     //fChain->SetBranchAddress("HPSTauHcalTotOverPLead", &HPSTauHcalTotOverPLead, &b_HPSTauHcalTotOverPLead);
     //fChain->SetBranchAddress("HPSTauIsolationMVArawDiscr", &HPSTauIsolationMVArawDiscr, &b_HPSTauIsolationMVArawDiscr);
     //fChain->SetBranchAddress("HPSTauIsolationPFChargedHadrCandsPtSum", &HPSTauIsolationPFChargedHadrCandsPtSum, &b_HPSTauIsolationPFChargedHadrCandsPtSum);
     //fChain->SetBranchAddress("HPSTauIsolationPFGammaCandsEtSum", &HPSTauIsolationPFGammaCandsEtSum, &b_HPSTauIsolationPFGammaCandsEtSum);
     //fChain->SetBranchAddress("HPSTauLeadPFChargedHadrCandsignedSipt", &HPSTauLeadPFChargedHadrCandsignedSipt, &b_HPSTauLeadPFChargedHadrCandsignedSipt);
     //fChain->SetBranchAddress("HPSTauLeadVtxDistXY", &HPSTauLeadVtxDistXY, &b_HPSTauLeadVtxDistXY);
     //fChain->SetBranchAddress("HPSTauLeadVtxDistZ", &HPSTauLeadVtxDistZ, &b_HPSTauLeadVtxDistZ);
     //fChain->SetBranchAddress("HPSTauLooseCombinedIsolationDeltaBetaCorr3HitsDiscr", &HPSTauLooseCombinedIsolationDeltaBetaCorr3HitsDiscr, &b_HPSTauLooseCombinedIsolationDeltaBetaCorr3HitsDiscr);
     //fChain->SetBranchAddress("HPSTauLooseCombinedIsolationDeltaBetaCorrDiscr", &HPSTauLooseCombinedIsolationDeltaBetaCorrDiscr, &b_HPSTauLooseCombinedIsolationDeltaBetaCorrDiscr);
     //fChain->SetBranchAddress("HPSTauLooseIsolationDeltaBetaCorrDiscr", &HPSTauLooseIsolationDeltaBetaCorrDiscr, &b_HPSTauLooseIsolationDeltaBetaCorrDiscr);
     //fChain->SetBranchAddress("HPSTauLooseIsolationDiscr", &HPSTauLooseIsolationDiscr, &b_HPSTauLooseIsolationDiscr);
     //fChain->SetBranchAddress("HPSTauLooseIsolationMVA2Discr", &HPSTauLooseIsolationMVA2Discr, &b_HPSTauLooseIsolationMVA2Discr);
     //fChain->SetBranchAddress("HPSTauLooseIsolationMVADiscr", &HPSTauLooseIsolationMVADiscr, &b_HPSTauLooseIsolationMVADiscr);
     //fChain->SetBranchAddress("HPSTauMatchedGenJetEta", &HPSTauMatchedGenJetEta, &b_HPSTauMatchedGenJetEta);
     //fChain->SetBranchAddress("HPSTauMatchedGenJetPhi", &HPSTauMatchedGenJetPhi, &b_HPSTauMatchedGenJetPhi);
     //fChain->SetBranchAddress("HPSTauMatchedGenJetPt", &HPSTauMatchedGenJetPt, &b_HPSTauMatchedGenJetPt);
     //fChain->SetBranchAddress("HPSTauMatchedGenParticleEta", &HPSTauMatchedGenParticleEta, &b_HPSTauMatchedGenParticleEta);
     //fChain->SetBranchAddress("HPSTauMatchedGenParticlePhi", &HPSTauMatchedGenParticlePhi, &b_HPSTauMatchedGenParticlePhi);
     //fChain->SetBranchAddress("HPSTauMatchedGenParticlePt", &HPSTauMatchedGenParticlePt, &b_HPSTauMatchedGenParticlePt);
     //fChain->SetBranchAddress("HPSTauMaximumHCALPFClusterEt", &HPSTauMaximumHCALPFClusterEt, &b_HPSTauMaximumHCALPFClusterEt);
     //fChain->SetBranchAddress("HPSTauMediumCombinedIsolationDeltaBetaCorr3HitsDiscr", &HPSTauMediumCombinedIsolationDeltaBetaCorr3HitsDiscr, &b_HPSTauMediumCombinedIsolationDeltaBetaCorr3HitsDiscr);
     //fChain->SetBranchAddress("HPSTauMediumCombinedIsolationDeltaBetaCorrDiscr", &HPSTauMediumCombinedIsolationDeltaBetaCorrDiscr, &b_HPSTauMediumCombinedIsolationDeltaBetaCorrDiscr);
     //fChain->SetBranchAddress("HPSTauMediumIsolationDeltaBetaCorrDiscr", &HPSTauMediumIsolationDeltaBetaCorrDiscr, &b_HPSTauMediumIsolationDeltaBetaCorrDiscr);
     //fChain->SetBranchAddress("HPSTauMediumIsolationDiscr", &HPSTauMediumIsolationDiscr, &b_HPSTauMediumIsolationDiscr);
     //fChain->SetBranchAddress("HPSTauMediumIsolationMVA2Discr", &HPSTauMediumIsolationMVA2Discr, &b_HPSTauMediumIsolationMVA2Discr);
     //fChain->SetBranchAddress("HPSTauMediumIsolationMVADiscr", &HPSTauMediumIsolationMVADiscr, &b_HPSTauMediumIsolationMVADiscr);
     fChain->SetBranchAddress("HPSTauPhi", &HPSTauPhi, &b_HPSTauPhi);
     //fChain->SetBranchAddress("HPSTauPhiLeadCharged", &HPSTauPhiLeadCharged, &b_HPSTauPhiLeadCharged);
     //fChain->SetBranchAddress("HPSTauPhiphiMoment", &HPSTauPhiphiMoment, &b_HPSTauPhiphiMoment);
     fChain->SetBranchAddress("HPSTauPt", &HPSTauPt, &b_HPSTauPt);
     //fChain->SetBranchAddress("HPSTauPtLeadCharged", &HPSTauPtLeadCharged, &b_HPSTauPtLeadCharged);
     //fChain->SetBranchAddress("HPSTauSignalPFChargedHadrCandsCount", &HPSTauSignalPFChargedHadrCandsCount, &b_HPSTauSignalPFChargedHadrCandsCount);
     //fChain->SetBranchAddress("HPSTauSignalPFChargedHadrCandsEta", &HPSTauSignalPFChargedHadrCandsEta, &b_HPSTauSignalPFChargedHadrCandsEta);
     //fChain->SetBranchAddress("HPSTauSignalPFChargedHadrCandsPhi", &HPSTauSignalPFChargedHadrCandsPhi, &b_HPSTauSignalPFChargedHadrCandsPhi);
     //fChain->SetBranchAddress("HPSTauSignalPFChargedHadrCandsPt", &HPSTauSignalPFChargedHadrCandsPt, &b_HPSTauSignalPFChargedHadrCandsPt);
     //fChain->SetBranchAddress("HPSTauSignalPFGammaCandsCount", &HPSTauSignalPFGammaCandsCount, &b_HPSTauSignalPFGammaCandsCount);
     //fChain->SetBranchAddress("HPSTauSignalPFGammaCandsEta", &HPSTauSignalPFGammaCandsEta, &b_HPSTauSignalPFGammaCandsEta);
     //fChain->SetBranchAddress("HPSTauSignalPFGammaCandsPhi", &HPSTauSignalPFGammaCandsPhi, &b_HPSTauSignalPFGammaCandsPhi);
     //fChain->SetBranchAddress("HPSTauSignalPFGammaCandsPt", &HPSTauSignalPFGammaCandsPt, &b_HPSTauSignalPFGammaCandsPt);
     //fChain->SetBranchAddress("HPSTauSignalPFNeutrHadrCandsCount", &HPSTauSignalPFNeutrHadrCandsCount, &b_HPSTauSignalPFNeutrHadrCandsCount);
     //fChain->SetBranchAddress("HPSTauSignalPFNeutrHadrCandsEta", &HPSTauSignalPFNeutrHadrCandsEta, &b_HPSTauSignalPFNeutrHadrCandsEta);
     //fChain->SetBranchAddress("HPSTauSignalPFNeutrHadrCandsPhi", &HPSTauSignalPFNeutrHadrCandsPhi, &b_HPSTauSignalPFNeutrHadrCandsPhi);
     //fChain->SetBranchAddress("HPSTauSignalPFNeutrHadrCandsPt", &HPSTauSignalPFNeutrHadrCandsPt, &b_HPSTauSignalPFNeutrHadrCandsPt);
     //fChain->SetBranchAddress("HPSTauTightCombinedIsolationDeltaBetaCorr3HitsDiscr", &HPSTauTightCombinedIsolationDeltaBetaCorr3HitsDiscr, &b_HPSTauTightCombinedIsolationDeltaBetaCorr3HitsDiscr);
     //fChain->SetBranchAddress("HPSTauTightCombinedIsolationDeltaBetaCorrDiscr", &HPSTauTightCombinedIsolationDeltaBetaCorrDiscr, &b_HPSTauTightCombinedIsolationDeltaBetaCorrDiscr);
     //fChain->SetBranchAddress("HPSTauTightIsolationDeltaBetaCorrDiscr", &HPSTauTightIsolationDeltaBetaCorrDiscr, &b_HPSTauTightIsolationDeltaBetaCorrDiscr);
     //fChain->SetBranchAddress("HPSTauTightIsolationDiscr", &HPSTauTightIsolationDiscr, &b_HPSTauTightIsolationDiscr);
     //fChain->SetBranchAddress("HPSTauTightIsolationMVA2Discr", &HPSTauTightIsolationMVA2Discr, &b_HPSTauTightIsolationMVA2Discr);
     //fChain->SetBranchAddress("HPSTauTightIsolationMVADiscr", &HPSTauTightIsolationMVADiscr, &b_HPSTauTightIsolationMVADiscr);
     //fChain->SetBranchAddress("HPSTauVLooseCombinedIsolationDeltaBetaCorrDiscr", &HPSTauVLooseCombinedIsolationDeltaBetaCorrDiscr, &b_HPSTauVLooseCombinedIsolationDeltaBetaCorrDiscr);
     //fChain->SetBranchAddress("HPSTauVLooseIsolationDeltaBetaCorrDiscr", &HPSTauVLooseIsolationDeltaBetaCorrDiscr, &b_HPSTauVLooseIsolationDeltaBetaCorrDiscr);
     //fChain->SetBranchAddress("HPSTauVLooseIsolationDiscr", &HPSTauVLooseIsolationDiscr, &b_HPSTauVLooseIsolationDiscr);
     //fChain->SetBranchAddress("HPSTauVtxDistXY", &HPSTauVtxDistXY, &b_HPSTauVtxDistXY);
     //fChain->SetBranchAddress("HPSTauVtxDistZ", &HPSTauVtxDistZ, &b_HPSTauVtxDistZ);
     //fChain->SetBranchAddress("MuonBackToBackCompatibility", &MuonBackToBackCompatibility, &b_MuonBackToBackCompatibility);
     //fChain->SetBranchAddress("MuonBeamSpotDXY", &MuonBeamSpotDXY, &b_MuonBeamSpotDXY);
     //fChain->SetBranchAddress("MuonBeamSpotDXYError", &MuonBeamSpotDXYError, &b_MuonBeamSpotDXYError);
     //fChain->SetBranchAddress("MuonBestTrackVtxDistXY", &MuonBestTrackVtxDistXY, &b_MuonBestTrackVtxDistXY);
     //fChain->SetBranchAddress("MuonBestTrackVtxDistZ", &MuonBestTrackVtxDistZ, &b_MuonBestTrackVtxDistZ);
     //fChain->SetBranchAddress("MuonCocktailEta", &MuonCocktailEta, &b_MuonCocktailEta);
     //fChain->SetBranchAddress("MuonCocktailEtaError", &MuonCocktailEtaError, &b_MuonCocktailEtaError);
     //fChain->SetBranchAddress("MuonCocktailGlobalChi2", &MuonCocktailGlobalChi2, &b_MuonCocktailGlobalChi2);
     //fChain->SetBranchAddress("MuonCocktailP", &MuonCocktailP, &b_MuonCocktailP);
     //fChain->SetBranchAddress("MuonCocktailPhi", &MuonCocktailPhi, &b_MuonCocktailPhi);
     //fChain->SetBranchAddress("MuonCocktailPhiError", &MuonCocktailPhiError, &b_MuonCocktailPhiError);
     //fChain->SetBranchAddress("MuonCocktailPt", &MuonCocktailPt, &b_MuonCocktailPt);
     //fChain->SetBranchAddress("MuonCocktailPtError", &MuonCocktailPtError, &b_MuonCocktailPtError);
     //fChain->SetBranchAddress("MuonCocktailQOverPError", &MuonCocktailQOverPError, &b_MuonCocktailQOverPError);
     //fChain->SetBranchAddress("MuonCocktailTrkD0", &MuonCocktailTrkD0, &b_MuonCocktailTrkD0);
     //fChain->SetBranchAddress("MuonCocktailTrkD0Error", &MuonCocktailTrkD0Error, &b_MuonCocktailTrkD0Error);
     //fChain->SetBranchAddress("MuonCocktailTrkDz", &MuonCocktailTrkDz, &b_MuonCocktailTrkDz);
     //fChain->SetBranchAddress("MuonCocktailTrkDzError", &MuonCocktailTrkDzError, &b_MuonCocktailTrkDzError);
     //fChain->SetBranchAddress("MuonCocktailTrkValidFractionOfHits", &MuonCocktailTrkValidFractionOfHits, &b_MuonCocktailTrkValidFractionOfHits);
     //fChain->SetBranchAddress("MuonCosmicCompatibility", &MuonCosmicCompatibility, &b_MuonCosmicCompatibility);
   }
   fChain->SetBranchAddress("MuonEcalIso", &MuonEcalIso, &b_MuonEcalIso);
   fChain->SetBranchAddress("MuonEcalVetoIso", &MuonEcalVetoIso, &b_MuonEcalVetoIso);
   fChain->SetBranchAddress("MuonEnergy", &MuonEnergy, &b_MuonEnergy);
   fChain->SetBranchAddress("MuonEta", &MuonEta, &b_MuonEta);
   fChain->SetBranchAddress("MuonEtaError", &MuonEtaError, &b_MuonEtaError);
   fChain->SetBranchAddress("MuonGlobalChi2", &MuonGlobalChi2, &b_MuonGlobalChi2);
   fChain->SetBranchAddress("MuonHLTSingleIsoMuonMatchEta", &MuonHLTSingleIsoMuonMatchEta, &b_MuonHLTSingleIsoMuonMatchEta);
   fChain->SetBranchAddress("MuonHLTSingleIsoMuonMatchPhi", &MuonHLTSingleIsoMuonMatchPhi, &b_MuonHLTSingleIsoMuonMatchPhi);
   fChain->SetBranchAddress("MuonHLTSingleIsoMuonMatchPt", &MuonHLTSingleIsoMuonMatchPt, &b_MuonHLTSingleIsoMuonMatchPt);
   fChain->SetBranchAddress("MuonHLTSingleMuonMatchEta", &MuonHLTSingleMuonMatchEta, &b_MuonHLTSingleMuonMatchEta);
   fChain->SetBranchAddress("MuonHLTSingleMuonMatchPhi", &MuonHLTSingleMuonMatchPhi, &b_MuonHLTSingleMuonMatchPhi);
   fChain->SetBranchAddress("MuonHLTSingleMuonMatchPt", &MuonHLTSingleMuonMatchPt, &b_MuonHLTSingleMuonMatchPt);
   fChain->SetBranchAddress("MuonHOIso", &MuonHOIso, &b_MuonHOIso);
   fChain->SetBranchAddress("MuonHcalIso", &MuonHcalIso, &b_MuonHcalIso);
   fChain->SetBranchAddress("MuonHcalVetoIso", &MuonHcalVetoIso, &b_MuonHcalVetoIso);
   fChain->SetBranchAddress("MuonMatchedGenParticleEta", &MuonMatchedGenParticleEta, &b_MuonMatchedGenParticleEta);
   fChain->SetBranchAddress("MuonMatchedGenParticlePhi", &MuonMatchedGenParticlePhi, &b_MuonMatchedGenParticlePhi);
   fChain->SetBranchAddress("MuonMatchedGenParticlePt", &MuonMatchedGenParticlePt, &b_MuonMatchedGenParticlePt);
   fChain->SetBranchAddress("MuonOverlapCompatibility", &MuonOverlapCompatibility, &b_MuonOverlapCompatibility);
   fChain->SetBranchAddress("MuonP", &MuonP, &b_MuonP);
   fChain->SetBranchAddress("MuonPFIsoR03ChargedHadron", &MuonPFIsoR03ChargedHadron, &b_MuonPFIsoR03ChargedHadron);
   fChain->SetBranchAddress("MuonPFIsoR03ChargedParticle", &MuonPFIsoR03ChargedParticle, &b_MuonPFIsoR03ChargedParticle);
   fChain->SetBranchAddress("MuonPFIsoR03NeutralHadron", &MuonPFIsoR03NeutralHadron, &b_MuonPFIsoR03NeutralHadron);
   fChain->SetBranchAddress("MuonPFIsoR03NeutralHadronHT", &MuonPFIsoR03NeutralHadronHT, &b_MuonPFIsoR03NeutralHadronHT);
   fChain->SetBranchAddress("MuonPFIsoR03PU", &MuonPFIsoR03PU, &b_MuonPFIsoR03PU);
   fChain->SetBranchAddress("MuonPFIsoR03Photon", &MuonPFIsoR03Photon, &b_MuonPFIsoR03Photon);
   fChain->SetBranchAddress("MuonPFIsoR03PhotonHT", &MuonPFIsoR03PhotonHT, &b_MuonPFIsoR03PhotonHT);
   fChain->SetBranchAddress("MuonPFIsoR04ChargedHadron", &MuonPFIsoR04ChargedHadron, &b_MuonPFIsoR04ChargedHadron);
   fChain->SetBranchAddress("MuonPFIsoR04ChargedParticle", &MuonPFIsoR04ChargedParticle, &b_MuonPFIsoR04ChargedParticle);
   fChain->SetBranchAddress("MuonPFIsoR04NeutralHadron", &MuonPFIsoR04NeutralHadron, &b_MuonPFIsoR04NeutralHadron);
   fChain->SetBranchAddress("MuonPFIsoR04NeutralHadronHT", &MuonPFIsoR04NeutralHadronHT, &b_MuonPFIsoR04NeutralHadronHT);
   fChain->SetBranchAddress("MuonPFIsoR04PU", &MuonPFIsoR04PU, &b_MuonPFIsoR04PU);
   fChain->SetBranchAddress("MuonPFIsoR04Photon", &MuonPFIsoR04Photon, &b_MuonPFIsoR04Photon);
   fChain->SetBranchAddress("MuonPFIsoR04PhotonHT", &MuonPFIsoR04PhotonHT, &b_MuonPFIsoR04PhotonHT);
   fChain->SetBranchAddress("MuonPhi", &MuonPhi, &b_MuonPhi);
   fChain->SetBranchAddress("MuonPhiError", &MuonPhiError, &b_MuonPhiError);
   fChain->SetBranchAddress("MuonPrimaryVertexDXY", &MuonPrimaryVertexDXY, &b_MuonPrimaryVertexDXY);
   fChain->SetBranchAddress("MuonPrimaryVertexDXYError", &MuonPrimaryVertexDXYError, &b_MuonPrimaryVertexDXYError);
   fChain->SetBranchAddress("MuonPt", &MuonPt, &b_MuonPt);
   fChain->SetBranchAddress("MuonPtError", &MuonPtError, &b_MuonPtError);
   fChain->SetBranchAddress("MuonQOverPError", &MuonQOverPError, &b_MuonQOverPError);
   fChain->SetBranchAddress("MuonTimeCompatibility", &MuonTimeCompatibility, &b_MuonTimeCompatibility);
   fChain->SetBranchAddress("MuonTrackChi2", &MuonTrackChi2, &b_MuonTrackChi2);
   fChain->SetBranchAddress("MuonTrackerIsoSumPT", &MuonTrackerIsoSumPT, &b_MuonTrackerIsoSumPT);
   fChain->SetBranchAddress("MuonTrkD0", &MuonTrkD0, &b_MuonTrkD0);
   fChain->SetBranchAddress("MuonTrkD0Error", &MuonTrkD0Error, &b_MuonTrkD0Error);
   fChain->SetBranchAddress("MuonTrkDz", &MuonTrkDz, &b_MuonTrkDz);
   fChain->SetBranchAddress("MuonTrkDzError", &MuonTrkDzError, &b_MuonTrkDzError);
   fChain->SetBranchAddress("MuonTrkEta", &MuonTrkEta, &b_MuonTrkEta);
   fChain->SetBranchAddress("MuonTrkEtaError", &MuonTrkEtaError, &b_MuonTrkEtaError);
   fChain->SetBranchAddress("MuonTrkIso", &MuonTrkIso, &b_MuonTrkIso);
   fChain->SetBranchAddress("MuonTrkPhi", &MuonTrkPhi, &b_MuonTrkPhi);
   fChain->SetBranchAddress("MuonTrkPhiError", &MuonTrkPhiError, &b_MuonTrkPhiError);
   fChain->SetBranchAddress("MuonTrkPt", &MuonTrkPt, &b_MuonTrkPt);
   fChain->SetBranchAddress("MuonTrkPtError", &MuonTrkPtError, &b_MuonTrkPtError);
   fChain->SetBranchAddress("MuonTrkValidFractionOfHits", &MuonTrkValidFractionOfHits, &b_MuonTrkValidFractionOfHits);
   fChain->SetBranchAddress("MuonTrkVx", &MuonTrkVx, &b_MuonTrkVx);
   fChain->SetBranchAddress("MuonTrkVy", &MuonTrkVy, &b_MuonTrkVy);
   fChain->SetBranchAddress("MuonTrkVz", &MuonTrkVz, &b_MuonTrkVz);
   fChain->SetBranchAddress("MuonVtxDistXY", &MuonVtxDistXY, &b_MuonVtxDistXY);
   fChain->SetBranchAddress("MuonVtxDistZ", &MuonVtxDistZ, &b_MuonVtxDistZ);
   if(setall){
     //fChain->SetBranchAddress("PFCandEnergyLeptLink", &PFCandEnergyLeptLink, &b_PFCandEnergyLeptLink);
     //fChain->SetBranchAddress("PFCandEtaLeptLink", &PFCandEtaLeptLink, &b_PFCandEtaLeptLink);
     //fChain->SetBranchAddress("PFCandPhiLeptLink", &PFCandPhiLeptLink, &b_PFCandPhiLeptLink);
     //fChain->SetBranchAddress("PFCandPtLeptLink", &PFCandPtLeptLink, &b_PFCandPtLeptLink);
     //fChain->SetBranchAddress("PFJetBestVertexTrackAssociationFactor", &PFJetBestVertexTrackAssociationFactor, &b_PFJetBestVertexTrackAssociationFactor);
     //fChain->SetBranchAddress("PFJetBeta", &PFJetBeta, &b_PFJetBeta);
     //fChain->SetBranchAddress("PFJetBetaClassic", &PFJetBetaClassic, &b_PFJetBetaClassic);
     //fChain->SetBranchAddress("PFJetBetaStar", &PFJetBetaStar, &b_PFJetBetaStar);
     //fChain->SetBranchAddress("PFJetBetaStarClassic", &PFJetBetaStarClassic, &b_PFJetBetaStarClassic);
   }
   
   fChain->SetBranchAddress("PFJetChargedEmEnergyFraction", &PFJetChargedEmEnergyFraction, &b_PFJetChargedEmEnergyFraction);
   fChain->SetBranchAddress("PFJetChargedHadronEnergyFraction", &PFJetChargedHadronEnergyFraction, &b_PFJetChargedHadronEnergyFraction);
   fChain->SetBranchAddress("PFJetChargedMuEnergyFraction", &PFJetChargedMuEnergyFraction, &b_PFJetChargedMuEnergyFraction);
   fChain->SetBranchAddress("PFJetClosestVertexWeighted3DSeparation", &PFJetClosestVertexWeighted3DSeparation, &b_PFJetClosestVertexWeighted3DSeparation);
   fChain->SetBranchAddress("PFJetClosestVertexWeightedXYSeparation", &PFJetClosestVertexWeightedXYSeparation, &b_PFJetClosestVertexWeightedXYSeparation);
   fChain->SetBranchAddress("PFJetClosestVertexWeightedZSeparation", &PFJetClosestVertexWeightedZSeparation, &b_PFJetClosestVertexWeightedZSeparation);   
   fChain->SetBranchAddress("PFJetCombinedSecondaryVertexBTag", &PFJetCombinedSecondaryVertexBTag, &b_PFJetCombinedSecondaryVertexBTag);
   if(setall){
     ///  fChain->SetBranchAddress("PFJetCombinedInclusiveSecondaryVertexBTag", &PFJetCombinedInclusiveSecondaryVertexBTag, &b_PFJetCombinedInclusiveSecondaryVertexBTag);
     ///fChain->SetBranchAddress("PFJetCombinedMVABTag", &PFJetCombinedMVABTag, &b_PFJetCombinedMVABTag);
     //fChain->SetBranchAddress("PFJetCombinedSecondaryVertexMVABTag", &PFJetCombinedSecondaryVertexMVABTag, &b_PFJetCombinedSecondaryVertexMVABTag);
   }

   fChain->SetBranchAddress("PFJetElectronEnergyFraction", &PFJetElectronEnergyFraction, &b_PFJetElectronEnergyFraction);
   fChain->SetBranchAddress("PFJetEnergy", &PFJetEnergy, &b_PFJetEnergy);
   fChain->SetBranchAddress("PFJetEnergyRaw", &PFJetEnergyRaw, &b_PFJetEnergyRaw);
   fChain->SetBranchAddress("PFJetEta", &PFJetEta, &b_PFJetEta);
   fChain->SetBranchAddress("PFJetHFEMEnergyFraction", &PFJetHFEMEnergyFraction, &b_PFJetHFEMEnergyFraction);
   fChain->SetBranchAddress("PFJetHFHadronEnergyFraction", &PFJetHFHadronEnergyFraction, &b_PFJetHFHadronEnergyFraction);
   fChain->SetBranchAddress("PFJetJECUnc", &PFJetJECUnc, &b_PFJetJECUnc);
   fChain->SetBranchAddress("PFJetJetBProbabilityBTag", &PFJetJetBProbabilityBTag, &b_PFJetJetBProbabilityBTag);
   fChain->SetBranchAddress("PFJetJetProbabilityBTag", &PFJetJetProbabilityBTag, &b_PFJetJetProbabilityBTag);
   if(setall){
     //fChain->SetBranchAddress("PFJetL1FastJetJEC", &PFJetL1FastJetJEC, &b_PFJetL1FastJetJEC);
     //fChain->SetBranchAddress("PFJetL2L3ResJEC", &PFJetL2L3ResJEC, &b_PFJetL2L3ResJEC);
     //fChain->SetBranchAddress("PFJetL2RelJEC", &PFJetL2RelJEC, &b_PFJetL2RelJEC);
     //fChain->SetBranchAddress("PFJetL3AbsJEC", &PFJetL3AbsJEC, &b_PFJetL3AbsJEC);
   }
   fChain->SetBranchAddress("PFJetMuonEnergyFraction", &PFJetMuonEnergyFraction, &b_PFJetMuonEnergyFraction);
   fChain->SetBranchAddress("PFJetNeutralEmEnergyFraction", &PFJetNeutralEmEnergyFraction, &b_PFJetNeutralEmEnergyFraction);
   fChain->SetBranchAddress("PFJetNeutralHadronEnergyFraction", &PFJetNeutralHadronEnergyFraction, &b_PFJetNeutralHadronEnergyFraction);
   fChain->SetBranchAddress("PFJetPhi", &PFJetPhi, &b_PFJetPhi);
   fChain->SetBranchAddress("PFJetPhotonEnergyFraction", &PFJetPhotonEnergyFraction, &b_PFJetPhotonEnergyFraction);
   fChain->SetBranchAddress("PFJetPt", &PFJetPt, &b_PFJetPt);
   if(setall){
     //fChain->SetBranchAddress("PFJetPtRaw", &PFJetPtRaw, &b_PFJetPtRaw);
     // //fChain->SetBranchAddress("PFJetSimpleSecondaryVertexHighEffBTag", &PFJetSimpleSecondaryVertexHighEffBTag, &b_PFJetSimpleSecondaryVertexHighEffBTag);
     ///fChain->SetBranchAddress("PFJetSimpleSecondaryVertexHighPurBTag", &PFJetSimpleSecondaryVertexHighPurBTag, &b_PFJetSimpleSecondaryVertexHighPurBTag);
     //fChain->SetBranchAddress("PFJetSoftElectronByIP3dBTag", &PFJetSoftElectronByIP3dBTag, &b_PFJetSoftElectronByIP3dBTag);
     //fChain->SetBranchAddress("PFJetSoftElectronByPtBTag", &PFJetSoftElectronByPtBTag, &b_PFJetSoftElectronByPtBTag);
   }
   fChain->SetBranchAddress("PFJetSoftMuonBTag", &PFJetSoftMuonBTag, &b_PFJetSoftMuonBTag);
   fChain->SetBranchAddress("PFJetSoftMuonByIP3dBTag", &PFJetSoftMuonByIP3dBTag, &b_PFJetSoftMuonByIP3dBTag);
   fChain->SetBranchAddress("PFJetSoftMuonByPtBTag", &PFJetSoftMuonByPtBTag, &b_PFJetSoftMuonByPtBTag);
   fChain->SetBranchAddress("PFJetTrackCountingHighPurBTag", &PFJetTrackCountingHighPurBTag, &b_PFJetTrackCountingHighPurBTag);
   if(setall){
     //fChain->SetBranchAddress("PFJetTrackCountingHighEffBTag", &PFJetTrackCountingHighEffBTag, &b_PFJetTrackCountingHighEffBTag);
     // fChain->SetBranchAddress("PFMET", &PFMET, &b_PFMET);
   }
   fChain->SetBranchAddress("PFMETPhi", &PFMETPhi, &b_PFMETPhi);     
   fChain->SetBranchAddress("PFSumET", &PFSumET, &b_PFSumET);
   fChain->SetBranchAddress("PFMETPhiType01Cor", &PFMETPhiType01Cor, &b_PFMETPhiType01Cor);     
   fChain->SetBranchAddress("PFMETType01Cor", &PFMETType01Cor, &b_PFMETType01Cor);
   fChain->SetBranchAddress("PFSumETType01Cor", &PFSumETType01Cor, &b_PFSumETType01Cor);
   fChain->SetBranchAddress("PFMETPhiType01XYCor", &PFMETPhiType01XYCor, &b_PFMETPhiType01XYCor);
   fChain->SetBranchAddress("PFMETType01XYCor", &PFMETType01XYCor, &b_PFMETType01XYCor);
   fChain->SetBranchAddress("PFSumETType01XYCor", &PFSumETType01XYCor, &b_PFSumETType01XYCor);
   fChain->SetBranchAddress("PFMETPhiType1Cor", &PFMETPhiType1Cor, &b_PFMETPhiType1Cor);
   fChain->SetBranchAddress("PFMETType1Cor", &PFMETType1Cor, &b_PFMETType1Cor);
   fChain->SetBranchAddress("PFSumETType1Cor", &PFSumETType1Cor, &b_PFSumETType1Cor);
   if(setall){
     //fChain->SetBranchAddress("PhotonAlpha", &PhotonAlpha, &b_PhotonAlpha);
     //fChain->SetBranchAddress("PhotonChi2ConvPhot", &PhotonChi2ConvPhot, &b_PhotonChi2ConvPhot);
     //fChain->SetBranchAddress("PhotonDPhiTracksAtVtxConvPhot", &PhotonDPhiTracksAtVtxConvPhot, &b_PhotonDPhiTracksAtVtxConvPhot);
     //fChain->SetBranchAddress("PhotonDistOfMinApproachConvPhot", &PhotonDistOfMinApproachConvPhot, &b_PhotonDistOfMinApproachConvPhot);
     //fChain->SetBranchAddress("PhotonE2OverE9", &PhotonE2OverE9, &b_PhotonE2OverE9);
     //fChain->SetBranchAddress("PhotonE3x3", &PhotonE3x3, &b_PhotonE3x3);
     //fChain->SetBranchAddress("PhotonE4SwissCross", &PhotonE4SwissCross, &b_PhotonE4SwissCross);
     //fChain->SetBranchAddress("PhotonE5x5", &PhotonE5x5, &b_PhotonE5x5);
     //fChain->SetBranchAddress("PhotonEOverPConvPhot", &PhotonEOverPConvPhot, &b_PhotonEOverPConvPhot);
     //fChain->SetBranchAddress("PhotonEcalIsoDR03", &PhotonEcalIsoDR03, &b_PhotonEcalIsoDR03);
     //fChain->SetBranchAddress("PhotonEcalIsoDR04", &PhotonEcalIsoDR04, &b_PhotonEcalIsoDR04);
     //fChain->SetBranchAddress("PhotonEnergy", &PhotonEnergy, &b_PhotonEnergy);
     //fChain->SetBranchAddress("PhotonEta", &PhotonEta, &b_PhotonEta);
     //fChain->SetBranchAddress("PhotonHcalIsoDR03", &PhotonHcalIsoDR03, &b_PhotonHcalIsoDR03);
     //fChain->SetBranchAddress("PhotonHcalIsoDR03FullCone", &PhotonHcalIsoDR03FullCone, &b_PhotonHcalIsoDR03FullCone);
     //fChain->SetBranchAddress("PhotonHcalIsoDR04", &PhotonHcalIsoDR04, &b_PhotonHcalIsoDR04);
     //fChain->SetBranchAddress("PhotonHcalIsoDR04FullCone", &PhotonHcalIsoDR04FullCone, &b_PhotonHcalIsoDR04FullCone);
     //fChain->SetBranchAddress("PhotonHoE", &PhotonHoE, &b_PhotonHoE);
     //fChain->SetBranchAddress("PhotonNDofConvPhot", &PhotonNDofConvPhot, &b_PhotonNDofConvPhot);
     //fChain->SetBranchAddress("PhotonPairCotThetaSeparationConvPhot", &PhotonPairCotThetaSeparationConvPhot, &b_PhotonPairCotThetaSeparationConvPhot);
     //fChain->SetBranchAddress("PhotonPairInvariantMassConvPhot", &PhotonPairInvariantMassConvPhot, &b_PhotonPairInvariantMassConvPhot);
     //fChain->SetBranchAddress("PhotonPairMomentumxConvPhot", &PhotonPairMomentumxConvPhot, &b_PhotonPairMomentumxConvPhot);
     //fChain->SetBranchAddress("PhotonPairMomentumyConvPhot", &PhotonPairMomentumyConvPhot, &b_PhotonPairMomentumyConvPhot);
     //fChain->SetBranchAddress("PhotonPairMomentumzConvPhot", &PhotonPairMomentumzConvPhot, &b_PhotonPairMomentumzConvPhot);
     //fChain->SetBranchAddress("PhotonPhi", &PhotonPhi, &b_PhotonPhi);
     //fChain->SetBranchAddress("PhotonPt", &PhotonPt, &b_PhotonPt);
     //fChain->SetBranchAddress("PhotonSCenergy", &PhotonSCenergy, &b_PhotonSCenergy);
     //fChain->SetBranchAddress("PhotonSCeta", &PhotonSCeta, &b_PhotonSCeta);
     //fChain->SetBranchAddress("PhotonSCphi", &PhotonSCphi, &b_PhotonSCphi);
     //fChain->SetBranchAddress("PhotonSCseedEnergy", &PhotonSCseedEnergy, &b_PhotonSCseedEnergy);
     //fChain->SetBranchAddress("PhotonSEtaEta", &PhotonSEtaEta, &b_PhotonSEtaEta);
     //fChain->SetBranchAddress("PhotonSEtaPhi", &PhotonSEtaPhi, &b_PhotonSEtaPhi);
     //fChain->SetBranchAddress("PhotonSMajMaj", &PhotonSMajMaj, &b_PhotonSMajMaj);
     //fChain->SetBranchAddress("PhotonSMinMin", &PhotonSMinMin, &b_PhotonSMinMin);
     //fChain->SetBranchAddress("PhotonSPhiPhi", &PhotonSPhiPhi, &b_PhotonSPhiPhi);
     //fChain->SetBranchAddress("PhotonSigmaIEtaIEta", &PhotonSigmaIEtaIEta, &b_PhotonSigmaIEtaIEta);
     //fChain->SetBranchAddress("PhotonTimeSeed", &PhotonTimeSeed, &b_PhotonTimeSeed);
     //fChain->SetBranchAddress("PhotonTrkIsoHollowDR03", &PhotonTrkIsoHollowDR03, &b_PhotonTrkIsoHollowDR03);
     //fChain->SetBranchAddress("PhotonTrkIsoHollowDR04", &PhotonTrkIsoHollowDR04, &b_PhotonTrkIsoHollowDR04);
     //fChain->SetBranchAddress("PhotonTrkIsoSolidDR03", &PhotonTrkIsoSolidDR03, &b_PhotonTrkIsoSolidDR03);
     //fChain->SetBranchAddress("PhotonTrkIsoSolidDR04", &PhotonTrkIsoSolidDR04, &b_PhotonTrkIsoSolidDR04);
     //fChain->SetBranchAddress("PhotonXVtxConvPhot", &PhotonXVtxConvPhot, &b_PhotonXVtxConvPhot);
     //fChain->SetBranchAddress("PhotonYVtxConvPhot", &PhotonYVtxConvPhot, &b_PhotonYVtxConvPhot);
     //fChain->SetBranchAddress("PhotonZVtxConvPhot", &PhotonZVtxConvPhot, &b_PhotonZVtxConvPhot);
     //fChain->SetBranchAddress("TCMET", &TCMET, &b_TCMET);
     //fChain->SetBranchAddress("TCMETPhi", &TCMETPhi, &b_TCMETPhi);
     //fChain->SetBranchAddress("TCSumET", &TCSumET, &b_TCSumET);
   }

   fChain->SetBranchAddress("VertexIsFake", &VertexIsFake, &b_VertexIsFake);
   fChain->SetBranchAddress("VertexChi2", &VertexChi2, &b_VertexChi2);
   fChain->SetBranchAddress("VertexNDF", &VertexNDF, &b_VertexNDF);
   fChain->SetBranchAddress("VertexRho", &VertexRho, &b_VertexRho);
   fChain->SetBranchAddress("VertexX", &VertexX, &b_VertexX);
   fChain->SetBranchAddress("VertexXErr", &VertexXErr, &b_VertexXErr);
   fChain->SetBranchAddress("VertexY", &VertexY, &b_VertexY);
   fChain->SetBranchAddress("VertexYErr", &VertexYErr, &b_VertexYErr);
   fChain->SetBranchAddress("VertexZ", &VertexZ, &b_VertexZ);
   fChain->SetBranchAddress("VertexZErr", &VertexZErr, &b_VertexZErr);
   fChain->SetBranchAddress("PileUpInteractionsTrue", &PileUpInteractionsTrue, &b_PileUpInteractionsTrue);
   if(setall){
     //fChain->SetBranchAddress("HLTFilterObjEta", &HLTFilterObjEta, &b_HLTFilterObjEta);
     //fChain->SetBranchAddress("HLTFilterObjPhi", &HLTFilterObjPhi, &b_HLTFilterObjPhi);
     //fChain->SetBranchAddress("HLTFilterObjPt", &HLTFilterObjPt, &b_HLTFilterObjPt);
   }
   fChain->SetBranchAddress("ElectronCharge", &ElectronCharge, &b_ElectronCharge);
   fChain->SetBranchAddress("ElectronClassif", &ElectronClassif, &b_ElectronClassif);
   fChain->SetBranchAddress("ElectronMissingHits", &ElectronMissingHits, &b_ElectronMissingHits);
   fChain->SetBranchAddress("ElectronMissingHitsEG", &ElectronMissingHitsEG, &b_ElectronMissingHitsEG);
   fChain->SetBranchAddress("ElectronNumberOfBrems", &ElectronNumberOfBrems, &b_ElectronNumberOfBrems);
   fChain->SetBranchAddress("ElectronOverlaps", &ElectronOverlaps, &b_ElectronOverlaps);
   fChain->SetBranchAddress("ElectronPassEGammaIDEoP", &ElectronPassEGammaIDEoP, &b_ElectronPassEGammaIDEoP);
   fChain->SetBranchAddress("ElectronPassEGammaIDLoose", &ElectronPassEGammaIDLoose, &b_ElectronPassEGammaIDLoose);
   fChain->SetBranchAddress("ElectronPassEGammaIDMedium", &ElectronPassEGammaIDMedium, &b_ElectronPassEGammaIDMedium);
   fChain->SetBranchAddress("ElectronPassEGammaIDTight", &ElectronPassEGammaIDTight, &b_ElectronPassEGammaIDTight);
   fChain->SetBranchAddress("ElectronPassEGammaIDTrigTight", &ElectronPassEGammaIDTrigTight, &b_ElectronPassEGammaIDTrigTight);
   fChain->SetBranchAddress("ElectronPassEGammaIDTrigWP70", &ElectronPassEGammaIDTrigWP70, &b_ElectronPassEGammaIDTrigWP70);
   fChain->SetBranchAddress("ElectronPassEGammaIDVeto", &ElectronPassEGammaIDVeto, &b_ElectronPassEGammaIDVeto);
   fChain->SetBranchAddress("ElectronPassId", &ElectronPassId, &b_ElectronPassId);
   fChain->SetBranchAddress("ElectronPassIsoPAT", &ElectronPassIsoPAT, &b_ElectronPassIsoPAT);
   if(setall){
     //fChain->SetBranchAddress("ElectronVtxIndex", &ElectronVtxIndex, &b_ElectronVtxIndex);
     //fChain->SetBranchAddress("GenWElectronMotherIndex", &GenWElectronMotherIndex, &b_GenWElectronMotherIndex);
     //fChain->SetBranchAddress("GenWElectronNumDaught", &GenWElectronNumDaught, &b_GenWElectronNumDaught);
     //fChain->SetBranchAddress("GenWElectronPdgId", &GenWElectronPdgId, &b_GenWElectronPdgId);
     //fChain->SetBranchAddress("GenWElectronStatus", &GenWElectronStatus, &b_GenWElectronStatus);
     //fChain->SetBranchAddress("GenWElectronTauDecayMode", &GenWElectronTauDecayMode, &b_GenWElectronTauDecayMode);
     //fChain->SetBranchAddress("GenZElectronMotherIndex", &GenZElectronMotherIndex, &b_GenZElectronMotherIndex);
     //fChain->SetBranchAddress("GenZElectronNumDaught", &GenZElectronNumDaught, &b_GenZElectronNumDaught);
     //fChain->SetBranchAddress("GenZElectronPdgId", &GenZElectronPdgId, &b_GenZElectronPdgId);
     //fChain->SetBranchAddress("GenZElectronStatus", &GenZElectronStatus, &b_GenZElectronStatus);
     //fChain->SetBranchAddress("GenZElectronTauDecayMode", &GenZElectronTauDecayMode, &b_GenZElectronTauDecayMode);
     //fChain->SetBranchAddress("PileUpInteractions", &PileUpInteractions, &b_PileUpInteractions);
     //fChain->SetBranchAddress("PileUpOriginBX", &PileUpOriginBX, &b_PileUpOriginBX);
     //fChain->SetBranchAddress("GenWMuMotherIndex", &GenWMuMotherIndex, &b_GenWMuMotherIndex);
     //fChain->SetBranchAddress("GenWMuNumDaught", &GenWMuNumDaught, &b_GenWMuNumDaught);
     //fChain->SetBranchAddress("GenWMuPdgId", &GenWMuPdgId, &b_GenWMuPdgId);
     //fChain->SetBranchAddress("GenWMuStatus", &GenWMuStatus, &b_GenWMuStatus);
     //fChain->SetBranchAddress("GenWMuTauDecayMode", &GenWMuTauDecayMode, &b_GenWMuTauDecayMode);
     //fChain->SetBranchAddress("GenZMuMotherIndex", &GenZMuMotherIndex, &b_GenZMuMotherIndex);
     //fChain->SetBranchAddress("GenZMuNumDaught", &GenZMuNumDaught, &b_GenZMuNumDaught);
     //fChain->SetBranchAddress("GenZMuPdgId", &GenZMuPdgId, &b_GenZMuPdgId);
     //fChain->SetBranchAddress("GenZMuStatus", &GenZMuStatus, &b_GenZMuStatus);
     //fChain->SetBranchAddress("GenZMuTauDecayMode", &GenZMuTauDecayMode, &b_GenZMuTauDecayMode);
   }
   fChain->SetBranchAddress("GenParticleMotherIndex", &GenParticleMotherIndex, &b_GenParticleMotherIndex);
   fChain->SetBranchAddress("GenParticleNumDaught", &GenParticleNumDaught, &b_GenParticleNumDaught);
   fChain->SetBranchAddress("GenParticlePdgId", &GenParticlePdgId, &b_GenParticlePdgId);
   fChain->SetBranchAddress("GenParticleStatus", &GenParticleStatus, &b_GenParticleStatus);
   fChain->SetBranchAddress("GenParticleTauDecayMode", &GenParticleTauDecayMode, &b_GenParticleTauDecayMode);
   if(setall){
     //fChain->SetBranchAddress("GenWTauMotherIndex", &GenWTauMotherIndex, &b_GenWTauMotherIndex);
     //fChain->SetBranchAddress("GenWTauNumDaught", &GenWTauNumDaught, &b_GenWTauNumDaught);
     //fChain->SetBranchAddress("GenWTauPdgId", &GenWTauPdgId, &b_GenWTauPdgId);
     //fChain->SetBranchAddress("GenWTauStatus", &GenWTauStatus, &b_GenWTauStatus);
     //fChain->SetBranchAddress("GenWTauTauDecayMode", &GenWTauTauDecayMode, &b_GenWTauTauDecayMode);
     //fChain->SetBranchAddress("GenZTauMotherIndex", &GenZTauMotherIndex, &b_GenZTauMotherIndex);
     //fChain->SetBranchAddress("GenZTauNumDaught", &GenZTauNumDaught, &b_GenZTauNumDaught);
     //fChain->SetBranchAddress("GenZTauPdgId", &GenZTauPdgId, &b_GenZTauPdgId);
     //fChain->SetBranchAddress("GenZTauStatus", &GenZTauStatus, &b_GenZTauStatus);
     //fChain->SetBranchAddress("GenZTauTauDecayMode", &GenZTauTauDecayMode, &b_GenZTauTauDecayMode);
     fChain->SetBranchAddress("HPSTauCharge", &HPSTauCharge, &b_HPSTauCharge);
     //fChain->SetBranchAddress("HPSTauDecayMode", &HPSTauDecayMode, &b_HPSTauDecayMode);
     //fChain->SetBranchAddress("HPSTauIsCaloTau", &HPSTauIsCaloTau, &b_HPSTauIsCaloTau);
     fChain->SetBranchAddress("HPSTauIsPFTau", &HPSTauIsPFTau, &b_HPSTauIsPFTau);
     //fChain->SetBranchAddress("HPSTauVtxIndex", &HPSTauVtxIndex, &b_HPSTauVtxIndex);
     //fChain->SetBranchAddress("MuonBestTrackVtxIndex", &MuonBestTrackVtxIndex, &b_MuonBestTrackVtxIndex);
     //fChain->SetBranchAddress("MuonCocktailCharge", &MuonCocktailCharge, &b_MuonCocktailCharge);
     //fChain->SetBranchAddress("MuonCocktailRefitID", &MuonCocktailRefitID, &b_MuonCocktailRefitID);
     //fChain->SetBranchAddress("MuonCocktailTrkHits", &MuonCocktailTrkHits, &b_MuonCocktailTrkHits);
   }
   
   fChain->SetBranchAddress("MuonCharge", &MuonCharge, &b_MuonCharge);
   fChain->SetBranchAddress("MuonGlobalTrkValidHits", &MuonGlobalTrkValidHits, &b_MuonGlobalTrkValidHits);
   fChain->SetBranchAddress("MuonIsGlobal", &MuonIsGlobal, &b_MuonIsGlobal);
   fChain->SetBranchAddress("MuonIsPF", &MuonIsPF, &b_MuonIsPF);
   fChain->SetBranchAddress("MuonIsTracker", &MuonIsTracker, &b_MuonIsTracker);
   fChain->SetBranchAddress("MuonPassID", &MuonPassID, &b_MuonPassID);
   fChain->SetBranchAddress("MuonPixelHits", &MuonPixelHits, &b_MuonPixelHits);
   fChain->SetBranchAddress("MuonSegmentMatches", &MuonSegmentMatches, &b_MuonSegmentMatches);
   fChain->SetBranchAddress("MuonStationMatches", &MuonStationMatches, &b_MuonStationMatches);
   fChain->SetBranchAddress("MuonTrackLayersWithMeasurement", &MuonTrackLayersWithMeasurement, &b_MuonTrackLayersWithMeasurement);
   fChain->SetBranchAddress("MuonTrkHits", &MuonTrkHits, &b_MuonTrkHits);
   fChain->SetBranchAddress("MuonTrkHitsTrackerOnly", &MuonTrkHitsTrackerOnly, &b_MuonTrkHitsTrackerOnly);
   fChain->SetBranchAddress("MuonTrkPixelHits", &MuonTrkPixelHits, &b_MuonTrkPixelHits);
   fChain->SetBranchAddress("MuonVtxIndex", &MuonVtxIndex, &b_MuonVtxIndex);
   fChain->SetBranchAddress("PFCandChargeLeptLink", &PFCandChargeLeptLink, &b_PFCandChargeLeptLink);
   fChain->SetBranchAddress("PFJetBestVertexTrackAssociationIndex", &PFJetBestVertexTrackAssociationIndex, &b_PFJetBestVertexTrackAssociationIndex);
   fChain->SetBranchAddress("PFJetChargedHadronMultiplicity", &PFJetChargedHadronMultiplicity, &b_PFJetChargedHadronMultiplicity);
   fChain->SetBranchAddress("PFJetChargedMultiplicity", &PFJetChargedMultiplicity, &b_PFJetChargedMultiplicity);
   fChain->SetBranchAddress("PFJetClosestVertex3DIndex", &PFJetClosestVertex3DIndex, &b_PFJetClosestVertex3DIndex);
   fChain->SetBranchAddress("PFJetClosestVertexXYIndex", &PFJetClosestVertexXYIndex, &b_PFJetClosestVertexXYIndex);
   fChain->SetBranchAddress("PFJetClosestVertexZIndex", &PFJetClosestVertexZIndex, &b_PFJetClosestVertexZIndex);
   fChain->SetBranchAddress("PFJetElectronMultiplicity", &PFJetElectronMultiplicity, &b_PFJetElectronMultiplicity);
   fChain->SetBranchAddress("PFJetHFEMMultiplicity", &PFJetHFEMMultiplicity, &b_PFJetHFEMMultiplicity);
   fChain->SetBranchAddress("PFJetHFHadronMultiplicity", &PFJetHFHadronMultiplicity, &b_PFJetHFHadronMultiplicity);
   fChain->SetBranchAddress("PFJetMuonMultiplicity", &PFJetMuonMultiplicity, &b_PFJetMuonMultiplicity);
   fChain->SetBranchAddress("PFJetNConstituents", &PFJetNConstituents, &b_PFJetNConstituents);
   fChain->SetBranchAddress("PFJetNeutralHadronMultiplicity", &PFJetNeutralHadronMultiplicity, &b_PFJetNeutralHadronMultiplicity);
   fChain->SetBranchAddress("PFJetNeutralMultiplicity", &PFJetNeutralMultiplicity, &b_PFJetNeutralMultiplicity);
   fChain->SetBranchAddress("PFJetPartonFlavour", &PFJetPartonFlavour, &b_PFJetPartonFlavour);
   fChain->SetBranchAddress("PFJetPassLooseID", &PFJetPassLooseID, &b_PFJetPassLooseID);
   fChain->SetBranchAddress("PFJetPassTightID", &PFJetPassTightID, &b_PFJetPassTightID);
   fChain->SetBranchAddress("PFJetPhotonMultiplicity", &PFJetPhotonMultiplicity, &b_PFJetPhotonMultiplicity);
   fChain->SetBranchAddress("PhotonNTracksConvPhot", &PhotonNTracksConvPhot, &b_PhotonNTracksConvPhot);
   fChain->SetBranchAddress("HLTInsideDatasetTriggerPrescales", &HLTInsideDatasetTriggerPrescales, &b_HLTInsideDatasetTriggerPrescales);
   fChain->SetBranchAddress("HLTOutsideDatasetTriggerPrescales", &HLTOutsideDatasetTriggerPrescales, &b_HLTOutsideDatasetTriggerPrescales);
   if(setall){
     //fChain->SetBranchAddress("L1PhysBits", &L1PhysBits, &b_L1PhysBits);
     //fChain->SetBranchAddress("L1TechBits", &L1TechBits, &b_L1TechBits);
     //fChain->SetBranchAddress("VertexNTracksW05", &VertexNTracksW05, &b_VertexNTracksW05);
     //fChain->SetBranchAddress("HLTFilterObjId", &HLTFilterObjId, &b_HLTFilterObjId);
     //fChain->SetBranchAddress("bunch", &bunch, &b_bunch);
     //fChain->SetBranchAddress("ls", &ls, &b_ls); 
     //fChain->SetBranchAddress("orbit", &orbit, &b_orbit);                                                                                                    
     //fChain->SetBranchAddress("ProcessID", &ProcessID, &b_ProcessID);     
   }
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("VertexNTracks", &VertexNTracks, &b_VertexNTracks);
   fChain->SetBranchAddress("event", &event, &b_event);

      
   nentries = fChain->GetEntries();
   
   Notify();
}

void Data::setBranchStatus(void) {

  fChain->SetBranchStatus("*",1);


}

Bool_t Data::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Data::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

Int_t Data::Cut(Long64_t entry)
{
  fChain->GetEntry(entry);
  
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}

#endif