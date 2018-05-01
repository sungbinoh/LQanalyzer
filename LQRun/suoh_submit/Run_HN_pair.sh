#signal
#sktree -a HN_pair_all -s FLATCAT -list hn_pair_mm -n 20 -SIG
sktree -a HN_pair_all -s FLATCAT -i HNpair_MuMu_WR5000_Zp4000_HN100_official -n 10 -SIG

#bkg MC
#sktree -a HN_pair_MM -s SKTree_DiLepSkim -list dilep_suoh -n 20
#sktree -a HN_pair_all -s SKTree_DiLepSkim -list fake_sub_big -n 40
#sktree -a HN_pair_all -s SKTree_DiLepSkim -list fake_sub_small -n 20
#sktree -a HN_pair_MM -s SKTree_LeptonSkim -i DYJets -n 60

#bkg fake
#sktree -a HN_pair_MM -S DoubleMuon -s SKTree_DiLepSkim -userflag RunFake -n 40


#Data
#sktree -a HN_pair_MM -S DoubleMuon -s SKTree_DiLepSkim -n 40
#sktree -a HN_pair_all -S SingleMuon -s SKTree_DiLepSkim -n 30
#sktree -a HN_pair_all -S ElectronMuon -s SKTree_DiLepSkim -n 30
#sktree -a HN_pair_all -S DoubleEG -s SKTree_DiLepSkim -n 30

#debug
#sktree -a HN_pair_MM -s FLATCAT -i HNpair_MuMu_WR5000_Zp3000_HN900_official -n 1 -SIG
#sktree -a HN_pair_all -s SKTree_DiLepSkim -i DYJets -n 1 -events 100000
#sktree -a HN_pair_MM -s FLATCAT -i DYJets -n 1 -events 1
#sktree -a HN_pair_MM -s SKTree_LeptonSkim -i DYJets -n 1 -events 1

#sktree -a HN_pair_MM -S SingleMuon -s SKTree_DiLepSkim -p B -n 1 -events 100000
#sktree -a HN_pair_MM -S SingleMuon -s SKTree_LeptonSkim -p B -n 1 -events 100000
#sktree -a HN_pair_MM -S DoubleMuon -s SKTree_LeptonSkim -p B -n 1 -events 1000
#sktree -a HN_pair_MM -s SKTree_DiLepSkim -i DYJets_10to50 -n 1 -events 100000