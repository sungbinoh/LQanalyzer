
#signal
#sktree -a HN_pair_MM -s FLATCAT -list hn_pair_mm -n 20 -SIG


#bkg MC
#sktree -a HN_pair_MM -s SKTree_DiLepSkim -list dilep_suoh -n 20
sktree -a HN_pair_MM -s SKTree_DiLepSkim -list fake_sub_big -n 40
sktree -a HN_pair_MM -s SKTree_DiLepSkim -list fake_sub_small -n 20


#bkg fake
#sktree -a HN_pair_MM -S DoubleMuon -s SKTree_DiLepSkim -userflag RunFake -n 40


#Data
#sktree -a HN_pair_MM -S DoubleMuon -s SKTree_DiLepSkim -n 40
sktree -a HN_pair_MM -S SingleMuon -s SKTree_LeptonSkim -n 30

#debug
#sktree -a HN_pair_MM -s FLATCAT -i HNpair_MuMu_WR5000_Zp3000_HN900_official -n 1 -SIG
#sktree -a HN_pair_MM -S SingleMuon -s SKTree_DiLepSkim -p B -n 1 -events 100000
#sktree -a HN_pair_MM -S SingleMuon -s SKTree_LeptonSkim -p B -n 1 -events 100000
#sktree -a HN_pair_MM -S DoubleMuon -s SKTree_LeptonSkim -p B -n 1 -events 1000
#sktree -a HN_pair_MM -s SKTree_DiLepSkim -i DYJets_10to50 -n 1 -events 100000