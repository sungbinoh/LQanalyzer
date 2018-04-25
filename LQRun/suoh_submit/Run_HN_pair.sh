
#signal
sktree -a HN_pair_MM -s FLATCAT -list hn_pair_mm -n 20 -SIG


#bkg MC
#sktree -a HN_pair_MM -s SKTree_DiLepSkim -list dilep_suoh -n 20



#bkg fake
#sktree -a HN_pair_MM -S DoubleMuon -s SKTree_DiLepSkim -userflag RunFake -n 40


#Data
#sktree -a HN_pair_MM -S DoubleMuon -s SKTree_DiLepSkim -n 40


#debug
#sktree -a HN_pair_MM -s FLATCAT -i HNpair_MuMu_WR5000_Zp3000_HN900_official -n 1 -SIG