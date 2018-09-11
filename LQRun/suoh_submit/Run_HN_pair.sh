#bkg MC
#sktree -a HN_pair_all -s SKTree_DiLepSkim -list fake_sub_big -n 100 #-F
#sktree -a HN_pair_all -s SKTree_DiLepSkim -list DY_amcnlo -userflag tau_veto -n 100 #-F
sktree -a HN_pair_all -s SKTree_DiLepSkim -list fake_sub_small -n 20 #-F
#sktree -a HN_pair_all -s SKTree_LeptonSkim -list DY_mass_binned -n 20

#bkg fake
#sktree -a HN_pair_all -S SingleMuon -s SKTree_DiLepSkim -userflag RunFake -o /data2/CAT_SKTreeOutput/JobOutPut/suoh/LQanalyzer/data/output/CAT/HN_pair_all/periodBtoH/Fake -n 40 -F
#sktree -a HN_pair_all -S DoubleEG -s SKTree_DiLepSkim -userflag RunFake -o /data2/CAT_SKTreeOutput/JobOutPut/suoh/LQanalyzer/data/output/CAT/HN_pair_all/periodBtoH/Fake -n 40 -F
#sktree -a HN_pair_all -S SingleMuon -s SKTree_DiLepSkim -p B -userflag RunFake -n 1 -events 1000000

#Data
sktree -a HN_pair_all -S SingleMuon -s SKTree_DiLepSkim -n 30 #-F
sktree -a HN_pair_all -S DoubleEG -s SKTree_DiLepSkim -n 30 #-F

#CF
#sktree -a CF_MC -s SKTree_LeptonSkim -list MC_CF -n 40
#sktree -a CF_MC -s SKTree_LeptonSkim -i DYJets -n 100 #-F
#sktree -a CF_MC -s SKTree_LeptonSkim -i DYJets_10to50 -n 50 #-F
#sktree -a CF_MC -s SKTree_LeptonSkim -i DYJets -n 1 -nskip 427000 -events 428000

#signal
#sktree -a HN_pair_all -s FLATCAT -list hn_pair_ee -userflag hn_pair_ee -n 20 -SIG
#sktree -a HN_pair_all -s FLATCAT -list hn_pair_mm -userflag hn_pair_mm -n 20 -SIG

#fake rate with QCD samples
#sktree -a HN_pair_fake -s SKTree_LeptonSkim -userflag DijetFake -list qcd_mm -n 30
#sktree -a HN_pair_fake -s SKTree_LeptonSkim -userflag DijetFake -list qcd_ee -n 30



#debug

#sktree -a HN_pair_MM -s FLATCAT -i HNpair_MuMu_WR5000_Zp3000_HN900_official -n 1 -SIG
#sktree -a HN_pair_all -s SKTree_DiLepSkim -i DYJets -n 1 -events 100000
#sktree -a HN_pair_all -S DoubleEG -s SKTree_DiLepSkim -p B -fake True -userflag RunFake -n 1 # -events 100000
#sktree -a HN_pair_all -S SingleMuon -s SKTree_LeptonSkim -p B -n 1 -events 1
#sktree -a HN_pair_all -S DoubleEG -s SKTree_LeptonSkim -p B -n 1 -events 1

#sktree -a HN_pair_MM -s FLATCAT -i DYJets -n 1 -events 1
#sktree -a HN_pair_MM -s SKTree_LeptonSkim -i DYJets -n 1 -events 1

#sktree -a HN_pair_MM -S SingleMuon -s SKTree_DiLepSkim -p B -n 1 -events 100000
#sktree -a HN_pair_MM -S SingleMuon -s SKTree_LeptonSkim -p B -n 1 -events 100000
#sktree -a HN_pair_MM -S DoubleMuon -s SKTree_LeptonSkim -p B -n 1 -events 1000
#sktree -a HN_pair_MM -s SKTree_DiLepSkim -i DYJets_10to50 -n 1 -events 100000