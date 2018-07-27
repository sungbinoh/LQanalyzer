#sktree -a HN_pair_fake -s SKTree_LeptonSkim -userflag DijetFake -i QCD_Pt-15to20_MuEnriched -n 1
#sktree -a HN_pair_fake -s SKTree_LeptonSkim -userflag DijetFake -list qcd_mm -n 30
#sktree -a HN_pair_fake -s SKTree_LeptonSkim -userflag DijetFake -list qcd_ee -n 30

#sktree -a FR_syst_cal -S DoubleMuon -s SKTree_LeptonSkim -userflag DijetPrompt -n 20
#sktree -a Muon_FR_cal_all -s SKTree_LeptonSkim -userflag DijetPrompt -list fake_sub_big -n 40 -F
#sktree -a Muon_FR_cal_all -s SKTree_LeptonSkim -userflag DijetPrompt -list fake_sub_small -n 30

sktree -a Muon_FR_cal_all -S SingleMuon -s SKTree_LeptonSkim -userflag DijetPrompt -n 20
sktree -a Muon_FR_cal_all -s SKTree_LeptonSkim -userflag DijetPrompt -list fake_sub_big -n 40
sktree -a Muon_FR_cal_all -s SKTree_LeptonSkim -userflag DijetPrompt -list fake_sub_small -n 30

#$sktree -a Muon_FR_cal_all -s FLATCAT -userflag DijetPrompt -list fake_sub_big -n 40 -F
#sktree -a Muon_FR_cal_all -s FLATCAT -userflag DijetPrompt -list fake_sub_small -n 30
#sktree -a Muon_FR_cal_all -S SingleMuon -s FLATCAT -userflag DijetPrompt -n 20


#debug
#sktree -a FR_syst_cal -S DoubleMuon -s SKTree_LeptonSkim -p B -userflag DijetPrompt -n 1 -events 1
#sktree -a Muon_FR_cal_all -S SingleMuon -s SKTree_LeptonSkim -p B -userflag DijetPrompt -n 1