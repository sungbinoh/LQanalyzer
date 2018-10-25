sktree -a Muon_FR_cal_all -S SingleMuon -s SKTree_LeptonSkim -userflag DijetPrompt -n 20
sktree -a Muon_FR_cal_all -s SKTree_LeptonSkim -userflag DijetPrompt -list fake_sub_big -n 40
sktree -a Muon_FR_cal_all -s SKTree_LeptonSkim -userflag DijetPrompt -list fake_sub_small -n 30

#sktree -a Muon_FR_cal_all -s FLATCAT -userflag DijetPrompt -i DY_pt_50to100 -n 40

#sktree -a Muon_FR_cal_all -s FLATCAT -userflag DijetPrompt -list DY_pt_binned -n 40
#sktree -a Muon_FR_cal_all -s SKTree_LeptonSkim -userflag DijetPrompt -i LowStat_DYJets -c v8-0-7 -n 40
#sktree -a Muon_FR_cal_all -S SingleMuon -s SKTree_LeptonSkim -userflag DijetPrompt -c v8-0-7 -n 20
#sktree -a Muon_FR_cal_all -s SKTree_LeptonSkim -userflag DijetPrompt -list fake_sub_big_807 -c v8-0-7 -n 40
#sktree -a Muon_FR_cal_all -s SKTree_LeptonSkim -userflag DijetPrompt -list fake_sub_small_807 -c v8-0-7 -n 30 -F

#sktree -a Electron_FR_cal_all -S SinglePhoton -s SKTree_LeptonSkim -userflag DijetPrompt -n 20

#sktree -a Muon_FR_cal_all -s FLATCAT -userflag DijetPrompt -list fake_sub_big -n 40
#sktree -a Muon_FR_cal_all -s FLATCAT -userflag DijetPrompt -list fake_sub_small -n 30
#sktree -a Muon_FR_cal_all -S SingleMuon -s FLATCAT -userflag DijetPrompt -n 20




#SinglePhoton 
sktree -a Electron_FR_cal_all -S SinglePhoton -s SKTree_LeptonSkim -userflag DijetPrompt -n 20 -F
sktree -a Electron_FR_cal_all -s SKTree_LeptonSkim -userflag DijetPrompt -list fake_sub_big -n 40
sktree -a Electron_FR_cal_all -s SKTree_LeptonSkim -userflag DijetPrompt -list fake_sub_small -n 30

#sktree -a Electron_FR_cal_all -S SinglePhoton -s FLATCAT -userflag DijetPrompt -n 20
#sktree -a Electron_FR_cal_all -s FLATCAT -userflag DijetPrompt -list fake_sub_big -n 40
#sktree -a Electron_FR_cal_all -s FLATCAT -userflag DijetPrompt -list fake_sub_small -n 30
#sktree -a Electron_FR_cal_all -s FLATCAT -userflag DijetPrompt -i DY50plus -n 40

#debug
#sktree -a Muon_FR_cal_all -s SKTree_LeptonSkim -userflag DijetPrompt -i DY50plus -n 1 # -events 1
#sktree -a FR_syst_cal -S DoubleMuon -s SKTree_LeptonSkim -p B -userflag DijetPrompt -n 1 -events 1
#sktree -a Muon_FR_cal_all -S SingleMuon -s SKTree_LeptonSkim -p B -userflag DijetPrompt -n 1
#sktree -a Electron_FR_cal_all -S SinglePhoton -s FLATCAT -p B -userflag DijetPrompt -n 1 #-events 10000
#sktree -a Electron_FR_cal_all -S SingleElectron -s SKTree_LeptonSkim -p B -userflag DijetPrompt -n 1 -events 100
