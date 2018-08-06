#sktree -a Muon_FR_cal_all -S SingleMuon -s SKTree_LeptonSkim -userflag DijetPrompt -n 20
#sktree -a Muon_FR_cal_all -s SKTree_LeptonSkim -userflag DijetPrompt -list fake_sub_big -n 40
#sktree -a Muon_FR_cal_all -s SKTree_LeptonSkim -userflag DijetPrompt -list fake_sub_small -n 30

#sktree -a Electron_FR_cal_all -S SinglePhoton -s SKTree_LeptonSkim -userflag DijetPrompt -n 20

#sktree -a Electron_FR_cal_all -S SinglePhoton -s FLATCAT -userflag DijetPrompt -n 20
#sktree -a Electron_FR_cal_all -s SKTree_LeptonSkim -userflag DijetPrompt -list fake_sub_big -n 40
#sktree -a Electron_FR_cal_all -s SKTree_LeptonSkim -userflag DijetPrompt -list fake_sub_small -n 30

#$sktree -a Muon_FR_cal_all -s FLATCAT -userflag DijetPrompt -list fake_sub_big -n 40 -F
#sktree -a Muon_FR_cal_all -s FLATCAT -userflag DijetPrompt -list fake_sub_small -n 30
#sktree -a Muon_FR_cal_all -S SingleMuon -s FLATCAT -userflag DijetPrompt -n 20



#SinglePhoton period by period
sktree -a Electron_FR_cal_all -S SinglePhoton -s SKTree_LeptonSkim -p B -userflag DijetPrompt -n 70 -F
sktree -a Electron_FR_cal_all -S SinglePhoton -s SKTree_LeptonSkim -p C -userflag DijetPrompt -n 70 -F
sktree -a Electron_FR_cal_all -S SinglePhoton -s SKTree_LeptonSkim -p D -userflag DijetPrompt -n 70 -F
sktree -a Electron_FR_cal_all -S SinglePhoton -s SKTree_LeptonSkim -p E -userflag DijetPrompt -n 70 -F
sktree -a Electron_FR_cal_all -S SinglePhoton -s SKTree_LeptonSkim -p F -userflag DijetPrompt -n 70 -F
sktree -a Electron_FR_cal_all -S SinglePhoton -s SKTree_LeptonSkim -p H_v2 -userflag DijetPrompt -n 70 -F
sktree -a Electron_FR_cal_all -S SinglePhoton -s SKTree_LeptonSkim -p H_v3 -userflag DijetPrompt -n 70 -F


#debug
#sktree -a FR_syst_cal -S DoubleMuon -s SKTree_LeptonSkim -p B -userflag DijetPrompt -n 1 -events 1
#sktree -a Muon_FR_cal_all -S SingleMuon -s SKTree_LeptonSkim -p B -userflag DijetPrompt -n 1
#sktree -a Electron_FR_cal_all -S SinglePhoton -s SKTree_LeptonSkim -p B -userflag DijetPrompt -n 1 -events 10000
#sktree -a Electron_FR_cal_all -S SingleElectron -s SKTree_LeptonSkim -p B -userflag DijetPrompt -n 1 -events 100
