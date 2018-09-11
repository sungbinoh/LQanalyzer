#sktree -a Muon_FR_cal_all -s SKTree_LeptonSkim -userflag DijetPrompt -i LowStat_DYJets -c v8-0-7 -nskip 494187 -n 1 #-events 494400
#sktree -a Muon_FR_cal_all -s SKTree_LeptonSkim -userflag DijetPrompt -i DY50plus -n 1

#sktree -a Muon_FR_cal_all -s SKTree_LeptonSkim -userflag DijetPrompt -i LowStat_DYJets -c v8-0-7 -n 20 #50 -F
#sktree -a Muon_FR_cal_all -s SKTree_LeptonSkim -userflag DijetPrompt -i DY50plus -n 20 #50 -F

sktree -a Muon_FR_cal_all -s SKTree_LeptonSkim -userflag DijetPrompt -i TT_powheg  -c v8-0-7 -n 20
sktree -a Muon_FR_cal_all -s SKTree_LeptonSkim -userflag DijetPrompt -i TT_powheg -n 20 #50 -F

#sktreemaker -a SKTreeMaker -i DY_pt_50to100 -c v8-0-8 -n 20
#sktreemaker -a SKTreeMaker -i DY_pt_100to250 -c v8-0-8 -n 20
#sktreemaker -a SKTreeMaker -i DY_pt_250to400 -c v8-0-8 -n 20