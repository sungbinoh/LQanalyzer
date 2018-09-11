#signal
#sktree -a HN_pair_all -s FLATCAT -list hn_pair_ee -userflag hn_pair_ee -n 20 -SIG
#sktree -a HN_pair_all -s FLATCAT -list hn_pair_mm -userflag hn_pair_mm -n 20 -SIG

sktree -a HN_pair_all -s FLATCAT -i HNpair_MuMu_WR5000_Zp3000_HN900_official -n 20 -c v8-0-7 -SIG
sktree -a HN_pair_all -s FLATCAT -i ZprimetoNN_WR_MuMu_Z3000_N900 -n 20 -SIG

#sktree -a HN_pair_all -s FLATCAT -i HNpair_MuMu_WR5000_Zp3000_HN900_official,ZprimetoNN_WR_MuMu_Z3000_N900 -n 20 -c v8-0-7 -SIG



#debug

#sktree -a HN_pair_all -s FLATCAT -i HNpair_MuMu_WR5000_Zp3000_HN900_official -n 1 -c v8-0-7 -events 10 -SIG
#sktree -a HN_pair_all -s FLATCAT -i ZprimetoNN_WR_MuMu_Z2400_N1000 -n 1 -events 10 -SIG
