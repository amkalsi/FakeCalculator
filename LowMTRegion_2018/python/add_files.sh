#!/bin/bash
hadd Data.root OutputFiles/Con*data_Data*.root
hadd TTinc.root   OutputFiles/Con*TT_inc_div_*root
rm OutputFiles/Con*TT_inc_div_*root
hadd DYJets.root  OutputFiles/Con*DY*inc*root
hadd TT2L2Nu.root OutputFiles/Con*TT_inc*.root
hadd WWinc.root      OutputFiles/Con*WW_inc_div*.root
rm OutputFiles/Con*WW_inc_div*.root

#hadd  TT_1200.root OutputFiles/Con*TT_1200.root
#hadd  TT_1800.root OutputFiles/Con*TT_1800.root
#hadd  TT_500.root  OutputFiles/Con*TT_500.root
#hadd  TT_800.root  OutputFiles/Con*TT_800.root
hadd WW.root      OutputFiles/Con*WW_inc*.root
hadd WZ.root      OutputFiles/Con*ZZ.root
hadd ZZ.root      OutputFiles/Con*WZ.root
hadd ZZTo2L2Nu.root      OutputFiles/Con*ZZTo2L2Nu*.root    
hadd ZZTo2L2Q.root      OutputFiles/Con*ZZTo2L2Q*.root
hadd ZZTo4L.root         OutputFiles/Con*ZZTo4L*.root 
hadd WZTo2L2Q.root       OutputFiles/Con*WZTo2L2Q*.root  
hadd WZTo3LNu.root      OutputFiles/Con*WZTo3LNu*.root
hadd ST_tW_top_5f.root     OutputFiles/Con*ST_tW_top_5f*.root
hadd ST_tW_antitop_5f.root  OutputFiles/Con*ST_tW_antitop_5f*.root
#hadd  WW_600.root      OutputFiles/Con*WW_600_*root
#hadd  WW_1200.root     OutputFiles/Con*WW_1200*root
#hadd  WW_2500.root      OutputFiles/Con*WW_2500*.root
#hadd  WW_200.root      OutputFiles/Con*WW_200*root
hadd DY_1000.root     OutputFiles/Con*DY_1000.root
hadd DY_1500.root     OutputFiles/Con*DY_1500.root
hadd DY_2000.root     OutputFiles/Con*DY_2000.root
hadd DY_3000.root OutputFiles/Con*DY_3000.root
hadd DY_400.root  OutputFiles/Con*DY_400.root
hadd DY_500.root  OutputFiles/Con*DY_500.root
hadd DY_700.root  OutputFiles/Con*DY_700.root
hadd DY_800.root  OutputFiles/Con*DY_800.root
hadd DY_200.root  OutputFiles/Con*DY_200.root
hadd DY_100.root  OutputFiles/Con*DY_100.root
hadd  TT_1200.root OutputFiles/Con*TT_1200*.root
hadd  TT_1800.root OutputFiles/Con*TT_1800*.root
hadd  TT_500.root  OutputFiles/Con*TT_500*.root
hadd  TT_800.root  OutputFiles/Con*TT_800*.root
hadd  WW_600.root      OutputFiles/Con*WW_600*root
hadd  WW_1200.root     OutputFiles/Con*WW_1200*root
hadd  WW_2500.root      OutputFiles/Con*WW_2500*.root
hadd  WW_200.root      OutputFiles/Con*WW_200*root
