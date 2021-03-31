#!/bin/bash          

##Some general shell commands
#voms-proxy-init --voms cms --valid 192:0

#cp $X509_USER_PROXY /user/$USER/
STR="Hello World!"
echo $STR    
for counter in `cat fileList_parent.txt`; do

cp /user/amkalsi/CMSSW_8_0_17/src/Plots/ETauAnalysis_Test/python/${counter}.txt .

done
