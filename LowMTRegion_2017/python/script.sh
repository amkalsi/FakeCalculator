#!/bin/bash          

##Some general shell commands
#voms-proxy-init --voms cms --valid 192:0

#cp $X509_USER_PROXY /user/$USER/
STR="Hello World!"
echo $STR    
for counter in `cat fileList.txt`; do

cat > file_${counter}.sh <<EOF
#!/bin/bash          

echo "step 1"
pwd=$PWD
echo $PWD
source $VO_CMS_SW_DIR/cmsset_default.sh
export X509_USER_PROXY=/user/$USER/x509up_u$(id -u $USER)

cd /user/amkalsi/CMSSW_8_0_17/src/Plots/ETauAnalysis/python/
eval \`scramv1 runtime -sh\`
cmsRun MuJetAnalysisBGZ.py inputFilename=${counter}

echo "COMPLETED -------"
EOF

qsub -q localgrid@cream02 -o log_${counter}.stdout  -e log_${counter}.stderr file_${counter}.sh 

done
