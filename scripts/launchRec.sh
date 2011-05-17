#!/usr/local/bin/zsh

# modify here ! ------------------------
nruns=445;  

inputDir="/afs/cern.ch/user/c/covarell/scratch0/jpsi-trig/CMSSW_2_1_12/src/HLTrigger/Configuration/python";
cmsswRelDir="/afs/cern.ch/user/c/covarell/scratch0/jpsi-trig/CMSSW_2_1_12/src";
castorDir="/castor/cern.ch/user/c/covarell/QCD_EMenriched";
sample="QCD_EMenr_3080_HLTskimmed";
# modify here ! ------------------------

writeScript() {

    jobNumber=`echo $1 | awk '{ printf("%0d",$1) }'`;
    
    # preparing the cfg file
    echo "preparing the cfg file" 
    echo "cp $inputDir/JPsieeRecoFromRaw_cfg.py ./JPsieeRecoFromRaw_1_${jobNumber}_cfg.py;"
    cp $inputDir/JPsieeRecoFromRaw_cfg.py ./JPsieeRecoFromRaw_1_${jobNumber}_cfg.py;
    
    sed -e "s/NUMBER/${jobNumber}/" JPsieeRecoFromRaw_1_${jobNumber}_cfg.py > JPsieeRecoFromRaw_${jobNumber}_cfg.py
    rm JPsieeRecoFromRaw_1_${jobNumber}_cfg.py
    
    # preparing the script to launch
    echo "#!/usr/local/bin/zsh

cd $cmsswRelDir;
eval \`scramv1 runtime -sh\`;
cd -;

cmsRun $inputDir/JPsieeRecoFromRaw_${jobNumber}_cfg.py &> log_$jobNumber.log;
mv ${sample}_HLTandReco.root ${sample}_HLTandReco_${jobNumber}.root

for x in *root 
do
   rfcp \$x $castorDir;
done
rm -f ${sample}_HLTandReco_${jobNumber}.root
" > myScript_${jobNumber}.sh;
    
    chmod 755 myScript_${jobNumber}.sh;
    echo "bsub -q 8nh -o log2_$jobNumber.log myScript_${jobNumber}.sh";
    bsub -q 8nh myScript_${jobNumber}.sh;    
}


for (( i=0; i<nruns; i++ )) {
   writeScript $i;
}
