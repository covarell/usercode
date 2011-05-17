#!/usr/local/bin/zsh

## GENERIC SCRIPT TO LAUNCH HLT JOBS

# modify here ------------------------
nruns=100;  

inputDir="Configuration/GenProduction/python";
cmsswRelDir="/afs/cern.ch/user/c/covarell/scratch0/quarkonia/CMSSW_2_2_10/src";
castorDir="/castor/cern.ch/user/c/covarell/Jpsi-nonp-mumu";
sample="EVTGEN_inclBtoJpsiMuMu_10TeV_GEN_DIGI_L1_RAW_RECO_HLT_IDEAL";
# modify here ------------------------

writeScript() {

    jobNumber=`echo $1 | awk '{ printf("%0d",$1) }'`;
    
    # preparing the cfg file
    echo "preparing the cfg file" 
    echo "cp $cmsswRelDir/$inputDir/run_HLT_RECO_IDEAL.py ./run_HLT_RECO_IDEAL_temp_${jobNumber}.py;"
    cp $cmsswRelDir/$inputDir/run_HLT_RECO_IDEAL.py ./run_HLT_RECO_IDEAL_temp_${jobNumber}.py;
    
    sed -e "s/NUMBER/${jobNumber}/" run_HLT_RECO_IDEAL_temp_${jobNumber}.py > run_HLT_RECO_IDEAL_${jobNumber}.py
    rm run_HLT_RECO_IDEAL_temp_${jobNumber}.py
    
    # preparing the script to launch
    echo "#!/usr/local/bin/zsh

cd $cmsswRelDir;
eval \`scramv1 runtime -sh\`;
cd - 

cmsRun $cmsswRelDir/$inputDir/run_HLT_RECO_IDEAL_${jobNumber}.py &> log_$jobNumber.log;
mv ${sample}.root ${sample}_${jobNumber}.root

for x in *root 
do
   rfcp \$x $castorDir;
done
rm -f ${sample}_${jobNumber}.root
" > myScript_${jobNumber}.sh;
    
    chmod 755 myScript_${jobNumber}.sh;
    echo "bsub -q 8nh -o log2_$jobNumber.log myScript_${jobNumber}.sh";
    bsub -q 8nh myScript_${jobNumber}.sh;    
}


for (( i=1; i<nruns+1; i++ )) {
   writeScript $i;
}
