#!/usr/bin/tcsh

mkdir testdir
cd testdir

# Change process name in Matthias' config files

foreach process ( "Dijet" "HJ" "HJJ" "HW" "HWJ" "HZ" "HZJ" "ST_sch" "ST_tch" "ST_tch_4f" "ST_wtch_DR" "ST_wtch_DS" "VBF_H" "VBF_W-Z" "VBF_Wp_Wm" "VBF_Wp_Wp" "VBF_Z" "W" "W2jet" "WW" "WZ" "W_ew-BMNNP" "W_ew-BW" "Wbb" "Wj" "Wp_Wp_J_J" "Z" "Z2jet" "ZZ" "Z_ew-BMNNPV" "Zj" "Zjj" "dislepton" "gg_H" "gg_H_MSSM" "gg_H_quark-mass-effects" "hvq" "ttJ" "ttb_dec" )

  mkdir -p ${process}
  cd ${process}
  cat > runJob.sh <<EOF
#!/bin/bash
mkdir -p /tmp/covarell/${process}
cd /tmp/covarell/${process}
export SCRAM_ARCH=slc6_amd64_gcc481 
scram p CMSSW CMSSW_7_2_0_pre5
cd CMSSW_7_2_0_pre5/src 
cp $1/eiko_template.py ./${process}_cfg.py
cp $1/create_powheg_tarball.sh .
cp $1/runcmsgrid_powheg.sh .
sed s/TEMPLATE/${process}/g ${process}_cfg.py > ${process}_newcfg.py
mv ${process}_newcfg.py ${process}_cfg.py
eval \`scram runtime -sh\`
cmsRun ${process}_cfg.py >& log_${process}_cfg.log
mv log_${process}_cfg.log $1/testdir/${process}
EOF

  chmod 755 runJob.sh
  bsub -R "type=SLC6_64" -J ${process} -q 2nd < runJob.sh
  cd ..

end

cd ..
