set thisdata = $1
set thisfile = $2

foreach ll ( 1 2 3 4 5 6 7 8 ) 
  touch $thisfile
  echo "${thisdata}_$ll.log" >> $thisfile
  echo "  " >> $thisfile
  grep "hltSiStripMatchedRecHits" log/${thisdata}_$ll.log >> $thisfile
  # grep "hltSiPixelRecHits" log/${thisdata}_$ll.log >> $thisfile
  # grep "hltL1NonIsoStartUpElectronPixelSeeds" log/${thisdata}_$ll.log >> $thisfile
  echo "  " >> $thisfile 
end
