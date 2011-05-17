set thisdata = $1
set thisfile = $2

foreach ll ( 1 2 3 4 5 6 7 8 ) 
  touch $thisfile
  echo "${thisdata}_$ll.log" >> $thisfile
  echo "  " >> $thisfile
  ## grep "hltL1NonIsoHLTNonIsoSingleElectronLWEt15PixelMatchFilter" log/${thisdata}_$ll.log >> $thisfile
  ## grep "hltL1NonIsoHLTNonIsoSingleElectronEt15PixelMatchFilter" log/${thisdata}_$ll.log >> $thisfile
  ## grep "hltL1NonIsoHLTNonIsoSingleElectronLWEt15DetaDphiFilter" log/${thisdata}_$ll.log >> $thisfile
  ## grep "hltL1NonIsoHLTNonIsoSingleElectronEt15DetaDphiFilter" log/${thisdata}_$ll.log >> $thisfile
  grep "hltL1NonIsoLargeWindowDoubleElectronTrackIsolFilter" log/${thisdata}_$ll.log >> $thisfile
  echo "  " >> $thisfile 
end
