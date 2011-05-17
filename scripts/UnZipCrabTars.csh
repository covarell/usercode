#!/usr/bin/tcsh

# Untar CMSSW results in case of getoutput problems

foreach file (`ls out_*.tar`)    
  # gunzip ${file}
  # newfile = ${file//'tgz'/'tar'}
  tar -xvf ${file}
  rm -f $file
end

echo OK
 
