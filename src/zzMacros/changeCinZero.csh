#!/usr/bin/tcsh

# Change process name in Matthias' config files

foreach file (`ls text/paramsPTOverMCJLST*.txt`)    
# echo sed 's/'\"'$1'\"'/'\"'$2'\"'/g' ${file} > ${file}new
 sed s/'C'/'+\/- 0.00'/g ${file} > ${file}new
 mv ${file}new ${file}
end

 
