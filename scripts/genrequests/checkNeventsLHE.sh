#!/usr/bin/tcsh

# Change process name in Matthias' config files

foreach ie (0 1 2 3 4 5 6 7 8 9)    
cmsStage -f /store/lhe/$1/lhefileMerged_$ie.lhe .
grep '<event>' lhefileMerged_$ie.lhe | wc -l 
rm -rf lhefileMerged_$ie.lhe
end

echo Process name changed from $1 
 
