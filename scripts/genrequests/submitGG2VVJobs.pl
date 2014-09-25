#!/usr/bin/perl
###############################################################################
# Launch GG2VV jobs
###############################################################################

# job steering ----------------------------------------------------------------

$jobName = $ARGV[0];
$nJobs = $ARGV[1];
$theSeed = $ARGV[2];
$onlyHist = $ARGV[3];
$copyInit = $ARGV[4];
$initName = $ARGV[5];
$farm = $ARGV[6];

$gg2VVvers = "3.1.5";

if (!($onlyHist == 1 || $onlyHist == 0)) { 
  print("Argument n. 3 must be: 0 (generate events) or 1 (only cross-section) \n");
  exit;
} 

if ($onlyHist == 1) { $nJobs = 1; }

if (!($copyInit == 1 || $copyInit == 0)) { 
  print("Argument n. 4 must be: 0 (create init-new) or 1 (copy init) \n");
  exit;
} 

# name of job
$castorarea="/store/user/covarell/GG2VV-lhe";

# name of job
$basedir="/afs/cern.ch/work/c/covarell/gg2VV/";

# interactive or lxbatch queue
# $farm="I";
# $farm="2nd";
# $farm="1nw";

# site-specific paths etc

# -----------------------------------------------------------------------------

print "\n";
print "-------------------------------------------------------------------------------\n";
print "L A U N C H I N G   G G 2 V V   J O B S\n";
print "-------------------------------------------------------------------------------\n";
print "Workdir: ${basedir}/gg2VV-${gg2VVvers}/gg2VV \n";
print "EOS dir: ${castorarea}/${jobName}\n";
# -----------------------------------------------------------------------------
# create N jobs

if ($onlyHist == 0) {
    system("cmsMkdir ${castorarea}/${jobName}");
}
system("cd ${basedir}/gg2VV-${gg2VVvers}/gg2VV/");

$ijob=1;
while ( $ijob <= $nJobs ) {

    system("mkdir ${jobName}_$ijob");
    my $seed = $theSeed + 876*$ijob;
    open(OUTFILE,">${jobName}_$ijob/gg2VV.rngseed");
    print OUTFILE $seed; 
    close(OUTFILE);

    system("cp gg2VV ${jobName}_$ijob"); 

    if ($copyInit == 1) { 
	system("cp ${initName}.init ${jobName}_$ijob/gg2VV.init");
        system("cp ${initName}.evt ${jobName}_$ijob/gg2VV.evt");
    }

    $dir="./${jobName}_$ijob";
    make_subjob($dir,$ijob);

    print "Submit job $ijob ... ";
    $test=0;
    while($test == 0) {
	$rc=system("bsub -o $dir/queue_$ijob.log -J GG2VV$ijob -q $farm < $dir/$dir.sh");
	# $rc=system("cd $dir ; bsub -q $farm < job$ijob/subjob");
	if ($rc == 0) { $test=1; }
	else {
	    print "ERROR in submitting job: $rc  .. retrying in 10s ...\n";
	    sleep(5);
	}
    }
    sleep(1);
    $ijob++;
}

exit 0;

###############################################################################

sub make_subjob {

$dir=@_[0];
$ijob=@_[1];

# for parallel
$subjob="#!/bin/bash 
export PATH=\$" . "PATH:.
export PATH=\$" . "PATH:${basedir}/lhapdf-5.9.1-out/bin
export LD_LIBRARY_PATH=${basedir}/lhapdf-5.9.1-out/lib:${basedir}/omniORB/lib

export LHAPATH=${basedir}/lhapdf-5.9.1-out/share/lhapdf

export TOPDIR=${basedir}/gg2VV-${gg2VVvers} 
export PATH=\$" . "PATH:${basedir}/omniORB/bin 

source /afs/cern.ch/sw/IntelSoftware/linux/all-setup.sh

cd ${basedir}/gg2VV-${gg2VVvers}/gg2VV/${jobName}_$ijob
./gg2VV > ${jobName}_$ijob.log
";

if ($onlyHist == 0) {
$subjob = $subjob . "/bin/mv gg2VV_unweightedEvents.dat ${jobName}_$ijob.lhe
/afs/cern.ch/cms/caf/scripts/cmsStage -f ${jobName}_$ijob.lhe ${castorarea}/${jobName}
/bin/rm -f ${jobName}_$ijob.lhe
";
}

open(FILE,">$dir/$dir.sh");
print FILE "$subjob";
close(FILE);
system("chmod u+x $dir/$dir.sh");

}

