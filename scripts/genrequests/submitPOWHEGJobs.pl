#!/usr/bin/perl
###############################################################################
# Launch POWHEG jobs
###############################################################################

# job steering ----------------------------------------------------------------

$pwgIn = $ARGV[0];
# $pwgProcess = $ARGV[1];
$jobName = $ARGV[1];
$nJobs = $ARGV[2];
$hMass = $ARGV[3];
$theSeed = $ARGV[4];
$farm = $ARGV[5];

# name of job
$castorarea="/store/user/covarell/POWHEG-lhe";

# interactive or lxbatch queue
# $farm="I";
# $farm="2nw";
# $farm="1nw";

# site-specific paths etc

#$basedir="/afs/cern.ch/work/c/covarell/powheg/POWHEG-BOX/${pwgProcess}";
system("pwd > pwdfile");
open(INFILE2,"pwdfile") or die "cannot open pwdfile!";;
@log2=<INFILE2>;
close(INFILE2);
foreach $line (@log2) {
    chomp($line);
    $basedir = $line;
}
system("rm -f pwdfile");

# -----------------------------------------------------------------------------

print "\n";
print "-------------------------------------------------------------------------------\n";
print "L A U N C H I N G   P O W H E G   J O B S\n";
print "-------------------------------------------------------------------------------\n";
print "Workdir: ${basedir}\n";
print "EOS dir: ${castorarea}/${jobName}\n";
print "powheg.input: ${pwgIn}\n";
print "Higgs Mass: ${hMass}\n";

system("cd ~/work/mela/CMSSW_5_2_5/src");
system("cmsenv");
system("cd ${basedir}");
system("cmsMkdir ${castorarea}/${jobName}");

# -----------------------------------------------------------------------------
# lookup Higgs width
open(INFILE,"../widthhiggs.txt") or die "cannot open widthhiggs.txt!";;
@log=<INFILE>;
close(INFILE);

$hWidth = 0.000270375;
foreach $line (@log) {
    @splitline = split(/ +/, $line);
    if ($splitline[0] == ${hMass}) {
	$hWidth = $splitline[3];	
    }
}
chomp($hWidth);
print "Higgs Width: ${hWidth}\n";

$hFact = $hMass/1.2;

# -----------------------------------------------------------------------------
# create N jobs

$ijob=1;
while ( $ijob <= $nJobs ) {

    system("mkdir ${jobName}_$ijob");
    # print "sed s/THMASS/${hMass}/g ${pwgIn} > ${basedir}/${jobName}_$ijob/powheg.input";
    system("sed s/THMASS/${hMass}/g ${pwgIn} > ${basedir}/${jobName}_$ijob/powheg.input");
    # print "sed s/THWIDTH/\'${hWidth}\'/g ${basedir}/${jobName}_$ijob/powheg.input > tmpfile";
    system("sed s/THWIDTH/\'${hWidth}\'/g ${basedir}/${jobName}_$ijob/powheg.input > tmpfile");
    system("mv tmpfile ${basedir}/${jobName}_$ijob/powheg.input");
    # print "sed s/THHFACT/\'${hFact}\'/g ${basedir}/${jobName}_$ijob/powheg.input > tmpfile";
    system("sed s/THHFACT/\'${hFact}\'/g ${basedir}/${jobName}_$ijob/powheg.input > tmpfile");
    system("mv tmpfile ${basedir}/${jobName}_$ijob/powheg.input");
    my $seed = $theSeed + 876*$ijob;
    system("sed s/THSEED/${seed}/g ${basedir}/${jobName}_$ijob/powheg.input > tmpfile");
    system("mv tmpfile ${basedir}/${jobName}_$ijob/powheg.input");

    $dir="/tmp/covarell/";
    make_subjob($dir,$ijob);

    print "Submit job $ijob ... ";
    $test=0;
    while($test == 0) {
	$rc=system("bsub -o /tmp/covarell/stafava_$ijob.log -J POWHEG$ijob -q $farm < $dir/subjob_$ijob");
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
$subjob="#!/bin/zsh 
cd ~/work/mela/CMSSW_5_2_5/src
eval \`scramv1 runtime -sh\`
export LHAPATH=/afs/cern.ch/work/c/covarell/powheg/lhapdf-5.8.5-out/share/lhapdf
cd ${basedir}/${jobName}_$ijob
../pwhg_main > ${basedir}/${jobName}_$ijob/${jobName}_$ijob.log
head -n -1 ${basedir}/${jobName}_$ijob/pwgevents.lhe > ${basedir}/${jobName}_$ijob/pwgevents.new.lhe
mv ${basedir}/${jobName}_$ijob/pwgevents.new.lhe ${basedir}/${jobName}_$ijob/pwgevents.lhe
cmsStage -f ${basedir}/${jobName}_$ijob/pwgevents.lhe ${castorarea}/${jobName}/pwgevents_$ijob.lhe
rm -f ${basedir}/${jobName}_$ijob/pwgevents.lhe
";

open(FILE,">$dir/subjob_$ijob");
print FILE "$subjob";
close(FILE);
system("chmod u+x $dir/subjob_$ijob");

}

