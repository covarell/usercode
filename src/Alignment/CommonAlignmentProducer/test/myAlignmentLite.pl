#!/usr/bin/perl
###############################################################################
# Simple processing for alignment (from the main directory itself) 
###############################################################################

# determine location
$host=`echo \$HOST`;
if    ( $host =~ /fpslife/ ) { $location="fpslife"; }
elsif ( $host =~ /lxplus/ )  { $location="lxplus"; }
elsif ( $host =~ /cnaf/ )    { $location="cnaf"; }
elsif ( $host =~ /lxcms/ )  { $location="lxcmsg1"; }

# job steering ----------------------------------------------------------------

# name of job
$jobname="McTIBTOBstr-xy-phiRange2";

# db output files to be deleted except for last iteration
$dbfiles="Alignments.db condbcatalog.xml";

# db authentication file to be copied to exec dir
$authfile="authentication.xml";

# db input file
# $sqlitefile="HIPScenario.db"; $condbcatalogfile="condbcatalog.xml";
$sqlitefile=""; $condbcatalogfile="";

# number of events per job
$nevent=2000;

# first event
$firstev=0;

# number of jobs
$njobs=8;

#starting iteration
$firstiteration=19;

# number of total iterations (excluding initial step)
$iterations=20;

# interactive or lxbatch queue
# $farm="I";
# $farm="8nm"; $resource="";
$farm="dedicated -R cmsalca";

$cmsswvers="CMSSW_1_3_1";
$scramarch="slc3_ia32_gcc323";
$scram=scramv1;

# sleep time in seconds between two check cycles
$sleeptime=30;

# site-specific paths etc

if ( $location eq "lxcmsg1" ) {
  $homedir="/afs/cern.ch/user/c/covarell/scratch0/alignment-mc";
  $basedir="${homedir}/${cmsswvers}";
  $outdir="/data/covarell/joboutput";
}
elsif ( $location eq "lxplus" ) {
  # die "ERROR: Root files for the analysis are on lxcmsg1!\n";
   $homedir="/afs/cern.ch/user/c/covarell/scratch0/alignment-mc";
   $basedir="${homedir}/${cmsswvers}";
   $outdir="/afs/cern.ch/user/c/covarell/scratch0/joboutput";
}
elsif ( $location eq "fpslife" ) {
   die "ERROR: location: $location unsupported!\n";
#  $homedir="/afs/cern.ch/user/f/fpschill";
#  $basedir="${homedir}/cms/${cmsswvers}";
#  $outdir="/x01/usr/fpschill/cms/joboutput";
}
elsif ( $location eq "cnaf" ) {
   die "ERROR: location: $location unsupported!\n";
#  $homedir="/home/CMS/fpschill";
#  $basedir="${homedir}/${cmsswvers}";
#  $outdir="/afs/infn.it/cnaf/user/fpschill/joboutput";
}
else {
  die "ERROR: location: $location unknown!\n";
}

# -----------------------------------------------------------------------------

$workdir="${basedir}/src/Alignment/CommonAlignmentProducer/test";

print "\n";
print "-------------------------------------------------------------------------------\n";
print "A l i g n m e n t \n";
print "-------------------------------------------------------------------------------\n";
print "Location: $location\n";
print "Workdir: ${workdir}\n";
print "Outdir: ${outdir}\n";
if ($njobs>1) { 
  print "Parallel Jobs: ${njobs} with $nevent events per job\n";
}
else {
  print "Events: $nevent \n";
}
print "Iterations: ${iterations}\n";

# $odir=".";

# -----------------------------------------------------------------------------
# create sandbox and set up local environment

system("
cd $outdir/$jobname/$cmsswvers;
");
# rm -rf $outdir/$jobname;
# mkdir $outdir/$jobname;
# cd $outdir/$jobname/$cmsswvers;
# $scram project CMSSW $cmsswvers > /dev/null;
# cd $cmsswvers;
# mkdir -p lib/${scramarch};
# mkdir -p bin/${scramarch};
# mkdir -p config/${scramarch};
# cp $basedir/lib/${scramarch}/* lib/${scramarch} > /dev/null;
# cp $basedir/module/${scramarch}/* module/${scramarch} > /dev/null;
# ");

# deal with modified cf? files
# chdir "$basedir/src";
# @cfgfilelist = `ls */*/*/*.cf?`;
# foreach $file ( @cfgfilelist ) {
#  system("tar rf cfg.tar $file");
#  }
#  system("cp cfg.tar $odir/$cmsswvers/src; rm cfg.tar; cd $odir/$cmsswvers/src; tar xf cfg.tar");


# if one cpu create subjob and submit #########################################

if ($njobs == 1) {

  # make subjob
  $dir=".";
  make_subjob($dir,1);

  # create cfg file
  # system("cp ${workdir}/$steering $dir/cfgfile");
  # $repl="replace PoolSource.maxEvents  = $nevent ";
  # replace("$dir/cfgfile",$repl);

  if ($farm eq "I") {
    print "Run interactively ...\n";
    system("cd $dir; ./subjob");
  }
  else {
    die "ERROR: myAlignmentLite is only for interactive running.";
  }

}

# PARALLEL ####################################################################

else {

$ijob=1;
$ifirst=$firstev;

while ( ${ijob} <= ${njobs} ) {
  print "Job ${ijob}/${njobs} starting at event $ifirst \n";
  ${ijob}++;
  ${ifirst}=${ifirst}+${nevent};
};


$dir="$outdir/$jobname/$cmsswvers";
$iteration=$firstiteration;
system("rm -f $dir/iteration.txt");
open(FILE,">$dir/iteration.txt");
print FILE "$firstiteration";
close(FILE);

# enter loop ... ==============================================================

LOOP: {

  # check if all are finished
  if ($iteration gt $firstiteration ) { $alldone=&check_finished(); }

  # need to run collector / resubmit for next iteration
  if ( $alldone eq 1 or $iteration eq $firstiteration ) {

    # run collector
    print "Run collector for iteration $iteration ...\n";
    $iret=&run_collector();
    if ($iret ne 0) { die "ERROR in collector!\n";}

    # submit jobs for next iteration
    if ($iteration < $iterations) {
      $iteration++;

      system("rm -f $dir/iteration.txt");
      open(FILE,">$dir/iteration.txt");
      print FILE "$iteration";
      close(FILE);

      print "New iteration $iteration --------------------------------------------------------------\n";
      print "Submit all jobs for iteration $iteration ...\n";
      $iret=&submit_jobs;
      if ($iret ne 0) { die "ERROR in submit jobs!\n";}

      redo LOOP;
    }
    # we are done
    else {
      system("
        cd $dir
        rm -rf tmp src logs lib config bin
     ");
      print "Finished!\n";
    }

  }

  # else, just wait ...
  else {
    # print "sleeping for $sleeptime seconds ...\n";
    sleep($sleeptime);
    redo LOOP;
  }

}

# save some disk space
system("cd $dir; rm -rf job*"); 

}

exit 0;

###############################################################################

sub check_finished  
{

$alldone=1;
$ijob=1;

$timestamp=`date +%H:%M:%S_%d/%m/%y`;

print "Job status: ";

while ( $ijob <= $njobs ) {
  if ( -e "$dir/job$ijob/DONE" ) { print "1 "; }
  else { print "0 "; $alldone=0; }
  $ijob++;
}

print "at $timestamp";
if ($alldone eq 1) { return 1;}
else { return 0; }


}

###############################################################################

sub run_collector 
{
  system(" 
    cd $dir/main
    subjob
  ");
  if ($iteration gt 0) {
    system("cp $dir/job1/HIPAlignmentEvents.root $dir/main");
  }

  return 0;
}

###############################################################################

sub submit_jobs 
{

$ijob=1;
while ( $ijob <= $njobs ) {

  system("
    rm -f $dir/job$ijob/DONE

    cp $dir/main/IOTruePositions.root       $dir/job$ijob
    cp $dir/main/IOMisalignedPositions.root $dir/job$ijob
    cp $dir/main/IOAlignedPositions.root    $dir/job$ijob
    cp $dir/main/IOIteration.root           $dir/job$ijob

  ");

  if ( $farm eq "I" ) {
    print "Run job $ijob interactively ...\n";
    system("$dir/job$ijob/subjob &");
  }
  else {
    print "Submit job $ijob ... ";
    $test=0;
    while($test == 0) {
      $rc=system("cd $dir ; bsub -o /tmp/junk -q $farm < job$ijob/subjob");
      if ($rc == 0) { $test=1; }
      else {
	print "ERROR in submitting job: $rc  .. retrying in 10s ...\n";
	sleep(5);
      }
    }
  }
  sleep(1);
  $ijob++;
}

  return 0;
}

###############################################################################

sub make_subjob {

$dir=@_[0];
$single=@_[1];

# system("mkdir $dir");
# system("cp ${workdir}/$authfile $dir");
if ($sqlitefile !="") { system("cp ${workdir}/$sqlitefile $dir"); }
if ($condbcatalogfile !="") { system("cp ${workdir}/$condbcatalogfile $dir"); }

if ($single == 1) {

# for single job
system("touch $dir/timing.log");
$subjob="#!/bin/zsh -f
#BSUB -J \"ALIGN\"
#BSUB -C 0
cd $dir
eval \`$scram runtime -sh\`
rehash
thisiter=1
until ((thisiter > $iterations )); do
rm -f ${dbfiles}
echo \"Iteration \${thisiter} ...\"
date | read starttime
cmsRun cfgfile >& log.\${thisiter}.log
date | read endtime
gzip log.\${thisiter}.log
mv alignment.log alignment.\${thisiter}.log
gzip alignment.\${thisiter}.log
echo \"Iteration \${thisiter}: \${starttime} ... \${endtime}\" >> timing.log
let thisiter++
done
";

}

else {

# for parallel
$subjob="#!/bin/zsh 
#BSUB -J \"ALIGN\" 
#BSUB -C 0
cd $dir
eval \`$scram runtime -sh\`
rehash
cat ../iteration.txt | read iter
rm -f $dbfiles
cmsRun cfgfile >& log.\$iter.log
gzip log.\$iter.log
mv alignment.log alignment.\$iter.log
gzip alignment.\$iter.log
touch DONE
";

}

open(FILE,">$dir/subjob");
print FILE "$subjob";
close(FILE);
system("chmod u+x $dir/subjob");

}

###############################################################################

sub replace {

$infile = @_[0];
$repl = @_[1];

open(INFILE,"$infile") or die "cannot open $infile";;
@log=<INFILE>;
close(INFILE);

system("rm -f tmp");
open(OUTFILE,">tmp");

foreach $line (@log) {
  if ($line =~ /REPLACEME/) { print OUTFILE $repl; }
  else { print OUTFILE $line; }
}

close(OUTFILE);
system("mv tmp $infile");

}
