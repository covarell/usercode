#!/usr/bin/perl
###############################################################################
# Launch a lot of jobs
###############################################################################

# determine location
$host=`echo \$HOST`;
if    ( $host =~ /fpslife/ ) { $location="fpslife"; }
elsif ( $host =~ /lxplus/ )  { $location="lxplus"; }
elsif ( $host =~ /cnaf/ )    { $location="cnaf"; }
elsif ( $host =~ /lxcms/ )   { $location="lxcmsg1"; }

# job steering ----------------------------------------------------------------

# name of job
$castorarea="/castor/cern.ch/user/c/covarell/egamma-trig/NEWRECO";

# cfg file
$steering="MinBias_210_RECO_all_cfg.py";

$datalist="MinBias_210_raw.list";
# number of files per job
$nfilesperjob=1;

$outfilebase="MinBias_217_RECO_startup";

# interactive or lxbatch queue
# $farm="I";
# $farm="8nm"; $resource="";
$farm="8nh";

$cmsswvers="CMSSW_2_1_7";
$scramarch="slc4_ia32_gcc345";
$scram=scramv1;

# site-specific paths etc

if ( $location eq "lxcmsg1" ) {
  $homedir="/afs/cern.ch/user/c/covarell/scratch0/egamma-trig";
  $basedir="${homedir}/${cmsswvers}";
  $outdir="/tmp";
}
elsif ( $location eq "lxplus" ) {
  # die "ERROR: Root files for the analysis are on lxcmsg1!\n";
  $homedir="/afs/cern.ch/user/c/covarell/scratch0/egamma-trig";
  $basedir="${homedir}/${cmsswvers}";
  $outdir="/tmp";
}
else {
  die "ERROR: location: $location unknown!\n";
}

# -----------------------------------------------------------------------------

$workdir="${basedir}/src/HLTriggerOffline/Egamma/python";
$datalistdir="${basedir}/src/HLTriggerOffline/Egamma/data/" . ${datalist};

print "\n";
print "-------------------------------------------------------------------------------\n";
print "L A U N C H I N G   J O B S\n";
print "-------------------------------------------------------------------------------\n";
print "Location: $location\n";
print "Workdir: ${workdir}\n";
print "Outdir: ${outdir}\n";
print "CASTOR dir: ${castorarea}\n";
print "Data list: ${datalistdir}\n";

# Calculate njobs (= nfiles / nfilesperjob, last files not taken)
system("more ${datalistdir} | grep root | wc -l >> nofiles.txt;");
open(INFILE,"nofiles.txt") or die "cannot open nofiles.txt";;
@nofiles=<INFILE>;
close(INFILE);
foreach $line1 (@nofiles) {$njobs = int($line1 / ${nfilesperjob});}
system("rm -f nofiles.txt;");

if ($njobs>1) { 
  print "Parallel Jobs: ${njobs} with ${nfilesperjob} file(s) per job\n";
}
else {
  # only a few files 
  print "Ma vaffanculo, vai!\n";
  exit;
}

# -----------------------------------------------------------------------------
# create N jobs

$ijob=1;
$ifile=0;

open(INFILE2,"$datalistdir") or die "cannot open $datalistdir";;
@listdata=<INFILE2>;
close(INFILE2);

foreach $line (@listdata) { 
    if ($line =~ /root/) {
	${ifile}++;
        if (${nfilesperjob} == 1) {   # 1 FILE PER JOB
            # make job dir and copy aux files
	    $dir="${workdir}/tmpscripts";
	    make_subjob($dir,$ijob);
	    
	    # create cfg file
	    system("cp ${workdir}/$steering ${workdir}/tmpcfg/cfgfile_$ijob.py");
	    $repl= "process.source = cms.Source(\"PoolSource\",
      fileNames = cms.untracked.vstring(
        $line
        )
      )
   
process.out = cms.OutputModule(\"PoolOutputModule\",
   fileName = cms.untracked.string('file:${outdir}/${outfilebase}_$ijob.root'),
                  ";
	    
	    replace("${workdir}/tmpcfg/cfgfile_$ijob.py",$repl);
	    ${ijob}++;
	} else {   # 2 OR MORE FILES PER JOB
	    if ($ifile % ${nfilesperjob} == 0) { 
		
		# make job dir and copy aux files
		$dir="${workdir}/tmpscripts";
		make_subjob($dir,$ijob);
		
		# create cfg file
		system("cp ${workdir}/$steering ${workdir}/tmpcfg/cfgfile_$ijob.py");
		
		$repl= $repl . "
       $line
       )
    )
   
process.out = cms.OutputModule(\"PoolOutputModule\",
   fileName = cms.untracked.string('file:${outdir}/${outfilebase}_$ijob.root'),
                                         
                  ";
	
		replace("${workdir}/tmpcfg/cfgfile_$ijob.py",$repl);
		${ijob}++;
	    
	    } elsif ($ifile % ${nfilesperjob} == 1) {
		$repl= "process.source = cms.Source(\"PoolSource\",
    fileNames = cms.untracked.vstring(
      $line,
                   ";	
	    } else {
		$repl= $repl . "$line,";
	    }
	}
    }
}

# system("cd ${workdir}") ;
$iret=&submit_jobs;
if ($iret ne 0) { die "ERROR in submit jobs!\n";}

exit 0;

###############################################################################

sub submit_jobs {

$ijob=1;
while ( $ijob <= $njobs ) {

#   if ($ijob ==  40 || $ijob == 8 || $ijob == 9) {  
       print "Submit job $ijob ... ";
       $test=0;
       while($test == 0) {
	   $rc=system("bsub -o log/log_$ijob.log -q $farm < tmpscripts/subjob_$ijob");
	   # $rc=system("cd $dir ; bsub -q $farm < job$ijob/subjob");
	   if ($rc == 0) { $test=1; }
	   else {
	       print "ERROR in submitting job: $rc  .. retrying in 10s ...\n";
	       sleep(5);
	   }
       }
       sleep(1);

#    }
    $ijob++;
}

return 0;
}

###############################################################################

sub make_subjob {

$dir=@_[0];
$ijob=@_[1];

# for parallel
$subjob="#!/bin/zsh 
#BSUB -J \"RECO\" 
#BSUB -C 0
cd ${workdir}
eval \`$scram runtime -sh\`
rehash
EdmPluginRefresh
cmsRun tmpcfg/cfgfile_$ijob.py > /dev/null
rfcp ${outdir}/${outfilebase}_$ijob.root $castorarea
";


open(FILE,">$dir/subjob_$ijob");
print FILE "$subjob";
close(FILE);
system("chmod u+x $dir/subjob_$ijob");

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
