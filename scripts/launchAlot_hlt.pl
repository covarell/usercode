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
$castorarea="/castor/cern.ch/user/c/covarell/egamma-trig/HLT";

# cfg file
$steering="RelVal_HLT2_8E29_noinput.py";

$datalist="WenuRelVal_340.list";
# $datalist="GammaJet_220_few.list";
# $datalist="Gamma1000_220_few.list";
# $datalist="MinBiasNoFilt_210_few.list";

# number of files per job
$nfilesperjob=1;
# $nfilesperjob=20;

$outfilebase="Wenu340_withNewAPE_HLT8E29";

# interactive or lxbatch queue
# $farm="I";
$farm="8nm";
# $farm="1nh";

$cmsswvers="CMSSW_3_4_0";
$scramarch="slc5_ia32_gcc434";
$scram=scramv1;

# site-specific paths etc

if ( $location eq "lxcmsg1" ) {
  $homedir="/afs/cern.ch/user/c/covarell/scratch0/tagprobe";
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

$workdir="${basedir}/src/HLTrigger/Configuration/test";
$datalistdir="${basedir}/src/HLTrigger/Configuration/data/" . ${datalist};

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

print "Parallel Jobs: ${njobs} with ${nfilesperjob} file(s) per job\n";

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
   
process.output = cms.OutputModule(\"PoolOutputModule\",
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
   
process.output = cms.OutputModule(\"PoolOutputModule\",
   fileName = cms.untracked.string('file:${outdir}/${outfilebase}_$ijob.root'),
                                         
                  ";
	
		replace("${workdir}/tmpcfg/cfgfile_$ijob.py",$repl);
		${ijob}++;
	    
	    } elsif ($ifile % ${nfilesperjob} == 1) {
		$repl= "process.source = cms.Source(\"PoolSource\",
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
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
	   $rc=system("bsub -o log/${outfilebase}_$ijob.log -q $farm < tmpscripts/subjob_$ijob");
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
#BSUB -J \"HLT\" 
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
