#!/usr/bin/perl

# Call CMSDriver in an intelligent way

$workdir = ".";
$sample = $ARGV[0];
$thisnev = $ARGV[1];
$thistyp = $ARGV[2];

if (!($thistyp == 1 || $thistyp == 2)) { 
    print("Argument n. 3 must be: 1 (GEN only) or 2 (GENSIM) \n");
    exit;
} 

if ($thistyp == 1) {
    system("cmsDriver.py Configuration/GenProduction/python/${sample}.py -s GEN --no_output --conditions START311_V2::All --datatier GEN-SIM-RAW  --eventcontent RAWSIM --customise=Configuration/GenProduction/customise_SilentMessageLogger.py --no_exec -n ${thisnev}");
} else {
    system("cmsDriver.py Configuration/GenProduction/python/${sample}.py -s GEN,SIM --conditions START311_V2::All --datatier GEN-SIM --eventcontent RAWSIM --beamspot Realistic7TeV2011Collision --customise=Configuration/GenProduction/customise_SilentMessageLogger.py --no_exec -n ${thisnev}");
} 

my $repl = "process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True)";
if ($thistyp == 1) {
    replace("${workdir}/${sample}_py_GEN.py",$repl);
} else {
    replace("${workdir}/${sample}_py_GEN_SIM.py",$repl);
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
  if ($line =~ /process/ && $line =~ /options/) { print OUTFILE $repl; }
  else { print OUTFILE $line; }
}

close(OUTFILE);
system("mv tmp $infile");

}
