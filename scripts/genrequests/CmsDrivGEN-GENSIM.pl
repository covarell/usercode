#!/usr/bin/perl

# Call CMSDriver in an intelligent way

$workdir = ".";
$sample = $ARGV[0];
$thisnev = $ARGV[1];
$thistyp = $ARGV[2];
$thiscamp = $ARGV[3];

if (!($thistyp == 1 || $thistyp == 2)) { 
    print("Argument n. 3 must be: 1 (GEN only) or 2 (GENSIM) \n");
    exit;
} 

if (!($thiscamp == 1 || $thiscamp == 2 || $thiscamp == 3)) { 
    print("Argument n. 4 must be: 1 (Summer11-7TeV) or 2 (Summer12-8TeV) or 3 (Upgradde-14TeV) \n");
    exit;
} 

if ($thistyp == 1) {
    if ($thiscamp == 1) {
	system("cmsDriver.py Configuration/GenProduction/python/${sample}.py -s GEN --no_output --conditions START311_V2::All --datatier GEN-SIM-RAW  --eventcontent RAWSIM --customise=Configuration/GenProduction/customise_SilentMessageLogger.py --no_exec -n ${thisnev}");
    } 
    if ($thiscamp == 2) {
	system("cmsDriver.py Configuration/GenProduction/python/EightTeV/${sample}.py -s GEN --no_output --conditions START50_V13::All --datatier GEN-SIM-RAW  --eventcontent RAWSIM --customise=Configuration/GenProduction/customise_SilentMessageLogger.py --no_exec -n ${thisnev}");
    } 
    if ($thiscamp == 3) {
	system("cmsDriver.py Configuration/GenProduction/python/FourteenTeV/${sample}.py -s GEN --no_output --step GEN,SIM --geometry Extended2017 --beamspot Gauss --conditions STAR17_61_V1A::All --customise SLHCUpgradeSimulations/Configuration/combinedCustoms.cust_2017+fixRPCConditions --eventcontent FEVTDEBUG --datatier GEN-SIM-RAW --no_exec -n ${thisnev}");
    } 
} else {
    if ($thiscamp == 1) {
	system("cmsDriver.py Configuration/GenProduction/python/${sample}.py -s GEN,SIM --conditions START311_V2::All --datatier GEN-SIM --eventcontent RAWSIM --beamspot Realistic7TeV2011Collision --pileup NoPileUp --datamix NODATAMIXER --customise=Configuration/GenProduction/customise_SilentMessageLogger.py --no_exec -n ${thisnev}");
    } 
    if ($thiscamp == 2) {
	system("cmsDriver.py Configuration/GenProduction/python/EightTeV/${sample}.py -s GEN,SIM --conditions START50_V13::All --datatier GEN-SIM --eventcontent RAWSIM --beamspot Realistic8TeVCollision  --pileup NoPileUp --datamix NODATAMIXER --customise=Configuration/GenProduction/customise_SilentMessageLogger.py --no_exec -n ${thisnev}");
    } 
    if ($thiscamp == 3) {
	system("cmsDriver.py Configuration/GenProduction/python/FourteenTeV/${sample}.py -s GEN,SIM --geometry Extended2017 --beamspot Gauss --conditions STAR17_61_V1A::All --pileup NoPileUp --datamix NODATAMIXER --customise SLHCUpgradeSimulations/Configuration/combinedCustoms.cust_2017+fixRPCConditions --eventcontent FEVTDEBUG --datatier GEN-SIM --no_exec -n ${thisnev}");
    }     
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
