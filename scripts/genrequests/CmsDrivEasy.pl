#!/usr/bin/perl

# Call CMSDriver in an intelligent way

$workdir = ".";
$thisnev = $ARGV[0];
$theCommandWithoutMinusNTen = "cmsDriver.py Configuration/GenProduction/python/EightTeV/POWHEG_PYTHIA6_Tauola_H_ZZ_2l2nu_8TeV_cff.py --filein lhe:8432   --mc --eventcontent RAWSIM            --datatier GEN-SIM    --conditions START53_V7C::All    --beamspot Realistic8TeVCollision     --step GEN,SIM --python_filename HIG-Summer12-01583_1_cfg.py --no_exec";

system("${theCommandWithoutMinusNTen} -n ${thisnev}");

$i = 0;
$sample = "";
$user = $ENV{'USER'};

@splitline = split(/ +/, $theCommandWithoutMinusNTen);
foreach $line (@splitline) {
    if ($line =~ /python_filename/) {
	$sample = $splitline[$i+1];
    } 
    $i++;
}
chomp($sample);

my $repl = "process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True)";
replace("${workdir}/${sample}",$repl);

my $repl2 = "    fileName = cms.untracked.string('/tmp/${user}/${sample}.root'), \n";
replacef("${workdir}/${sample}",$repl2);

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

sub replacef {

$infile = @_[0];
$repl = @_[1];

open(INFILE,"$infile") or die "cannot open $infile";;
@log=<INFILE>;
close(INFILE);

system("rm -f tmp");
open(OUTFILE,">tmp");

foreach $line (@log) {
  if ($line =~ /fileName/ && $line =~ /untracked.string/) { print OUTFILE $repl; }
  else { print OUTFILE $line; }
}

close(OUTFILE);
system("mv tmp $infile");

}
