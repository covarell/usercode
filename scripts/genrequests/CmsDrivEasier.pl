#!/usr/bin/perl

# Call CMSDriver in an even more intelligent way

$workdir = ".";
$request = $ARGV[0];
$thisnev = $ARGV[1];
$gensim = $ARGV[2];

if (!($gensim == 1 || $gensim == 2)) { 
  print("Argument n. 3 must be: 1 (GEN step) or 2 (GEN-SIM step) \n");
  exit;
} 

system("wget --no-check-certificate https://cms-pdmv.cern.ch/mcm/public/restapi/requests/get_setup/${request}");

open(REQUEST,"${request}") or die "cannot open $request";;
@logr=<REQUEST>;
close(REQUEST);

$ir = 0;
$theCommand = "";

foreach $liner (@logr) {
  if ($liner =~ /curl/) { 
      chomp($liner); 
      $liner =~ s/github.com/githubusercontent.com/g; 
      print "\n DOING NOW ${liner} \n";
      system("${liner}");      
  }

  if ($liner =~ /scram b/) { 
      chomp($liner); 
      print "\n DOING NOW ${liner} \n";
      system("${liner}");
  }
  if ($liner =~ /cmsDriver/ && $ir == 0) { 
      $theCommand = $liner; 
      $ir++;
  }
} 
system("rm -f ${request}");

my $rightmarker = rindex($theCommand,'n') - length($theCommand) - 2;
my $theCommandWithoutMinusNTen = substr($theCommand,0,$rightmarker);

# $theCommandWithoutMinusNTen =~ s/\"/\\\"/g;
if ($gensim == 1) {$theCommandWithoutMinusNTen =~ s/GEN,SIM/GEN/g;}

print "\n DOING NOW $theCommandWithoutMinusNTen -n ${thisnev} \n";

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
