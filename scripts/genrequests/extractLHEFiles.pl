#!/usr/bin/perl

###############################################################################
# Extract many files from GEN
###############################################################################

$datasetIn = $ARGV[0];
$dir = $ARGV[1];

print "$datasetIn \n";

# job steering ----------------------------------------------------------------
system('cd ~/scratch0/mcrequests/CMSSW_5_3_8/src/');
system('cmsenv');
system('cd /tmp/covarell');

# print("dbs search --query=\'find file where dataset like ${datasetIn}\' > tmpfile \n");
system("dbs search --query=\'find file where dataset like ${datasetIn}\' > tmpfile");

open(INDATA,"tmpfile") or die "cannot open tmpfile";;
@log2=<INDATA>;
close(INDATA);

$ijob=0;
my $str = "";

foreach $line2 (@log2) {
    if ($line2 =~ /store/) {
	${ijob}++;
	$str = 'cp /afs/cern.ch/user/c/covarell/scratch0/mcrequests/CMSSW_5_3_8/src/GeneratorInterface/LHEInterface/test/testExternalLHEAsciiDumper_tpl.py ./launch_' . $ijob . '.py';
	system($str); 
	chomp($line2);
        $str = 'launch_' . $ijob . '.py';
	$line2 = "\'" . $line2 . "\'"; 
	replace($str,$line2);
	$str = 'cmsRun launch_' . $ijob .'.py' ; 
	system($str);
        $str = 'mv ascii_dump.lhe ' . $dir . '/lhefile_' . $ijob . '.lhe';
        system($str);
    }
}

# system('rm -f tmpfile');

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
  if ($line =~ /STOCAZZO/) { print OUTFILE $repl; }
  else { print OUTFILE $line; }
}

close(OUTFILE);
system("mv tmp $infile");

}
