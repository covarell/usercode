#!/usr/bin/perl

# inspect triggers in a list of runs

$startrun = $ARGV[0];
$endrun = $ARGV[1];
my $thisTrig = $ARGV[2];

system("rm -f trigfile.txt");

for ($count = $startrun; $count <= $endrun; $count = $count + 5) {
   system("echo ${count}");
   system("echo ${count} >> trigfile.txt");
   system("edmConfigFromDB --runNumber ${count} | hltDumpStream | grep \'\\<${thisTrig}\\>\' | head -n1 >> trigfile.txt");
}

# send staging requests

print "Results in trigfile.txt. \n";
 
