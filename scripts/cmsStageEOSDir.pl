#!/usr/bin/perl

# CmsStage a EOS dir or part of it

$workdir = ".";
$thisdata = $ARGV[0];
my $thispart = $ARGV[1];
my $newName = $ARGV[2];
my $total = $ARGV[3];

# retrieve cfi file
system("cmsLs ${thisdata} >> tmpfile");
open(INFILE,"tmpfile") or die "cannot open tmpfile";;
@log=<INFILE>;
close(INFILE);

# send staging requests
foreach $line (@log) {
    if ($line =~ m/(${thispart})/) {
        my $filename = substr($line,43);
	# my $rightmarker = rindex($partline,'e') - length($partline) + 1;
        #my $filename = substr($partline,0,$rightmarker);
       #  $filename = $thisdata . "/" . $filename;
        chomp($filename);
	print "cmsStage -f $filename . \n";
        system("cmsStage -f $filename .");
    }
}

$jj = 0;

while ($jj <= $total) {
    $jj++;
    print "mv pwgevents_$jj.lhe ${newName}_$jj.lhe \n";
    system("mv pwgevents_$jj.lhe ${newName}_$jj.lhe");
}

system("rm -f tmpfile");
print "Directory clean! Bye bye. \n";
 
