#!/usr/bin/perl

# Prestage all files in a CASTOR dir or part of them

$workdir = ".";
$thisdata = $ARGV[0];
my $thispart = $ARGV[1];

# retrieve cfi file
system("rfdir ${thisdata} >> tmpfile");
open(INFILE,"tmpfile") or die "cannot open tmpfile";;
@log=<INFILE>;
close(INFILE);

# send staging requests
foreach $line (@log) {
    if ($line =~ m/(${thispart})/) {   
        my $partline = substr($line,67);
	my $rightmarker = rindex($partline,'t') - length($partline) + 1;
        my $filename = substr($partline,0,$rightmarker);
        $filename = $thisdata . "/" . $filename;
        
	print "stager_get -M $filename \n";
        system("stager_get -M $filename");
    }
}

system("rm -f tmpfile");
print "Files staging... Be VERY patient... \n";
 
