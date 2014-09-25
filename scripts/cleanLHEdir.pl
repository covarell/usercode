#!/usr/bin/perl

# Clean a EOS dir or part of it

$workdir = ".";
$thisdata = $ARGV[0];
my $thispart = $ARGV[1];

# retrieve cfi file
system("cmsLs ${thisdata} >> tmpfile");
open(INFILE,"tmpfile") or die "cannot open tmpfile";;
@log=<INFILE>;
close(INFILE);

# send staging requests
foreach $line (@log) {
    if ($line =~ m/(${thispart})/) {   
        my $partline = substr($line,43);
	my $rightmarker = rindex($partline,'e') - length($partline) + 1;
        my $filename = substr($partline,0,$rightmarker);
       #  $filename = $thisdata . "/" . $filename;
        
	print "cmsRm $filename \n";
        system("cmsRm $filename");
    }
}

system("rm -f tmpfile");
print "Directory clean! Bye bye. \n";
 
