#!/usr/bin/perl

# Clean a CASTOR dir or part of it

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
        
	print "rfrm $filename \n";
        system("rfrm $filename");
    }
}

system("rm -f tmpfile");
print "Directory clean! Bye bye. \n";
 
