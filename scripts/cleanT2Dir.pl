#!/usr/bin/perl

# Clean a CASTOR dir or part of it

$workdir = ".";
$thisdata = $ARGV[0];
my $thispart = $ARGV[1];

# retrieve cfi file
system("unsetenv SRM_PATH");
system("srmls ${thisdata} >> tmpfile");
open(INFILE,"tmpfile") or die "cannot open tmpfile";;
@log=<INFILE>;
close(INFILE);

# send staging requests
foreach $line (@log) {
    if ($line =~ m/(${thispart})/) {
        my $leftmarker = rindex($line,'/');
        my $filename = substr($line,$leftmarker+1,length($line));
        chomp($filename);  
        $filename = $thisdata . "/" . $filename;
        
	print "srmrm $filename \n";
        system("srmrm $filename");
    }
}

system("rm -f tmpfile");
print "Directory clean! Bye bye. \n";
 
