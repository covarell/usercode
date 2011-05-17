#!/usr/bin/perl

# Copy a CASTOR dir or part of it

$workdir = ".";
$thisdata = $ARGV[0];
$newdir = $ARGV[1];
my $thispart = $ARGV[2];

my $beginline = 0; 
system("rfdir ${thisdata} >> tmpfile");
$beginline = 67;

open(INFILE,"tmpfile") or die "cannot open tmpfile";;
@log=<INFILE>;
close(INFILE);

# copy files
foreach $line (@log) {
    if ($line =~ m/(${thispart})/) {   
        my $partline = substr($line,$beginline);
	my $rightmarker = rindex($partline,'t') - length($partline) + 1;
        my $filename = substr($partline,0,$rightmarker);
        my $newrightmarker = rindex($filename,'.') - length($filename) - 6;
        my $newfilename = substr($filename,0,$newrightmarker);
        $newfilename = $newfilename . ".root";
        
	print "rfcp ${thisdata}/${filename} ${newdir}/${newfilename} \n";
	system("rfcp ${thisdata}/${filename} ${newdir}/${newfilename}");
    }
}

system("rm -f tmpfile");
print "Directory copied! Bye bye. \n";
 
