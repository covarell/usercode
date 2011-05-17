#!/usr/bin/perl

# 'hadd' a lot of files 

$workdir = ".";
$thisdata = $ARGV[0];
$thispart = $ARGV[1];
$resultfilename = $ARGV[2];
                 
# retrieve cfi file
system("cmsenv");
system("rfdir ${thisdata} | grep ${thispart} >> tmpfile");
open(INFILE,"tmpfile") or die "cannot open tmpfile";;
@log=<INFILE>;
close(INFILE);
my $totalstring = '';

foreach $line (@log) {
   my $partline = substr($line,67);
   my $rightmarker = rindex($partline,'t') - length($partline) + 1;
   my $filename = substr($partline,0,$rightmarker);
   chomp($filename);
   $totalstring = $totalstring . " rfio:" . $thisdata . "/" . $filename;   
}

system("hadd ${resultfilename} ${totalstring}");
system("rm -f tmpfile");
print "Result in ${resultfilename}\n";
 
