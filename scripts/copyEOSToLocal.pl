#!/usr/bin/perl

# make copy many EOS files

$workdir = ".";
$thisdata = $ARGV[0];
$thispart = $ARGV[1];
                 
# retrieve cfi file
system("cmsLs ${thisdata} | grep ${thispart} >> tmpfile");
open(INFILE,"tmpfile") or die "cannot open tmpfile";;
@log=<INFILE>;
close(INFILE);

foreach $line (@log) {
   my $filename = substr($line,43);
   #my $rightmarker = rindex($partline,'t') - length($partline) + 1;
   #my $filename = substr($partline,0,$rightmarker);
   chomp($filename);  
   print "cmsStage -f ${filename} . \n";
   system("cmsStage -f ${filename} . \n");  
}

system("rm -f tmpfile tmpfile2");
 
