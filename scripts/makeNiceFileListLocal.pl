#!/usr/bin/perl

# make a nice list of files 

$workdir = ".";
$thisdata = $ARGV[0];
$thispart = $ARGV[1];
$textfilename = $ARGV[2];
$forpython = $ARGV[3];

if (!($forpython == 1 || $forpython == 2)) { 
  print("Argument n. 3 must be: 1 (a simple list) or 2 (a list for the PoolSource in a python cfg) \n");
  exit;
} 
                 
# retrieve cfi file
system("ls -latr ${thisdata} | grep ${thispart} >> tmpfile");
open(INFILE,"tmpfile") or die "cannot open tmpfile";;
@log=<INFILE>;
close(INFILE);

# lines of cfi file
system("ls -latr ${thisdata} | grep ${thispart} | wc -l >> tmpfile2");
open(INFILE2,"tmpfile2") or die "cannot open tmpfile2";;
@log2=<INFILE2>;
close(INFILE2);

open(MYFILE,">${textfilename}");;

my $nlines = 0;
foreach $line2 (@log2) {
    $nlines = $line2;
}

my $ilines = 0;
foreach $line (@log) {
   $ilines = $ilines + 1;
   my $partline = substr($line,52);
   my $rightmarker = rindex($partline,'t') - length($partline) + 1;
   my $filename = substr($partline,0,$rightmarker);
   chomp($filename);
   if ($forpython == 1) {
       $filename = "file:" . $thisdata . "/" . $filename;
   } else {
       if ($ilines == $nlines) {
	   $filename = "\'file:" . $thisdata . "/" . $filename . "\'";
       } else {
	   $filename = "\'file:" . $thisdata . "/" . $filename . "\',"; 
       }
   }
   print MYFILE "$filename \n";   
}

system("rm -f tmpfile tmpfile2");
print "Result in ${textfilename}\n";
 
