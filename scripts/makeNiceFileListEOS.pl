#!/usr/bin/perl

# make a nice list of files from EOS

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
system("cmsLs ${thisdata} | grep ${thispart} >> tmpfile");
open(INFILE,"tmpfile") or die "cannot open tmpfile";;
@log=<INFILE>;
close(INFILE);

# lines of cfi file
system("cmsLs ${thisdata} | grep ${thispart} | wc -l >> tmpfile2");
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
   my $filename = substr($line,43);
   #my $rightmarker = rindex($partline,'t') - length($partline) + 1;
   #my $filename = substr($partline,0,$rightmarker);
   chomp($filename);  
   print "${filename} \n";
   system("edmFileUtil -d ${filename} >> utilfile");
   open(UTFILE,"utilfile") or die "cannot open utilfile";;
   @util=<UTFILE>;
   close(UTFILE);
   foreach $utline (@util) {
       chomp($utline);
       # print "utline $utline \n";
       $correctpath = $utline;
   }
   if ($forpython == 2) {
       if ($ilines == $nlines) {
	   $correctpath = "\'" . $correctpath . "\'";
       } else {
	   $correctpath = "\'" . $correctpath . "\',"; 
       }
   }
   print MYFILE "$correctpath \n"; 
   if ($forpython == 2) { 
       if ($ilines % 250 == 249) {  # every 250 lines
          print MYFILE "]) \n";
	  print MYFILE "readFiles.extend([ \n";
      }
   }

   system("rm -f utilfile");
}

system("rm -f tmpfile tmpfile2");
print "Result in ${textfilename}\n";
 
