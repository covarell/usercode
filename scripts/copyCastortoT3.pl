#!/usr/bin/perl

# Copy to T3 a list of castor files

$workdir = ".";
$thissource = $ARGV[0];
$thisdest = $ARGV[1];
$thispart = $ARGV[2];
   
system("rfdir ${thissource} >> tmpfile2");
open(INFILE2,"tmpfile2") or die "cannot open tmpfile2";;
@log2=<INFILE2>;
close(INFILE2);
              
system("rfdir /castor/cern.ch/user/c/covarell/${thisdest} >> tmpfile3");
open(INFILE3,"tmpfile3") or die "cannot open tmpfile3";;
@log3=<INFILE3>;
close(INFILE3);

my $ok = 0;

system("setenv STAGE_HOST castorcms");
system("setenv STAGE_SVCCLASS cmst3");

# find missing files
foreach $line2 (@log2) {
  $ok = 0;
  if ($line2 =~  m/(${thispart})/) {
      my $partline2 = substr($line2,67);
      my $rightmarker2 = rindex($partline2,'t') - length($partline2) + 1;
      my $filename2 = substr($partline2,0,$rightmarker2);
      chomp($filename2);      
   
      foreach $line3 (@log3) {
	  if ($line3 =~  m/(${filename2})/) {$ok = 1;};
      }
      # print("$ok \n");
      if ($ok == 0) {
	  print("rfcp ${thissource}/${filename2} /castor/cern.ch/user/c/covarell/${thisdest} \n");
	  system("rfcp ${thissource}/${filename2} /castor/cern.ch/user/c/covarell/${thisdest}");
          ## UNCOMMENT TO DELETE ORIGINAL FILE (NOT RECOMMENDED...)
	  # print("rfrm /castor/cern.ch/user/c/covarell/${thissource}/${filename2} \n");
	  # system("rfrm /castor/cern.ch/user/c/covarell/${thissource}/${filename2}");
      }
  }
};
system("rm -f tmpfile2 tmpfile3");  

print "Done... \n";
 
