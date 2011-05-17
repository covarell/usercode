#!/usr/bin/perl

# Copy to T3 a list of T2 files

$workdir = ".";
$thissource = $ARGV[0];
$thisdest = $ARGV[1];
$thispart = $ARGV[2];
                 
system("srmls ${thissource} >> tmpfile");
open(INFILE,"tmpfile") or die "cannot open tmpfile";;
@log=<INFILE>;
close(INFILE);

system("rfdir /castor/cern.ch/user/c/covarell/${thisdest} >> tmpfile3");
open(INFILE3,"tmpfile3") or die "cannot open tmpfile3";;
@log3=<INFILE3>;
close(INFILE3);

my $ok = 0;

# find missing files
foreach $line (@log) {
  $ok = 0;
  if ($line =~  m/(${thispart})/) {
      my $leftmarker = rindex($line,'/');
      my $filename = substr($line,$leftmarker+1,length($line));
      chomp($filename);
      # print("$filename \n");
      foreach $line3 (@log3) {
	  if ($line3 =~  m/(${filename})/) {$ok = 1;};
      }
      # print("$ok \n");
      if ($ok == 0) {
	  print("lcg-cp -b -D srmv2 -t 4800 --verbose --vo cms ${thissource}/${filename} srm://srm-cms.cern.ch:8443/srm/managerv2\?SFN=/castor/cern.ch/user/c/covarell/temp/${filename} \n");
	  system("lcg-cp -b -D srmv2 -t 4800 --verbose --vo cms ${thissource}/${filename} srm://srm-cms.cern.ch:8443/srm/managerv2\?SFN=/castor/cern.ch/user/c/covarell/temp/${filename}");
      };
  }
};

system("rm -f tmpfile");
system("rm -f tmpfile3");

system("T3Castor");

system("rfdir /castor/cern.ch/user/c/covarell/temp/ >> tmpfile2");
open(INFILE2,"tmpfile2") or die "cannot open tmpfile2";;
@log2=<INFILE2>;
close(INFILE2);

foreach $line2 (@log2) {
  if ($line2 =~  m/(${thispart})/) {
      my $partline2 = substr($line2,67);
      my $rightmarker2 = rindex($partline2,'t') - length($partline2) + 1;
      my $filename2 = substr($partline2,0,$rightmarker2);
        
      print("rfcp /castor/cern.ch/user/c/covarell/temp/${filename2} /castor/cern.ch/user/c/covarell/${thisdest} \n");
      system("rfcp /castor/cern.ch/user/c/covarell/temp/${filename2} /castor/cern.ch/user/c/covarell/${thisdest}");
      print("rfrm /castor/cern.ch/user/c/covarell/temp/${filename2} \n");
      system("rfrm /castor/cern.ch/user/c/covarell/temp/${filename2}");
      
  }
};
system("rm -f tmpfile2");  

print "Done... sigh... \n";
 
