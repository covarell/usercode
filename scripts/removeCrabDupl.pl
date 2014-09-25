#!/usr/bin/perl

# remove strange CRAB duplicates

$workdir = ".";
$thisdata = $ARGV[0];
$thispart = $ARGV[1];
                 
system("cmsLs ${thisdata} >> tmpfile");
open(INFILE,"tmpfile") or die "cannot open tmpfile";;
@log=<INFILE>;
close(INFILE);

$fileno = 0;

# find extra files
foreach $line (@log) {
  if ($line =~  m/(${thispart})/) {
      my $filenobefore = $fileno;
      my $partline = substr($line,43);
      my $rightmarker = rindex($partline,'t') - length($partline) + 1;
      my $filename = substr($partline,0,$rightmarker);
      chomp($filename);
      my $leftmarker2 = rindex($filename,'_',length($filename)-12);
      my $rightmarker2 = rindex($filename,'_',length($filename)-10) - length($filename);
      $fileno = substr($filename,$leftmarker2+1,$rightmarker2);
      print("File number = $fileno, Next... \n");
      # print("$ok \n");
      if ($fileno == $filenobefore) {
         print("A crab duplicate! \n");
	 print("cmsRm ${filename} \n");
	 system("cmsRm ${filename} \n");
      }
  };
};

system("rm -f tmpfile");  
print "Done \n";
 
