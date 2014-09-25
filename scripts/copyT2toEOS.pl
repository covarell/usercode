#!/usr/bin/perl

# Copy to T3 a list of T2 files

$workdir = ".";
$thissource = $ARGV[0];
$thisdest = $ARGV[1];
$thispart = $ARGV[2];
                 
# system("srmls ${thissource} >> tmpfile");
open(INFILE,"tmpfile") or die "cannot open tmpfile";;
@log=<INFILE>;
close(INFILE);

system("cmsLs /store/user/covarell/${thisdest} >> tmpfile3");
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
	  print("lcg-cp -b -D srmv2 -t 4800 --verbose --vo cms ${thissource}/${filename} /tmp/covarell/${filename} \n");
	  system("lcg-cp -b -D srmv2 -t 4800 --verbose --vo cms ${thissource}/${filename} /tmp/covarell/${filename}");
          print("cmsStage /tmp/covarell/${filename} /store/user/covarell/${thisdest}"); 
          system("cmsStage /tmp/covarell/${filename} /store/user/covarell/${thisdest}"); 
	  print("rm -f /tmp/covarell/${filename}");
          system("rm -f /tmp/covarell/${filename}"); 
      };
  }
};

print "Done... sigh... \n";
 
