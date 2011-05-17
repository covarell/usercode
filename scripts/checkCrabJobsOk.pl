#!/usr/bin/perl

# Find missing successful jobs in a list of CASTOR files

$workdir = ".";
$thisdata = $ARGV[0];
$thispart = $ARGV[1];
$thisnumber = $ARGV[2];
$deepcheck = $ARGV[3]; 
if (!($deepcheck == 1 || $deepcheck == 2)) { 
  print("Argument n. 4 must be: 1 (check only missing files) or 2 (check also files with abnormal size) \n");
  exit;
} 
                 
my $thispart0 = "zh                  35";
my $thispart1 = "zh                  32";  # smart strings
my $thispart2 = "zh                  36";  # indicating successful job
my $thispart3 = "zh                  34";  # and file in CASTOR ok
my $thispart4 = "zh                  25";  # (see example)
my $thispart5 = "zh                  24";

# retrieve cfi file
system("rfdir ${thisdata} | grep '${thispart}' >> tmpfile");
# print("rfdir ${thisdata} | grep ${thispart} >> tmpfile");
open(INFILE,"tmpfile") or die "cannot open tmpfile";;
@log=<INFILE>;
close(INFILE);

open(MYFILE,">toresubmit.txt");;

my $result = "";
my $ok = 0;

# find missing files
for($i = 1; $i <= $thisnumber; $i++) {
  system("grep '_$i.root' tmpfile > tmpfile2");
#  print("grep '_$i.root' tmpfile > tmpfile2 \n");
  open(INFILE2,"tmpfile2") or die "cannot open tmpfile2";;
  @log2=<INFILE2>;
  close(INFILE2);
  $ok = 0;
  foreach $line2 (@log2) {
    if ($line2 =~ "root" ) {$ok = 1;};
  }
  if ($ok == 0)  {$result = $result . $i . ",";}; 
  system("rm -f tmpfile2");
}

# find other jobs to resubmit
if  ($deepcheck == 2) {
foreach $line (@log) {
    if (!($line =~ m/(${thispart0})/ || $line =~ m/(${thispart1})/ || $line =~ m/(${thispart2})/ || $line =~ m/(${thispart3})/ || $line =~ m/(${thispart4})/ || $line =~ m/(${thispart5})/)) {   
        my $leftmarker = rindex($line,'_');
	my $rightmarker = rindex($line,'.');
        my $number = substr($line,$leftmarker+1,$rightmarker-$leftmarker-1);
        $result = $result . $number . ",";

    }
}
}

my $realresult = substr($result,0,length($result)-1);
print MYFILE "$realresult \n";
# system("rm -f tmpfile");
print "Job to be resubmitted in toresubmit.txt \n";
 
