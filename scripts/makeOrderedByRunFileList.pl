#!/usr/bin/perl

# make a nice list of files 

$workdir = ".";
$thisdata = $ARGV[0];
$textfilename = $ARGV[1];
$forpython = $ARGV[2];

if (!($forpython == 1 || $forpython == 2)) { 
  print("Argument n. 3 must be: 1 (a simple list) or 2 (a list for the PoolSource in a python cfg) \n");
  exit;
} 

if ($thisdata =~ /castor/) { 
  print("Argument n. 1 must be the directory WITHOUT /castor/cern.ch/cms/, i.e. starting with store/ \n");
  exit;
} 
                 
# retrieve cfi file
system("rfdir /castor/cern.ch/cms/${thisdata} >> tmpfile");
open(INFILE,"tmpfile") or die "cannot open tmpfile";;
@log=<INFILE>;
close(INFILE);

open(MYFILE,">${textfilename}");;

my $ilines = 0;
foreach $line (@log) {
   my $dirname = find_dir_in_line($line);
   # print "dirname $dirname \n";
   system("rfdir /castor/cern.ch/cms/${thisdata}/${dirname} >> tmpfile2");
   open(INFILE2,"tmpfile2") or die "cannot open tmpfile2";;
   @log2=<INFILE2>;
   close(INFILE2);
   foreach $line2 (@log2) {
       my $dirname2 = find_dir_in_line($line2);
       # print "dirname2 $dirname2 \n";
       system("rfdir /castor/cern.ch/cms/${thisdata}/${dirname}/${dirname2} >> tmpfile3");
       open(INFILE3,"tmpfile3") or die "cannot open tmpfile3";;
       @log3=<INFILE3>;
       close(INFILE3);
       foreach $line3 (@log3) {
	   my $dirname3 = find_dir_in_line($line3);
	   # print "dirname3 $dirname3 \n";
	   system("rfdir /castor/cern.ch/cms/${thisdata}/${dirname}/${dirname2}/${dirname3} >> tmpfile4");
	   open(INFILE4,"tmpfile4") or die "cannot open tmpfile4";;
	   @log4=<INFILE4>;
	   close(INFILE4);
	   system("edmFileUtil -d ${thisdata}/${dirname}/${dirname2}/${dirname3} >> utilfile");
	   open(UTFILE,"utilfile") or die "cannot open utilfile";;
	   @util=<UTFILE>;
	   close(UTFILE);
	   foreach $utline (@util) {
	       chomp($utline);
	       my $correctpath = $utline;
	   }
	   foreach $line4 (@log4) {
	       my $filename = find_file_in_line($line4);
	       # print "filename $filename \n";
	       if ($forpython == 1) {
		   $filename = $correctpath . "/" . $filename;
	       } else {
		   $filename = "\'" . $correctpath . "/" . $filename . "\',"; 
	       }
	       print MYFILE "$filename \n";    
	   }
	   system("rm -f tmpfile4");
       }
       system("rm -f tmpfile3 utilfile");
   }
   system("rm -f tmpfile2");
}

system("rm -f tmpfile");
print "Result in ${textfilename}\n";
 
## No fixed lenght, look for the t of root in the name
sub find_file_in_line
{
    my $line=shift;
    my $partline = substr($line,67);
    my $rightmarker = rindex($partline,'t') - length($partline) + 1;
    my $filename = substr($partline,0,$rightmarker);
    chomp($filename);
    return $filename;
}

## Fixed lenght, 3 chars
sub find_dir_in_line
{
    my $line=shift;
    my $partline = substr($line,67);
    my $dirname = substr($partline,0,4);
    chomp($dirname);
    return $dirname;
}
