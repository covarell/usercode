#!/usr/bin/perl

###############################################################################
# Compute filter efficiency from logs of automatic production
###############################################################################

$datasetIn = $ARGV[0];

print "$datasetIn \n";

# job steering ----------------------------------------------------------------
system("rfdir $datasetIn | grep LogCollect-1 > tmpfile");

$filename = "";
findfilename("tmpfile",67,$filename);
print "$filename \n";

system('rm -f tmpfile');
$mycommand = "rfcp $datasetIn/" . $filename . " ."; 
system($mycommand);
$mycommand = "tar -xvf " . $filename; 
system($mycommand);

system("ls WMTaskSpace/logCollect1/*.tar.gz > tmpfile2");

open(INFILE2,"tmpfile2") or die "cannot open tmpfile2";;
@log2=<INFILE2>;
close(INFILE2);

$run = 0;
$passed = 0;

foreach $line (@log2) {
   chomp($line);
   print "$line \n";
   $mycommand = "tar -zxvf " . $line; 
   system($mycommand);
   system('more cmsRun1/cmsRun1-stdout.log | grep Nev > tmpfile3');
   
   open(INFILE3,"tmpfile3") or die "cannot open tmpfile3";;
   @log3=<INFILE3>;
   close(INFILE3);
    
   foreach $line3 (@log3) {
       $thisrun = substr($line3,14,5);
       print "Run: $thisrun \n";
       $run += $thisrun;
   }
   system('rm -f tmpfile3');

   system('more cmsRun1/FrameworkJobReport.xml | grep TotalEvents > tmpfile4');
   
   open(INFILE4,"tmpfile4") or die "cannot open tmpfile4";;
   @log4=<INFILE4>;
   close(INFILE4);
    
   foreach $line4 (@log4) {
       $tempor = substr($line4,13);
       my $rightmarker = rindex($tempor,'<') - length($tempor);
       my $thispassed = substr($tempor,0,$rightmarker);
       print "Passed: $thispassed \n";
       $passed += $thispassed;
   }
   system('rm -f tmpfile4');
}

$ratio = $passed/$run;
$error = sqrt($ratio*(1-$ratio)/$run);
print "\nTotal run: $run \n";
print "Total passed: $passed \n";
print "Efficiency: ";
print &restrict_num_decimal_digits($ratio,7);
print " +/- ";
print &restrict_num_decimal_digits($error,7);
print " \n";

system('rm -f tmpfile2');

###############################################################################

sub findfilename {

$infile = @_[0];
$beginline = @_[1];
$filename = @_[2];

open(INFILE,"$infile") or die "cannot open $infile";;
@log=<INFILE>;
close(INFILE);

foreach $line (@log) {
    if (!($line =~ /Production/)) {
	chomp($line);
	$filename = substr($line,$beginline);
    }
}

}

sub restrict_num_decimal_digits  {

  my $num=shift;#the number to work on
  my $digs_to_cut=shift;# the number of digits after 
		        # the decimal point to cut 
		        #(eg: $digs_to_cut=3 will leave 
	         	# two digits after the decimal point)

  if ($num=~/\d+\.(\d){$digs_to_cut,}/)
  {
    # there are $digs_to_cut or 
    # more digits after the decimal point
    $num=sprintf("%.".($digs_to_cut-1)."f", $num);
  }
  return $num;
}
