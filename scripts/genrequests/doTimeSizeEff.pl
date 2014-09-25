#!/usr/bin/perl

# Calculate filter efficiency and time per event

$workdir = ".";
$thislog = $ARGV[0];

open(INFILE,"$workdir/${thislog}.log") or die "cannot open ${thislog}.log!";;
@log=<INFILE>;
close(INFILE);
 
my $passedevents = 0; 
my $efffilter = 0;
my $errefffilter = 0;
my $effmatch = 0;
my $erreffmatch = 0;
my $timing = 0;
my $eventsize = 0;
my $xsec = 0;

$user = $ENV{'USER'};
$xsecfound = 0;

# Calculate the cross section (PYTHIA6)
foreach $line (@log) {
    if ($line =~ /subprocesses/) {  
	@splitline = split(/ +/, $line);
        # print "$splitline[10]\n";
        $integ = $splitline[10];
        @resplitline = split('D', $integ);
        $xsec = $resplitline[0] * (10 ** $resplitline[1]) * (10 ** 9) ;	
        $xsecfound = 1; 
    }
}

# Calculate the cross section (PYTHIA8)
if ($xsecfound == 0) { 
    foreach $line (@log) {
	if ($line =~ /sum /) {  
	    @splitline = split(/ +/, $line);
	    # print "$splitline[8]\n";
	    $integ = $splitline[8];
	    @resplitline = split('e', $integ);
	    $xsec = $resplitline[0] * (10 ** $resplitline[1]) * (10 ** 9) ;	
	}
    }
}

# Calculate the filter eff
foreach $line (@log) {
    if ($line =~ /TrigReport/ && $line =~ /Events total/) {  
	@splitline = split(/ +/, $line);
        # print "$splitline[7]\n";
        $passedevents = $splitline[7];	
	$efffilter = $splitline[7]/$splitline[4];
        $errefffilter = sqrt($efffilter*(1-$efffilter)/$splitline[4]);
    }
}

my $i=0;
# Calculate the matching eff
foreach $line (@log) {
    if ($line =~ /TrigReport/ && $line =~ /generator/ && $line =~ /  1  / && $i==0) {  
	@splitline = split(/ +/, $line);
  #      print "$splitline[3] $splitline[4]\n";	
	$effmatch = $splitline[4]/$splitline[3];
        $erreffmatch = sqrt($effmatch*(1-$effmatch)/$splitline[3]);
        $i++;
    }
}

# Divide!
# $xsec = $xsec/$effmatch;
$efffilter = $efffilter/$effmatch;
$errefffilter = $errefffilter/$effmatch;
if ($efffilter == 1) {$errefffilter = 0};

# Calculate the timing
foreach $line (@log) {
    if ($line =~ /TimeReport/ && $line =~ /CPU/ && $line =~ /event/) {  
	@splitline = split(/ +/, $line);	
	$timing = $splitline[3];
    }
}

# Calculate the event size
system("ls -l /tmp/${user}/${thislog}.py.root > dummyfile");

open(INFILE2,"$workdir/dummyfile") or die "cannot open dummyfile!";;
@log2=<INFILE2>;
close(INFILE2);
foreach $line (@log2) {  
    @splitline = split(/ +/, $line);	
    my $filesize = $splitline[4];
    $eventsize = $filesize/(1000*$passedevents); 
}
system("rm -f dummyfile");

print "Cross section = ";
print $xsec;
print " pb \n";

print "Eff filter = ";
print &restrict_num_decimal_digits($efffilter,7);
print " +/- ";
print &restrict_num_decimal_digits($errefffilter,7);
print "\n";

print "Eff matching = ";
print &restrict_num_decimal_digits($effmatch,7);
print " +/- ";
print &restrict_num_decimal_digits($erreffmatch,7);
print "\n";

print "Timing = ";
print &restrict_num_decimal_digits($timing,2);
print " sec/event \n"; 

print "Event size = ";
print &restrict_num_decimal_digits($eventsize,2);
print " kB/event \n";

sub restrict_num_decimal_digits
{
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
