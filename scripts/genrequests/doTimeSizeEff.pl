#!/usr/bin/perl

# Calculate filter efficiency and time per event

$workdir = ".";
$thislog = $ARGV[0];

open(INFILE,"$workdir/${thislog}.log") or die "cannot open ${thislog}.log!";;
@log=<INFILE>;
close(INFILE);
 
my $passedevents = 0; 
my $efffilter = 0;
my $timing = 0;
my $eventsize = 0;

# Calculate the filter eff
foreach $line (@log) {
    if ($line =~ /TrigReport/ && $line =~ /Events total/) {  
	@splitline = split(/ +/, $line);
        # print "$splitline[7]\n";
        $passedevents = $splitline[7];	
	$efffilter = $splitline[7]/$splitline[4];
    }
}
print "Eff filter = ";
print &restrict_num_decimal_digits($efffilter,7);
print "\n";

# Calculate the timing
foreach $line (@log) {
    if ($line =~ /TimeReport/ && $line =~ /CPU/ && $line =~ /event/) {  
	@splitline = split(/ +/, $line);	
	$timing = $splitline[3]/$efffilter;
    }
}
print "Timing = ";
print &restrict_num_decimal_digits($timing,2);
print " sec/event \n"; 

# Calculate the event size
system("ls -l ${thislog}.root > dummyfile");

open(INFILE2,"$workdir/dummyfile") or die "cannot open dummyfile!";;
@log2=<INFILE2>;
close(INFILE2);
foreach $line (@log2) {  
    @splitline = split(/ +/, $line);	
    my $filesize = $splitline[4];
    $eventsize = $filesize/(1000*$passedevents); 
}
system("rm -f dummyfile");

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
