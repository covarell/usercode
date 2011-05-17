#!/usr/bin/perl

# Sum Rates from many files

$workdir = "log";
$thislogbase = $ARGV[0];
$nfiles = $ARGV[1];

my @triggernames;
my @total;
my @passed;

my $itrig = 0;

for($i = 1; $i < $nfiles+1; $i++) {

  open(INFILE,"$workdir/${thislogbase}_$i.log") or die "cannot open any $thislogbase";;
  @log=<INFILE>;
  close(INFILE);

  $itrig = 0;

  # Sum the passing events
  foreach $line (@log) {
    if ($line =~ /HLT-Report/ && $line =~ /0 HLT/) {  
	if ($line =~ /EG/ || $line =~ /Photon/ || $line =~ /Ele/ ) {

            $itrig++;
            @splitline = split(/ +/, $line);

            if ($i == 1) {
		$total[$itrig] = 0;
                $passed[$itrig] = 0; 
                chomp($splitline[5]);
		$triggernames[$itrig] = $splitline[5];
	    }

	    $total[$itrig] = $total[$itrig] + $splitline[2];
            $passed[$itrig] = $passed[$itrig] + $splitline[3]; 
	    # print "$itrig $triggernames[$itrig]  $total[$itrig] $passed[$itrig] \n";
	  }
      }
  }
}

# print "$itrig \n";

for($n = 1; $n < $itrig+1; $n++) {
   my $ratio = $passed[$n]/$total[$n];
   print "$triggernames[$n] : $passed[$n] / $total[$n] = ";
   print &restrict_num_decimal_digits($ratio,5);
   print "\n";
}

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
