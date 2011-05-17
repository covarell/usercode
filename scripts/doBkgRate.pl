#!/usr/bin/perl

# Calculate all MB rates

#### Lumi = 8e29 / 1e31
#### sigma(pp) = 51.5 mb
#### e1 = e_skim(L1SingleEG5) = 0.553%
#### e2 = e_skim(L1SingleEG5 or DoubleEG5) = 0.879%

#### R = Lumi * sigma(pp) * e_skim(L1SingleEG5) * e_hlt
####   = (227.8 / 2848) * e_hlt

#### R = Lumi * sigma(pp) * e_skim(L1SingleEG5 or DoubleEG3) * e_hlt
####   = (362.1 / 4527) * e_hlt

$workdir = ".";
$thislog = $ARGV[0];
my $thislumi = $ARGV[1];
my $skim = $ARGV[2];

if (!($thislumi eq "8E29") && !($thislumi eq "1E31")) {
    die "Argument 2 must be either 8E29 or 1E31!"
}

if (!($skim eq "L") && !($skim eq "T")) {
    die "Argument 3 must be either T (tight skim: SingleEG5) or L (loose skim: SingleEG5 or DoubleEG3)!"
} 

open(INFILE,"$workdir/$thislog") or die "cannot open $thislog";;
@log=<INFILE>;
close(INFILE);


print "Rates for the $thislumi menu \n \n";

# Calculate the rates
foreach $line (@log) {
    if ($line =~ /HLT-Report/ && $line =~ /0 HLT/) {  
	if ($line =~ /EG/ || $line =~ /Photon/ || $line =~ /Ele/ ) {

	    @splitline = split(/ +/, $line);

	    my $effhlt = $splitline[3]/$splitline[2];
            my $errhlt = sqrt($splitline[3])/$splitline[2];
            my $rate = 0.;  my $errrate = 0.;
            if ($thislumi eq "8E29") {
		$rate = $effhlt * 227.8; 
		$errrate = $errhlt * 227.8;
		if ($skim eq "L") {
		    $rate = $effhlt * 362.1; 
		    $errrate = $errhlt * 362.1;
		}
	    } else {
         	$rate = $effhlt * 2848;
                $errrate = $errhlt * 2848;
		if ($skim eq "L") {
		    $rate = $effhlt * 4527; 
		    $errrate = $errhlt * 4527;
		}
	    }
            
	    chomp($splitline[5]);
            print "$splitline[5] : (";
	    print &restrict_num_decimal_digits($rate,3);
	    print " +/- ";
	    print &restrict_num_decimal_digits($errrate,3);
	    print ") Hz \n";
	}
    }
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
