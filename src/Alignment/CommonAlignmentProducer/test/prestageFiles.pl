#!/usr/bin/perl

# Pre-stage CASTOR files of TIF data

$workdir = ".";
$thisdata = $ARGV[0];

# retrieve cfi file
open(INFILE,"../data/data-FNAL-$thisdata.cfi") or die "cannot open ../data/data-FNAL-$thisdata.cfi";;
@log=<INFILE>;
close(INFILE);

# send staging requests
foreach $line (@log) {
    if ($line =~ /store/) {
	my $leftmarker = index($line,'s');     
        my $partline = substr($line,$leftmarker);
	my $rightmarker = rindex($partline,'t') - length($partline) + 1;
        my $filename = substr($partline,0,$rightmarker);
        $filename = "/castor/cern.ch/cms/" . $filename;
        
	print "stager_get -M $filename \n";
        system("stager_get -M $filename");
    }
}

# wait till all are staged
sleep(15);
my $okstage = 0;

LOOP: {

foreach $line (@log) {
    if ($line =~ /store/) {
	my $leftmarker = index($line,'s');     
        my $partline = substr($line,$leftmarker);
	my $rightmarker = rindex($partline,'t') - length($partline) + 1;
        my $filename = substr($partline,0,$rightmarker);
        $filename = "/castor/cern.ch/cms/" . $filename;

        # print "stager_qry -M $filename \n";
        system("stager_qry -M $filename > tmpstat");

        open(TMPSTAT,"tmpstat") or die "cannot open tmpstat";;
	@stat=<TMPSTAT>;
	close(TMPSTAT);

        $okstage = 0;
        foreach $line1 (@stat) {
	    if ($line1 =~ /STAGED/) {$okstage = 1;}
	}
	system("rm -f tmpstat");

        if ($okstage == 0) {
	    print "Still staging in $filename ... \n";
            sleep(60);
            redo LOOP;
	}   
    }
}
}

print "Ok, all done! Bye bye. \n";
 
