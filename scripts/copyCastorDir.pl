#!/usr/bin/perl

# Copy a CASTOR dir or part of it

$workdir = ".";
$thisdata = $ARGV[0];
my $thispart = $ARGV[1];
my $direction = $ARGV[2];
if (!($direction == 1 || $direction == 2)) { 
  print("Argument n. 2 must be: 1 (copy FROM castor) or 2 (copy TO castor) \n");
  exit;
}

my $beginline = 0;
if ($direction == 1) { 
  system("rfdir ${thisdata} >> tmpfile");
  $beginline = 67;
} else {
  system("ls -l >> tmpfile");
  $beginline = 47;
}
open(INFILE,"tmpfile") or die "cannot open tmpfile";;
@log=<INFILE>;
close(INFILE);

# copy files
foreach $line (@log) {
    if ($line =~ m/(${thispart})/) {   
        my $partline = substr($line,$beginline);
	my $rightmarker = rindex($partline,'t') - length($partline) + 1;
        my $filename = substr($partline,0,$rightmarker);
        
        if ($direction == 1) { 
	    print "rfcp ${thisdata}/${filename} .\n";
	    system("rfcp ${thisdata}/${filename} .");
	} else {
            print "rfcp ${filename} ${thisdata} \n";
	    system("rfcp ${filename} ${thisdata}");
	}
    }
}

system("rm -f tmpfile");
print "Directory copied! Bye bye. \n";
 
