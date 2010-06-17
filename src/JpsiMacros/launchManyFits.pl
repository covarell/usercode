#!/usr/bin/perl

# Launch many fits

$fitcommand = "FitMassDataRange -f DatiDMuOpen.root";
# $fitcommand = "Fit2DRange -f totalDataSet_allTriggers_05pb.root -g 1 -t 0 -b -u -c";
$ptfile = "prangesdata.txt";
$etafile = "etaranges.txt";

open(INFILE,${ptfile}) or die "cannot open ${ptfile}";;
@log=<INFILE>;
close(INFILE);

open(INFILE2,${etafile}) or die "cannot open ${etafile}";;
@log2=<INFILE2>;
close(INFILE2);

# loop on files
foreach $line2 (@log2) {
    foreach $line (@log) {
	chomp($line);
	chomp($line2);
        my $currentfit = $fitcommand . ' -p ' . $line . ' -e ' . $line2;
        print "${currentfit} \n";
	system("${currentfit}");
    }
}
print "All fits done! Bye bye. \n";
 
