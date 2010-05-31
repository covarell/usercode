#!/usr/bin/perl

# Copy a CASTOR dir or part of it

$fitcommand = "FitMassRange -f totalDataSet_Mu3_05pb.root -c";
$ptfile = "pranges.txt";
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
 
