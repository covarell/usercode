#!/usr/bin/perl

# Launch many fits

# $fitcommand = "FindMeans -f datasets/DataSet_314nb_PV.root -t 0";
# $fitcommand = "Fit2DDataSyst -f datasets/DataSet_314nb_PV.root -m datasets/totalDataSet_allTriggers_05pb.root -t 0 -u -b -c -s 1";
$fitcommand = "Fit2DSimRange -m datasets/MCFall10_jpsi_DMu0.root -c datasets/MCFall10_psip_DMu0.root -u -b -s";
$ptfile = $ARGV[0];
$etafile = $ARGV[1];

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
        @splitline = split(/ +/, $line);
	chomp($line2);
        @splitline2 = split(/ +/, $line2);

        my $currentfit = $fitcommand . ' -p ' . $splitline[0] . ' -l ' . $splitline[1] .' -y ' . $splitline2[0] .' -f datasets/Data2010_rap' . $splitline2[1] . '.root';
        # my $currentfit = $fitcommand . ' -p ' . $splitline[0] .' -e ' . $line2;
        print "${currentfit} \n";
	system("${currentfit}");
    }
}
print "All fits done! Bye bye. \n";
 
