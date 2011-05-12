#!/usr/bin/perl

# Launch many fits

# $fitcommand = "./Fit2DSimMC -m datasets/MCFall10_jpsi_DMu0.root -c datasets/MCFall10_psip_DMu0.root -f datasets/MCFall10_QCD_DMu0.root -u -b -s";
# $fitcommand = "./Fit2DSimRange -m datasets/MCFall10_jpsi_DMu0.root -c datasets/MCFall10_psip_DMu0.root -u -b -s";
# $fitcommand = "./Fit2DSimPEE -m datasets/MCFall10_jpsi_DMu0.root -c datasets/MCFall10_psip_DMu0.root -u -b -s";
# $fitcommand = "./Fit2DJpsiPEE -m datasets/MCFall10_jpsi_DMu0.root -u -b -s";
$fitcommand = "./FitMassSim";

$ptfile = $ARGV[0];
$etafile = $ARGV[1];
$thequeue = "8nh";

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

	# my $currentfit = $fitcommand . ' -p ' . $splitline[0] . ' -l ' . $splitline[1] . ' -r ' . $splitline[2] . ' -y ' . $splitline2[0] .' -d datasets/Data2010_rap' . $splitline2[1] . '.root';
        my $currentfit = $fitcommand . ' -p ' . $splitline[0] . ' -y ' . $splitline2[0] .' -d datasets/DataMergedWide2010_rap' . $splitline2[1] . '.root';

        print "SUBMITTING ${currentfit} ... \n";
	$subjob="#!/bin/zsh 
#BSUB -J \"FIT\" 
#BSUB -C 0
cd /afs/cern.ch/user/c/covarell/scratch0/quarkonia/CMSSW_3_8_6/roofit/UserCode/Covarell/src/JpsiMacros
export SCRAM_ARCH=slc5_ia32_gcc434
cmsenv
${currentfit}
";

        my $scriptname = 'submit_' . $splitline[0] . '_' . $splitline2[0];
        open(FILE,">scripts/${scriptname}");
        print FILE "$subjob";
        close(FILE);
        system("chmod u+x scripts/${scriptname}");

	$test=0;
	while($test == 0) {
	    $rc=system("bsub -o log/${scriptname}.log -q $thequeue < scripts/${scriptname}");
	    if ($rc == 0) { $test=1; }
	    else {
		print "ERROR in submitting job: $rc  .. retrying in 10s ...\n";
		sleep(5);
	    }
	}
	sleep(1);
    }
}
print "All fits submitted! Bye bye. \n";
 
