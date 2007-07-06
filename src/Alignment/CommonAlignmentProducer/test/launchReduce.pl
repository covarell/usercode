#!/usr/bin/perl

# Launch "skimming jobs" on TIF data

$workdir = ".";
$thisdata = $ARGV[0];
$istest = $ARGV[1];
$originalcfg = "ReduceEventsCTF.cfg";

# create cfg file
system("cp ${workdir}/$originalcfg ${workdir}/tmp.cfg");

$repl1= "             untracked vstring destinations = { \"cout\", \"reduce-" . $thisdata . "\" } \n           untracked vstring statistics = { \"cout\", \"reduce-" . $thisdata . "\" }";

$repl2= "             include \"Alignment/CommonAlignmentProducer/data/data-FNAL-". $thisdata . ".cfi\" \n             replace PoolSource.maxEvents = " ;
if ($istest =~ /test/ ) { $repl2= $repl2 . "500" ;} 
else { $repl2= $repl2 . "-1" ;}

$repl3=  "             untracked string fileName = \"/data/covarell/inputfiles/dataRED-" . $thisdata . "-CTF.root\"" ;
 
open(INFILE,"tmp.cfg") or die "cannot open tmp.cfg";;
@log=<INFILE>;
close(INFILE);

system("rm -f tmp");
open(OUTFILE,">tmp");

foreach $line (@log) {
  if ($line =~ /REPLACELOG/) { print OUTFILE $repl1; }
  elsif ($line =~ /REPLACEIN/) { print OUTFILE $repl2; }
  elsif ($line =~ /REPLACEOUT/) { print OUTFILE $repl3; }
  else { print OUTFILE $line; }
}

close(OUTFILE);
system("mv tmp tmp.cfg");

system("cmsRun tmp.cfg");

