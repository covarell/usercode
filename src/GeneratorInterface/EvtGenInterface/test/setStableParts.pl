#!/usr/bin/perl

# Replace stable particles in Pythia setup

$workdir = ".";

$originalcfg = "urs.cfg";
$finalcfg = "pythiabb.cfg";
$partlist = "/afs/cern.ch/user/c/covarell/scratch0/evtgen/CMSSW_1_3_1/src/GeneratorInterface/EvtGenInterface/test/partlist.txt";

# create cfg file
system("cp ${workdir}/$originalcfg ${workdir}/$finalcfg");

open(PARTLIST,"$partlist") or die "cannot open $partlist";;
@part=<PARTLIST>;
close(PARTLIST);
$repl="";

foreach $number (@part) {
  chop($number);
  if ($number > 0) {
      $repl=$repl . "       'MDCY(PYCOMP(" . $number . "),1) = 0', \n ";
  }
}

open(INFILE,"$finalcfg") or die "cannot open $finalcfg";;
@log=<INFILE>;
close(INFILE);

system("rm -f tmp");
open(OUTFILE,">tmp");

foreach $line (@log) {
  if ($line =~ /REPLACEME/) { print OUTFILE $repl; }
  else { print OUTFILE $line; }
}

close(OUTFILE);
system("mv tmp $finalcfg");

