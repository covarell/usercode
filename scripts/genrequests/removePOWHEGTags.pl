#!/usr/bin/perl

# Call CMSDriver in an intelligent way

$workdir = ".";

system("ls $workdir > lista");

open(INFILE,"lista") or die "cannot open lista";;
@log=<INFILE>;
close(INFILE);

foreach $line (@log) {
  replace("${line}");
}

system("rm -f lista");

###############################################################################

sub replace {

$infile = @_[0];
$repl = @_[1];

open(INFILE,"$infile") or die "cannot open $infile";;
@log2=<INFILE>;
close(INFILE);

system("rm -f tmp");
open(OUTFILE,">tmp");

foreach $line (@log2) {
  if ($line =~ /xit values/) { print "${line} \n"}
  else { print OUTFILE $line; }
}

close(OUTFILE);
system("mv tmp $infile");

}
