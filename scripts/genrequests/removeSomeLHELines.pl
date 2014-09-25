#!/usr/bin/perl

# Call CMSDriver in an intelligent way

$filer = $ARGV[0];

system("ls *${filer}* > lista");

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
  @splitline = split(/ +/, $line);
  # print "$splitline[1] \n";
  if ($splitline[1] == "5" || $splitline[1] == "-5" || $splitline[1] == "15" || $splitline[1] == "-15") {}
  else { print OUTFILE $line; }
}

close(OUTFILE);
system("mv tmp $infile");

}
