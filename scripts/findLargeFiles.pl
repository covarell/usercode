#!/usr/bin/perl

# Find which files are causing AFS quota full

$workdir = ".";
$thispart = "/";

# produce a long file
system("du -h >> tmpfile");
open(INFILE,"tmpfile") or die "cannot open tmpfile";;
@log=<INFILE>;
close(INFILE);

# send staging requests
foreach $line (@log) {
    my $mystring = $line;
    if ($mystring =~ m/(${thispart})/) {
        my $vindex = index($mystring, $thispart);   
        my $partstring = substr($mystring,$vindex + 1);
        if (!($partstring =~ m/(${thispart})/)) {    
             print "$line"; 
        }
    }
}

system("rm -f tmpfile");
 
