#!/usr/bin/perl

system("du *.lhe > pwdfile");
open(INFILE2,"pwdfile") or die "cannot open pwdfile!";;
@log2=<INFILE2>;
close(INFILE2);
system("rm -f pwdfile");

foreach $line (@log2) {
    chomp($line);
    @splitline = split('h', $line);
    my $space = substr($splitline[0],2,1);
    # print "$splitline[0] $splitline[1] $space \n";
    if ($space == 'l') {
        my $filename = 'lh' . $splitline[1] . 'he';
        print "$filename \n";
	system("rm -f $filename");
    }
}
