#!/usr/bin/perl

# make a nice list of files 

$workdir = ".";
$thisdata = $ARGV[0];
$thispart = $ARGV[1];
$textfilename = $ARGV[2];
$forpython = $ARGV[3];
$filesperline = $ARGV[4];

if (!($forpython == 1 || $forpython == 2)) { 
  print("Argument n. 4 must be: 1 (a simple list) or 2 (a list for the PoolSource in a python cfg) \n");
  exit;
} 

if ($thisdata =~ /eos/) { 
  print("Argument n. 1 must be the directory WITHOUT /eoscms/, i.e. starting with /store/ \n");
  exit;
} 
                 
# retrieve cfi file
system("cmsLs ${thisdata} >> tmpfile");
open(INFILE,"tmpfile") or die "cannot open tmpfile";;
@log=<INFILE>;
close(INFILE);

my $correctpath = '';
open(MYFILE,">${textfilename}");;

my $ilines = 0;
my $totfilename = '';

foreach $line (@log) {
   if ($line =~ /dr/) {
      my $dirname = find_dir_in_line($line);
      print "dirname $dirname \n";
      system("cmsLs ${dirname} >> tmpfile2");
      # system("cmsLs ${thisdata}/${dirname} >> tmpfile2");
      open(INFILE2,"tmpfile2") or die "cannot open tmpfile2";;
      @log2=<INFILE2>;
      close(INFILE2);
      foreach $line2 (@log2) {
	 if ($line2 =~ /dr/) { 
	    my $dirname2 = find_dir_in_line($line2);
	    print "dirname2 $dirname2 \n";
	    # system("cmsLs ${thisdata}/${dirname}/${dirname2} >> tmpfile3");
	    system("cmsLs ${dirname2} >> tmpfile3");
	    open(INFILE3,"tmpfile3") or die "cannot open tmpfile3";;
	    @log3=<INFILE3>;
	    close(INFILE3);
	    foreach $line3 (@log3) {
	       if ($line3 =~ /dr/ && $line3 =~ /${thispart}/) {	
		  my $dirname3 = find_dir_in_line($line3);
		  print "dirname3 $dirname3 \n";
		  # system("cmsLs ${thisdata}/${dirname}/${dirname2}/${dirname3} >> tmpfile4");
		  system("cmsLs ${dirname3} >> tmpfile4");
		  open(INFILE4,"tmpfile4") or die "cannot open tmpfile4";;
		  @log4=<INFILE4>;
		  close(INFILE4);
		  # system("edmFileUtil -d ${thisdata}/${dirname}/${dirname2}/${dirname3} >> utilfile");
		  foreach $line4 (@log4) {
		     if ($line4 =~ /root/) { 
			my $filename = find_dir_in_line($line4);
			# my $filename = find_file_in_line($line4);
			# print "filename $filename \n";
			if ($forpython == 1) {
			    $totfilename = "/" . $filename;
			} else { 
			    system("edmFileUtil -d ${filename} >> utilfile");
			    open(UTFILE,"utilfile") or die "cannot open utilfile";;
			    @util=<UTFILE>;
			    close(UTFILE);
			    foreach $utline (@util) {
				chomp($utline);
				# print "utline $utline \n";
				$correctpath = $utline;
			    }
			    if ($correctpath =~ /root/) { 
				$totfilename = $totfilename . "\'" . $correctpath . "\',";
				$ilines = $ilines + 1;
			    }
			} 
		   
			if ($ilines % $filesperline == 0) {
			    print MYFILE "$totfilename \n";
			    print "$totfilename \n";    
			    $totfilename = '';
			}
			system("rm -f utilfile");
		    }
		 }
		 system("rm -f tmpfile4");
	      }
	   }
	   system("rm -f tmpfile3");
	}
     }
     system("rm -f tmpfile2");
  }
}

system("rm -f tmpfile");
print "Result in ${textfilename}\n";
 
## No fixed length, look for the t of root in the name
sub find_file_in_line
{
    my $line=shift;
    my $partline = substr($line,43);
    my $rightmarker = rindex($partline,'t') - length($partline) + 1;
    my $filename = substr($partline,0,$rightmarker);
    chomp($filename);
    return $filename;
}

## Fixed lenght, 3 chars
sub find_dir_in_line
{
    my $line=shift;
    my $dirname = substr($line,43);
    # my $partline = substr($line,43);
    # my $dirname = substr($partline,0,4);
    chomp($dirname);
    return $dirname;
}
