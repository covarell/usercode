#!/usr/bin/perl

# CmsStage a EOS dir or part of it

$workdir = ".";
my $newName = $ARGV[0];
$upload = $ARGV[1];
$artNum = $ARGV[2];
$ii = 0;

while ($ii < 10) {
     system("ls ${newName}*$ii.lhe > lista$ii");
     system("~/myscripts/genRequests/mergeLheFiles lista$ii");
     system("mv out.lhe ${newName}Merged_$ii.lhe");
     # if ($upload == 1) system("cmsStage -f ${newName}Merged_$ii.lhe /store/lhe/$artNum/");
     $ii++;
}

if ($upload == 1) {
    print("~/scratch0/mcrequests/CMSSW_5_3_3/src/GeneratorInterface/LHEInterface/scripts/cmsLHEtoEOSManager.py -n $artNum -f ${newName}Merged_0.lhe,${newName}Merged_1.lhe,${newName}Merged_2.lhe,${newName}Merged_3.lhe,${newName}Merged_4.lhe,${newName}Merged_5.lhe,${newName}Merged_6.lhe,${newName}Merged_7.lhe,${newName}Merged_8.lhe,${newName}Merged_9.lhe \n");
    system("~/scratch0/mcrequests/CMSSW_5_3_3/src/GeneratorInterface/LHEInterface/scripts/cmsLHEtoEOSManager.py -n $artNum -f ${newName}Merged_0.lhe,${newName}Merged_1.lhe,${newName}Merged_2.lhe,${newName}Merged_3.lhe,${newName}Merged_4.lhe,${newName}Merged_5.lhe,${newName}Merged_6.lhe,${newName}Merged_7.lhe,${newName}Merged_8.lhe,${newName}Merged_9.lhe"); 
}
