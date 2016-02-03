#!/usr/bin/perl -w
#

die "copyToWindowsShare.pl <miseq folder>\n" if (!(@ARGV));
die "copyToWindowsShare.pl <miseq_folder>\n" if ( $#ARGV != 0 );

$folder = $ARGV[0];
$winshare = "/media/windowsshare/GML/NGS/NGS_calls";
$rootfolder = $folder . "/calling/final";

@samplefolders = <$rootfolder/*>;

foreach $samplef (@samplefolders) {
    if($samplef !~ /project/) {
        @samplefolderyears = <$samplef/*>;
        foreach $sfy (@samplefolderyears){
            if (-d $sfy) {
                @tmp = split(/\//, $sfy);
                if($tmp[-1] =~ /^\d\d\d\d$/){
                    print "$sfy\n";
                	system("rsync -r --ignore-existing $sfy $winshare");
                }
            }
        }
    }
}

#system("winshare_date_no.pl 1");
