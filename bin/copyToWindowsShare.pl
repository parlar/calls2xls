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
	@spl = split(/\//, $samplef);
	push @spl, $spl[-1];
	
	$from = join('/', @spl);
	$from =~ s/\/\//\//g;
	print "$from\n";
	system("rsync -r --ignore-existing $from $winshare");
    }
}

system("winshare_date_no.pl 1");
