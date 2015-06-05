#!/usr/bin/perl -w
#
use Spreadsheet::WriteExcel;
use Encode;

die "winshare_date_no.pl \n" if (!(@ARGV));

$winshare = "/media/windowsshare/GML/NGS/NGS_calls";

@samples = <$winshare/*>;

foreach $sample (@samples) {
    unless ($sample =~ /perdatum/ || $sample =~ /Planering/){
	@rundates = <$sample/*>;
	foreach $rundate (@rundates){
	    if(-d $rundate) {
		$samp = (split(/\//, $sample))[-1];
		$run1 = (split(/\//, $rundate))[-1];
		$run2 = (split(/\./, $run1))[0];
		$r2s{$run2}{$samp} = 1;
	    }
	}
    }
}

$xls_name = $winshare . "/" . "remissnummer_kördatum.xls";
my $wb = Spreadsheet::WriteExcel->new($xls_name);
$ws = $wb->add_worksheet('remissnummer_rundate');

$ws->set_column('A:A', 18);
$ws->set_column('B:B', 18);


$zero = $wb->add_format(bg_color=> 9,border   => 1, );
$zerobold = $wb->add_format(bg_color=> 9,border   => 1, );
$kd =   decode("utf-8", "kördatum");
@head = ($kd,"remissnummer");
$i = 0;
$ws->write_row($i,0, \@head, $zerobold);
$i+=1;

foreach $r (sort { $b <=> $a } keys(%r2s)) {
    foreach $s (sort { $a cmp $b } keys(%{$r2s{$r}})) {
	@outrow = ($r,$s);
	$ws->write_row($i,0, \@outrow, $zero);
	$i++;
    }
}
$i++;
$cells = "A1:B" . $i;

$ws->autofilter($cells);

$wb->close();

