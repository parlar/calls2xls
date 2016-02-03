#!/usr/bin/perl -w
#
# $Id$
# $Revision$
# $Date$

use Spreadsheet::WriteExcel;
use List::Util qw(max min);
use Config::General;

die "cov2xls.pl <sample.perase.coverage.bed> <disease file> <disease> <output file xls>\n" if (!(@ARGV));
die "cov2xls.pl <sample.perase.coverage.bed> <disease file> <disease> <output file xls>\n" if ( $#ARGV != 3 );

chomp($perbase = $ARGV[0]);
chomp($diseasef = $ARGV[1]);
chomp($disease = $ARGV[2]);
chomp($xls_name = $ARGV[3]);


## Read configs
my $conf = Config::General->new("$ENV{'CALLS2XLS'}/calls2xls.cfg");
my %c = $conf->getall;

open(FH0, "<$diseasef");
chomp(@GeneL = <FH0>);
close(FH0);

my $curr_d = "none";
my %disease2genes;
foreach my $row (@GeneL) {
	if($row =~ /^>(\S+)/) {
		$curr_d = $1;
		next;
	}
	if($row !~ /^>/){
#		print "$curr_d $row\n";
		my ($gene,$core,$vartype) = split(/\s+/, $row);
#		print "$curr_d $gene,$core,$vartype\n";
		$disease2genes{$curr_d}{$gene}{'core'} = $core;
		$disease2genes{$curr_d}{$gene}{'vartype'} = $vartype;
	}
}


open(FH0, "<$perbase");
chomp(@RCF = <FH0>);
close(FH0);



%genecov = ();
%allcov = ();
%genecov_30 = ();	

foreach $row (@RCF){

	@tmp = split(/\t/, $row);
	$genecov{$tmp[3]}{addcov} += $tmp[5];
	$genecov{$tmp[3]}{rows} += 1;
	$genecov{$tmp[3]}{chromosome} = $tmp[0];
	
	if($tmp[5]>=30){
		$genecov_30{$tmp[3]}{addcov} += $tmp[5];
		$genecov_30{$tmp[3]}{rows} += 1;
		$genecov_30{$tmp[3]}{chromosome} = $tmp[0];
	}
	push( @{ $genecov{$tmp[3]}{positions}}, $tmp[1]);
	push( @{ $genecov{$tmp[3]}{positions}}, $tmp[2]);
	
	$allcov{addcov} += $tmp[5];
	$allcov{rows} += 1;
}

my $wb = Spreadsheet::WriteExcel->new($xls_name);

$format{'5'} = $wb->add_format(bg_color=>$clr{'red'},border   => 1);
$format{'4'} = $wb->add_format(bg_color=>$clr{'pink'},border   => 1);
$format{'3'} = $wb->add_format(bg_color=>$clr{'yellow'},border   => 1);
$format{'2'} = $wb->add_format(bg_color=>$clr{'lgreen'},border   => 1);
$format{'1'} = $wb->add_format(bg_color=>$clr{'green'},border   => 1);
$format{'0'} = $wb->add_format(bg_color=> 9,border   => 1, );
$format{'0b'} = $wb->add_format(bg_color=> 9,border   => 1, );
$format{'0b'} -> set_bold();


$ws2 = $wb->add_worksheet('selected genes');
$ws2->write(0,0, 'genes in disease hypothesis file', $format{'0b'});
$ws2->write(0,1, 'disease hypothesis', $format{'0b'});
$ws2->write(0,2, 'exists on panel', $format{'0b'});
$ws2->write(0,3, 'core', $format{'0b'});
$ws2->write(0,4, 'vartype', $format{'0b'});    
	
$wsrowno = 1;

foreach $gene (sort { $a cmp $b } keys(%{$disease2genes{$disease}})){
#	print "$gene\n";
	
	if($disease2genes{$disease}{$gene}{'core'} eq "+") {
		$ft_tmp = $format{'0b'};
	}
	else{
		$ft_tmp = $format{'0'};
	}
		
	$gene =~ s/\s+//g;
	$ws2->write($wsrowno,0, $gene, $ft_tmp);
	$ws2->write($wsrowno,1, $disease, $ft_tmp);
	
	if(exists($genecov{$gene})) {
		$ws2->write($wsrowno,2, "yes", $ft_tmp);
	}
	else{
		$ws2->write($wsrowno,2, "no", $ft_tmp);
	}
	if($disease2genes{$disease}{$gene}{'core'} eq "+") {
		$ws2->write($wsrowno,3, "+", $ft_tmp);
	}
	else{
		$ws2->write($wsrowno,3, "-", $ft_tmp);
	}
	$ws2->write($wsrowno,4, $disease2genes{$disease}{$gene}{'vartype'}, $ft_tmp);
	$wsrowno++;
}

$ws4 = $wb->add_worksheet('genes coverage & completeness');
$ws4->write(0,0, 'chrom',$format{'0b'});
$ws4->write(0,1, 'start',$format{'0b'});
$ws4->write(0,2, 'end',$format{'0b'});
$ws4->write(0,3, 'gene',$format{'0b'});
$ws4->write(0,4, 'core',$format{'0b'});
$ws4->write(0,5, 'av_cov',$format{'0b'});
$ws4->write(0,6, 'completeness_30x',$format{'0b'});	

$gci = 1;
foreach $gene (sort { $a cmp $b } keys (%genecov)) {
	print "$gene\n";
	$gene =~ s/\s+//g;
	if (exists($disease2genes{$disease}{$gene})) {
		if($disease2genes{$disease}{$gene}{'core'} eq "+") {
			$ft_tmp = $format{'0b'};
		}
		else{
			$ft_tmp = $format{'0'};
		}

		my $chr = $genecov{$gene}{chromosome};
		my $max = max @{$genecov{$gene}{positions}};
		my $min = min @{$genecov{$gene}{positions}};
		my $cov = sprintf("%.2f", $genecov{$gene}{addcov}/$genecov{$gene}{rows});
		
		if(!exists($genecov_30{$gene}{rows})){
			$genecov_30{$gene}{rows} = 0;
		}

		$completeness = sprintf("%.2f", $genecov_30{$gene}{rows}/$genecov{$gene}{rows}*100) . "%";
		
		$ws4->write($gci,0, $chr,$ft_tmp);
		$ws4->write($gci,1, $max,$ft_tmp);
		$ws4->write($gci,2, $min,$ft_tmp);
		$ws4->write($gci,3, $gene,$ft_tmp);
		$ws4->write($gci,4, $disease2genes{$disease}{$gene}{'core'},$ft_tmp);
		$ws4->write($gci,5, $cov,$ft_tmp);
		$ws4->write($gci,6, $completeness,$ft_tmp);		
		
		$gci++;
	}
}
