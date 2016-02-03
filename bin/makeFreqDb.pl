#!/usr/bin/perl -w
#

use Sort::Key::Natural qw(natsort);
use Sort::Key::Maker chr_pos_sort => qw(natural integer);
use Config::General;
use MCE::Map
      max_workers =>16,             ## Default 'auto'
      chunk_size => 10              ## Default 'auto'
;

die "makeFreqDb.pl <in.vardb> <in.covdb> <out.freqdb>\n" if (!(@ARGV));
die "makeFreqDb.pl <in.vardb> <in.covdb> <out.freqdb>\n" if ( $#ARGV != 2 );

chomp($vardb = $ARGV[0]);
chomp($covdb = $ARGV[1]);
chomp($freqdb = $ARGV[2]);

$conf = Config::General->new("$ENV{'CALLS2XLS'}/calls2xls.cfg");
%c = $conf->getall;

open FH, "gunzip -c $vardb | ";
while(<FH>){
	push @VCF,$_ unless /^#/;
}
close(FH);


sub getcov{

    @vcfrow = split(/\s/, $_);

    print "$_\n";
    
	%covSamples = ();
	@covout = ();
	$all_sample_found_alleles = 0;
	$all_sample_tot_alleles = 0;
    $position = $vcfrow[0] . ":" . $vcfrow[1] . "-" . $vcfrow[1];
    print "$position\n";
    
    unless (open(GFH, "tabix $covdb $position 2>&1 |")) {
		die "Can't spawn external command!";
    }
    chomp(@covout = <GFH>);
    close(GFH);
    
    foreach $row (@covout){
		@cov = split(/\t/, $row);
		$id = $cov[5] . ":" . $cov[6] . ":" . $cov[7];
		if($cov[3] >= $c{'min_cov'}) {
			if($cov[0] eq "chrM"){
				$all_sample_tot_alleles += 1;
				$covSamples{$id}{'tot_alleles'} = 1;
			}
			elsif($cov[0] eq "chrY" && $cov[8] eq "male"){
				$all_sample_tot_alleles += 1;
				$covSamples{$id}{'tot_alleles'} = 1;
			}
			elsif($cov[0] eq "chrX" && $cov[8] eq "female"){
				$all_sample_tot_alleles += 2;
				$covSamples{$id}{'tot_alleles'} = 2;
			}
			elsif($cov[0] eq "chrX" && $cov[8] eq "male"){
				$all_sample_tot_alleles += 1;
				$covSamples{$id}{'tot_alleles'} = 1;
			}
			else{
				$all_sample_tot_alleles += 2;
				$covSamples{$id}{'tot_alleles'} = 2;
			}
		}
	}

	@entries = split(/;/, $vcfrow[7]);
	$all_sample_found_alleles = 0;
	foreach $row (@entries){
		if($row =~ /=/){
			($id, $value) = split(/=/, $row);
			if(exists($covSamples{$id})){
				$value =~ /(\d+)\/(\d+)/;
				$all_sample_found_alleles += $covSamples{$id}{'tot_alleles'} * $1 / $2;
			}
		}
	}
	if($all_sample_found_alleles > 0){
		$freq = sprintf("%.4f",  $all_sample_found_alleles/$all_sample_tot_alleles);
		$vcfrow[7] = "AF=$freq;AC=$all_sample_found_alleles;AN=$all_sample_tot_alleles";
		$vcfrow[2] = ".";
		$vcfrow[5] = "100";  

		return join("\t", @vcfrow);
	}
	else{
		$vcfrow[7] = "AF=0;AC=$all_sample_found_alleles;AN=$all_sample_tot_alleles";
		$vcfrow[2] = ".";
		$vcfrow[5] = "100";

		return join("\t", @vcfrow);
	}
}
@results =  mce_map { getcov } @VCF;
chr_pos_sort_inplace {/^(\S+)\t(\d+)/; $1, $2} @results;

open FH2, "| bgzip -c >$freqdb ";
print FH2 "##fileformat=VCFv4.1\n";
print FH2 "##INFO=<ID=LOC_AF,Number=A,Type=Float,Description=\"Alternate allele frequency\">\n";
print FH2 "##INFO=<ID=LOC_AC,Number=1,Type=Integer,Description=\"Alternate allele count\">\n";
print FH2 "##INFO=<ID=LOC_AN,Number=A,Type=Integer,Description=\"Chromosome count\">\n";
print FH2 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";


foreach $row (@results){	
	print FH2 "$row\n";
}
close(FH2);
system("tabix -p vcf $freqdb");
