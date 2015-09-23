#!/usr/bin/perl -w
#

use Time::localtime;
use List::Util qw(max min);

die "mut2vcf.pl <final_path> \n" if (!(@ARGV));
die "mut2vcf.pl <final_path> \n" if ( $#ARGV != 0 );

chomp($final = $ARGV[0]);
$final =~ s/\/$//g;
$of = $final . "/alamut_mut.vcf";

open (FH, "</home/data_in/calls2xls/calls2xls.cfg");
while ( $line = <FH> ) {
    if($line !~ /#/ && $line =~ /\S/){
	@tmp = split(/=/, $line);
	foreach $row (@tmp) {
	    $row =~ s/\s//g;
	}
	$c{$tmp[0]}=$tmp[1];
	if($tmp[0] eq "info_dp"){
	    @info = split(/,/,$tmp[1]);
	    foreach $i (@info){
		$c{$tmp[0]}{$i} = 1;
	    }
	}
	if($tmp[0] eq "info_ao"){
	    @info = split(/,/,$tmp[1]);
	    foreach $i (@info){
		$c{$tmp[0]}{$i} = 1;
	    }
	}
    }
}



$path = $c{'alamutpath'};
chomp(@files = <$path/*.mut>);

$hg19 = $c{'genome'};

$tm = localtime;
$date = sprintf("%04d-%02d-%02d", $tm->year+1900, ($tm->mon)+1, $tm->mday);

foreach $mutfile (@files) {
    open(FH, "<$mutfile");
    chomp(@data = <FH>);
    close(FH);

    $gene = "";
    $assembly = "";
    $baseFrom_g = "";
    $baseTo_g = "";
    $inserted_g = "";

    $gene = ".";
    $assembly = ".";
    $baseFrom_g = ".";
    $baseTo_g = ".";
    $inserted_g = ".";
    $type = ".";
    $posi_g = ".";
    $chr = ".";
    $class_val = ".";
    $class_index = ".";
    $created = ".";
    $updated = ".";
    %crh = ();
    %uph = ();

    foreach $line (@data) {
	if( $line =~ /\s+<Mutation id=/ ) {

	    $line =~ /\s+<Mutation id="(.*?)"/;
	    $id = $1;
	    $line =~ /\s+chr="(.*?)"/;
	    $chr = $1;
	    $line =~ /\s+geneSym="(.*?)"/;
	    $gene = $1;
	    $line =~ /\s+refAssembly="(.*?)"/;
	    $assembly = $1;	    
	    $assembly =~ s/\s//;

	    $baseFrom_g = " ";
	    $baseTo_g = " ";
	    $inserted_g = " ";
	    $class_val = " ";
	    $class_index = " ";

	}

	if($assembly =~ /NCBI/ || $assembly =~ /GRCh/){
	    if($line =~ /<gNomen val="g.(.*)(\w)&gt;(\w)"\/>$/){
		$posi_g = $1;
		$baseFrom_g = $2;
		$baseTo_g = $3;
		$type = "subst";
	    }
	    elsif($line =~ /<gNomen val="g.([0-9_]+)ins(\w+)"\/>$/){
		$posi_g = $1;
		$inserted_g = $2;
		$type = "ins";
	    }
	    elsif($line =~ /<gNomen val="g.([0-9_]+)delins(\w+)"\/>$/){	
		$posi_g = $1;
		$inserted_g = $2;
		$type = "delins";
	    }
	    elsif($line =~ /<gNomen val="g.([0-9_]+)del"\/>$/){
		$posi_g = $1;
		$type = "del";
	    }
	    elsif($line =~ /<gNomen val="g.([0-9_]+)dup"\/>$/){
		$posi_g = $1;
		$type = "dup";
	    }
	    if($line =~ /<Classification\s+/){
		$line =~ /val="(.*?)"/;
		$class_val = $1;
		$line =~ /index="(.*?)"/;
		$class_index = $1;		
	    }
	    if($line =~ /<Created\s+/){
		$line =~ /date="(.*?)"/;
		$created = $1;
		$created =~ s/-//g;
		substr($created, 0, 2, "");
		$uph{$created} = 1;
	    }
	    if($line =~ /<Updated\s+/){
		$line =~ /date="(.*?)"/;
		$updated = $1;	
		$updated =~ s/-//g;
		substr($updated, 0, 2, "");
		$crh{$updated} = 1;
	
	    }
	}
	if( $line =~ /<\/Mutation>/ ) {
	    
	    $sizecrh =  keys %crh;
	    $sizeuph =  keys %uph;

#	    print "$sizecrh $sizeuph\n";

	    if($sizecrh > 0){
		$created = min keys %crh;
	    }
	    else{
		$created = ".";
	    }
	    if($sizeuph > 0){
		$updated = min keys %uph;
	    }
	    else{
		$updated = ".";
	    }
	    


#	    print "$id $assembly $chr $gene $posi_g $baseFrom_g $baseTo_g $inserted_g $type $class_val $class_index\n";
	    $alldat{$assembly}{$id}{chr} = $chr;
	    $alldat{$assembly}{$id}{posi} = $posi_g;
	    $alldat{$assembly}{$id}{type} = $type;
	    $alldat{$assembly}{$id}{created} = $updated;
	    $alldat{$assembly}{$id}{updated} = $updated;
	    if($type eq "subst") {
		$alldat{$assembly}{$id}{baseFrom} = $baseFrom_g;
		$alldat{$assembly}{$id}{baseTo} = $baseTo_g;
	    }
	    elsif($type eq "ins" || $type eq "delins") {
		$alldat{$assembly}{$id}{insert} = $inserted_g;
	    }
	    if($class_val =~ /\S/) {
		$alldat{$assembly}{$id}{clval} = $class_val;
		$alldat{$assembly}{$id}{clind} = $class_index;
	    }
	    %crh = ();
	    %uph = ();

	    $gene = ".";
	    $assembly = ".";
	    $baseFrom_g = ".";
	    $baseTo_g = ".";
	    $inserted_g = ".";
	    $type = ".";
	    $posi_g = ".";
	    $chr = ".";
	    $class_val = ".";
	    $class_index = ".";
	    $created = ".";
	    $updated = ".";

	}
    }
}
	
foreach $ass (sort { $a cmp $b } keys(%alldat)) {
    $bedfile = "tmp." . $ass . ".bed";
    open(FH, ">/tmp/$bedfile");
    if($ass eq "NCBI36") {
	$chain = "/home/data_in/calls2xls/chainfiles/hg18ToHg19.over.chain.gz";
    }
    elsif($ass eq "GRCh38"){
	$chain = "/home/data_in/calls2xls/chainfiles/hg38ToHg19.over.chain.gz";
    }
#    print "$ass $chain\n";
    if($ass eq "GRCh38" || $ass eq "NCBI36"){
	foreach $id (sort { $a cmp $b } keys(%{$alldat{$ass}})) {
	    if($alldat{$ass}{$id}{posi} =~ /(\d+)_(\d+)/){
		$p1 = $1;
		$p2 = $2;
		print FH "chr$alldat{$ass}{$id}{chr}\t$p1\t$p2\t$id\t10\t+\n";
		
	    }
	    elsif($alldat{$ass}{$id}{posi} =~ /(\d+)/){
		$p1 = $1;
#		print "chr$alldat{$ass}{$id}{chr}\t$p1\t$p1\t$id\t10\t+\n";
		print FH "chr$alldat{$ass}{$id}{chr}\t$p1\t$p1\t$id\t10\t+\n";
	    }
	}
	close(FH);
	if($ass eq "NCBI36" || $ass eq "GRCh38"){
	    open(FH2, "CrossMap.py bed $chain /tmp/$bedfile  2>&1 |");
	    chomp(@CM = <FH2>);
	    unless (close(FH2)) {  
		die "External command failed: $?";
	    }
	    foreach $row (@CM) {
#		print "$row\n";
		unless($row =~ /^@/){
#		    print "$row\n";
		    @tmp = split(/\t/, $row);
		    $conv{$tmp[3]} = $row;
		}
	    }
	}
    }
}

foreach $ass (sort { $a cmp $b } keys(%alldat)) {
    foreach $id (sort { $a cmp $b } keys(%{$alldat{$ass}})) {
	
	$id2 = $id;
	$id2 =~ s/\{//g;
	$id2 =~ s/\}//g;
	$mutid = "ALMUTID=" . $id2;
#	print "$ass $id\n";
	$alc = $mutid;
	if(exists($alldat{$ass}{$id}{clind})){
	    $alc_index = ";ALCINDEX=" . $alldat{$ass}{$id}{clind};
	    $alc_val = ";ALCVAL=" . $alldat{$ass}{$id}{clval};
	    $alc .= ";ALAMUT=1" .  $alc_index . $alc_val;
	}
	if(exists($alldat{$ass}{$id}{created})){
	    $alc .= ";ALCREATE=" .  $alldat{$ass}{$id}{created};
	}
	if(exists($alldat{$ass}{$id}{created})){
	    $alc .= ";ALUPDATE=" .  $alldat{$ass}{$id}{updated};
	}
	if(($ass eq "NCBI36" || $ass eq "GRCh38") || $ass eq "GRCh37"){
	    $chr = $alldat{$ass}{$id}{'chr'};
	    $chr_ch = $chr;
	    if(exists($conv{$id})){
		@tmp = split(/\t/, $conv{$id});
		$first = $tmp[8];
		$second = $tmp[9];
#		$dir = $tmp[12]
	    }
	    if($chr eq "MT"){
		$chr = "M";
	    }
	    if($chr eq "X"){
		$chr_ch = 30;
	    }
	    if($chr eq "Y"){
		$chr_ch = 31;
	    }
	    if($chr eq "MT"){
		$chr_ch = 32;
	    }
	    if($chr eq "M"){
		$chr_ch = 33;
	    }
	    if($alldat{$ass}{$id}{type} eq "dup"){
#		print "dup\n";
		if($alldat{$ass}{$id}{posi} =~ /(\d+)_(\d+)/){
		    if(exists($conv{$id})){
			$reg = "chr$chr:$first-$second";
			$start = $first;
		    }
		    else{
			$reg = "chr$chr:$1-$2";
			$start = $1;
		    }
		}
		else{
		    if(exists($conv{$id})){
			$reg = "chr$chr:$first-$second";
			$start = $first;
		    }
		    else{
			$reg = "chr$chr:$alldat{$ass}{$id}{posi}-$alldat{$ass}{$id}{posi}";
			$start = $alldat{$ass}{$id}{posi};
		    }
		}
		
		open(FH2, "samtools faidx $hg19 $reg  2>&1 |");
		chomp(@REG = <FH2>);
		close(FH2);

		shift(@REG);
		$ref = join('', @REG);
		$alt = $ref . $ref;
#		print "$reg\n$ref\n$alt\n$start\n";
		if(length($ref) < 30) {
		    $outdat{$chr_ch}{$start} = "chr$chr\t$start\t.\t$ref\t$alt\t100\t.\t$alc\tGT\t1/1";
		}
	    }
	    
	    if($alldat{$ass}{$id}{type} eq "del"){
#		print "del\n$alldat{$ass}{$id}{posi}\n";
		if($alldat{$ass}{$id}{posi} =~ /(\d+)_(\d+)/){
		    if(exists($conv{$id})){
			$start = $first - 1;
			$end = $second;
			$reg = "chr$chr:$start-$end";
		    }
		    else{
			$start = $1 - 1;
			$end = $2;
			$reg = "chr$chr:$start-$end";
		    }
		}
		else{
		    if(exists($conv{$id})){
			$start = $first - 1;
			$end = $second;
			
		    }
		    else{
			 $start = $alldat{$ass}{$id}{posi} -1;
			 $end = $alldat{$ass}{$id}{posi};
		    }
		    $reg = "chr$chr:$start-$end";
		}
#		print "$reg\n";
		open(FH2, "samtools faidx $hg19 $reg  2>&1 |");

		chomp(@REG = <FH2>);
		close(FH2);
		shift(@REG);

		$ref = join('', @REG);
		$alt = substr($ref, 0, 1);
#		print "$ref $alt\n";
		if(length($ref) < 30) {
		    $outdat{$chr_ch}{$start} = "chr$chr\t$start\t.\t$ref\t$alt\t100\t.\t$alc\tGT\t1/1";
		}
	    }
	    if($alldat{$ass}{$id}{type} eq "ins"){
#		print "ins\n";
		if(exists($conv{$id})){
		    $start = $first;
		    $end = $second;
		}
		else{
		    ($start, $end) = split(/_/, $alldat{$ass}{$id}{posi});
		}
		
		$reg = "chr$chr:$start-$end";
		$insert = $alldat{$ass}{$id}{insert};
#		print "$reg\n";
       
		open(FH2, "samtools faidx $hg19 $reg  2>&1 |");
		chomp(@REG = <FH2>);
		close(FH2);
		shift(@REG);
		
		$ref = join('', @REG);
		$alt = $ref;

#		print "$ref $alt\n";
		substr($alt, 1, 0, $insert);

		if(length($ref) < 20) {
#		    print "chr$chr\t$start\t.\t$ref\t$alt\t.\t$alc\tGT\t1/1\n";
		    $outdat{$chr_ch}{$start} = "chr$chr\t$start\t.\t$ref\t$alt\t100\t.\t$alc\tGT\t1/1";
		}
	    }
	    if($alldat{$ass}{$id}{type} eq "delins"){
#		print "delins\n";
		$insert = $alldat{$ass}{$id}{insert};

		if($alldat{$ass}{$id}{posi} =~ /(\d+)_(\d+)/){
		    if(exists($conv{$id})){
			$start = $first - 1;
			$end = $second;
			$reg = "chr$chr:$start-$end";
		    }
		    else{
			$start = $1 - 1;
			$end = $2;
			$reg = "chr$chr:$start-$end";
		    } 
		}
		else{
		    if(exists($conv{$id})){
			$start = $first - 1;
			$end = $second;
		    }
		    else{
			$start = $alldat{$ass}{$id}{posi} -1;
			$end = $alldat{$ass}{$id}{posi};
		    } 
		}		

		$reg = "chr$chr:$start-$end";
		open(FH2, "samtools faidx $hg19 $reg 2>&1 |");
		chomp(@REG = <FH2>);
		close(FH2);
		shift(@REG);
		$ref = join('', @REG);
		$alt = substr($ref, 0, 1);
		$alt .= $insert;

#		print "$reg\n$ref $alt\n";
		if(length($ref) < 30) {
		    $outdat{$chr_ch}{$start} = "chr$chr\t$start\t.\t$ref\t$alt\t100\t.\t$alc\tGT\t1/1";
		}
	    }
	    if($alldat{$ass}{$id}{type} eq "subst"){
#		print "$alldat{$ass}{$id}{type}\n";
		if(exists($conv{$id})){
		    $start = $first;
		    $end = $second;
		}
		else{
		    $start = $alldat{$ass}{$id}{posi};
		    $end = $alldat{$ass}{$id}{posi};
		} 

		$outdat{$chr_ch}{$start} = "chr$chr\t$start\t.\t$alldat{$ass}{$id}{baseFrom}\t$alldat{$ass}{$id}{baseTo}\t100\t.\t$alc\tGT\t1/1";
#		print "$outdat{$chr_ch}{$start}\n";
#		print "$chr\n";
#		print "$start\n";
#		print "$alldat{$ass}{$id}{baseFrom}\n";
#		print "$alldat{$ass}{$id}{baseTo}\n";
#		print "$alc\n";

	    }
#	    print "$ass $start $id\n---\n";
	}
    }
}

$of2 = $of;
$of = $of . ".tmp.vcf";
open(OF, ">$of");

print OF "##fileformat=VCFv4.1\n";
print OF "##fileDate=$date\n";
print OF "##reference=hg19\n";
print OF "##phasing=none\n";
print OF "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";

print OF "##contig=<ID=chrM,length=16571,assembly=hg19>\n";
print OF "##contig=<ID=chr1,length=249250621,assembly=hg19>\n";
print OF "##contig=<ID=chr2,length=243199373,assembly=hg19>\n";
print OF "##contig=<ID=chr3,length=198022430,assembly=hg19>\n";
print OF "##contig=<ID=chr4,length=191154276,assembly=hg19>\n";
print OF "##contig=<ID=chr5,length=180915260,assembly=hg19>\n";
print OF "##contig=<ID=chr6,length=171115067,assembly=hg19>\n";
print OF "##contig=<ID=chr7,length=159138663,assembly=hg19>\n";
print OF "##contig=<ID=chr8,length=146364022,assembly=hg19>\n";
print OF "##contig=<ID=chr9,length=141213431,assembly=hg19>\n";
print OF "##contig=<ID=chr10,length=135534747,assembly=hg19>\n";
print OF "##contig=<ID=chr11,length=135006516,assembly=hg19>\n";
print OF "##contig=<ID=chr12,length=133851895,assembly=hg19>\n";
print OF "##contig=<ID=chr13,length=115169878,assembly=hg19>\n";
print OF "##contig=<ID=chr14,length=107349540,assembly=hg19>\n";
print OF "##contig=<ID=chr15,length=102531392,assembly=hg19>\n";
print OF "##contig=<ID=chr16,length=90354753,assembly=hg19>\n";
print OF "##contig=<ID=chr17,length=81195210,assembly=hg19>\n";
print OF "##contig=<ID=chr18,length=78077248,assembly=hg19>\n";
print OF "##contig=<ID=chr19,length=59128983,assembly=hg19>\n";
print OF "##contig=<ID=chr20,length=63025520,assembly=hg19>\n";
print OF "##contig=<ID=chr21,length=48129895,assembly=hg19>\n";
print OF "##contig=<ID=chr22,length=51304566,assembly=hg19>\n";
print OF "##contig=<ID=chrX,length=155270560,assembly=hg19>\n";
print OF "##contig=<ID=chrY,length=59373566,assembly=hg19>\n";
print OF "##contig=<ID=chr1_gl000191_random,length=106433,assembly=hg19>\n";
print OF "##contig=<ID=chr1_gl000192_random,length=547496,assembly=hg19>\n";
print OF "##contig=<ID=chr4_ctg9_hap1,length=590426,assembly=hg19>\n";
print OF "##contig=<ID=chr4_gl000193_random,length=189789,assembly=hg19>\n";
print OF "##contig=<ID=chr4_gl000194_random,length=191469,assembly=hg19>\n";
print OF "##contig=<ID=chr6_apd_hap1,length=4622290,assembly=hg19>\n";
print OF "##contig=<ID=chr6_cox_hap2,length=4795371,assembly=hg19>\n";
print OF "##contig=<ID=chr6_dbb_hap3,length=4610396,assembly=hg19>\n";
print OF "##contig=<ID=chr6_mann_hap4,length=4683263,assembly=hg19>\n";
print OF "##contig=<ID=chr6_mcf_hap5,length=4833398,assembly=hg19>\n";
print OF "##contig=<ID=chr6_qbl_hap6,length=4611984,assembly=hg19>\n";
print OF "##contig=<ID=chr6_ssto_hap7,length=4928567,assembly=hg19>\n";
print OF "##contig=<ID=chr7_gl000195_random,length=182896,assembly=hg19>\n";
print OF "##contig=<ID=chr8_gl000196_random,length=38914,assembly=hg19>\n";
print OF "##contig=<ID=chr8_gl000197_random,length=37175,assembly=hg19>\n";
print OF "##contig=<ID=chr9_gl000198_random,length=90085,assembly=hg19>\n";
print OF "##contig=<ID=chr9_gl000199_random,length=169874,assembly=hg19>\n";
print OF "##contig=<ID=chr9_gl000200_random,length=187035,assembly=hg19>\n";
print OF "##contig=<ID=chr9_gl000201_random,length=36148,assembly=hg19>\n";
print OF "##contig=<ID=chr11_gl000202_random,length=40103,assembly=hg19>\n";
print OF "##contig=<ID=chr17_ctg5_hap1,length=1680828,assembly=hg19>\n";
print OF "##contig=<ID=chr17_gl000203_random,length=37498,assembly=hg19>\n";
print OF "##contig=<ID=chr17_gl000204_random,length=81310,assembly=hg19>\n";
print OF "##contig=<ID=chr17_gl000205_random,length=174588,assembly=hg19>\n";
print OF "##contig=<ID=chr17_gl000206_random,length=41001,assembly=hg19>\n";
print OF "##contig=<ID=chr18_gl000207_random,length=4262,assembly=hg19>\n";
print OF "##contig=<ID=chr19_gl000208_random,length=92689,assembly=hg19>\n";
print OF "##contig=<ID=chr19_gl000209_random,length=159169,assembly=hg19>\n";
print OF "##contig=<ID=chr21_gl000210_random,length=27682,assembly=hg19>\n";
print OF "##contig=<ID=chrUn_gl000211,length=166566,assembly=hg19>\n";
print OF "##contig=<ID=chrUn_gl000212,length=186858,assembly=hg19>\n";
print OF "##contig=<ID=chrUn_gl000213,length=164239,assembly=hg19>\n";
print OF "##contig=<ID=chrUn_gl000214,length=137718,assembly=hg19>\n";
print OF "##contig=<ID=chrUn_gl000215,length=172545,assembly=hg19>\n";
print OF "##contig=<ID=chrUn_gl000216,length=172294,assembly=hg19>\n";
print OF "##contig=<ID=chrUn_gl000217,length=172149,assembly=hg19>\n";
print OF "##contig=<ID=chrUn_gl000218,length=161147,assembly=hg19>\n";
print OF "##contig=<ID=chrUn_gl000219,length=179198,assembly=hg19>\n";
print OF "##contig=<ID=chrUn_gl000220,length=161802,assembly=hg19>\n";
print OF "##contig=<ID=chrUn_gl000221,length=155397,assembly=hg19>\n";
print OF "##contig=<ID=chrUn_gl000222,length=186861,assembly=hg19>\n";
print OF "##contig=<ID=chrUn_gl000223,length=180455,assembly=hg19>\n";
print OF "##contig=<ID=chrUn_gl000224,length=179693,assembly=hg19>\n";
print OF "##contig=<ID=chrUn_gl000225,length=211173,assembly=hg19>\n";
print OF "##contig=<ID=chrUn_gl000226,length=15008,assembly=hg19>\n";
print OF "##contig=<ID=chrUn_gl000227,length=128374,assembly=hg19>\n";
print OF "##contig=<ID=chrUn_gl000228,length=129120,assembly=hg19>\n";
print OF "##contig=<ID=chrUn_gl000229,length=19913,assembly=hg19>\n";
print OF "##contig=<ID=chrUn_gl000230,length=43691,assembly=hg19>\n";
print OF "##contig=<ID=chrUn_gl000231,length=27386,assembly=hg19>\n";
print OF "##contig=<ID=chrUn_gl000232,length=40652,assembly=hg19>\n";
print OF "##contig=<ID=chrUn_gl000233,length=45941,assembly=hg19>\n";
print OF "##contig=<ID=chrUn_gl000234,length=40531,assembly=hg19>\n";
print OF "##contig=<ID=chrUn_gl000235,length=34474,assembly=hg19>\n";
print OF "##contig=<ID=chrUn_gl000236,length=41934,assembly=hg19>\n";
print OF "##contig=<ID=chrUn_gl000237,length=45867,assembly=hg19>\n";
print OF "##contig=<ID=chrUn_gl000238,length=39939,assembly=hg19>\n";
print OF "##contig=<ID=chrUn_gl000239,length=33824,assembly=hg19>\n";
print OF "##contig=<ID=chrUn_gl000240,length=41933,assembly=hg19>\n";
print OF "##contig=<ID=chrUn_gl000241,length=42152,assembly=hg19>\n";
print OF "##contig=<ID=chrUn_gl000242,length=43523,assembly=hg19>\n";
print OF "##contig=<ID=chrUn_gl000243,length=43341,assembly=hg19>\n";
print OF "##contig=<ID=chrUn_gl000244,length=39929,assembly=hg19>\n";
print OF "##contig=<ID=chrUn_gl000245,length=36651,assembly=hg19>\n";
print OF "##contig=<ID=chrUn_gl000246,length=38154,assembly=hg19>\n";
print OF "##contig=<ID=chrUn_gl000247,length=36422,assembly=hg19>\n";
print OF "##contig=<ID=chrUn_gl000248,length=39786,assembly=hg19>\n";
print OF "##contig=<ID=chrUn_gl000249,length=38502,assembly=hg19>\n";
print OF "##INFO=<ID=ALCINDEX,Number=1,Type=String,Description=\"ALAMUT class index type\">\n";
print OF "##INFO=<ID=ALCVAL,Number=1,Type=String,Description=\"ALAMUT class value\">\n";
print OF "##INFO=<ID=ALAMUT,Number=1,Type=Integer,Description=\"Variant present in local Alamut mut repository\">\n";
print OF "##INFO=<ID=ALMUTID,Number=1,Type=String,Description=\"ALAMUT mutation id\">\n";
print OF "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\n";
foreach $chr (sort { $a <=> $b } keys(%outdat)) {
    foreach $posi (sort { $a <=> $b } keys(%{$outdat{$chr}})) {
	print OF "$outdat{$chr}{$posi}\n";
    }
}
close(OF);

system("vt normalize $of -r $hg19 -o $of2");

system("rm $of");

sub rc {
  my $dna = shift;
  my $revcomp = reverse($dna);
  $revcomp =~ tr/ACGTacgt/TGCAtgca/;
  return $revcomp;
}
