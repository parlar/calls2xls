#!/usr/bin/perl -w
#

use IPC::Run;
use List::Util qw(max min);
use Spreadsheet::WriteExcel;
use YAML::XS 'LoadFile';
use Time::localtime;
use Getopt::Long;

$o{'in_silico'} = 0;
$o{'help'} = 0;
$o{'version'} = "v0.42";

## Read configs
open (FH, "</home/data_in/calls2xls/calls2xls.cfg");
while ( $line = <FH> ) {
    if($line !~ /#/ && $line =~ /\S/){
	@tmp = split(/=/, $line);
	foreach $row (@tmp) {
	    $row =~ s/\s//g;
	}
	$c{$tmp[0]}=$tmp[1];
    }
}
close(FH);

GetOptions (\%o, 'in_silico', "help", "data_root=s", "old_sample=s", "new_sample=s", "disease_indication=s", "disease_gene_file=s");
if(!exists($o{'data_root'})){
    &helpmess();
    exit 0;
}
if($o{'help'}){
    &helpmess();
    exit 0;
}
if($o{'help'}){
    &helpmess();
    exit 0;
}
if($o{'in_silico'} && ((!exists($o{'old_sample'}) || !exists($o{'new_sample'})) || !exists($o{'disease_indication'})) || !exists($o{'disease_gene_file'})){
    print "\nin silico run requires that options --old_sample --new_sample --disease_indication --disease_gene_file are set.\n";
    &helpmess();
    exit 0;
}


sub helpmess{
    print "\nversion: $o{'version'}\n";
    print "\nusage: calls2xls.pl [options]\n\n";
    print "  ---data_root=path (required)\n";
    print "    path to data root folder containing raw data and variant calls \n";
    print "    from bcbio-nextgen.\n";    
    print "  --in_silico\n";
    print "    run in silico.\n";
    print "  --old_sample=sample number, e.g. 15-1234\n";
    print "    with --in_silico: old sample number to transfer calling data from.\n";
    print "  --new_sample=sample_number, e.g. 15-1235\n";
    print "    with --in_silico: new sample number to transfer calling data to.\n";
    print "  --disease_indication=indication, e.g. HCM\n";
    print "    with --in_silico: new disease indication.\n";    
    print "  --disease_gene_file=file\n";
    print "    disease gene file defining gene sets for disease indications.\n";

    chomp(@disgenefiles = <$c{'designAndGeneListsRootPath'}/*.txt>);
    print "\n";
    if(exists($disgenefiles[0])) {
	print "disease gene files and disease indications:\n";
    }
    foreach $file (@disgenefiles){
        @tmp = split(/\//, $file);
        print "  $tmp[-1]:\n";
        open (FH, "<$file");
	@diseases = ();
        while ( $line = <FH> ) {
            if($line =~ /^>(\S+)/){
                push @diseases, $1;
            }
        }
	close(FH);
	$dstr = "";
        foreach $row (@diseases){
	    $dstr .= "$row ";
	    if(length($dstr) > 60) {
		print "     $dstr\n";
		$dstr = "";
	    }
	}
	print "     $dstr\n";
    }
    print "\n";
} 


chomp($miseqroot = $o{'data_root'});

#colors
%clr=(
    'red' => 10,
    'pink' => 45,
    'yellow' => 13,
    'lgreen'=> 42,
    'green' => 11
    );

$rundate = (split /_/, $miseqroot)[0];

$tm = localtime;
$nowtime = sprintf("%04d.%02d.%02d-%02d.%02d.%02d", $tm->year+1900, ($tm->mon)+1, $tm->mday, $tm->hour, $tm->min, $tm->sec );

$finalpath = $miseqroot . "/calling/final";

$samplesheetfile = $miseqroot . "/SampleSheet.csv";
open(FH, "<$samplesheetfile");
chomp(@samplesheet = <FH>);
close(FH);

foreach $row (@samplesheet){
    $row =~ s/\s+$//g;
    if($row =~ /^\[(.*?)\]/) {
	$field = $1;
    }
    elsif($row =~ /\S/ && $row =~ /(.*)/){
	@tmp = split(/,/, $1);
	if($#tmp == 1){
#	    print "1: @tmp\n";
	    $ssheet{$field}{$tmp[0]} = $tmp[1];
	}
	elsif($#tmp == 0){
#	    print "0: @tmp\n";
	    $ssheet{$field} = $tmp[0];
	}
	elsif($#tmp > 1){
	    if($tmp[0] eq "Sample_ID"){
		@header = @tmp;
	    }
	    else{
		for ($i = 0; $i <= $#tmp; $i++){
		    if($header[$i] eq "Description"){
			@desc_header = ("design","gender","inquiry","genelist");
			@tmp2 = split(/\$/, $tmp[$i]);
			for ($j = 0; $j <= $#tmp2; $j++){
			    print "$tmp2[$j]\n";
			    $ssheet{$field}{$header[$i]}{$tmp[0]}{$desc_header[$j]} = $tmp2[$j];
			}
		    }
		}
	    }
	}
    }
}

#move to sampleC
foreach $sample ( sort { $b cmp $a } keys(%{$ssheet{'Data'}{'Description'}})) {
    $sampleC{$sample}{sex} = $ssheet{'Data'}{'Description'}{$sample}{'gender'};
    $sampleC{$sample}{regions} = $ssheet{'Data'}{'Description'}{$sample}{'design'};
    $sampleC{$sample}{disease} = $ssheet{'Data'}{'Description'}{$sample}{'inquiry'};
    $sampleC{$sample}{genelist} = $ssheet{'Data'}{'Description'}{$sample}{'genelist'};
    $sampleC{$sample}{disease} =~ s/\s//g;
}

if($o{'in_silico'}){
    foreach $field ( sort { $b cmp $a } keys(%{$sampleC{$oldsample}})) {
	$sampleC{$newsample}{$field} = $sampleC{$oldsample}{$field};
    }
    
    $sampleC{$newsample}{genelist} = $c{'-disease_gene_file'};
    $sampleC{$newsample}{disease} = $c{'disease_indication'};
    
    foreach $sample ( sort { $b cmp $a } keys(%sampleC)) {
	if($sample ne $newsample){
	    delete $sampleC{$sample};
	}
    }
    
    $newsamppath = $finalpath . "/" . $newsample;
    
    if(!-d $newsamppath){
	system("mkdir $newsamppath");
    }
#    system("rm $finalpath/$newsample/*");
#    system("cp -r $finalpath/$oldsample/* $finalpath/$newsample");
#    system("rename s/$oldsample-/$newsample-/g $newsamppath/*");
#    system("rm $newsamppath/sample*");
#    system("rm -r $newsamppath/$oldsample");
    
    print "rm $finalpath/$newsample/*\n";
    print "cp -r $finalpath/$oldsample/* $finalpath/$newsample\n";
    print "rename s/$oldsample-/$newsample-/g $newsamppath/*\n";
    print "rm $newsamppath/sample*\n";
    print "rm -r $newsamppath/$oldsample\n";
}
=pod
$projectpath = <$finalpath/*project*>;
foreach $sample ( sort { $b cmp $a } keys(%sampleC)) {
    $dir = $finalpath . "/" . $sample;
    if($o{'in_silico'}){
	if(($dir !~ /project/ && -d $dir) && $dir =~ /$newsample/) {
	    push @sampledirs, $dir;
	}
    }
    else{
	if($dir !~ /project/ && -d $dir) {
	    push @sampledirs, $dir;
	}
    }
}

system("mut2vcf_v3.pl $projectpath");

#get calldate
$summary_file = $projectpath . "/project-summary.yaml";
$config = LoadFile( $summary_file );
$dt = $config -> {date};
$dt =~ /\d\d(\d\d)-(\d\d)-(\d\d)\s/;
$calldate = $1 . $2 . $3;

#get programs
$programs = $projectpath . "/programs.txt";
open(FH, "<$programs") or die $!;
chomp(@PROGS = <FH>);
close(FH) or die $!;

#get calling.yaml
$bcbio_config = $miseqroot . $c{'bcbio_config'};
open(FH, "<$bcbio_config") or die $!;
chomp(@cf = <FH>);
close(FH) or die $!;
$cs = "";

for ($i = 0;  $i <= $#cf; $i++){
    if($cf[$i] =~ /  - files: /) {
	$cf[$i+1] =~ /    description: (\S+)/;
	$cs = $1;
	$configHash{$cs}{$i} = $cf[$i];
    }
    elsif($cs =~ /\S/) {
	$configHash{$cs}{$i} = $cf[$i];
    }
}

#get diseasegenesdef info
foreach $sample ( sort { $b cmp $a } keys(%sampleC)) {
    $diseasegenesdef = $c{'designAndGeneListsRootPath'} . "/" . $sampleC{$sample}{genelist};
    $diseasegenesdef_file = $sampleC{$sample}{genelist};
    open(FH0, "<$diseasegenesdef");
    chomp(@GeneL = <FH0>);
    close(FH0);

    foreach $row (@GeneL) {
	$row =~ s/\s+//g;
	if($row =~ />(\S+)/) {
	    $dise = $1;
	}
	elsif($row =~ /(\S+)/){
	    $disease2genes{$sample}{$dise}{$1} = 1;
#	    print "$sample $dise $1\n"; 
	}
    }
}

foreach $sample (sort { $a cmp $b } keys(%sampleC)) {
    $ccnf = $c{'designAndGeneListsRootPath'} . "/" .  $sampleC{$sample}{regions};
    open(FH0, "<$ccnf");
    chomp(@REG = <FH0>);
    close(FH0);
    foreach $row (@REG) {
	@tmp = split(/\t/, $row);
	$tmp[0] =~ s/\s//g; #chr
	$tmp[1] =~ s/\s//g; #p1
	$tmp[2] =~ s/\s//g; #p2
	$tmp[3] =~ s/\s//g; #gene
	
	if($tmp[0] !~ /chr/){
	    if($tmp[0] eq "MT") {
		$tmp[0] = "M";
	    }
	    $chr = "chr". $tmp[0];
	}
	else{
	    $chr = $tmp[0];
	}
	
	if(exists($disease2genes{$sample}{$sampleC{$sample}{disease}}{$tmp[3]})){
	    for($i = $tmp[1] ; $i<= $tmp[2]; $i++) {
		$onTarget{$sample}{$chr}{$i} = $tmp[3];
	    }
	}
    }
}

#$c{'coreRoiRootPath = $c{'designAndGeneListsRootPath'} . "/core_rois";

@cores = <$c{'coreRoiRootPath'}/*.bed>;
foreach $corebed (@cores) {
#    print "$corebed\n";
    @tmp1 = split(/\//, $corebed);
    @tmp2 = split(/\./, $tmp1[-1]);
#    print "$tmp2[0] $corebed\n";
    $coreH{$tmp2[0]} = $corebed;
}

system("addLocAfDb.pl $miseqroot");


$command0 = "comm." . $miseqroot . ".0.txt";
$command1 = "comm." . $miseqroot . ".1.txt";
$command2 = "comm." . $miseqroot . ".2.txt";
$command3 = "comm." . $miseqroot . ".3.txt";
$command4 = "comm." . $miseqroot . ".4.txt";

$command0 =~ s/\///g;
$command1 =~ s/\///g;
$command2 =~ s/\///g;
$command3 =~ s/\///g;
$command4 =~ s/\///g;

open(FH0, ">$command0");
open(FH1, ">$command1");
open(FH2, ">$command2");
open(FH3, ">$command3");
open(FH4, ">$command4");


foreach $sampledir (@sampledirs) {
    print "$sampledir\n";
    $sample = (split /\//, $sampledir )[-1];
    @tmp = split(/\//, $sampledir);
    $bam = $sampledir . "/" . $tmp[-1] . "-ready.bam";

    $reg =  $c{'designAndGeneListsRootPath'} . $sampleC{$sample}{regions};

    print "$reg\n";

    print "$sampledir\n";
    print "$bam\n";

    $qcdir = $sampledir . "/qc";

    if ( !-d $qcdir ) {
	system("mkdir $qcdir");
    }

    $destdir = $sampledir . "/qc/teqc";
    if ( !-d $destdir ) {
	system("mkdir $destdir");
    }

    $outfile = $sampledir . "/" . "sample.coverage.bed";
    $outfile2 = $sampledir . "/" . "sample.perbase.coverage.bed";
    $outfile_le30merge = $sampledir . "/" . "sample.panel.bg.coverage.le30merge.bed";
    $coreregs = $sampledir . "/" . "sample.core.bg.coverage.le30merge.bed";

    $bam = $sampledir . "/" . $sample . "-ready.bam";

    
    print "bedtools genomecov -ibam $bam -bga -g $c{'genome'} -dz \n";
    print "bedtools intersect -wb -a stdin -b $reg\n";
    print "cut -f 1-4,8 | sed \'s/.*/&\n";
    print "$sample\n";
    print "$rundate\n";
    print "$calldate\n";
    print "$sampleC{$sample}{sex}\n";
    print "$outfile\n";
    print FH0 "bedtools coverage -d -abam $bam -b $reg > $outfile2\n";
    print FH1 "bedtools genomecov -ibam $bam -bga -g $c{'genome'} | bedtools intersect -wb -a stdin -b $reg | cut -f 1-4,8 | sed \'s/.*/&\t$sample\t$rundate\t$calldate\t$sampleC{$sample}{sex}/\' >$outfile\n";
    print FH2 "filter_coverage.pl $outfile | sort -k1,1 -k2,2n | bedtools merge -i stdin -c 5,6,7,8 -o distinct,distinct,distinct,distinct > $outfile_le30merge\n";
    print FH3 "samtools index $bam\n";
    if(exists($coreH{$sampleC{$sample}{disease}})){
	print FH3 "bedtools intersect -a $outfile_le30merge -b $coreH{$sampleC{$sample}{disease}} >$coreregs\n";
    }
    else{
	open(FLOC, ">$coreregs");
	print FLOC  "Inga core-regioner definierade";
	close(FLOC);
    }
    print FH4 "run_TEQC.pl $bam $reg $sample $destdir\n";
}

close(FH0);
close(FH1);
close(FH2);
close(FH3);
close(FH4);

system("parallel -j 10 -- < $command0");
system("parallel -j 10 -- < $command1");
system("parallel -j 10 -- < $command2");
system("parallel -j 10 -- < $command3");
system("parallel -j 10 -- < $command4");

system("rm $command0");
system("rm $command1");
system("rm $command2");
system("rm $command3");
system("rm $command4");

#unpack ensemble vcf
$ensvcfgz = $projectpath . "/batch1-ensemble.vcf.gz";
$ensvcf = $projectpath . "/batch1-ensemble.vcf";
$ensdecompvcf = $projectpath . "/batch1-ensemble.d.vcf";
$ann1vcf = $projectpath . "/batch1-ensemble.d.1kg.esp.exac.vcf";
$ann2vcf = $projectpath . "/batch1-ensemble.d.1kg.esp.exac.ala.vcf";
$ann3vcf = $projectpath . "/batch1-ensemble.d.1kg.esp.exac.ala.frq.vcf";
$ann4vcf = $projectpath . "/batch1-ensemble.d.1kg.esp.exac.ala.frq.ann.vcf";
$ann4avcf = $projectpath . "/batch1-ensemble.d.1kg.esp.exac.ala.frq.ann." . $newsample . ".vcf";

system("pigz -f -d -k $ensvcfgz");
system("vt decompose $ensvcf -o $ensdecompvcf");
system("snpsift_annotate_esp_1kg_exac.pl $ensdecompvcf");
system("snpsift_annotate_alamut.pl $ann1vcf");
system("getfreq.pl $ann2vcf > $ann3vcf");

system("java -Xmx4g -jar $c{'snpeff_path'} hg19 -c /home/data_in/databases/snpeff/snpEff.config $ann3vcf > $ann4vcf");

if($o{'in_silico'}){
    open(FH, "<$ann4vcf") or die $!;
    open(FH2, ">$ann4avcf") or die $!;
    while ($line = <FH>) {
	if($line =~ /^\S/) {
	    chomp($line);
	    if($line =~ /^#CHROM/){
		@tmp = split(/\t/, $line);
		for ($i = 0; $i <= $#tmp; $i++){
		    if($i >= 9){
			if ($tmp[$i] =~ /^(\d\d)-(\d\d\d)$/){
			    $tmp[$i] = $1 . "-0" . $2;
			}
		    }
		    
		    if($tmp[$i] eq $oldsample){
			$s2c{$newsample} = $i;
		    }
		}
		$line .= "\t" . $newsample;
		print FH2 "$line\n";
	    }
	    elsif($line =~ /^#/){
		print FH2 "$line\n";
	    }
	    else{
		@tmp = split(/\t/, $line);
		push @tmp, $tmp[$s2c{$newsample}];
		$line = join("\t", @tmp);
		print FH2 "$line\n";
	    }
	}
    }
    close(FH);
    close(FH2);
    $ann4vcf = $ann4avcf;
}




open(FH, "<$ann4vcf") or die $!;
chomp(@gout = <FH>);
close(FH) or die $!;

%sample2column = ();
%column2sample = ();
%calls = ();

$i = 0;
foreach $row (@gout) {
    if($row =~ /^#CHROM/) {
#	print "$row\n";
	@tmp = split(/\s+/, $row);
	for ($j = 8; $j <= $#tmp; $j++) {
	    $tmp[$j] =~ s/\s+//g;
	    $column2sample{$j} = $tmp[$j];
	    $sample2column{$tmp[$j]} = $j;
#	    print "$j $tmp[$j]\n";
	}
    }

    elsif($row !~ /^#/ && $row =~ /\S/) {
	$i++;
	@tmp = split(/\s+/, $row);
	$chr = $tmp[0];
	$p = $tmp[1];
	$ref = $tmp[3];
	$alt = $tmp[4];

	$exac_af = 0;
	$esp_af_ea = 0;
	$esp_af_aa = 0;
	$esp_af_all = 0;
	$onekg_af_all = 0;
	$onekg_af_afr = 0;
	$onekg_af_amr = 0;
	$onekg_af_eas = 0;
	$onekg_af_eur = 0;
	$onekg_af_sas = 0;

	$alamut = ".";
	$alcindex = ".";
	$alcval = ".";
	$alcreate = ".";
	$alupdate = ".";
	
	%annhash = ();

	@INFO = split(/;/, $tmp[7]);
#	print "@INFO\n";
	%infohash = ();
	foreach $in (@INFO) {
	    ($index, $val) = split(/=/, $in);
	    $infohash{$index} = $val;
#	    print "$index  $infohash{$index}\n";
	}
	if(exists($infohash{'EXAC_AF'})) {
	    $exac_af = sprintf("%.3f", $infohash{'EXAC_AF'});
	}
	if(exists($infohash{'ESP_MAF'})) {
	    @esps = split(/,/, $infohash{'ESP_MAF'});
	    foreach $x (@esps) { 
		$x = sprintf("%.1f", $x ); 
	    }
	    $esp_af_ea = $esps[0];
	    $esp_af_aa = $esps[1];
	    $esp_af_all = $esps[2];
	}
	if(exists($infohash{'1KG_AF'})) {
	    $onekg_af_all = sprintf("%.1f", $infohash{'1KG_AF'} * 100); 
	}
	if(exists($infohash{'1KG_AFR_AF'})) {
	    $onekg_af_afr = sprintf("%.1f", $infohash{'1KG_AFR_AF'}* 100); 
	}
	if(exists($infohash{'1KG_AMR_AF'})) {
	    $onekg_af_amr = sprintf("%.1f", $infohash{'1KG_AMR_AF'}* 100); 
	}
	if(exists($infohash{'1KG_EAS_AF'})) {
	    $onekg_af_eas = sprintf("%.1f", $infohash{'1KG_EAS_AF'}* 100); 
	}
	if(exists($infohash{'1KG_EUR_AF'})) {
	    $onekg_af_eur = sprintf("%.1f", $infohash{'1KG_EUR_AF'}* 100); 
	}
	if(exists($infohash{'1KG_SAS_AF'})) {
	    $onekg_af_sas = sprintf("%.1f", $infohash{'1KG_SAS_AF'}* 100); 
	}
	if(exists($infohash{'ALAMUT'})) {
	    $alamut = $infohash{'ALAMUT'}; 
	}
	if(exists($infohash{'ALCINDEX'})) {
	    $alcindex = $infohash{'ALCINDEX'};
	}
	if(exists($infohash{'ALCVAL'})) {
	    $alcval = $infohash{'ALCVAL'}; 
	}
	if(exists($infohash{'ALCREATE'})) {
	    $alcreate = $infohash{'ALCREATE'}; 
	}
	if(exists($infohash{'ALUPDATE'})) {
	    $alupdate = $infohash{'ALUPDATE'}; 
	}
	if(exists($infohash{'LOCFREQ'})) {
	    $locfreq = sprintf("%.1f", $infohash{'LOCFREQ'}* 100);  
#	    print "$locfreq\n";
	}
	if(exists($infohash{'LOCALS'})) {
	    $locals = $infohash{'LOCALS'};
#	    print "$locals\n";	    
	}
	
	if(exists($infohash{'ANN'})){
	    @anns = split(/,/, $infohash{'ANN'});
	    for($ano = 0; $ano <= $#anns; $ano++) {
		$anns[$ano] =~ s/\|\|/\| \|/g;
		@anntmp = split(/\|/, $anns[$ano]);
		if($anntmp[0] eq $alt) {
		    if(exists($anntmp[1])){
			$annhash{$ano}{type} = $anntmp[1];
		    }
		    else{
			$annhash{$ano}{type} = " ";
		    }
		    if(exists($anntmp[2])){
			$annhash{$ano}{impact} = $anntmp[2];
		    }
		    else{
			$annhash{$ano}{impact} = " ";
		    }
		    if(exists($anntmp[4])){
			$annhash{$ano}{genename} = $anntmp[4];
		    }
		    else{
			$annhash{$ano}{genename} = " ";
		    }
		    if(exists($anntmp[5])){
			$annhash{$ano}{fid} = $anntmp[5];
		    }
		    else{
			$annhash{$ano}{fid} = " ";
		    }
		    if(exists($anntmp[6])){
			$annhash{$ano}{btype} = $anntmp[6];
		    }
		    else{
			$annhash{$ano}{btype} = " ";
		    }
		    if(exists($anntmp[8])){
			$annhash{$ano}{hgvsc} = $anntmp[9];
		    }
		    else{
			$annhash{$ano}{hgvsc} = " ";
		    }
		    if(exists($anntmp[9])){
			$annhash{$ano}{hgvsp} = $anntmp[10];
		    }
		    else{
			$annhash{$ano}{hgvsp} = " ";
		    }		    
		}
	    }
	}
	
	@format = split(/:/, $tmp[8]);
	
	$dp_no = "NA";
	$ao_no = "NA";
	
	foreach ($j = 0; $j <= $#format; $j++){
	    if($format[$j] eq "NR") {
		$dp_no = $j;
	    }
	    elsif($format[$j] eq "DP") {
		$dp_no = $j;
	    }
	    if($format[$j] eq "AO") {
		$ao_no = $j;
	    }
	    elsif($format[$j] eq "NV") {
		$ao_no = $j;
	    }
	}

	foreach $sample (sort { $a cmp $b } keys(%sample2column)) {
	    @tmp2 = split(/:/, $tmp[$sample2column{$sample}]);
	    if($tmp2[0] =~ /[1-9]/ && exists($onTarget{$sample}{$chr}{$p})) {
		foreach $ano (sort { $a <=> $b } keys(%annhash)) {
		    if(!exists($calls{$sample})){
			$outno = 0;
		    }
		    else{
			$outno = max keys(%{$calls{$sample}}) + 1;
		    }
#		    print "$sample $outno\n";
		    $calls{$sample}{$outno}{'gt'} =  $tmp2[0];
		    if($dp_no ne "NA") {
			$calls{$sample}{$outno}{'dp'} =  $tmp2[$dp_no];
		    }
		    else{
			$calls{$sample}{$outno}{'dp'} =  " ";
		    }
		    if($ao_no ne "NA") {
			$calls{$sample}{$outno}{'ao'} =  $tmp2[$ao_no];
		    }
		    else{
			$calls{$sample}{$outno}{'ao'} =  " ";
		    }
		    $calls{$sample}{$outno}{'varno'} =  $varcount{$sample}{$chr}{$p}{$ref}{$alt};
		    $calls{$sample}{$outno}{'chr'} = $chr;
		    $calls{$sample}{$outno}{'p'} = $p;
		    $calls{$sample}{$outno}{'ref'} = $ref;
		    $calls{$sample}{$outno}{'alt'} = $alt;
		    $calls{$sample}{$outno}{'exac_af'} = $exac_af;
		    $calls{$sample}{$outno}{'esp_af_ea'} = $esp_af_ea;
		    $calls{$sample}{$outno}{'esp_af_aa'} = $esp_af_aa;
		    $calls{$sample}{$outno}{'1kg_af_all'} = $onekg_af_all;
		    $calls{$sample}{$outno}{'1kg_af_afr'} = $onekg_af_afr;
		    $calls{$sample}{$outno}{'1kg_af_amr'} = $onekg_af_amr;
		    $calls{$sample}{$outno}{'1kg_af_eas'} = $onekg_af_eas;
		    $calls{$sample}{$outno}{'1kg_af_eur'} = $onekg_af_eur;
		    $calls{$sample}{$outno}{'1kg_af_sas'} = $onekg_af_sas;
		    $calls{$sample}{$outno}{'type'} = $annhash{$ano}{type};
		    $calls{$sample}{$outno}{'impact'} = $annhash{$ano}{impact};
		    $calls{$sample}{$outno}{'genename'} = $annhash{$ano}{genename};
		    $calls{$sample}{$outno}{'fid'} = $annhash{$ano}{fid};
		    $calls{$sample}{$outno}{'btype'} = $annhash{$ano}{btype};
		    $calls{$sample}{$outno}{'hgvsc'} = $annhash{$ano}{hgvsc};
		    $calls{$sample}{$outno}{'hgvsp'} = $annhash{$ano}{hgvsp};
		    $calls{$sample}{$outno}{'alamut'} = $alamut;
		    $calls{$sample}{$outno}{'alcindex'} = $alcindex;
		    $calls{$sample}{$outno}{'alcval'} = $alcval;
		    $calls{$sample}{$outno}{'alcreated'} = $alcreate;
		    $calls{$sample}{$outno}{'alupdated'} = $alupdate;
		    $calls{$sample}{$outno}{'locfreq'} = $locfreq;
		    $calls{$sample}{$outno}{'locals'} = $locals;

		}
	    }
	}
    }
}


foreach $sampledir (@sampledirs) {
    $sample = (split /\//, $sampledir )[-1];
    $bam_file = $sample . "-ready.bam";
    $bam = $sampledir . "/" . $sample . "-ready.bam";

    %genesHash = ();
    
    foreach $gene (sort { $a cmp $b } keys(%{$disease2genes{$sample}{$sampleC{$sample}{disease}}})) {
	$gene =~ s/\s+//g;
	$genesHash{$gene} = 1;
	print "$sample $sampleC{$sample}{disease} $gene\n";
    }

    
    
    $sampleout_folder = $sampledir . "/" . $sample;
    $sampleout_rundate_folder = $sampleout_folder . "/" . $rundate . "." . $nowtime;
    
    $winpath = $sample . "\\" . $rundate . "." . $nowtime;
    
    if(!-d $sampleout_folder) {
	system("mkdir $sampleout_folder");
    }
    else{
    }
    if(!-d $sampleout_rundate_folder) {
	system("mkdir $sampleout_rundate_folder");
    }
    else{
    }
    
    $xls_name = $sampleout_rundate_folder . "/" . $sample . "_calldata.xls";
    

    
    my $wb = Spreadsheet::WriteExcel->new($xls_name);
    $ws = $wb->add_worksheet('calldata');
    
    $format = $wb->add_format();

    $ws->set_column('A:A', 5);
    $ws->set_column('B:B', 5);
    $ws->set_column('C:C', 11);
    $ws->set_column('D:D', 4);
    $ws->set_column('E:E', 4);
    $ws->set_column('F:F', 5);
    $ws->set_column('G:G', 5);
    $ws->set_column('H:H', 5);
    $ws->set_column('I:I', 5);
    $ws->set_column('J:J', 5);
    $ws->set_column('K:K', 17);
    $ws->set_column('L:L', 13);
    $ws->set_column('O:O', 15);
    $ws->set_column('P:P', 13);
    $ws->set_column('Q:Q', 13);
    $ws->set_column('R:R', 4);
    $ws->set_column('S:S', 4);
    $ws->set_column('T:T', 4);
    $ws->set_column('U:U', 6);
    $ws->set_column('V:V', 7);
    $ws->set_column('X:X', 11);
    $ws->set_column('Y:Y', 11);
    
    $ws->split_panes(12.75, 0, 1, 0);
    
    $five = $wb->add_format(bg_color=>$clr{'red'},border   => 1);
    $four = $wb->add_format(bg_color=>$clr{'pink'},border   => 1);
    $three = $wb->add_format(bg_color=>$clr{'yellow'},border   => 1);
    $two = $wb->add_format(bg_color=>$clr{'lgreen'},border   => 1);
    $one = $wb->add_format(bg_color=>$clr{'green'},border   => 1);
    $zero = $wb->add_format(bg_color=> 9,border   => 1, );
    $zerobold = $wb->add_format(bg_color=> 9,border   => 1, );
    $zerobold->set_bold();
    

    $i = 0;

    @ROWS =();
    @ONEROWS = ();
    @TWOROWS = ();
    @THREEROWS = ();
    @FOURROWS = ();
    @FIVEROWS = ();
    
    foreach $n (sort { $a <=> $b } keys(%{$calls{$sample}})) {
#	print "$sample $calls{$sample}{$n}{'chr'} $calls{$sample}{$n}{'p'} $onTarget{$sample}{$calls{$sample}{$n}{'chr'}}{$calls{$sample}{$n}{'p'}} \n";
	@outrow = ();
	$targetGene = $onTarget{$sample}{$calls{$sample}{$n}{'chr'}}{$calls{$sample}{$n}{'p'}};
#	push @outrow, $calls{$sample}{$n}{'varno'};
	push @outrow, $calls{$sample}{$n}{'chr'};
	push @outrow, $calls{$sample}{$n}{'p'};
	push @outrow, $calls{$sample}{$n}{'ref'};
	push @outrow, $calls{$sample}{$n}{'alt'}; 
#	push @outrow, $calls{$sample}{$n}{'exac_af'};
	push @outrow, $calls{$sample}{$n}{'esp_af_ea'};
	push @outrow, $calls{$sample}{$n}{'esp_af_aa'};
	push @outrow, $calls{$sample}{$n}{'1kg_af_all'};
#	push @outrow, $calls{$sample}{$n}{'1kg_af_afr'}; 
#	push @outrow, $calls{$sample}{$n}{'1kg_af_amr'}; 
#	push @outrow, $calls{$sample}{$n}{'1kg_af_eas'}; 
#	push @outrow, $calls{$sample}{$n}{'1kg_af_eur'}; 
#	push @outrow, $calls{$sample}{$n}{'1kg_af_sas'};
	push @outrow, $calls{$sample}{$n}{'locfreq'}; 
	push @outrow, $calls{$sample}{$n}{'locals'};
	push @outrow, $calls{$sample}{$n}{'type'}; 
	push @outrow, $calls{$sample}{$n}{'impact'}; 
	push @outrow, $calls{$sample}{$n}{'genename'};
	push @outrow, $targetGene;
#	push @outrow, $calls{$sample}{$n}{'fid'}; 
	push @outrow, $calls{$sample}{$n}{'btype'}; 
	push @outrow, $calls{$sample}{$n}{'hgvsc'}; 
	push @outrow, $calls{$sample}{$n}{'hgvsp'}; 
	push @outrow, $calls{$sample}{$n}{'gt'};
	push @outrow, $calls{$sample}{$n}{'dp'};
	push @outrow, $calls{$sample}{$n}{'ao'};
#	push @outrow, $calls{$sample}{$n}{'alamut'}; 
	push @outrow, $calls{$sample}{$n}{'alcindex'}; 
	push @outrow, $calls{$sample}{$n}{'alcval'}; 
	push @outrow, $calls{$sample}{$n}{'alcreated'}; 
	push @outrow, $calls{$sample}{$n}{'alupdated'};

#	print "@outrow\n";

	if($calls{$sample}{$n}{'alcindex'} eq "1" && $calls{$sample}{$n}{'alcval'} eq "CMGS_VGKL_5") {
	    push @ONEROWS, join(';', @outrow);
	}
	elsif($calls{$sample}{$n}{'alcindex'} eq "2" && $calls{$sample}{$n}{'alcval'} eq "CMGS_VGKL_5") {
	    push @TWOROWS, join(';', @outrow);
	}
	elsif($calls{$sample}{$n}{'alcindex'} eq "3" && $calls{$sample}{$n}{'alcval'} eq "CMGS_VGKL_5") {
	    push @THREEROWS, join(';', @outrow);
	}
	elsif($calls{$sample}{$n}{'alcindex'} eq "4" && $calls{$sample}{$n}{'alcval'} eq "CMGS_VGKL_5") {
	    push @FOURROWS, join(';', @outrow);
	}
	elsif($calls{$sample}{$n}{'alcindex'} eq "5" && $calls{$sample}{$n}{'alcval'} eq "CMGS_VGKL_5") {
	    push @FIVEROWS, join(';', @outrow);
	}
	else{
	    push @ROWS, join(';', @outrow);
	}
	
    }

    @houtrow = ("varno", "chr","position","ref","alt","espEA","espAA","1kg", "locAF","locAC","type","impact","genename","design_gene","btype","hgvsc","hgvsp","gt","dp","ao","klass","alind","skapad","uppdaterad","mut_link", "bam_link");
    
    $varcount = 0;
    %printedvars = ();

    $ws->write_row($i,0, \@houtrow, $zerobold);
    $i+=1;
    $start = $i + 1;
    $afilt = "A" . $start;
    $ws->write_row($i,0, \@houtrow, $zerobold);
    $i+=1;
    

    foreach $r (@FIVEROWS) {
	@outrow = split(/;/, $r);
	if(!exists($printedvars{$outrow[0]}{$outrow[1]}{$outrow[2]}{$outrow[3]})){
	    $varcount++;
	    $printedvars{$outrow[0]}{$outrow[1]}{$outrow[2]}{$outrow[3]} = 1;
	}
	unshift(@outrow, $varcount);

	if($outrow[20] eq "5" && $outrow[21] eq "CMGS_VGKL_5"){
	    $format = $five;
	}
	elsif($outrow[20] eq "4" && $outrow[21] eq "CMGS_VGKL_5"){
	    $format = $four;
	}
	elsif($outrow[20] eq "3" && $outrow[21] eq "CMGS_VGKL_5"){
	    $format = $three;
	}	
	elsif($outrow[20] eq "2" && $outrow[21] eq "CMGS_VGKL_5"){
	    $format = $two;
	}
	elsif($outrow[20] eq "1" && $outrow[21] eq "CMGS_VGKL_5"){
	    $format = $one;
	}
	else{
	    $format = $zero;
	}
	$mlink = "http://localhost:10000/show?request=" . $outrow[1] . ":" . $outrow[2] . "-" . $outrow[2];
	$blink = "http://localhost:10000/show?request=BAM<$c{'winshare_root'}" . $sample . "\\" . $rundate . ".$nowtime" . "\\" . $bam_file;
	$ws->write_row($i,0, \@outrow, $format);
	$ws->write_url($i, $#outrow+1, $mlink, $format, 'mut_link');
	$ws->write_url($i, $#outrow+2, $blink, $format, 'bam_link');
	$i++;
    }
    $i+=1;
    
    foreach $r (@FOURROWS) {
	@outrow = split(/;/, $r);
	if(!exists($printedvars{$outrow[0]}{$outrow[1]}{$outrow[2]}{$outrow[3]})){
	    $varcount++;
	    $printedvars{$outrow[0]}{$outrow[1]}{$outrow[2]}{$outrow[3]} = 1;
	}
	unshift(@outrow, $varcount);

	if($outrow[20] eq "5" && $outrow[21] eq "CMGS_VGKL_5"){
	    $format = $five;
	}
	elsif($outrow[20] eq "4" && $outrow[21] eq "CMGS_VGKL_5"){
	    $format = $four;
	}
	elsif($outrow[20] eq "3" && $outrow[21] eq "CMGS_VGKL_5"){
	    $format = $three;
	}	
	elsif($outrow[20] eq "2" && $outrow[21] eq "CMGS_VGKL_5"){
	    $format = $two;
	}
	elsif($outrow[20] eq "1" && $outrow[21] eq "CMGS_VGKL_5"){
	    $format = $one;
	}
	else{
	    $format = $zero;
	}
	$mlink = "http://localhost:10000/show?request=" . $outrow[1] . ":" . $outrow[2] . "-" . $outrow[2];
	$blink = "http://localhost:10000/show?request=BAM<$c{'winshare_root'}" . $sample . "\\" . $rundate . ".$nowtime" . "\\" . $bam_file;
	$ws->write_row($i,0, \@outrow, $format);
	$ws->write_url($i, $#outrow+1, $mlink, $format, 'mut_link');
	$ws->write_url($i, $#outrow+2, $blink, $format, 'bam_link');
	$i++;
    }
    $i+=1;
    
    foreach $r (@ROWS) {
	@outrow = split(/;/, $r);
	if(!exists($printedvars{$outrow[0]}{$outrow[1]}{$outrow[2]}{$outrow[3]})){
	    $varcount++;
	    $printedvars{$outrow[0]}{$outrow[1]}{$outrow[2]}{$outrow[3]} = 1;
	}
	unshift(@outrow, $varcount);

	if($outrow[20] eq "5" && $outrow[21] eq "CMGS_VGKL_5"){
	    $format = $five;
	}
	elsif($outrow[20] eq "4" && $outrow[21] eq "CMGS_VGKL_5"){
	    $format = $four;
	}
	elsif($outrow[20] eq "3" && $outrow[21] eq "CMGS_VGKL_5"){
	    $format = $three;
	}	
	elsif($outrow[20] eq "2" && $outrow[21] eq "CMGS_VGKL_5"){
	    $format = $two;
	}
	elsif($outrow[20] eq "1" && $outrow[21] eq "CMGS_VGKL_5"){
	    $format = $one;
	}
	else{
	    $format = $zero;
	}
	$mlink = "http://localhost:10000/show?request=" . $outrow[1] . ":" . $outrow[2] . "-" . $outrow[2];
	$blink = "http://localhost:10000/show?request=BAM<$c{'winshare_root'}" . $sample . "\\" . $rundate . ".$nowtime" . "\\" . $bam_file;
	$ws->write_row($i,0, \@outrow, $format);
	$ws->write_url($i, $#outrow+1, $mlink, $format, 'mut_link');
	$ws->write_url($i, $#outrow+2, $blink, $format, 'bam_link');
	$i++;
    }
    $i+=1;
    
    foreach $r (@THREEROWS) {
	@outrow = split(/;/, $r);
	if(!exists($printedvars{$outrow[0]}{$outrow[1]}{$outrow[2]}{$outrow[3]})){
	    $varcount++;
	    $printedvars{$outrow[0]}{$outrow[1]}{$outrow[2]}{$outrow[3]} = 1;
	}
	unshift(@outrow, $varcount);

	if($outrow[20] eq "5" && $outrow[21] eq "CMGS_VGKL_5"){
	    $format = $five;
	}
	elsif($outrow[20] eq "4" && $outrow[21] eq "CMGS_VGKL_5"){
	    $format = $four;
	}
	elsif($outrow[20] eq "3" && $outrow[21] eq "CMGS_VGKL_5"){
	    $format = $three;
	}	
	elsif($outrow[20] eq "2" && $outrow[21] eq "CMGS_VGKL_5"){
	    $format = $two;
	}
	elsif($outrow[20] eq "1" && $outrow[21] eq "CMGS_VGKL_5"){
	    $format = $one;
	}
	else{
	    $format = $zero;
	}
	$mlink = "http://localhost:10000/show?request=" . $outrow[1] . ":" . $outrow[2] . "-" . $outrow[2];
	$blink = "http://localhost:10000/show?request=BAM<$c{'winshare_root'}" . $sample . "\\" . $rundate . ".$nowtime" . "\\" . $bam_file;
	$ws->write_row($i,0, \@outrow, $format);
	$ws->write_url($i, $#outrow+1, $mlink, $format, 'mut_link');
	$ws->write_url($i, $#outrow+2, $blink, $format, 'bam_link');
	$i++;
    }
    $i+=1;

    foreach $r (@TWOROWS) {
	@outrow = split(/;/, $r);
	if(!exists($printedvars{$outrow[0]}{$outrow[1]}{$outrow[2]}{$outrow[3]})){
	    $varcount++;
	    $printedvars{$outrow[0]}{$outrow[1]}{$outrow[2]}{$outrow[3]} = 1;
	    
	}
	unshift(@outrow, $varcount);
	if($outrow[20] eq "5" && $outrow[21] eq "CMGS_VGKL_5"){
	    $format = $five;
	}
	elsif($outrow[20] eq "4" && $outrow[21] eq "CMGS_VGKL_5"){
	    $format = $four;
	}
	elsif($outrow[20] eq "3" && $outrow[21] eq "CMGS_VGKL_5"){
	    $format = $three;
	}	
	elsif($outrow[20] eq "2" && $outrow[21] eq "CMGS_VGKL_5"){
	    $format = $two;
	}
	elsif($outrow[20] eq "1" && $outrow[21] eq "CMGS_VGKL_5"){
	    $format = $one;
	}
	else{
	    $format = $zero;
	}
	$mlink = "http://localhost:10000/show?request=" . $outrow[1] . ":" . $outrow[2] . "-" . $outrow[2];
	$blink = "http://localhost:10000/show?request=BAM<$c{'winshare_root'}" . $sample . "\\" . $rundate . ".$nowtime" . "\\" . $bam_file;
	$ws->write_row($i,0, \@outrow, $format);
	$ws->write_url($i, $#outrow+1, $mlink, $format, 'mut_link');
	$ws->write_url($i, $#outrow+2, $blink, $format, 'bam_link');
	$i++;
    }
    
    $i+=1;

    foreach $r (@ONEROWS) {
	@outrow = split(/;/, $r);
	if(!exists($printedvars{$outrow[0]}{$outrow[1]}{$outrow[2]}{$outrow[3]})){
	    $varcount++;
	    $printedvars{$outrow[0]}{$outrow[1]}{$outrow[2]}{$outrow[3]} = 1;
	}
	unshift(@outrow, $varcount);
	if($outrow[20] eq "5" && $outrow[21] eq "CMGS_VGKL_5"){
	    $format = $five;
	}
	elsif($outrow[20] eq "4" && $outrow[21] eq "CMGS_VGKL_5"){
	    $format = $four;
	}
	elsif($outrow[20] eq "3" && $outrow[21] eq "CMGS_VGKL_5"){
	    $format = $three;
	}	
	elsif($outrow[20] eq "2" && $outrow[21] eq "CMGS_VGKL_5"){
	    $format = $two;
	}
	elsif($outrow[20] eq "1" && $outrow[21] eq "CMGS_VGKL_5"){
	    $format = $one;
	}
	else{
	    $format = $zero;
	}
	$mlink = "http://localhost:10000/show?request=" . $outrow[1] . ":" . $outrow[2] . "-" . $outrow[2];
	$blink = "http://localhost:10000/show?request=BAM<$c{'winshare_root'}" . $sample . "\\" . $rundate . ".$nowtime" . "\\" . $bam_file;
	$ws->write_row($i,0, \@outrow, $format);
	$ws->write_url($i, $#outrow+1, $mlink, $format, 'mut_link');
	$ws->write_url($i, $#outrow+2, $blink, $format, 'bam_link');
	$i++;
    }
    $afilt .= ":AQ" . $i;
    $ws->autofilter($afilt);

    $i+=3;
    
    $format = $wb->add_format(bg_color=>$clr{'red'},border   => 1);
    $ws->write_blank($i, 1, $format);
    $format = $wb->add_format(bg_color=> 9,border   => 1);
    $ws->write($i, 2, 'class5 sannolikt_patogen');
    $ws->write($i, 3, 'Ska rapporteras i Alamut');

    $i++;

    $format= $wb->add_format(bg_color=>$clr{'pink'},border   => 1);
    $ws->write_blank($i, 1, $format);
    $format = $wb->add_format(bg_color=> 9,border   => 1);
    $ws->write($i, 2, 'class4 sannolikt_patogen');
    $ws->write($i, 3, 'Ska rapporteras i Alamut');

    $i++;
    
    $format= $wb->add_format(bg_color=>$clr{'yellow'},border   => 1);
    $ws->write_blank($i, 1, $format);
    $format = $wb->add_format(bg_color=> 9,border   => 1);
    $ws->write($i, 2, 'class3 VUS');
    $ws->write($i, 3, 'Ska rapporteras i Alamut');
 
    $i++;
 
    $format= $wb->add_format(bg_color=>$clr{'lgreen'},border   => 1);
    $ws->write_blank($i, 1, $format);
    $format = $wb->add_format(bg_color=> 9,border   => 1);
    $ws->write($i, 2, 'class2 sannolikt_benign');
    $ws->write($i, 3, 'Ska rapporteras i Alamut missense');
 
    $i++;

    $format= $wb->add_format(bg_color=>$clr{'green'},border   => 1);
    $ws->write_blank($i, 1, $format);
    $format = $wb->add_format(bg_color=> 9,border   => 1);
    $ws->write($i, 2, 'class1 benign');
    $ws->write($i, 3, 'Rapporteras ej i Alamut: allelfrekvens 5-95%, ingen spliceverkan');
    $format = $wb->add_format(bg_color=> 9,border   => 1);


    $outfile_le30merge = $sampledir . "/" . "sample.panel.bg.coverage.le30merge.bed";

    open(FH0, "<$outfile_le30merge");
    chomp(@CLE30 = <FH0>);
    close(FH0);
    
    $ws3 = $wb->add_worksheet('cov<30');
    $ws3->write(0,0, 'chrom');
    $ws3->write(0,1, 'start');
    $ws3->write(0,2, 'end');
    $ws3->write(0,3, 'size');	
    $ws3->write(0,4, 'gene');
    $ws3->write(0,5, 'reg_link');
    $ws3->write(0,6, 'bam_link');
    
    $i = 0;
    foreach $row (@CLE30) {
	@tmp = split(/\t/, $row);
	splice(@tmp, 4, 3);
	$tmp[3] =~ s/\s+//g;
	$size = $tmp[2] - $tmp[1];
	splice(@tmp, 3, 0, $size);
	if (exists($genesHash{$tmp[4]})) {
	    $ws3->write_row($i+1,0, \@tmp);
	    
	    $start = $tmp[1] +1;
	    $end = $tmp[2];
	    
	    $mlink = "http://localhost:10000/show?request=" . $tmp[0] . ":" . $start . "-" . $end;
	    $blink = "http://localhost:10000/show?request=BAM<$c{'winshare_root'}" . $winpath . "\\". $bam_file;	
	    
	    $ws3->write_url($i+1, 5, $mlink, $format, 'reg_link');
	    $ws3->write_url($i+1, 6, $blink, $format, 'bam_link');
	    $i++;
	}
    }
    
    $corefile = $sampledir . "/" . "sample.core.bg.coverage.le30merge.bed";
    open(FH0, "<$corefile");
    chomp(@CORE = <FH0>);
    close(FH0);
    
    $ws3b = $wb->add_worksheet('core-regioner');
    $ws3b->write(0,0, 'chrom');
    $ws3b->write(0,1, 'start');
    $ws3b->write(0,2, 'end');
    $ws3b->write(0,3, 'size');	
    $ws3b->write(0,4, 'gene');
    $ws3b->write(0,5, 'reg_link');
    $ws3b->write(0,6, 'bam_link');
    
    $i = 0;
    foreach $row (@CORE) {
	@tmp = ();
	if($row =~ /core-regioner/){
	    print "$row\n";
	    $tmp[0] = $row;
	    $ws3b->write_row($i+1,0, \@tmp);
	    $i++;
	}
	else{
	    @tmp = split(/\t/, $row);
	    splice(@tmp, 4, 3);
	    $tmp[3] =~ s/\s+//g;
	    $size = $tmp[2] - $tmp[1];
	    splice(@tmp, 3, 0, $size);
	    pop(@tmp);
	    pop(@tmp);
	    $ws3b->write_row($i+1,0, \@tmp);


	    $start = $tmp[1] +1;
	    $end = $tmp[2];
	    
	    $mlink = "http://localhost:10000/show?request=" . $tmp[0] . ":" . $start . "-" . $end;
	    $blink = "http://localhost:10000/show?request=BAM<$c{'winshare_root'}" . $winpath . "\\". $bam_file;	
	    
	    $ws3b->write_url($i+1, 5, $mlink, $format, 'reg_link');
	    $ws3b->write_url($i+1, 6, $blink, $format, 'bam_link');
	    $i++;
	}
    }

    
    $regcovfile =  $sampledir . "/" . "sample.perbase.coverage.bed";

    open(FH0, "<$regcovfile");
    chomp(@RCF = <FH0>);
    close(FH0);

    %genecov = ();
    %allcov = ();
    %genecov_30 = ();	

    foreach $row (@RCF){
#	print "$row\n";
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
    $ws2 = $wb->add_worksheet('selected genes');
    
    $ws2->write(0,0, 'Genes in selection file');
    $ws2->write(0,1, 'Disease group');
    $ws2->write(0,2, 'Exists on panel');
    
    print "$sampleC{$sample}{disease}\n";

    $i = 0;
    foreach $gene (sort { $a cmp $b } keys(%genesHash)){
	print "$gene\n";
	$gene =~ s/\s+//g;
	$ws2->write($i+1,0, $gene);
	$ws2->write($i+1,1, $sampleC{$sample}{disease});
	
	if(exists($genecov{$gene})) {
	    $ws2->write($i+1,2, "yes");
	}
	else{
	    $ws2->write($i+1,2, "no");
	}  
	$i++;
    }
    
    $ws4 = $wb->add_worksheet('genes coverage & completeness');
    $ws4->write(0,0, 'chrom');
    $ws4->write(0,1, 'start');
    $ws4->write(0,2, 'end');
    $ws4->write(0,3, 'gene');
    $ws4->write(0,4, 'av_cov');
    $ws4->write(0,5, 'completeness_30x');	
    
    $gci = 0;
    
    foreach $gene (sort { $a cmp $b } keys (%genecov)) {
	$gene =~ s/\s+//g;
	if (exists($genesHash{$gene})) {
	    $chr = $genecov{$gene}{chromosome};
	    $max = max @{$genecov{$gene}{positions}};
	    $min = min @{$genecov{$gene}{positions}};
	    $cov = sprintf("%.2f", $genecov{$gene}{addcov}/$genecov{$gene}{rows});
	    
	    if(!exists($genecov_30{$gene}{rows})){
		$genecov_30{$gene}{rows} = 0;
	    }
	    $completeness = sprintf("%.2f", $genecov_30{$gene}{rows}/$genecov{$gene}{rows}*100) . "%";
	    
	    $ws4->write($gci+1,0, $chr);
	    $ws4->write($gci+1,1, $max);
	    $ws4->write($gci+1,2, $min);
	    $ws4->write($gci+1,3, $gene);
	    $ws4->write($gci+1,4, $cov);
	    $ws4->write($gci+1,5, $completeness);		
	    
	    $gci++;
	}
    }
    
    $destdir = $sampledir . "/qc/teqc";
    
    if ( !-d $destdir ) {
	system("mkdir $destdir");
    }
    
    $ws5 = $wb->add_worksheet('quality');
    $ws5->write(1,0, 'TEQC:');
    $teqclink = $c{'winshare_root'} . $winpath . "\\" . "qc\\teqc\\index.html";
    $ws5->write_url(1, 1, $teqclink, $format, 'teqc_link');
    
    $ws5->write(2,0, 'FASTQC:');
    $fqclink = $c{'winshare_root'} . $winpath . "\\" . "qc\\fastqc\\fastqc_report.html";
    $ws5->write_url(2, 1, $fqclink, $format, 'fastqc_link');
    
    $ws5->write(3,0, 'Commands log:');
    $fqclink = $c{'winshare_root'} . $winpath . "\\" . "qc\\bcbio-nextgen-commands.log";
    $ws5->write_url(3, 1, $fqclink, $format, 'commands_link');
    $diseasegenesdef_link = $diseasegenesdef_file . "_link";
    $ws5->write(4,0, 'Indications/genes:');
    $fqclink = $c{'winshare_root'} . $winpath . "\\" . "qc\\" . $diseasegenesdef_link;
    $ws5->write_url(4, 1, $fqclink, $format, $diseasegenesdef_file);
    
    
    $bamtoolsfile = $sampledir . "/qc/bamtools/bamtools_stats.txt";
    
    open (FH1,"<$bamtoolsfile");
    chomp(@BT = <FH1>);
    close(FH1);
    
    for ($bti = 0; $bti <= $#BT; $bti++) {
	$ws5->write($bti+5,0,$BT[$bti]);
    }
    $ws6 = $wb->add_worksheet('bcbio config parameters');
    $crow = 0;
    foreach $row (sort { $a <=> $b } keys(%{$configHash{$sample}})) {
	$ws6->write($crow,0,$configHash{$sample}{$row});
	$crow++;
    }
    $ws7 = $wb->add_worksheet('programs');
    $crow = 0;
    foreach $row (@PROGS) {
	$ws7->write($crow,0,$row);
	$crow++;
    }
    $ws8 = $wb->add_worksheet('analyzed regions');
    $crow = 0;
    foreach $row (@REG) {
	@tmp = split(/\t/, $row);
	if(exists($genesHash{$tmp[3]})){
#	    print "write_row @tmp\n";
	    $ws8->write_row($crow,0,\@tmp);
	    $crow++;
	}
    }
    $configdatfile = $sampledir . "/qc/sample_config_file.csv";
    open(CFG, ">$configdatfile");
    print CFG "#sample,gender,regions,disease\n";
    print CFG "$sample,$sampleC{$sample}{sex},$sampleC{$sample}{regions},$sampleC{$sample}{disease}\n";
    close(CFG);

    $qc_from = $sampledir . "/qc";
    $qc_to = $sampledir . "/" . $sample . "/" . $rundate . "." . $nowtime;
    $bamindex = $bam . ".bai";

    system("cp $diseasegenesdef $qc_from");
    system("cp $summary_file $qc_from"); 
    system("cp $programs $qc_from");
    
    system("cp -r $qc_from $qc_to");
    system("cp $bam $qc_to");
    system("cp $bamindex $qc_to");

}

