#!/usr/bin/perl -w
#

#vep: variant_effect_predictor.pl --cache --vcf -i batch1-ensemble.vcf -o batch1-ensemble-vep.vcf --refseq --offline --force_overwrite --everything --hgvs --fork 16 --fasta /home/data_in/.vep/homo_sapiens/80_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa


use IPC::Run;
use List::Util qw(max min);
use Spreadsheet::WriteExcel;
use Time::localtime;
use Getopt::Long;
use Sort::Key::Natural;
use Sort::Key::Maker int_nat_sort => qw(integer natural);

$o{'in_silico'} = 0;
$o{'help'} = 0;
$o{'version'} = "v0.45";

## Read configs
open (FH, "</home/data_in/calls2xls/calls2xls.cfg");
while ( $line = <FH> ) {
    if($line !~ /#/ && $line =~ /\S/){
	@tmp = split(/\s=\s/, $line);
	foreach $row (@tmp) {
	    $row =~ s/\s//g;
	}
	if($tmp[0] eq "info_dp"){
	    @info = split(/,/,$tmp[1]);
	    foreach $i (@info){
		$c{$tmp[0]}{$i} = 1;
	    }
	}
	elsif($tmp[0] eq "info_ao"){
	    @info = split(/,/,$tmp[1]);
	    foreach $i (@info){
		$c{$tmp[0]}{$i} = 1;
	    }
	}
	else{
	    $c{$tmp[0]}=$tmp[1];
#	    print "$tmp[0]=$tmp[1]\n";
	}
    }
}

close(FH);

########## get options and define run mode ##########

GetOptions (\%o, "help", "runfolder=s", "old_sample=s", "new_sample=s", "disease_hypothesis=s", "disease-gene-file=s", "ped-file=s");

if(!(@ARGV)){
    &helpmess0();
    exit 0;  
}

chomp($subcommand = $ARGV[0]);

if(($subcommand eq "single" || $subcommand eq "insilico") || $subcommand eq "family"){
    if($subcommand eq "single" && !exists($o{'runfolder'})){
	&helpmess1();
	exit 0;
    }
    if($subcommand eq "insilico" && ((((!exists($o{'runfolder'}) || !exists($o{'old_sample'})) || !exists($o{'new_sample'})) || !exists($o{'disese_hypothesis'})) || !exists($o{'disese-gene-file'}))){
	&helpmess2();
	exit 0;
    }
    if($subcommand eq "family" && ((!exists($o{'runfolder'}) || !exists($o{'disese_hypothesis'})) || !exists($o{'disese-gene-file'}) || !exists($o{'ped-file'}))){
	&helpmess3();
	exit 0;
    }
}
elsif($subcommand eq "pdga"){
    &pdga();
}
elsif($subcommand eq "pped"){
    &ppedfiles();
}
elsif($subcommand eq "pcfg"){
    &pcfgfiles();
}
else{
    &helpmess0();
    exit 0;
}
if($o{'help'}){
    &helpmess0();
    exit 0;
}
########## help messge subroutines  ##########

sub helpmess0{
    print "\nversion: $o{'version'}\n";
    print "usage: calls2xls.pl subcommand [options]\n\n";
    print "subcommands:     \n";    
    print "  single     Annotate single single samples\n";
    print "  insilico   Transfer data to new sample number and annotate\n";
    print "  family     Transfer data to new sample number and annotate\n"; 
    print "             using family information\n\n";
    print "  pdga       Print disease-genes-association files\n";
    print "  pped       Print available ped files\n";
    print "  pcfg       Print available config files\n\n";
}
sub helpmess1{
    print "\nversion: $o{'version'}\n";
    print "\nusage: calls2xls.pl single [options]\n\n";
    print "  --runfolder=path (required)\n";
    print "      Path to runfolder containing raw data and variant calls \n";
    print "      from bcbio-nextgen.\n";    
    print "  --config_file=file (optional)\n";
    print "      Path to config file with rows for sample IDs, disease hypothesis, \n"; 
    print "      gender, disease-gene file. Use if configurations are not defined\n";
    print "      in SampleSheet.csv.\n";
    print "      Example row: 15-1234 male HCM disease-gene-file.txt\n\n"; 
} 


sub helpmess2{
    print "\nversion: $o{'version'}\n";
    print "usage: calls2xls.pl insilico [options]\n\n";
    print "  --runfolder=path (required)\n";
    print "    path to data root folder containing raw data and variant calls \n";
    print "    from bcbio-nextgen.\n";    
    print "  --old_sample=sample number, e.g. 15-1234\n";
    print "      old sample number to transfer calling data from.\n";
    print "  --new_sample=sample_number, e.g. 15-1235\n";
    print "      new sample number to transfer calling data to.\n";
    print "  --disease_hypothesis=hypothesis, e.g. HCM\n";
    print "      new disease indication.\n";    
    print "  --disease-gene-file=file\n";
    print "      path to disease-gene associations.\n\n";
} 
sub helpmess3{
    print "\nversion: $o{'version'}\n";
    print "usage: calls2xls.pl family [options]\n\n";
    print "  --runfolder=path (required)\n";
    print "    path to funfolder containing raw data and variant calls \n";
    print "    from bcbio-nextgen.\n";    
    print "  --ped-file=file\n";
    print "      ped file for inheritance annotations.\n";
    print "  --disease-gene-file=file\n";
    print "      disease-gene associations file\n\n";
    print "  --disease-gene-file=file\n";
} 

sub pdga{
    print "\nversion: $o{'version'}\n";
    print "usage: calls2xls.pl subcommand [options]\n";
    chomp(@dga = <$c{'diseaseGeneAssocPath'}/*.txt>);
    print "\n";
    if(exists($dga[0])) {
	print "disease gene association files and disease hypotheses:\n";
    }
    foreach $file (@dga){
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
sub ppedfiles{
    print "\nversion: $o{'version'}\n";
    print "usage: calls2xls.pl subcommand [options]\n";
    chomp(@pedf = <$c{'pedFilesPath'}/*.ped>);
    print "\n";
    if(exists($pedf[0])) {
	print "ped files:\n";
    }
    else{
	print "no ped files available\n";
    }
    if(exists($pedf[0])) {
	foreach $file (@pedf){
	    @tmp = split(/\//, $file);
	    print "  $tmp[-1]\n";
	}
    }
}
sub pcfgfiles{
    print "\nversion: $o{'version'}\n";
    print "usage: calls2xls.pl subcommand [options]\n";
    chomp(@pcfgf= <$c{'configFilesPath'}/*.csv>);
    print "\n";
    if(exists($pcfgf[0])) {
	print "config files:\n";
    }
    foreach $file (@pcfgf){
	@tmp = split(/\//, $file);
	print "  $tmp[-1]\n";
    }
    
    print "\n";
}


chomp($rundir = $o{'runfolder'});

#print "$rundir\n";

#open(PRG, "java -jar ~/bin/SnpSift.jar 2>&1 |");
#chomp(@snpsift_version = <PRG>);
#close(PRG);

#open(PRG, "vt normalize 2>&1 |");
#chomp(@vt_version = <PRG>);
#close(PRG);


#colors
%clr=(
    'red' => 10,
    'pink' => 45,
    'yellow' => 13,
    'lgreen'=> 42,
    'green' => 11
    );

$finalpath = $rundir . "/" . $c{'call_output_root'};

opendir(DIR, $finalpath);
@files = grep { /^project_/ } readdir(DIR);
closedir(DIR);
$projectpath = $finalpath . "/" . $files[0];

$summaryfile = $projectpath . '/' . $c{'callRunSummaryFile'};
$samplesheetfile = $rundir . "/SampleSheet.csv";

$rundate = (split /_/, $rundir)[0];
open(RD, "read_bcbio_summary_yaml.pl --summary_file $summaryfile --output_type date  2>&1 |");
chomp($calldate = <RD>);
close(RD);

#print "$calldate\n";

$tm = localtime;
$nowtime = sprintf("%04d.%02d.%02d-%02d.%02d.%02d", $tm->year+1900, ($tm->mon)+1, $tm->mday, $tm->hour, $tm->min, $tm->sec );

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
	    $ssheet{$field}{$tmp[0]} = $tmp[1];
	}
	elsif($#tmp == 0){
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
#			    print "$tmp2[$j]\n";
			    $ssheet{$field}{$header[$i]}{$tmp[0]}{$desc_header[$j]} = $tmp2[$j];
			}
		    }
		}
	    }
	}
    }
}

#move to sampleConfig
foreach $sample ( sort { $b cmp $a } keys(%{$ssheet{'Data'}{'Description'}})) {
    $sampleConfig{$sample}{gender} = $ssheet{'Data'}{'Description'}{$sample}{'gender'};
    $sampleConfig{$sample}{regions} = $ssheet{'Data'}{'Description'}{$sample}{'design'};
    $sampleConfig{$sample}{disease} = $ssheet{'Data'}{'Description'}{$sample}{'inquiry'};
    $sampleConfig{$sample}{genelist} = $ssheet{'Data'}{'Description'}{$sample}{'genelist'};
    $sampleConfig{$sample}{disease} =~ s/\s//g;
}

@variant_fields = split(/,/, $c{'variant_field_order'});

@vfhtmp = split(/,/, $c{'variant_header'});
foreach $row (@vfhtmp){
#    print "vcftmp $row\n";
    @kv = split(/:/, $row);
    $vfh{$kv[0]} = $kv[1];
#    print "$vfh{$kv[0]} = $kv[1]\n";
}

@csztmp = split(/,/, $c{'var_cols_widths'});
@kv = split(/:/, $csztmp[-1]);
$variant_max_col = $kv[0];

foreach $row (@csztmp){
    @kv = split(/:/, $row);
    $k = $kv[0] . ":" . $kv[0];
    $csz{$k}=$kv[1];
}

if($o{'in_silico'}){
    foreach $field ( sort { $b cmp $a } keys(%{$sampleConfig{$o{'old_sample'}}})) {
	$sampleConfig{$o{'new_sample'}}{$field} = $sampleConfig{$o{'old_sample'}}{$field};
#	print "$o{'new_sample'} $sampleConfig{$o{'new_sample'}}{$field}\n";
    }
    
    $sampleConfig{$o{'new_sample'}}{genelist} = $o{'disease_gene_file'};
    $sampleConfig{$o{'new_sample'}}{disease} = $o{'disease_indication'};
#    print "$sampleConfig{$o{'new_sample'}}{genelist}\n";
#    print "$sampleConfig{$o{'new_sample'}}{disease}\n";
    
    foreach $sample ( sort { $b cmp $a } keys(%sampleConfig)) {
	if($sample ne $o{'new_sample'}){
	    delete $sampleConfig{$sample};
	}
    }
    
    $newsamppath = $finalpath . "/" . $o{'new_sample'};

#    print "$newsamppath\n";

    if(! -d $newsamppath){
	system "mkdir $newsamppath\n";
	system "cp -r $finalpath/$o{'old_sample'}/* $finalpath/$o{'new_sample'}";
    	system "rm -r $finalpath/$o{'new_sample'}/$o{'old_sample'}";
	system "rename s/$o{'old_sample'}-/$o{'new_sample'}-/g $newsamppath/*";
    }
}

foreach $sample ( sort { $b cmp $a } keys(%sampleConfig)) {
    $dir = $finalpath . "/" . $sample;
    if($o{'in_silico'}){
	if(($dir !~ /project/ && -d $dir) && $dir =~ /$o{'new_sample'}/) {
	    push @sampledirs, $dir;
	}
    }
    else{
	if($dir !~ /project/ && -d $dir) {
	    push @sampledirs, $dir;
	}
    }
}

system("mut2vcf.pl $projectpath");
print "mut2vcf.pl $projectpath\n";


#get programs
$programs = $projectpath . "/programs.txt";
open(FH, "<$programs") or die $!;
chomp(@PROGS = <FH>);
close(FH) or die $!;

#get calling.yaml
$bcbio_config = $rundir . "/" . $c{'bcbio_config'};
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
foreach $sample ( sort { $b cmp $a } keys(%sampleConfig)) {
#    print "$sample\n";
    $diseasegenesdef = $c{'diseaseGeneAssocPath'} . "/" . $sampleConfig{$sample}{genelist};
    $diseasegenesdef_file = $sampleConfig{$sample}{genelist};

#    print "$diseasegenesdef $diseasegenesdef_file\n";

    open(FH0, "<$diseasegenesdef");
    chomp(@GeneL = <FH0>);
    close(FH0);

    foreach $row (@GeneL) {
	#$row =~ s/\s+//g;
	if($row =~ /^>(\S+)/) {
	    $dise = $1;
	}
	elsif($row =~ /\S/){
	    @tmp = split(/\s+/, $row);
#	    print "$dise $sampleConfig{$sample}{genelist} @tmp\n";
	    $disease2genes{$sample}{$dise}{$tmp[0]}{'core'} = $tmp[1];
	    $disease2genes{$sample}{$dise}{$tmp[0]}{'vartype'} = $tmp[2];

#	    print "$sample $dise $tmp[0] $disease2genes{$sample}{$dise}{$tmp[0]}{'vartype'}  $disease2genes{$sample}{$dise}{$tmp[0]}{'core'}\n";

	}
    }
}

foreach $sample (sort { $a cmp $b } keys(%sampleConfig)) {
#    print "$sample\n";
    $ccnf = $c{'designAndGeneListsRootPath'} . "/" .  $sampleConfig{$sample}{regions};
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
	
	if(exists($disease2genes{$sample}{$sampleConfig{$sample}{disease}}{$tmp[3]})){
	    for($i = $tmp[1] ; $i<= $tmp[2]; $i++) {
		$onTarget{$sample}{$chr}{$i} = $tmp[3];
	    }
	}
    }
}

$command0 = "comm." . $rundir . ".0.txt";
$command1 = "comm." . $rundir . ".1.txt";
$command2 = "comm." . $rundir . ".2.txt";
$command3 = "comm." . $rundir . ".3.txt";
$command4 = "comm." . $rundir . ".4.txt";

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
#    print "parallel $sampledir\n";
    $sample = (split /\//, $sampledir )[-1];
    @tmp = split(/\//, $sampledir);
    $bam = $sampledir . "/" . $tmp[-1] . "-ready.bam";

    $reg =  $c{'designAndGeneListsRootPath'} . "/" .  $sampleConfig{$sample}{regions};

    $qcdir = $sampledir . "/qc";

    if ( !-d $qcdir ) {
	system("mkdir $qcdir");
    }

    $destdir = $sampledir . "/qc/teqc";
    if ( !-d $destdir ) {
	system("mkdir -p $destdir");
    }

    $sc = $sampledir . "/" . "sample.coverage.bed";
    $perbasec = $sampledir . "/" . "sample.perbase.coverage.bed";
    $le30merge = $sampledir . "/" . "sample.panel.bg.coverage.le30merge.bed";

    $bam = $sampledir . "/" . $sample . "-ready.bam";
    
#    print "$sc $perbasec $le30merge $bam\n";

    print FH0 "bedtools coverage -d -b $bam -a $reg > $perbasec\n";
    print FH1 "bedtools genomecov -ibam $bam -bga  | bedtools intersect -wb -a stdin -b $reg | cut -f 1-4,8 | sed \'s/.*/&\t$sample\t$rundate\t$calldate\t$sampleConfig{$sample}{gender}/\' >$sc\n";
    print FH2 "filterCoverage.pl $sc | sort -k1,1 -k2,2n | bedtools merge -i stdin -c 5,6,7,8 -o distinct,distinct,distinct,distinct > $le30merge\n";
    print FH3 "samtools index $bam\n";
    print FH4 "runTEQC.pl $bam $reg $sample $destdir\n";
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


if(!$o{'in_silico'}){ 
    system("putLocAfDb.pl $rundir");
}

#unpack ensemble vcf
$ensvcfgz = $projectpath . "/batch1-ensemble.vcf.gz";
$ensvcf = $projectpath . "/batch1-ensemble.vcf";
$ensdecompvcf = $projectpath . "/batch1-ensemble.d.vcf";
$ann1vcf = $projectpath . "/batch1-ensemble.d.1kg.esp.exac.vcf";
$ann2vcf = $projectpath . "/batch1-ensemble.d.1kg.esp.exac.ala.vcf";
$ann3vcf = $projectpath . "/batch1-ensemble.d.1kg.esp.exac.ala.frq.vcf";
$ann4vcf = $projectpath . "/batch1-ensemble.d.1kg.esp.exac.ala.frq.ann.vcf";
if($o{'in_silico'}) {
    $ann4avcf = $projectpath . "/batch1-ensemble.d.1kg.esp.exac.ala.frq.ann." . $o{'new_sample'} . ".vcf";
}

system("pigz -f -d -k $ensvcfgz");
system("vt decompose $ensvcf -o $ensdecompvcf");
system("snpsift_annotate_esp_1kg_exac.pl $ensdecompvcf");
system("snpsift_annotate_alamut.pl $ann1vcf");
system("getLocAf.pl $ann2vcf > $ann3vcf");

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
		    
		    if($tmp[$i] eq $o{'old_sample'}){
			$s2c{$o{'new_sample'}} = $i;
		    }
		}
		$line .= "\t" . $o{'new_sample'};
		print FH2 "$line\n";
	    }
	    elsif($line =~ /^#/){
		print FH2 "$line\n";
	    }
	    else{
		@tmp = split(/\t/, $line);
		push @tmp, $tmp[$s2c{$o{'new_sample'}}];
		$line = join("\t", @tmp);
		print FH2 "$line\n";
	    }
	}
    }
    close(FH);
    close(FH2);
    $ann4vcf = $ann4avcf;
}
#print "$ann4vcf\n";
open(FH, "<$ann4vcf") or die $!;
chomp(@gout = <FH>);
close(FH) or die $!;

%sample2column = ();
%column2sample = ();
%calls = ();
#print "$c{'alamut_index'}\n";

$i = 0;
foreach $row (@gout) {
#    print "$row\n";
    if($row =~ /^#CHROM/) {

	@tmp = split(/\s+/, $row);

	for ($j = 8; $j <= $#tmp; $j++) {
	    $tmp[$j] =~ s/\s+//g;
	    $column2sample{$j} = $tmp[$j];
	    $sample2column{$tmp[$j]} = $j;
	}
    }

    elsif($row !~ /^#/ && $row =~ /\S/) {
#	print "$row\n";
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
	$alcindex = "0";
	$alcval = ".";
	$alcreate = ".";
	$alupdate = ".";
	
	%annhash = ();

	@INFO = split(/;/, $tmp[7]);

	%infohash = ();
	foreach $in (@INFO) {
	    if($in =~ /=/){
		($index, $val) = split(/=/, $in);
		$infohash{$index} = $val;
#		print "$index  $infohash{$index}\n";
	    }
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
	if(exists($infohash{'ALCINDEX'}) && exists($infohash{'ALCVAL'})){
#	    print "c" . $infohash{'ALCVAL'} . "x" . $c{'alamut_index'} . "c\n";
	    if($infohash{'ALCVAL'} eq $c{'alamut_index'}){
		$alcindex = $infohash{'ALCINDEX'};
		$alcval = $infohash{'ALCVAL'};
#		print "$alcindex $alcval\n";
	    }
	}
	    
	if(exists($infohash{'ALAMUT'})) {
	    $alamut = $infohash{'ALAMUT'}; 
	}
	if(exists($infohash{'ALCREATE'})) {
	    $alcreate = $infohash{'ALCREATE'}; 
	}
	if(exists($infohash{'ALUPDATE'})) {
	    $alupdate = $infohash{'ALUPDATE'}; 
	}
	if(exists($infohash{'LOC_AF'})) {
	    $locfreq = sprintf("%.1f", $infohash{'LOC_AF'}* 100);  
	}
	if(exists($infohash{'LOC_AC'})) {
	    $locals = $infohash{'LOC_AC'};
#	    print "LOC found\n";
	}
	else{
#	    print "LOC not found\n";
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
	
	for ($j = 0; $j <= $#format; $j++){
	    if(exists($c{'info_dp'}{$format[$j]})) {
		$dp_no = $j;
	    }
	    elsif(exists($c{'info_ao'}{$format[$j]})) {
		$ao_no = $j;
	    }
	}

	foreach $sample (sort { $a cmp $b } keys(%sample2column)) {
#	    print "$sample\n";
	    @tmp2 = split(/:/, $tmp[$sample2column{$sample}]);
	    if($tmp2[0] =~ /[1-9]/ && exists($onTarget{$sample}{$chr}{$p})) {
		foreach $ano (sort { $a <=> $b } keys(%annhash)) {
		    
		    $variant = $chr . "_" . $p . "_" . $ref . "_" . $alt . "_" . $ano;

		    $calls{$sample}{$variant}{'gt'} =  $tmp2[0];

		    if($dp_no ne "NA") {
			$calls{$sample}{$variant}{'dp'} =  $tmp2[$dp_no];
		    }
		    else{
			$calls{$sample}{$variant}{'dp'} =  "NA";
		    }
		    if($ao_no ne "NA") {
			$calls{$sample}{$variant}{'ao'} =  $tmp2[$ao_no];
		    }
		    else{
			$calls{$sample}{$variant}{'ao'} =  "NA";
		    }
		    $calls{$sample}{$variant}{'chr'} = $chr;
		    $calls{$sample}{$variant}{'p'} = $p;
		    $calls{$sample}{$variant}{'ref'} = $ref;
		    $calls{$sample}{$variant}{'alt'} = $alt;
		    $calls{$sample}{$variant}{'exac_af'} = $exac_af;
		    $calls{$sample}{$variant}{'esp_af_ea'} = $esp_af_ea;
		    $calls{$sample}{$variant}{'esp_af_aa'} = $esp_af_aa;
		    $calls{$sample}{$variant}{'1kg_af_all'} = $onekg_af_all;
		    $calls{$sample}{$variant}{'1kg_af_afr'} = $onekg_af_afr;
		    $calls{$sample}{$variant}{'1kg_af_amr'} = $onekg_af_amr;
		    $calls{$sample}{$variant}{'1kg_af_eas'} = $onekg_af_eas;
		    $calls{$sample}{$variant}{'1kg_af_eur'} = $onekg_af_eur;
		    $calls{$sample}{$variant}{'1kg_af_sas'} = $onekg_af_sas;
		    $calls{$sample}{$variant}{'type'} = $annhash{$ano}{type};
		    $calls{$sample}{$variant}{'impact'} = $annhash{$ano}{impact};
		    $calls{$sample}{$variant}{'genename'} = $annhash{$ano}{genename};
		    $calls{$sample}{$variant}{'fid'} = $annhash{$ano}{fid};
		    $calls{$sample}{$variant}{'btype'} = $annhash{$ano}{btype};
		    $calls{$sample}{$variant}{'hgvsc'} = $annhash{$ano}{hgvsc};
		    $calls{$sample}{$variant}{'hgvsp'} = $annhash{$ano}{hgvsp};
		    $calls{$sample}{$variant}{'alamut'} = $alamut;
		    $calls{$sample}{$variant}{'alcindex'} = $alcindex;
		    $calls{$sample}{$variant}{'alcval'} = $alcval;
		    $calls{$sample}{$variant}{'alcreated'} = $alcreate;
		    $calls{$sample}{$variant}{'alupdated'} = $alupdate;
		    $calls{$sample}{$variant}{'loc_AF'} = $locfreq;
		    $calls{$sample}{$variant}{'loc_AC'} = $locals;
		    $calls{$sample}{$variant}{'targetGene'} = $onTarget{$sample}{$calls{$sample}{$variant}{'chr'}}{$calls{$sample}{$variant}{'p'}};
#		    print "$sample $variant $calls{$sample}{$variant}{'targetGene'}\n";
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

    foreach $gene (sort { $a cmp $b } keys(%{$disease2genes{$sample}{$sampleConfig{$sample}{'disease'}}})) {
#	print "1. $gene $sampleConfig{$sample}{'disease'} \n";
	$gene =~ s/\s+//g;
	$genesHash{$gene}{'core'} = $disease2genes{$sample}{$sampleConfig{$sample}{'disease'}}{$gene}{'core'};
	$genesHash{$gene}{'vartype'} = $disease2genes{$sample}{$sampleConfig{$sample}{'disease'}}{$gene}{'vartype'};

#	print "2. $sample $sampleConfig{$sample}{disease} $gene $genesHash{$gene}{'core'} $genesHash{$gene}{'vartype'}\n";
	
    }

    $sampleout_rundate_folder = $sampledir . "/" . $sample . "/" . $rundate . "." . $nowtime;
    
    $winpath = $sample . "\\" . $rundate . "." . $nowtime;
    
    if(!-d $sampleout_rundate_folder) {
	system("mkdir -p $sampleout_rundate_folder");
    }
    
    $xls_name = $sampleout_rundate_folder . "/" . $sample . "_calldata.xls";
    
    my $wb = Spreadsheet::WriteExcel->new($xls_name);
    $ws = $wb->add_worksheet('calldata');
    
    $format = $wb->add_format();
    foreach $key (keys %csz){
	$ws->set_column($key, $csz{$key});
    }

    $ws->freeze_panes(1, 0);
    
    $format{'5'} = $wb->add_format(bg_color=>$clr{'red'},border   => 1);
    $format{'4'} = $wb->add_format(bg_color=>$clr{'pink'},border   => 1);
    $format{'3'} = $wb->add_format(bg_color=>$clr{'yellow'},border   => 1);
    $format{'2'} = $wb->add_format(bg_color=>$clr{'lgreen'},border   => 1);
    $format{'1'} = $wb->add_format(bg_color=>$clr{'green'},border   => 1);
    $format{'0'} = $wb->add_format(bg_color=> 9,border   => 1, );
    $format{'0b'} = $wb->add_format(bg_color=> 9,border   => 1, );
    $format{'0b'} -> set_bold();
    

    @indexorder = split(/,/, $c{'alamut_index_sort_order'});
    %sortmap = map { $indexorder[$_] => $_ } 0 .. $#indexorder;

    $wsrowno = 0;

    #populate header
    @houtrow = ('varno');
    foreach $v (@variant_fields){
#	print "$v\n";
	push @houtrow, $vfh{$v};
    }
    push @houtrow, "link";
    push @houtrow, "link";

    $ws->write_row($wsrowno,0, \@houtrow, $format{'0b'});
    $wsrowno++;

    $ws->write_row($wsrowno,0, \@houtrow, $format{'0'});
    $autofstart = $wsrowno + 1;
    $autofilt = "A" . $autofstart;
    $wsrowno++;

    %uniqvariants = ();

    $varno = 0;

    foreach $variant (int_nat_sort {   $sortmap{$calls{$sample}{$_}{'alcindex'}}, $_   } keys %{$calls{$sample}}){
	$genpos = $variant;
	$genpos =~ s/_\d+$//g;
	
#	print "$variant $calls{$sample}{$variant}{'alcindex'}\n";
	if (!exists($uniqvariants{$genpos})){
	    $uniqvariants{$genpos} = 1;
	    $varno++;
	}

	@outrow = ($varno);
	foreach $v (@variant_fields){
	    push @outrow, $calls{$sample}{$variant}{$v};
	}

	$mlink = $c{'mlink'} . $calls{$sample}{$variant}{'chr'} . ":" . $calls{$sample}{$variant}{'p'} . "-" . $calls{$sample}{$variant}{'p'};
	$blink = $c{'blink'} . $c{'winshare_root'} . $sample . "\\" . $rundate . ".$nowtime" . "\\" . $bam_file;

	$ws->write_row($wsrowno, 0, \@outrow, $format{$calls{$sample}{$variant}{'alcindex'}});
	$ws->write_url($wsrowno, $#outrow+1, $mlink, $format{$calls{$sample}{$variant}{'alcindex'}}, 'mut');
	$ws->write_url($wsrowno, $#outrow+2, $blink, $format{$calls{$sample}{$variant}{'alcindex'}}, 'bam');
	$wsrowno++;
#	print "$blink\n";
#	print "mlink\n";
#	print $c{'blink'} . "\n";
#	print $c{'mlink'} . "\n";
    }
        
    $autofilt .= ":" . $variant_max_col . $#houtrow;
    $ws->autofilter($autofilt);

    $wsrowno+=3;
        
    $ws->write_blank($wsrowno, 1, $format{"5"});
    $ws->write($wsrowno, 2, 'Class 5, Pathogenic', $format{"0"});

    $wsrowno++;

    $ws->write_blank($wsrowno, 1, $format{"4"});
    $ws->write($wsrowno, 2, 'Class 4, Likely pathogenic');

    $wsrowno++;
    
    $ws->write_blank($wsrowno, 1, $format{"3"});
    $ws->write($wsrowno, 2, 'Class 3, VUS');
 
    $wsrowno++;

    $ws->write_blank($wsrowno, 1, $format{"2"});
    $ws->write($wsrowno, 2, 'Class 2, Likely benign');
 
    $wsrowno++;

    $ws->write_blank($wsrowno, 1, $format{"1"});
    $ws->write($wsrowno, 2, 'Class 1, Benign');

    $le30merge = $sampledir . "/" . "sample.panel.bg.coverage.le30merge.bed";

    open(FH0, "<$le30merge");
    chomp(@CLE30 = <FH0>);
    close(FH0);

    $ws3 = $wb->add_worksheet('cov<30');
    $ws3->write(0,0, 'chrom', $format{'0b'});
    $ws3->write(0,1, 'start', $format{'0b'});
    $ws3->write(0,2, 'end', $format{'0b'});
    $ws3->write(0,3, 'size', $format{'0b'});	
    $ws3->write(0,4, 'gene', $format{'0b'});
    $ws3->write(0,5, 'core', $format{'0b'});
    $ws3->write(0,6, 'reg_link', $format{'0b'});
    $ws3->write(0,7, 'bam_link', $format{'0b'});
    
    $wsrowno = 0;
    foreach $row (@CLE30) {
	@tmp = split(/\t/, $row);
	splice(@tmp, 4, 3);
	$tmp[3] =~ s/\s+//g;
	$size = $tmp[2] - $tmp[1];
	splice(@tmp, 3, 0, $size);
#	print "hej @tmp\n";

	if (exists($genesHash{$tmp[4]})) {
#	    print "finns\n";
#	    foreach $urk (sort { $a cmp $b } keys(%{$genesHash{$tmp[4]}})) {
#		print "$urk $genesHash{$tmp[4]}{$urk}\n";
#	    }
#	    print "$genesHash{$tmp[4]}{'core'}\n";
	    push @tmp, $genesHash{$tmp[4]}{'core'};

	    if($genesHash{$tmp[4]}{'core'}  eq "+") {
		$ft_tmp = $format{'0b'};
	    }
	    else{
		$ft_tmp = $format{'0'};
	    }

	    $ws3->write_row($wsrowno+1,0, \@tmp, $ft_tmp);
	    
	    $start = $tmp[1];
	    $end = $tmp[2];
	    
	    $mlink = $c{'mlink'} . $tmp[0] . ":" . $start . "-" . $end;
	    $blink = $c{'blink'} . $c{'winshare_root'} . $winpath . "\\". $bam_file;	
	    
	    $ws3->write_url($wsrowno+1, 6, $mlink, $ft_tmp, 'reg_link');
	    $ws3->write_url($wsrowno+1, 7, $blink, $ft_tmp, 'bam_link');
	    $wsrowno++;
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
    
    $ws2->write(0,0, 'genes in disease hypothesis file', $format{'0b'});
    $ws2->write(0,1, 'disease hypothesis', $format{'0b'});
    $ws2->write(0,2, 'exists on panel', $format{'0b'});
    $ws2->write(0,3, 'core', $format{'0b'});
    $ws2->write(0,4, 'vartype', $format{'0b'});    
    
#    print "$sampleConfig{$sample}{disease}\n";

    $wsrowno = 1;
    foreach $gene (sort { $a cmp $b } keys(%genesHash)){
	print "$gene\n";
	
	if($genesHash{$gene}{'core'} eq "+") {
	    $ft_tmp = $format{'0b'};
	}
	else{
	    $ft_tmp = $format{'0'};
	}
	    
	$gene =~ s/\s+//g;
	$ws2->write($wsrowno,0, $gene, $ft_tmp);
	$ws2->write($wsrowno,1, $sampleConfig{$sample}{disease}, $ft_tmp);
	
	if(exists($genecov{$gene})) {
	    $ws2->write($wsrowno,2, "yes", $ft_tmp);
	}
	else{
	    $ws2->write($wsrowno,2, "no", $ft_tmp);
	}
	if($genesHash{$gene}{'core'} eq "+") {
	    $ws2->write($wsrowno,3, "+", $ft_tmp);
	}
	else{
	    $ws2->write($wsrowno,3, "-", $ft_tmp);
	}
	$ws2->write($wsrowno,4, $genesHash{$gene}{'vartype'}, $ft_tmp);
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
	$gene =~ s/\s+//g;
	if (exists($genesHash{$gene})) {
	    
	    if($genesHash{$gene}{'core'} eq "+") {
		$ft_tmp = $format{'0b'};
	    }
	    else{
		$ft_tmp = $format{'0'};
	    }

	    $chr = $genecov{$gene}{chromosome};
	    $max = max @{$genecov{$gene}{positions}};
	    $min = min @{$genecov{$gene}{positions}};
	    $cov = sprintf("%.2f", $genecov{$gene}{addcov}/$genecov{$gene}{rows});
	    
	    if(!exists($genecov_30{$gene}{rows})){
		$genecov_30{$gene}{rows} = 0;
	    }
	    $completeness = sprintf("%.2f", $genecov_30{$gene}{rows}/$genecov{$gene}{rows}*100) . "%";
	    
	    $ws4->write($gci,0, $chr,$ft_tmp);
	    $ws4->write($gci,1, $max,$ft_tmp);
	    $ws4->write($gci,2, $min,$ft_tmp);
	    $ws4->write($gci,3, $gene,$ft_tmp);
	    $ws4->write($gci,4, $genesHash{$gene}{'core'},$ft_tmp);
	    $ws4->write($gci,5, $cov,$ft_tmp);
	    $ws4->write($gci,6, $completeness,$ft_tmp);		
	    
	    $gci++;
	}
    }
    
    $destdir = $sampledir . "/qc/teqc";
    
    if ( !-d $destdir ) {
	system("mkdir $destdir");
    }
    
    $ws5 = $wb->add_worksheet('quality');
    $ws5->write(1,0, 'TEQC:', $format{'0b'});
    $teqclink = $c{'winshare_root'} . $winpath . "\\" . "qc\\teqc\\index.html";
    $ws5->write_url(1, 1, $teqclink, $format, 'teqc_link');
    
    $ws5->write(2,0, 'FASTQC:',$format{'0b'});
    $fqclink = $c{'winshare_root'} . $winpath . "\\" . "qc\\fastqc\\fastqc_report.html";
    $ws5->write_url(2, 1, $fqclink, $format, 'fastqc_link');
    
    $ws5->write(3,0, 'Commands log:',$format{'0b'});
    $fqclink = $c{'winshare_root'} . $winpath . "\\" . "qc\\bcbio-nextgen-commands.log";
    $ws5->write_url(3, 1, $fqclink, $format, 'commands_link');
    $diseasegenesdef_link = $diseasegenesdef_file . "_link";
    $ws5->write(4,0, 'Indications/genes:',$format{'0b'});
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
    print CFG "$sample,$sampleConfig{$sample}{gender},$sampleConfig{$sample}{regions},$sampleConfig{$sample}{disease}\n";
    close(CFG);

    $qc_from = $sampledir . "/qc";
    $qc_to = $sampledir . "/" . $sample . "/" . $rundate . "." . $nowtime;
    $bamindex = $bam . ".bai";

    system("cp $diseasegenesdef $qc_from");
    system("cp $summaryfile $qc_from"); 
    system("cp $programs $qc_from");
    
    system("cp -r $qc_from $qc_to");
    system("cp $bam $qc_to");
    system("cp $bamindex $qc_to");

}



