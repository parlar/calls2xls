#!/usr/bin/perl -w
#
# $Id
# $Revision$
# $Date$

use strict;
use warnings;
use vars qw/ $VERSION /;

use IPC::Run;
use List::Util qw(max min);
use Spreadsheet::WriteExcel;
use Time::localtime;
use Getopt::Long;
use Sort::Key::Natural qw(natsort);;
use Sort::Key::Maker int_nat_sort => qw(integer natural);
use Sort::Key::Maker chr_nat_sort => qw(natural);
use Config::General;
use YAML::XS 'LoadFile';
use MCE::Map
      max_workers =>16,             ## Default 'auto'
      chunk_size => 1              ## Default 'auto'
;


$VERSION = '0.1';

my %o;

$o{'in_silico'} = 0;
$o{'help'} = 0;

#mce_map sub, must be declared here

sub run_parallel{
	system($_);
}
	

## Read configs
my $conf = Config::General->new("$ENV{'CALLS2XLS'}/calls2xls.cfg");
my %c = $conf->getall;



########## get options and define run mode ##########

GetOptions (\%o, "help", "runfolder=s", "old_sample=s", "new_sample=s", "disease_hypothesis=s", "disease-gene-file=s", "ped-file=s");

if(!(@ARGV)){
    &usage0();
    exit 0;  
}

chomp(my $subcommand = $ARGV[0]);

if(($subcommand eq "single" || $subcommand eq "insilico") || $subcommand eq "family"){
    if($subcommand eq "single" && !exists($o{'runfolder'})){
		usage1();
		exit 0;
    }
    if($subcommand eq "insilico" && ((((!exists($o{'runfolder'}) || !exists($o{'old_sample'})) || !exists($o{'new_sample'})) || !exists($o{'disease_hypothesis'})) || !exists($o{'disease-gene-file'}))){

#		print "$subcommand $o{'runfolder'} $o{'old_sample'} $o{'new_sample'} $o{'disease_hypothesis'} $o{'disease-gene-file'}\n";

		usage2();
		exit 0;
    }
    if($subcommand eq "family" && ((!exists($o{'runfolder'}) || !exists($o{'disese_hypothesis'})) || !exists($o{'disese-gene-file'}) || !exists($o{'ped-file'}))){
		usage3();
		exit 0;
    }
}
else{
    usage0();
    exit 0;
}
if($o{'help'}){
    usage0();
    exit 0;
}


########## help messge subroutines  ##########

sub usage0{
	print <<"::_USAGE0_BLOCK_END_::";
	
usage: calls2xls.pl subcommand [options]
subcommands:         
  single     Annotate single single samples
  insilico   Transfer data to new sample number and annotate
  family     Analyse and annotate multiple samples using ped file 

::_USAGE0_BLOCK_END_::
	
	exit 1
}

sub usage1{
	print <<"::_USAGE1_BLOCK_END_::";

usage: calls2xls.pl single [options]
  --runfolder=path (required)
      Path to runfolder containing raw data and variant calls
      from bcbio-nextgen.    
  --config_file=file (optional)
      Path to config file with rows for sample IDs, disease hypothesis, 
      gender, disease-gene file. Use if configurations are not defined
      in SampleSheet.csv.
      Example row: 15-1234 male HCM disease-gene-file.txt 

::_USAGE1_BLOCK_END_::
	
	exit 1
} 


sub usage2{
	print <<"::_USAGE2_BLOCK_END_::";

usage: calls2xls.pl insilico [options]
  --runfolder=path (required)
    path to data root folder containing raw data and variant calls
    from bcbio-nextgen.    
  --old_sample=sample number, e.g. 15-1234
      old sample number to transfer calling data from.
  --new_sample=sample_number, e.g. 15-1235
      new sample number to transfer calling data to.
  --disease_hypothesis=hypothesis, e.g. HCM
      new disease indication.    
  --disease-gene-file=file
      path to disease-gene associations.

::_USAGE2_BLOCK_END_::
	
	exit 1

} 
sub usage3{
	print <<"::_USAGE3_BLOCK_END_::";

usage: calls2xls.pl family [options]\n\n";
  --runfolder=path (required)\n";
    path to funfolder containing raw data and variant calls \n";
    from bcbio-nextgen.\n";    
  --ped-file=file\n";
      ped file for inheritance annotations.\n";
  --disease-gene-file=file\n";
      disease-gene associations file\n\n";

::_USAGE3_BLOCK_END_::
	
	exit 1
} 

#colors
my %clr=(
    'red' => 10,
    'pink' => 45,
    'yellow' => 13,
    'lgreen'=> 42,
    'green' => 11
    );

my $finalpath       = $o{'runfolder'} . "/" . $c{'call_output_root'};
my $projectpath 	= &get_projectpath($finalpath);
my $summaryfile 	= $projectpath . '/' . $c{'callRunSummaryFile'};
my $programsfile 	= $projectpath . '/' . $c{'callRunProgramsFile'};
my $samplesheetfile = $o{'runfolder'} . "/SampleSheet.csv";
my $rundate 		= (split /_/, $o{'runfolder'})[0];
my $calldate 		= &get_calldate($summaryfile);
my $nowtime 		= &get_nowtime();
my @csq_heada 		= split /,/, $c{'CSQ'};
my %csq_headh;

#print "1. $c{'CSQ'}\n";
#print "2. @csq_heada\n";

for(my $csq_no = 0; $csq_no <= $#csq_heada; $csq_no++){
	$csq_headh{$csq_heada[$csq_no]} = $csq_no;

#    print "$csq_heada[$csq_no] $csq_no\n";
}





open(FH, "<$samplesheetfile");
chomp(my @samplesheet = <FH>);
close(FH);

my %sampleConfig;

for(my $i = 0; $i <= $#samplesheet; $i++){ 
	if($samplesheet[$i] =~ /^\[Data\]/){
		$i = $i+2; 						#actual samples start two rows below [Data]
		while($i <= $#samplesheet){
			my @sample_info = split(/,/, $samplesheet[$i]);
			my @description = split(/\$/, $sample_info[-1]);
			#description is last column: design$gender$hypothesis$genelist

			if(!exists($description[3])){
				die "no sample-specific data in samplesheet\n";
			}
			
			$sampleConfig{$sample_info[0]}{gender} = $description[1];
			$sampleConfig{$sample_info[0]}{design} = $description[0];
			$sampleConfig{$sample_info[0]}{disease} = $description[2];
			$sampleConfig{$sample_info[0]}{genelist} = $description[3];
			$sampleConfig{$sample_info[0]}{disease} =~ s/\s//g;
			$i++;
			
			
		}
	}
}

my @variant_fields 	= split(/,/, $c{'variant_field_order'});
my @vfhtmp 			= split(/,/, $c{'variant_header'});
my %vfh 			= map {split /:/} @vfhtmp; #variant variable name => column name  


my $variant_max_col = (split(/:/, (split(/,/, $c{'var_cols_widths'}))[-1]))[0];
my @csztmp 			= split(/,/, $c{'var_cols_widths'});

my %formdp;
my %formao;

foreach my $dp_tmp (split(/,/, $c{'dp'})){
	$formdp{$dp_tmp} = 1;
}

foreach my $ao_tmp (split(/,/, $c{'ao'})){
	$formao{$ao_tmp} = 1;
}

my %csz;
foreach my $row (@csztmp){ 			# variant sheet column widths
    my @kv = split(/:/, $row);
    my $k = $kv[0] . ":" . $kv[0];
    $csz{$k}=$kv[1];
}


my @samples_remove;

if($subcommand eq "insilico"){
#	print "hej\n";
    foreach my $field ( sort { $b cmp $a } keys(%{$sampleConfig{$o{'old_sample'}}})) {
		$sampleConfig{$o{'new_sample'}}{$field} = $sampleConfig{$o{'old_sample'}}{$field};
    }
    
    $sampleConfig{$o{'new_sample'}}{genelist}   = $o{'disease-gene-file'};
    $sampleConfig{$o{'new_sample'}}{disease}    = $o{'disease_hypothesis'};
    
    foreach my $sample ( sort { $b cmp $a } keys(%sampleConfig)) {
		if($sample ne $o{'old_sample'} && $sample ne $o{'new_sample'}){
			push @samples_remove, $sample;
		}
		if($sample ne $o{'new_sample'}){
			delete $sampleConfig{$sample};
		}
	}
    
    my $newsamppath = $finalpath . "/" . $o{'new_sample'};

#    print "$newsamppath\n";

    if(! -d $newsamppath){
		system "mkdir $newsamppath\n";
		system "cp -r $finalpath/$o{'old_sample'}/* $finalpath/$o{'new_sample'}";
		system "rm -r $finalpath/$o{'new_sample'}/$o{'old_sample'}";

    }
    my @sfiles = &get_sample_files_in_path($newsamppath,$o{'old_sample'});

	foreach my $file (@sfiles){
		my $nfile = $file;
		$nfile =~ s/^$o{'old_sample'}/$o{'new_sample'}/g;  

		$nfile  = $newsamppath . "/" . $nfile;
		$file   = $newsamppath . "/" . $file;

		rename($file,$nfile);
	}
}

# make sample dirs array

my @sampledirs;

foreach my $sample ( sort { $b cmp $a } keys(%sampleConfig)) {
    my $dir = $finalpath . "/" . $sample;
    if($subcommand eq "insilico"){
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
#print "@sampledirs\n";

# get alamut data
system("mut2vcf.pl $projectpath");
my $alamutfile = $projectpath . "/" . "alamut_mut.vcf.gz.vcf.gz";


#get programs
my @PROGS = &get_programs_used($projectpath);


#get calling.yaml config data
my $bcbio_config_path = $o{'runfolder'} . "/" . $c{'bcbio_config'};

open(FH, "<$bcbio_config_path") or die $!;
chomp(my @bcbio_config = <FH>);
close(FH) or die $!;

my %configHash;
my $cs = "";
for (my $i = 0; $i <= $#bcbio_config; $i++){
    if($bcbio_config[$i] =~ /  - files: /) {
		if($bcbio_config[$i+1] =~ /    description: (\S+)/){
			$cs = $1;
			$configHash{$cs}{$i} = $bcbio_config[$i];
		}
    }
    elsif($cs =~ /\S/) {
		$configHash{$cs}{$i} = $bcbio_config[$i];
    }
}

my %disease2genes;
my %sample2genes;


#get diseasegenesdef info
foreach my $sample ( sort { $b cmp $a } keys(%sampleConfig)) {
	my $dise;
#	print "$sample $sampleConfig{$sample}{disease}\n";

    my $diseasegenesdef = $c{'diseaseGeneAssocPath'} . "/" . $sampleConfig{$sample}{genelist};

#    print "$diseasegenesdef\n";


    open(FH0, "<$ENV{'CALLS2XLS'}/$diseasegenesdef");
    chomp(my @GeneL = <FH0>);
    close(FH0);
#	print "$sampleConfig{$sample}{disease}\n";
    foreach my $row (@GeneL) {
#	print "$row\n";
		if($row =~ /^>(\S+)/) {
			$dise = $1;
#			print "GeneL: $dise\n";
		}
		elsif($row =~ /\S/){
			my ($gene,$core,$vartype) = split(/\s+/, $row);
			$gene = uc($gene);
			if($dise eq $sampleConfig{$sample}{disease}){
#			    print "$dise\n";
				$disease2genes{$sample}{$sampleConfig{$sample}{disease}}{$gene}{'core'} = $core;
				$disease2genes{$sample}{$sampleConfig{$sample}{disease}}{$gene}{'vartype'} = $vartype;
				$sample2genes{$sample}{$gene} = 1;
			}	
		}
	}
}

#populate onTarget
my %onTarget;
my %onTargetReg;
my %onTargetGene;

foreach my $sample (sort { $a cmp $b } keys(%sampleConfig)) {
#	print "$sample\n";

    my $ccnf =  $ENV{'CALLS2XLS'} . "/" . $c{'designAndGeneListsRootPath'} . "/" .  $sampleConfig{$sample}{design};
    open(FH0, "<$ccnf");
    chomp(my @REG = <FH0>);
    close(FH0);

    foreach my $row (@REG) {
#	print "$sample $row\n";
		my ($chr, $p1, $p2, $gene) = split(/\t/, $row);
		$chr =~ s/\s//g; #chr
		$p1 =~ s/\s//g; #p1
		$p2 =~ s/\s//g; #p2
		$gene =~ s/\s//g; #gene
		$gene = uc($gene);
	
		if($gene !~ /chr/){
			if($chr eq "MT") {
				$chr = "chrM";
			}
			if($chr eq "X") {
				$chr = "chrX";
			}
			if($chr eq "Y") {
				$chr = "chrY";
			}
		}

		$onTargetGene{$sample}{$gene} = 1;
		
		if(exists($sample2genes{$sample}{$gene})){
#			print "$gene $sample $chr $p1 $p2 \n";
			$onTargetReg{$sample}{$chr}{$p1}{$p2} = $gene;
			$onTargetGene{$sample}{$gene} = 1;
			for(my $i = $p1 ; $i<= $p2; $i++) {
			    $onTarget{$sample}{$chr}{$i} = $gene;
			}
		}
    }
}

# run programs in parallel

my @cmd1;
my @cmd2;
my @cmd3;
my @cmd4;
my @cmd5;

foreach my $sampledir (@sampledirs) {
    my $sample 	= (split /\//, $sampledir )[-1];
    my $bam 	= $sampledir . "/" . $sample . "-ready.bam";
    my $reg 	= $ENV{'CALLS2XLS'} . "/" . $c{'designAndGeneListsRootPath'} . "/" .  $sampleConfig{$sample}{design};
    my $qcdir 	= $sampledir . "/qc";
    my $destdir = $sampledir . "/qc/teqc";
    my $sc      = $sampledir . "/" . "sample.coverage.bed";
    my $perbasec = $sampledir . "/" . "sample.perbase.coverage.bed";
    my $le30merge = $sampledir . "/" . "sample.panel.bg.coverage.le30merge.bed";

    if ( !-d $qcdir ) {
		system("mkdir $qcdir");
    }
    if ( !-d $destdir ) {
		system("mkdir -p $destdir");
    }

    push @cmd1, "$c{'bt'} coverage -d -b $bam -a $reg > $perbasec";
    push @cmd2, "$c{'bt'} genomecov -ibam $bam -bga  | $c{'bt'} intersect -wb -a stdin -b $reg | cut -f 1-4,8 | sed \'s/.*/&\t$sample\t$rundate\t$calldate\t$sampleConfig{$sample}{gender}/\' >$sc\n";
    push @cmd3, "filterCoverage.pl $sc | $c{'bt'} sort -i stdin | $c{'bt'} merge -i stdin -c 5,6,7,8 -o distinct,distinct,distinct,distinct > $le30merge\n";
    push @cmd4, "samtools index $bam\n";
    push @cmd5, "runTEQC.pl $bam $reg $sample $destdir\n";
}
print "running bamtools coverage\n";
#my @out1 = mce_map { run_parallel } @cmd1;

print "running bamtools genomecovrage\n";
#my @out2 = mce_map { run_parallel } @cmd2;

print "running filterCoverage\n";
#my @out3 = mce_map { run_parallel } @cmd3;

print "running samtools index\n";
#my @out4 = mce_map { run_parallel } @cmd4;

print "running TEQC\n";
#my @out5 = mce_map { run_parallel } @cmd5;

#annotation

my $vcfgz 			= $projectpath . "/batch1-ensemble.vcf.gz";
my $vcf 			= $projectpath . "/batch1-ensemble.vcf";
my $decompvcf 		= $projectpath . "/batch1-ensemble.d.vcf";
my $decompvcfanno 	= $decompvcf . ".ann.vcf";
my @vcf_annotations;

#create local annotation freqs db

system("putLocVarCovDBs.pl $o{'runfolder'}");
system("makeFreqDb.pl $c{'vardb'} $c{'covdb'} $c{'freqdb'}");

#annotate ensemble vcf

system("pigz -f -d -k $vcfgz");
system("vt decompose -s $vcfgz | vt normalize -r $c{'genome'} - | vt sort - | vt uniq - | pbgzip -c >$decompvcf");
system("$c{'vep'} --cache --cache_version $c{'ens-cache-version'} --vcf -i $decompvcf -o $decompvcfanno --refseq --offline --force_overwrite --everything --hgvs --fork $c{'noCores'} --fasta $c{'vep-fasta'} --assembly $c{'vep-assembly'} -custom $c{'freqdb'},LOC,vcf,exact,0,AF,AC,AN -custom $alamutfile,ALA,vcf,exact,0,ALAMUT,CTYPE,CLEVEL,CREATED,UPDATED");

if($subcommand eq "insilico") {
	
	my $oldsamps = join(' ', @samples_remove);

	my $decompvcfannonewfile = $decompvcfanno . "." . $o{'new_sample'} . ".vcf";
	my $decompvcfannonewtmp = $decompvcfanno . "." . $o{'new_sample'} . ".tmp";
	system("vcfremovesamples $decompvcfanno $oldsamps > $decompvcfannonewtmp");
	
    open(FH, "<$decompvcfannonewtmp") or die $!;
    open(FH2, ">$decompvcfannonewfile") or die $!;
    while (my $line = <FH>) {
		if($line =~ /^#CHROM/){
			$line =~ s/(\s)$o{'old_sample'}(\s)/$1$o{'new_sample'}$2/g;
			$line =~ s/\t+/\t/g;
		}
		print FH2 $line;
		push @vcf_annotations, $line;
	}
	system("rm $decompvcfannonewtmp");
}
else{
    print "decompvcfanno $decompvcfanno\n";
	open(FH, "<$decompvcfanno") or die $!;
	chomp(@vcf_annotations = <FH>);
	close(FH) or die $!;
}

my %field2col;
my %col2field;
my %calls;




foreach my $row (@vcf_annotations) {
	if($row =~ /^#CHROM/){
		print "$row\n";
		my @vcfrow = split(/\s+/, $row);
		for (my $fieldno = 0; $fieldno <= $#vcfrow; $fieldno++){
			$field2col{$vcfrow[$fieldno]} = $fieldno;
			$col2field{$fieldno} = $vcfrow[$fieldno];
			
#			print " $fieldno $vcfrow[$fieldno] = $fieldno;\n";
			
		}
	}
	elsif($row !~ /^#/ && $row =~ /\S/){
		my ($chr,$p,$id,$ref,$alt,$qual,$filter,$info,$form) = split(/\s+/, $row);

#		print "$info\n";

		my @rowa = split(/\s+/, $row);
		
		my @infoa = split(/;/, $info);
		my %infoh = ();
		my %csqh = ();
		my %formh = ();

		foreach my $in (@infoa) {
#		    print "$in\n";
			if($in =~ /=/){
				my ($key, $value) = split /=/, $in;
				$infoh{$key} = $value;
			}
		}
		
		# process CSQ fields

		$infoh{'CSQ'} =~ s/\|\|/\|.\|/g;
		$infoh{'CSQ'} =~ s/\|\|/\|.\|/g; #need to do twice..

		#structure format field
		
		my @f_tmp = split(/:/, $form);
		for (my $f_no = 0; $f_no <= $#f_tmp; $f_no++){
			$formh{$f_tmp[$f_no]} = $f_no;
		}
		

		my $al_ctype = "NA";
		my $al_clevel = 0;
		my $al_create = "NA";
		my $al_update = "NA";
		my $loc_freq = "NA";
		my $loc_an = "NA";
		my $loc_ac = "NA";
		my $al_alamut = "NA";

		if( exists $infoh{'ALA_ALAMUT'}){
			$al_alamut = $infoh{'ALA_ALAMUT'};
		}
		
		if( exists $infoh{'ALA_CTYPE'}){
			if($infoh{'ALA_CTYPE'} eq $c{'alamut_ctype'}){
				$al_ctype = $infoh{'ALA_CTYPE'};
			}
		}
		
		if( exists $infoh{'ALA_CLEVEL'}){
			if($infoh{'ALA_CTYPE'} eq $c{'alamut_ctype'}){
				$al_clevel = $infoh{'ALA_CLEVEL'};
			}
		}
		
		if( exists $infoh{'ALA_CREATED'}){
			$al_create = $infoh{'ALA_CREATED'};
		}
		
		if( exists $infoh{'ALA_UPDATED'}){
			$al_update = $infoh{'ALA_UPDATED'};
		}
		
		if( exists $infoh{'LOC_AF'}){
			$loc_freq = $infoh{'LOC_AF'};
		}
		
		if( exists $infoh{'LOC_AN'}){
			$loc_an = $infoh{'LOC_AN'};
		}
		
		if( exists $infoh{'LOC_AC'}){
			$loc_ac = $infoh{'LOC_AC'};
		}
		
		
		foreach my $sample (sort { $a cmp $b } keys(%sampleConfig)) {
		    print "$sample\t";
			my @sample_gt_data = split(/:/, $rowa[$field2col{$sample}]);

			if($sample_gt_data[$formh{GT}] =~ /[1-9]/ && exists($onTarget{$sample}{$chr}{$p})) {
				foreach my $csq (split(/,/,$infoh{'CSQ'})) {
#				    print "$csq\n";
					my @csq_fields = split(/\|/, $csq);
					
					my $variant = $chr . "_" . $p . "_" . $ref . "_" . $alt . "_" . $csq;

					print "$variant\n";

					$calls{$sample}{$variant}{'gt'} =  $sample_gt_data[$formh{GT}];
					$calls{$sample}{$variant}{'dp'} = "NA";
					$calls{$sample}{$variant}{'ao'} = "NA";

					foreach my $dp (keys(%formdp)){
						if(exists ($formh{$dp})){
							$calls{$sample}{$variant}{'dp'} = $sample_gt_data[ $formh{$dp} ];
						}
					}
					foreach my $ao (keys(%formao)){
						if(exists ($formh{$ao})){
							$calls{$sample}{$variant}{'ao'} = $sample_gt_data[ $formh{$ao} ];
						}
					}	 

					
					$csq_fields[$csq_headh{EA_MAF}] =~ s/[ACGT:-]+//g;
					$csq_fields[$csq_headh{AA_MAF}] =~ s/[ACGT:-]+//g;
					$csq_fields[$csq_headh{GMAF}] =~ s/[ACGT:-]+//g;

					$calls{$sample}{$variant}{'chr'} = $chr;
					$calls{$sample}{$variant}{'pos'} = $p;
					$calls{$sample}{$variant}{'ref'} = $ref;
					$calls{$sample}{$variant}{'alt'} = $alt;
					$calls{$sample}{$variant}{'e.ea'} = $csq_fields[$csq_headh{EA_MAF}];
					$calls{$sample}{$variant}{'e.aa'} = $csq_fields[$csq_headh{AA_MAF}];
					$calls{$sample}{$variant}{'gmaf'} = $csq_fields[$csq_headh{GMAF}];
					$calls{$sample}{$variant}{'type'} = $csq_fields[$csq_headh{Consequence}];
					$calls{$sample}{$variant}{'impact'} = $csq_fields[$csq_headh{IMPACT}];
					$calls{$sample}{$variant}{'symbol'} = $csq_fields[$csq_headh{SYMBOL}];
					$calls{$sample}{$variant}{'feature'} = $csq_fields[$csq_headh{Feature}];
					$calls{$sample}{$variant}{'hgvsc'} = $csq_fields[$csq_headh{HGVSc}];
					$calls{$sample}{$variant}{'hgvsp'} = $csq_fields[$csq_headh{HGVSp}];
					$calls{$sample}{$variant}{'a.ctype'} = $al_ctype;
					$calls{$sample}{$variant}{'a.clevel'} = $al_clevel;
					$calls{$sample}{$variant}{'a.create'} = $al_create;
					$calls{$sample}{$variant}{'a.update'} = $al_update;
					$calls{$sample}{$variant}{'l.af'} = $loc_freq;
					$calls{$sample}{$variant}{'l.ac'} = $loc_ac;
					$calls{$sample}{$variant}{'l.an'} = $loc_an;					
					$calls{$sample}{$variant}{'tsymbol'} = $onTarget{$sample}{$calls{$sample}{$variant}{'chr'}}{$calls{$sample}{$variant}{'pos'}};

#					print "$sample $variant\n";

				}
			}
		}
    }
}

foreach my $sampledir (@sampledirs) {
    my $sample 	= (split /\//, $sampledir )[-1];
    my $sampleyear = "20" . substr($sample, 0, 2);
    my $disease 	= $sampleConfig{$sample}{disease};
    my $bam_file 	= $sample . "-ready.bam";
    my $bam 		= $sampledir . "/" . $sample . "-ready.bam";
    my $bamindex 	= $bam . ".bai";
    my $sampleout_sampleyear_rundate_folder = $sampledir . "/" . $sampleyear . "/" . $sample . "/" . $rundate . "." . $nowtime;
    my $winbampath                          = $c{'winshare_root'} . $sampleyear . "\\" . $sample . "\\" .  $rundate . ".$nowtime";
    my $xls_name 	= $sampleout_sampleyear_rundate_folder . "/" . $sample . "_calldata.xls";      
    my $winpath 	= $sampleyear . "\\" . $sample . "\\" . $rundate . "." . $nowtime;
	my $disease_genes_file = $ENV{'CALLS2XLS'} . "/diseasegenelists/" . $sampleConfig{$sample}{genelist};
    my $qc_from 	= $sampledir . "/qc";
    my $qc_to 		= $sampledir . "/" . $sampleyear . "/" . $sample . "/" . $rundate . "." . $nowtime;

#    print "$sample\n";

    if(!-d $sampleout_sampleyear_rundate_folder ){
		system("mkdir -p $sampleout_sampleyear_rundate_folder");
    }

	# create workbook!
 
    my $wb 	= Spreadsheet::WriteExcel->new($xls_name);

	# make bg formats for workbook

	my %format;

    $format{'5'} 	= $wb->add_format(bg_color=>$clr{'red'},border   => 1);
    $format{'4'} 	= $wb->add_format(bg_color=>$clr{'pink'},border   => 1);
    $format{'3'} 	= $wb->add_format(bg_color=>$clr{'yellow'},border   => 1);
    $format{'2'} 	= $wb->add_format(bg_color=>$clr{'lgreen'},border   => 1);
    $format{'1'} 	= $wb->add_format(bg_color=>$clr{'green'},border   => 1);
    $format{'0'} 	= $wb->add_format(bg_color=> 9,border   => 1, );
    $format{'0b'} 	= $wb->add_format(bg_color=> 9,border   => 1, );
    $format{'0b'} 	-> set_bold();

	# create variants worksheet

    my $ws_variant 	= $wb->add_worksheet('calldata');
    
	# set column widths for variants worksheet
    
	foreach my $key (keys %csz){
		$ws_variant->set_column($key, $csz{$key});
    }

    $ws_variant->freeze_panes(1, 0);
    
    my @indexorder = split(/,/, $c{'alamut_clevel_sort_order'});
    my %sortmap 	= map { $indexorder[$_] => $_ } 0 .. $#indexorder;

    my $wsrowno = 0;

    #populate header

    $ws_variant->write_row($wsrowno,0, \@variant_fields, $format{'0b'});
    $wsrowno++;
    $ws_variant->write_row($wsrowno,0, \@variant_fields, $format{'0'});

    my $autofstart = $wsrowno + 1;
    my $autofilt 	= "A" . $autofstart;
    $wsrowno++;

    my %uniqvariants = ();

    my $varno = 0;
	my $genpos;
    foreach my $variant (int_nat_sort {   $sortmap{ $calls{$sample}{$_}{'a.clevel'} }, $_   } keys %{$calls{$sample}}){
        print "$sample $variant\n";
		$genpos = (split(/\|/, $variant))[0];
			
		if (!exists($uniqvariants{$genpos})){
			$uniqvariants{$genpos} = 1;
			$varno++;
		}

		my @outrow = ($varno);
		foreach my $v (@variant_fields){
            if ($v ne "vno") {
               push @outrow, $calls{$sample}{$variant}{$v}; 
            }
 		}

		my $mlink = $c{'mlink'} . $calls{$sample}{$variant}{'chr'} . ":" . $calls{$sample}{$variant}{'pos'} . "-" . $calls{$sample}{$variant}{'pos'};
		my $blink = $c{'blink'} . $winbampath . "\\" . $bam_file;

		$ws_variant->write_row($wsrowno, 0, \@outrow, $format{$calls{$sample}{$variant}{'a.clevel'}});
		$ws_variant->write_url($wsrowno, $#outrow+1, $mlink, $format{$calls{$sample}{$variant}{'a.clevel'}}, 'mut');
		$ws_variant->write_url($wsrowno, $#outrow+2, $blink, $format{$calls{$sample}{$variant}{'a.clevel'}}, 'bam');
		$wsrowno++;
    }        
	$autofilt .= ":" . $variant_max_col . $#variant_fields;
	$ws_variant->autofilter($autofilt);

	$wsrowno+=3;
		
	$ws_variant->write_blank($wsrowno, 1, $format{"5"});
	$ws_variant->write($wsrowno, 2, 'Class 5, Pathogenic', $format{"0"});

	$wsrowno++;

	$ws_variant->write_blank($wsrowno, 1, $format{"4"});
	$ws_variant->write($wsrowno, 2, 'Class 4, Likely pathogenic');

	$wsrowno++;
	
	$ws_variant->write_blank($wsrowno, 1, $format{"3"});
	$ws_variant->write($wsrowno, 2, 'Class 3, VUS');
 
	$wsrowno++;

    $ws_variant->write_blank($wsrowno, 1, $format{"2"});
    $ws_variant->write($wsrowno, 2, 'Class 2, Likely benign');
 
    $wsrowno++;

    $ws_variant->write_blank($wsrowno, 1, $format{"1"});
    $ws_variant->write($wsrowno, 2, 'Class 1, Benign');

    my $le30merge = $sampledir . "/" . "sample.panel.bg.coverage.le30merge.bed";


    open(FH0, "<$le30merge");
    chomp(my @CLE30 = <FH0>);
    close(FH0);
    
    @CLE30 = natsort(@CLE30);
    
    print "$le30merge\n";

    my $ws_lowcov = $wb->add_worksheet('cov<30');
    $ws_lowcov->write(0,0, 'chrom', $format{'0b'});
    $ws_lowcov->write(0,1, 'start', $format{'0b'});
    $ws_lowcov->write(0,2, 'end', $format{'0b'});
    $ws_lowcov->write(0,3, 'size', $format{'0b'});	
    $ws_lowcov->write(0,4, 'gene', $format{'0b'});
    $ws_lowcov->write(0,5, 'core', $format{'0b'});
    $ws_lowcov->write(0,6, 'reg_link', $format{'0b'});
    $ws_lowcov->write(0,7, 'bam_link', $format{'0b'});
    
    my $lowcovno = 0;

	

    foreach my $row (@CLE30) {
		my ($chr, $start, $end, $gene, $sampl) = split(/\t/, $row);
		my $size = $end - $start;
		my @covprint = ();

		if (exists($disease2genes{$sample}{$disease}{$gene})) {
			@covprint = ($chr,$start,$end,$size,$gene);
			push @covprint, $disease2genes{$sample}{$disease}{$gene}{'core'} ;

			my $covformat = $format{'0'};
			if($disease2genes{$sample}{$disease}{$gene}{'core'} eq "+") {
				$covformat = $format{'0b'};
			}
			
			$ws_lowcov->write_row($lowcovno+1,0, \@covprint, $covformat);
	    
			my $mlink = $c{'mlink'} . $covprint[0] . ":" . $start . "-" . $end;
			my $blink = $c{'blink'} . $c{'winshare_root'} . $winpath . "\\". $bam_file;	
			
			$ws_lowcov->write_url($lowcovno+1, 6, $mlink, $covformat, 'reg_link');
			$ws_lowcov->write_url($lowcovno+1, 7, $blink, $covformat, 'bam_link');
			$lowcovno++;
		}
    }
    
    my $regcovfile =  $sampledir . "/" . "sample.perbase.coverage.bed";

    open(FH0, "<$regcovfile");
    chomp(my @RCF = <FH0>);
    close(FH0);
    

    my %genecov = ();
    my %allcov = ();
    my %genecov_30 = ();	

    foreach my $row (@RCF){
		my @tmp = split(/\t/, $row);
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
    my $ws_genes = $wb->add_worksheet('selected genes');
    
    $ws_genes->write(0,0, 'genes in disease hypothesis file', $format{'0b'});
    $ws_genes->write(0,1, 'disease hypothesis', $format{'0b'});
    $ws_genes->write(0,2, 'exists on panel', $format{'0b'});
    $ws_genes->write(0,3, 'core', $format{'0b'});
    $ws_genes->write(0,4, 'vartype', $format{'0b'});    
    
    my $selgenesno = 1;
    foreach my  $gene (sort { $a cmp $b } keys(%{$disease2genes{$sample}{$disease}})){
		print "$gene\n";
		my $selectformat = $format{'0'};


		if($disease2genes{$sample}{$disease}{$gene}{core} eq "+") {
			$selectformat = $format{'0b'};
		}
					
		$gene =~ s/\s+//g;
		$ws_genes->write($selgenesno,0, $gene, $selectformat);
		$ws_genes->write($selgenesno,1, $sampleConfig{$sample}{disease}, $selectformat);
		
		if(exists($onTargetGene{$sample}{$gene})) {
			$ws_genes->write($selgenesno,2, "yes", $selectformat);
		}
		else{
			$ws_genes->write($selgenesno,2, "no", $selectformat);
		}
		if($disease2genes{$sample}{$disease}{$gene}{'core'} eq "+") {
			$ws_genes->write($selgenesno,3, "+", $selectformat);
		}
		else{
			$ws_genes->write($selgenesno,3, "-", $selectformat);
		}
		$ws_genes->write($selgenesno,4, $disease2genes{$sample}{$disease}{$gene}{'vartype'}, $selectformat);
		$selgenesno++;
    }
    
    my $ws_completeness = $wb->add_worksheet('genes coverage & completeness');

    $ws_completeness->write(0,0, 'chrom',$format{'0b'});
    $ws_completeness->write(0,1, 'start',$format{'0b'});
    $ws_completeness->write(0,2, 'end',$format{'0b'});
    $ws_completeness->write(0,3, 'gene',$format{'0b'});
    $ws_completeness->write(0,4, 'core',$format{'0b'});
    $ws_completeness->write(0,5, 'av_cov',$format{'0b'});
    $ws_completeness->write(0,6, 'completeness_30x',$format{'0b'});	
    
    my $completenessno = 1;
    
    foreach my $gene (sort { $a cmp $b } keys (%genecov)) {
		$gene =~ s/\s+//g;
		if (exists($disease2genes{$sample}{$disease}{$gene})) {
			
			my $completenessformat = $format{'0'};
			
			if($disease2genes{$sample}{$disease}{$gene}{'core'} eq "+") {
				$completenessformat = $format{'0b'};
			}

			my $chr = $genecov{$gene}{chromosome};
			my $max = max @{$genecov{$gene}{positions}};
			my $min = min @{$genecov{$gene}{positions}};
			my $avcov = sprintf("%.2f", $genecov{$gene}{addcov}/$genecov{$gene}{rows});
			
			if(!exists($genecov_30{$gene}{rows})){
				$genecov_30{$gene}{rows} = 0;
			}
			my $completeness = sprintf("%.2f", $genecov_30{$gene}{rows}/$genecov{$gene}{rows}*100) . "%";
			
			my @completenessprint = ($chr, $max, $min, $gene, $disease2genes{$sample}{$disease}{$gene}{'core'}, $avcov, $completeness);
			
			$ws_completeness->write_row($completenessno, 0, \@completenessprint, $completenessformat);
			
			$completenessno++;
		}
    }
    
    my $destdir = $sampledir . "/qc/teqc";
    
    if ( !-d $destdir ) {
		system("mkdir $destdir");
    }
    
    my $ws_quality = $wb->add_worksheet('quality');
    $ws_quality->write(1,0, 'TEQC:', $format{'0b'});

    my $teqc_link = $c{'winshare_root'} . $winpath . "\\" . "qc\\teqc\\index.html";
    $ws_quality->write_url(1, 1, $teqc_link, $format{'0b'}, 'teqc_link');
    
    $ws_quality->write(2,0, 'FASTQC:',$format{'0b'});
    my $fqc_report = $c{'winshare_root'} . $winpath . "\\" . "qc\\fastqc\\fastqc_report.html";
    $ws_quality->write_url(2, 1, $fqc_report, $format{'0b'}, 'fastqc_link');
    
    $ws_quality->write(3,0, 'Commands log:',$format{'0b'});

    my $bcbio_link = $c{'winshare_root'} . $winpath . "\\" . "qc\\bcbio-nextgen-commands.log";
    $ws_quality->write_url(3, 1, $bcbio_link, $format{'0b'}, 'commands_link');

    my $diseasegenesdef_link = $sampleConfig{$sample}{genelist} . "_link";
    $ws_quality->write(4,0, 'Indications/genes:',$format{'0b'});

    my $disease_link = $c{'winshare_root'} . $winpath . "\\" . "qc\\" . $diseasegenesdef_link;
    $ws_quality->write_url(4, 1, $disease_link, $format{'0b'}, $sampleConfig{$sample}{genelist});
    
    my $bamtoolsfile = $sampledir . "/qc/bamtools/bamtools_stats.txt";
    
    open (FH1,"<$bamtoolsfile");
    chomp(my @BT = <FH1>);
    close(FH1);
    
    for (my $bti = 0; $bti <= $#BT; $bti++) {
		$ws_quality->write($bti+5,0,$BT[$bti]);
    }

    my $ws_bcbio_config = $wb->add_worksheet('bcbio config parameters');
    my $configrow = 0;
    foreach my $row (sort { $a <=> $b } keys(%{$configHash{$sample}})) {
		$ws_bcbio_config->write($configrow,0,$configHash{$sample}{$row});
		$configrow++;
    }

    my $ws_programs = $wb->add_worksheet('programs');
#	print @PROGS;	
    my $progrow = 0;
    foreach my $row (@PROGS) {
		$row =~ s/\s//g;
		$ws_programs->write($progrow,0,$row);
		$progrow++;
    }

    my $ws_design_regions = $wb->add_worksheet('analyzed regions');
	my $design_reg_row = 0;
#	print "$sample\n";
    foreach my $chr (chr_nat_sort {$_} keys(%{$onTargetReg{$sample}})) {
#		print "$chr\n";
		foreach my $p1 (sort { $a <=> $b } keys(%{$onTargetReg{$sample}{$chr}})) {
			foreach my $p2 (sort { $a <=> $b } keys(%{$onTargetReg{$sample}{$chr}{$p1}})) {
				my @regprint = "$chr,$p1,$p2,$onTargetReg{$sample}{$chr}{$p1}{$p2}";
				$ws_design_regions->write_row($design_reg_row,0,\@regprint);
				$design_reg_row++;
			}
		}
	}
    my $configdatfile = $sampledir . "/qc/sample_config_file.csv";

    open(CFG, ">$configdatfile");
		print CFG "#sample,gender,regions,disease\n";
		print CFG "$sample,$sampleConfig{$sample}{gender},$sampleConfig{$sample}{design},$sampleConfig{$sample}{disease}\n";
    close(CFG);

    system("cp $disease_genes_file $qc_from");
    system("cp $summaryfile $qc_from"); 
    system("cp $programsfile $qc_from");
    
    system("cp -r $qc_from $qc_to");
    system("cp $bam $qc_to");
    system("cp $bamindex $qc_to");

}

# subroutines

sub get_projectpath(){
	my $finalpath = shift;
	my @projectdirs;
	opendir(DIR, $finalpath);
	@projectdirs = grep { /project_/ } readdir(DIR);
	close(DIR);
	
#	print "@projectdirs\n";
	return $finalpath . "/" . $projectdirs[0];
}

sub get_sample_files_in_path(){
	my ($path,$sample) = @_;

	print "$path\n";
	print "$sample\n";

	my @files;
	opendir(DIR, $path);
	@files = grep {/$sample/} readdir(DIR);
	close(DIR);
	
#	print "@projectdirs\n";
	return @files;
}

sub get_calldate(){
	my $sfile = shift;
	my $config = LoadFile($sfile);
	my $dt_tmp = $config -> {date};
	$dt_tmp =~ /\d\d(\d\d)-(\d\d)-(\d\d)\s/;
	my $calldate = $1 . $2 . $3;
	return $calldate;
}

sub get_nowtime(){
	my $tm = localtime;
	return sprintf("%04d.%02d.%02d-%02d.%02d.%02d", $tm->year+1900, ($tm->mon)+1, $tm->mday, $tm->hour, $tm->min, $tm->sec );
}

sub get_programs_used(){
	my $in = shift;
	my $bcbio_programs = $in . "/programs.txt";
#	print "$bcbio_programs\n";
	my @PROGS;
	push @PROGS, "bcbio_nextgen:";
	open(FH, "<$bcbio_programs") or die $!;
	while(<FH>){
		push @PROGS, $_;
	}
	close(FH) or die $!;

	push @PROGS, "not part of bcbio:";

	return @PROGS;

}
