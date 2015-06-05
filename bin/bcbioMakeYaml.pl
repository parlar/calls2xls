#!/usr/bin/perl -w
#
#       split a multifasta file to separate fasta files
#
#

use Cwd;

die "usage: bcbio_make_yaml_germline_v3.pl <fastq folder>\n" if (!exists($ARGV[0]));
die "usage: bcbio_make_yaml_germline_v3.pl <fastq folder>\n" if ( $#ARGV < 0 );

chomp($folder = $ARGV[0]);
$folder =~ s/\/$//g;
$root = "/home/data_in/calling/haloplex";
$design_genelistsRoot = $root . "/design_genelists";

@fpath = split(/\//, $folder);
@folder_explode = split(/_/, $fpath[-1]);

#$rundate = $folder_explode[0];

$sample_file = $folder . "/SampleSheet.csv";

open(FH, "<$sample_file");
chomp(@SAMPLESHEET = <FH>);
close(FH);

$getdata = 0;
$getsample = 0;

for ($i = 0; $i <= $#SAMPLESHEET; $i++) {
    if($SAMPLESHEET[$i] =~ /^\[Data\]/ ){
	$getdata = 1;
	next;
    }
    elsif($getdata && $SAMPLESHEET[$i] =~ /^Sample_ID/){
	$getsample = 1;
	next;
    }
    if($getsample && $getdata) {
	if($SAMPLESHEET[$i] =~ /\S/){
	    @tmp = split(/,/, $SAMPLESHEET[$i]);
	    @tmp2 = split(/\$/, $tmp[7]);
	    
	    $sdat{$tmp[0]}{'regions_design'} = $design_genelistsRoot . "/" . $tmp2[0];
	    $sdat{$tmp[0]}{'sex'} = $tmp2[1];
	    $sdat{$tmp[0]}{'indication'} = $tmp2[2];
	    $sdat{$tmp[0]}{'genedisease'} = $design_genelistsRoot . "/" . $tmp2[3];
	    
	    print "$tmp[0]\t$sdat{$tmp[0]}{'sex'}\t$sdat{$tmp[0]}{'regions_design'}\t$sdat{$tmp[0]}{'indication'}\n";
	}
    }
}

system("make_fastq_links.pl $folder");

$folder =~ s/\/$//;
$fastq_lns = $folder . "/" . "fastqs-lns";
$calldir = $folder . "/" . "calling";
$confdir = $calldir . "/" . "config";
$workdir = $calldir . "/" . "work";

system("mkdir $calldir");
system("mkdir $confdir");
system("mkdir $workdir");


unless (open(FH0, "find $fastq_lns -name \"*.fq.gz\" 2>&1 |")) {
        die "Can't spawn external command!";
}       
chomp(@GZ = <FH0>);
unless (close(FH0)) {  
        die "External command failed: $?";
}

$pwd = getcwd;
foreach $file (@GZ) {
    $file =~ s/^\.\///g;
    
    print "$file\n";
    @tmp = split(/\//, $file);
    @sample = split(/_/, $tmp[-1]);
    $file = $pwd . "/" . $file;
    $file =~ s/\/\///g;
    $data{$sample[0]}{$file} = 1;
}

$outfile = $confdir . "/" . "calling.yaml";

## Make yaml
open(FH, ">$outfile");
print FH "---\n";
print FH "upload:\n";
print FH "  dir: ../final\n";
print FH "details:\n";

foreach $sample (sort { $a cmp $b } keys(%sdat) ) {
#    print "$sample\n";
    
    $sex = $sdat{$sample}{sex};
    $reg = $sdat{$sample}{regions_design};
    $f1 = "";
    $f2 = "";
    foreach $file (sort { $a cmp $b } keys(%{$data{$sample}})) {
	if($file =~ /_1.fq.gz/) {
	    $f1 = $file;
	}
	if($file =~ /_2.fq.gz/) {
	    $f2 = $file;
	}
    }
    print FH "  - files: [$f1, $f2]\n";
    print FH "    description: $sample\n";
    print FH "    metadata:\n";
    print FH "      batch: batch1\n";    
    print FH "      sex: $sex\n";
    print FH "    analysis: variant2\n";
    print FH "    genome_build: hg19\n";
    print FH "    algorithm:\n";
    print FH "      aligner: bwa\n";
#    print FH "      align_split_size: 5000000\n";
    print FH "      mark_duplicates: false\n";
    print FH "      variantcaller: [freebayes, gatk-haplotype, platypus, samtools]\n";
#    print FH "      jointcaller: [freebayes-joint, gatk-haplotype-joint, platypus-joint, samtools-joint]\n";
    print FH "      quality_format: Standard\n";
    print FH "      coverage_interval: amplicon\n";
    print FH "      recalibrate: false\n";
    print FH "      realign: gatk\n";
    print FH "      merge_bamprep: true\n";
    print FH "      variant_regions: $reg\n";
    print FH "      clinical_reporting: true\n";
    print FH "      effects: false\n";
    print FH "      ensemble:\n";
    print FH "        numpass: 1\n";
#    print FH "        format-filters: [DP < 4]\n";
#    print FH "        classifiers:\n"; 
#    print FH "          balance: [AD, FS, Entropy]\n";
#    print FH "          calling: [ReadPosEndDist, PL, PLratio, Entropy, NBQ]\n";
#    print FH "        classifier-params:\n";
#    print FH "          type: svm\n";
#    print FH "        trusted-pct: 0.65\n";
}

close(FH);

chdir "$workdir";
system("bcbio_nextgen.py -t local -n 30 ../config/calling.yaml");
chdir "../../";
