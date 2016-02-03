#!/usr/bin/perl -w
#
#       Perform fastq trimming, create bcbio-nextgen folder structure for variant calling
#
#

use Config::General;
use Cwd 'abs_path';
use Time::localtime;

die "usage: bcbioMakeYaml_validate.pl <runfolder>\n" if (!exists($ARGV[0]));
die "usage: bcbioMakeYaml_validate.pl <runfolder>\n" if ( $#ARGV != 0 );

chomp($runfolder = $ARGV[0]);
$runfolder =~ s/\/$//g;

my $nowtime = &get_nowtime();

## Read configs
$conf = Config::General->new("$ENV{'CALLS2XLS'}/calls2xls.cfg");
%c = $conf->getall;

$ssfile = $runfolder . "/SampleSheet.csv";

print "$ssfile\n";

open(FH, "<$ssfile");
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
	    @sample_desc = split(/,/, $SAMPLESHEET[$i]);
	    @callspec = split(/\$/, $sample_desc[7]);
	    $sdat{$sample_desc[0]}{'designtarget'} = $ENV{'CALLS2XLS'} . "/" . $c{'designFilesPath'} . "/" . $callspec[0];
	    $sdat{$sample_desc[0]}{'gender'} = $callspec[1];
	}
    }
}

$runfolder =~ s/\/$//;
$fastq_lns = $runfolder . "/" . "fastqs-lns_" . $nowtime;
$calldir = $runfolder . "/" . "calling_" . $nowtime;
$confdir = $calldir . "/" . "config";
$workdir = $calldir . "/" . "work";

system("mkdir -p $calldir");
system("mkdir -p $confdir");
system("mkdir -p $workdir");

system("makeFastqLinks_validate.pl $runfolder $nowtime");

chomp(@GZ = <$fastq_lns/*.fq.gz>);

#print "@GZ\n";

foreach $fqgz (@GZ) {
    $fqgz =~ s/^\.\///g;
    $fqgz_file = (split(/\//, $fqgz))[-1];
    $sample = (split(/_/, $fqgz_file))[0];
    $data{$sample}{$fqgz} = 1;
}

$outfile = $confdir . "/" . "calling.yaml";

$vcfdir = $runfolder . "/vcfs";

@refvcfs =  <$vcfdir/*.vcf>;

## Make yaml
open(FH, ">$outfile");
print FH "---\n";
print FH "upload:\n";
print FH "  dir: ../final\n";
print FH "details:\n";


foreach $sample (sort { $a cmp $b } keys(%sdat) ) {
    $f1 = "";
    $f2 = "";
    $validate_set = "";

    foreach $file (sort { $a cmp $b } keys(%{$data{$sample}})) {
	if($file =~ /_1.fq.gz/) {
	    $f1 = abs_path($file);
	}
	if($file =~ /_2.fq.gz/) {
	    $f2 = abs_path($file);
	}
    }
    foreach $vcf (@refvcfs){
	print "$vcf\n";
	if($vcf =~ /\/$sample-/){
	    $validate_set = abs_path($vcf);
	}
    }
    
    print FH "  - files: [$f1, $f2]\n";
    print FH "    description: $sample\n";
    print FH "    metadata:\n";
    print FH "      batch: batch1\n";    
    print FH "      sex: $sdat{$sample}{'gender'}\n";
    print FH "    analysis: variant2\n";
    print FH "    genome_build: hg19\n";
    print FH "    algorithm:\n";
    print FH "      aligner: bwa\n";
    print FH "      mark_duplicates: false\n";
    print FH "      variantcaller: [freebayes, gatk-haplotype, platypus, samtools]\n";
    print FH "      svcaller: [lumpy, manta, cnvkit]\n";
    print FH "      quality_format: Standard\n";
    print FH "      coverage_interval: amplicon\n";
    print FH "      recalibrate: false\n";
    print FH "      realign: gatk\n";
    print FH "      merge_bamprep: true\n";
    print FH "      variant_regions: $sdat{$sample}{designtarget}\n";
    print FH "      clinical_reporting: true\n";
    print FH "      effects: false\n";
    print FH "      ensemble:\n";
    print FH "        numpass: 1\n";
    print FH "      validate: $validate_set\n";
    print FH "      validate_regions: $sdat{$sample}{designtarget}\n";
    print FH "      validate_method: rtg\n";

}

close(FH);

sub get_nowtime(){
	my $tm = localtime;
	return sprintf("%04d.%02d.%02d-%02d.%02d.%02d", $tm->year+1900, ($tm->mon)+1, $tm->mday, $tm->hour, $tm->min, $tm->sec );
}
