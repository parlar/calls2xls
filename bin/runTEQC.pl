#!/usr/bin/perl -w
#
use Statistics::R;
use IPC::Run;

die "run_TEQC.pl <bam file> <target file> <sample name> <out dir>\n" if (!(@ARGV));
die "run_TEQC.pl <bam file> <target file> <sample name> <out dir>\n" if ( $#ARGV != 3 );

chomp($b = $ARGV[0]);
chomp($t = $ARGV[1]);
chomp($s = $ARGV[2]);
chomp($destdir = $ARGV[3]);

@bamex = split(/\//, $b);
$bam_file = pop(@bamex);
$bam_path = join('/', @bamex);

@tarex = split(/\//, $t);
$target_file = pop(@tarex);
$target_path = join('/', @tarex);

$R = Statistics::R->new();
$R->run('library(TEQC)');
$R->run("readsfile <- file.path(\"$bam_path\", \"$bam_file\")");
$R->run("reads <- get.reads(readsfile, filetype=\"bam\")");
$R->run("targetfile <- file.path(\"$target_path\", \"$target_file\")");
$R->run('targets <- get.reads(targetfile, chrcol = 1, startcol = 2, endcol = 3)');
$R->run("TEQCreport(sampleName = \"$s\", targetsName = \"targets\", referenceName = \"hg19\", destDir = \"$destdir\",reads, targets, Offset = 0, pairedend = FALSE, genome = \"hg19\", k=c(1, 5, 10, 20, 30, 50, 100))");
$R->clean_up();
$R->stopR();
