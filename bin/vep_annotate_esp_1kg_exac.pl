#!/usr/bin/perl -w

die "vep_annotate_esp_1kg_exac.pl <vcf> \n" if (!(@ARGV));
die "vep_annotate_esp_1kg_exac.pl <vcf> \n" if ( $#ARGV != 0 );

chomp($vcf = $ARGV[0]);

$invcf = $vcf;
$outvcf =~ s/.vcf$/.d.1kg.esp.exac.vcf/g;

system("variant_effect_predictor.pl --cache --vcf -i $invcf -o $outvcf --refseq --offline --force_overwrite --everything --hgvs --fork 16 --fasta /home/data_in/.vep/homo_sapiens/80_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa");
